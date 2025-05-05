import torch

def avg_dist_mat(grid):
    """
    Calculates the average distance between each vector in the grid and all of its 8 neighbors (including diagonals).
    Uses a 3x3 kernel with circular padding.

    WARNING! Restriction:
        This function assumes a symmetric distance measure. It computes distances using a normalization trick that is only valid if the distance from x to y is the same as from y to x.

    Args:
        grid (torch.Tensor): 4D tensor with shape (batch, channels, height, width)

    Returns:
        torch.Tensor: 3D tensor with shape (batch, height, width) containing the average neighbor distances for each vector.
    """
    b, c, h, w = grid.shape
    padded = torch.nn.functional.pad(grid, (1, 1, 1, 1), mode="circular")
    unf = torch.nn.functional.unfold(
        padded, kernel_size=(3, 3), stride=1, padding=0
    ).reshape(b, c, 3, 3, h, w)
    # Compute L2 distance between center and all neighbors
    dist = torch.linalg.norm(
        unf - unf[:, :, 1, 1, :, :].unsqueeze(2).unsqueeze(2), dim=1
    )
    # Average over 8 neighbors (excluding the center)
    avg_dist = dist.reshape(b, 3 * 3 * 1, h * w) / 8.0
    return torch.nn.functional.fold(
        avg_dist, output_size=(h + 2, w + 2), kernel_size=(3, 3), stride=1
    )[:, 0, :, :]

def sum_avg_dist(grid):
    """
    Computes the sum of the average neighbor distances for each vector in the grid.

    Args:
        grid (torch.Tensor): 4D tensor with shape (batch, channels, height, width)

    Returns:
        torch.Tensor: 1D tensor with shape (batch,) containing the sum of average neighbor distances for each batch.
    """
    avg_dist = avg_dist_mat(grid)
    sum_avg_dist = avg_dist.sum(dim=(1, 2)).reshape(-1,)
    return sum_avg_dist

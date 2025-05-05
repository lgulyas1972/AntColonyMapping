import torch
import re
import pandas as pd
from io import StringIO
from glob import glob
import math
import json

def load_netlogo_state(path):
    """
    Load a NetLogo state file and extract the node vectors as a torch tensor.

    The function expects the data between "TURTLES" and "EXTENSIONS" to be a CSV.
    Only rows with breed '{breed nodes}' are used. The 'vector' column is parsed
    into a 2D tensor of shape (vector_size, grid_size*grid_size).

    Args:
        path (str): Path to the NetLogo state file.

    Returns:
        torch.Tensor: 2D tensor (vector_size, grid_size*grid_size).
    """
    with open(path, "r") as f:
        lines = f.readlines()

    start = False
    content = ""
    for line in lines:
        if start:
            if "EXTENSIONS" in line:
                break
            content += line
        if "TURTLES" in line:
            start = True

    strio = StringIO(content)
    df = pd.read_csv(strio, encoding="utf-8")

    data = df[df["breed"] == "{breed nodes}"]["vector"]
    data = (
        torch.tensor([list(map(int, re.findall(r"\d+", x))) for x in data])
        .float()
        .transpose(0, 1)
    )
    return data

def load_experiments_from_folder(glob_str_path="./data/*.state", max_batch_size=40):
    """
    Load and batch NetLogo experiments from files matching a glob pattern.

    Each experiment is grouped by vector and grid size, and batched up to max_batch_size.
    Returns a list of dicts with metadata and data tensors.

    Args:
        glob_str_path (str, optional): Glob string to find files. Defaults to "./data/*.state".
        max_batch_size (int, optional): Maximum batch size. Defaults to 40.

    Returns:
        List[dict]: Each dict contains 'vecsize', 'gridsize', 'files', and 'data'.
    """
    files = glob(glob_str_path)
    print("Files found in", glob_str_path, " - ", len(files))

    expdf = pd.DataFrame(files, columns=["file"])
    expdf["data"] = expdf["file"].apply(load_netlogo_state)
    expdf["vecsize"] = expdf["data"].apply(lambda x: x.shape[0])
    expdf["gridsize"] = expdf["data"].apply(lambda x: math.sqrt(x.shape[1]))

    # Find unique vecsize-gridsize pairs
    unique_pairs = (
        expdf.groupby(["vecsize", "gridsize"])
        .size()
        .reset_index()
        .rename(columns={0: "count"})
    )

    # Stack data to batched torch tensors if they have the same vector size and grid size
    data = []
    for idx, row in unique_pairs.iterrows():
        vecsize = row["vecsize"]
        gridsize = row["gridsize"]
        subset = expdf[(expdf["vecsize"] == vecsize) & (expdf["gridsize"] == gridsize)]
        stacked = torch.stack(subset["data"].to_list())
        # Group the data according to the max batch size
        for i in range(0, len(stacked), max_batch_size):
            data.append(
                {
                    "vecsize": vecsize,
                    "gridsize": gridsize,
                    "files": subset["file"][i : i + max_batch_size].to_list(),
                    "data": stacked[i : i + max_batch_size],
                }
            )
    return data

class BatchLogger:
    """
    Logger for batched experiments. Stores headers and data for each batch and saves/loads JSON logs.
    """
    def __init__(self, batchsize, paths):
        self.headers = [{} for _ in range(batchsize)]
        self.data = [[] for _ in range(batchsize)]
        self.paths = paths

    def log_same_header_for_all(self, **kwargs):
        """
        Update all batch headers with the same key-value pairs.
        """
        for i in range(len(self.headers)):
            self.headers[i].update(kwargs)

    def log_data_list(self, **kwargs):
        """
        Log a list of data points for each batch. Each kwarg should be a list of length batchsize.
        """
        for i in range(len(self.data)):
            current_data = {}
            for key, value in kwargs.items():
                current_data[key] = value[i]
            self.data[i].append(current_data)

    def save_all(self):
        """
        Save all batch logs to their respective files.
        """
        for i in range(len(self.headers)):
            with open(self.paths[i], "w") as f:
                json.dump({"header": self.headers[i], "data": self.data[i]}, f)

    def load_all(self):
        """
        Load all batch logs from their respective files.
        """
        for i in range(len(self.headers)):
            with open(self.paths[i], "r") as f:
                data = json.load(f)
                self.headers[i] = data["header"]
                self.data[i] = data["data"]

    def add_extra_header(self, batch_idx, **kwargs):
        """
        Add extra key-value pairs to a specific batch header.
        """
        self.headers[batch_idx].update(kwargs)

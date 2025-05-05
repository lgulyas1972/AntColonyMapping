import time
import torch
from util import BatchLogger
from tqdm import tqdm

class BatchSSGA_LS:
    """
    Implements a batched Steady-State Genetic Algorithm (SSGA) with optional local search.
    Supports both Baldwinian and Lamarckian evolution strategies.
    """

    def __init__(
        self,
        batch_size,
        grid_size,
        vector_size,
        popsize,
        mutprob,
        lsmutprob,
        lssize,
        evoltype,
        costfunc,
        device,
        grid=None,
        logpaths=None,
        save_frequency=1000,
    ):
        """
        Initialize the batched SSGA with local search.

        Args:
            batch_size (int): Number of experiments in the batch.
            grid_size (int): Size of the grid (one dimension).
            vector_size (int): Size of the feature vector.
            popsize (int): Number of individuals per batch.
            mutprob (float): Mutation probability per gene.
            lsmutprob (float): Local search mutation probability per gene.
            lssize (int): Number of local search iterations.
            evoltype (str): "Baldwinian" or "Lamarckian" evolution.
            costfunc (function): Fitness function.
            device (str): "cuda" or "cpu".
            grid (torch.Tensor, optional): Initial grid state (batch_size, vector_size, grid_size*grid_size).
            logpaths (list, optional): List of log file paths for each batch.
            save_frequency (int, optional): Save logs every N steps.
        """
        self.grid_size = grid_size
        self.vector_size = vector_size
        self.batch_size = batch_size
        self.device = device
        self.evoltype = evoltype
        self.lssize = lssize
        self.mutprob = torch.tensor(mutprob, device=device)
        self.lsmutprob = torch.tensor(lsmutprob, device=device)
        self.save_frequency = save_frequency
        self.starttime = time.time()

        assert self.evoltype in ["Baldwinian", "Lamarckian"], "Invalid evolution type"

        if grid is not None:
            self.gridvals = grid.to(device)
        else:
            self.gridvals = torch.normal(
                0, 20, (batch_size, vector_size, grid_size * grid_size)
            ).to(device)
        self.gridvals.requires_grad = False

        # Initialize population: each individual is a permutation of grid indices
        self.pop = (
            torch.vstack(
                [
                    torch.randperm(grid_size * grid_size)
                    for _ in range(batch_size * popsize)
                ]
            )
            .to(device)
            .reshape(batch_size, popsize, grid_size * grid_size)
        )
        self.popsize = popsize
        self.costfunc = costfunc

        # Compute initial fitness and apply local search
        expanded = self.gather_values(self.pop).reshape(
            batch_size * popsize, vector_size, grid_size, grid_size
        )
        self.fitness = self.costfunc(expanded).reshape(batch_size, popsize)
        self.pheno, self.fitness = self.local_search_sequential(self.pop, self.fitness)

        # Initialize logger
        if logpaths is not None:
            self.logger = BatchLogger(batchsize=batch_size, paths=logpaths)
        else:
            self.logger = BatchLogger(
                batchsize=batch_size,
                paths=[
                    "SSGA_LS_log_"
                    + time.strftime("%Y%m%d-%H%M%S")
                    + "b_"
                    + str(i)
                    + ".json"
                    for i in range(batch_size)
                ],
            )

        self.logger.log_same_header_for_all(
            finished=False,
            batch_size=batch_size,
            grid_size=grid_size,
            vector_size=vector_size,
            popsize=popsize,
            mutprob=mutprob,
            lsmutprob=lsmutprob,
            lssize=lssize,
            evoltype=evoltype,
            costfunc=costfunc.__name__,
            device=device,
            grid=grid.tolist() if grid is not None else None,
            start_time=time.strftime("%Y%m%d-%H%M%S"),
        )

    def log_config_state_paths(self, statelist):
        """
        Add state file paths to the logger headers for each batch.
        """
        for i in range(self.batch_size):
            self.logger.add_extra_header(i, statefile=statelist[i])

    def gather_values(self, popslice):
        """
        Gather grid values according to the population's permutations.

        Args:
            popslice (torch.Tensor): (batch_size, slice_size, grid_size*grid_size)

        Returns:
            torch.Tensor: (batch_size, slice_size, vector_size, grid_size, grid_size)
        """
        b, p, gg = popslice.shape
        v = self.vector_size
        return (
            torch.index_select(self.gridvals, -1, popslice.flatten())
            .reshape(b, v, b, p, gg)[0, ...]
            .transpose(0, 1)
            .transpose(1, 2)
            .reshape(b, p, v, self.grid_size, self.grid_size)
        )

    def find_missing_scalars(self, source, removed):
        """
        Find values in 'source' not present in 'removed'.

        Args:
            source (torch.Tensor): 1D tensor.
            removed (torch.Tensor): 1D tensor.

        Returns:
            torch.Tensor: 1D tensor of missing values.
        """
        return source[~torch.isin(source, removed)]

    def copy_with_crossover(self, bestidx, otheridx):
        """
        Perform one-point crossover between best and another parent.

        Args:
            bestidx (torch.Tensor): Indices of best individuals.
            otheridx (torch.Tensor): Indices of other parents.

        Returns:
            torch.Tensor: New genomes after crossover.
        """
        if self.evoltype == "Baldwinian":
            data = self.pop.clone()
        elif self.evoltype == "Lamarckian":
            data = self.pheno.clone()
        crossover = torch.randint(
            low=1,
            high=self.grid_size * self.grid_size - 1,
            size=(self.batch_size,),
            device=self.device,
        ).long()
        newgenome = torch.ones(
            (self.batch_size, self.grid_size * self.grid_size), device=self.device
        )
        for i in range(self.batch_size):
            selected = data[i, bestidx[i], : crossover[i]].clone()
            complement = self.find_missing_scalars(
                data[i, otheridx[i], :], selected
            ).clone()
            newgenome[i, :] = torch.cat((selected, complement))
        return newgenome.clone()

    def mutate_scramble(self, vector, mutprob=None):
        """
        Scramble-mutate the vector by randomly shuffling a subset of indices.

        Args:
            vector (torch.Tensor): (batch_size, grid_size*grid_size)
            mutprob (float, optional): Mutation probability per gene.

        Returns:
            torch.Tensor: Mutated vector.
        """
        batch, seq = vector.shape
        if mutprob is None:
            mutprob = self.mutprob

        shifted_perm = torch.stack(
            [torch.randperm(seq, device=self.device) + i * seq for i in range(batch)]
        )
        rands = torch.rand_like(shifted_perm, dtype=torch.float, device=self.device)
        mask = rands < mutprob
        randnums = mask.sum(axis=1)
        prevrandnums = torch.cat(
            [torch.cumsum(randnums, 0), torch.tensor([0], device=self.device)]
        )
        scramble_perm = torch.cat(
            [
                torch.randperm(randnums[i], device=self.device) + prevrandnums[i - 1]
                for i in range(batch)
            ]
        )
        selected = vector[mask].clone()
        selected = selected[scramble_perm].clone()
        vector[mask] = selected.clone()
        return vector.clone()

    def local_search_parallel(self, popslice):
        """
        Perform local search in parallel for each individual.

        Args:
            popslice (torch.Tensor): (batch_size, slice_size, grid_size*grid_size)

        Returns:
            tuple: (phenotype, fitness)
        """
        b, p, gg = popslice.shape
        expanded = (
            popslice.unsqueeze(2)
            .repeat(1, 1, self.lssize, 1)
            .reshape(b * p * self.lssize, gg)
        )
        expanded = self.mutate_scramble(expanded, self.lsmutprob).reshape(
            b, p * self.lssize, gg
        )
        expanded_batch = self.gather_values(expanded).reshape(
            b * p * self.lssize, self.vector_size, self.grid_size, self.grid_size
        )
        fits = self.costfunc(expanded_batch).reshape(b, p, self.lssize)
        minidx = torch.argmin(fits, axis=2)
        minfits = torch.gather(fits, 2, minidx.unsqueeze(2)).squeeze(2)
        minpops = torch.gather(
            expanded.reshape(b, p, self.lssize, gg),
            2,
            minidx.unsqueeze(2).unsqueeze(2).repeat(1, 1, self.lssize, gg),
        )[:, :, 0, :]
        return minpops, minfits

    def local_search_sequential(self, popslice, fitslice):
        """
        Perform sequential local search for each individual.

        Args:
            popslice (torch.Tensor): (batch_size, slice_size, grid_size*grid_size)
            fitslice (torch.Tensor): (batch_size, slice_size)

        Returns:
            tuple: (phenotype, fitness)
        """
        pop = popslice.clone()
        fit = fitslice.clone()
        b, p, gg = pop.shape

        for _ in range(self.lssize):
            mutated = self.mutate_scramble(
                pop.reshape(b * p, gg), self.lsmutprob
            ).reshape(b, p, gg)
            expanded_batch = self.gather_values(mutated).reshape(
                b * p, self.vector_size, self.grid_size, self.grid_size
            )
            fits = self.costfunc(expanded_batch).reshape(b, p)
            mask = fit < fits
            pop[mask, :] = mutated[mask, :].clone()
            fit[mask] = fits[mask].clone()
        return pop, fit

    def step(self):
        """
        Perform a single SSGA step: selection, crossover, mutation, local search, and replacement.
        """
        worstidx = torch.argmax(self.fitness, axis=1)
        bestidx = torch.argmin(self.fitness, axis=1)
        otheridx = torch.randint(
            0, self.popsize, (self.batch_size,), device=self.device
        )
        while (otheridx == bestidx).sum() > 0 or (otheridx == worstidx).sum() > 0:
            otheridx = torch.randint(
                0, self.popsize, (self.batch_size,), device=self.device
            )
        newgenome = self.copy_with_crossover(bestidx, otheridx).long()
        newgenome = self.mutate_scramble(newgenome).long()
        expanded_batch = self.gather_values(
            newgenome.reshape(self.batch_size, 1, -1)
        ).squeeze(1)
        newfitness = self.costfunc(expanded_batch)
        newpheno, newfitness = self.local_search_sequential(
            newgenome.reshape(self.batch_size, 1, -1),
            newfitness.reshape(self.batch_size, 1),
        )
        newpheno = newpheno.squeeze(1)
        newfitness = newfitness.squeeze(1)
        self.pop[torch.arange(self.batch_size, device=self.device), worstidx, :] = (
            newgenome.clone()
        )
        self.fitness[torch.arange(self.batch_size, device=self.device), worstidx] = (
            newfitness.clone()
        )
        self.pheno[torch.arange(self.batch_size, device=self.device), worstidx, :] = (
            newpheno.clone()
        )

    def optimize(self, steps):
        """
        Run the SSGA optimization for a given number of steps.

        Args:
            steps (int): Number of optimization steps.
        """
        self.logger.log_data_list(
            min_fitness=self.fitness.min(dim=1).values.tolist(),
            step=[-1] * self.batch_size,
        )
        self.logger.save_all()

        for i in tqdm(range(steps)):
            self.step()
            self.logger.log_data_list(
                min_fitness=self.fitness.min(dim=1).values.tolist(),
                step=[i] * self.batch_size,
            )
            if i % self.save_frequency == 0:
                self.logger.save_all()

        self.logger.log_same_header_for_all(
            runtime_s=time.time() - self.starttime, finished=True
        )
        self.logger.save_all()

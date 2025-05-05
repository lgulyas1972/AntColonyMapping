from ssma.SSMA_batched_local import BatchSSMA
from grid import sum_avg_dist
from util import load_experiments_from_folder
import torch
import os
import gc


# -------------------- Experiment Config --------------------
# Folder containing experiment state files from NetLogo runs

experiment_state_folder = "../data/acom"

# -------------------- Baldwinian Evolution Config --------------------
max_batch_size_baldwin = 10  # Maximum batch size for experiments (depends on available memory)
pop_size_baldwin = 100       # Population size per experiment
mut_rate_baldwin = 0.005     # Mutation rate per gene
device_baldwin = "cuda"      # "cuda" for GPU, "cpu" for CPU
generations_baldwin = 20000  # Number of generations to run
evol_type_baldwin = "Baldwinian"  # Evolution type: "Baldwinian" or "Lamarckian"
local_random_rate_baldwin = 0.005 # Local search mutation rate per gene
local_search_iters_baldwin = 2    # Number of local search iterations
logfolder_name_baldwin = "./logs/log_baldwin_ssma"
save_frequency_baldwin = 10000    # Save logs every N generations

# -------------------- Lamarckian Evolution Config --------------------
max_batch_size_lamarck = 10
pop_size_lamarck = 100
mut_rate_lamarck = 0.005
device_lamarck = "cuda"
generations_lamarck = 20000
evol_type_lamarck = "Lamarckian"
local_random_rate_lamarck = 0.005
local_search_iters_lamarck = 2
logfolder_name_lamarck = "./logs/log_lamarck_ssma"
save_frequency_lamarck = 10000

# Set torch random seed for reproducibility
torch.manual_seed(42)

# --- Baldwinian Evolution Run ---
if not os.path.exists(logfolder_name_baldwin):
    os.makedirs(logfolder_name_baldwin)

# Load and batch experiments from data files
experiments = load_experiments_from_folder(
    f"{experiment_state_folder}/*.state", max_batch_size=max_batch_size_baldwin
)

for i, exp in enumerate(experiments):
    grid_size = int(exp["gridsize"])
    vector_size = int(exp["vecsize"])
    grid = exp["data"]
    files = exp["files"]
    print(files)
    actual_batch_size = grid.shape[0]
    print(
        "Running experiment",
        i + 1,
        "/",
        len(experiments),
        "gsize:",
        grid_size,
        "vsize:",
        vector_size,
        "batch:",
        actual_batch_size,
    )

    optimizer = BatchSSMA(
        actual_batch_size,
        grid_size,
        vector_size,
        pop_size_baldwin,
        mut_rate_baldwin,
        local_random_rate_baldwin,
        local_search_iters_baldwin,
        evol_type_baldwin,
        sum_avg_dist,
        device_baldwin,
        grid=grid,
        logpaths=[f.replace(experiment_state_folder, logfolder_name_baldwin) for f in files],
        save_frequency=save_frequency_baldwin,
    )

    optimizer.optimize(generations_baldwin)

    del optimizer
    gc.collect()

# --- Lamarckian Evolution Run ---
if not os.path.exists(logfolder_name_lamarck):
    os.makedirs(logfolder_name_lamarck)

experiments = load_experiments_from_folder(
    f"{experiment_state_folder}/*.state", max_batch_size=max_batch_size_lamarck
)

for i, exp in enumerate(experiments):
    grid_size = int(exp["gridsize"])
    vector_size = int(exp["vecsize"])
    grid = exp["data"]
    files = exp["files"]
    print(files)
    actual_batch_size = grid.shape[0]
    print(
        "Running experiment",
        i + 1,
        "/",
        len(experiments),
        "gsize:",
        grid_size,
        "vsize:",
        vector_size,
        "batch:",
        actual_batch_size,
    )

    optimizer = BatchSSMA(
        actual_batch_size,
        grid_size,
        vector_size,
        pop_size_lamarck,
        mut_rate_lamarck,
        local_random_rate_lamarck,
        local_search_iters_lamarck,
        evol_type_lamarck,
        sum_avg_dist,
        device_lamarck,
        grid=grid,
        logpaths=[f.replace(experiment_state_folder, logfolder_name_lamarck) for f in files],
        save_frequency=save_frequency_lamarck,
    )

    optimizer.optimize(generations_lamarck)

    del optimizer
    gc.collect()

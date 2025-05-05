# Steady-state Memetic Algorithm (SSMA) Pytorch Implementation

## SSMA for Mapping Problems

This implementation provides a steady-state memetic algorithm using PyTorch for optimization tasks. It leverages GPU acceleration for improved performance. It parallelizes fitness evaluation using fold and unfold to evaluate multiple solutions simultaneously, enhancing the overall efficiency of the algorithm. The algorithm uses parallelized local search and a scramble mutation operator that is also parallelized. Both Baldwinian and Lamarckian approaches are implemented.

The code is designed to copy existing NetLogo runs with the same parameters and starting state. It needs the NetLogo state files to load experimental parameters from them to be placed in `./data/acom/*.state`.

## Installation

The code is implemented in Python 3.10 and PyTorch 2.1.0 (with CUDA 11.8). The code was tested on Windows 11, but no OS-specific restrictions are expected.

Install required packages using pip:

```bash
pip install -r requirements.txt
```

## Usage

Configure the parameters of the algorithm in `main_local.py`'s header.
The default parameters are set to run experiments needed for our paper using a 24GB VRAM GPU. Adjust batch size if you have less VRAM available.

To run the algorithm, use the following command:

```bash
python main_local.py
```

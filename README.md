# AntColonyMapping

This is the source code for the Ant Colony Mapping (ACOM) algorithm, mapping feature vectors to discrete topologies (networks), such that similar vectors occupy neighbouring nodes. The problem is inspired by Kohonen’s Self-Organising Maps (SOM).

## Requirements

The ACOM code is implemented in the NetLogo (6.2.0) framework available at [NetLogo](https://ccl.northwestern.edu/netlogo/).

The SSGA baseline's GPU implementation is based on Python 3.10 and PyTorch 2.1.0 (with CUDA 11.8). The code was tested on Windows 11, but no OS-specific restrictions are expected.

## Data

The NetLogo generated data files are available in the `data/acom` directory. The baseline SSMA data files can be generated using the source code provided in the `ssga` directory.

## Publication

The paper detailing the algorithm is under review. Until it is published, please refer to this repository using the following citation:

```bibtex
@misc{gulyas2025antcolony,
  author       = {L\'{a}szl\'{o} Guly\'{a}s and Natabara M\'{a}t\'{e} Gy\"{o}ngy\"{o}ssy and J\'{a}nos Botzheim},
  title        = {Ant Colony Mapping for Dimensionality Reduction to Discrete Spaces: Code Repository},
  year         = {2025},
  howpublished = {\url{https://github.com/lgulyas1972/AntColonyMapping}},
  note         = {Accessed: ADD DATE HERE},
  institution  = {ELTE E\"{o}tv\"{o}s Lor\'{a}nd University, Department of Artificial Intelligence},
}
```

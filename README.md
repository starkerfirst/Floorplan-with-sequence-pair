# Floorplan with sequence pair
## Overview

Chip floorplanning is a critical step in the physical synthesis process that involves determining the positions and shape of circuit modules to minimize chip area while satisfying various constraints. I present my floorplan implementation of **Simulated Annealing algorithm using sequence pair representation**. The sequence pair representation encodes both slicing floorplan and non-slicing floorplan into sequence pair space, while simulated annealing helps navigate the complex solution space to find near-optimal floorplans.

My implementation uses Python 3.10.16 and only imports minimal libraries like pyplot and math. The seed for random generator is 754. The iteration limit is set to 15000 due to the time limit. I tested the algorithm on four benchmark problems of increasing size: 10, 50, 200, and 500 modules. The final results are listed as follows (all chip widths meet constraints):

| #Modules | Chip Width | Chip Height | Utilization | Runtime (s) |
|----------|------------|-------------|-------------|-------------|
| 10       | 49         | 16          | 94.52%      | 1.31        |
| 50       | 97         | 56          | 73.90%      | 8.59        |
| 200      | 298        | 353         | 51.10%      | 88.44       |
| 500      | 998        | 1447        | 42.69%      | 566.95      |

## Usage

After changing to `floorplan` directory, simply run the command `python floorplan.py` to execute the program. The program will generate output files and visualization figs in `./output` and `./fig`. The input traces have already been provided in the directory `./2-floorplan-v0`.
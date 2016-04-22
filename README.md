# PredefinedKMC
Adam's predefined transition KMC code

Requires LAKMC graphs module

This script runs KMC on a hexagonal ZnO(0001) surface (to be generalised later)
Assumes PBC in x and z directions and deposition in y
Grid positions start at (0,0,0)

Setup:
KMC.py          - main source code (python)
lattice.dat     - initial lattice file
Volumes         - directory containing saved volume lattices
Transitions     - directory containing transitions for each volume
Output          - directory containing lattices after KMC steps


# LatAKMC
##Adam's Lattice based Adaptive KMC code

### Requirements
LAKMC Graphs, NEB, Lattice, Minimise, Input, Vectors modules <br>
lkmcInput.IN    - input file from LAKMC <br>
md directory    - directory full of information for MD <br>

### About to code
This script runs KMC on a hexagonal <b>ZnO(0001)</b> surface (to be generalised later) <br>
Assumes PBC in <i>x</i> and <i>z</i> directions and deposition in <i>y</i> <br>
Grid positions start at (0,0,0) <br> 


### Setup:
KMC.py            - main source code (python)<br>
lattice.dat      - initial lattice file<br>
Volumes           - directory containing saved volume lattices<br>
Transitions       - directory containing transitions for each volume<br>
Output            - directory containing lattices after KMC steps<br>

#### Volumes
Each volume file is named by the volume hashkey. <br>
Format:<br>
Number of atoms in lattice<br>
List of atoms in lattice (species, x, y, z coordinates)<br>
List of Atom indices in hashkey<br>

#### Transitions
Each Transition file is named by the initial volume hashkey. <br>
Format:<br>
Displacement vectors in integer lattice units (x y z) <br>
List of final Hashkeys and Barrier heights (eV)

#### Output
Each ouput file is called KMC + KMC step.<br>
Format:<br>
Number of atoms in lattice<br>
Simulation time<br>
Barrier height (or equivalent for deposition)<br>
Simulation cell dimensions<br>
List of atoms in lattice (species, x, y, z coordinates, charge)<br>

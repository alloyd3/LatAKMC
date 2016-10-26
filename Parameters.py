# -*- coding: utf-8 -*-
"""
Input module.

"""

import os
import sys
from LKMC import Utilities

class Input(object):
    def __init__(self):
        # default input
        self.jobStatus = 'CNTIN'            # BEGIN or CNTIN run
        self.atom_species = 'Ag'         # species to deposit
        self.numberDepos = 35  	        # number of initial depositions
        self.total_steps = 100000           # total number of steps to run
        self.latticeOutEvery = 10         # write output lattice every n steps
        self.volumesOutEvery = 20        # write out volumes to file every n steps
        self.temperature = 300           # system temperature in Kelvin
        self.prefactor = 1.00E+13        # fixed prefactor for Arrhenius eq. (typically 1E+12 or 1E+13)
        self.boltzmann = 8.62E-05        # Boltzmann constant (8.62E-05)
        self.graphRad = 5.9                # graph radius of defect volumes (Angstroms)
        self.depoRate = 5184              # deposition rate
        self.maxMoveCriteria = 0.87        # maximum distance an atom can move after relaxation (pre NEB)
        self.maxHeight = 30              # Dimension of cell in y direction (A)
        self.includeUpTrans = 0          # Booleon: Include transitions up step edges (turning off speeds up simulation)
        self.includeDownTrans = 1        # Booleon: Include transitions down step edges
        self.statsOut = 0              # Recieve extra information from your run
        self.useBasin = 1           # Booleon: use the basin method or not
        self.basinBarrierTol = 0.25      # barriers below this are considered in a basin (eV)
        self.basinBarrierSubTol = 0.40   # if one barrier is above this, it is considered an escaping transition not internal
        self.basinDistTol = 0.6         # distance between states to be considered the same state (A)
        self.checkMoveDist = 2          # distance used in checkMoveDist. Do not allow an atom to move within this distance of another atom. Needed for basin method
        self.reverseBarrierTol = 0.03   # Tolerance to allow transitions with reverse barriers greater than this only
        self.maxCoordNum = 9            # max coordination to be considered a Defects
        self.bondDist = 3.4             # bond distance between atoms

        # for (0001) ZnO only
        self.x_grid_dist = 0.9497411251   # distance in x direction between each atom in lattice (A)
        self.y_grid_dist = 1.55            # distance in y direction between surface and first layer (eg O layer - Ag layer)
        self.y_grid_dist2 = 2.1             # distance between deposited layers (eg. Ag layer - Ag layer)
        self.z_grid_dist = 1.6449999809   # distance in z direction between each atom in lattice

def createInputObject():
    return Input();

def getInput(inputFile = "input.IN"):

    try:
        f = open(inputFile, "r")
    except:
        sys.exit(__name__+": ERROR: could not open file: " + inputFile)

    inputObject = createInputObject()

    lookForName = 1
    lookForValue = 0

    for line in f:
        paramValueStr = None
        if line[0] == '!':
            continue
        if (lookForName == 1):
            if line[0] != "%":
                print line
                sys.exit(__name__+": ERROR: something is wrong with the input file: " + inputFile)
            else:
                lineLen = len(line)
                paramName = line[1:lineLen-1]

                if hasattr(inputObject , paramName):
                    paramType = type(getattr(inputObject , paramName)).__name__
                else:
                    sys.exit(__name__+": ERROR: unexpected parameter:" + paramName + ": " + inputFile)

                lookForName = 0
                lookForValue = 1

        elif (lookForValue == 1):
            lineLen = len(line)
            paramValueStr = line[0:lineLen-1]

            paramValue = Utilities.convertStrToType(paramValueStr, paramType)

            setattr(inputObject , paramName, paramValue)

            lookForName = 1
            lookForValue = 0
        else:
            sys.exit(__name__+": ERROR: undefined action while reading the input file: " + inputFile)

    f.close()

    return inputObject

# def main():
#     params = getInput()
#     print params.depoRate
#
# if __name__ == "__main__":
#     main()

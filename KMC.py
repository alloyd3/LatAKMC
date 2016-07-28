# -*- coding: utf-8 -*-

# This script runs KMC on a hexagonal ZnO(0001) surface
# TODO: make more general for any hexagonal surface
# Assumes PBC in x and z directions and deposition in y
# grid positions start at (0,0,0)


# copyright Adam Lloyd 2016
import time
import sys
import os
import shutil
import random
import math
import copy
from decimal import Decimal
import numpy as np
from LKMC import Graphs, NEB, Lattice, Minimise, Input, Vectors

#------------------------------------------------------------------------------
#- User inputs hard coded in script
jobStatus = 'BEGIN'            # BEGIN or CNTIN run
atom_species = 'Ag'         # species to deposit
numberDepos = 2  	        # number of initial depositions
total_steps = 3           # total number of steps to run
latticeOutEvery = 1         # write output lattice every n steps
volumesOutEvery = 10        # write out volumes to file every n steps
temperature = 300           # system temperature in Kelvin
prefactor = 1.00E+13        # fixed prefactor for Arrhenius eq. (typically 1E+12 or 1E+13)
boltzmann = 8.62E-05        # Boltzmann constant (8.62E-05)
graphRad = 5.9                # graph radius of defect volumes (Angstroms)
depoRate = 1            # deposition rate
maxMoveCriteria = 0.6        # maximum distance an atom can move after relaxation (pre NEB)
MaxHeight = 30              # Dimension of cell in y direction (A)
IncludeUpTrans = 0          # Booleon: Include transitions up step edges (turning off speeds up simulation)
IncludeDownTrans = 1        # Booleon: Include transitions down step edges
StatsOut = False               # Recieve extra information from your run
useBasin = True            # Booleon: use the basin method or not
basinBarrierTol = 0.3      # barriers below this are considered in a basin (eV)
basinDistTol = 0.3         # distance between states to be considered the same state (A)


# for (0001) ZnO only
x_grid_dist = 0.9497411251   # distance in x direction between each atom in lattice (A)
y_grid_dist = 1.55            # distance in y direction between surface and first layer (eg O layer - Ag layer)
y_grid_dist2 = 2.1             # distance between deposited layers (eg. Ag layer - Ag layer)
z_grid_dist = 1.6449999809   # distance in z direction between each atom in lattice
#------------------------------------------------------------------------------

# Defined some useful functions

# classes used in hashkey calculation from LKMC
class lattice(object):
    def __init__(self):
        self.pos = []
        self.specie = []
        self.cellDims = [box_x,0,0,0,box_y,0,0,0,box_z]
        self.specieList = ['O_','Zn','Ag']
        self.NAtoms = 0
        self.charge = 0

class params(object):
    def __init__(self):
        self.graphRadius = 5.9

class volume(object):
    def __init__(self):
        self.hashkey = None
        self.directions = []
        self.finalKeys = {}
        self.pos = []
        self.specie = []
        self.volumeAtoms = []

    def addTrans(self, direction, finalKey, barrier, rate, reverseBarrier):
        if direction not in self.directions:
            self.directions.append(direction)

        if finalKey not in self.finalKeys:
            newKey = key()
            newKey.barrier = barrier
            newKey.rate = rate
            newKey.reverseBarrier = reverseBarrier
            # newKey.hashkey = finalKey
            self.finalKeys[finalKey] = newKey

    def addDirection(self, direction):
        if direction not in self.directions:
            self.directions.append(direction)

    # store new volume atoms in memory
    def addVolumeAtoms(self, volumeAtoms,lattice_positions,specie_list):
        if len(self.volumeAtoms) < 1:
            for i in xrange(len(lattice_positions)):
                self.pos.append(lattice_positions[i])
            for i in xrange(len(specie_list)):
                self.specie.append(specie_list[i])
            for i in xrange(len(volumeAtoms)):
                self.volumeAtoms.append(volumeAtoms)

class key(object):
    def __init__(self):
        self.barrier = None
        self.rate = None
        self.reverseBarrier = None
        # self.hashkey = None

# transition object for the basin
class basinTransition(object):
    def __init__(self,finPos,rate,barrier,reverseBarrier):
        self.finPos = finPos
        self.rate = rate
        self.barrier = barrier
        self.reverseBarrier = reverseBarrier
        self.finRef = None

# basinPosition object for basinPos list (similar to basin DV in LAKMC)
class basinPosition(object):
    def __init__(self):
        self.transitionList = []
        self.explored = 0
        self.iniPos = None

# basin type object
class basin(object):
    def __init__(self):
        self.atomNum = None
        self.currentPos = None
        self.positions = []
        self.basinPos = []
        self.exploredList = []
        self.transitionList = []
        self.connectivity = None

    # add transition to basin
    def addTransition(self,iniPos,finPos,rate,barrier,reverseBarrier):
        # is this an internal event or escaping?
        if barrier < basinBarrierTol or reverseBarrier < basinBarrierTol:
            flag = 0
        else:
            flag = 1

        createFlagF = 1
        i=0
        # check if initial position exists in the basin
        for i in range(len(self.basinPos)):
            pos = self.basinPos[i].iniPos
            if PBC_distance(iniPos[0],iniPos[1],iniPos[2],pos[0],pos[1],pos[2]) < basinDistTol:
                createFlagF = 0
                break
        # if does not exist, add
        if createFlagF:
            i = len(self.basinPos)
            newTrans = basinTransition(finPos,rate,barrier,reverseBarrier)
            newPos = basinPosition()
            newPos.explored = 1
            newPos.iniPos = iniPos
            newPos.transitionList.append(newTrans)
            self.basinPos.append(newPos)
        else:
            # check if this transition already exists in the basin
            createFlagF = 1
            for trans in self.basinPos[i].transitionList:
                if PBC_distance(finPos[0],finPos[1],finPos[2],trans.finPos[0],trans.finPos[1],trans.finPos[2]) < basinDistTol:
                    createFlagF = 0
            # if not, add transition
            if createFlagF:
                newTrans = basinTransition(finPos,rate,barrier,reverseBarrier)
                self.basinPos[i].transitionList.append(newTrans)
                self.basinPos[i].explored = 1


        # if non escaping transition, add reverse transition
        if not flag:
            # locate final position in basin
            j=0
            createFlag = 1
            # check final position exists in the basin
            for j in range(len(self.basinPos)):
                pos = self.basinPos[j].iniPos
                if PBC_distance(finPos[0],finPos[1],finPos[2],pos[0],pos[1],pos[2]) < basinDistTol:
                    createFlag = 0
                    break

            # if does not exist, add
            if createFlag:
                j = len(self.basinPos)
                revRate = calc_rate(reverseBarrier)
                newTransR = basinTransition(iniPos,revRate,reverseBarrier,barrier)
                newTransR.finRef = i
                newPosR = basinPosition()
                newPosR.explored = 0
                newPosR.iniPos = finPos
                newPosR.transitionList.append(newTransR)
                self.basinPos.append(newPosR)
            else:
                # check transition exists in the basin
                createFlag = 1
                for trans in self.basinPos[j].transitionList:
                    if PBC_distance(iniPos[0],iniPos[1],iniPos[2],trans.finPos[0],trans.finPos[1],trans.finPos[2]) < basinDistTol:
                        createFlag = 0
                # add transition if not
                if createFlag:
                    newTransR = basinTransition(finPos,rate,barrier,reverseBarrier)
                    newTransR.finRef = i
                    self.basinPos[j].transitionList.append(newTransR)

            # assign ref to positions
            if createFlagF:
                self.basinPos[i].transitionList[-1].finRef = j
                # print "Adding transition: ",i,j

    # build connectivity matrix. All elements are transition numbers
    def buildConnectivity(self):
        N = len(self.basinPos)
        self.connectivity = [[[] for i in range(N)] for j in range(N)]

        # create connectivity matrix
        for i in range(len(self.basinPos)):
            pos = self.basinPos[i]
            for j in range(len(pos.transitionList)):
                trans = pos.transitionList[j]
                if trans.finRef is not None:
                    self.connectivity[i][trans.finRef].append(j)

        print "Connectivity matrix:"
        for i in range(N):
            print self.connectivity[i]

        explor = []
        for i in range(N):
            explor.append(self.basinPos[i].explored)
        print "Explored states: ", explor


    # check if position is in this basin
    def thisBasin(self, pos):
        cPos = self.currentPos
        if PBC_distance(cPos[0],cPos[1],cPos[2],pos[0],pos[1],pos[2]) < basinDistTol:
            return True

        for basPos in self.basinPos:
            bPos = basPos.iniPos
            if PBC_distance(bPos[0],bPos[1],bPos[2],pos[0],pos[1],pos[2]) < basinDistTol:
                return True

        return False

    # calculate mean rates within the basin
    def meanRate(self):
        result = []
        stateNum = 0
        stateMapping = []
        for i in range(len(self.basinPos)):
            if self.basinPos[i].explored:
                stateNum +=1
                stateMapping.append(i)

        if stateNum == 1:
            return result

        transMatrix = np.zeros([stateNum, stateNum])
        tao1 = np.zeros(stateNum, np.float64)
        transNum = np.zeros(stateNum,dtype=np.int)
        for j in range(stateNum):
            transNum[j] = len(self.basinPos[stateMapping[j]].transitionList)
            sum = 0.0
            for i in range(transNum[j]):
                sum += self.basinPos[stateMapping[j]].transitionList[i].rate
            if sum == 0.0:
                print "Error: Sum of rates for State %d (%d trans) is zero!"%(j, transNum[j])
                return result
            tao1[j] = 1 / sum

        for i in range(stateNum):
            for j in range(stateNum):
                for k in self.connectivity[stateMapping[i]][stateMapping[j]]:
                    transMatrix[j][i] += tao1[i] * self.basinPos[stateMapping[i]].transitionList[k].rate



        for k in range(len(self.basinPos)):
            q = self.basinPos[k]
            if PBC_distance(q.iniPos[0],q.iniPos[1],q.iniPos[2],self.currentPos[0],self.currentPos[1],self.currentPos[2]) < basinDistTol:
                startDV = k
                break

        occupVect0 = np.zeros(stateNum)
        occupVect0[startDV] = 1.0

        matrix2bInv = np.matrix(np.identity(stateNum)) - transMatrix
        try:
            occupVect = np.linalg.inv(matrix2bInv)
        except np.linalg.LinAlgError:
            print "Error: Transition Matrix not invertible."
            return None
        occupVect = np.inner(occupVect, occupVect0)
        occupVect = np.squeeze(np.asarray(occupVect))

        tao = np.zeros(stateNum, np.float64)
        taoSum = 0.0
        for j in range(stateNum):
            tao[j] = tao1[j] * occupVect[j]
            taoSum += tao[j]


        negRate = False
        for i in range(stateNum):
            for j in range(transNum[i]):
                finalDV = self.basinPos[stateMapping[i]].transitionList[j].finRef
                if finalDV is not None:
                    if self.basinPos[finalDV].explored:
                        result.append(0.0)
                    else:
                        result.append(self.basinPos[stateMapping[i]].transitionList[j].rate)
                else:
                    localRate = tao[i]/taoSum*self.basinPos[stateMapping[i]].transitionList[j].rate
                    result.append(localRate)
                    if localRate < 0.0:
                        print "WARNING: get a negative mean rate!! %r"%localRate
                        negRate = True

        print "Mean rate results: ", result

        if not negRate:
            return result
        else:
            return None

    # if basin is being deleted, create normal event list
    def addUnchangedEvents(self,atomNum):
        event_list = []
        for trans in self.basinPos[0].transitionList:
            event = [trans.rate,atomNum,trans.finPos,trans.barrier]
            event_list.append(event)

        return event_list

    # update rates within the basin and create new events
    def addChangedEvents(self,atomNum):
        event_list = []
        result = self.meanRate()

        if len(result):
            k=0
            for i in range(len(self.basinPos)):
                if self.basinPos[i].explored:
                    for j in range(len(self.basinPos[i].transitionList)):
                        trans = self.basinPos[i].transitionList[j]
                        newRate = result[k]
                        event = [newRate,atomNum,trans.finPos,trans.barrier]
                        event_list.append(event)
                        k += 1
                else:
                    for j in range(len(self.basinPos[i].transitionList)):
                        trans = self.basinPos[i].transitionList[j]
                        event = [trans.rate,atomNum,trans.finPos,trans.barrier]
                        event_list.append(event)
        else:
            for i in range(len(self.basinPos)):
                if self.basinPos[i].explored:
                    for j in range(len(self.basinPos[i].transitionList)):
                        trans = self.basinPos[i].transitionList[j]
                        event = [trans.rate,atomNum,trans.finPos,trans.barrier]
                        event_list.append(event)

        return event_list

# calculate the rate of an event given barrier height (Arrhenius eq.)
def calc_rate(barrier):
    rate = prefactor * math.exp(- barrier / (boltzmann * temperature))
    return rate

# find barrier height given rate
def find_barrier_height(rate):
    barrier = round(- math.log(rate / prefactor) * (boltzmann * temperature),6)
    return barrier

# read lattice header at input filename
def read_lattice_header(input_lattice_path):
    file = open(input_lattice_path, 'r')
    latline = file.readline()
    line = latline.split()
    atoms = int(line[0])
    latline = file.readline()
    line = latline.split()
    box_x = float(line[0])
    box_y = float(line[1])
    box_z = float(line[2])
    file.close()
    return atoms,box_x,box_y,box_z

# read in lattice file
def read_lattice(lattice,m):
    if (os.path.isfile(lattice)):
        input_file = open(lattice, 'r')
        for i in range(0,m):
            line = input_file.readline()
        latticeLines = []
        while 1:
            latline = input_file.readline()
            line = latline.split()
            if len(line) < 5:
                break
            latticeLines.append([str(line[0]),float(line[1]),float(line[2]),float(line[3]),float(line[4])])
            #latticeLines.append(line)
        return latticeLines
    else:
        print "Cannot read lattice"
        print input_lattice_path
        sys.exit()
    input_file.close()

# function reads of lattice, and finds the maximum y-coordinate
def find_max_height(input_lattice_path):
    max_height = 0.0
    if (os.path.isfile(input_lattice_path)):
        input_file = open(input_lattice_path, 'r')
        for i in range(0,2):
            line = input_file.readline()
            #print i, line
        while 1:
            line = input_file.readline()
            if not line: break
            #print line
            line = line.split()
            if(float(line[2]) > max_height):
                max_height = float(line[2])
        input_file.close()
        return max_height
    else:
        print " Input file not found in find_max_height(input_filename) function"
        print input_lattice_path
        sys.exit()

# find max height and species of max atom at a point in the x,z plane
def find_max_height_at_points(surface_lattice, full_depo_index, x, z):
    max_height = 0.0
    max_height_atom = None
    full_lattice = surface_lattice + full_depo_index

    # check each x value in lattice against given x
    for i in xrange(len(full_lattice)):
        line = full_lattice[i]
        new_x = round(PBC_pos(float(line[1])-x,box_x),2)
        if new_x == 0 or new_x == round(box_x,2):
            new_z = round(PBC_pos(float(line[3])-z,box_z),2)

            # if x value is the same, check z value
            if new_z == 0 or new_z == round(box_z,2):
                if(float(line[2]) > max_height):
                    max_height = float(line[2])
                    max_height_atom = str(line[0])

    del full_lattice
    return max_height, max_height_atom

# function reads off lattice, and finds atom below a point
def find_atom_below(surface_lattice, full_depo_index,x,z):
    max_height = 0.0
    max_height_atom = None
    atom_below = None

    # check each x value in lattice against given x
    full_lattice = surface_lattice + full_depo_index
    for i in xrange(len(full_lattice)):
        line = full_lattice[i]
        new_x = round(PBC_pos(float(line[1])-x,box_x),2)
        if new_x == 0 or new_x == round(box_x,2):
            new_z = round(PBC_pos(float(line[3])-z,box_z),2)

            # if x value is the same, check z value
            if new_z == 0 or new_z == round(box_z,2):
                atom_below = max_height_atom
                max_height = float(line[2])
                max_height_atom = str(line[0])

    del full_lattice
    return max_height, atom_below

# find how many grid points in each direction
def grid_size(box_x,box_y,box_z):
    #assuming first atom starts at 0,0,0
    x_grid_points = round(box_x/x_grid_dist)
    y_grid_points = round(box_y/y_grid_dist)
    z_grid_points = round(box_z/z_grid_dist)

    box_x = x_grid_points * x_grid_dist
    box_z = z_grid_points * z_grid_dist

    return x_grid_points,y_grid_points,z_grid_points, box_x, box_z

# find random position in x,z plane to deposit
def deposition_xz(box_x,box_z,x_grid_dist,z_grid_dist):
	x_random = random.random()
	z_random = random.random()

    # check if row is even or odd
	x_coordinate = round((x_random * box_x)/x_grid_dist)
	even = x_coordinate % 2

	# case where deposit on far edge
	if x_coordinate == x_grid_points:
		even = 0

	x_2 = x_coordinate * x_grid_dist
	x_coordinate = PBC_pos(x_coordinate * x_grid_dist, box_x)


	z_coordinate = round((z_random * box_z * 0.5)/z_grid_dist)
	z_coordinate = PBC_pos(z_coordinate * z_grid_dist * 2 + even * z_grid_dist, box_z)

	return x_coordinate ,z_coordinate

# returns height of and species of neighbour atoms
def deposition_y(full_depo_index,x_coord,z_coord):
    # max height at single point
    #y_max_0, atom_below = find_max_height_at_point(lattice_path,x_coord,z_coord)
    y_max_0, atom_below = find_max_height_at_points(surface_lattice, full_depo_index,x_coord, z_coord)

    # check height at surrounding points
    neighbour_pos, neighbour_species = find_neighbours(x_coord,z_coord,atom_below,y_max_0,full_depo_index)
    neighbour_heights = []

    i = 0
    while i < 7:
    	neighbour_heights.append(round(neighbour_pos[i*3+1],6))
    	i += 1
#     print neighbour_heights


    # find max height of local points
    y_max = max(neighbour_heights)

    # add y dist to max height found
    y_coordinate = round(y_max,6)
    return y_coordinate, neighbour_species, neighbour_heights

# do deposition
def deposition(box_x,box_z,x_grid_dist,z_grid_dist,full_depo_index,natoms):
    x_coord, z_coord = deposition_xz(box_x,box_z,x_grid_dist,z_grid_dist)
    y_coord, nlist, hlist = deposition_y(full_depo_index,x_coord,z_coord)

    # Potential Deposition erros for ZnO-Ag system
    maxAg = []
    nAg = nlist.count('Ag')
    nH = hlist.count(y_coord)
    if nlist[0] != 'Ag':
        for i in xrange(len(nlist)):
            if nlist[i] == 'Ag' and hlist[i] == y_coord:
                maxAg.append(i)
                if len(maxAg) == 3:
                    print "landed on top of 3 Ag atoms"

    # check if y_coord is at least surface height
    if (y_coord < initial_surface_height):
        print "ERROR: y height is less than initial surface height", y_coord
        sys.exit()
        return

    # if y_coord is initial surface height
    if round(y_coord - initial_surface_height,2) == 0:
        y_coord += y_grid_dist

    # if y_coord > initial surface height
    else:
        y_coord += y_grid_dist2


    # print "Trying to deposit %s atom at %f, %f, %f" % (atom_species, x_coord, y_coord, z_coord)
    #print nlist
    if len(maxAg) != 3:
        for x in nlist:
            if x == 'Ag':
                # print "deposited near Ag, redo deposition"
                return None
        if nlist[0] == 'O_':
        	# print "deposited on top of surface Oxygen: unstable position"
        	return None
        # if nlist[0] == 'Zn':
        # 	print "deposited on top of surface Oxygen: unstable position"
        # 	return None


    natoms += 1

    print "SUCCESS: Number of atoms: ", natoms

    deposition_list = [atom_species, x_coord, y_coord, z_coord, natoms]
    return deposition_list

# find PBC position in one direction
def PBC_pos(x,box_x):
	x = round(x,6)
	box_x = round(box_x,6)
	if x >= box_x:
		while x >= box_x:
			x = x - box_x
	if x < 0:
		while x < 0:
			x = x + box_x
	return x

def PBC_distance(x1,y1,z1,x2,y2,z2):
    rx = x2-x1
    ry = y2-y1
    rz = z2-z1

    # assume PBC in x and z
    rx = rx - round( rx / box_x) * box_x
    #ry = ry - round( ry / box_y) * box_y
    rz = rz - round( rz / box_z) * box_z

    # square of distance
    r2 = rx * rx + ry * ry + rz * rz
    return r2

# find x and z of the 6 positions surrounding points (1-6)
def find_neighbours(x,z,atom_below,y_max_0,full_depo_index):
    n_x = PBC_pos(x + 2 * x_grid_dist,box_x)
    n_z = z

    nw_x = PBC_pos(x + 1 * x_grid_dist,box_x)
    nw_z = PBC_pos(z - 1 * z_grid_dist,box_z)

    sw_x = PBC_pos(x - 1 * x_grid_dist,box_x)
    sw_z = PBC_pos(z - 1 * z_grid_dist,box_z)

    s_x = PBC_pos(x - 2 * x_grid_dist,box_x)
    s_z = z

    se_x = PBC_pos(x - 1 * x_grid_dist,box_x)
    se_z = PBC_pos(z + 1 * z_grid_dist,box_z)

    ne_x = PBC_pos(x + 1 * x_grid_dist,box_x)
    ne_z = PBC_pos(z + 1 * z_grid_dist,box_z)

    n_y, n = find_max_height_at_points(surface_lattice,full_depo_index,n_x,n_z)
    nw_y, nw = find_max_height_at_points(surface_lattice,full_depo_index,nw_x,nw_z)
    sw_y, sw = find_max_height_at_points(surface_lattice,full_depo_index,sw_x,sw_z)
    s_y, s = find_max_height_at_points(surface_lattice,full_depo_index,s_x,s_z)
    se_y, se = find_max_height_at_points(surface_lattice,full_depo_index,se_x,se_z)
    ne_y, ne = find_max_height_at_points(surface_lattice,full_depo_index,ne_x,ne_z)

    neighbour_species = [atom_below,n,nw,sw,ne,se,s]
    #print neighbour_species
    neighbour_pos = [x,y_max_0,z,n_x,n_y,n_z,nw_x,nw_y,nw_z,sw_x,sw_y,sw_z,s_x,s_y,s_z,se_x,se_y,se_z,ne_x,ne_y,ne_z]

    return neighbour_pos, neighbour_species

# find list of second neighbours (1-12)
# returns coordinates and species
def find_second_neighbours(x,z,full_depo_index):
    nb2 = []
    nb2_species = []

    # find x and z coordinates
    nx = PBC_pos(x + 4 * x_grid_dist,box_x)
    nz = z
    nb2.append([nx, nz])
    nx = PBC_pos(x + 3 * x_grid_dist,box_x)
    nz = PBC_pos(z - 1 * z_grid_dist,box_z)
    nb2.append([nx, nz])
    nx = PBC_pos(x + 2 * x_grid_dist,box_x)
    nz = PBC_pos(z - 2 * z_grid_dist,box_z)
    nb2.append([nx, nz])
    nx = x
    nz = PBC_pos(z - 2 * z_grid_dist,box_z)
    nb2.append([nx, nz])
    nx = PBC_pos(x - 2 * x_grid_dist,box_x)
    nz = PBC_pos(z - 2 * z_grid_dist,box_z)
    nb2.append([nx, nz])
    nx = PBC_pos(x - 3 * x_grid_dist,box_x)
    nz = PBC_pos(z - 1 * z_grid_dist,box_z)
    nb2.append([nx, nz])

    nx = PBC_pos(x - 4 * x_grid_dist,box_x)
    nz = z
    nb2.append([nx, nz])
    nx = PBC_pos(x - 3 * x_grid_dist,box_x)
    nz = PBC_pos(z + 1 * z_grid_dist,box_z)
    nb2.append([nx, nz])
    nx = PBC_pos(x - 2 * x_grid_dist,box_x)
    nz = PBC_pos(z + 2 * z_grid_dist,box_z)
    nb2.append([nx, nz])
    nx = x
    nz = PBC_pos(z + 2 * z_grid_dist,box_z)
    nb2.append([nx, nz])
    nx = PBC_pos(x + 2 * x_grid_dist,box_x)
    nz = PBC_pos(z + 2 * z_grid_dist,box_z)
    nb2.append([nx, nz])
    nx = PBC_pos(x + 3 * x_grid_dist,box_x)
    nz = PBC_pos(z + 1 * z_grid_dist,box_z)
    nb2.append([nx, nz])

    # find y and species
    for k in xrange(len(nb2)):
        nx, nz = nb2[k]
        ny, species = find_max_height_at_points(surface_lattice,full_depo_index,nx,nz)
        nb2[k] = [nx, ny, nz]
        nb2_species.append(species)

    return nb2, nb2_species

# write out new lattice
def write_lattice(index,full_depo_index,surface_lattice,natoms,time,barrier):
    new_lattice = output_dir_name_prefac + '/KMC' + str(index) + '.dat'
    outfile = open(new_lattice, 'w')
    line = str(natoms)
    outfile.write(line + '\n')
    line = str(time)
    outfile.write(line + '\n')
    line = str(barrier)
    outfile.write(line + '\n')

    outfile.write(str(box_x)+'  '+str(box_y)+'  '+str(box_z)+'  ' + '\n')

    # write surface lattice first (does not change)
    for j in xrange(len(surface_lattice)):
        #outfile.write(str(latticeLines[x]))
        outfile.write(str(surface_lattice[j][0]) + '   ' + str(surface_lattice[j][1])+ '    '+str(surface_lattice[j][2])+ '   '+ str(surface_lattice[j][3])+'  '+str(surface_lattice[j][4])+'\n')

    # write all deposited atom positions
    for i in xrange(len(full_depo_index)):
    	line = full_depo_index[i]
    	outfile.write(str(atom_species) + '	' + str(line[1]) + '   ' + str(line[2])+ '    '+ str(line[3]) + '   ' +  '0' + '\n')
    outfile.close

    return

# write temp lattice.dat file
def write_lattice_LKMC(index,full_depo_index,surface_lattice,natoms):
    new_lattice = NEB_dir_name_prefac + str(index) + '.dat'
    outfile = open(new_lattice, 'w')
    line = str(natoms)
    outfile.write(line + '\n')
    outfile.write(str(box_x)+'  '+str(30)+'  '+str(box_z)+'  ' + '\n')

    # write surface lattice first (does not change)
    for j in xrange(len(surface_lattice)):
        #outfile.write(str(latticeLines[x]))
        outfile.write(str(surface_lattice[j][0]) + '   ' + str(surface_lattice[j][1])+ '    '+str(surface_lattice[j][2])+ '   '+ str(surface_lattice[j][3])+'  '+str(surface_lattice[j][4])+'\n')

    # write all deposited atom positions
    for i in xrange(len(full_depo_index)):
    	line = full_depo_index[i]
    	outfile.write(str(atom_species) + '	' + str(line[1]) + '   ' + str(line[2])+ '    '+ str(line[3]) + '   ' +  '0' + '\n')
    outfile.close

    return

# move event
def move_atom(depo_list, dir_vector ,full_depo_index):
    moved_list = None
    x = depo_list[1]
    y = depo_list[2]
    z = depo_list[3]

    x = round(PBC_pos(x+dir_vector[0]*x_grid_dist,box_x),6)
    y = round(y+dir_vector[1]*y_grid_dist2,6)
    z = round(PBC_pos(z+dir_vector[2]*z_grid_dist,box_z),6)

    y2, neighbour_species, neighbour_heights = deposition_y(full_depo_index,x,z)
    #print neighbour_species

    # check if large up/down move has taken place. Then check for tripod of atoms
    if (np.abs(dir_vector[0])+np.abs(dir_vector[1])+np.abs(dir_vector[2])) > 3:
        if round(y-y2,2) > (y_grid_dist2*1.1):
            nH = neighbour_heights.count(round(y2,6))
            if nH < 3:
                #print "Moved to 'floating' unstable position", nH
                return None
    else:
        nH = neighbour_heights.count(round(y-y_grid_dist2,6))
        if nH < 3:
            nH = neighbour_heights.count(round(y-y_grid_dist,6))
            if nH < 3:
                #print "Moved to 'floating' unstable position", nH
                return None



    # check for move fails
    if round(y-y_grid_dist-neighbour_heights[0],2) == 0:
        #print "Moved to unstable position", neighbour_species[0]
        print x,y,z
        return None
    if round(y-y_grid_dist2-neighbour_heights[0],2) == 0:
        #print "Moved to unstable position", neighbour_species[0]
        print x,y,z
        return None

    if round(y-neighbour_heights[0],2) == 0:
        #print "Moved into existing atom!"
        print x,y,z
        return None

    AdNeighbours = 0
    for i in xrange(len(neighbour_species)):
            if round(neighbour_heights[0]-y,2) == 0:
                AdNeighbours += 1

    # if near neighbours exist in plane fail move
    if dir_vector[1] == 0:
        # will always be at least 1 as it includes self pre-move
        if AdNeighbours > 1:
            #print "Moved too close to existing atoms"
            print x,y,z
            return None
    else:
        if AdNeighbours > 0:
            #print "Moved too close to existing atoms"
            print x,y,z
            return None

    #print "Moved atom"

    moved_list = [atom_species,x,y,z,depo_list[4]]

    return moved_list


# pick an event from an event list
def choose_event(event_list,Time):
    # add on deposition event
    event_list.append([depoRate,0,['Depo'],0])

    # print "Choose event from event list:"
    # print event_list

    # create cumulative function
    R = []
    TotalRate = 0
    TotalBarrier = 0
    num = len(event_list)
    for i in xrange(num):
        if event_list[i][0] > 0 and event_list[i][0] != 'None':
            try:
                TotalRate += event_list[i][0]
                TotalBarrier += float(event_list[i][3])
                R.append([TotalRate,event_list[i][1],event_list[i][2],event_list[i][3]])
            except TypeError:
                print "WARNING! Type error in choose_event: "
                print "Rate: ",  event_list[i][0], "\tBarrier: ", event_list[i][3]
                sys.exit()
            except ValueError:
                print "WARNING! Value error in choose_event: "
                print "Rate: ",  event_list[i][0], "\tBarrier: ", event_list[i][3]
                sys.exit()
    TotalRate

    numEvents = len(event_list)


    # Debug
    #print R

    # find random number
    u = random.random()
    Q = u*TotalRate
    #print "DEBUG: random number: ", Q

    # choose event
    for i in xrange(len(R)):
        if Q < R[i][0]:
            chosenRate = R[i][0]
            chosenEvent = R[i][2]
            chosenAtom = R[i][1]
            chosenBarrier = R[i][3]
            print "Chosen event:",R[i][2],"on atom:",R[i][1]
            print "Rate:", R[i][0]
            break

    # increase time
    u = random.random()
    Time += (np.log(1/u)/TotalRate)*1E15

    return chosenRate, chosenEvent, chosenAtom, Time, chosenBarrier

# calcuate list of atoms with graph radius of defect
def find_volume_atoms(lattice_pos,x,y,z):
    volume_atoms = []

    for i in xrange(len(lattice_pos)/3):
        # find distance squared between 2 atoms
        dist = PBC_distance(lattice_pos[3*i],lattice_pos[3*i+1],lattice_pos[3*i+2],x,y,z)
        if dist < (graphRad * graphRad):
            volume_atoms.append(i)

    return volume_atoms

# calculate hashkey for a defect
def hashkey(lattice_positions,specie_list,volumeAtoms):
    # set up parameters for hashkey calculation

    #lattice.pos = np.asarray(lattice_positions,dtype=np.float64)
    Lattice1 = lattice()
    Lattice1.pos = []
    species = []
    for i in volumeAtoms:
        Lattice1.pos.append(lattice_positions[3*i])
        Lattice1.pos.append(lattice_positions[3*i+1])
        Lattice1.pos.append(lattice_positions[3*i+2])
        species.append(specie_list[i])
    Lattice1.pos = np.asarray(Lattice1.pos,dtype=np.float64)

    for i in xrange(len(species)):
        if species[i] == 'O_':
            species[i] = 0
        elif species[i] == 'Zn':
            species[i] = 1
        elif species[i] == 'Ag':
            species[i] = 2

    for i in xrange(len(specie_list)):
        if specie_list[i] == 'O_':
            specie_list[i] = 0
        elif specie_list[i] == 'Zn':
            specie_list[i] = 1
        elif specie_list[i] == 'Ag':
            specie_list[i] = 2
        # else:
        #     specie_list[i] = 3

    # set lattice object values for hashkey
    Lattice1.specie = np.asarray(species,np.int32)
    Lattice1.cellDims = np.asarray([box_x,0,0,0,MaxHeight,0,0,0,box_z],dtype=np.float64)
    Lattice1.specieList = np.asarray(['O_','Zn','Ag'],dtype=np.character)
    Lattice1.pos = np.around(Lattice1.pos,decimals = 5)

    volumeAtoms = np.arange(len(volumeAtoms),dtype=np.int32)

    params.graphRadius = graphRad

    # get hashkey
    hashkey = Graphs.getHashKeyForAVolume(params,volumeAtoms,Lattice1)

    del Lattice1
    return hashkey


# save defect volume to compare against
# - stores lattice + volume atom indices
def SaveVolume(hashkey,volumeAtoms,lattice_positions,specie_list):
    VolumeFile = Volumes_dir + '/'+str(hashkey)+'.txt'

    # check if file already exist for hashkey
    if (os.path.isfile(VolumeFile)):
        return

    else:
        outfile = open(VolumeFile, 'w')
        outfile.write(str(len(specie_list))+'\n')
        for i in xrange(len(specie_list)):
            outfile.write(str(specie_list[i]) + '   ' + str(lattice_positions[3*i])+ '    '+str(lattice_positions[3*i+1])+ '   '+ str(lattice_positions[3*i+2])+'\n')
        # print "length of volume atoms is", len(volumeAtoms)
        for j in xrange(len(volumeAtoms)):
            outfile.write(str(volumeAtoms[j]) +'\n')
        outfile.close()
        return



# find the final hashkey for a given defect and transition
def findFinal(dir_vector,atom_index,full_depo_index,surface_positions):
    adatom_positions = []
    adatom_specie = []
    full_depo = copy.deepcopy(full_depo_index)
    depo_list = full_depo[atom_index]

    # move atom to final position
    moved_list = move_atom(depo_list, dir_vector ,full_depo)

    if moved_list:
        full_depo[atom_index] = moved_list
        depo_list = full_depo[atom_index]

        for i in xrange(len(full_depo)):
            adatom_specie.append(full_depo[i][0])
            adatom_positions.append(full_depo[i][1])
            adatom_positions.append(full_depo[i][2])
            adatom_positions.append(full_depo[i][3])
        lattice_positions = surface_positions + adatom_positions
        specie_list = surface_specie + adatom_specie

        # find atoms in defect volume
        volumeAtoms = find_volume_atoms(lattice_positions,depo_list[1],depo_list[2],depo_list[3])

        # create hashkey
        final_key = hashkey(lattice_positions,specie_list,volumeAtoms)

        #write_lattice(1000,full_depo_index,surface_lattice,401,0,0)
        del full_depo
        return final_key, [float(depo_list[1]),float(depo_list[2]),float(depo_list[3])]
    else:
        return None, None

# create the list of possible events
def create_events_list(full_depo_index,surface_lattice, volumes):

    event_list = []
    adatom_positions = []
    adatom_specie = []

    # move to global parameters
    trans_dir = initial_dir + '/Transitions/'

    # add all deposited atoms to list of positions + species
    for i in xrange(len(full_depo_list)):
        adatom_specie.append(full_depo_list[i][0])
        adatom_positions.append(full_depo_list[i][1])
        adatom_positions.append(full_depo_list[i][2])
        adatom_positions.append(full_depo_list[i][3])
    lattice_positions = surface_positions + adatom_positions
    specie_list = surface_specie + adatom_specie

    # find transitions for each adatom
    for j in xrange(len(full_depo_list)):
        final_keys = []
        directions = []
        depo_list = full_depo_list[j]
        hashkeyExists = []

        # find atoms in volume
        volumeAtoms = find_volume_atoms(lattice_positions,depo_list[1],depo_list[2],depo_list[3])

        # create hashkey for each adatom + store volume
        vol_key = hashkey(lattice_positions,specie_list,volumeAtoms)
        # print vol_key
        try:
            vol = volumes[vol_key]
        except KeyError:
            volumes[vol_key] = volume()
            vol = volumes[vol_key]

        if useBasin:
            basinExists = False
            iniPos = [depo_list[1],depo_list[2],depo_list[3]]
            # check if a basin exists for this state
            for basId in basinList:
                if basId.atomNum == j:
                    if basId.thisBasin(iniPos):
                        bas = basId
                        basinExists = True
                        keepBasin = True
                        break

            if not basinExists:
                bas = basin()
                bas.atomNum = j
                bas.currentPos = iniPos
                basinList.append(bas)
                keepBasin = False



        if len(vol.directions) != 0:
            vol = volumes[vol_key]
            # vol.hashkey = vol_key
            for direc in vol.directions:
                final_key, final_pos = findFinal(direc,j,full_depo_index,surface_positions)
                if final_key is not None:
                    try:
                        trans = vol.finalKeys[final_key]

                        if useBasin:
                            bas.addTransition(iniPos,final_pos,trans.rate,trans.barrier,trans.reverseBarrier)
                            if trans.barrier < basinBarrierTol or trans.reverseBarrier < basinBarrierTol:
                                keepBasin = True

                        # trans.hashkey = final_key
                        # event_list.append([trans.rate,j,final_pos,trans.barrier])
                    except KeyError:
                        result, vol = singleNEB(direc,full_depo_index,surface_lattice,j,vol_key,final_key,natoms,vol)
                        if result:
                            if result[2] != "None":
                                rate = calc_rate(float(result[2]))
                                bas.addTransition(iniPos,final_pos,rate,float(result[2]),vol.finalKeys[final_key].reverseBarrier)
                                if float(result[2]) < basinBarrierTol or vol.finalKeys[final_key].reverseBarrier < basinBarrierTol:
                                    keepBasin = True
                                # event_list.append([rate,j,final_pos,float(result[2])])


        else:
            print "Cannot find volume transitions. Doing searches now ", vol_key
            # do searches on volume and save to new trans file
            if useBasin:
                status, result, vol, keepBasin = autoNEB(full_depo_index,surface_lattice,j,vol_key,natoms,vol,bas)
            else:
                status, result, vol, _ = autoNEB(full_depo_index,surface_lattice,j,vol_key,natoms,vol,None)
            if status:
                break
            else:
                event_list = event_list + result
                volumes[vol_key] = vol

        del volumeAtoms

        if useBasin:
            if not keepBasin:
                events = bas.addUnchangedEvents(j)
                event_list = event_list + events
                basinList.pop()
            else:
                bas.buildConnectivity()
                events = bas.addChangedEvents(j)
                event_list = event_list + events


    del lattice_positions
    del adatom_positions

    return event_list, volumes

# check if event exists
def add_to_events(final_keys,hashkey,barrier,index,directions,hashkeyExists):
    events = []

    for k in xrange(len(final_keys)):
        if hashkey == final_keys[k]:
            if barrier != "None":
                rate = calc_rate(float(barrier))
                events.append([rate,index,directions[k],barrier])
            hashkeyExists[k] = 1

    return events, hashkeyExists


# run NEB to find barriers that are not known
def autoNEB(full_depo_index,surface_lattice,atom_index,hashkey,natoms,vol,bas):
    print "AUTO NEB", "="*60

    barrier = []
    final_keys = []
    results = []

    keepBasin = False

    # create initial lattice
    write_lattice_LKMC('/initial',full_depo_index,surface_lattice,natoms)
    ini = Lattice.readLattice(NEB_dir_name_prefac+"/initial.dat")
    iniMin = copy.deepcopy(ini)
    iniMin.calcForce(correctTE=1)
    print "ini energy:", iniMin.totalEnergy

    # create cell dimensions
    cellDims = np.asarray([box_x,0,0,0,MaxHeight,0,0,0,box_z],dtype=np.float64)

    # minimise lattice
    minimiser = Minimise.getMinimiser(params)
    status = minimiser.run(iniMin)
    if status:
        print " Warning: failed to minimise initial lattice"
        sys.exit()

    # check max movement
    Index, maxMove, avgMove, Sep = Vectors.maxMovement(ini.pos, iniMin.pos, cellDims)

    if maxMove < maxMoveCriteria:
        #ini.writeLattice("initialMin.dat")
        full_depo = copy.deepcopy(full_depo_index)
        depo_list = full_depo[atom_index]

        # check 6 initial directions
        dir_vector = []
        dir_vector.append([2,0,0])
        dir_vector.append([1,0,-1])
        dir_vector.append([-1,0,-1])
        dir_vector.append([1,0,1])
        dir_vector.append([-1,0,1])
        dir_vector.append([-2,0,0])

        # include down transitions
        if IncludeDownTrans:
            # if atom if above first layer
            atom_height = full_depo_index[atom_index][2]
            if atom_height > (initial_surface_height + y_grid_dist*1.1):
                print "Adding move down transitions"
                dir_vector.append([4,-1,0])
                dir_vector.append([2,-1,-2])
                dir_vector.append([-2,-1,-2])
                dir_vector.append([2,-1,2])
                dir_vector.append([-2,-1,2])
                dir_vector.append([-4,-1,0])

        # include up transitions
        if IncludeUpTrans:
            #check if surrounds atom
            AdNeighbours = 0
            nb_pos, nb_species = find_second_neighbours(full_depo_index[atom_index][1],full_depo_index[atom_index][3],full_depo_index)
            for j in xrange(len(nb_pos)):
                if round(nb_pos[j][1] - atom_height,2) == 0:
                    if nb_species[j] == atom_species:
                        AdNeighbours += 1

            # if 2 or more atoms surround current atom, look at up moves
            if AdNeighbours > 1:
                print "Adding move up transitions"
                dir_vector.append([4,1,0])
                dir_vector.append([2,1,-2])
                dir_vector.append([-2,1,-2])
                dir_vector.append([2,1,2])
                dir_vector.append([-2,1,2])
                dir_vector.append([-4,1,0])


        for i in xrange(len(dir_vector)):
            # move atom
            moved_list = move_atom(depo_list, dir_vector[i] ,full_depo_index)
            print "Trying direction: ", dir_vector[i]
            vol.addDirection(dir_vector[i])

            if moved_list:
                full_depo[atom_index] = moved_list

                # create initial lattice and Minimise
                write_lattice_LKMC('/'+str(i),full_depo,surface_lattice,natoms)
                fin = Lattice.readLattice(NEB_dir_name_prefac+"/" + str(i) +".dat")
                finMin = copy.deepcopy(fin)

                finMin.calcForce(correctTE=1)
                print "fin energy: ", finMin.totalEnergy

                # minimise lattice
                mini_fin = Minimise.getMinimiser(params)
                status = mini_fin.run(finMin)
                if status:
                    continue

                final_key, final_pos = findFinal(dir_vector[i],atom_index,full_depo_index,surface_positions)

                # check that initial and final are different
                Index, maxMove, avgMove, Sep = Vectors.maxMovement(iniMin.pos, finMin.pos, cellDims)
                if maxMove < 0.4:
                    print " difference between ini and fin is too small:", maxMove
                    barrier = str("None")
                    results.append([0,atom_index, final_pos, barrier])
                    vol.addTrans(dir_vector[i], final_key, barrier, 0, str("None"))
                    continue

                # check max movement
                Index, maxMove, avgMove, Sep = Vectors.maxMovement(fin.pos, finMin.pos, cellDims)
                if maxMove < maxMoveCriteria:
                    # run NEB on initial and final lattices
                    neb = NEB.NEB(params)
                    status = neb.run(iniMin, finMin)

                    if status:
                        print "WARNING: NEB failed to converge"
                        print "Try changing parameters in lkmcInput.IN"
                        barrier = str("None")

                        results.append([0,atom_index, final_pos, barrier])
                        vol.addTrans(dir_vector[i], final_key, barrier, str("None"))
                        continue

                    neb.barrier = round(neb.barrier,6)
                    reverseBarrier = round((iniMin.totalEnergy-finMin.totalEnergy)+neb.barrier,6)
                    print "Reverse barrier: ", reverseBarrier
                    # find final hashkey
                    final_key, final_pos = findFinal(dir_vector[i],atom_index,full_depo_index,surface_positions)
                    rate = calc_rate(neb.barrier)
                    results.append([rate, atom_index, final_pos, neb.barrier])
                    vol.addTrans(dir_vector[i], final_key, neb.barrier, rate, reverseBarrier)

                    # add result to basin
                    if useBasin:
                        iniPos = copy.copy(full_depo_index[atom_index])
                        iniPos.pop(0)
                        iniPos.pop()
                        bas.addTransition(iniPos,final_pos,rate,neb.barrier,reverseBarrier)
                        if neb.barrier < basinBarrierTol or reverseBarrier < basinBarrierTol :
                            keepBasin = True

                else:
                    print "WARNING: maxMove too large in final lattice:", maxMove
                    barrier = str("None")
                    results.append([0,atom_index, final_pos, barrier])
                    vol.addTrans(dir_vector[i], final_key, barrier, 0, str("None"))

        if results:
            print results
            #write_trans_file(hashkey,results)
            return 0, results, vol, keepBasin;
    else:
        print "WARNING: maxMove too large in initial lattice"
        sys.exit()

    del ini, iniMin
    del fin, finMin
    del Sep

    return 1, results, vol, keepBasin;

# do a single NEB and add transition to trans files
def singleNEB(direction,full_depo_index,surface_lattice,atom_index,hashkey,final_key,natoms, vol):
    print "SINGLE NEB", "="*60

    barrier = []
    final_keys = []
    results = None

    print natoms
    # create initial lattice
    write_lattice_LKMC('/initial',full_depo_index,surface_lattice,natoms)
    ini = Lattice.readLattice(NEB_dir_name_prefac+"/initial.dat")
    iniMin = copy.deepcopy(ini)

    # create cell dimensions
    cellDims = np.asarray([box_x,0,0,0,MaxHeight,0,0,0,box_z],dtype=np.float64)

    # minimise lattice
    minimiser = Minimise.getMinimiser(params)
    status = minimiser.run(iniMin)
    if status:
        print " WARNING! failed to minimise initial lattice"
        sys.exit()

    # check max movement
    Index, maxMove, avgMove, Sep = Vectors.maxMovement(ini.pos, iniMin.pos, cellDims)

    if maxMove < maxMoveCriteria:
        #ini.writeLattice("initialMin.dat")
        full_depo = copy.deepcopy(full_depo_index)
        depo_list = full_depo[atom_index]

        # move atom
        moved_list = move_atom(depo_list, direction ,full_depo_index)

        if moved_list:
            full_depo[atom_index] = moved_list

            # create initial lattice and Minimise
            write_lattice_LKMC('/6',full_depo,surface_lattice,natoms)
            fin = Lattice.readLattice(NEB_dir_name_prefac+"/" + '6' +".dat")
            finMin = copy.deepcopy(fin)
            Index, maxMove, avgMove, Sep = Vectors.maxMovement(iniMin.pos, finMin.pos, cellDims)
            if maxMove < 0.4:
                print " difference between ini and fin is too small", maxMove
                barrier = str("None")
                results = [direction, final_key, barrier]
                results[0] = map(int,results[0])
                del ini, iniMin
                del fin, finMin
                vol.addTrans(results[0], final_key, barrier, 0, str("None"))
                # add_to_trans_file(hashkey,results)
                return results, vol

            # minimise lattice

            mini_fin = Minimise.getMinimiser(params)
            status = mini_fin.run(finMin)
            if status:
                print " Failed to minimise final"
                barrier = str("None")
                results = [direction, final_key, barrier]
                results[0] = map(int,results[0])
                del ini, iniMin
                del fin, finMin
                vol.addTrans(results[0], final_key, barrier, 0, str("None"))
                # add_to_trans_file(hashkey,results)
                return results, vol

            # check max movement
            Index, maxMove, avgMove, Sep = Vectors.maxMovement(fin.pos, finMin.pos, cellDims)
            if maxMove < maxMoveCriteria:
                # run NEB on initial and final lattices
                neb = NEB.NEB(params)
                status = neb.run(iniMin, finMin)

                if status:
                    print "WARNING: NEB failed to converge"
                    print "Try changing parameters in lkmcInput.IN"
                    barrier = str("None")
                    results = [direction, final_key, barrier]
                    results[0] = map(int,results[0])
                    del ini, iniMin
                    del fin, finMin
                    vol.addTrans(results[0], final_key, barrier, 0, str("None"))
                    #add_to_trans_file(hashkey,results)
                    return results, vol

                print neb.barrier
                barrier = round(neb.barrier,6)
                reverseBarrier = round((iniMin.totalEnergy-finMin.totalEnergy)+neb.barrier,6)
                print "Reverse barrier: ", reverseBarrier

                results = [direction, final_key, barrier]
                results[0] = map(int,results[0])
                rate = calc_rate(barrier)
                vol.addTrans(results[0], final_key, barrier, rate, reverseBarrier)
                print direction
            else:
                print "WARNING: maxMove too large in final lattice"
                barrier = str("None")
                results = [direction, final_key, barrier]
                results[0] = map(int,results[0])
                del ini, iniMin
                del fin, finMin
                vol.addTrans(results[0], final_key, barrier, 0, str("None"))
                #add_to_trans_file(hashkey,results)
                return results, vol
        else:
            print "WARNING: move atom failed"
            barrier = str("None")
            results = [direction, final_key, barrier]
            results[0] = map(int,results[0])
            del ini, iniMin
            vol.addTrans(results[0], final_key, barrier, 0 , str("None"))
            #add_to_trans_file(hashkey,results)
            return results, vol

        # if results:
        #     add_to_trans_file(hashkey,results)
    else:
        print "WARNING: maxMove too large in initial lattice"
        sys.exit()

    del ini, iniMin
    del Sep
    print "Finished SINGlE NEB", "="*60
    return results, vol

# write out all volumes and transitions to a file
def writeVolumes(volumes):
    volfile = initial_dir + '/Volumes.txt'

    out = open(volfile, 'w')
    for i, vol in enumerate(volumes):
        # print "vol:", vol
        numDir = len(volumes[vol].directions)
        numTrans = len(volumes[vol].finalKeys)
        out.write(str(vol)+'\t'+str(numDir)+'\t'+str(numTrans)+'\n')
        for j in range(numDir):
            direc = volumes[vol].directions[j]
            out.write(str(direc[0])+'\t'+str(direc[1])+'\t'+str(direc[2])+'\n')
        for j, trans in enumerate(volumes[vol].finalKeys):
            out.write(str(trans)+'\t'+str(volumes[vol].finalKeys[trans].barrier)+'\t'+str(volumes[vol].finalKeys[trans].rate)+'\t'+str(volumes[vol].finalKeys[trans].reverseBarrier)+'\n')

    # writeVolAtoms(volumes)
    return

# write out volume atom positions to a file
def writeVolAtoms(volumes):
    volfile = initial_dir + '/VolumeAtoms.txt'

    out = open(volfile, 'w')
    for i, vol in enumerate(volumes):
        numAtoms = len(volumes[vol].volumeAtoms)
        out.write(str(vol)+'\t'+str(numAtoms)+'\n')
        cV = volumes[vol]
        for j in range(numAtoms):
            out.write(str(cV.specie[j]) + '   ' + str(cV.pos[3*j])+ '    '+str(cV.pos[3*j+1])+ '   '+ str(cV.pos[3*j+2])+'\n')
    return

# read volumes from file
def readVolumes(volumes):

    # read in transitions
    volPath = initial_dir + '/Volumes.txt'
    if (os.path.isfile(volPath)):
        input_file = open(volPath, 'r')
        while 1:
            latline = input_file.readline()
            line = latline.split()
            if len(line) < 3:
                break
            key = str(line[0])
            numDir = int(line[1])
            numTrans = int(line[2])
            vol = volume()
            for i in range(numDir):
                latline = input_file.readline()
                line = latline.split()
                if len(line) < 3:
                    break
                direc = [int(line[0]),int(line[1]),int(line[2])]
                vol.directions.append(direc)
            for i in range(numTrans):
                latline = input_file.readline()
                line = latline.split()
                if len(line) < 3:
                    break
                if line[1]!= "None":
                    vol.addTrans(direc, str(line[0]), float(line[1]), float(line[2]), float(line[3]))
                else:
                    vol.addTrans(direc, str(line[0]), str(line[1]), float(0.0), str(line[3]))
            volumes[key]=vol

    return volumes

# write stats to a file
def StatsOutput(event_list,CurrentStep,numAdatoms):
    statsFile = Stats_dir + '/Stats.txt'
    if (os.path.isfile(statsFile)):
        outfile = open(statsFile, 'a')
        TotalRate = 0
        TotalBarrier = 0
        num = len(event_list)
        for i in xrange(num):
            TotalRate += float(event_list[i][0])
            if event_list[i][3] != None:
                try:
                    TotalBarrier += float(event_list[i][3])
                except ValueError:
                    continue

        AveRate = TotalRate/num
        AveBarrier = TotalBarrier/num
        outfile.write(str(AveRate)+','+str(AveBarrier)+','+str(len(event_list))+','+str(numAdatoms)+','+str(CurrentStep)+'\n')
        outfile.close()
    else:
        print "Warning! Could not find Stats.txt"
    return






# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ====================== MAIN SCRIPT =========================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================


# set initial values
Time = 0
full_depo_list = []
surface_specie= []
surface_positions = []
startTimeSub = time.time()
CurrentStep = 0
volumes = {}
basinList = []

print "="*80
print "~~~~~~~~ Starting lattice KMC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "="*80


# directory setup
initial_dir = os.getcwd()
output_dir_name_prefac = initial_dir + '/Output'
Trans_dir = initial_dir + '/Transitions'
Volumes_dir = initial_dir + '/Volumes'
NEB_dir_name_prefac = initial_dir + '/Temp'
Stats_dir = initial_dir + '/Stats'

print " Current directory           : ", initial_dir
print " Output directory prefactor  : ", output_dir_name_prefac

# determine if all input files are present
# check we have lattice.dat
if (os.path.isdir(initial_dir)):
    # check for required files
    # lattice.dat
    check_path = initial_dir + '/lattice.dat'
    input_lattice_path = check_path
    if not (os.path.isfile(check_path)):
        print "lattice.dat: ", check_path, " not found, exiting ..."
        sys.exit()

# check if volumes and output directories exists, if not: creates them
if not os.path.exists(Volumes_dir):
    os.makedirs(Volumes_dir)
if not os.path.exists(output_dir_name_prefac):
    os.makedirs(output_dir_name_prefac)
if not os.path.exists(NEB_dir_name_prefac):
    os.makedirs(NEB_dir_name_prefac)
if not os.path.exists(Trans_dir):
    os.makedirs(Trans_dir)
if StatsOut:
    if not os.path.exists(Stats_dir):
        os.makedirs(Stats_dir)
    statsFile = Stats_dir + '/Stats.txt'
    outfile = open(statsFile, 'w')
    outfile.write('Average Rate'+', Average Barrier'+', No. Events'+', No. Adatoms'+', Step'+'\n')
    outfile.close()
print "~"*80

# read lattice header
natoms,box_x,box_y,box_z = read_lattice_header(input_lattice_path)
print "Initial atoms  : ",natoms
print "Lattice size: ",box_x,box_y,box_z, " Angstroms"
sys.stdout.flush()

# Y-height of initial lattice
initial_surface_height = find_max_height(input_lattice_path)
print "initial_surface_height: ", initial_surface_height

# store whole surface lattice first
surface_lattice = read_lattice(input_lattice_path,2)
for i in xrange(len(surface_lattice)):
    surface_specie.append(surface_lattice[i][0])
    surface_positions.append(round(surface_lattice[i][1],6))
    surface_positions.append(round(surface_lattice[i][2],6))
    surface_positions.append(round(surface_lattice[i][3],6))

# set up temp initial and final lattices.dat
params = Input.getLKMCParams(1, "", "lkmcInput.IN")
Input.readGlobals("lkmcInput.IN")
volumes = readVolumes(volumes)

# find size of grid_size
x_grid_points, y_grid_points, z_grid_points, box_x, box_z = grid_size(box_x,initial_surface_height,box_z)
print "grid size: %d * %d * %d" % (x_grid_points,y_grid_points,z_grid_points)
if (x_grid_points%6):
    print " WARNING WARNING: x grid points must be a multiple of 6 for PBC to work"
    sys.exit()
if (z_grid_points%2):
    print " WARNING WARNING: z grid points must be a multiple of 2 for PBC to work"
    sys.exit()
print "New lattice size: ",box_x,box_y,box_z, " Angstroms"
print "-" * 80

# check if continue or begin run
if jobStatus == 'CNTIN':
    num = 0
    orig_len = len(surface_lattice)
    while 1:
        kmcFile = output_dir_name_prefac + '/KMC' + str(num) + '.dat'

        # read in last KMC file
        if not os.path.exists(kmcFile):
            input_file = open(output_dir_name_prefac + '/KMC' + str(num-1) + '.dat', 'r')
            # skip first line
            latline = input_file.readline()
            latline = input_file.readline()
            Time = float(latline)
            input_file.close()

            full_depo_list = read_lattice(output_dir_name_prefac + '/KMC' + str(num-1) + '.dat',orig_len+4)

            # add adatoms to full_depo_list
            for i in xrange(len(full_depo_list)):
                full_depo_list[i][4] = orig_len + 1 + i

            natoms += len(full_depo_list)
            CurrentStep = (num-1) * latticeOutEvery
            break
        num += 1

# do initial consecutive depositions
while CurrentStep < (numberDepos):
    depo_list = []
    depo_list = deposition(box_x,box_z,x_grid_dist,z_grid_dist,full_depo_list,natoms)
    if depo_list:
        print "Current Step: ", CurrentStep
        natoms = depo_list[4]
        full_depo_list.append(depo_list)
        write_lattice(CurrentStep,full_depo_list,surface_lattice,natoms,0,0)
        print "Writing lattice: KMC 0"
        CurrentStep += 1

print "-" * 80
print "Number of initial adatoms: " , len(full_depo_list)
index = CurrentStep

# store previous 2 list of positions and events
full_depo_list1 = [0,0]
full_depo_list2 = [0,0]


# ============================================================================
# ================== do KMC run ==============================================
# ============================================================================


while CurrentStep < (total_steps + 1):
    # TODO: include a verbosity level
    print "Current Step: ", CurrentStep

    # check if in same position as 2 steps ago
    if full_depo_list != full_depo_list2[1]:
        event_list, volumes = create_events_list(full_depo_list,surface_lattice, volumes)
    else:
        # do quick reuse
        event_list = full_depo_list2[0]

    full_depo_list2 = copy.deepcopy(full_depo_list1)
    full_depo_list1[1] = copy.deepcopy(full_depo_list)
    full_depo_list1[0] = copy.deepcopy(event_list)

    # write out volumes file
    if CurrentStep%volumesOutEvery == 0 or CurrentStep == total_steps:
        writeVolumes(volumes)

    # write out stats
    if StatsOut:
        StatsOutput(event_list,CurrentStep,len(full_depo_list))

    # choose event
    chosenRate, chosenEvent, chosenAtom, Time, chosenBarrier = choose_event(event_list, Time)

    # do deposition
    if chosenEvent[0] == 'Depo':
        while index < (CurrentStep+1):
            depo_list = []
            depo_list = deposition(box_x,box_z,x_grid_dist,z_grid_dist,full_depo_list,natoms)
            if depo_list:
                natoms = depo_list[4]
                full_depo_list.append(depo_list)
                index += 1
        CurrentStep += 1

    # do move
    else:
        while index < (CurrentStep+1):
            full_depo_list[chosenAtom] = [full_depo_list[chosenAtom][0], chosenEvent[0],chosenEvent[1],chosenEvent[2],full_depo_list[chosenAtom][4]]
            index += 1
        CurrentStep += 1


    # write out lattice
    if (CurrentStep-1)%latticeOutEvery == 0:
        latticeNo = (CurrentStep-1)/latticeOutEvery
        print "Writing lattice: KMC", latticeNo

        Barrier = chosenBarrier
        write_lattice(latticeNo,full_depo_list,surface_lattice,natoms,Time,Barrier)


    print "-" * 80


# last lines of output
print "====== Finished KMC run ========================================================"
print "Number of steps completed:	", (CurrentStep-1)
print "Number of depositions:		", len(full_depo_list)
print "="*80
print "Lattice KMC        Copyright (c) Adam Lloyd 2016 "
print "="*80

# Timer
FinalTimeSub = time.time() - startTimeSub
print "Time: ", FinalTimeSub
print "Average Time per step: ", FinalTimeSub/CurrentStep

if StatsOut:
    if (os.path.isfile(statsFile)):
        input_file = open(statsFile, 'r')
        # skip first line
        line = input_file.readline()
        AveRate = 0
        AveBarrier = 0
        AveEvents = 0
        while 1:
            line = input_file.readline()
            if not line: break
            #print line
            line = line.split(',')
            AveRate += float(line[0])
            AveBarrier += float(line[1])
            AveEvents += float(line[2])
        input_file.close()
        AveRate = AveRate/(CurrentStep-numberDepos)
        AveBarrier = AveBarrier/(CurrentStep-numberDepos)
        AveEvents = AveEvents/(CurrentStep-numberDepos)
        print "Average Rate: ", AveRate, "\tAverage Barrier: ", AveBarrier, "\tAverage Number of events: ", AveEvents


del volumes
del full_depo_list
del surface_lattice

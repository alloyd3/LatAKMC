# -*- coding: utf-8 -*-

# This script runs KMC on a hexagonal ZnO(0001) surface
# TODO: make more general for any hexagonal surface
# Assumes PBC in x and z directions and deposition in y
# grid positions start at (0,0,0)


# copyright Adam Lloyd 2016

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
jobStatus = 'BEGIN'            # BEGIN or CNTIN (not implemented yet) run
atom_species = 'Ag'         # species to deposit
numberDepos = 2		        # number of initial depositions
total_steps = 3          # total number of steps to run
latticeOutEvery = 1         # write output lattice every n steps
temperature = 300           # system temperature in Kelvin
prefactor = 1.00E+13        # fixed prefactor for Arrhenius eq. (typically 1E+12 or 1E+13)
boltzmann = 8.62E-05        # Boltzmann constant (8.62E-05)
graphRad = 5.9                # graph radius of defect volumes (Angstroms)
depoRate = 1.2e3            # deposition rate
maxMoveCriteria = 0.6        # maximum distance an atom can move after relaxation (pre NEB)

# for (0001) ZnO only
x_grid_dist = 0.9497411251   # distance in x direction between each atom in lattice
# TODO: second layer increase to 2.1
y_grid_dist = 1.55            # distance in y direction between each atom in lattice
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

# calculate the rate of an event given barrier height (Arrhenius eq.)
def calc_rate(barrier):
    rate = prefactor * np.exp(-barrier/(boltzmann*temperature))
    return rate

# find barrier height given rate
def find_barrier_height(rate):
    barrier = -(boltzmann*temperature)*np.log(rate/prefactor)
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
    for i in xrange(len(full_lattice)):
        line = full_lattice[i]
        if(float(line[1]) > (x-0.1) and float(line[1]) < (x+0.1)):
            if(float(line[3]) > (z-0.1) and float(line[3]) < (z+0.1)):
                if(float(line[2]) > max_height):
                    max_height = float(line[2])
                    max_height_atom = str(line[0])
            elif((float(line[3])+box_z) > (z-0.1) and (float(line[3])+box_z) < (z+0.1)):
                if(float(line[2]) > max_height):
                    max_height = float(line[2])
                    max_height_atom = str(line[0])
        elif((float(line[1])+box_x) > (x-0.1) and (float(line[1])+box_x) < (x+0.1)):
            if(float(line[3]) > (z-0.1) and float(line[3]) < (z+0.1)):
                if(float(line[2]) > max_height):
                    max_height = float(line[2])
                    max_height_atom = str(line[0])
            elif((float(line[3])+box_z) > (z-0.1) and (float(line[3])+box_z) < (z+0.1)):
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
    full_lattice = surface_lattice + full_depo_index
    for i in xrange(len(full_lattice)):
        line = full_lattice[i]
        if(float(line[1]) > (x-0.1) and float(line[1]) < (x+0.1)):
            if(float(line[3]) > (z-0.1) and float(line[3]) < (z+0.1)):
                if(float(line[2]) > max_height):
                    atom_below = max_height_atom
                    max_height = float(line[2])
                    max_height_atom = str(line[0])
            elif((float(line[3])+box_z) > (z-0.1) and (float(line[3])+box_z) < (z+0.1)):
                if(float(line[2]) > max_height):
                    atom_below = max_height_atom
                    max_height = float(line[2])
                    max_height_atom = str(line[0])
        elif((float(line[1])+box_x) > (x-0.1) and (float(line[1])+box_x) < (x+0.1)):
            if(float(line[3]) > (z-0.1) and float(line[3]) < (z+0.1)):
                if(float(line[2]) > max_height):
                    atom_below = max_height_atom
                    max_height = float(line[2])
                    max_height_atom = str(line[0])
            elif((float(line[3])+box_z) > (z-0.1) and (float(line[3])+box_z) < (z+0.1)):
                if(float(line[2]) > max_height):
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
    	neighbour_heights.append(round(neighbour_pos[i*3+1],7))
    	i += 1
#     print neighbour_heights


    # find max height of local points
    y_max = max(neighbour_heights)

    # add y dist to max height found
    y_coordinate = round(y_max,7)
    return y_coordinate, neighbour_species, neighbour_heights

# do deposition
def deposition(box_x,box_z,x_grid_dist,z_grid_dist,full_depo_index,natoms):
    x_coord, z_coord = deposition_xz(box_x,box_z,x_grid_dist,z_grid_dist)
    y_coord, nlist, hlist = deposition_y(full_depo_index,x_coord,z_coord)

    maxAg = []
    nAg = nlist.count('Ag')
    nH = hlist.count(y_coord)
    if nlist[0] != 'Ag':
        for i in xrange(len(nlist)):
            if nlist[i] == 'Ag' and hlist[i] == y_coord:
                maxAg.append(i)
                if len(maxAg) == 3:
                    print "landed on top of 3 Ag atoms"

    y_coord += y_grid_dist
    if (y_coord < initial_surface_height):
        print "ERROR y height is less than initial surface height", y_coord
        return


    print "Trying to deposit %s atom at %f, %f, %f" % (atom_species, x_coord, y_coord, z_coord)
    #print nlist
    if len(maxAg) != 3:
        for x in nlist:
            if x == 'Ag':
            	# TODO: take out later
                print "deposited near Ag, redo deposition"
                return None
        if nlist[0] == 'O_':
        	print "deposited on top of surface Oxygen: unstable position"
        	return None

    natoms += 1

    print "SUCCESS: Number of atoms: ", natoms

    deposition_list = [atom_species, x_coord, y_coord, z_coord, natoms]
    return deposition_list

# find PBC position in one direction
def PBC_pos(x,box_x):
	x = round(x,7)
	box_x = round(box_x,7)
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
    #neighbour_pos = [x,y_max_0,z,n_x,n_y,n_z,nw_x,nw_y,nw_z,sw_x,sw_y,sw_z,ne_x,ne_y,ne_z,se_x,se_y,se_z,s_x,s_y,s_z]

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
    x = depo_list[1]
    y = depo_list[2]
    z = depo_list[3]

    x = round(PBC_pos(x+dir_vector[0]*x_grid_dist,box_x),7)
    y = round(y+dir_vector[1]*y_grid_dist,7)
    z = round(PBC_pos(z+dir_vector[2]*z_grid_dist,box_z),7)


    y2, neighbour_species, neighbour_heights = deposition_y(full_depo_index,x,z)
    #print neighbour_species


    # check if Ag moved to unstable sites
    # TODO: check if on top of another Ag
    if neighbour_species[0] == 'O_':
    	print "moved on top of surface Oxygen: unstable position"
    	return None
    #print "Moved atom"

    moved_list = [atom_species,x,y,z,depo_list[4]]
    return moved_list

#def move_atom(depo_list, direction_index,full_depo_index):

# delete
# create a roulette table
def create_roulette(AtomIndex,x_coord,z_coord,second_nb_pos, second_nb_species, full_depo_index):
    # find probability for atom to move in any of 6 directions
    # max height at single point
    y_max_0, atom_below = find_atom_below(surface_lattice,full_depo_index,x_coord,z_coord)
    #print atom_below

    # check height at surrounding points
    nb_pos, nb_species = find_neighbours(x_coord,z_coord,atom_below,y_max_0,full_depo_index)

    # define an event list : rate, move_direction
    event_list = []
    for i in xrange(6):
        if nb_species[i+1] == 'O_':
            if nb_pos[3*(i+1)+1] == initial_surface_height:
                event_list.append([0,i,AtomIndex])
            else:
                # TODO: figure out prob if on silver
                event_list.append([0,i,AtomIndex])


        elif 'Ag' not in second_nb_species:
            if nb_species[0] == 'Zn':
                event_list.append([3e10,i,AtomIndex])
            elif nb_species[0] == None:
                event_list.append([2.7e4,i,AtomIndex])
            else:
                print "Error in roulette: Atom below must be Ag"
        else:
            event_list.append([0,i,AtomIndex])

            # near_second_neighbours = [second_nb_species[(2*i-2)%12],second_nb_species[(2*i-1)%12],second_nb_species[(2*i)%12]]
            # for j  in xrange(len(near_second_neighbours)):
            #     if 'Ag' in near_second_neighbours:
            #         # TODO: could climb or drop with positive rate
            #         event_list.append([0,i])
            #         break

    return event_list

# pick an event from an event list
def choose_event(event_list):
    # add on deposition event
    event_list.append([depoRate,0,['Depo']])

    # create cumulative function
    R = []
    TotalRate = 0
    for i in xrange(len(event_list)):
        TotalRate += event_list[i][0]
        R.append([TotalRate,event_list[i][1],event_list[i][2]])
    TotalRate

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
            print "Chosen event:",R[i][2],"on atom:",R[i][1]
            print "Rate:", event_list[i][0]
            break

    # TODO: increase time

    return chosenRate, chosenEvent, chosenAtom

# calcuate list of atoms with graph radius of defect
def find_volume_atoms(lattice_pos,x,y,z):
    volume_atoms = []
    for i in xrange(len(lattice_pos)/3):
        # find distance squared between 2 atoms
        dist = PBC_distance(lattice_pos[3*i],lattice_pos[3*i+1],lattice_pos[3*i+2],x,y,z)
        if dist < (graphRad * graphRad):
            # might be i+1
            volume_atoms.append(i)
    return volume_atoms

# calculate hashkey for a defect
def hashkey(lattice_positions,specie_list,volumeAtoms):
    # set up parameters for hashkey calculation

    lattice.pos = np.asarray(lattice_positions,dtype=np.float64)
    for i in xrange(len(specie_list)):
        if specie_list[i] == 'O_':
            specie_list[i] = 0
        elif specie_list[i] == 'Zn':
            specie_list[i] = 1
        elif specie_list[i] == 'Ag':
            specie_list[i] = 2
        # else:
        #     specie_list[i] = 3

    lattice.specie = np.asarray(specie_list,np.int32)
    lattice.cellDims = np.asarray([box_x,0,0,0,box_y,0,0,0,box_z],dtype=np.float64)
    lattice.specieList = np.asarray(['O_','Zn','Ag'],dtype=np.character)
    lattice.pos = np.around(lattice.pos,decimals = 5)

    volumeAtoms = np.asarray(volumeAtoms,dtype=np.int32)
    # Globals.PBC = [1,0,1]
    params.graphRadius = graphRad

    hashkey = Graphs.getHashKeyForAVolume(params,volumeAtoms,lattice)

    return hashkey, lattice

# find the transform matrix between 2 defect volumes
def TransformMatrix(hashkey,volumeAtoms,Lattice1,lattice_positions):
    lattice2_positions = []
    lattice2_species = []
    volume2Atoms = []

    #Lattice1 = lattice()
    Lattice2 = lattice()

    # read in lattice from file to compare
    VolumeFile = Volumes_dir + '/'+str(hashkey)+'.txt'
    if (os.path.isfile(VolumeFile)):
        input_file = open(VolumeFile, 'r')
        line = input_file.readline()
        LatLength = int(line)
        for i in xrange(LatLength):
            latline = input_file.readline()
            lat = latline.split()
            lattice2_species.append(lat[0])
            lattice2_positions.append(lat[1])
            lattice2_positions.append(lat[2])
            lattice2_positions.append(lat[3])
        # save positions to a lattice object
        Lattice2.specie = np.asarray(lattice2_species,np.int32)
        Lattice2.pos = np.asarray(lattice2_positions,dtype=np.float64)
        for j in xrange(len(volumeAtoms)):
            volline = input_file.readline()
            volume2Atoms.append(int(volline))
        volume2Atoms = np.asarray(volume2Atoms,dtype=np.int32)
        #print volume2Atoms
    else:
        sys.exit()

    Lattice1.cellDims = np.asarray([box_x,0,0,0,30,0,0,0,box_z],dtype=np.float64)
    Lattice2.cellDims = np.asarray([box_x,0,0,0,30,0,0,0,box_z],dtype=np.float64)
    volumeAtoms = np.asarray(volumeAtoms,dtype=np.int32)
    params.graphRadius = graphRad



    # -------------------------------------------------------------------------------------------
    Lattice2.charge = np.asarray([0]*LatLength,dtype=np.float64)
    Lattice1.charge = np.asarray([0]*LatLength,dtype=np.float64)
    Lattice2.specieList = np.asarray(['O_','Zn','Ag'],dtype=np.character)
    Lattice2.NAtoms = LatLength
    Lattice._writeLattice("Lattice2.dat",Lattice2,append=False, zipfile=False)
    Lattice._writeLattice("Lattice1.dat",Lattice1,append=False, zipfile=False)
    # -------------------------------------------------------------------------------------------



    print "HASHKEY:" , hashkey
    # find the transform matrix
    trnsfMatrix, lattice1, atLst1, lattice2, atLst2, cntr1, cntr2 = Graphs.prepareTheMatrix(params,volume2Atoms,Lattice2, True, volumeAtoms, Lattice1, True)

    trnsfMatrix = np.around(trnsfMatrix, decimals=5)

    #print trnsfMatrix
    del Lattice2
    del lattice2,lattice1,atLst1,atLst2,cntr1,cntr2

    # multiply transform matrix by north direction vector
    result = np.dot(trnsfMatrix,[2*x_grid_dist,0,0])
    result = np.around(result, decimals=3)
    comparex = round(x_grid_dist,3)
    comparez = round(z_grid_dist,3)

    # compare to x and z grid distances to find direction after rotation
    if comparex-0.05 < result[0] < comparex+0.05:
        if comparez-0.05 < result[2] < comparez+0.05:
            initial_direction = 3
        if -comparez-0.05 < result[2] < -comparez+0.05:
            initial_direction = 1
    elif -comparex-0.05 < result[0] < -comparex+0.05:
        if comparez-0.05 < result[2] < comparez+0.05:
            initial_direction = 4
        if -comparez-0.05 < result[2] < -comparez+0.05:
            initial_direction = 2
    elif -2*comparex-0.1 < result[0] < -2*comparex+0.1:
        initial_direction = 5
    elif 2*comparex-0.1 < result[0] < 2*comparex+0.1:
        initial_direction = 0
    else:
        print "Error in initital direction method"
        sys.exit()

    print "Initial direction", initial_direction

    return trnsfMatrix

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

        #new_list = [x+1 for x in volumeAtoms]


        # create hashkey
        final_key, Lattice1 = hashkey(lattice_positions,specie_list,volumeAtoms)
        del Lattice1

        #write_lattice(1000,full_depo_index,surface_lattice,401,0,0)
        del full_depo
        return final_key
    else:
        sys.exit()

# create the list of possible events
def create_events_list(full_depo_index,surface_lattice):
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

        #new_list = [x+1 for x in volumeAtoms]
        #print new_list

        # create hashkey for each adatom + store volume
        vol_key, Lattice1 = hashkey(lattice_positions,specie_list,volumeAtoms)
        Lattice1.NAtoms = len(specie_list)
        SaveVolume(vol_key,volumeAtoms,lattice_positions,specie_list)
        trnsfMatrix = TransformMatrix(vol_key,volumeAtoms,Lattice1,lattice_positions)
        del Lattice1

        trans_file = trans_dir + str(vol_key)+'.txt'

        # read in trans file
        if (os.path.isfile(trans_file)):
            input_file = open(trans_file, 'r')
            # skip first line
            while 1:
                latline = input_file.readline()
                barriers = latline.split()
                if len(barriers) > 1:
                    if len(barriers) == 3:
                        # multiply vector by transformation matrix
                        dis_x = float(barriers[0])*x_grid_dist
                        dis_y = float(barriers[1])*y_grid_dist
                        dis_z = float(barriers[2])*z_grid_dist
                        dir_vector = np.dot(trnsfMatrix,[dis_x,dis_y,dis_z])
                        dir_vector[0] = int(round(dir_vector[0]/x_grid_dist))
                        dir_vector[1] = int(round(dir_vector[1]/y_grid_dist))
                        dir_vector[2] = int(round(dir_vector[2]/z_grid_dist))
                        # find final hashkeys and save
                        final_key = findFinal(dir_vector,j,full_depo_index,surface_positions)
                        final_keys.append(final_key)
                        directions.append(dir_vector)
                        hashkeyExists.append(0)

                    if len(barriers) == 2:
                        # compare final hashkey to trans file
                        for k in xrange(len(final_keys)):
                            if barriers[0] == final_keys[k]:
                                rate = calc_rate(float(barriers[1]))
                                event_list.append([rate,j,directions[k]])
                                hashkeyExists[k] = 1


                else:
                    #TODO: does not seem to be working correctly!
                    for k in xrange(len(final_keys)):
                        if not hashkeyExists[k]:
                            result = singleNEB(directions[k],full_depo_index,surface_lattice,j,vol_key,final_keys[k],natoms)
                            rate = calc_rate(float(result[2]))
                            event_list.append([rate,j,directions[k]])
                    else:
                        break
        else:
            print "WARNING: Cannot find transitions file. Do searches "
            # do searches on volume and save to new trans file
            autoNEB(full_depo_index,surface_lattice,j,vol_key,natoms)
            input_file = open(trans_file, 'r')
            # skip first line
            while 1:
                latline = input_file.readline()
                barriers = latline.split()
                if len(barriers) > 1:
                    if len(barriers) == 3:
                        # multiply vector by transformation matrix
                        dis_x = float(barriers[0])*x_grid_dist
                        dis_y = float(barriers[1])*y_grid_dist
                        dis_z = float(barriers[2])*z_grid_dist
                        dir_vector = np.dot(trnsfMatrix,[dis_x,dis_y,dis_z])
                        dir_vector[0] = int(round(dir_vector[0]/x_grid_dist))
                        dir_vector[1] = int(round(dir_vector[1]/y_grid_dist))
                        dir_vector[2] = int(round(dir_vector[2]/z_grid_dist))
                        # find final hashkeys and save
                        final_key = findFinal(dir_vector,j,full_depo_index,surface_positions)
                        final_keys.append(final_key)
                        directions.append(dir_vector)

                    if len(barriers) == 2:
                        # compare final hashkey to trans file
                        for k in xrange(len(final_keys)):
                            if barriers[0] == final_keys[k]:
                                rate = calc_rate(float(barriers[1]))
                                event_list.append([rate,j,directions[k]])
                                hashkeyExists = 1
                else:
                    if hashkeyExists == None:
                        print "WARNING: Final Hashkey not found, more info needed"
                        sys.exit()
                    else:
                        break
        input_file.close()
        del volumeAtoms


    del lattice_positions
    del adatom_positions

    #print event_list
    return event_list

# run NEB to find barriers that are not known
def autoNEB(full_depo_index,surface_lattice,atom_index,hashkey,natoms):
    barrier = []
    final_keys = []
    results = []


    # create initial lattice
    write_lattice_LKMC('/initial',full_depo_index,surface_lattice,natoms)
    ini = Lattice.readLattice(NEB_dir_name_prefac+"/initial.dat")
    iniMin = copy.deepcopy(ini)

    # create cell dimensions
    cellDims = np.asarray([box_x,0,0,0,30,0,0,0,box_z],dtype=np.float64)

    # minimise lattice
    minimiser = Minimise.getMinimiser(params)
    status = minimiser.run(iniMin)

    # check max movement
    Index, maxMove, avgMove, Sep = Vectors.maxMovement(ini.pos, iniMin.pos, cellDims)
    print maxMove
    if maxMove < maxMoveCriteria:
        #ini.writeLattice("initialMin.dat")
        full_depo = copy.deepcopy(full_depo_index)
        depo_list = full_depo[atom_index]

        # just look at 6 initial directions
        dir_vector = []
        dir_vector.append([2,0,0])
        dir_vector.append([1,0,-1])
        dir_vector.append([-1,0,-1])
        dir_vector.append([1,0,1])
        dir_vector.append([-1,0,1])
        dir_vector.append([-2,0,0])

        for i in xrange(len(dir_vector)):
            # move atom
            moved_list = move_atom(depo_list, dir_vector[i] ,full_depo_index)
            if moved_list:
                full_depo[atom_index] = moved_list

                # create initial lattice and Minimise
                write_lattice_LKMC('/'+str(i),full_depo,surface_lattice,natoms)
                fin = Lattice.readLattice(NEB_dir_name_prefac+"/" + str(i) +".dat")
                finMin = copy.deepcopy(fin)

                # minimise lattice
                mini_fin = Minimise.getMinimiser(params)
                status = mini_fin.run(finMin)

                # check max movement
                Index, maxMove, avgMove, Sep = Vectors.maxMovement(fin.pos, finMin.pos, cellDims)
                if maxMove < maxMoveCriteria:
                    # run NEB on initial and final lattices
                    neb = NEB.NEB(params)
                    status = neb.run(iniMin, finMin)

                    if status:
                        print "WARNING: NEB failed to converge"
                        print "Try changing parameters in lkmcInput.IN"
                        sys.exit()

                    neb.barrier = round(neb.barrier,6)

                    # find final hashkey
                    final_key = findFinal(dir_vector[i],atom_index,full_depo_index,surface_positions)
                    results.append([dir_vector[i], final_key, neb.barrier])
                else:
                    print "WARNING: maxMove too large in final lattice:", maxMove

        print results
        write_trans_file(hashkey,results)
    else:
        print "WARNING: maxMove too large in initial lattice"

    del ini, iniMin
    del fin, finMin
    del Sep

    return

# do a single NEB and add transition to trans files
def singleNEB(direction,full_depo_index,surface_lattice,atom_index,hashkey,final_key,natoms):
    print "SINGE NEB", "="*60

    barrier = []
    final_keys = []

    # set up temp initial and final lattices.dat
    params = Input.getLKMCParams(1, "", "lkmcInput.IN")
    Input.readGlobals("lkmcInput.IN")

    print natoms
    # create initial lattice
    write_lattice_LKMC('/initial',full_depo_index,surface_lattice,natoms)
    ini = Lattice.readLattice(NEB_dir_name_prefac+"/initial.dat")
    iniMin = copy.deepcopy(ini)

    # create cell dimensions
    cellDims = np.asarray([box_x,0,0,0,30,0,0,0,box_z],dtype=np.float64)

    # minimise lattice
    minimiser = Minimise.getMinimiser(params)
    status = minimiser.run(iniMin)

    # check max movement
    Index, maxMove, avgMove, Sep = Vectors.maxMovement(ini.pos, iniMin.pos, cellDims)
    print maxMove
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
                print " difference between ini and fin is too small"
                sys.exit()

            # minimise lattice
            mini_fin = Minimise.getMinimiser(params)
            status = mini_fin.run(finMin)

            # check max movement
            Index, maxMove, avgMove, Sep = Vectors.maxMovement(fin.pos, finMin.pos, cellDims)
            if maxMove < maxMoveCriteria:
                # run NEB on initial and final lattices
                neb = NEB.NEB(params)
                status = neb.run(iniMin, finMin)

                if status:
                    print "WARNING: NEB failed to converge"
                    print "Try changing parameters in lkmcInput.IN"
                    sys.exit()

                print neb.barrier
                neb.barrier = round(neb.barrier,6)

                results = [direction, final_key, neb.barrier]
                results[0] = map(int,results[0])
                print direction
            else:
                print "WARNING: maxMove too large in final lattice"
                sys.exit()

        print results
        add_to_trans_file(hashkey,results)
    else:
        print "WARNING: maxMove too large in initial lattice"
        sys.exit()

    del ini, iniMin
    del fin, finMin
    del Sep
    print "Finished SINGlE NEB", "="*60
    return results

# after autoNEB, create new transition file with results
def write_trans_file(hashkey,results):
    trans_file = initial_dir + '/Transitions/'+str(hashkey)+'.txt'
    vectors = []
    keys =[]
    barr =[]

    if os.path.exists(trans_file):
        while 1:
            # read in existing transition file
            input_file = open(trans_file, 'r')
            line = input_file.readline()
            line = line.split()
            if len(line) > 1:
                if len(line) == 3:
                    vectors.append([int(line[0]),int(line[1]),int(line[2])])
                if len(line) == 2:
                    keys.append(line[0])
                    barr.append(float(line[1]))
            else:
                break
            input_file.close()

    # remove duplicates from list
    new_keys = []
    new_barr = []
    new_keys.append(results[0][1])
    new_barr.append(results[0][2])
    for i in xrange(len(results)):
        if results[i][1] not in new_keys:
            new_keys.append(results[i][1])
            new_barr.append(results[i][2])
    print new_keys

    outfile = open(trans_file, 'w')
    # write vectors first
    for i in xrange(len(vectors)):
        outfile.write(str(vectors[i][0])+'  '+str(vectors[i][1])+'  '+str(vectors[i][2]) + '\n')
    for i in xrange(len(results)):
        result = results[i]
        if result[0] not in vectors:
            outfile.write(str(result[0][0])+'  '+str(result[0][1])+'  '+str(result[0][2])+'\n')

    # write hashkeys + barriers
    for i in xrange(len(keys)):
        outfile.write(str(keys[i])+'    '+str(barr[i])+'\n')
    for i in xrange(len(new_keys)):
        if new_keys[i] not in keys:
            outfile.write(str(new_keys[i])+'  '+str(new_barr[i])+'\n')
    outfile.close()

    return

def add_to_trans_file(hashkey,result):
    trans_file = initial_dir + '/Transitions/'+str(hashkey)+'.txt'
    vectors = []
    keys =[]
    barr =[]

    if os.path.exists(trans_file):
        input_file = open(trans_file, 'r')
        while 1:
            # read in existing transition file
            line = input_file.readline()
            line = line.split()
            if len(line) > 1:
                if len(line) == 3:
                    vectors.append([int(line[0]),int(line[1]),int(line[2])])
                if len(line) == 2:
                    keys.append(line[0])
                    print keys
                    barr.append(float(line[1]))
            else:
                break
        input_file.close()

    outfile = open(trans_file, 'w')
    # write vectors first
    for i in xrange(len(vectors)):
        outfile.write(str(vectors[i][0])+'  '+str(vectors[i][1])+'  '+str(vectors[i][2]) + '\n')
    if result[0] not in vectors:
        outfile.write(str(result[0][0])+'  '+str(result[0][1])+'  '+str(result[0][2])+'\n')

    # write hashkeys + barriers
    for i in xrange(len(keys)):
        outfile.write(str(keys[i]) +'   '+str(barr[i])+'\n')
    if result[1] not in keys:
            outfile.write(str(result[1])+'  '+str(result[2])+'\n')
    print "Current vol", hashkey
    print result[1], keys
    outfile.close()

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
Time = 0
full_depo_list = []
surface_specie= []
surface_positions = []

print "="*80
print "~~~~~~~~ Starting predefined transitions KMC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
print "="*80
# directory setup
initial_dir = os.getcwd()
output_dir_name_prefac = initial_dir + '/Output'
Volumes_dir = initial_dir + '/Volumes'
NEB_dir_name_prefac = initial_dir + '/Temp'


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

print "~"*80

# read lattice header (in md_input dir)
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
    surface_positions.append(round(surface_lattice[i][1],7))
    surface_positions.append(round(surface_lattice[i][2],7))
    surface_positions.append(round(surface_lattice[i][3],7))

# set up temp initial and final lattices.dat
params = Input.getLKMCParams(1, "", "lkmcInput.IN")
Input.readGlobals("lkmcInput.IN")

CurrentStep = 0
if jobStatus == 'CNTIN':
    num = 0
    orig_len = len(surface_lattice)
    while 1:
        kmcFile = output_dir_name_prefac + '/KMC' + str(num) + '.dat'
        if not os.path.exists(kmcFile):
            full_depo_list = read_lattice(output_dir_name_prefac + '/KMC' + str(num-1) + '.dat',orig_len+4)
            for i in xrange(len(full_depo_list)):
                full_depo_list[i][4] = orig_len + 1 + i
            print full_depo_list
            natoms += len(full_depo_list)
            CurrentStep = len(full_depo_list)
            break
        num += 1

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

# do initial consecutive depos

#New_lattice_path = input_lattice_path
print "Current Step: 0"
while CurrentStep < (numberDepos):
    depo_list = []
    depo_list = deposition(box_x,box_z,x_grid_dist,z_grid_dist,full_depo_list,natoms)
    if depo_list:
        natoms = depo_list[4]
        full_depo_list.append(depo_list)
        write_lattice(CurrentStep,full_depo_list,surface_lattice,natoms,0,0)
        print "Writing lattice: KMC 0"
        #write_lattice_post_depo(index, natoms, depo_list[0], depo_list[1], depo_list[2], atom_species, New_lattice_path)
        CurrentStep += 1
#New_lattice_path = initial_dir + '/KMC' + str(index-1) + '.dat'
#write_lattice(0,full_depo_list,surface_lattice,natoms,0,0)
print "-" * 80
print full_depo_list
index = CurrentStep

# do KMC run
while CurrentStep < (total_steps + 1):
    # TODO: include a verbosity level
    print "Current Step: ", CurrentStep+1

    event_list = create_events_list(full_depo_list,surface_lattice)
    chosenRate, chosenEvent, chosenAtom = choose_event(event_list)

    # do deposition
    if chosenEvent[0] == 'Depo':
        while index < (CurrentStep+1):
            depo_list = []
            depo_list = deposition(box_x,box_z,x_grid_dist,z_grid_dist,full_depo_list,natoms)
            if depo_list:
                natoms = depo_list[4]
                full_depo_list.append(depo_list)
                #write_lattice(index,full_depo_list,surface_lattice,natoms,0,0)
                #write_lattice_post_depo(index, natoms, depo_list[0], depo_list[1], depo_list[2], atom_species, New_lattice_path)

                index += 1
        	#New_lattice_path = initial_dir + '/KMC' + str(index-1) + '.dat'
        	full_depo_list.append(depo_list)
        CurrentStep += 1

    # do move
    else:
        while index < (CurrentStep+1):
            moved_list =[]
            #pick random atom to move and random direction
            #move_atom_index = random.randint(0,len(full_depo_list)-1)
            move_atom_index = chosenAtom
            direction_index = chosenEvent
            moved_list = move_atom(full_depo_list[move_atom_index], direction_index, full_depo_list)
            if moved_list:
                print "pre move position: ", full_depo_list[move_atom_index]
                print "post move position: ", moved_list
                full_depo_list[move_atom_index] = moved_list
                #write_lattice(index,full_depo_list,surface_lattice,natoms,0,0)
                #write_lattice_post_move(index,natoms,full_depo_list,New_lattice_path)
                index += 1
            #New_lattice_path = initial_dir + '/KMC' + str(index-1) + '.dat'
        CurrentStep += 1
    if (CurrentStep-1)%latticeOutEvery == 0:
        latticeNo = CurrentStep/latticeOutEvery
        print "Writing lattice: KMC", latticeNo
        Time = Time + (1/chosenRate)*1E15
        Barrier = find_barrier_height(chosenRate)
        write_lattice(latticeNo,full_depo_list,surface_lattice,natoms,Time,Barrier)
    print "-" * 80

#autoNEB(full_depo_list,surface_lattice,0,5,401)

# for i in xrange(len(full_depo_list)):
#     depo_list = full_depo_list[i]
#     create_events_list(full_depo_list,surface_lattice,depo_list[1],depo_list[3],depo_list[4])



# last lines of output
print "====== Finished KMC run ========================================================"
print "Number of steps completed:	", (CurrentStep-1)
print "Number of depositions:		", len(full_depo_list)
print "="*80
print "Predefined KMC        Copyright Adam Lloyd 2016 "
print "="*80

del full_depo_list
del surface_lattice
#print full_depo_list[0][0]

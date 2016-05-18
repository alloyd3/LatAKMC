#!/usr/bin/env python

# script used to find exit rates from islands
# lattice consists of all possible basin states
# assumes hexagonal lattice (to be generalised)
# assumes 2D island
# written to include 2 types of basin states (called O_, Zn)


import os
import sys
import numpy as np
from LKMC import Lattice

# input parameters
x_grid_dist = 0.9497411251   # distance in x direction between each atom in lattice
z_grid_dist = 1.6449999809   # distance in z direction between each atom in lattice
type1_rate = 8.77E7          # rate for type1 basin state to move to type2 basin state
type2_rate = 1.56E10         # rate for type2 basin state to move to type1 basin state
type1_escapeRate = 8.02      # rate to escape basin from type1 basin state
type2_escapeRate = 1.37E7    # rate to escape basin from type2 basin state
start = 2                    # length on hexagon edge lengths for islands to consider
finish = 4
nbdist = 2                   # distance(Ang)




# Make ZnO (000-1) island (edit for your system)
def make_lattice(length,x,z):
    height = 2*length-2

    x_values = []
    z_values = []
    lattice = []

    row_len = length-1
    step = 0

    # above oxygen sites
    for i in xrange(height):
        for j in xrange(row_len):
            z_values.append(z+2*j-step)
            x_values.append(x+3*i)
            lattice.append(['O_',x+3*i,0,z+2*j-step])
        if i > length-2:
            row_len -= 1
            step -= 1
        else:
            row_len += 1
            step += 1

    row_len = length
    step = 1

    # above zn sites
    for i in xrange(height):
        for j in xrange(row_len):
            z_values.append(z+2*j-step)
            x_values.append(x+3*i+1)
            lattice.append(['Zn',x+3*i+1,0,z+2*j-step])
        if i > length-3:
            row_len -= 1
            step -= 1
        else:
            row_len += 1
            step += 1

    for i in xrange(len(lattice)):
        lattice[i][1] = lattice[i][1]*x_grid_dist
        lattice[i][3] = lattice[i][3]*z_grid_dist
        lattice[i].append(i)



    return x_values, z_values, lattice

# write out final island as a lattice
def write_lattice(lattice):
    new_lattice = 'meanLat.dat'
    outfile = open(new_lattice, 'w')
    line = str(len(lattice))
    outfile.write(line + '\n')

    # dimensions of lattice box (make larger than island)
    outfile.write(str(1000)+'  '+str(30)+'  '+str(1000)+'  ' + '\n')

    # write surface lattice first (does not change)
    for j in xrange(len(lattice)):
        #outfile.write(str(latticeLines[x]))
        outfile.write(str(lattice[j][0]) + '   ' + str(lattice[j][1])+ '    '+str(lattice[j][2])+ '   '+ str(lattice[j][3])+'  '+str(0)+'\n')

    outfile.close()
    return

# find neighbours indices (could be optimised)
def find_neighbours(lattice,n):
    nblist = []
    pos = lattice[n]
    dist2 = nbdist*nbdist

    for i in lattice:

        dist = 0
        dist += (float(i[1])-float(pos[1]))*(float(i[1])-float(pos[1]))
        dist += (float(i[2])-float(pos[2]))*(float(i[2])-float(pos[2]))
        dist += (float(i[3])-float(pos[3]))*(float(i[3])-float(pos[3]))

        if dist < dist2 and dist > 0:
            nblist.append(i[4])

    return nblist

# find T matrix (connectivity + associated rates)
def connectivity(lattice):
    #conMat = [0]*len(lattice)
    Matrix = np.zeros([len(lattice), len(lattice)])
    tao1 = np.zeros(len(lattice))
    boundaryAtoms = []

    for i in xrange(len(lattice)):
        #conMat2 = [0]*len(lattice)
        nblist = find_neighbours(lattice,i)
        for j in nblist:

            # check if boundary basin state
            if len(nblist) == 2:

                if lattice[i][0] == 'O_':
                    escape = type1_escapeRate
                    Matrix[i][j] += np.float64(type1_rate/(2*type1_rate+escape))

                elif lattice[i][0] == 'Zn':
                    escape = type2_escapeRate
                    Matrix[i][j] += np.float64(type2_rate/(2*type2_rate+escape))
                else:
                    print "Error, atom not recognised 1"
                    print lattice[i][0]
                    sys.exit()

            # assume 3 fold symmetry
            else:
                Matrix[i][j] = np.float64(1.0)/3

        # check if boundary basin state
        if len(nblist) == 2:

            if lattice[i][0] == 'O_':
                escape = type1_escapeRate
                tao1[i] = np.float64(1/(2*type1_rate+escape))
                boundaryAtoms.append(['O_',i])

            elif lattice[i][0] == 'Zn':
                escape = type2_escapeRate
                tao1[i] = np.float64(1/(2*type2_rate+escape))
                boundaryAtoms.append(['Zn',i])

            else:
                print "Error, atom not recognised 1"
                print lattice[i][0]
                sys.exit()
        else:
            if lattice[i][0] == 'O_':
                tao1[i] = np.float64(1/(3*type1_rate))

            elif lattice[i][0] == 'Zn':
                tao1[i] = np.float64(1/(3*type2_rate))

        if len(nblist) > 3:
            print "Error, more than 3 neighbours"
            sys.exit()

        if len(nblist) < 2:
            print "Error, less than 2 neighbours"
            sys.exit()

    return Matrix, tao1, boundaryAtoms



# ========================================================================
# ============= MAIN SCRIPT ==============================================
# ========================================================================


# compute mean exit rates for varying island sizes
for t in range(start,finish):

    # creat lattice
    x_values, z_values, lattice = make_lattice(t,2,300)

    # find number of atoms in island (not basin states)
    numAg = 0
    for q in range(t-1):
        numAg += t+q
    numAg = 2*numAg
    numAg += 2*t-1

    print numAg

    # find T matrix
    conMat, tao1, boundAtoms = connectivity(lattice)

    # prepare matrix to be inversed
    matrix2bInv = np.matrix(np.identity(len(lattice))) - conMat

    # attempt inverse
    try:
        occupVect = np.linalg.inv(matrix2bInv)
    except np.linalg.LinAlgError:
        print "Error: Transition Matrix not invertible."

    # assume current state is 0 (try others)
    occupVect0 = np.zeros(len(lattice))
    occupVect0[0] = 1.0
    occupVect = np.inner(occupVect, occupVect0)
    occupVect = np.squeeze(np.asarray(occupVect))

    #print occupVect
    #print "Basin states: " , len(lattice)

    # mean residence time in basin state i before leaving the basin
    tao = np.zeros(len(lattice), np.float64)
    taoSum = 0.0
    for j in range(len(lattice)):
        tao[j] = tao1[j] * occupVect[j]
        taoSum += tao[j]

    print len(lattice),taoSum

    # optional write out lattice
    write_lattice(lattice)

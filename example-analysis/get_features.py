#!/usr/bin/python

import argparse
import matplotlib.pyplot as plt
import sys
import numpy as np
from ase.io import read, write
from ase import Atoms
from ase.neighborlist import NeighborList
from scipy.spatial.distance import cdist
from tqdm import tqdm

def read_frame(filedesc):
    # 1-3
    comment = filedesc.readline()
    comment = filedesc.readline()
    comment = filedesc.readline()
    # 4
    natoms = int(filedesc.readline())
    #print(natoms)
    # 5
    comment = filedesc.readline()
    # 6,7,8
    cell = np.zeros((3,3),float)
    # line 1
    size = filedesc.readline().split()
    cell[0,0] = abs(float(size[1])-float(size[0]))
    # line 2
    size = filedesc.readline().split()
    cell[1,1] = abs(float(size[1])-float(size[0]))
    # line 1
    size = filedesc.readline().split()
    cell[2,2] = abs(float(size[1])-float(size[0]))

    #print(cell)
    # 9
    comment = filedesc.readline()
    # ITEM: ATOMS type x y z
    names = np.zeros(natoms, int)
    q = np.zeros((natoms,3),float)
    for i in range(natoms):
        line = filedesc.readline().split()
        if line[0] == 'Ti':
            names[i] = 22
        elif line[0] == 'O':
            names[i] = 8
        else:
            names[i] = 1
        q[i,0] = float(line[1])
        q[i,1] = float(line[2])
        q[i,2] = float(line[3])
    #print([natoms, cell, names, q])
    return [natoms, cell, names, q]


def main(system_now, s_start=1000, s_stop=10000, s_stride=10):

    # read the location of the slab
    sys_info_dict = {}
    sys_info_list = np.genfromtxt('Ti-min-max.dat',str)
    for sys_info in sys_info_list:
        sys_info_dict[sys_info[0]] = [float(sys_info[1]), float(sys_info[2]),float(sys_info[3])]

    # read the trajectory
    fxyz = system_now+'/out.lammpstrj'
    nwater = 128
    ixyz = open(fxyz,"r")
    frames = []
    while True:
        try:
            [ na, cell, names, pos] = read_frame(ixyz)
            if np.sum(cell) > 0:
                pbc = [True, True, True]
            else:
                pbc = [False, False, False]
            #print(names, cell, pos, pbc)
            frtemp = Atoms(numbers=names,cell=cell,positions=pos,pbc=pbc)
            #frtemp.wrap()
            print(frtemp)
            frames.append(frtemp)

        except:
            break

    # number of the hydrogen atoms has to be the same
    nhydrogen = len(np.where(frames[0].get_atomic_numbers()==1)[0])
    natom = len(frames[0].get_positions())
    #print(nhydrogen)
    nframe = len(frames)
    print("Tot frame:", nframe)

    # the index of the hydrogen atoms, has to be the same across all frames
    h_array = np.where(frames[0].get_atomic_numbers()==1)
    h_index = h_array[0]
    print("H: ", h_index)

    o_array = np.where(frames[0].get_atomic_numbers()==8)
    o_index = o_array[0]
    print("O: ", o_index)

    # last nwater*3 atoms are from water molecules
    o_tio2_index = [ i for i in o_index if i < natom-nwater*3 ]
    print("O in TiO2: ", o_tio2_index)

    # there are 64 water molecules
    o_water_index = [ i for i in o_index if i >= natom-nwater*3 ]
    print("O in H2O: ", o_water_index)

    ti_array = np.where(frames[0].get_atomic_numbers()==22)
    ti_index = ti_array[0]
    print("Ti: ", ti_index)

    slab_min = sys_info_dict[system_now][0]
    slab_max = sys_info_dict[system_now][1]
    lz = sys_info_dict[system_now][2]
    print("slab: ", slab_min, slab_max, lz)



    # O-H bond in water molecule < 1.3
    # O-H hydrogen bond < 2.4
    # H-H distance in water molecule < 1.8
    # H in first layer h2o/oh -surface Ti < 3.6
    cutoff = 4.5

    selected_frames = frames[int(s_start):int(s_stop):int(s_stride)]

    h_dis_all = np.zeros((len(selected_frames), nhydrogen, 13))
    h_env_all = np.zeros((len(selected_frames), nhydrogen, 7, 4))

    # build neighborlist
    r_cut_list =  np.ones(natom)*cutoff/2.
    r_cut_list[ti_array] = cutoff
    nl = NeighborList(r_cut_list, skin=0., sorted=False, self_interaction=False,
                     bothways=True)

    for num_frame,frame in tqdm(enumerate(selected_frames)):
        nl.update(frame)

        # the distance vector array: 
        # central_atom, rHTi, rHH, rHO_TiO2, rOTi, rOwTi, v, mu, rOO, r_O_1_z
        # $0: index_of_H;   
        # $1: H-Ti;  $2: H-H; $3: rHO_TiO2; $4:rOTi; $5:rOwTi; $6: r_O_1_z
        # $7: v; $8: mu; $9: O-O; 
        h_dis = np.zeros((nhydrogen, 13))
        for h_i,central_atom in enumerate(h_index):

            # this is to identify whether the H is closer to top or bottom slab
            H_z_now = frame.get_positions()[central_atom,2]
            if H_z_now < 0: H_z_now += lz
            if H_z_now > lz: H_z_now -= lz
            h_top = H_z_now-slab_max
            h_bottom = lz+slab_min-H_z_now
            if abs(h_top) > lz/2:
                h_top -= lz*(h_top)/abs(h_top)
            if abs(h_bottom) > lz/2:
                h_bottom -= lz*(h_bottom)/abs(h_bottom)
            if abs(h_top) < abs(h_bottom): # closer to the top slab
                topslab = True
                dis_2_slab = np.amin([np.abs(h_top),10])
            else: # closer to the bottom slab
                topslab = False
                dis_2_slab = np.amin([np.abs(h_bottom),10])
            #if dis_2_slab<0: print(dis_2_slab, lz, 'hz:', H_z_now, "h_top:", h_top, "h_bottom:", h_bottom)


            indices, offsets = nl.get_neighbors(central_atom)

            # compute displacements r_ij
            displacements = np.zeros((len(indices),5))
            j=0
            for i, offset in zip(indices, offsets):
                displacements[j,0] = i
                rij = frame.positions[i] + np.dot(offset, frame.get_cell()) - frame.positions[central_atom]
                displacements[j,1:4] = rij
                displacements[j,4] = np.linalg.norm(rij) # scalar distance
                j+=1

            # build sorted lists
            rO_list = np.array([d for d in displacements if d[0] in o_index])
            rO_list = rO_list[rO_list[:, 4].argsort()]
            rH_list = np.array([d for d in displacements if d[0] in h_index])
            rH_list = rH_list[rH_list[:, 4].argsort()]
            rTi_list = np.array([d for d in displacements if d[0] in ti_index])
            if len(rTi_list) > 1:
                r_Ti = rTi_list[rTi_list[:, 4].argsort()][0]
            else:
                r_Ti = rTi_list
            #print(r_Ti)

            r_vec = []
            # collect the displacement of neiboring atoms
            r_O_1, r_O_2 = rO_list[0], rO_list[1]
            r_vec.append([8, r_O_1[1], r_O_1[2], r_O_1[3]])
            r_vec.append([8, r_O_2[1], r_O_2[2], r_O_2[3]])
            # find the Hs that are bonded to the Os
            nh = 0
            for hh in rH_list:
                if np.linalg.norm(hh[1:4]-r_O_1[1:4]) < 1.5 or np.linalg.norm(hh[1:4]-r_O_2[1:4]) < 1.5:
                    nh +=1
                    if nh <= 4:
                        r_vec.append([1, hh[1], hh[2], hh[3]])
            # Ti
            if len(rTi_list) > 1:
                r_Ti = rTi_list[rTi_list[:, 4].argsort()][0]
                r_vec.append([22, r_Ti[1], r_Ti[2], r_Ti[3]])

            r_vec = np.reshape(np.asarray(r_vec),(-1,4))
            if not topslab:
                r_vec[:, 3] *= -1        
            # store
            h_env_all[num_frame,h_i,:len(r_vec),:] = r_vec

            # compute the (scalar) |rHTi| between H and closest Ti
            try: 
                rHTi = np.amin([d[4] for d in displacements if int(d[0]) in ti_index])
                #print(rHTi)
            except:
                rHTi = 10    

            # compute the (scalar) |rHH| between H and closest H
            rHH = rH_list[0,4]
            # compute the (scalar) |rHH| between H and second closest H
            rHH2 = rH_list[1,4]
            # take the z component of the vector H->H
            if topslab: # closer to the top slab
                rHH_z = -rH_list[0,3]/rH_list[0,4]
            else: # closer to the bottom slab
                rHH_z = rH_list[0,3]/rH_list[0,4]

            # compute the (scalar) |rHOti| between H and closest O in TiO2
            try: 
                rHO_TiO2 = np.amin([d[4] for d in displacements if int(d[0]) in o_tio2_index])
                #print(rHTi)
            except:
                rHO_TiO2 = 10            

            # compute the (scalar) |rOti| between the closest O and Ti
            try:
                rOTi = np.linalg.norm((r_O[1:4]-r_Ti[1:4]))
            except:
                rOTi = 10

            # compute the (scalar) |rOti| between the closest Ow and Ti
            try:
                rOw_list = np.array([d for d in displacements if d[0] in o_water_index]) ###
                r_O = rOw_list[rOw_list[:, 4].argsort()][0]
                rOwTi_old = np.linalg.norm(r_O[1:4]-r_Ti[1:4])
                rOwTi = np.amin(np.array([ np.linalg.norm(r_O[1:4]-r_Ti_now[1:4]) for r_Ti_now in rTi_list ]))
                if (rOwTi_old - rOwTi)**2.>1: print(rOwTi_old,rOwTi)
            except:
                rOwTi = 10

            # compute the (scalar) |rOwOt| between closest O in water and closest O in TiO2
            try:
                rOt_list = np.array([d for d in displacements if d[0] in o_tio2_index])
                r_O = rOw_list[rOw_list[:, 4].argsort()][0]
                rOwOt = np.amin(np.array([ np.linalg.norm(r_O[1:4]-r_Ot_now[1:4]) for r_Ot_now in rOt_list ]))
            except:
                rOwOt = 10


            # compute (vector) displacements between H and closest 2 oxygen atoms (r_O_1, r_O_2)
            r_O_1, r_O_2 = rO_list[0], rO_list[1]
            rHO1, rHO2 = r_O_1[4], r_O_2[4]
            rOO = np.linalg.norm((r_O_1[1:4]-r_O_2[1:4]))
            #proton-transfer coordinate ν = d(D-H) − d(A-H), 
            v = rHO1-rHO2
            #the symmetric stretch coordinate μ = d(D-H) + d(A-H)
            mu = rHO1+rHO2
            #print(rHO1, rHO2, rOO)

            # take the z component of the vector O->H
            if topslab: # closer to the top slab
                r_O_1_z = -r_O_1[3]
            else: # closer to the bottom slab
                r_O_1_z = r_O_1[3]


            h_dis[h_i] = [ central_atom, rHTi, rHH, rHO_TiO2, rOTi, rOwTi, r_O_1_z, rHH_z, v, mu, rOO, rHH2, rOwOt]
        h_dis_all[num_frame,:,:] = h_dis

    with open(system_now+'-h-dis-env.npy', 'wb') as f:
        np.save(f, h_dis_all)
        np.save(f, h_env_all)

if __name__ == '__main__':
    main(*sys.argv[1:])

#!usr/bin/env python
import random
import numpy as np
def metropolis(lattice,kT,J,h):
    "metropolis algorithm for updating lattice"
    N=len(lattice)
    for i in range(N):
        for j in range(N):
            neighbours= lattice[i,(j+1)%N] + lattice[i,(j-1)%N] + lattice[(i+1)%N,j] + lattice[(i-1)%N,j]
            delta_H = (-2*lattice[i,j])*(J*neighbours + h)
            if delta_H>=0:
                lattice[i,j]=-lattice[i,j]
            elif random.random()<np.exp((delta_H)/kT):
                lattice[i,j]=-lattice[i,j]
return lattice

def convert_to_pbm(lattice):

    for i in range(len(lattice)):
        for j in range(len(lattice)):
            if lattice[i,j]==-1:
                lattice[i,j]=0



#!usr/bin/env python


"importing relevant packages"

import random
import numpy as np




"metropolis algorithms for lattices of different dimensions"

def metropolis(lattice,kT,J,h):
    "Performs one sweep of the metropolis algorithm for the 2D lattice"
    N=len(lattice)
    for i in range(N):
        for j in range(N):
            neighbours= lattice[i,(j+1)%N] + lattice[i,(j-1)%N] + lattice[(i+1)%N,j] + lattice[(i-1)%N,j] # modular gives periodicity
            delta_H = (-2*lattice[i,j])*(J*neighbours + h)
            if delta_H>=0:
                lattice[i,j]=-lattice[i,j]
            elif random.random()<np.exp((delta_H)/kT):
                lattice[i,j]=-lattice[i,j]
return lattice


def metropolis_3D(lattice,kT,J,h):
    "Performs one sweep of the metropolis algorithm for the 3D lattice"
    N=len(lattice)
    for i in range(N):
        for j in range(N):
            for k in range(N):
                neighbours= lattice[i,(j+1)%N,k] + lattice[i,(j-1)%N,k] + lattice[(i+1)%N,j,k] + lattice[(i-1)%N,j,k]+ lattice[i,j,(k+1)%N] + lattice[i,j,(k-1)%N]
                delta_H=(-2*lattice[i,j,k])*(J*neighbours + h)
                if delta_H>=0:
                    lattice[i,j,k]=-lattice[i,j,k]
                elif random.random()<np.exp((delta_H)/kT):
                    lattice[i,j,k]=-lattice[i,j,k]
    return lattice


def metropolis_1D(lattice,kT,J,h):
    "Performs one sweep of the metropolis algorithm for the 1D lattice"
    N=len(lattice)
    for i in range(N):
                neighbours=lattice[(i+1)%N]+lattice[(i-1)%N]
                delta_H=(-2*lattice[i])*(J*neighbours + h)
                if delta_H>=0:
                    lattice[i]=-lattice[i]
                elif random.random()<np.exp((delta_H)/kT):
                    lattice[i]=-lattice[i]
    return lattice




"Functions that find energy per lattice point of lattices of different dimensions"

def find_energy(lattice,J,h):
    "finding the energy per lattice point for 2D lattice"
    N=len(lattice)
    energy=0.0
    for i in range(N):
        for j in range(N):
            s=lattice[i,j]
            neighbours= lattice[i,(j+1)%N] + lattice[i,(j-1)%N] + lattice[(i+1)%N,j] + lattice[(i-1)%N,j]
            energy=energy-(J*neighbours*s)-4*h*s
    return float(energy/(2*N*N)) #divide by 2 to overcome double counting


def find_energy_3D(lattice,J,h):
    "finding the energy per lattice point for 3D lattice"
    N=len(lattice)
    energy=0.0
    for i in range(N):
        for j in range(N):
            for k in range(N):
                s=lattice[i,j,k]
                neighbours= lattice[i,(j+1)%N,k] + lattice[i,(j-1)%N,k] + lattice[(i+1)%N,j,k] + lattice[(i-1)%N,j,k]+ lattice[i,j,(k+1)%N] + lattice[i,j,(k-1)%N]
                energy=energy-(J*neighbours*s)-4*h*s
    return float((energy/(2*N*N*N))) #divide by 2 to overcome double counting

   
def find_energy_1D(lattice,J,h):
    "finding the energy per lattice point for the 1D lattice
    N=len(lattice)
    energy=0.0
    for i in range(N):
        s=lattice[i]
        neighbours=lattice[(i+1)%N]+lattice[(i-1)%N]
        energy=energy-(J*neighbours*s)-h*s
    return float((energy/(2*N))) #divide by 2 to overcome double counting




"Functions that find the magnetization"

def find_magnetization(lattice):
    "finding magnetization for 2D lattice"
    return abs(float(np.sum(lattice))/(len(lattice)*len(lattice)))


def find_magnetization_3D(lattice):
    "finding magnetization for 3D lattice"
    return abs(float(np.sum(lattice))/(len(lattice)*len(lattice)*len(lattice)))


def find_magnetization_1D(lattice):
    "finding magnetization for 1D lattice"
    return abs(float(np.sum(lattice))/(len(lattice)))



"Miscellaneous functions"

def convert_to_pbm(lattice):
    "Puts matrix into pbm format, spin down=white cell, spin up=black cell"
    for i in range(len(lattice)):
        for j in range(len(lattice)):
            if lattice[i,j]==-1:
                lattice[i,j]=0

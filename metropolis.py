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


def metropolis_3D(lattice,kT,J,h):
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
    N=len(lattice)
    for i in range(N):
                neighbours=lattice[(i+1)%N]+lattice[(i-1)%N]
                delta_H=(-2*lattice[i])*(J*neighbours + h)
                if delta_H>=0:
                    lattice[i]=-lattice[i]
                elif random.random()<np.exp((delta_H)/kT):
                    lattice[i]=-lattice[i]
    return lattice


def convert_to_pbm(lattice):

    for i in range(len(lattice)):
        for j in range(len(lattice)):
            if lattice[i,j]==-1:
                lattice[i,j]=0


def find_energy(lattice,J,h):
    "finding the energy per lattice point"
    N=len(lattice)
    energy=0.0
    for i in range(N):
        for j in range(N):
            s=lattice[i,j]
            neighbours= lattice[i,(j+1)%N] + lattice[i,(j-1)%N] + lattice[(i+1)%N,j] + lattice[(i-1)%N,j]
            energy=energy-(J*neighbours*s)-4*h*s
    return float(energy/(2*N*N))


def find_energy_3D(lattice,J,h):
    N=len(lattice)
    energy=0.0
    for i in range(N):
        for j in range(N):
            for k in range(N):
                s=lattice[i,j,k]
                neighbours= lattice[i,(j+1)%N,k] + lattice[i,(j-1)%N,k] + lattice[(i+1)%N,j,k] + lattice[(i-1)%N,j,k]+ lattice[i,j,(k+1)%N] + lattice[i,j,(k-1)%N]
                energy=energy-(J*neighbours*s)-4*h*s
    return float((energy/(2*N*N*N)))

    
def find_energy_1D(lattice,J,h):
    N=len(lattice)
    energy=0.0
    for i in range(N):
        s=lattice[i]
        neighbours=lattice[(i+1)%N]+lattice[(i-1)%N]
        energy=energy-(J*neighbours*s)-h*s
    return float((energy/(2*N)))



def find_magnetization(lattice):
    "finding average magnetization per site"
    return abs(float(np.sum(lattice))/(len(lattice)*len(lattice)))

def find_magnetization_3D(lattice):
    "finding average magnetization per site"
    return abs(float(np.sum(lattice))/(len(lattice)*len(lattice)*len(lattice)))



def find_magnetization_1D(lattice):
    "finding average magnetization per site"
    return abs(float(np.sum(lattice))/(len(lattice)))






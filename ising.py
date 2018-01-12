#!usr/bin/env python



"importing packages"

import random
import functions
import matplotlib.pyplot as plt
import numpy as np
import mypbm
import plotting



"defining relevant parameters"

N=600
kT=0.5
J=1.
h=0.0
steps=2000
steps_to_equilibrium=2000
steps_to_calculate=1000
initial_state = 2*np.random.randint(2, size=(N,N))-1
initial_state_3D=2*np.random.randint(2, size=(N,N,N))-1
initial_state_1D = 2*np.random.randint(2, size=(N))-1
mag=np.zeros(steps)
step_array=np.arange(0,steps,1)
lattice_1D=initial_state_1D
lattice=initial_state
lattice_3D=initial_state_3D




"extra functions section"

def plot_lattice():
	for i in range(steps):
		functions.metropolis(lattice,kT,J,h)

	functions.convert_to_pbm(lattice)
	mypbm.myplot(lattice,"plot.pbm")


def find_number_of_steps_to_equilibrium():
	lattice=initial_state
	for i in range(steps):
		functions.metropolis(lattice,kT,J,h)
		mag[i]=functions.find_magnetization(lattice)
	fig=plt.figure()
	plt.plot(step_array,mag)
	fig.suptitle('Average magnetization per site versus number of sweeps',fontsize=15)
	plt.xlabel('Number of sweeps', fontsize=13)
	plt.ylabel('Average magnetization per site [$\\mu$]', fontsize=13)
	plt.show()



"Here is where the functions required were called, below is an example"
plot_lattice()


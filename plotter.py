#!usr/bin/env python


"importing relevant packages"

import random
import functions
import matplotlib.pyplot as plt
import numpy as np




"defining relevant parameters"

N=10					#define your lattice size
J=30.0					#interaction term
h=0.0					#external magnetic field
steps=1000
steps_to_equilibrium=1500
steps_to_calculate=1000
initial_state = 2*np.random.randint(2, size=(N,N))-1
initial_state_3D=2*np.random.randint(2, size=(N,N,N))-1			#initial random configurations of different dimension lattices
initial_state_1D = 2*np.random.randint(2, size=(N))-1




"The following three functions plot the physical quantities in the case of 1D, 2D and 3D"

def plot_physical_quantities_2D():
	intervals=200
	T=np.linspace(0.8,4.0,intervals)
	energy=np.zeros(intervals)
	magnetization=np.zeros(intervals)
	susceptibility=np.zeros(intervals)
	specific_heat=np.zeros(intervals)
	for j in range(len(T)):
		t=T[j]
		lattice=initial_state
		E1=0.0
		B1=0.0
		B2=0.0
		E2=0.0
		for i in range(steps_to_equilibrium):
			functions.metropolis(lattice,t,J,h)

		for i in range(steps_to_calculate):
			functions.metropolis(lattice,t,J,h)
			E=functions.find_energy(lattice,J,h)
			B=functions.find_magnetization(lattice)
			E1=E1+E
			B1=B1+B
			B2=B2+(B*B)
			E2=E2+(E*E)
		E1_av=E1/steps_to_calculate
		E2_av=E2/steps_to_calculate
		B1_av=B1/steps_to_calculate
		B2_av=B2/steps_to_calculate

		energy[j]=E1_av
		magnetization[j]=B1_av
		susceptibility[j]=float((B2_av-(B1_av)**2))/t
		specific_heat[j]=float((E2_av-(E1_av)**2))/(t**2)
	
	def plot_energy():
		fig=plt.figure()
		plt.plot(T,energy,'ro')
		fig.suptitle('Average energy per site versus temperature\n2D lattice', fontsize=15)
		plt.xlabel('Temperature [J/k]',fontsize=13)
		plt.ylabel('Average energy per site [J]',fontsize=13)
		plt.show()
	def plot_mag():
		fig=plt.figure()
		plt.plot(T,magnetization,'ro')
		fig.suptitle('Average magnetization per site versus temperature\n2D lattice',fontsize=15)
		plt.xlabel('Temperature [J/k]',fontsize=13)
		plt.ylabel('Absolute average magnetization per site [$\\mu$]',fontsize=13)
		plt.show()
	def plot_sus():
		fig=plt.figure()
		plt.plot(T,susceptibility,'ro')
		fig.suptitle('Average magnetic susceptibility per site versus temperature\n2D lattice',fontsize=15)
		plt.xlabel('Temperature [J/k]', fontsize=13)
		plt.ylabel('Average magnetic susceptibility per site [$\\mu$/k]', fontsize=13)
		plt.show()
	def plot_sc():
		fig=plt.figure()
		plt.plot(T,specific_heat,'ro')
		fig.suptitle('Average specific heat per site versus temperature\n2D lattice',fontsize=15)
		plt.xlabel('Temperature [J/k]',fontsize=13)
		plt.ylabel('Average specific heat per site [J/k^2]', fontsize=13)
		plt.show()

	plot_energy()
	plot_mag()
	plot_sus()
	plot_sc()


def plot_physical_quantities_3D():
	intervals=600
	T2=np.linspace(2.3,6.3,intervals)
	energy_3D=np.zeros(intervals)
	magnetization_3D=np.zeros(intervals)
	susceptibility_3D=np.zeros(intervals)
	specific_heat_3D=np.zeros(intervals)

	for j in range(len(T2)):
		t=T2[j]
		lattice=initial_state_3D
		E1=0.0
		B1=0.0
		B2=0.0
		E2=0.0
		for i in range(steps_to_equilibrium):
			functions.metropolis_3D(lattice,t,J,h)

		for i in range(steps_to_calculate):
			functions.metropolis_3D(lattice,t,J,h)
			E=functions.find_energy_3D(lattice,J,h)
			B=functions.find_magnetization_3D(lattice)
			E1=E1+E
			B1=B1+B
			B2=B2+(B*B)
			E2=E2+(E*E)
		E1_av=E1/steps_to_calculate
		E2_av=E2/steps_to_calculate
		B1_av=B1/steps_to_calculate
		B2_av=B2/steps_to_calculate

		energy_3D[j]=E1_av
		magnetization_3D[j]=B1_av
		susceptibility_3D[j]=float((B2_av-(B1_av)**2))/t
		specific_heat_3D[j]=float((E2_av-(E1_av)**2))/(t**2)
	
	def plot_energy():
		fig=plt.figure()
		plt.plot(T2,energy_3D,'ro')
		fig.suptitle('Average energy per site versus temperature\n3D lattice', fontsize=15)
		plt.xlabel('Temperature [J/k]',fontsize=13)
		plt.ylabel('Average energy per site [J]',fontsize=13)
		plt.show()
	def plot_mag():
		fig=plt.figure()
		plt.plot(T2,magnetization_3D,'ro')
		fig.suptitle('Average magnetization per site versus temperature\n3D lattice',fontsize=15)
		plt.xlabel('Temperature [J/k]',fontsize=13)
		plt.ylabel('Absolute average magnetization per site [$\\mu$]',fontsize=13)
		plt.show()
	def plot_sus():
		fig=plt.figure()
		plt.plot(T2,susceptibility_3D,'ro')
		fig.suptitle('Average magnetic susceptibility per site versus temperature\n3D lattice',fontsize=15)
		plt.xlabel('Temperature [J/k]', fontsize=13)
		plt.ylabel('Average magnetic susceptibility per site [$\\mu$/k]', fontsize=13)
		plt.show()
	def plot_sc():
		fig=plt.figure()
		plt.plot(T2,specific_heat_3D,'ro')
		fig.suptitle('Average specific heat per site versus temperature\n3D lattice',fontsize=15)
		plt.xlabel('Temperature [J/k]',fontsize=13)
		plt.ylabel('Average specific heat per site [J/k^2]', fontsize=13)
		plt.show()

	plot_energy()
	plot_mag()
	plot_sus()
	plot_sc()


def plot_physical_quantities_1D():
	intervals=600
	T1=np.linspace(0.1,7.0,intervals)
	energy_1D=np.zeros(intervals)
	magnetization_1D=np.zeros(intervals)
	susceptibility_1D=np.zeros(intervals)
	specific_heat_1D=np.zeros(intervals)
	for j in range(len(T1)):
		t=T1[j]
		lattice=initial_state_1D
		E1=0.0
		B1=0.0
		B2=0.0
		E2=0.0
		for i in range(steps_to_equilibrium):
			functions.metropolis_1D(lattice,t,J,h)

		for i in range(steps_to_calculate):
			functions.metropolis_1D(lattice,t,J,h)
			E=functions.find_energy_1D(lattice,J,h)
			B=functions.find_magnetization_1D(lattice)
			E1=E1+E
			B1=B1+B
			B2=B2+(B*B)
			E2=E2+(E*E)
		E1_av=E1/steps_to_calculate
		E2_av=E2/steps_to_calculate
		B1_av=B1/steps_to_calculate
		B2_av=B2/steps_to_calculate

		energy_1D[j]=E1_av
		magnetization_1D[j]=B1_av
		susceptibility_1D[j]=float((B2_av-(B1_av)**2))/t
		specific_heat_1D[j]=float((E2_av-(E1_av)**2))/(t**2)
    	def plot_energy():
		fig=plt.figure()
		plt.plot(T1,energy_1D,'ro')
		fig.suptitle('Average energy per site versus temperature\n1D lattice', fontsize=15)
		plt.xlabel('Temperature [J/k]',fontsize=13)
		plt.ylabel('Average energy per site [J]',fontsize=13)
		plt.show()
	def plot_mag():
		fig=plt.figure()
		plt.plot(T1,magnetization_1D,'ro')
		fig.suptitle('Average magnetization per site versus temperature\n1D lattice',fontsize=15)
		plt.xlabel('Temperature [J/k]',fontsize=13)
		plt.ylabel('Absolute average magnetization per site [$\\mu$]',fontsize=13)
		plt.show()
	def plot_sus():
		fig=plt.figure()
		plt.plot(T1,susceptibility_1D,'ro')
		fig.suptitle('Average magnetic susceptibility per site versus temperature\n1D lattice',fontsize=15)
		plt.xlabel('Temperature [J/k]', fontsize=13)
		plt.ylabel('Average magnetic susceptibility per site [$\\mu$/k]', fontsize=13)
		plt.show()
	def plot_sc():
		fig=plt.figure()
		plt.plot(T1,specific_heat_1D,'ro')
		fig.suptitle('Average specific heat per site versus temperature\n1D lattice',fontsize=15)
		plt.xlabel('Temperature [J/k]',fontsize=13)
		plt.ylabel('Average specific heat per site [J/k^2]', fontsize=13)
		plt.show()

	plot_energy()
	plot_mag()
	plot_sus()
	plot_sc()




"This function plots the recorded critical temperatures as a function of lattice size to extrapolate the critical temperature of an infinite lattice"

def find_critical_temperature():
	x=[0.0001,0.0004,0.0025,0.01]
	y=[2.27,2.28,2.33,2.52]
	fig=plt.figure()
	plt.plot(x,y,'ko')
	fig.suptitle('Critical temperature for different sized lattices',fontsize=15)
	plt.xlabel('1/N^2',fontsize=13)
	plt.ylabel('Critical temperature [J/k]', fontsize=13)
	p=np.polyfit(x,y,1)
	print p[0]
	print p[1]
	y_model=[]
	for i in x:
		y_model.append(i*p[0]+p[1])
	plt.plot(x,y_model)
	plt.show()




"These functions find the relationship between the interaction term J and the ground energy of the system, in the 1, 2 and 3 dimensional cases"

def find_ground_energy_2D():
	energy_model=[]
	intervals=10
	J1=np.linspace(1.0,10.0,intervals)
	energy=np.zeros(intervals)
	for j in range(len(J1)):
		j_part=J1[j]
		T=0.001
		lattice=initial_state
		E1=0.0
		for i in range(steps_to_equilibrium):
			functions.metropolis(lattice,T,j_part,h)

		for i in range(steps_to_calculate):
			functions.metropolis(lattice,T,j_part,h)
			E=functions.find_energy(lattice,j_part,h)
			E1=E1+E
		E1_av=E1/steps_to_calculate
		energy[j]=E1_av
	
	fig=plt.figure()
	plt.plot(J1,energy,'ro')
	fig.suptitle('2D ground state energy as function of J',fontsize=15)
	plt.xlabel('J-interaction term',fontsize=13)
	plt.ylabel('ground state energy [J]', fontsize=13)
	p=np.polyfit(J1,energy,1)
	for i in J1:
		energy_model.append(i*p[0]+p[1])
	plt.plot(J1,energy_model)
	print p[0]
	plt.show()

def find_ground_energy_3D():
	energy_model=[]
	intervals=10
	J1=np.linspace(1.0,10.0,intervals)
	energy=np.zeros(intervals)
	for j in range(len(J1)):
		j_part=J1[j]
		T=0.001
		lattice=initial_state_3D
		E1=0.0
		for i in range(steps_to_equilibrium):
			functions.metropolis_3D(lattice,T,j_part,h)

		for i in range(steps_to_calculate):
			functions.metropolis_3D(lattice,T,j_part,h)
			E=functions.find_energy_3D(lattice,j_part,h)
			E1=E1+E
		E1_av=E1/steps_to_calculate
		energy[j]=E1_av
	
	fig=plt.figure()
	plt.plot(J1,energy,'ro')
	fig.suptitle('3D ground state energy as function of J',fontsize=15)
	plt.xlabel('J-interaction term',fontsize=13)
	plt.ylabel('ground state energy [J]', fontsize=13)
	p=np.polyfit(J1,energy,1)
	for i in J1:
		energy_model.append(i*p[0]+p[1])
	plt.plot(J1,energy_model)
	print p[0]
	plt.show()


def find_ground_energy_1D():
	energy_model=[]
	intervals=10
	J1=np.linspace(1.0,10.0,intervals)
	energy=np.zeros(intervals)
	for j in range(len(J1)):
		j_part=J1[j]
		T=0.001
		lattice=initial_state_1D
		E1=0.0
		for i in range(steps_to_equilibrium):
			functions.metropolis_1D(lattice,T,j_part,h)

		for i in range(steps_to_calculate):
			functions.metropolis_1D(lattice,T,j_part,h)
			E=functions.find_energy_1D(lattice,j_part,h)
			E1=E1+E
		E1_av=E1/steps_to_calculate
		energy[j]=E1_av
	
	fig=plt.figure()
	plt.plot(J1,energy,'ro')
	fig.suptitle('1D ground state energy as function of J',fontsize=15)
	plt.xlabel('J-interaction term',fontsize=13)
	plt.ylabel('ground state energy [J]', fontsize=13)
	p=np.polyfit(J1,energy,1)
	for i in J1:
		energy_model.append(i*p[0]+p[1])
	plt.plot(J1,energy_model)
	print p[0]
	plt.show()




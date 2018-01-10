N=500
J=1.0
h=-1.0
steps=1000
steps_to_equilibrium=1000
steps_to_calculate=500
initial_state = 2*np.random.randint(2, size=(N,N))-1

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

#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt

##################################################################################################
# Physical Constants

M = 39.948
angstrom = 10e-10
avogadro = 6.022e-23
kb = 1.38064852e-23
femto = 10e-15
uma =  1.660538921e-27 
J = uma*angstrom**2/femto**2
kb_rescaled = kb/J

##################################################################################################
# System Parameters
N = 864
m =10
T = 94.4
L = 34.1580
nsteps = 10000
V = L**3
rho = N/V
sigma = 3.4

radii = np.linspace(0,20,100, endpoint=True)

##################################################################################################
# Simulation File names

file_dists = './data_analysis/distances.dat'
file_temps = './data_analysis/temps.dat'
file_energies = './data_analysis/energies.dat'
file_vels = './data_analysis/sq_velocities.dat'

##################################################################################################
# Figure names

rdf_fig = './documents/rdf.jpg'
temps_fig = './documents/T_vs_t.jpg'
energy_fig = './documents/E_vs_t.jpg'
vels_fig = './documents/vels_dist.jpg'

##################################################################################################
# File Opening

distances =  np.genfromtxt(file_dists, delimiter=" ", usecols=range(0,1), unpack=False)
temp_series =  np.genfromtxt(file_temps, delimiter=" ", usecols=range(0,1), unpack=False)
ener_series =  np.genfromtxt(file_energies, delimiter=" ", usecols=range(0,1), unpack=False)
vels_series =  np.sqrt(np.genfromtxt(file_vels, delimiter=" ", usecols=range(0,1), unpack=False))

##################################################################################################
# Radial Distribution Function Calculation

distances = distances[distances != 0]
histogram_distances = np.histogram(distances,radii, density=False)[0]

dr = radii[1]
rdf = histogram_distances
rdf = rdf/(N*rho)

for i, r in enumerate(radii[1:]):
    dr = radii[i+1] - radii[i]
    rdf[i] = rdf[i]/(4 * np.pi * r**2 * dr)

max_num = 10
local_max = []
max_rdf_idxs = np.argsort(rdf)[-max_num:]

for idx in max_rdf_idxs:
    if rdf[idx - 1] < rdf[idx] and rdf[idx] > rdf[idx + 1]: local_max.append(idx)
local_max

##################################################################################################
# Theoretical Maxwell Bolzmann Velocity Distribution

x = np.linspace(0,0.007,1000)
y = (4*np.pi*(np.sqrt(M/(2*np.pi*kb_rescaled*T))**3) * (x**2) *np.exp(-x**2*M/(2*kb_rescaled*T)))


##################################################################################################
# Velocity Distribution Plotting

plt.title('Maxwell Boltzmann Distribution T={T} K'.format(T=T))
n3, bins3, patches3 = plt.hist(vels_series, 200, density=True, facecolor='g', alpha=0.75, \
                              label='Generated Histogram')
plt.plot(x,y, c='r', label='Theoretical Distribution for T={T} K'.format(T=T))
plt.grid(True)
plt.legend()
plt.xlabel(r'Velocity [$\AA / fs$]')
plt.ylabel('F(v)')
plt.savefig(vels_fig)
plt.show()

##################################################################################################
# Energy vs time Plotting

plt.title('Total Energy vs Time')
plt.scatter(np.arange(len(ener_series))*100, ener_series*0.004184, label='Total Energy [MJ/mol]', c='r', s=4)
plt.plot(np.arange(len(ener_series))*100, ener_series*0.004184, c='r', lw=0.5 )
plt.grid(True)
plt.xlabel('Time [fs]')
plt.ylabel('Energy [MJ/Mol]')
plt.legend()
plt.grid(True)
plt.ylim(-3.8,-3.2)
plt.savefig(energy_fig)
plt.show()

##################################################################################################
# RDF Plotting

plt.title('Radial Distribution Functions')
plt.plot(radii[1:]/sigma, rdf, label='g(r)')
for i in local_max:
    x=radii[1:][i]/sigma
    y=rdf[i]
    plt.scatter(x, y)
    plt.annotate(r'{x} $\AA$'.format(x=round(x*3.4,1)), (x, y))
plt.grid(True)
plt.legend()
plt.xlabel(r'$r/\sigma$')
plt.ylabel('g(r)')
plt.xlim(0,4)
plt.savefig(rdf_fig)
plt.show()


##################################################################################################
# Temperature vs time plotting

start=5000
l=400
plt.title('Temperature vs Time')
plt.scatter(np.arange(start,start+l)*5 , temp_series[start:start+l], label='Temperature [k]', c='b', s=4)
plt.plot(np.arange(start,start+l)*5 , temp_series[start:start+l], c='b', lw=0.5)
plt.xlabel('Time [fs]')
plt.ylabel('Temperature [K]')
plt.ylim(85,105)
plt.legend()
plt.grid(True)
plt.savefig(temps_fig)
plt.show()
plt.close()

##################################################################################################







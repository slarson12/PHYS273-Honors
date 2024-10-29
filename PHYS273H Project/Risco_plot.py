# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 13:12:20 2024

@author: sierr
"""
import matplotlib.pyplot as plt
import numpy as np

# Constants 
G = 6.67 * 10**-11 # Nm^2 / kg^2
c = 2.99 * 10**8 #m/s
Ms = 2 * 10**30 #kg
M = 4.297 * 10**6 * Ms # mass of Sag A* in solar masses (Ms) 
R_g = G*M/(c**2)
print("gravitational radius of Sag A*:", R_g)

# plot the ISCO, MB orbit, and Photon circular orbit as functions of the spin parameter 

#define spin values from -1 to 1
a_values = np.linspace(-1, 1, 100) 

# Innermost Stable Circular Orbit as a function of a
def ISCO_func(a):
    Z1 =(1+(1-a**2)**(1/3) * ((1+a)**(1/3) + (1-a)**(1/3)))
    Z2 = ((3*a**2 + (Z1**2))**(1/2))
    if a > 0:
        R_isco = 3+Z2 - ((3-Z1)*(3+Z1+2*Z2))**(1/2) # for a prograde spin
    else:
        R_isco = 3+Z2 + ((3-Z1)*(3+Z1+2*Z2))**(1/2) # for a retrograde spin
    return R_isco
Risco = []
for a in a_values:
    Risco.append(ISCO_func(a))
    
# marginally bound orbit as a function of a
def mb_func(a):
    R_mb = (2- a +2*(np.sqrt(1-a))) 
    return R_mb

R_margbound = []
for a in a_values:
    R_margbound.append(mb_func(a))
    
# photon circular orbit as a function of a
def ph_func(a):
    R_ph = 2*(1+np.cos((2/3)*np.arccos(-a))) 
    return R_ph

R_phsphere = []
for a in a_values:
    R_phsphere.append(ph_func(a))
    

## How is the radius of the event horizen related to a?

def eventhorizon_func(a):
    R_evt = (1+np.sqrt(1-a**2))
    return R_evt

Event_Horizon = []
for a in a_values:
    Event_Horizon.append(eventhorizon_func(a))
    
# Plot all of the orbits on the same plot

plt.plot(a_values, Risco, linestyle='-', color='mediumvioletred', label='ISCO')
plt.plot(a_values, R_margbound, linestyle='-', color='indigo', label='Marginally Bound Orbit')
plt.plot(a_values, R_phsphere, linestyle='-', color='blue', label='Photon Circular Orbit')
plt.plot(a_values, Event_Horizon, linestyle='-', color='black', label='Event Horizon')
plt.title('Orbital Radius vs. Spin of Black Holes', fontsize=11, color='black') 
plt.xlabel('Spin Parameter, $a$', fontsize=10, color='black') 
plt.ylabel('Orbital Radius($R_g$)', fontsize=10, color='black')
plt.grid(True)
plt.gcf().set_dpi(300) 
plt.savefig('plot.png', dpi=300)  
plt.legend(fontsize = '10')
plt.show()
    


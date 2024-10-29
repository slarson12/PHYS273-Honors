# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 19:39:48 2024

@author: sierr
"""
import matplotlib.pyplot as plt
import numpy as np

# Constants 
G = 6.67 * 10**-11 # Nm^2 / kg^2
c = 2.99 * 10**8 #m/s
Ms = 2 * 10**30 #kg
M = 4.297 * 10**6 * Ms
grav_radius = G*M/(c**2)
R_g = G*M/(c**2)
Rs = 2*G*M/(c**2)

# define range of Radius values
r_values = np.linspace(1*grav_radius, 9*grav_radius, 100)

# define Newtonian orbital velocity fuction
def orbital_velocity_func(r):
    v = np.sqrt((G*M)/r)
    return v

obital_velocity = []
for r in r_values:
    obital_velocity.append(orbital_velocity_func(r))
    
#plot orbital velocity vs radius
plt.plot(r_values/grav_radius, obital_velocity, linestyle='-', color='purple', label='orbital velocity')
plt.title('Orbital Velocity vs Radius',fontsize= 12)
plt.xlabel('Radius from Black hole ($R_g$)',fontsize= 12)
plt.ylabel('Orbital Velocity (m/s)',fontsize= 12)
plt.grid(True)
plt.gcf().set_dpi(300) 
plt.savefig('plot.png', dpi=300)  
plt.legend(fontsize = '9')
plt.show() 

## what does the angular momentum look like as a function of radius??

R_values = np.linspace(1*R_g, 9*R_g, 100)

def angular_momentum_func(R):
    L = np.sqrt(G*M*R)
    return L

angular_momentum = []
for R in R_values:
    angular_momentum.append(angular_momentum_func(R))
    
# now plot L vs radius

plt.plot(R_values/grav_radius, angular_momentum, linestyle='-', color='indigo', label='angular momentum')
plt.title('Angular Momentum vs Radius',fontsize= 12)
plt.xlabel('Radius from Black hole ($R_g$)',fontsize= 12)
plt.ylabel('Angular momentum per unit mass',fontsize= 12)
plt.grid(True)
plt.gcf().set_dpi(300) 
plt.savefig('plot.png', dpi=300)  
plt.legend(fontsize = '10')
plt.show() 

## ******** plot angular momentum for different a values********
def angular_momentum_func1(R):
    L = (0+R**(2/3))**-1
    return L

angular_momentum1 = []
for R in R_values:
    angular_momentum1.append(angular_momentum_func1(R))
    
def angular_momentum_func2(R):
    L = (1+R**(2/3))**-1
    return L

angular_momentum2 = []
for R in R_values:
    angular_momentum2.append(angular_momentum_func2(R))
    
def angular_momentum_func3(R):
    L = (-1+R**(2/3))**-1
    return L

angular_momentum3 = []
for R in R_values:
    angular_momentum3.append(angular_momentum_func3(R))

def angular_momentum_func4(R):
    L = (.5+R**(2/3))**-1
    return L

angular_momentum4 = []
for R in R_values:
    angular_momentum4.append(angular_momentum_func4(R))
    
def angular_momentum_func5(R):
    L = (-.5+R**(2/3))**-1
    return L

angular_momentum5 = []
for R in R_values:
    angular_momentum5.append(angular_momentum_func5(R))

#plot the angular momentum for different a values
plt.plot(R_values, angular_momentum1, linestyle='-', label='a=0')
plt.plot(R_values, angular_momentum2, linestyle='-', label='a=1')
#plt.plot(R_values, angular_momentum3, linestyle='-', label='a=-1')
plt.plot(R_values, angular_momentum4, linestyle='-',label='a=0.5')
plt.plot(R_values, angular_momentum5, linestyle='-', label='a=-0.5')
plt.title('Angular Momentum vs Radius',fontsize= 12)
plt.xlabel('Radius from Black hole ($R_g$)',fontsize= 12)
plt.ylabel('Angular momentum',fontsize= 12)
plt.grid(True)
plt.gcf().set_dpi(300) 
plt.savefig('plot.png', dpi=300)  
plt.legend(fontsize = '10')
plt.show() 

a_values = np.linspace(-1,1,100)

def angular_momentum_func1(a):
    r = 2*R_g
    L = (a+r**(2/3))**-1
    return L

angular_momentum1 = []
for a in a_values:
    angular_momentum1.append(angular_momentum_func1(a))
    
def angular_momentum_func2(a):
    r = 5*R_g
    L = (a+r**(2/3))**-1
    return L

angular_momentum2 = []
for a in a_values:
    angular_momentum2.append(angular_momentum_func2(a))
    
def angular_momentum_func3(a):
    r = 9*R_g
    L = (a+r**(2/3))**-1
    return L

angular_momentum3 = []
for a in a_values:
    angular_momentum3.append(angular_momentum_func3(a))

plt.plot(a_values, angular_momentum1, linestyle='-', label='radius = 2$R_g$')
plt.plot(a_values, angular_momentum2, linestyle='-', label='radius = 5$R_g$')
plt.plot(a_values, angular_momentum3, linestyle='-', label='radius = 9$R_g$')
plt.title('Angular Momentum vs Spin Parameter',fontsize= 12)
plt.xlabel('Spin parameter a',fontsize= 12)
plt.ylabel('Angular Momentum',fontsize= 12)
plt.grid(True)
plt.gcf().set_dpi(300) 
plt.savefig('plot.png', dpi=300)  
plt.legend(fontsize = '11')
plt.show() 

# plot angular velocity proportional to R^(-3/2)

def angular_velocity_func(r):
    omega = r**(-(3/2))
    return omega

angular_velocity = []
for r in r_values:
    angular_velocity.append(angular_velocity_func(r))
    
# plotting
plt.plot(r_values/grav_radius, angular_velocity, linestyle='-', color='darkblue', label='angular velocity')
plt.title('Angular velocity vs Radius',fontsize= 12)
plt.xlabel('Radius from Black hole ($R_g$)',fontsize= 12)
plt.ylabel('Angular velocity',fontsize= 12)
plt.gcf().set_dpi(300) 
plt.savefig('plot.png', dpi=300)  
plt.grid(True)
plt.legend(fontsize = '9')
plt.show() 

## plot the doppler effect of a particle at a fixed radius over 1 period of its oribit

# define theta values over 2 periods
theta_values = np.linspace(0,4*np.pi,100)

# define a doppler shift function
def doppler_func1(theta):
    v_tan = .1*c   #use arbritary velocity
    B = v_tan/c
    gamma = 1/(np.sqrt(1-B**2))
    z = 1/(gamma*(np.sqrt(1+B*np.cos(theta))))
    return z 

doppler1 = []
for theta in theta_values:
    doppler1.append(doppler_func1(theta)) 
    
# use a different v_tan now
    
def doppler_func2(theta):
    v_tan = .2*c   #use different v_tan
    B = v_tan/c
    gamma = 1/(np.sqrt(1-B**2))
    z = 1/(gamma*(np.sqrt(1+B*np.cos(theta))))
    return z 

doppler2 = []
for theta in theta_values:
    doppler2.append(doppler_func2(theta)) 
    
# im sure there are better ways but make a 3rd function for a different v_tan value
def doppler_func3(theta):
    v_tan = .3*c   #new v_tan !
    B = v_tan/c
    gamma = 1/(np.sqrt(1-B**2))
    z = 1/(gamma*(np.sqrt(1+B*np.cos(theta))))
    return z 

doppler3 = []
for theta in theta_values:
    doppler3.append(doppler_func3(theta)) 
    
#what happens if v_tan is very big?? 
def doppler_func8(theta):
    v_tan = .8*c   #new v_tan !
    B = v_tan/c
    gamma = 1/(np.sqrt(1-B**2))
    z = 1/(gamma*(np.sqrt(1+B*np.cos(theta))))
    return z 

doppler8 = []
for theta in theta_values:
    doppler8.append(doppler_func8(theta)) 
    
# plotting the dopplershift versus theta for the different values of velocity
plt.plot(theta_values, doppler1, color='darkblue', label='v= 0.1c')
plt.plot(theta_values, doppler2, color='teal', label='v = 0.2c')
plt.plot(theta_values, doppler3, color='dodgerblue', label='v = 0.3c')
#plt.plot(theta_values/np.pi, doppler8, color='black', label='v = 0.8c')
plt.title('Doppler shift over two Periods of Rotation',fontsize= 12)
plt.xlabel('$\\alpha$ (rad)',fontsize= 12)
plt.ylabel('Z (Femitt/Fobs)',fontsize= 12)
plt.xticks(np.linspace(0, 4*np.pi, 5), ['0', r'$\pi$', r'$2\pi$', r'$3\pi$', r'$4\pi$'])
plt.grid(True)
plt.gcf().set_dpi(300) 
plt.savefig('plot.png', dpi=300)  
plt.legend(fontsize = '8')
plt.show() 

# what would happen if Vtan was sqrt(GM/r)

R1 = 2*R_g 
v = np.sqrt((G*M)/R1) 
print(v) 
#function using this new vtan
def doppler_func10(theta):
    v_tan = np.sqrt((G*M)/R1)   
    B = v_tan/c
    gamma = 1/(np.sqrt(1-B**2))
    z = 1/(gamma*(np.sqrt(1+B*np.cos(theta))))
    return z 
doppler10 = []
for theta in theta_values:
    doppler10.append(doppler_func10(theta)) 

#okay, now do a few other R values
R5 = 6*R_g 
def doppler_func11(theta):
    v_tan = np.sqrt((G*M)/R5)   
    B = v_tan/c
    gamma = 1/(np.sqrt(1-B**2))
    z = 1/(gamma*(np.sqrt(1+B*np.cos(theta))))
    return z 
doppler11 = []
for theta in theta_values:
    doppler11.append(doppler_func11(theta)) 
    
# another R value to plot
R9 = 9*R_g 
def doppler_func12(theta):
    v_tan = np.sqrt((G*M)/R9)   
    B = v_tan/c
    gamma = 1/(np.sqrt(1-B**2))
    z = 1/(gamma*(np.sqrt(1+B*np.cos(theta))))
    return z 
doppler12 = []
for theta in theta_values:
    doppler12.append(doppler_func12(theta)) 
    
plt.plot(theta_values, doppler10, color='blueviolet', label=' R = 2$R_g$')
plt.plot(theta_values, doppler11, color='darkgreen', label=' R = 6$R_g$')
plt.plot(theta_values, doppler12, color='navy', label=' R = 9$R_g$')
plt.title('Doppler shift Parameter Around a Black Hole',fontsize= 12)
plt.xlabel('$\\alpha$ (rad)',fontsize= 12)
plt.ylabel('Z',fontsize= 10)
plt.xticks(np.linspace(0, 4*np.pi, 5), ['0', r'$\pi$', r'$2\pi$', r'$3\pi$', r'$4\pi$'])
plt.grid(True)
plt.gcf().set_dpi(300) 
plt.savefig('plot.png', dpi=300)  
plt.legend(fontsize = '10')
plt.show() 


## GRAVITATIONAL REDSHIFT

def grav_redshift_func(r):
    z_g = (1-(Rs/r))**(-1/2) - 1
    return z_g

grav_redshift = []
for r in r_values:
    grav_redshift.append(grav_redshift_func(r))
    
# plot graviational redhsift vs radius
plt.plot(r_values/grav_radius,grav_redshift, linestyle='-', color='darkgreen', label='gravitational redshift z')
plt.title('Gravitational Redshift vs Radius',fontsize= 13)
plt.xlabel('Radius from Black hole ($R_g$)',fontsize= 12)
plt.ylabel('Z',fontsize= 12)
plt.gcf().set_dpi(300) 
plt.savefig('plot.png', dpi=300)  
plt.legend(fontsize = '10')
plt.show() 


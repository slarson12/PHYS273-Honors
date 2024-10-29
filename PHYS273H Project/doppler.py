# -*- coding: utf-8 -*-
"""
Created on Wed May  1 17:57:10 2024

@author: sierr
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter
from matplotlib.cm import get_cmap


## Doppler shift for inclination angles- eek

#constants
G = 6.67 * 10**-11 # gravitational constant Nm^2 / kg^2
c = 2.99 * 10**8 # speed of light( m/s )
Ms = 2 * 10**30 #solar mass (kg)
M = 4.297 * 10**6 * Ms # mass of Sag A*
R_g = G*M/(c**2) # Gravitational radius of a black hole
Rs = 2*R_g # Schwarzschild radius  (m)
h = 6.626e-34  # Planck's constant (J*s)
k = 1.38e-23  # Boltzmann constant (J/K)

#define numpy symbols
pi = np.pi
sin = np.sin
cos= np.cos

#define alpha, the angle around the circular orbit, ranges from 0 to 2pi
alpha_values = np.linspace(0, 4*pi, 1000)

# define the x,y,z hat values as vectors
x_hat = np.array([1,0,0])
y_hat = np.array([0,0,1])
z_hat = np.array([0,1,0])

#function for the dopper shift that takes into accound the alpha angle around the orbit and the inclination angle theta

def doppler_func(alpha, theta):

    # Define omega as a vector
    omega_vector = np.array([0 , sin(theta), cos(theta)])  # Angular velocity vector
    
    #omega hat unit vector
    omega_magnitude = np.linalg.norm(omega_vector)
    omega_hat = omega_vector / omega_magnitude

    # r = cos(alpha)*x_hat + sin(alpha)*(omega_hat cross x_hat)
    r_vector = cos(alpha)*x_hat + sin(alpha)*np.cross(omega_hat, x_hat) 
    
    # v = omega cross r
    v_vector = np.cross(omega_vector, r_vector)

    #magnitude of v
    v_magnitude = np.linalg.norm(v_vector)

    #cos(phi)
    cos_phi = (np.dot(v_vector, z_hat)) / v_magnitude

    #calcualte velocity v
    v = np.sqrt((G*M) / (2.1*R_g))

    # calculate the Lorentz factor to plug into doppler shift
    B = v/c
    gamma = 1/(np.sqrt(1 - B**2))

    # doppler shift
    z = 1 / (gamma*(np.sqrt(1 + B*cos_phi)))
    return z

# define inclindation angle values 
theta_0 = 0
theta_1 = pi/6
theta_2= pi/4
theta_3 = pi/3
theta_4 = pi/2
theta_5 = 2*pi/3
theta_6 = 3*pi/4
theta_7 = 5*pi/6
theta_8 = pi

# Calculate Doppler shift for each alpha value for the first inclination angle
doppler_values_0 = [doppler_func(alpha, theta_0) for alpha in alpha_values]
doppler_values_1 = [doppler_func(alpha, theta_1) for alpha in alpha_values]
doppler_values_2 = [doppler_func(alpha, theta_2) for alpha in alpha_values]
doppler_values_3 = [doppler_func(alpha, theta_3) for alpha in alpha_values]
doppler_values_4 = [doppler_func(alpha, theta_4) for alpha in alpha_values]
doppler_values_5 = [doppler_func(alpha, theta_5) for alpha in alpha_values]
doppler_values_6 = [doppler_func(alpha, theta_6) for alpha in alpha_values]
doppler_values_7 = [doppler_func(alpha, theta_7) for alpha in alpha_values]
doppler_values_8 = [doppler_func(alpha, theta_8) for alpha in alpha_values]


# Plot theta values 0 to pi/2
plt.plot(alpha_values, doppler_values_0, color='darkviolet', label='inc. angle = 0')
plt.plot(alpha_values, doppler_values_1, color='red', label=r'inc. angle = $\frac{\pi}{6}$')
plt.plot(alpha_values, doppler_values_2, color='mediumvioletred', label=r'inc. angle = $\frac{\pi}{4}$')
plt.plot(alpha_values, doppler_values_3, color='gold', label=r'inc. angle = $\frac{\pi}{3}$')
plt.plot(alpha_values, doppler_values_4, color='slategrey', label=r'inc. angle = $\frac{\pi}{2}$')
plt.xlabel('$\\alpha$ (rad)',fontsize= 12)
plt.ylabel('Doppler Shift ($\\delta$)',fontsize= 12)
plt.title('Doppler Shift as a function of $\\alpha$',fontsize= 13)
plt.xticks(np.linspace(0, 4*np.pi, 5), ['0', r'$\pi$', r'$2\pi$', r'$3\pi$', r'$4\pi$'])
plt.grid(True)
plt.gcf().set_dpi(300) 
plt.savefig('plot.png', dpi=300)  
plt.legend(loc="right", bbox_to_anchor=(1.4, .5), borderaxespad=0,fontsize = 12)
plt.show()

# Plot theta values from the second quadrent of the circle - shows a phase shift
plt.plot(alpha_values, doppler_values_8, color='darkolivegreen', label=r'inc. angle = $\pi$')
plt.plot(alpha_values, doppler_values_7, color='lightseagreen', label=r'inc. angle = $\frac{5\pi}{6}$')
plt.plot(alpha_values, doppler_values_6, color='steelblue', label=r'inc. angle = $\frac{3\pi}{4}$')
plt.plot(alpha_values, doppler_values_5, color='midnightblue', label=r'inc. angle = $\frac{2\pi}{3}$')
plt.xlabel('$\\alpha$ (rad)',fontsize= 12)
plt.ylabel('Doppler Shift ($\\delta$)',fontsize= 12)
plt.title('Doppler Shift as a function of $\\alpha$',fontsize= 13)
plt.xticks(np.linspace(0,4*pi, 5), ['0', r'$\pi$', r'$2\pi$', r'$3\pi$', r'$4\pi$'])
plt.gcf().set_dpi(300) 
plt.savefig('plot.png', dpi=300)  
plt.grid(True)
plt.legend(loc="right", bbox_to_anchor=(1.4,.5), borderaxespad=0,fontsize= 12)
plt.show()


# Okay - now try to simulate an emission spectrum , then apply doppler shifts 

# ------------------Calculate blackbody spectrum------------------------

# Equation for plotting Plank's law

def plancks_law_equation(wavelengths, temp):
    intensity_PLE = (2*h*c**2) / (wavelengths**5 *(np.exp((h*c) / (wavelengths*k*temp)) - 1))
    return intensity_PLE

# Wavelengths in the X-ray range (0.03 nm to 3 nm)
wavelengths_emit = np.linspace(0.03e-9, 3e-9, 1000)  # Wavelengths in meters

#Temp of the accretion disk (in K)
temperature = 1e7 # 10^7 K

#Find intensity using planck's law
intensity_emit = plancks_law_equation(wavelengths_emit, temperature)

# Plot  of the blackbody spectrum
plt.plot(wavelengths_emit*1e9, intensity_emit, label='X-ray Spectrum', color='blue')  # Convert to nm for plotting
plt.xlabel('Wavelength (nm)',fontsize= 12)
plt.ylabel('Intensity (W sr$^{-1}$ m$^{-3}$)',fontsize= 12)
plt.title('Simulated X-ray Spectrum (Blackbody)',fontsize= 12)
plt.legend()
plt.gcf().set_dpi(300) 
plt.savefig('plot.png', dpi=300)  
plt.grid(True)
plt.show()


# define function that calculates the wavelengths observed due to the doppler shift
# then plug those wavelengths into planks law to plot a shifted spectrum

#Define the alpha values we will use
alpha_values = np.linspace(0,2*pi,1000)

#define the wavelength emitted array
wavelengths_emit = np.linspace(0.03e-9, 3e-9, 1000)  # lamda in meters

#define the frequency emitted
freq_emit = c/wavelengths_emit

#define new coordinate system
x_hat = np.array([1,0,0])
y_hat = np.array([0,1,0])
z_hat = np.array([0,0,1])

def shifted_emission_func(alpha, theta, temp, wavelengths_emit):
    # define omega as a vector
    omega_vector = np.array([0, sin(theta), cos(theta)])  # Angular velocity vector

    #compute the omega hat unit vector
    omega_magnitude = np.linalg.norm(omega_vector)
    omega_hat = omega_vector / omega_magnitude

    # r is a vector with r = cos(alpha)*x_hat + sin(alpha)*(omega_hat cross x_hat)
    r_vector = cos(alpha) * x_hat + sin(alpha) * np.cross(omega_hat, x_hat)

    #  v = omega cross r
    v_vector = np.cross(omega_vector, r_vector)

    #magnitude of v
    v_magnitude = np.linalg.norm(v_vector)

    #cos(phi), phi is the angle between v_vector and z axis
    cos_phi = (np.dot(v_vector, z_hat)) / v_magnitude

    # Calculate velocity v
    v = np.sqrt((G*M) / (2*R_g))

    #Lorentz factor for doppler shift
    B = v/c
    gamma = 1/(np.sqrt(1 - B**2))

    # Calculate the Doppler shift
    z = 1/(gamma*(np.sqrt(1 + B*cos_phi)))

    # freq observed using the doppler shift and freq emitted
    freq_obs = freq_emit*z

    # find the wavelength observed from the freq observed
    wavelengths_obs = c / freq_obs

    #plancks function with the observed wavelengths
    intensity_obs = (2*h*c**2) / (
                wavelengths_obs**5 * (np.exp((h*c) / (wavelengths_obs*k*temp)) - 1))
    return intensity_obs


# define the shifted emission function for different theta values
shifted_emission_func0 = shifted_emission_func(0, theta_0, temperature, wavelengths_emit)
shifted_emission_func1 = shifted_emission_func(0, theta_1, temperature, wavelengths_emit)
shifted_emission_func2 = shifted_emission_func(0, theta_2, temperature, wavelengths_emit)
shifted_emission_func3 = shifted_emission_func(0, theta_3, temperature, wavelengths_emit)
shifted_emission_func4 = shifted_emission_func(0, theta_4, temperature, wavelengths_emit)
shifted_emission_func5 = shifted_emission_func(0, theta_5, temperature, wavelengths_emit)
shifted_emission_func6 = shifted_emission_func(0, theta_6, temperature, wavelengths_emit)
shifted_emission_func7 = shifted_emission_func(0, theta_7, temperature, wavelengths_emit)
shifted_emission_func8 = shifted_emission_func(0, theta_8, temperature, wavelengths_emit)

#plot emission spectrum of theta values 0 to pi/2
plt.plot(wavelengths_emit * 1e9, shifted_emission_func0, color='blue', label='Original Spectrum')
#plt.plot(wavelengths_emit * 1e9, shifted_emission_func1, color='maroon', label=r'inc. angle = $\frac{\pi}{6}$')
plt.plot(wavelengths_emit * 1e9, shifted_emission_func2, color='magenta', label=r'inc. angle = $\frac{\pi}{4}$')
#plt.plot(wavelengths_emit * 1e9, shifted_emission_func3, color='teal', label=r'inc. angle = $\frac{\pi}{3}$')
plt.plot(wavelengths_emit * 1e9, shifted_emission_func4, color='slategrey', label=r'inc. angle = 0')

plt.xlabel('Wavelength (nm)',fontsize= 12)
plt.ylabel('Intensity (W sr$^{-1}$ m$^{-3}$)',fontsize= 12)
plt.title('Simulated Blackbody X-ray Spectrum',fontsize= 12)
plt.gcf().set_dpi(300) 
plt.savefig('plot.png', dpi=300)  
plt.grid(True)
plt.legend()
plt.show()


# Now apply the gravitational redshift function to the emitted wavelengths

# look at the gravitational redshift at the Risco R = 2.1 Rg

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter

# Constants
G = 6.67 * 10**-11  # Gravitational constant (Nm^2 / kg^2)
c = 2.99 * 10**8  # Speed of light (m/s)
Ms = 2 * 10**30  # Solar mass (kg)
M = 4.297 * 10**6 * Ms  # Mass of Sag A* (kg)
h = 6.626e-34  # Planck's constant (J*s)
k = 1.38e-23  # Boltzmann constant (J/K)

# Define Planck's function
def plancks_law(wavelengths, temp):
    intensity = (2 * h * c**2) / (wavelengths**5 * (np.exp((h * c) / (wavelengths * k * temp)) - 1))
    return intensity

# Define gravitational redshift function
def gravitational_redshift(wavelengths, R):
    Ve = np.sqrt((2 * G * M) / R)  # Escape velocity
    z_g = 1 / (np.sqrt(1 - (Ve / c)**2)) - 1  # Gravitational redshift
    return wavelengths * (1 + z_g)

# Define Gaussian function for broadening
def gaussian(x, mu, sigma, intensity):
    return intensity * np.exp(-0.5 * ((x - mu) / sigma)**2)

# Parameters
temperature = 1e7  # Temperature of the black hole corona (K)
central_wavelength = 1.5e-9  # Central wavelength for Gaussian broadening (m)
sigma = 0.03e-9  # Standard deviation for Gaussian broadening (m)

# Define radii for different gravitational redshifts
radii = [2.1 * (G * M / c**2), 6.0 * (G * M / c**2), 9.0 * (G * M / c**2)]

# Get original intensity data from Planck's function
wavelengths_emit = np.linspace(0.03e-9, 3e-9, 1000)
intensity_emit = plancks_law(wavelengths_emit, temperature)

# Define colors
colors = ['indigo', 'darkred', 'green']  # Setting the first color to indigo

# Plot the original spectrum
plt.figure(figsize=(10, 6))
plt.plot(wavelengths_emit * 1e9, intensity_emit, label='Original Spectrum', color='blue')

# Calculate and plot the gravitational redshift spectra for each radius
for i, R in enumerate(radii):
    # Calculate the gravitational redshift for the given radius
    wavelengths_observed = gravitational_redshift(wavelengths_emit, R)

    # Use Gaussian broadening to simulate the broadening of intensity data
    broadened_spectrum = gaussian_filter(intensity_emit, sigma)

    # Plot the broadened spectrum
    plt.plot(wavelengths_observed * 1e9, broadened_spectrum, label=f'Gravitational Redshift at R = {R/(G*M/c**2):.1f} $R_g$', color=colors[i], linewidth=1)

# Customize the plot
plt.xlabel('Wavelength (nm)', fontsize=12)
plt.ylabel('Intensity', fontsize=12)
plt.title('Blackbody Spectrum with Gravitational Redshift at Different Radii', fontsize=12)
plt.legend(fontsize= 14)
plt.grid(True)
plt.xlim(-0.03,10)
plt.gcf().set_dpi(300) 
plt.savefig('plot.png', dpi=300)
plt.show()



# Comparison of the effect of the doppler shift and the gravitational redshift on the emission spectrum

# get original intensity data using planck's function
wavelengths_emit = np.linspace(0.03e-9, 3e-9, 1000)  # redefine wavelegnths emitted
intensity_emit = plancks_law(wavelengths_emit, temperature) # plug into planck's law

R = 2.1* R_g  # Radius for gravitational redshift
# wavelengths from gravitational redshift
Ve = np.sqrt((2 * G * M) / R)  # escape velocity
z_g = 1 / (np.sqrt(1 - (Ve / c)**2)) - 1  # gravredshift
wavelengths_gravitational = wavelengths_emit  * (1 + z_g)

# apply both doppler shift and gravitational redshift  
intensity_gravitational = gaussian_filter(intensity_emit, sigma=0.01e-9)  

# plot the spectra with different effects
# assume we are looking from the side(theta = 0) so we have the max possible dopplershift
plt.figure(figsize=(10, 6))
plt.plot(wavelengths_emit * 1e9, intensity_emit, label='Original Spectrum',color='blue')
plt.plot(wavelengths_emit * 1e9, shifted_emission_func4, color='slategrey', label=r'Spectrum with Doppler Shift')
plt.plot(wavelengths_gravitational * 1e9, intensity_gravitational, label='Spectrum with Gravitational Redshift at 2.1 $R_g$',color='indigo')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Intensity')
plt.title('Blackbody Spectrum with Doppler Shift and Gravitational Redshift')
plt.legend(fontsize= 14)
plt.xlim(-.03,10)
plt.gcf().set_dpi(300) 
plt.savefig('plot.png', dpi=300)  
plt.grid(True)
plt.show()




















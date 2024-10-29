# -*- coding: utf-8 -*-
"""
Created on Sat May  4 19:38:15 2024

@author: sierr
"""
import matplotlib.pyplot as plt
import numpy as np

# Constants
G = 6.67 * 10**-11  # Gravitational constant Nm^2 / kg^2
c = 2.99 * 10**8  # Speed of light (m/s)
Ms = 2 * 10**30  # Solar mass (kg)
M = 4.297 * 10**6 * Ms  # Mass of Sag A*
R_g = G * M / (c**2)  # Gravitational radius of a black hole
h = 6.626e-34  # Planck's constant (J*s)
k = 1.38e-23  # Boltzmann constant (J/K)

# Define numpy symbols
pi = np.pi
sin = np.sin
cos = np.cos

# Define the x, y, z hat values as vectors
x_hat = np.array([1, 0, 0])
y_hat = np.array([0, 0, 1])
z_hat = np.array([0, 1, 0])

# another way of making an emission spectrum with different inc angles

# Define emitted wavelengths
wavelengths_emit = np.linspace(0.03e-9, 3e-8, 100)  # meters

# Define the frequency emitted
freq_emit = c / wavelengths_emit

# Function for the Doppler shift
def doppler_func(alpha, theta, wavelength):
    # Define omega as a vector
    omega_vector = np.array([0, sin(theta), cos(theta)])  # Angular velocity vector

    # Omega hat unit vector
    omega_magnitude = np.linalg.norm(omega_vector)
    omega_hat = omega_vector / omega_magnitude

    # Calculate r_vector
    r_vector = cos(alpha) * x_hat + sin(alpha) * np.cross(omega_hat, x_hat)

    # Calculate v_vector
    v_vector = np.cross(omega_vector, r_vector)

    # Calculate the magnitude of v_vector
    v_magnitude = np.linalg.norm(v_vector)

    # Calculate cos(phi)
    cos_phi = np.dot(v_vector, z_hat) / v_magnitude

    # Calculate velocity v
    v = np.sqrt((G * M) / (2.1 * R_g))

    # Calculate the Lorentz factor to plug into Doppler shift
    B = v / c
    gamma = 1 / np.sqrt(1 - B ** 2)

    # Calculate the Doppler shift
    z = 1 / (gamma * np.sqrt(1 + B * cos_phi))

    # Find the frequency received based on the Doppler shift
    freq_doppler_shifted = (1 + z) * freq_emit

    # Find the wavelength observed from the frequency observed
    wavelengths_shifted = c / freq_doppler_shifted
    return wavelengths_shifted

# Function to calculate the intensity using Planck's law
def plancks_law(wavelength, temperature):
    intensity = (2 * h * c ** 2) / (wavelength ** 5 * (np.exp((h * c) / (wavelength * k * temperature)) - 1))
    return intensity

def observed_intensity(alpha, theta, wavelengths_emit):
    # Calculate Doppler-shifted wavelengths for the specific alpha and theta
    wavelengths_shifted = np.array([doppler_func(alpha, theta, wavelength) for wavelength in wavelengths_emit])

    # Calculate the intensity using Planck's law for each shifted wavelength
    intensities = np.array([plancks_law(wavelength, temperature) for wavelength in wavelengths_shifted])

    return intensities

# Define the temperature of the blackbody (in Kelvin)
temperature = 1e6  # K

# Define alpha values
alpha_values = np.linspace(0, 2 * pi, 100)
theta_values = [0,pi/4,pi/2,pi]

# plot the blackbody spectrum
plt.figure(figsize=(10, 6))

# calculate the total intensity but with the unshifted spectrum
total_intensity = np.zeros_like(wavelengths_emit)
for alpha in alpha_values:
    intensity = observed_intensity(0, 0, wavelengths_emit)  
    average_intensity = np.mean(intensity, axis=0)  
    total_intensity += average_intensity


plt.plot(wavelengths_emit * 1e9, total_intensity, color='blue', label='Unshifted Spectrum')
# calculate the total intensity for the given alpha values
total_intensity = np.zeros_like(wavelengths_emit)
for alpha in alpha_values:
    intensity = observed_intensity(alpha, theta_values[0], wavelengths_emit)  # calculate intensity for fixed theta
    average_intensity = np.mean(intensity, axis=0)  # calculate average intensity across different angles
    total_intensity += average_intensity

plt.plot(wavelengths_emit * 1e9, total_intensity, color='red', label='inc.angle = 0')


total_intensity = np.zeros_like(wavelengths_emit)
for alpha in alpha_values:
    intensity = observed_intensity(alpha, theta_values[1], wavelengths_emit)  
    average_intensity = np.mean(intensity, axis=0)  
    total_intensity += average_intensity

plt.plot(wavelengths_emit * 1e9, total_intensity, color='orange', label=r'inc.angle = $\frac{\pi}{4}$')

total_intensity = np.zeros_like(wavelengths_emit)
for alpha in alpha_values:
    intensity = observed_intensity(alpha, theta_values[2], wavelengths_emit)  
    average_intensity = np.mean(intensity, axis=0) 
    total_intensity += average_intensity

plt.plot(wavelengths_emit * 1e9, total_intensity, color='green', label=r'inc. angle = $\frac{\pi}{2}$')

plt.xlabel('Wavelength (nm)')
plt.ylabel('Intensity (W sr$^{-1}$ m$^{-3}$)')
plt.title('Blackbody Spectrum with Doppler Shift')
plt.grid(True)
plt.gcf().set_dpi(300) 
plt.savefig('plot.png', dpi=300) 
plt.legend()
plt.show()

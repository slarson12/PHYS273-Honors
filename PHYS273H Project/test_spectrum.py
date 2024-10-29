# -*- coding: utf-8 -*-
"""
Created on Mon May  6 16:17:04 2024

@author: sierr
"""
'''
import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.67 * 10**-11  # Gravitational constant Nm^2 / kg^2
c = 2.99 * 10**8  # Speed of light (m/s)
Ms = 2 * 10**30  # Solar mass (kg)
M = 4.297 * 10**6 * Ms  # Mass of Sag A*
R_g = G * M / (c**2)  # Gravitational radius of the black hole
h = 6.626e-34  # Planck's constant (J*s)
k = 1.38e-23  # Boltzmann constant (J/K)

keV_to_J = 1.60218e-16  # Conversion factor from keV to Joules

# Define emitted energies (in keV)
energies_emit_keV = np.linspace(6.0, 6.8, 1000)  # keV

# Function to calculate the intensity using Planck's law
def plancks_law(energy_keV, temperature_keV):
    wavelength = (h * c) / (energy_keV * keV_to_J)  # Convert energy to wavelength (m)
    intensity = (2 * h * c**2) / (wavelength**5 * (np.exp((h * c) / (wavelength * temperature_keV * keV_to_J)) - 1))
    return intensity

# Function to generate the Iron K-alpha line profile
def iron_line_profile(energies_emit_keV, temperature_keV):
    # Parameters for iron K-alpha line
    iron_line_center_keV = 6.4  # keV
    iron_line_width_keV = 0.1   # keV
    iron_line_amplitude = 1e-13  # Intensity amplitude of Iron K-alpha line

    # generate the line profile (Gaussian)
    line_profile = iron_line_amplitude * np.exp(-(energies_emit_keV - iron_line_center_keV)**2 / (2 * iron_line_width_keV**2))
    return line_profile

# Generate the iron K-alpha line profile
temperature_keV = 1.0  # keV (temperature of the accretion disk)
line_profile = iron_line_profile(energies_emit_keV, temperature_keV)

# Plot the Iion K-alpha line 
plt.figure(figsize=(10, 6))
plt.plot(energies_emit_keV, line_profile, color='blue')
plt.xlabel('Energy (keV)')
plt.ylabel('Intensity (arbitrary units)')
plt.title('Iron K-alpha Line Profile')
plt.grid(True)
plt.show()


#***********maybe try another way***********

# Define emitted energies (in keV)
energies_emit_keV = np.linspace(6.0, 6.8, 1000)  # keV

# Function to generate the Iron K-alpha line profile
def custom_iron_line_profile(energies_emit_keV):
    # Parameters for Iron K-alpha line
    iron_line_center_keV = 6.4  # keV
    iron_line_width_keV = 0.01  # keV
    iron_line_amplitude = 100 # Intensity amplitude of Iron K-alpha line

    # Generate the Gaussian part for the peak at 6.4 keV
    gaussian_part = iron_line_amplitude * np.exp(-(energies_emit_keV - iron_line_center_keV)**2 / (2 * iron_line_width_keV**2))

    # Generate the step function to make it zero after 6.4 keV
    step_function_part = np.where(energies_emit_keV <= iron_line_center_keV, gaussian_part, 0)

    return step_function_part

# Generate the custom Iron K-alpha line profile
line_profile = custom_iron_line_profile(energies_emit_keV)

# Plot the custom Iron K-alpha line profile
plt.figure(figsize=(10, 6))
plt.plot(energies_emit_keV, line_profile, color='blue')
plt.xlabel('Energy (keV)')
plt.ylabel('Intensity (arbitrary units)')
plt.title('Custom Iron K-alpha Line Profile')
plt.grid(True)
plt.show()

# Find peak wavelength
peak_energy_keV = 6.4

# Convert peak energy from keV to joules
peak_energy_J = peak_energy_keV * 1.60218e-16  # 1 keV = 1.60218e-16 J

# Calculate the corresponding peak wavelength in meters
peak_wavelength_m = h * c / peak_energy_J

print("Peak wavelength:", peak_wavelength_m, "m")

# Define the x, y, z hat values as vectors
x_hat = np.array([1, 0, 0])
y_hat = np.array([0,1 , 0])
z_hat = np.array([0, 0, 1])

# Define emitted wavelengths
wavelengths_emit = np.linspace(0.03e-9, 3e-9, 1000)  # meters

# Define the frequency emitted
freq_emit = c / wavelengths_emit

# Function for the Doppler shift
def doppler_func(alpha, theta, line_profile):
    # Define omega as a vector
    omega_vector = np.array([0, np.sin(theta), np.cos(theta)])  # Angular velocity vector

    # Omega hat unit vector
    omega_magnitude = np.linalg.norm(omega_vector)
    omega_hat = omega_vector / omega_magnitude

    # Calculate r_vector
    r_vector = np.cos(alpha) * x_hat + np.sin(alpha) * np.cross(omega_hat, x_hat)

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

    # Calculate the Doppler shift for each emitted energy
    z = 1 / (gamma * np.sqrt(1 + B * cos_phi))

    # Calculate Observed wavelengths
    wavelengths_observed = wavelengths_emit * (1 + z)

    # Convert observed wavelengths back to energies in keV
    energies_observed_keV = (h * c) / wavelengths_observed / 1.60218e-16

    # Apply Doppler shift to the line profile
    shifted_line_profile = line_profile(energies_observed_keV)
    return shifted_line_profile

# Function to calculate the intensity using the new line profile
def observed_intensity(alpha, theta, line_profile):
    # Apply Doppler-shifted line profile for the specific alpha and theta
    shifted_line_profile = doppler_func(alpha, theta, line_profile)

    # Intensity is just the shifted line profile
    return shifted_line_profile

# Define the new line profile function
def new_line_profile(energies_keV):
    # Create a line profile with a sharp peak at 6.4 keV and zero elsewhere
    peak_energy_keV = 6.4
    peak_intensity = 1.0
    sigma = 0.8  # Standard deviation of the Gaussian peak

    # Calculate the Gaussian profile
    gaussian_profile = peak_intensity * np.exp(-((energies_keV - peak_energy_keV) / sigma) ** 2 / 2)

    return gaussian_profile

# Define alpha values
alpha_values = np.linspace(0, 2 * np.pi, 100)
theta_values = [0, np.pi / 4, np.pi / 2]

# Plot the blackbody spectrum with the new line profile
plt.figure(figsize=(10, 6))

# Plot the unshifted spectrum
total_intensity = np.zeros_like(wavelengths_emit)
for alpha in alpha_values:
    intensity = observed_intensity(np.pi/2, np.pi/2, custom_iron_line_profile)  # Calculate intensity with the new line profile
    total_intensity += intensity

plt.plot((h * c) / wavelengths_emit / 1.60218e-16, total_intensity, color='blue', label='Unshifted Spectrum')

# Calculate the total intensity for the given alpha values

for theta in theta_values:
    total_intensity = np.zeros_like(wavelengths_emit)
    for alpha in alpha_values:
        intensity = observed_intensity(alpha, theta, new_line_profile)  # Calculate intensity with the new line profile
        total_intensity += intensity

    plt.plot((h * c) / wavelengths_emit / 1.60218e-16, total_intensity, label=f'inc. angle = {theta / np.pi:.2f}π')
    
plt.xlabel('Energy (keV)')
plt.ylabel('Intensity')
plt.title('Blackbody Spectrum with Doppler-Shifted Line Profile')
plt.grid(True)
plt.xlim(0,20)
plt.legend()
plt.show()
'''

import matplotlib.pyplot as plt
import numpy as np

# Constants
G = 6.67 * 10**-11  # Gravitational constant (Nm^2 / kg^2)
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
    z = 1 / (gamma * (1 + B * cos_phi)) - 1

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
temperature = 2.8e6  # K

# Define alpha values
alpha_values = np.linspace(0, 2 * pi, 100)
theta_values = [0, pi/4,pi/3]

# Plot the blackbody spectrum
plt.figure(figsize=(10, 6))

# Calculate the total intensity but with the unshifted spectrum
total_intensity = np.zeros_like(wavelengths_emit)
for alpha in alpha_values:
    intensity = plancks_law(wavelengths_emit, temperature)
    total_intensity += intensity

plt.plot(wavelengths_emit * 1e9, total_intensity, color='blue', label='Unshifted Spectrum')

# Calculate the total intensity for the given alpha values
for theta in theta_values:
    total_intensity = np.zeros_like(wavelengths_emit)
    for alpha in alpha_values:
        intensity = observed_intensity(alpha, theta, wavelengths_emit)
        average_intensity = np.mean(intensity, axis=0)
        total_intensity += average_intensity

    plt.plot(wavelengths_emit * 1e9, total_intensity, label=f'inc.angle = {theta / pi:.2f}π')

plt.xlabel('Wavelength (nm)')
plt.ylabel('Intensity (W sr$^{-1}$ m$^{-3}$)')
plt.xlim(-.3,15)
plt.title('Blackbody Spectrum with Doppler Shift')
plt.grid(True)
plt.legend()
plt.show()

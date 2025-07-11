# Will be using this python script to generate a trend line for the Short Hydrophobe surfactant data (IFT vs surfactant concentration)


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

##Exponential Regression
# Data from Paper (doi) -> https://dx.doi.org/10.1021/acs.energyfuels.0c02720
surfactant_concentration_wt_frac = np.array([0, 0.000008, 0.000016, 0.000031, 0.000062, 0.00013, 0.00025, 0.0005, 0.001, 0.0025, 0.005, 0.01])
wt_percent = []
for x in surfactant_concentration_wt_frac:
    wt_percent.append(x*100)
surfactant_concentration_wt_percent = np.array(wt_percent)
ift_dyne_per_cm = np.array([15.8, 6.19, 5.79, 3.3, 1.8, 0.098, 0.025, 0.023, 0.023, 0.023, 0.025, 0.025])

# # Define exponential function
# def exp_func(x, a, b):
#     return a * np.exp(b * x)

# Define exponential function
def custom_func(x, a, b, c):
    return (a/(x+b)) + c

# Fit the curve (exponential regression)
params, _ = curve_fit(custom_func, surfactant_concentration_wt_frac, ift_dyne_per_cm, bounds=(0, np.inf))
a, b, c = params

print(f"a = {a}")
print(f"b = {b}")
print(f"c = {c}")

# Generate values for the fitted curve
x_fit = np.linspace(min(surfactant_concentration_wt_frac), max(surfactant_concentration_wt_frac), 100)
y_fit = custom_func(x_fit, a, b, c)
y_fit_r_squared_calc = custom_func(surfactant_concentration_wt_frac, a, b, c)

for i in range(len(surfactant_concentration_wt_frac)):
    yf_val = custom_func(surfactant_concentration_wt_frac[i], a, b, c)
    print(f"For x = {surfactant_concentration_wt_frac[i]}, y = {ift_dyne_per_cm[i]} and y_fit_value = {yf_val}")
    print(f"error= {ift_dyne_per_cm[i] - yf_val}")

# R^2 calculation
ss_res = np.sum((ift_dyne_per_cm - y_fit_r_squared_calc) ** 2)
ss_tot = np.sum((ift_dyne_per_cm - np.mean(ift_dyne_per_cm)) ** 2)
r_squared = 1 - (ss_res / ss_tot)

print(f"R² = {r_squared:.4f}")

# Plot original data and exponential trendline
plt.scatter(surfactant_concentration_wt_frac, ift_dyne_per_cm, label='Data')
plt.plot(x_fit, y_fit, color='red', label='Trendline')
plt.legend()
plt.xlabel('surfactant concentration (weight fraction)')
plt.ylabel('Interfacial tension (dyne/cm)')
plt.title('IFT (dyne/cm) vs Surfactant Concentration (weight fraction)')
plt.grid(True)
plt.show()


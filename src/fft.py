import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from scipy.signal.windows import hamming, blackman

# Folder settings (as in your script)
folder_name = 'output_TC2'
file_name = f'{folder_name}/coeff_2.csv'

# --- PHYSICAL PARAMETERS OF BENCHMARK (Test Case 2D-2) ---
U_m = 1.5              # Maximum velocity at center [m/s]
U_mean = 2 * U_m / 3   # Mean velocity (parabolic profile) [m/s] -> 1.0 m/s
D = 0.1                # Cylinder diameter [m]
timestep = 0.005       # Simulation time step [s]

# Benchmark bounds for vortex shedding frequency
lower_bound = 0.2950
upper_bound = 0.3050
lower_bound_freq = lower_bound * U_mean / D
upper_bound_freq = upper_bound * U_mean / D

# Read data
col_names = ['Step', 'Drag_Coefficient', 'Lift_Coefficient']
df = pd.read_csv(file_name, header=None, names=col_names)

# Time vector
time = df['Step'] * timestep

# --- 1. TRANSIENT REMOVAL ---
# Isolate only the periodic regime part of the simulation.
# Modify t_start if you see that your wake stabilizes before or after.
t_start = 3.0 
mask = time >= t_start
time_regime = time[mask].values
lift_regime = df['Lift_Coefficient'][mask].values

# Subtract the mean from the signal to eliminate the "DC component" (useless peak at 0 Hz)
lift_centered = lift_regime - np.mean(lift_regime)

# --- 2. FFT ANALYSIS ---
N = len(lift_centered)         # Total number of samples analyzed
delta_step = df['Step'].iloc[1] - df['Step'].iloc[0] 
T = delta_step * timestep                             # 10 * 0.005 = 0.05 s

window = hamming(N)
lift_windowed = lift_centered * window

# Compute Fourier transform and associated frequencies
yf = fft(lift_centered)
xf = fftfreq(N, T)[:N//2]      # Take only positive frequencies
    
# Compute spectrum amplitude
amplitude = (3.0/N) * np.abs(yf[0:N//2])

# --- 3. FREQUENCY AND STROUHAL EXTRACTION ---
# Find the index of maximum amplitude to identify the dominant frequency
idx_max = np.argmax(amplitude)
f_vortex = xf[idx_max]

# Strouhal number calculation
St = (f_vortex * D) / U_mean

# --- 4. SCREEN OUTPUT ---
print("\n" + "="*50)
print(" FFT ANALYSIS - TEST CASE 2D-2 (Unsteady)")
print("="*50)
print(f"Vortex shedding frequency (f) : {f_vortex:.4f} Hz")
print(f"Mean velocity (U_mean)        : {U_mean:.4f} m/s")
print(f"Cylinder diameter (D)         : {D:.4f} m")
print("-" * 50)
print(f"STROUHAL NUMBER (St)          : {St:.4f}")
print(f"Benchmark Range (Table 4)     : [0.2950 - 0.3050]")
print("="*50 + "\n")

# --- 5. SPECTRUM PLOT (Visual Verification) ---
plt.figure(figsize=(10, 5))
plt.plot(xf, amplitude, color='black', linewidth=1.5)
plt.axvline(f_vortex, color='black', linestyle='--', label=f'Peak frequency: {f_vortex:.3f} Hz')
plt.axvspan(lower_bound_freq, upper_bound_freq, color='green', alpha=0.3, label='Benchmark Bounds (0.2950 - 0.3050)')


plt.title('Frequency Spectrum (FFT) of Lift Coefficient - Steady State Phase')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude')
plt.xlim(0, 8)  # Limit visualization to first 10 Hz to see the peak better
plt.minorticks_on()
plt.grid(which='major', color='gray', linestyle='-', alpha=0.5)
plt.grid(which='minor', color='gray', linestyle=':', alpha=0)

plt.savefig('fft_spectrum.png')


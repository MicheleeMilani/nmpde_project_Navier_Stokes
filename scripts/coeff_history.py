import pandas as pd
import matplotlib.pyplot as plt

# Folder name
folder_name = 'output_TC1'

# File name
file_name = f'{folder_name}/coeff_2.csv'

# Read the CSV file
# Since the file has no header row, we use header=None and assign column names ourselves
col_names = ['Step', 'Drag_Coefficient', 'Lift_Coefficient']
df = pd.read_csv(file_name, header=None, names=col_names)

# Create a figure and two plots (subplots) arranged one above the other (2 rows, 1 column)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
timestep = 0.005

if file_name == 'output_TC2/coeff_2.csv':
    cd_lower, cd_upper = 3.22, 3.24
    cl_lower, cl_upper = 0.99, 1.01
elif file_name == 'output_TC3/coeff_2.csv':
    cd_lower, cd_upper = 2.93, 2.97
    cl_lower, cl_upper = 0.47, 0.49
elif file_name == 'output_TC1/coeff_2.csv':
    cd_lower, cd_upper = 5.57, 5.59
    cl_lower, cl_upper = 0.0104, 0.0110
else:
    raise ValueError("Unknown file name. Please check the file name and update the bounds accordingly.")

# --- First plot: Drag Coefficient ---
ax1.plot(df['Step']*timestep, df['Drag_Coefficient'], color='black', linestyle='-')
ax1.axhspan(cd_lower, cd_upper, color='green', alpha=0.3, label='Benchmark Bounds (5.57 - 5.59)')
ax1.set_title('Drag Coefficient Time History')
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Drag Coefficient')

ax1.set_xlim(0, 3)
#ax1.set_ylim(0, 8.1)

ax1.minorticks_on() # Enable minor ticks
ax1.grid(which='major', color='gray', linestyle='-', linewidth=0.8) # Major grid
ax1.grid(which='minor', color='gray', linestyle=':', linewidth=0.5, alpha=0.7) # Minor grid

# --- Second plot: Lift Coefficient ---
ax2.plot(df['Step']*timestep, df['Lift_Coefficient'], color='black', linestyle='-')
ax2.axhspan(cl_lower, cl_upper, color='green', alpha=0.3, label='Benchmark Bounds (0.0104 - 0.0110)')
ax2.set_title('Lift Coefficient Time History')
ax2.set_xlabel('Time [s]')
ax2.set_ylabel('Lift Coefficient')

ax2.set_xlim(0, 10)
#ax2.set_ylim(-0.01, 0.039)

# ACTIVATE MINOR GRID FOR ax2
ax2.minorticks_on() # Enable minor ticks
ax2.grid(which='major', color='gray', linestyle='-', linewidth=0.8) # Major grid
ax2.grid(which='minor', color='gray', linestyle=':', linewidth=0.5, alpha=0.7) # Minor grid

# Properly space the plots to avoid overlapping text
plt.tight_layout()

# Save the image (optional)
plt.savefig('coefficienti_plot.png')

# Display the plots on screen
plt.show()
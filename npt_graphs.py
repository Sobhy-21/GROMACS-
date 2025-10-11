#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import os

# Configuration
files = {
    "Temperature": "npt_temperature.xvg",
    "Pressure": "npt_pressure.xvg",
    "Potential": "npt_potential.xvg"
}

# Recommended thresholds
target_temp = 300.0         # K
temp_tol = 2.0               # K
pressure_tol = 200.0         # bar, pressure is noisy
potential_tol = 1e4          # arbitrary unit for detecting drift

summary = []

def read_xvg(filename):
    """Read XVG file, skip comments and return time and values"""
    data = []
    with open(filename) as f:
        for line in f:
            if line.startswith(("#", "@")):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                data.append([float(parts[0]), float(parts[1])])
    return np.array(data)

for key, filename in files.items():
    if not os.path.isfile(filename):
        print(f"File {filename} not found! Skipping {key}.")
        continue

    data = read_xvg(filename)
    time = data[:,0]
    values = data[:,1]

    mean = np.mean(values)
    std = np.std(values)

    # Plot
    plt.figure(figsize=(8,4))
    plt.plot(time, values, label=f'{key}', color='tab:blue')
    plt.xlabel('Time (ps)')
    plt.ylabel(key)
    plt.title(f'{key} vs Time')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"{key.lower()}_plot.png")
    plt.show()  # <-- Display graph interactively

    # Comment / recommendation
    comment = ""
    ok = True

    if key == "Temperature":
        if abs(mean - target_temp) > temp_tol:
            comment = f"Temperature mean {mean:.2f} K deviates from target {target_temp} K."
            ok = False
        else:
            comment = f"Temperature is stable around target: mean={mean:.2f} ± {std:.2f} K."
    elif key == "Pressure":
        if std > pressure_tol:
            comment = f"Pressure fluctuations are large: mean={mean:.2f} ± {std:.2f} bar."
            ok = False
        else:
            comment = f"Pressure is stable: mean={mean:.2f} ± {std:.2f} bar."
    elif key == "Potential":
        drift = values[-1] - values[0]
        if abs(drift) > potential_tol:
            comment = f"Potential energy shows drift of {drift:.2f} units."
            ok = False
        else:
            comment = f"Potential energy is stable: mean={mean:.2f} ± {std:.2f}."

    summary.append((key, mean, std, ok, comment))
    print(f"{key}: mean={mean:.2f}, std={std:.2f}, ok={ok}")
    print("Comment:", comment)
    print("-"*50)

# Final recommendation
all_ok = all([s[3] for s in summary])
if all_ok:
    final_comment = "NPT equilibration looks stable. Safe to proceed to production run."
else:
    final_comment = ("NPT equilibration shows instabilities. "
                     "Check comments above and consider longer equilibration, "
                     "adjusting thermostat/barostat or timestep.")

print("="*60)
print("FINAL RECOMMENDATION:")
print(final_comment)
print("="*60)

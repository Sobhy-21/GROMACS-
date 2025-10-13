# analyze_md_detailed.py
# Works inside JupyterLab on Puhti

import numpy as np
import matplotlib.pyplot as plt
import os

%matplotlib inline

# === Helper to read .xvg ===
def read_xvg(filename):
    data = []
    with open(filename) as f:
        for line in f:
            if line.startswith(('#', '@')) or len(line.strip()) == 0:
                continue
            parts = line.split()
            data.append([float(x) for x in parts[:2]])
    return np.array(data)

# === Files to process ===
files = {
    "Temperature": ("run_temp.xvg", "Temperature (K)", "md_temp.png"),
    "Pressure": ("run_press.xvg", "Pressure (bar)", "md_press.png"),
    "Potential": ("run_poten.xvg", "Potential Energy (kJ/mol)", "md_poten.png"),
}

# === Check existence ===
missing = [f for f, _, _ in files.values() if not os.path.exists(f)]
if missing:
    print(f"Missing files: {', '.join(missing)}")
    print("Run the Bash extraction script first (extract_props.sh).")
else:
    # Stores overall stability flags and detailed comments
    stability_flags = {}
    detailed_comments = {}

    for key, (fname, ylabel, outpng) in files.items():
        data = read_xvg(fname)
        time = data[:, 0]
        values = data[:, 1]

        mean = np.mean(values)
        std = np.std(values)
        drift = values[-1] - values[0]

        # --- Stability evaluation ---
        if key == "Temperature":
            if 290 <= mean <= 310 and abs(drift) < 5:
                comment_short = "Temperature is stable and well-controlled."
                comment_detail = (
                    "The system maintained a steady temperature near 300 K, "
                    "indicating the thermostat is working properly. "
                    "This suggests the system has equilibrated thermally."
                )
                stable = True
            else:
                comment_short = "Temperature instability detected."
                comment_detail = (
                    "The temperature shows noticeable drift or deviation from the target (300 K). "
                    "Consider extending the NVT equilibration, checking the thermostat coupling constant, "
                    "or ensuring the initial minimization was sufficient."
                )
                stable = False

        elif key == "Pressure":
            if abs(mean) < 50 or (0 <= mean <= 10):
                comment_short = "Pressure is within the acceptable range."
                comment_detail = (
                    "Pressure fluctuates within a reasonable range, which is normal for MD simulations. "
                    "The barostat appears to be functioning correctly."
                )
                stable = True
            else:
                comment_short = "Pressure shows high fluctuations."
                comment_detail = (
                    "Pressure deviates significantly from the expected range. "
                    "This can happen if equilibration is too short or the barostat coupling is too tight. "
                    "Consider a longer NPT equilibration or adjusting the barostat parameters."
                )
                stable = False

        elif key == "Potential":
            rel_drift = abs(drift / abs(mean))
            if rel_drift < 0.001:
                comment_short = "Potential energy has plateaued."
                comment_detail = (
                    "The potential energy reached a stable plateau, indicating structural stability "
                    "and that the system is equilibrated energetically."
                )
                stable = True
            else:
                comment_short = "Potential energy still drifting."
                comment_detail = (
                    "The potential energy continues to change over time, suggesting the system "
                    "has not reached full equilibration. Extend the NVT or NPT run before production."
                )
                stable = False

        stability_flags[key] = stable
        detailed_comments[key] = (comment_short, comment_detail)

        # --- Plot ---
        plt.figure(figsize=(8, 4))
        plt.plot(time, values, lw=1.5)
        plt.title(f"{key} vs Time", fontsize=14)
        plt.xlabel("Time (ps)")
        plt.ylabel(ylabel)
        plt.grid(True, alpha=0.3)
        plt.figtext(0.5, -0.08, f"{comment_short}\n{comment_detail}", ha="center", fontsize=10, wrap=True)
        plt.tight_layout()
        plt.savefig(outpng, dpi=300, bbox_inches="tight")
        plt.show()

    # --- Overall analysis summary ---
    print("\n--- Overall MD Stability Summary ---")

    if all(stability_flags.values()):
        overall_comment = (
            "The simulation appears stable across all metrics (temperature, pressure, potential energy). "
            "You can confidently proceed with production MD or analysis. "
            "Still, you may verify structural RMSD/RMSF to confirm structural stability."
        )
    elif sum(stability_flags.values()) >= 2:
        overall_comment = (
            "The system is mostly stable, but one property (e.g., pressure or potential energy) "
            "shows signs of minor instability. "
            "Consider a short additional equilibration (~100–200 ps) before starting production runs."
        )
    else:
        overall_comment = (
            "Multiple properties show significant instability. "
            "It’s recommended to repeat or extend the equilibration phases (energy minimization, NVT, NPT). "
            "Re-examine thermostat and barostat coupling parameters, and ensure the system is well-solvated."
        )

    print(overall_comment)

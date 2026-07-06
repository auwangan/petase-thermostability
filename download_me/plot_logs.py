import glob, os
import pandas as pd
import matplotlib.pyplot as plt

logs = sorted(glob.glob("*.log"))
data = {}
for log in logs:
    name = log.replace(".log", "")
    rows = []
    with open(log) as f:
        for line in f:
            line = line.strip()
            # progress data lines start with a number and have commas
            if line and line[0].isdigit() and "," in line:
                parts = line.split(",")
                if len(parts) >= 5:
                    try:
                        step = int(parts[1])
                        pe   = float(parts[2])
                        temp = float(parts[3])
                        rows.append((step, pe, temp))
                    except ValueError:
                        continue
    if rows:
        df = pd.DataFrame(rows, columns=["step", "pe", "temp"])
        df["ns"] = df["step"] * 0.002 / 1000  # steps -> ns
        data[name] = df

# --- Plot 1: potential energy, WT vs mutant at 400K ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for name, df in data.items():
    if "400K" in name:
        color = "red" if name.startswith("mut") else "blue"
        axes[0].plot(df["ns"], df["pe"], color=color, alpha=0.6,
                     label=name if "rep1" in name else None)
axes[0].set_title("Potential Energy at 400K (blue=WT, red=mutant)")
axes[0].set_xlabel("Time (ns)"); axes[0].set_ylabel("PE (kJ/mol)")
axes[0].legend()

# --- Plot 2: potential energy, WT vs mutant at 300K ---
for name, df in data.items():
    if "300K" in name:
        color = "red" if name.startswith("mut") else "blue"
        axes[1].plot(df["ns"], df["pe"], color=color, alpha=0.6,
                     label=name if "rep1" in name else None)
axes[1].set_title("Potential Energy at 300K (blue=WT, red=mutant)")
axes[1].set_xlabel("Time (ns)"); axes[1].set_ylabel("PE (kJ/mol)")
axes[1].legend()

plt.tight_layout()
plt.savefig("energy_comparison.png", dpi=120)
print("saved energy_comparison.png")

# also print summary stats
print("\nMean potential energy (last half of each run):")
for name, df in sorted(data.items()):
    tail = df[df["ns"] > df["ns"].max()/2]
    print(f"  {name}: PE {tail['pe'].mean():.0f} kJ/mol, T {tail['temp'].mean():.1f} K")
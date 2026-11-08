# Step 4: intervals between consecutive clicks of the same detector -> dead time, afterpulses
import numpy as np
import matplotlib.pyplot as plt
from common import load_data

det, T_tot = load_data()

fig, axs = plt.subplots(1, 3, figsize=(14, 4))
for ax, c, name in zip(axs, (1, 2, 3), ("ch1 transmitted", "ch2 reflected", "ch3 herald")):
    dt = np.diff(det[c])
    rate = len(det[c]) / T_tot
    bw = 10e-9
    bins = np.arange(0, 2e-6, bw)
    h, e = np.histogram(dt, bins=bins)
    x = 0.5 * (e[1:] + e[:-1])
    expo = len(dt) * rate * np.exp(-rate * x) * bw     # exponential expected for a Poisson process

    print(f"{name}: dead time (min interval) = {dt.min()*1e9:.1f} ns, "
          f"intervals < 1 us: {(dt<1e-6).mean()*100:.2f}% (Poisson expectation {100*(1-np.exp(-rate*1e-6)):.2f}%)")

    ax.step(x * 1e9, h, where="mid", label="data")
    ax.plot(x * 1e9, expo, "r--", label=f"exponential, {rate:.0f} Hz")
    ax.set_yscale("log"); ax.set_ylim(0.5, h.max() * 3)
    ax.set_xlabel(r"$\Delta t$ [ns]"); ax.set_ylabel("counts / 10 ns")
    ax.set_title(name); ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig("step4_interarrival.png", dpi=150)
plt.show()
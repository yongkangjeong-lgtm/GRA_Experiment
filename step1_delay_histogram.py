# Step 1: herald-detector delay histograms + Gaussian fits
import numpy as np
import matplotlib.pyplot as plt
from common import load_data, delays, fit_peak, gauss

det, T_tot = load_data()
print(f"Acquisition time {T_tot:.2f} s")
for c, name in ((1, "transmitted"), (2, "reflected"), (3, "herald")):
    print(f"ch{c} {name}: {len(det[c])} counts, {len(det[c])/T_tot:.0f} Hz")

fig, axs = plt.subplots(1, 2, figsize=(11, 4))
for ax, c, name in ((axs[0], 1, "transmitted (ch1)"), (axs[1], 2, "reflected (ch2)")):
    dt = delays(det[3], det[c])
    popt, perr, x, h = fit_peak(dt)
    A, mu, sig, B = popt
    print(f"{name}: mu = {mu*1e9:.3f} +- {perr[1]*1e9:.3f} ns, "
          f"sigma = {sig*1e9:.3f} +- {perr[2]*1e9:.3f} ns, background = {B:.2f}/bin")

    ax.step(x * 1e9, h, where="mid", label="data")
    xx = np.linspace(x[0], x[-1], 2000)
    ax.plot(xx * 1e9, gauss(xx, *popt), "r-",
            label=rf"fit: $\mu$={mu*1e9:.2f} ns, $\sigma$={sig*1e9:.2f} ns")
    ax.set_yscale("log"); ax.set_ylim(0.5, h.max() * 2)
    ax.set_xlabel(r"$t_{det} - t_{herald}$ [ns]"); ax.set_ylabel("counts / 0.2 ns")
    ax.set_title("herald - " + name); ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig("step1_delay_histogram.png", dpi=150)
plt.show()
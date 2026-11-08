# Step 3: alpha as a function of the window half-width w
import numpy as np
import matplotlib.pyplot as plt
from common import load_data, get_peaks, count_coinc, alpha_from

det, _ = load_data()
muT, sT, muR, sR = get_peaks(det)

ws = np.logspace(-10, -5, 60)          # 0.1 ns ~ 10 us
al, dal = [], []
for w in ws:
    a, da = alpha_from(*count_coinc(det, (muT - w, muT + w), (muR - w, muR + w)))
    al.append(a); dal.append(da)

for w, a, da in zip(ws[::6], al[::6], dal[::6]):
    print(f"w = {w*1e9:9.1f} ns : alpha = {a:.3f} +- {da:.3f}")

plt.figure(figsize=(7, 4.2))
plt.errorbar(ws * 1e9, al, dal, fmt="o-", ms=3, capsize=2, label=r"$\alpha(w)$")
plt.axhline(1, color="r", ls="--", label=r"classical limit $\alpha=1$")
plt.axvline(3 * max(sT, sR) * 1e9, color="gray", ls=":", label=r"$3\sigma$")
plt.xscale("log"); plt.xlabel("window half-width $w$ [ns]"); plt.ylabel(r"$\alpha$")
plt.legend(); plt.grid(alpha=.3); plt.tight_layout()
plt.savefig("step3_alpha_vs_window.png", dpi=150)
plt.show()
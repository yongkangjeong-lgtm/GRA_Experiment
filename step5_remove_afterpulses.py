# Step 5: remove afterpulses (threshold 300 ns, chosen from the step 4 plot) and repeat steps 2-3
import numpy as np
import matplotlib.pyplot as plt
from common import load_data, get_peaks, count_coinc, alpha_from, remove_afterpulses

det, _ = load_data()
clean = remove_afterpulses(det, thr=300e-9)

for label, d in (("raw", det), ("afterpulses removed", clean)):
    muT, sT, muR, sR = get_peaks(d)
    N1, Nt, Nr, Nc = count_coinc(d, (muT - 3*sT, muT + 3*sT), (muR - 3*sR, muR + 3*sR))
    a, da = alpha_from(N1, Nt, Nr, Nc)
    print(f"[{label}] N1={N1}, Nt={Nt}, Nr={Nr}, Nc={Nc}, alpha = {a:.3f} +- {da:.3f}", end="")
    print(f", upper limit alpha < {3.0*N1/(Nt*Nr):.3f}" if Nc == 0 else "")

# alpha(w) comparison
ws = np.logspace(-10, -5, 60)
plt.figure(figsize=(7, 4.2))
for label, d, m in (("raw", det, "o"), ("afterpulses removed", clean, "s")):
    muT, sT, muR, sR = get_peaks(d)
    al, dal = zip(*[alpha_from(*count_coinc(d, (muT - w, muT + w), (muR - w, muR + w))) for w in ws])
    plt.errorbar(ws * 1e9, al, dal, fmt=m + "-", ms=3, capsize=2, label=label, alpha=.8)
plt.axhline(1, color="r", ls="--", label=r"$\alpha=1$")
plt.xscale("log"); plt.xlabel("window half-width $w$ [ns]"); plt.ylabel(r"$\alpha$")
plt.legend(); plt.grid(alpha=.3); plt.tight_layout()
plt.savefig("step5_alpha_comparison.png", dpi=150)
plt.show()
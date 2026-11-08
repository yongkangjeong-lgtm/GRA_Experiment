# Step 2: count N1, Nt, Nr, Nc in +-3 sigma windows and compute alpha
from common import load_data, get_peaks, count_coinc, alpha_from

det, _ = load_data()
muT, sT, muR, sR = get_peaks(det)

N1, Nt, Nr, Nc = count_coinc(det, (muT - 3*sT, muT + 3*sT), (muR - 3*sR, muR + 3*sR))
pt, pr, pc = Nt / N1, Nr / N1, Nc / N1
a, da = alpha_from(N1, Nt, Nr, Nc)

print(f"N1 = {N1}, Nt = {Nt}, Nr = {Nr}, Nc = {Nc}")
print(f"p_t = {pt:.4e}, p_r = {pr:.4e}, p_c = {pc:.4e}")
print(f"alpha = {a:.3f} +- {da:.3f}")
print(f"Expected Nc if independent (classical) = N1 p_t p_r = {N1*pt*pr:.1f}")
if Nc == 0:
    # 0 observed -> Poisson mean < 3.0 at 95% CL
    print(f"Nc = 0 -> 95% CL upper limit: alpha < {3.0*N1/(Nt*Nr):.3f}")
# common.py : shared functions used by all steps (do not run directly)
import numpy as np
from scipy.optimize import curve_fit

UNIT = 80.955e-12            # 1 time-tag unit = 80.955 ps
FILE = "TimeTags_2.txt"


def load_data(file=FILE):
    """Read the file and return click times (in s) per channel: det[1], det[2], det[3]"""
    raw = np.loadtxt(file, delimiter=";", comments="#")
    t = raw[:, 0] * UNIT
    ch = raw[:, 1].astype(int)
    det = {c: t[ch == c] for c in (1, 2, 3)}
    T_tot = t[-1] - t[0]
    return det, T_tot


def delays(t_her, t_det, win=50e-9):
    """For every herald click, collect all time differences with detector clicks within +-win"""
    lo = np.searchsorted(t_det, t_her - win)
    hi = np.searchsorted(t_det, t_her + win)
    idx_h = np.repeat(np.arange(len(t_her)), hi - lo)
    idx_d = np.concatenate([np.arange(a, b) for a, b in zip(lo, hi)])
    return t_det[idx_d] - t_her[idx_h]


def gauss(x, A, mu, sig, B):
    return A * np.exp(-0.5 * ((x - mu) / sig) ** 2) + B


def fit_peak(dt, win=20e-9, bw=0.2e-9):
    """Fit the delay histogram with a Gaussian plus constant background"""
    bins = np.arange(-win, win + bw, bw)
    h, e = np.histogram(dt, bins=bins)
    x = 0.5 * (e[1:] + e[:-1])
    p0 = [h.max(), x[np.argmax(h)], 0.5e-9, np.median(h)]
    popt, pcov = curve_fit(gauss, x, h, p0=p0, sigma=np.sqrt(h + 1), absolute_sigma=True)
    perr = np.sqrt(np.diag(pcov))
    popt[2] = abs(popt[2])
    return popt, perr, x, h


def get_peaks(det):
    """mu and sigma of the transmitted (ch1) and reflected (ch2) peaks"""
    pT, _, _, _ = fit_peak(delays(det[3], det[1]))
    pR, _, _, _ = fit_peak(delays(det[3], det[2]))
    return pT[1], pT[2], pR[1], pR[2]


def count_coinc(det, win_T, win_R):
    """win = (start, end) relative to the herald. Returns N1, Nt, Nr, Nc"""
    h = det[3]

    def hits(td, w):
        lo = np.searchsorted(td, h + w[0])
        hi = np.searchsorted(td, h + w[1])
        return hi > lo                       # True/False for each herald

    ct = hits(det[1], win_T)
    cr = hits(det[2], win_R)
    return len(h), ct.sum(), cr.sum(), (ct & cr).sum()


def alpha_from(N1, Nt, Nr, Nc):
    """alpha = N1 Nc / (Nt Nr), with Poisson error propagation"""
    a = N1 * Nc / (Nt * Nr)
    rel2 = (1 / Nc if Nc > 0 else 0) + 1 / Nt + 1 / Nr + 1 / N1
    return a, a * np.sqrt(rel2)


def remove_afterpulses(det, thr=300e-9):
    """Remove clicks arriving within thr after the previous click of the same channel"""
    out = {}
    for c in (1, 2, 3):
        t = det[c]
        keep = np.ones(len(t), bool)
        last = -np.inf
        for i, ti in enumerate(t):
            if ti - last < thr:
                keep[i] = False
            else:
                last = ti
        out[c] = t[keep]
        print(f"ch{c}: {(~keep).sum()} / {len(t)} removed ({100*(~keep).mean():.2f}%)")
    return out
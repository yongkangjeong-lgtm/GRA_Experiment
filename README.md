# Grangier–Roger–Aspect Experiment – Data Analysis

Analysis of heralded single-photon anticorrelation data (Quantum Optics and Lasers).

## Files
- `common.py` – shared functions (loading, fitting, coincidence counting)
- `step1_delay_histogram.py` – herald–detector delay histograms and Gaussian fits
- `step2_alpha.py` – coincidences in ±3σ windows and anticorrelation parameter α
- `step3_alpha_vs_window.py` – α as a function of the window width
- `step4_interarrival.py` – time between consecutive clicks (dead time, afterpulses)
- `step5_remove_afterpulses.py` – afterpulse removal and comparison

## How to run
Put `TimeTags_2.txt` in the same folder, then:
    pip install numpy scipy matplotlib
    python step1_delay_histogram.py

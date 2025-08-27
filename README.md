# Theta-Chirp Detector

A fast detector for periodic **theta-chirp** structures in 1-D signals. It combines:

* **FrFT-swept gamma search** (fractional Fourier transform) for efficient chirp-rate estimation
* **AR(1) prewhitening** (Yule–Walker) for colored-noise suppression
* **Jacobi-comb multi-harmonic pooling** across harmonic centers
* **Per-q β–γ heatmaps** for quick debugging and model inspection
* **Runtime benchmark harness** to compare brute gamma grids vs. FrFT sweeps

> The design mirrors an asymptotic/Rademacher viewpoint: capture the dominant term accurately and control the tail exponentially—implemented here via compact atoms, FrFT acceleration, and exponential truncation of harmonics.

---

##  Features

* **Theta-chirp atoms** with Gaussian harmonic window and quadratic phase (`beta`, `gamma`)
  – `theta_chirp_atom`, `theta_chirp_atom_centered_m`
* **Brute-force grid** over `(q, β, γ)` or **FrFT-swept** gamma search
  – `periodogram_grid`, `theta_chirp_periodogram_frft`
* **Pooling across harmonics** (`none` / `avg` with weights `m^{-p}` / `max`)
* **AR(1) prewhitening** (`--ar1`) using Yule–Walker (`prewhiten_ar1`)
* **β–γ heatmaps per q** for diagnostics (`beta_gamma_heatmap_for_q`)
* **Benchmarks & plots**: periodogram, heatmaps, runtime curves (`runtime_benchmark`)

---

##  Quickstart (CLI)

Run the detector on the built-in demo signal (a noisy theta-chirp around `q≈32`, `γ≈0.012`) and save plots:

```bash
python theta_chirp_detector.py --frft --ar1 --save-plots --verbose
```

You’ll see a summary like:

```
Top q=32 with (beta,gamma)=(0.6,0.012)
```

Saved figures:

* `periodogram_q.png` – best energy vs. period q
* `heatmap_q<top>.png` – β–γ energy map for the top q
* (if `--bench --save-plots`) `runtime_scaling.png` – runtime vs. N

---

##  Command-Line Options

```text
--seed <int>                 RNG seed (default 42)
--Qmax <int>                 Max period to test (default 64)
--betas "0.3,0.6,1.0"        Comma list of β values
--gammas "-0.02,-0.01,0,0.01,0.02"  γ grid (brute mode)
--alphas "-0.20,...,0.20"    FrFT orders (FrFT mode)
--M <int>                    Harmonic half-width (default 5)
--zero-mean                  Zero-mean the input
--frft                       Use FrFT sweep for γ (else brute γ grid)
--pool {none,avg,max}        Harmonic pooling mode
--mset "1,2,3"               Harmonic centers for pooling
--p <float>                  Weights m^{-p} if pool=avg (default 1.0)
--ar1                        AR(1) prewhitening (Yule–Walker)
--save-plots                 Save periodogram/heatmaps/bench plots
--heatmap-topk <int>         Save β–γ maps for top-K q (default 1)
--bench                      Run runtime benchmarks
--verbose                    Print details
```

**Notes**

* In brute mode, γ is taken from `--gammas`.
* In FrFT mode, γ is **mapped** from best FrFT order and then **locally refined** around that value (`refine_gamma_span`, `refine_steps` inside the function).

---

##  Examples

**1) Brute γ grid with multi-harmonic pooling**

```bash
python theta_chirp_detector.py \
  --pool avg --mset 1,2,3,4 --p 0.5 \
  --gammas -0.03,-0.02,-0.01,0,0.01,0.02,0.03 \
  --save-plots --verbose
```

**2) FrFT-swept γ (fast) with AR(1) prewhitening**

```bash
python theta_chirp_detector.py --frft --ar1 --save-plots --verbose
```

**3) Benchmarks: brute vs. FrFT**

```bash
python theta_chirp_detector.py --bench --frft --save-plots
```

---

##  Library Usage (API)

```python
import numpy as np
from theta_chirp_detector import (
    periodogram_grid, theta_chirp_periodogram_frft,
    beta_gamma_heatmap_for_q, prewhiten_ar1
)

x = ...  # your complex/real 1-D signal (numpy array)
x_pw, rho = prewhiten_ar1(x)  # optional AR(1) prewhitening

# 1) Brute grid over (q, β, γ)
S, picks = periodogram_grid(x_pw, Qmax=64, betas=(0.3,0.6,1.0),
                            gammas=(-0.02,-0.01,0,0.01,0.02), M=5)
q_star = 1 + np.argmax(S.max(axis=(1,2)))

# 2) FrFT-swept gamma (faster)
best, picks, scores = theta_chirp_periodogram_frft(x_pw, Qmax=64,
                      betas=(0.3,0.6,1.0),
                      alphas=(-0.20,-0.12,-0.06,0,0.06,0.12,0.20),
                      M=5)
q_star = 1 + np.argmax(best)
beta_star, gamma_star = picks[q_star - 1]
```

API functions are defined in `theta_chirp_detector.py` and mirror the CLI behavior.

---

##  Method (Short)

* **Atoms.** A theta-chirp atom is a Gaussian-weighted *Jacobi comb* in frequency, modulated by a quadratic phase `exp(i π γ t²)`; we subtract the mean and L2-normalize.
* **Search.** We scan over periods `q`, β (harmonic window), and γ (chirp rate).
  – Brute grid evaluates γ directly.
  – **FrFT sweep** maps fractional orders → candidate γ and refines locally.
* **Pooling.** To be robust to harmonic content, we pool energies across selected harmonic centers `m0 ∈ mset` using `max` or weighted `avg` with `m^{-p}`.
* **Noise.** Optional **AR(1) prewhitening** removes colored-noise bias before detection.
* **Diagnostics.** β–γ heatmaps visualize how energy concentrates for each `q`.

This strategy is inspired by asymptotic expansions (dominant term + exponentially small tail), letting us **truncate harmonics** and **accelerate γ** without losing detection power.

---

##  Benchmarks

The built-in benchmark compares **brute γ grid** vs. **FrFT sweep** for `N ∈ {256,512,1024}` and reports runtimes (optionally saving `runtime_scaling.png`). Expect FrFT to **win** as N grows, since it reduces the γ search to a few FrFT orders plus a short local refinement.

---

##  Tips

* **Gamma range.** If your true chirp rate is outside the tested band, expand `--gammas` (brute) or `--alphas` / refinement span (FrFT). The mapping from FrFT order to γ (`frft_gamma_map`) is intentionally simple; adjust `kappa` inside the function if needed.
* **Harmonic span.** `--M` controls the half-width in harmonic index; for narrow-band signals, smaller `M` is fine. The internal `auto_M(beta)` is used when `M=None` and keeps the Gaussian tail below a tiny ε.
* **Zero-mean.** Use `--zero-mean` if your signal has a DC offset.
* **Real vs. complex.** Inputs can be real or complex; atoms are complex.

---

##  Roadmap

* File I/O (load signal from `.npy` / `.wav` / `.csv`)
* Batched / streaming inference
* Adaptive FrFT-γ mapping and learning-based refinement
* Confidence metrics and false-alarm estimation
* GPU acceleration for large N

---

##  License

MIT License

---


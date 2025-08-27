
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
theta_chirp_detector.py
-----------------------
Theta-chirp periodogram with:
 - FrFT sweep for gamma (chirp-rate) acceleration
 - AR(1) prewhitening (Yule-Walker)
 - Jacobi-comb multi-harmonic pooling
 - Per-q beta–gamma heatmaps
 - Runtime benchmark harness

Dependencies: numpy, matplotlib (for plotting).

Author: ChatGPT
"""
import argparse
import math
import time
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Utilities
# -----------------------------

def set_seed(seed=42):
    rng = np.random.default_rng(int(seed))
    return rng

def yule_walker_ar1(x):
    """
    Estimate AR(1) coefficient rho via Yule-Walker:
      rho = Cov(x_t, x_{t-1}) / Var(x_t)
    Returns rho in (-1,1).
    """
    x = np.asarray(x, dtype=float)
    x = x - x.mean()
    num = np.dot(x[1:], x[:-1])
    den = np.dot(x, x) + 1e-12
    rho = num / den
    rho = np.clip(rho, -0.99, 0.99)
    return float(rho)

def prewhiten_ar1(x, rho=None):
    """
    Prewhiten series with AR(1) parameter rho.
    If rho is None, estimate via Yule-Walker.
    Returns y where y[t] = x[t] - rho*x[t-1], y[0]=x[0].
    """
    x = np.asarray(x, dtype=float)
    if rho is None:
        rho = yule_walker_ar1(x)
    y = np.empty_like(x, dtype=float)
    y[0] = x[0]
    y[1:] = x[1:] - rho * x[:-1]
    return y, rho

def auto_M(beta, eps=1e-10):
    """
    Choose harmonic truncation M for theta envelope weights exp(-beta m^2)
    such that tail energy < eps.
    Simple bound: exp(-beta (M+1)^2) ~ eps -> M ~ sqrt(log(1/eps)/beta).
    """
    if beta <= 0: 
        return 3
    M = int(np.ceil(np.sqrt(max(1e-16, np.log(1.0/eps))) / max(1e-8, beta)))
    return max(3, min(M, 64))

# -----------------------------
# Atoms and periodograms
# -----------------------------

def theta_chirp_atom(N, q, beta=0.5, gamma=0.0, M=5):
    """
    A theta-chirp atom:
      phi[n] = (sum_{m=-M..M} exp(-beta m^2) * exp(2πi m n / q)) * exp(iπ gamma t^2)
    where t = (n - (N-1)/2)/N. Zero-mean, L2-normalized.
    """
    n = np.arange(N)
    t = (n - (N - 1) / 2) / N
    chirp = np.exp(1j * np.pi * gamma * (t ** 2))
    m = np.arange(-M, M + 1)
    weights = np.exp(-beta * (m ** 2))
    phase = np.exp(2j * np.pi * np.outer(n, m) / q)
    theta_env = phase @ weights
    atom = theta_env * chirp
    atom = atom - atom.mean()
    norm = np.linalg.norm(atom)
    return atom / (norm + 1e-12)

def theta_chirp_atom_centered_m(N, q, beta=0.5, gamma=0.0, M=3, m0=1):
    """
    Theta-chirp atom with harmonic center m0:
      frequencies at (m0 + u)/q for u=-M..M with Gaussian weights in u.
    """
    n = np.arange(N)
    t = (n - (N - 1) / 2) / N
    chirp = np.exp(1j * np.pi * gamma * (t ** 2))
    u = np.arange(-M, M + 1)
    m = m0 + u
    weights = np.exp(-beta * (u ** 2))
    phase = np.exp(2j * np.pi * np.outer(n, m) / q)
    theta_env = phase @ weights
    atom = theta_env * chirp
    atom = atom - atom.mean()
    norm = np.linalg.norm(atom)
    return atom / (norm + 1e-12)

def periodogram_grid(x, Qmax=64, betas=(0.3, 0.6, 1.0), gammas=(-0.02, -0.01, 0.0, 0.01, 0.02),
                     M=5, zero_mean=True, pool_mode='none', mset=(1,2,3), p=1.0):
    """
    Brute-force grid over q, beta, gamma using explicit atoms.
    Returns energies S[Q, B, G] and picks[(q-1)] = (beta*, gamma*).
    """
    x = np.asarray(x, dtype=complex)
    if zero_mean:
        x = x - x.mean()
    N = x.size
    S = np.zeros((Qmax, len(betas), len(gammas)), dtype=float)
    picks = []
    for qi in range(Qmax):
        q = qi + 1
        for bi, beta in enumerate(betas):
            M_eff = auto_M(beta) if M is None else M
            for gi, gamma in enumerate(gammas):
                if pool_mode == 'none':
                    atom = theta_chirp_atom(N, q, beta, gamma, M_eff)
                    coeff = np.vdot(atom, x)
                    S[qi, bi, gi] = float(np.abs(coeff) ** 2)
                else:
                    # Jacobi-comb pooling across harmonic centers in mset
                    energies = []
                    etas = np.array([m ** (-p) for m in mset], dtype=float) if pool_mode == 'avg' else None
                    for m0 in mset:
                        atom = theta_chirp_atom_centered_m(N, q, beta, gamma, max(1, M_eff//2), m0=m0)
                        energies.append(float(np.abs(np.vdot(atom, x)) ** 2))
                    energies = np.array(energies, dtype=float)
                    if pool_mode == 'avg':
                        S[qi, bi, gi] = float(np.sum(etas * energies) / (np.sum(etas) + 1e-12))
                    else:  # 'max'
                        S[qi, bi, gi] = float(np.max(energies))
    for qi in range(Qmax):
        bi, gi = np.unravel_index(np.argmax(S[qi]), (len(betas), len(gammas)))
        picks.append((float(betas[bi]), float(gammas[gi])))
    return S, picks

# -----------------------------
# FrFT implementation and FrFT-based gamma sweep
# -----------------------------

def frft(x, a):
    """
    Discrete Fractional Fourier Transform (order a).
    a=0 -> identity; a=1 -> FFT (up to normalization/shift).
    Implementation via chirp-FFT-chirp (O(N log N)).
    """
    x = np.asarray(x, dtype=complex)
    N = x.size

    # Handle integer-like orders quickly
    if np.isclose(a, 0):
        return x.copy()
    if np.isclose(a, 1):
        return np.fft.fftshift(np.fft.fft(np.fft.ifftshift(x))) / np.sqrt(N)
    if np.isclose(a, -1):
        return np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(x))) * np.sqrt(N)

    alpha = a * np.pi / 2
    tana2 = np.tan(alpha / 2.0)
    sina = np.sin(alpha)

    # Avoid singular sin(alpha)
    if np.isclose(sina, 0.0):
        # fall back to near-identity / FFT cases
        eps = 1e-6
        sina = np.sign(sina) * eps if np.sign(sina) != 0 else eps

    n = np.arange(N) - (N - 1) / 2.0
    c = np.exp(-1j * np.pi * (n ** 2) * tana2 / N)
    X = c * x

    k = np.arange(N) - (N - 1) / 2.0
    h = np.exp(1j * np.pi * (k ** 2) / (N * sina))
    H = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(h)))

    Y = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(
        np.fft.fftshift(np.fft.fft(np.fft.ifftshift(X))) * H
    )))
    y = c * Y * np.sqrt(1j / (sina + 1e-12))
    return y

def frft_gamma_map(alpha, N, kappa=1.0):
    """
    Map FrFT order alpha to approximate chirp-rate gamma.
    Simple near-identity linearization; kappa is a scale knob.
    """
    # Small-angle approx: gamma ≈ tan(π alpha / 2) / kappa scaled by N
    return float(np.tan(0.5 * np.pi * alpha) / (kappa + 1e-12))

def theta_chirp_periodogram_frft(x, Qmax=64, betas=(0.3,0.6,1.0),
                                 alphas=(-0.20,-0.12,-0.06,0.0,0.06,0.12,0.20),
                                 M=5, zero_mean=True, refine_K=5, refine_gamma_span=0.02,
                                 refine_steps=5, kappa=1.0, pool_mode='none', mset=(1,2,3), p=1.0):
    """
    FrFT-swept approximation of gamma grid.
    - Coarse: sweep alpha orders, get per-(q,beta) best alpha
    - Fine: for top-K q, refine gamma in a small window around mapped gamma
    Returns best_over_bg[q] and picks[q]=(beta*, gamma*).
    """
    x = np.asarray(x, dtype=complex)
    if zero_mean:
        x = x - x.mean()
    N = x.size
    n = np.arange(N)
    t = (n - (N - 1) / 2) / N

    # Precompute theta envelopes for all (q, beta)
    m = np.arange(-M, M + 1)
    theta_env = {}
    for q in range(1, Qmax + 1):
        P = np.exp(2j * np.pi * np.outer(n, m) / q)
        for beta in betas:
            w = np.exp(-beta * (m ** 2))
            theta_env[(q, beta)] = (P @ w)

    # Coarse sweep over alpha
    coarse_scores = np.zeros((Qmax, len(betas)), dtype=float)
    alpha_star = np.zeros((Qmax, len(betas)), dtype=float)

    for a in alphas:
        y = frft(x, a)
        for qi in range(Qmax):
            q = qi + 1
            for bi, beta in enumerate(betas):
                phi = theta_env[(q, beta)]
                phi0 = phi - phi.mean()
                if pool_mode == 'none':
                    s = np.vdot(phi0, y) / (np.linalg.norm(phi0) + 1e-12)
                    e = (np.abs(s) ** 2)
                else:
                    # Pooling across harmonic centers
                    energies = []
                    etas = np.array([m0 ** (-p) for m0 in mset], dtype=float) if pool_mode == 'avg' else None
                    for m0 in mset:
                        atom = theta_chirp_atom_centered_m(N, q, beta, 0.0, max(1, M//2), m0=m0)
                        atom = atom - atom.mean()
                        s = np.vdot(atom, y) / (np.linalg.norm(atom) + 1e-12)
                        energies.append(float(np.abs(s) ** 2))
                    energies = np.array(energies, dtype=float)
                    if pool_mode == 'avg':
                        e = float(np.sum(etas * energies) / (np.sum(etas) + 1e-12))
                    else:
                        e = float(np.max(energies))

                if e > coarse_scores[qi, bi]:
                    coarse_scores[qi, bi] = e
                    alpha_star[qi, bi] = a

    best_beta_idx = np.argmax(coarse_scores, axis=1)
    q_ranking = np.argsort(coarse_scores[np.arange(Qmax), best_beta_idx])[::-1]
    q_top = q_ranking[:max(1, refine_K)]

    # Fine refine gamma around mapped gamma from alpha*
    gammas_refine = np.linspace(-refine_gamma_span, refine_gamma_span, max(3, refine_steps))
    best_scores = coarse_scores.copy()
    gamma_pick = np.zeros((Qmax, len(betas)), dtype=float)

    for qi in q_top:
        q = qi + 1
        for bi, beta in enumerate(betas):
            a0 = alpha_star[qi, bi]
            g0 = frft_gamma_map(a0, N, kappa=kappa)
            for dg in gammas_refine:
                g = g0 + dg
                chirp = np.exp(1j * np.pi * g * (t ** 2))
                phi = theta_env[(q, beta)] * chirp
                phi0 = phi - phi.mean()

                if pool_mode == 'none':
                    s = np.vdot(phi0, x) / (np.linalg.norm(phi0) + 1e-12)
                    e = (np.abs(s) ** 2)
                else:
                    energies = []
                    etas = np.array([m0 ** (-p) for m0 in mset], dtype=float) if pool_mode == 'avg' else None
                    for m0 in mset:
                        atom = theta_chirp_atom_centered_m(N, q, beta, g, max(1, M//2), m0=m0)
                        atom = atom - atom.mean()
                        s = np.vdot(atom, x) / (np.linalg.norm(atom) + 1e-12)
                        energies.append(float(np.abs(s) ** 2))
                    energies = np.array(energies, dtype=float)
                    if pool_mode == 'avg':
                        e = float(np.sum(etas * energies) / (np.sum(etas) + 1e-12))
                    else:
                        e = float(np.max(energies))

                if e > best_scores[qi, bi]:
                    best_scores[qi, bi] = e
                    gamma_pick[qi, bi] = g

    best_over_bg = best_scores.max(axis=1)
    picks = []
    for qi in range(Qmax):
        bi = np.argmax(best_scores[qi])
        picks.append((float(betas[bi]), float(gamma_pick[qi, bi])))
    return best_over_bg, picks, best_scores

# -----------------------------
# Heatmaps for debugging
# -----------------------------

def beta_gamma_heatmap_for_q(x, q, betas=(0.3,0.6,1.0), gammas=(-0.02,-0.01,0.0,0.01,0.02),
                             M=5, zero_mean=True, pool_mode='none', mset=(1,2,3), p=1.0):
    """
    Compute a per-q beta–gamma heatmap (energies) for quick debugging.
    """
    x = np.asarray(x, dtype=complex)
    if zero_mean:
        x = x - x.mean()
    N = x.size
    n = np.arange(N)
    t = (n - (N - 1) / 2) / N
    S = np.zeros((len(betas), len(gammas)), dtype=float)

    # Precompute theta envelope for q, each beta
    m = np.arange(-M, M + 1)
    P = np.exp(2j * np.pi * np.outer(n, m) / q)
    envs = {}
    for beta in betas:
        w = np.exp(-beta * (m ** 2))
        envs[beta] = (P @ w)

    for bi, beta in enumerate(betas):
        for gi, gamma in enumerate(gammas):
            chirp = np.exp(1j * np.pi * gamma * (t ** 2))
            phi = envs[beta] * chirp
            phi0 = phi - phi.mean()

            if pool_mode == 'none':
                s = np.vdot(phi0, x) / (np.linalg.norm(phi0) + 1e-12)
                S[bi, gi] = float(np.abs(s) ** 2)
            else:
                energies = []
                etas = np.array([m0 ** (-p) for m0 in mset], dtype=float) if pool_mode == 'avg' else None
                for m0 in mset:
                    atom = theta_chirp_atom_centered_m(N, q, beta, gamma, max(1, M//2), m0=m0)
                    atom = atom - atom.mean()
                    s = np.vdot(atom, x) / (np.linalg.norm(atom) + 1e-12)
                    energies.append(float(np.abs(s) ** 2))
                energies = np.array(energies, dtype=float)
                if pool_mode == 'avg':
                    S[bi, gi] = float(np.sum(etas * energies) / (np.sum(etas) + 1e-12))
                else:
                    S[bi, gi] = float(np.max(energies))
    return S

def save_heatmap(S, betas, gammas, q, out_png="heatmap.png", title_prefix="β–γ map at q="):
    plt.figure(figsize=(6, 3.8))
    plt.imshow(S, aspect='auto', origin='lower', 
               extent=[0, len(gammas)-1, 0, len(betas)-1])
    plt.colorbar(label="Energy")
    plt.yticks(range(len(betas)), [f"{b:g}" for b in betas])
    plt.xticks(range(len(gammas)), [f"{g:g}" for g in gammas], rotation=0)
    plt.xlabel("γ index")
    plt.ylabel("β index")
    plt.title(f"{title_prefix}{q}")
    plt.tight_layout()
    plt.savefig(out_png, dpi=140)
    plt.close()

# -----------------------------
# Benchmarks
# -----------------------------

def runtime_benchmark(seed=1, Ns=(256, 512, 1024), Qmax=64,
                      betas=(0.3,0.6,1.0), gammas=(-0.02,-0.01,0.0,0.01,0.02),
                      alphas=(-0.20,-0.12,-0.06,0.0,0.06,0.12,0.20), M=5,
                      pool_mode='none', save_plots=False):
    """
    Compare runtime of brute gamma grid vs FrFT sweep.
    Generates a synthetic chirp at q=32, gamma=0.012.
    """
    rng = set_seed(seed)
    t_frft = []
    t_grid = []
    for N in Ns:
        n = np.arange(N)
        t = (n - (N - 1) / 2) / N
        sig = np.cos(2*np.pi*n/32 + np.pi*0.012*(t**2))
        x = sig + 0.6 * rng.standard_normal(N)

        # Grid
        t0 = time.time()
        S, _ = periodogram_grid(x, Qmax=Qmax, betas=betas, gammas=gammas, M=M, pool_mode=pool_mode)
        t1 = time.time()
        t_grid.append(t1 - t0)

        # FrFT
        t0 = time.time()
        best, picks, _ = theta_chirp_periodogram_frft(x, Qmax=Qmax, betas=betas, alphas=alphas, M=M, pool_mode=pool_mode)
        t1 = time.time()
        t_frft.append(t1 - t0)

    if save_plots:
        plt.figure(figsize=(6,4))
        plt.plot(Ns, t_grid, 'o-', label='Gamma grid')
        plt.plot(Ns, t_frft, 's--', label='FrFT sweep')
        plt.xlabel("N")
        plt.ylabel("Runtime (s)")
        plt.title("Runtime scaling")
        plt.legend()
        plt.tight_layout()
        plt.savefig("runtime_scaling.png", dpi=140)
        plt.close()

    return {'N': list(Ns), 'grid_s': t_grid, 'frft_s': t_frft}

# -----------------------------
# CLI
# -----------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Theta-Chirp Detector with FrFT sweep, AR(1) prewhitening, pooling, heatmaps, and benchmark.")
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--Qmax', type=int, default=64)
    p.add_argument('--betas', type=str, default="0.3,0.6,1.0", help="Comma-separated list")
    p.add_argument('--gammas', type=str, default="-0.02,-0.01,0.0,0.01,0.02", help="Comma-separated list")
    p.add_argument('--alphas', type=str, default="-0.20,-0.12,-0.06,0.0,0.06,0.12,0.20", help="Comma-separated list for FrFT")
    p.add_argument('--M', type=int, default=5)
    p.add_argument('--zero-mean', action='store_true', help='Zero-mean the input')
    p.add_argument('--frft', action='store_true', help='Use FrFT sweep for gamma')
    p.add_argument('--pool', default='none', choices=['none', 'avg', 'max'])
    p.add_argument('--mset', type=str, default="1,2,3", help='Harmonic centers for pooling')
    p.add_argument('--p', type=float, default=1.0, help='Pooling weights m^{-p} for avg')
    p.add_argument('--ar1', action='store_true', help='Apply AR(1) prewhitening')
    p.add_argument('--save-plots', action='store_true', help='Save periodogram and heatmap plots')
    p.add_argument('--heatmap-topk', type=int, default=1, help='Save heatmaps for top-K q')
    p.add_argument('--bench', action='store_true', help='Run runtime benchmark')
    p.add_argument('--verbose', action='store_true', help='Print details')
    return p.parse_args()

def main():
    args = parse_args()
    rng = set_seed(args.seed)

    betas = tuple(float(x) for x in args.betas.split(','))
    gammas = tuple(float(x) for x in args.gammas.split(','))
    alphas = tuple(float(x) for x in args.alphas.split(','))
    mset = tuple(int(x) for x in args.mset.split(','))

    # Demo signal (user can replace with I/O later)
    N = 256
    n = np.arange(N)
    t = (n - (N - 1) / 2) / N
    signal = np.cos(2 * np.pi * n / 32 + np.pi * 0.012 * (t ** 2))
    noise = 0.6 * rng.standard_normal(N)
    x = signal + noise

    if args.ar1:
        x, rho = prewhiten_ar1(x)
        if args.verbose:
            print(f"AR(1) prewhitening applied: rho={rho:.3f}")

    # Periodogram
    if args.frft:
        if args.verbose: print("Running FrFT-swept detector...")
        best, picks, best_scores = theta_chirp_periodogram_frft(
            x, Qmax=args.Qmax, betas=betas, alphas=alphas, M=args.M,
            zero_mean=args.zero_mean, refine_K=5, refine_gamma_span=0.02, refine_steps=7,
            kappa=1.0, pool_mode=args.pool, mset=mset, p=args.p
        )
    else:
        if args.verbose: print("Running brute gamma-grid detector...")
        S, picks = periodogram_grid(
            x, Qmax=args.Qmax, betas=betas, gammas=gammas, M=args.M,
            zero_mean=args.zero_mean, pool_mode=args.pool, mset=mset, p=args.p
        )
        best = S.max(axis=(1,2))

    q_star = int(np.argmax(best) + 1)
    beta_star, gamma_star = picks[q_star - 1]
    print(f"Top q={q_star} with (beta,gamma)=({beta_star:.3g},{gamma_star:.3g})")

    # Save plots
    if args.save_plots:
        # Periodogram over q (best over beta,gamma)
        plt.figure(figsize=(7,3.4))
        plt.plot(np.arange(1, args.Qmax + 1), best, 'o-')
        plt.xlabel("Period q")
        plt.ylabel("Best energy")
        plt.title("Theta-chirp periodogram (best over β,γ)")
        plt.tight_layout()
        plt.savefig("periodogram_q.png", dpi=150)
        plt.close()

        # Heatmaps for top-K q
        top_idx = list(np.argsort(best)[::-1][:max(1, args.heatmap_topk)])
        for idx in top_idx:
            q = int(idx + 1)
            H = beta_gamma_heatmap_for_q(
                x, q, betas=betas, gammas=gammas, M=args.M,
                zero_mean=args.zero_mean, pool_mode=args.pool, mset=mset, p=args.p
            )
            save_heatmap(H, betas, gammas, q, out_png=f"heatmap_q{q}.png")

    # Benchmark
    if args.bench:
        if args.verbose: print("Running runtime benchmark...")
        res = runtime_benchmark(seed=args.seed, Ns=(256,512,1024), Qmax=args.Qmax,
                                betas=betas, gammas=gammas, alphas=alphas, M=args.M,
                                pool_mode=args.pool, save_plots=args.save_plots)
        print("Benchmark (s):")
        for Nval, tg, tf in zip(res['N'], res['grid_s'], res['frft_s']):
            print(f"  N={Nval}: grid={tg:.3f}s, frft={tf:.3f}s")

if __name__ == "__main__":
    main()

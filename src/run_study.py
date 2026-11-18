"""Runs the timing study and writes figures/ and results.json."""
import json, itertools, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cfd_core import *
from matplotlib.ticker import ScalarFormatter, NullFormatter
from matplotlib.colors import LogNorm

plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3, "figure.dpi": 130})
OUT = "../figures/"
N = 4000
FR = [k / 8 for k in range(1, 8)]                 # CFD fractions K/8 (hardware-friendly)
LE_THR = [0.03, 0.05, 0.08, 0.12, 0.18, 0.25]     # leading-edge thresholds (fraction of Amax)
D_NS = [1, 2, 3, 4, 6, 8, 12]

def sweep(fs, snr, seed):
    """Evaluate every LE and CFD configuration on the same events. Returns dicts of (sigma_ps, eff)."""
    rng = np.random.default_rng(seed)
    x, t0, A = make_events(N, fs, snr, rng)
    ts, sig = 1 / fs, 1 / snr
    arm = max(6 * sig, 0.05)
    le, cf = {}, {}
    for thr in LE_THR:
        if thr < 3.5 * sig:
            continue
        t, ok = leading_edge(x, thr, ts)
        le[thr] = stats(t, ok, t0) + (t, ok)
    dset = sorted({max(1, int(round(d * fs))) for d in D_NS})
    for f, d in itertools.product(FR, dset):
        t, ok = cfd(x, f, d, arm, ts)
        cf[(f, d)] = stats(t, ok, t0) + (t, ok)
    return le, cf, (x, t0, A)

def best(cfgs, min_eff=0.99):
    cand = {k: v for k, v in cfgs.items() if v[1] >= min_eff and not np.isnan(v[0])}
    k = min(cand, key=lambda k: cand[k][0])
    return k, cand[k]

def reeval(fs, snr, kind, cfg, seed):
    """Evaluate one chosen configuration on fresh (held-out) events."""
    rng = np.random.default_rng(seed)
    x, t0, A = make_events(N, fs, snr, rng)
    ts, sig = 1 / fs, 1 / snr
    if kind == "le":
        t, ok = leading_edge(x, cfg, ts)
    else:
        t, ok = cfd(x, cfg[0], cfg[1], max(6 * sig, 0.05), ts)
    return stats(t, ok, t0)[0]

res = {}

# ---- Fig 0: illustration of the two methods on one noiseless pulse -----------------------------
fs0 = 4.0
tk = np.arange(int(T_WIN * fs0)) / fs0
t0_ = 40.0
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
for A_, col in [(0.35, "tab:blue"), (1.0, "tab:red")]:
    s = pulse(tk - t0_, A_)
    d = int(round(3.0 * fs0)); f = 0.25
    sd = np.zeros_like(s); sd[d:] = s[:-d]
    y = f * s - sd
    ax[0].plot(tk, s, color=col, label=f"A = {A_}")
    ax[1].plot(tk, y, color=col)
    i = np.argmax((y[:-1] > 0) & (y[1:] <= 0) & (tk[1:] > t0_))
    zc = tk[i] + y[i] / (y[i] - y[i + 1]) / fs0
    ax[1].plot(zc, 0, "o", color=col)
    ax[0].axvline(zc, color=col, ls=":", lw=0.8)
ax[0].axhline(0.1, color="k", ls="--", lw=0.8, label="fixed threshold 0.1")
ax[0].set_xlim(35, 80); ax[1].set_xlim(35, 80)
ax[0].set_title("Input pulses (two amplitudes)"); ax[1].set_title("CFD signal  $y=f\\,x(t)-x(t-d)$, f=0.25, d=3 ns")
ax[0].set_xlabel("t [ns]"); ax[1].set_xlabel("t [ns]"); ax[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(OUT + "fig0_principle.png"); plt.close(fig)

# ---- E1: time walk vs amplitude (high SNR, fine sampling) ----------------------------------------
fs, snr = 4.0, 1000
le, cf, (x, t0, A) = sweep(fs, snr, 11)
kc, vc = best(cf)
thr_ref = 0.12
tl, okl = le[thr_ref][2], le[thr_ref][3]
tc, okc = cf[kc][2], cf[kc][3]
bins = np.linspace(A_MIN, 1, 12)
cen = 0.5 * (bins[1:] + bins[:-1])
def walk(t, ok):
    e = (t - t0)
    m = np.array([e[(A >= lo) & (A < hi) & ok].mean() for lo, hi in zip(bins[:-1], bins[1:])])
    return 1000 * (m - np.nanmean(e[ok]))
wl, wc = walk(tl, okl), walk(tc, okc)
fig, ax = plt.subplots(figsize=(4.6, 3.0))
ax.plot(cen, wl, "o-", label=f"leading edge, thr = {thr_ref}")
ax.plot(cen, wc, "s-", label=f"CFD  f = {kc[0]:.3f}, d = {kc[1]/fs:.2f} ns")
ax.set_xlabel("pulse amplitude (fraction of Amax)"); ax.set_ylabel("mean timing error [ps]")
ax.set_title("Time walk (SNR = 1000, 4 GS/s)"); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(OUT + "fig1_walk.png"); plt.close(fig)
res["walk_peak_to_peak_ps"] = {"leading_edge": float(np.ptp(wl)), "cfd": float(np.ptp(wc)),
                               "cfd_params": {"f": kc[0], "d_ns": kc[1] / fs}}

# ---- E2: resolution vs SNR at 1 GS/s ---------------------------------------------------------------
fs = 1.0
snrs = [30, 50, 100, 200, 500]
r_le, r_cf, p_cf, p_le = [], [], [], []
for s_ in snrs:
    le, cf, _ = sweep(fs, s_, 100 + s_)
    kl, vl = best(le); kc, vc = best(cf)
    r_le.append(reeval(fs, s_, 'le', kl, 5000 + s_)); r_cf.append(reeval(fs, s_, 'cfd', kc, 5000 + s_)); p_le.append(kl); p_cf.append((kc[0], kc[1]))
fig, ax = plt.subplots(figsize=(4.6, 3.0))
ax.loglog(snrs, r_le, "o-", label="leading edge (best threshold)")
ax.loglog(snrs, r_cf, "s-", label="CFD (best f, d)")
ax.set_xlabel("SNR = Amax / noise rms"); ax.set_ylabel("single-channel timing rms [ps]")
ax.set_title("Timing resolution vs SNR (1 GS/s)"); ax.legend(fontsize=7)
ax.set_xticks(snrs); ax.xaxis.set_major_formatter(ScalarFormatter()); ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_yticks([20, 50, 100, 200, 300]); ax.yaxis.set_major_formatter(ScalarFormatter()); ax.yaxis.set_minor_formatter(NullFormatter())
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout(); fig.savefig(OUT + "fig2_snr.png"); plt.close(fig)
res["vs_snr_1GSps"] = {"snr": snrs, "le_ps": [float(v) for v in r_le], "cfd_ps": [float(v) for v in r_cf],
                       "le_best_thr": [float(p) for p in p_le], "cfd_best_f_d": [list(map(float, p)) for p in p_cf]}

# ---- E3: heatmap f x d at SNR 100, 1 GS/s ----------------------------------------------------------
le, cf, _ = sweep(1.0, 100, 7)
ds = sorted({k[1] for k in cf})
H = np.array([[cf[(f, d)][0] if cf[(f, d)][1] >= 0.99 else np.nan for d in ds] for f in FR])
fig, ax = plt.subplots(figsize=(4.6, 3.0))
im = ax.imshow(H, origin="lower", aspect="auto", cmap="viridis_r", norm=LogNorm(vmin=40, vmax=5000))
ax.set_xticks(range(len(ds))); ax.set_xticklabels([f"{d}" for d in ds])
ax.set_yticks(range(len(FR))); ax.set_yticklabels([f"{f:.3f}" for f in FR])
ax.set_xlabel("delay d [samples at 1 GS/s = ns]"); ax.set_ylabel("fraction f")
ax.set_title("CFD timing rms [ps] (log colour), SNR = 100"); ax.grid(False)
fig.colorbar(im, ax=ax, shrink=0.85)
fig.tight_layout(); fig.savefig(OUT + "fig3_heatmap.png"); plt.close(fig)
kc, vc = best(cf); kl, vl = best(le)
res["heatmap_best"] = {"cfd": {"f": kc[0], "d_samples": kc[1], "rms_ps_heldout": float(reeval(1.0, 100, 'cfd', kc, 9001))},
                       "le": {"thr": kl, "rms_ps_heldout": float(reeval(1.0, 100, 'le', kl, 9001))}}

# ---- E4: resolution vs sampling rate at SNR 100 ----------------------------------------------------
fss = [0.25, 0.5, 1.0, 2.0, 4.0]
a_le, a_cf = [], []
for f_ in fss:
    le, cf, _ = sweep(f_, 100, 300 + int(10 * f_))
    kl, kc = best(le)[0], best(cf)[0]
    a_le.append(reeval(f_, 100, 'le', kl, 7000 + int(10 * f_))); a_cf.append(reeval(f_, 100, 'cfd', kc, 7000 + int(10 * f_)))
fig, ax = plt.subplots(figsize=(4.6, 3.0))
ax.loglog([f * 1000 for f in fss], a_le, "o-", label="leading edge (best threshold)")
ax.loglog([f * 1000 for f in fss], a_cf, "s-", label="CFD (best f, d)")
ax.set_xlabel("sampling rate [MS/s]"); ax.set_ylabel("single-channel timing rms [ps]")
ax.set_title("Timing resolution vs sampling rate (SNR = 100)"); ax.legend(fontsize=7)
ax.set_xticks([f * 1000 for f in fss]); ax.xaxis.set_major_formatter(ScalarFormatter()); ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_yticks([30, 50, 100, 200, 500, 1000]); ax.yaxis.set_major_formatter(ScalarFormatter()); ax.yaxis.set_minor_formatter(NullFormatter())
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout(); fig.savefig(OUT + "fig4_fs.png"); plt.close(fig)
res["vs_fs_snr100"] = {"fs_MSps": [f * 1000 for f in fss], "le_ps": [float(v) for v in a_le], "cfd_ps": [float(v) for v in a_cf]}

json.dump(res, open("../results.json", "w"), indent=1)
print(json.dumps(res, indent=1))

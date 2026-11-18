"""Generates a quantised pulse stream, runs the Verilog CFD in Icarus Verilog, and compares the
RTL output to the bit-exact Python model; then reports the timing resolution of the RTL output."""
import subprocess, os, sys, numpy as np, json
from cfd_core import *

FS, SNR, NEV = 1.0, 100, 400
HOLD = 16
D, K = 6, 3                      # best fixed-point config from the sweep (f = 3/8, d = 6 ns)
AMAX = 1500                      # LSB (12-bit ADC, full scale = 2047)
THR = int(round(6 * AMAX / SNR)) # 6 sigma arming threshold = 90 LSB
L = int(T_WIN * FS)
ts = 1 / FS

rng = np.random.default_rng(2024)
x, t0, A = make_events(NEV, FS, SNR, rng)
xi = np.clip(np.round(x * AMAX), -2048, 2047).astype(int)
stream = xi.reshape(-1)

os.makedirs("../sim", exist_ok=True)
with open("../sim/input.hex", "w") as f:
    for v in stream:
        f.write("%03x\n" % (int(v) & 0xFFF))

cmd = f"iverilog -g2005 -P tb_cfd.NS={len(stream)} -P tb_cfd.D={D} -P tb_cfd.K={K} -P tb_cfd.THR={THR} -P tb_cfd.HOLD={HOLD} -o tb_cfd ../tb/tb_cfd.v ../rtl/cfd.v && vvp tb_cfd"
r = subprocess.run(cmd, shell=True, cwd="../sim", capture_output=True, text=True)
print(r.stdout[-300:], r.stderr[-500:])

hw = [tuple(map(int, l.split())) for l in open("../sim/hw_out.txt")]
ref = cfd_fixed_stream(stream, D, K, THR, HOLD)
match = hw == ref
print(f"RTL hits: {len(hw)}   Python fixed-point model hits: {len(ref)}   bit-exact match: {match}")

# timing resolution from the RTL output (interpolation done in software)
t_hw = np.full(NEV, np.nan)
cnt = np.zeros(NEV, int)
for n, yp, yc in hw:
    e = n // L
    cnt[e] += 1
    if cnt[e] == 1:
        t_hw[e] = ((n - e * L) - 1 + yp / (yp - yc)) * ts
ok = ~np.isnan(t_hw) & (cnt == 1)
sig_hw = 1000 * (t_hw - t0)[ok].std()
# float CFD on the unquantised and on the quantised data (same f, d)
t_f, ok_f = cfd(x, K / 8, D, THR / AMAX, ts)
t_q, ok_q = cfd(xi / AMAX, K / 8, D, THR / AMAX, ts)
sig_f = stats(t_f, ok_f, t0)[0]; sig_q = stats(t_q, ok_q, t0)[0]
out = {"events": NEV, "hits_rtl": len(hw), "bit_exact_vs_python_model": bool(match),
       "events_with_exactly_one_hit": int(ok.sum()),
       "rms_ps": {"float_unquantised": float(sig_f), "float_12bit_quantised": float(sig_q), "rtl_fixed_point": float(sig_hw)},
       "params": {"fs_GSps": FS, "snr": SNR, "D": D, "K": K, "THR_LSB": THR, "HOLD": HOLD, "amax_LSB": AMAX}}
json.dump(out, open("../hw_results.json", "w"), indent=1)
print(json.dumps(out, indent=1))

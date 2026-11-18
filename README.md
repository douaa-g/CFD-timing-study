# Digital constant-fraction discriminator (CFD): timing study + Verilog implementation

Preliminary study of the *principle* of CFD timing for digitised detector pulses (time-of-flight use case).
Pulses are simulated (illustrative bi-exponential shape, 2 ns rise / 20 ns decay), so the numbers are not
those of any real detector. Author: Douaa Guermassi.

## Results (single-channel timing rms; the best parameters are picked on one event set and reported on held-out events)

| Condition | Best leading-edge threshold | Best CFD |
|---|---|---|
| Time walk, amplitudes 0.3-1, SNR 1000, 4 GS/s (peak to peak) | 469 ps | < 1 ps |
| 1 GS/s, SNR 100 | 194 ps | 47 ps (f = 3/8, d = 6 ns) |
| 1 GS/s, SNR 30 / 500 | 324 / 184 ps | 135 / 23 ps |
| SNR 100, 250 MS/s / 4 GS/s | 859 / 93 ps | 372 / 31 ps |

RTL check (400 events, 12-bit quantised, 1 GS/s, SNR 100): the Verilog output is bit-exact against the Python
fixed-point model; 399/400 events give exactly one hit; timing rms 47.8 ps (RTL) vs 47.6 ps (float).

## Layout
```
src/cfd_core.py    pulse generator, leading edge, float CFD, bit-exact fixed-point model
src/run_study.py   sweeps (SNR, sampling rate, f, d) -> figures/ and results.json
src/hw_check.py    stream -> Icarus Verilog -> compare with the Python model -> hw_results.json
rtl/cfd.v          fixed-point CFD (12-bit, y = (K*x >>> 3) - x[n-D], arming, zero crossing, hold-off)
tb/tb_cfd.v        file-based testbench (input.hex -> hw_out.txt)
rapport_cfd.pdf    2-page report (French)
```

## Reproduce
```
pip install numpy matplotlib scipy      # and: apt-get install iverilog
cd src && python3 run_study.py && python3 hw_check.py
```

## Limitations
Identical pulse shapes (a real detector adds rise-time spread, hence residual walk), white noise only, no sampling-clock
jitter, linear interpolation, illustrative sampling rate and pulse shape. The sub-sample interpolation is done in software;
the RTL outputs the sample index and y[n-1], y[n]. Bug found and fixed during verification: without the hold-off the unit
re-triggered on tail noise (738 hits for 400 events).

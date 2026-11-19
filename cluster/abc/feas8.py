import sys, time, numpy as np
sys.path.insert(0, "/home/ecastron/araucaria/scripts/abc")
import sim8 as S
for mi, mod in enumerate(("H1", "H1b", "H2", "H3")):
    for k in range(4):
        seed = 5_000_000 + mi*100 + k
        t = time.time(); r = S.one(mod, seed); dt = time.time() - t
        if r is None:
            print(f"{mod} seed {seed}: FAILED in {dt:.0f}s {getattr(S.one, 'last_error', 'None')}", flush=True); continue
        _, v, p = r; c = 10**v[-3:]
        print(f"{mod} {dt:5.0f}s  Tyr {p['T_years']:>11,.0f}  mu_yr {p['mu_yr']:.1e}  g {p['g']:4.0f}  "
              f"Ne {p['Ne']:7.0f}  Nanc {p['Nanc']:9.0f}  piC {c[0]:.5f} piA {c[1]:.5f} dxy {c[2]:.5f}", flush=True)

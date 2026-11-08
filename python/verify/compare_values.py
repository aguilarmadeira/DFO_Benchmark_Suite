"""MATLAB-Python check of the DFO Benchmark Suite (python/verify).

    python compare_values.py            compare values_matlab.csv (from eval_points.m in MATLAB; if absent, the
                                        reference file values_matlab_reference.csv) with values_python.csv
    python compare_values.py --python   first recompute values_python.csv with the Python wrappers of ../core

Fails (exit code 1) when the two files do not hold exactly the same (folder, instance, point) keys, when a value is
not finite in one file and finite in the other, or when a relative difference |p - m| / max(1, |m|) exceeds 1e-12."""
import os, sys, math, importlib.util
TOL = 1e-12
here = os.path.dirname(os.path.abspath(__file__)); core = os.path.join(here, '..', 'core')

def load(f):
    out = {}
    for n, l in enumerate(open(os.path.join(here, f)), 1):
        c = l.rstrip('\n').split(',')
        key = tuple(c[:3])
        if len(c) != 4 or key in out: sys.exit(f'{f}, line {n}: malformed or repeated row')
        out[key] = float(c[3])
    return out

if '--python' in sys.argv:
    out = open(os.path.join(here, 'values_python.csv'), 'w'); cur = (None, None); PB = None
    for l in open(os.path.join(here, 'points.csv')):
        fo, name, j, z = l.rstrip('\n').split(',')
        if (fo, name) != cur:
            spec = importlib.util.spec_from_file_location(name, os.path.join(core, fo, name + '.py'))
            mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); PB = getattr(mod, name)(); cur = (fo, name)
        out.write(f"{fo},{name},{j},{PB.func_idx([float(v) for v in z.split(';')])!r}\n")
    out.close()

mfile = 'values_matlab.csv' if os.path.exists(os.path.join(here, 'values_matlab.csv')) else 'values_matlab_reference.csv'
m = load(mfile); p = load("values_python.csv")
keys_points = {tuple(l.split(',')[:3]) for l in open(os.path.join(here, 'points.csv'))}
fail = []
if set(m) != keys_points or set(p) != keys_points:
    fail.append(f'key sets differ: points {len(keys_points)}, {mfile} {len(m)}, values_python.csv {len(p)} '
                f'(missing in MATLAB {len(keys_points - set(m))}, in Python {len(keys_points - set(p))})')
rows = []
for k in sorted(keys_points & set(m) & set(p)):
    a, b = p[k], m[k]
    if math.isfinite(a) != math.isfinite(b) or (math.isnan(a) != math.isnan(b)):
        rows.append((math.inf, k, a, b)); continue
    if not math.isfinite(a):
        d = 0.0 if a == b else math.inf
    else:
        d = 0.0 if a == b else abs(a - b) / max(1.0, abs(b))
    rows.append((d, k, a, b))
same = sum(1 for d, *_ in rows if d == 0); worst = sorted(rows, key=lambda r: -r[0])[:10]; bad = [r for r in rows if r[0] > TOL]
print(f'{mfile} vs values_python.csv: {len(rows)} points, identical {same} ({100 * same / max(1, len(rows)):.1f}%), '
      f'max relative difference {worst[0][0]:.3g}, above {TOL:g}: {len(bad)}')
print('largest differences:')
for d, k, a, b in worst: print(f'   {d:.3g}  {k[0]}/{k[1]} point {k[2]}: python {a!r}, matlab {b!r}')
if bad: fail.append(f'{len(bad)} points above the tolerance')
if fail:
    print('FAIL: ' + '; '.join(fail)); sys.exit(1)
print('PASS')

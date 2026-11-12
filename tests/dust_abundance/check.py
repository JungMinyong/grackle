#!/usr/bin/env python3
"""Test the physical fit and real solver against explicit dust-density cells."""
import argparse
import ctypes
import math
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--library", type=Path, required=True)
parser.add_argument("--cell", type=Path, required=True)
parser.add_argument("--data", type=Path, required=True)
args = parser.parse_args()
library = ctypes.CDLL(str(args.library.resolve()))
functions = {}
for name in ["remy_ruyer_dust_to_gas_ratio"]:
    fn = getattr(library, name + "_")
    fn.argtypes = [ctypes.POINTER(ctypes.c_double)] * 2
    functions[name] = fn


def call(name, metallicity):
    z, result = ctypes.c_double(metallicity), ctypes.c_double()
    functions[name](ctypes.byref(z), ctypes.byref(result))
    return result.value


def fit(z):
    if z <= 0:
        return 0.
    x = math.log10(z / .014)
    return 10. ** (-(.68 - 3.08*x if x <= -.73 else 2.21-x))


def cell(model, h2_switch, dust_field, z, ratio):
    result = subprocess.run([
        str(args.cell.resolve()), str(args.data.resolve()), str(model),
        str(h2_switch), str(dust_field), repr(z), repr(ratio)],
        text=True, capture_output=True)
    if result.returncode:
        raise AssertionError(result.stdout + result.stderr)
    line = next(line for line in result.stdout.splitlines()
                if line.startswith("RESULT "))
    fields = line.split()[1:]
    return tuple(map(int, fields[:2])), tuple(map(float, fields[2:]))


def equal(a, b):
    assert all(math.isclose(x, y, rel_tol=1.e-9, abs_tol=1.e-20)
               for x, y in zip(a, b)), (a, b)


transition = .014 * 10.**(-.73)
metallicities = [0., 1.e-8, 5.4e-4, transition*(1-1.e-6),
                 transition, transition*(1+1.e-6), .014, .028]
for z in [-1.e-10] + metallicities:
    assert math.isclose(call("remy_ruyer_dust_to_gas_ratio", z),
                        fit(z), rel_tol=1.e-12, abs_tol=1.e-30)

for z in metallicities:
    physical = fit(z)
    corrected_flags, corrected = cell(5, 1, 0, z, 0.)
    assert corrected_flags == (5, 1)
    equal(corrected, cell(1, 1, 1, z, physical)[1])
    equal(cell(5, 1, 1, z, physical*2)[1],
          cell(1, 1, 1, z, physical*2)[1])
    equal(cell(1, 1, 0, z, 0.)[1],
          cell(1, 1, 1, z, .009387*z/.01295)[1])
    # Dust chemistry enables the internal switch without a public H2 request.
    equal(corrected, cell(5, 0, 0, z, 0.)[1])
    print(f"PASS Z={z:.9g}: modes 5/1, field override, automatic H2 activation")

# Reject invalid selectors rather than using an uninitialized dust2gas value.
for model in [2, 3, 4, 6, -1]:
    bad = subprocess.run([str(args.cell.resolve()), str(args.data.resolve()),
                          str(model), "1", "0", ".00054", "0"],
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    assert bad.returncode == 3
print("PASS: unsupported dust modes rejected, including mode 4")

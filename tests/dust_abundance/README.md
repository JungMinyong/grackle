# Corrected Rémy-Ruyer dust abundance

On `aeos-v2`, `dust_chemistry=5` enables the standard dust processes with an
absolute broken-power-law dust/gas mass ratio. Mode 1 retains its linear
abundance, mode 0 disables dust chemistry, and all other modes (including 4)
are rejected. No chemistry structure member or LW prescription is added.
`h2_on_dust` retains its existing boolean meaning.

For metal mass fraction `Z > 0`, let `x = log10(Z / 0.014)`:

```
D/G = 10**(-0.68 + 3.08*x)  if x <= -0.73
D/G = 10**(-2.21 + x)      otherwise
```

For nonpositive Z the ratio is zero. An explicit dust-density field overrides
the inferred abundance. `local_dust_to_gas_ratio` still normalizes the existing
dust rate coefficients; it does not rescale the empirical mode-5 abundance.
The existing fit coefficients are retained, including their small discontinuity
at the break.

Build a double-precision shared Grackle library, then compile `cell.c` against
its matching public and generated headers:

```bash
cc -std=c99 -Isrc/include -I/path/to/build/generated_public_headers \
  tests/dust_abundance/cell.c /path/to/libgrackle.so \
  -Wl,-rpath,/path/to/library -lm -o /tmp/grackle-dust-cell
python3 tests/dust_abundance/check.py \
  --library /path/to/libgrackle.so \
  --cell /tmp/grackle-dust-cell \
  --data /path/to/CloudyData_HM2012_highz_shielded.h5
```

The regression checks the analytic helper and compares real cooling times,
dust temperatures, H2 abundances, and energy after one year against cells
with independently specified dust densities. It covers eight metallicities,
including zero, the fiducial Z=0.00054, both sides of the break, the break,
and solar and supersolar values. It also checks mode-1 behavior, explicit
field precedence, automatic H2 activation, and rejection of unsupported modes.
These controlled cells do not predict a galaxy's evolved H2 or stellar mass.

# Dust-normalization regression

The corrected model passes a physical dust/gas mass ratio to the existing dust
rates. The historical normalized helper is unchanged and remains available
through `dust_chemistry=4`; corrected abundance uses `dust_chemistry=5`.
The standard path and explicit dust field
retain their conventions. No LW or Enzo photoelectric-heating code is changed.

`check.py` loads the actual library's Fortran helpers and compares the real
single-cell cooling, dust temperature, H2 abundance, and energy after one year
with results obtained by supplying an independent dust-density field.
It covers zero and low metallicity, the fiducial metallicity, both sides of the
break, the exact break, solar and supersolar metallicity; it also checks
automatic H2 activation, field precedence, and invalid-selector rejection.
The old helper's default-REAL rounding is retained in the legacy reference.

Build Grackle with matching 64-bit headers and library, then:

```bash
cc -std=c99 -DCONFIG_BFLOAT_8 -Isrc/clib \
  tests/dust_abundance/cell.c \
  -L/path/to/library -Wl,-rpath,/path/to/library \
  -lgrackle -lm -o /tmp/grackle-dust-cell
python3 tests/dust_abundance/check.py \
  --library /path/to/library/libgrackle.so \
  --cell /tmp/grackle-dust-cell \
  --data /path/to/CloudyData_HM2012_highz_shielded.h5
```

Tests use non-cosmological, initially neutral gas near 100 K with density
1000 in units of `1.6735575e-24 g/cm^3`, molecular chemistry, metal cooling,
and dust chemistry enabled, with UV background and photoelectric heating off.
These controlled cells test the abundance/rate interface and finite chemistry
evolution. They do not measure the change in a galaxy's H2 mass or star
formation.

The C/Fortran parameter layout and defaults are checked separately:

```bash
cc -DCONFIG_BFLOAT_8 -Isrc/clib -c tests/dust_abundance/abi.c -o /tmp/dust-abi.o
gfortran -DCONFIG_BFLOAT_8 -Isrc/clib tests/dust_abundance/abi.F /tmp/dust-abi.o \
  -L/path/to/library -Wl,-rpath,/path/to/library -lgrackle -o /tmp/dust-abi
/tmp/dust-abi
```

For aeos-v2, the updated adapter recognizes `GRACKLE_HAS_REMY_RUYER_DUST`:
`dust_chemistry=5` selects the corrected absolute fit and `4` selects historical
behavior after rebuilding against this library. The separate
`GrackleDustAbundanceModel` and Grackle `dust_abundance_model` parameters have
been removed. No field is added to the C `chemistry_data` structure.
Unsupported library builds must reject modes 4/5.

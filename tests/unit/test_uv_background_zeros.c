/* Exercise the production UV reader with synthetic two-point tables.
 * No external data files or chemistry initialization are needed. */
#include <math.h>
#include <stdio.h>
#include "grackle.h"

int update_UVbackground_rates(chemistry_data *, chemistry_data_storage *,
                             photo_rate_storage *, code_units *);

static int check(double lo, double hi, double redshift)
{
  chemistry_data c = {0};
  chemistry_data_storage r = {0};
  photo_rate_storage uv = {0};
  code_units u = {0};
  double z[] = {0.0, 50.0}, rate[] = {lo, hi};
  c.UVbackground = 1;
  c.primordial_chemistry = 2;
  c.self_shielding_method = 3;
  c.UVbackground_redshift_on = c.UVbackground_redshift_fullon = 50;
  c.UVbackground_redshift_off = c.UVbackground_redshift_drop = 0;
  u.a_units = 1;
  u.a_value = 1 / (1 + redshift);
  u.length_units = u.time_units = u.velocity_units = 1;
  u.density_units = 1.67262171e-24;
  r.UVbackground_table.Nz = 2;
  r.UVbackground_table.z = z;
#define ASSIGN(name) r.UVbackground_table.name = rate
  ASSIGN(k24); ASSIGN(k25); ASSIGN(k26); ASSIGN(k27);
  ASSIGN(k28); ASSIGN(k29); ASSIGN(k30); ASSIGN(k31);
  ASSIGN(piHI); ASSIGN(piHeI); ASSIGN(piHeII);
  ASSIGN(crsHI); ASSIGN(crsHeI); ASSIGN(crsHeII);
#undef ASSIGN
  if (!update_UVbackground_rates(&c, &r, &uv, &u)) return 1;
  double actual_z = 1 / u.a_value - 1;
  double w = fmax(0, fmin(1, log1p(actual_z) / log(51.0)));
  double expected = lo > 0 && hi > 0 ?
    exp((1-w)*log(lo) + w*log(hi)) : (1-w)*lo + w*hi;
  if (actual_z > 50 || actual_z < 0) expected = 0;
  double results[] = {
    uv.k24, uv.k25, uv.k26, uv.k27, uv.k28, uv.k29, uv.k30, uv.k31,
    uv.piHI * 1.67262171e-24 / 1.60217653e-12,
    uv.piHeI * 1.67262171e-24 / 1.60217653e-12,
    uv.piHeII * 1.67262171e-24 / 1.60217653e-12,
    uv.crsHI, uv.crsHeI, uv.crsHeII
  };
  /* At z=on the ramp is zero for rates/heating, but not cross sections. */
  for (int i=0; i<14; ++i) {
    double ref = actual_z == 50 && i < 11 ? 0 : expected;
    if (!isfinite(results[i]) || results[i] < 0 ||
        fabs(results[i]-ref) > 1e-12*fabs(ref)) {
      fprintf(stderr, "lo=%g hi=%g z=%.17g channel=%d: %.17g != %.17g\n",
              lo, hi, redshift, i, results[i], ref);
      return 1;
    }
  }
  return 0;
}

int main(void)
{
  const double endpoints[][2] = {
    {0,0}, {0,1e-12}, {1e-12,0}, {1e-12,1e-14}, {1e-20,1e-20}
  };
  const double z[] = {0, 0.01, 10, 25, 49.98987025196391, 50, 50.001, 70};
  for (int i=0; i<5; ++i)
    for (int j=0; j<8; ++j)
      if (check(endpoints[i][0], endpoints[i][1], z[j])) return 1;
  puts("PASS: 40 UV-table cases, 14 channels: zeros, transitions, positive rates, endpoints and off-range controls");
  return 0;
}

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "grackle.h"

/* One fully specified cell. Compare inferred abundance with an independently
 * supplied dust-density field through the real cooling and chemistry APIs. */
int main(int argc, char **argv)
{
  if (argc != 7) return 2;
  chemistry_data chemistry = _set_default_chemistry_parameters();
  chemistry_data_storage rates = {0};
  code_units units = {0};
  grackle_field_data fields = {0};
  chemistry.use_grackle = 1;
  chemistry.primordial_chemistry = 2;
  chemistry.metal_cooling = 1;
  chemistry.dust_chemistry = atoi(argv[2]);
  chemistry.photoelectric_heating = 0;
  chemistry.cmb_temperature_floor = 0;
  chemistry.UVbackground = 0;
  chemistry.grackle_data_file = argv[1];
  chemistry.h2_on_dust = atoi(argv[3]);
  chemistry.use_dust_density_field = atoi(argv[4]);
  const double metallicity = atof(argv[5]), dust_ratio = atof(argv[6]);
  units.density_units = 1.6735575e-24;
  units.length_units = 1.0;
  units.time_units = 1.0;
  units.velocity_units = 1.0;
  units.a_units = units.a_value = 1.0;
  if (!_initialize_chemistry_data(&chemistry, &rates, &units)) return 3;
  int dim[3] = {1,1,1}, begin[3] = {0,0,0}, end[3] = {0,0,0};
  fields.grid_rank = 1;
  fields.grid_dimension = dim; fields.grid_start = begin; fields.grid_end = end;
  fields.grid_dx = 1.0;
  gr_float rho = 1000., hi = 760., hii = 1.e-17, hm = 1.e-17;
  gr_float hei = 240., heii = 1.e-17, heiii = 1.e-17;
  gr_float h2 = 1.e-17, h2ii = 1.e-17, electron = 1.e-17;
  gr_float metal = rho * metallicity, dust = rho * dust_ratio;
  gr_float energy = 100. * 1.380649e-16 /
    (1.6735575e-24 * (1. / (0.76 + 0.24/4.)) * (chemistry.Gamma-1.));
  gr_float vx = 0., vy = 0., vz = 0.;
  fields.density = &rho; fields.HI_density = &hi; fields.HII_density = &hii;
  fields.HM_density = &hm; fields.HeI_density = &hei;
  fields.HeII_density = &heii; fields.HeIII_density = &heiii;
  fields.H2I_density = &h2; fields.H2II_density = &h2ii;
  fields.e_density = &electron; fields.metal_density = &metal;
  fields.dust_density = &dust; fields.internal_energy = &energy;
  fields.x_velocity = &vx; fields.y_velocity = &vy; fields.z_velocity = &vz;
  gr_float cooling_time = 0., dust_temperature = 0.;
  if (!local_calculate_cooling_time(&chemistry, &rates, &units,
                                    &fields, &cooling_time)) return 4;
  if (!local_calculate_dust_temperature(&chemistry, &rates, &units,
                                        &fields, &dust_temperature)) return 5;
  if (!local_solve_chemistry(&chemistry, &rates, &units, &fields, 31557600.)) return 6;
  if (!isfinite(cooling_time) || !isfinite(dust_temperature) ||
      !isfinite(h2) || !isfinite(energy) || h2 < 0. || energy <= 0.) return 7;
  printf("RESULT %d %d %.17g %.17g %.17g %.17g\n",
         chemistry.dust_chemistry, chemistry.h2_on_dust,
         (double)cooling_time, (double)dust_temperature, (double)h2, (double)energy);
  /* The legacy free routine has no defined return value. */
  _free_chemistry_data(&chemistry, &rates);
  return 0;
}

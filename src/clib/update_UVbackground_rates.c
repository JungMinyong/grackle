/***********************************************************************
/
/ Update UV background rates to current redshift
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"

/* function prototypes */

/* UV tables may contain physical zeros (e.g. before reionization).
 * Preserve logarithmic rate interpolation for positive endpoints. Otherwise
 * interpolate rates in log(1+z), preserving zero intervals exactly. Clamp
 * the fallback weight to avoid negative rates from endpoint roundoff. */
static double interpolate_uv_rate(double y0, double y1,
                                  double dz, double interval)
{
  if (y0 == 0.0 && y1 == 0.0) return 0.0;
  if (y0 <= 0.0 || y1 <= 0.0) {
    double w = fmax(0.0, fmin(1.0, dz / interval));
    return (1.0 - w) * y0 + w * y1;
  }
  return exp(dz * (log(y1 / y0) / interval) + log(y0));
}

int update_UVbackground_rates(chemistry_data *my_chemistry,
                              chemistry_data_storage *my_rates,
                              photo_rate_storage *my_uvb_rates,
                              code_units *my_units)
{
  /* Return if there is no radiation (rates should be all zero). */

  if (my_chemistry->UVbackground == 0 ||
      my_chemistry->primordial_chemistry == 0)
    return SUCCESS;

  /* Return if redshift is outside of on/off redshifts. */

  double Redshift = 1.0 / (my_units->a_value * my_units->a_units) - 1;
  if ( (Redshift < my_chemistry->UVbackground_redshift_off) ||
       (Redshift > my_chemistry->UVbackground_redshift_on) )
    return SUCCESS;

  /* ------------------------------------------------------------------ */
  /* First, calculate the ramp value, a number between 0 and 1 which
     is used as an external control to the radiation. */

  double Ramp = 0;

  if (Redshift < my_chemistry->UVbackground_redshift_on && 
      Redshift >= my_chemistry->UVbackground_redshift_off) {

    if (Redshift > my_chemistry->UVbackground_redshift_fullon)
      Ramp = 0.5 - 0.5 * tanh(15.0 * (Redshift - 0.5 *
        (my_chemistry->UVbackground_redshift_on +
         my_chemistry->UVbackground_redshift_fullon)));
    else if (Redshift < my_chemistry->UVbackground_redshift_drop)
      Ramp = 0.5 - 0.5 * tanh(15.0 * (0.5 *
        (my_chemistry->UVbackground_redshift_drop +
         my_chemistry->UVbackground_redshift_off) - Redshift));
    else
      Ramp = 1.0;

  }


  /* Interpolate the UV background table. */
  
  // find interpolation index
  double *zvec = my_rates->UVbackground_table.z;
  int index=0;
  while (Redshift > zvec[index])
    index++;
  if(index == 0) index=1;
  if(index == my_rates->UVbackground_table.Nz) index--;
  double zvec_grad = log((1+zvec[index]) / (1+zvec[index-1]));
  double redshift_grad = log((1+Redshift) / (1+zvec[index-1]));
  // printf("index = %d, %.3f <= %.3f <= %.3f\n",index,zvec[index-1],Redshift,zvec[index]);

  // *** k24 ***
  my_uvb_rates->k24 = interpolate_uv_rate(
    my_rates->UVbackground_table.k24[index-1],
    my_rates->UVbackground_table.k24[index], redshift_grad, zvec_grad);

  // *** k25 ***
  my_uvb_rates->k25 = interpolate_uv_rate(
    my_rates->UVbackground_table.k25[index-1],
    my_rates->UVbackground_table.k25[index], redshift_grad, zvec_grad);

  // *** k26 ***
  my_uvb_rates->k26 = interpolate_uv_rate(
    my_rates->UVbackground_table.k26[index-1],
    my_rates->UVbackground_table.k26[index], redshift_grad, zvec_grad);

  if (my_chemistry->primordial_chemistry > 1) {

    // *** k27 ***
    my_uvb_rates->k27 = interpolate_uv_rate(
    my_rates->UVbackground_table.k27[index-1],
    my_rates->UVbackground_table.k27[index], redshift_grad, zvec_grad);

    // *** k28 ***
    my_uvb_rates->k28 = interpolate_uv_rate(
    my_rates->UVbackground_table.k28[index-1],
    my_rates->UVbackground_table.k28[index], redshift_grad, zvec_grad);

    // *** k29 ***
    my_uvb_rates->k29 = interpolate_uv_rate(
    my_rates->UVbackground_table.k29[index-1],
    my_rates->UVbackground_table.k29[index], redshift_grad, zvec_grad);

    // *** k30 ***
    my_uvb_rates->k30 = interpolate_uv_rate(
    my_rates->UVbackground_table.k30[index-1],
    my_rates->UVbackground_table.k30[index], redshift_grad, zvec_grad);

    // *** k31 ***
    my_uvb_rates->k31 = interpolate_uv_rate(
    my_rates->UVbackground_table.k31[index-1],
    my_rates->UVbackground_table.k31[index], redshift_grad, zvec_grad);

  }

  // *** piHI ***
  my_uvb_rates->piHI = interpolate_uv_rate(
    my_rates->UVbackground_table.piHI[index-1],
    my_rates->UVbackground_table.piHI[index], redshift_grad, zvec_grad);

  // *** piHeII ***
  my_uvb_rates->piHeII = interpolate_uv_rate(
    my_rates->UVbackground_table.piHeII[index-1],
    my_rates->UVbackground_table.piHeII[index], redshift_grad, zvec_grad);

  // *** piHeI ***
  my_uvb_rates->piHeI = interpolate_uv_rate(
    my_rates->UVbackground_table.piHeI[index-1],
    my_rates->UVbackground_table.piHeI[index], redshift_grad, zvec_grad);

  //
  // Cross sections are read from table in cgs, and remain in cgs.
  // Analytic self-shielding method assumes cross sections in cgs units
  //

  // *** crsHI ***
  if (my_chemistry->self_shielding_method > 0){
    my_uvb_rates->crsHI = interpolate_uv_rate(
    my_rates->UVbackground_table.crsHI[index-1],
    my_rates->UVbackground_table.crsHI[index], redshift_grad, zvec_grad);

    // *** crsHeI ***
    my_uvb_rates->crsHeI = interpolate_uv_rate(
    my_rates->UVbackground_table.crsHeI[index-1],
    my_rates->UVbackground_table.crsHeI[index], redshift_grad, zvec_grad);

    // *** crsHeII ***
    my_uvb_rates->crsHeII = interpolate_uv_rate(
    my_rates->UVbackground_table.crsHeII[index-1],
    my_rates->UVbackground_table.crsHeII[index], redshift_grad, zvec_grad);
  }

  // Now convert the rates to code units.

  /* Get conversion units. */

  double co_length_units, co_density_units;
  if (my_units->comoving_coordinates == TRUE) {
    co_length_units = my_units->length_units;
    co_density_units = my_units->density_units;
  }
  else {
    co_length_units = my_units->length_units *
      my_units->a_value * my_units->a_units;
    co_density_units = my_units->density_units /
      POW(my_units->a_value * my_units->a_units, 3);
  }

  double tbase1 = my_units->time_units;
  double xbase1 = co_length_units /
    (my_units->a_value * my_units->a_units);
  double dbase1 = co_density_units *
    POW(my_units->a_value * my_units->a_units, 3);
  double mh     = 1.67262171e-24;
  double ev2erg = 1.60217653e-12;
  /* compared to Enzo source, there's an additional factor of
     1/ev2erg here, because the heating rates are stored as eV/s. */
  double CoolingUnits = (POW(my_units->a_units, 5) * xbase1*xbase1 * mh*mh) /
    (POW(tbase1, 3) * dbase1) / ev2erg;


  my_uvb_rates->k24 *= my_units->time_units;
  my_uvb_rates->k25 *= my_units->time_units;
  my_uvb_rates->k26 *= my_units->time_units;
  
  if (my_chemistry->primordial_chemistry > 1) {
    my_uvb_rates->k27 *= my_units->time_units;
    my_uvb_rates->k28 *= my_units->time_units;
    my_uvb_rates->k29 *= my_units->time_units;
    my_uvb_rates->k30 *= my_units->time_units;
    my_uvb_rates->k31 *= my_units->time_units;
  }
    
  my_uvb_rates->piHI /= CoolingUnits;
  my_uvb_rates->piHeII /= CoolingUnits;
  my_uvb_rates->piHeI /= CoolingUnits;

  // Now apply the Ramp factor
  my_uvb_rates->k24 *= Ramp;
  my_uvb_rates->k25 *= Ramp;
  my_uvb_rates->k26 *= Ramp;
  if (my_chemistry->primordial_chemistry > 1) {
    my_uvb_rates->k27 *= Ramp;
    my_uvb_rates->k28 *= Ramp;
    my_uvb_rates->k29 *= Ramp;
    my_uvb_rates->k30 *= Ramp;
    my_uvb_rates->k31 *= Ramp;
  }
  my_uvb_rates->piHI *= Ramp;
  my_uvb_rates->piHeII *= Ramp;
  my_uvb_rates->piHeI *= Ramp;

  /* Molecular hydrogen constant photo-dissociation */

  /* Note that k31 can be set by above by the UV background table, in
     which case it is overwritten here if (LWbackground_intensity >
     0.0). */

  if (my_chemistry->LWbackground_intensity > 0.0) 
    my_uvb_rates->k31 = 1.38e-12 * my_chemistry->LWbackground_intensity *
      my_units->time_units;
  
  /* LWbackground_sawtooth_suppression is supposed to account for the
     suppression of LW flux due to Lyman-series absorption (giving a
     sawtooth pattern), a la Haiman & Abel, & Rees (2000).  This is
     really only applicable when there is plenty of neutral hydrogen
     is around, so I'm scaling it with Ramp: when Ramp==1 (fully
     re-ionized universe) there's no suppression. */

  if (my_chemistry->LWbackground_sawtooth_suppression) {
    double LymanSawtoothSuppressionFactor = 0.1 + 0.9 * Ramp;
  
    my_uvb_rates->k31 *= LymanSawtoothSuppressionFactor;
  }

  /* Compton X-ray heating */

  if (my_chemistry->Compton_xray_heating) {

    double RedshiftXrayCutoff = 5.0;

    /* This is sigma_thompson * c * (effective <h \nu>/<m_e c^2>) *
       U_xray * 1eV.  U_xray is the energy density of XRB in , <h \nu>
       is the average photon energy in keV, corrected for relativistic
       effects.  Eq.(4) and Eq.(11) of Madau & Efstathiou (1999) */
    
    my_uvb_rates->comp_xray = 4.15e-13 * 3.0e10 *
      (31.8*POW(1.0+Redshift, 0.3333)/511.0) * 
      (6.3e-5 * 1.6e-12) * 
      POW(1.0 + Redshift, 4) * 
      exp(-POW(Redshift/RedshiftXrayCutoff, 2)) / 
      CoolingUnits; 
  
    /* The effective temperature (in K).  Eq.(10) of Madau &
       Efstathiou (1999) with U_xray(z=0) = 6.3e-5 and U_cmb(z=0) =
       0.256 eV/cm3 */
    
    my_uvb_rates->temp_xray = 31.8e3*POW(1.0+Redshift, 0.3333)*1.6e-12/
      (4.0*1.38e-16) *
      6.3e-5 * POW(1.0 + Redshift, 4) * 
      exp(-POW(Redshift/RedshiftXrayCutoff, 2)) /
      (0.256 * (1+Redshift));  

  }

  return SUCCESS;
}

/***********************************************************************
/
/  READ AND INTERPOLATE DISK GALAXY AND HALO DATA
/
/  written by: Wolfram Schmidt
/  date:       June 2023
/  modified1:  
/
/  PURPOSE: reads and interpolates data for disk galaxies
/           size of table limited by MAX_DISK_ZONES
/
************************************************************************/

#include <cmath> 
#include <fstream>

#include "preincludes.h" 
#include "EquilibriumGalaxyDisk.h"
#include "global_data.h"
#include "phys_constants.h"
 
void EquilibriumGalaxyDisk::ReadInData(char *fname)
{
  //
  // Read in gas density and circular velocity data from file. Units
  // are assumed to be pc (radius and height), kg/m^3 (density)
  // and m/s (velocity)
  //
  std::fstream inFile;

  inFile.open(fname, std::ios::in);

  if (!(inFile.is_open()))
    ENZO_VFAIL("ReadInData: failed to open disk data file  %s\n", fname);

  double r, z, rho, vcirc;
  double r_prev = -1.0;

  int i = 0; 
  int j = 0;
  int n = 0; // used for debugging output below

  if (debug && (MyProcessorNumber == ROOT_PROCESSOR))
    printf("ReadInData: %s\n", fname);

  while(inFile >> r >> z >> rho >> vcirc)
  {
    if (r > r_prev) {
      // record radial zones if midplane (z = 0)
      if (j == 0)
	      this->gas_disk_zones_r[i] = r * pc_cm;
    } else {
      this->gas_disk_nr = i;
      // next vertical zone
      i = 0;
      this->gas_disk_zones_z[j++] = z * pc_cm;
    }

    this->gas_disk_log_rho[i][j] = log(1e-3*rho); // SI to cgs
    this->gas_disk_vcirc[i][j] = 1e2*vcirc; // SI to cgs
    
    /* uncomment to check if table is read in correctly
    if (debug && (MyProcessorNumber == ROOT_PROCESSOR)) {
      n++;
      printf("%d, i=%d, r=%f, j=%d, z=%f, %f, %e\n",
	           n, i, this->gas_disk_zones_r[i], j, z, this->gas_disk_log_rho[i][j], this->gas_disk_vcirc[i][j]);
    } */

    i++;
    r_prev = r;

    if (i >= MAX_DISK_ZONES || j >= MAX_DISK_ZONES)
      ENZO_VFAIL("ReadInData: table size exceeds maximum of %"ISYM" rows/columns\n", MAX_DISK_ZONES);
  }
  this->gas_disk_nz = j+1;

  inFile.close();

  if (debug && (MyProcessorNumber == ROOT_PROCESSOR))
    printf("ReadInData: number of zones n_r=%d, n_z=%d\n", this->gas_disk_nr, this->gas_disk_nz);

  return;
}

double EquilibriumGalaxyDisk::InterpolateEquilibriumDensityTable(double r_cyl, double z)
{
  //
  // Interpolate the density from the read-in table. 
  //
  double rho;

  double w_r =   r_cyl/(this->gas_disk_zones_r[1] - this->gas_disk_zones_r[0]);
  double w_z = fabs(z)/(this->gas_disk_zones_z[1] - this->gas_disk_zones_z[0]);

  // integer parts of coordinates divided by zone widths
  int i = w_r;
  int j = w_z;

  // interpolation weights are given by remainders
  w_r = w_r - i;
  w_z = w_z - j;

  if (i+1 < this->gas_disk_nr && j+1 < this->gas_disk_nz) { 
    rho = exp((1.0 - w_r)*(1.0 - w_z) * this->gas_disk_log_rho[i][j] +
	            (1.0 - w_r)* w_z        * this->gas_disk_log_rho[i][j+1] +
	             w_r       *(1.0 - w_z) * this->gas_disk_log_rho[i+1][j] +
	             w_r       * w_z        * this->gas_disk_log_rho[i+1][j+1]);
  } else {
    rho = 0.0;
  }

  return rho;
}

double EquilibriumGalaxyDisk::InterpolateEquilibriumVcircTable(double r_cyl, double z)
{
  //
  // Interpolate the circular velocity from the read-in table. 
  //
  double vcirc;

  double w_r =   r_cyl/(this->gas_disk_zones_r[1] - this->gas_disk_zones_r[0]);
  double w_z = fabs(z)/(this->gas_disk_zones_z[1] - this->gas_disk_zones_z[0]);

  // integer parts of coordinates divided by zone widths
  int i = w_r;
  int j = w_z;

  // interpolation weights are given by remainders
  w_r = w_r - i;
  w_z = w_z - j;

  if (i+1 < this->gas_disk_nr && j+1 < this->gas_disk_nz) { 
    vcirc = (1.0 - w_r)*(1.0 - w_z) * this->gas_disk_vcirc[i][j] +
	          (1.0 - w_r)* w_z        * this->gas_disk_vcirc[i][j+1] +
	           w_r       *(1.0 - w_z) * this->gas_disk_vcirc[i+1][j] +
	           w_r       * w_z        * this->gas_disk_vcirc[i+1][j+1];

  } else {
    vcirc = 0.0;
  }

  return vcirc;
}

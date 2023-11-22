/**************************************************************************
 *
 *  EQUILIBRIUM GALAXY DISK CLASS
 *
 *  written by: Wolfram Schmidt
 *  date:       June, 2023
 *  modified1:
 *
 *  PURPOSE: reads and interpolates data for disk galaxies
 *           size of table limited by MAX_DISK_ZONES
 *
 **************************************************************************/

#ifndef EQUILIBRIUM_GALAXY_DISK_DEFINED__
#define EQUILIBRIUM_GALAXY_DISK_DEFINED__

#define MAX_DISK_ZONES 500

#include "macros_and_parameters.h"
#include "typedefs.h"

class EquilibriumGalaxyDisk {

 private:

  int gas_disk_nr, gas_disk_nz;

  double gas_disk_zones_r[MAX_DISK_ZONES];
  double gas_disk_zones_z[MAX_DISK_ZONES];

  double gas_disk_log_rho[MAX_DISK_ZONES][MAX_DISK_ZONES];
  double gas_disk_vcirc[MAX_DISK_ZONES][MAX_DISK_ZONES];

 public:

  //
  // Constructor
  //
  EquilibriumGalaxyDisk() {
    gas_disk_nr = 0;
    gas_disk_nz = 0;
  };

  // 
  // read in data file
  //
  void ReadInData(char *fname);

  // 
  // interpolate data
  //
  double InterpolateEquilibriumDensityTable(double r_cyl, double z);
  double InterpolateEquilibriumVcircTable(double r_cyl, double z);
};

#endif

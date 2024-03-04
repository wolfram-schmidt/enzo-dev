/***********************************************************************
/
/  INITIALIZE DISK GALAXIES WITH LIVE HALOS
/
/  written by: Simon Selg and Wolfram Schmidt 
/              (based on code by Kai Rodenbeck)
/  date:       July, 2020
/  modified:   October, 2023
/
/  PURPOSE:
/    Reads initial conditions from files and
/    sets up one or more disk galaxies and particle halos.
/    Use https://github.com/wolfram-schmidt/isolated-galaxy
/    or equivalent tool to generate initial conditions.
/
************************************************************************/

#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>

#include "preincludes.h"
#include <stdlib.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "phys_constants.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "EquilibriumGalaxyDisk.h"

/********************* PROTOTYPES *********************/

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, float Time);

/*******************************************************/

static float MetallicityFloor = 1e-12;
static float CosmologySimulationInitialFractionHII   = 1.2e-5;
static float CosmologySimulationInitialFractionHeII  = 1.0e-14;
static float CosmologySimulationInitialFractionHeIII = 1.0e-17;
static float CosmologySimulationInitialFractionHM    = 2.0e-9;
static float CosmologySimulationInitialFractionH2I   = 2.0e-20;
static float CosmologySimulationInitialFractionH2II  = 3.0e-14;

// auxiliary function for computing the temperature profile of CGM

/*
double F0(double r, double r_c, double a) {
  double asqrp = a*a + r_c*r_c;
  double asqrm = a*a - r_c*r_c;
  return (pi * asqrm / 2 + r_c * asqrp / (a + r) - asqrm * atan(r/r_c)
          - r_c * a * log((a + r)*(a + r) / (r*r + r_c*r_c))) 
         / (r_c * asqrp);
}

double F1(double r, double r_c) {
  double atanr = atan(r/r_c);
  return pi*pi / (8*r_c) - atanr*atanr / (2*r_c) - atanr / r;
}
*/

double compute_tempr_cgm(double x, double a, double b, double c) {
  return (b * (1 + a*a - (a*a - 1) * (a + x) * atan(x) 
               + a*(a + x) * (log(x*x + 1) - 2 * log(a + x))) 
              / (pow((a*a + 1), 2) * (a + x)) 
          + c) * (x*x + 1);
}


int grid::GalaxyLiveHaloInitializeGrid(int NumberOfSpheres,
                                       char* HaloDataFile[MAX_SPHERES],
                                       EquilibriumGalaxyDisk DiskTable[MAX_SPHERES],
                                       float SpherePosition[MAX_SPHERES][MAX_DIMENSION],
                                       float SphereRotAxis[MAX_SPHERES][MAX_DIMENSION],
                                       float SphereVelocity[MAX_SPHERES][MAX_DIMENSION],
                                       float SphereRadius[MAX_SPHERES],
                                       float SphereTemperature[MAX_SPHERES],
                                       float SphereBeta[MAX_SPHERES],
                                       float SphereMetallicity[MAX_SPHERES],
                                       float SphereMetallicityScale[MAX_SPHERES],
                                       float CGMCentralDensity[MAX_SPHERES],
                                       float CGMCoreRadius[MAX_SPHERES],
                                       float CGMMetallicity[MAX_SPHERES],
                                       float LiveHaloMass[MAX_SPHERES],
                                       float LiveHaloScaleLength[MAX_SPHERES],
                                       float LiveHaloVirialRadius[MAX_SPHERES],
                                       float InitialTemperature,
                                       float InitialDensity,
                                       float InitialMagnField,
                                       int UseParticles,
                                       int UseGas,
                                       int UseMetals,
                                       int UseCGM,
                                       int level,
                                       int SetBaryonFields,
                                       int partitioned)
{
    /* declarations */

    int dim, field, sphere, size;
    int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
        DINum, DIINum, HDINum, MetalNum;

    /* create fields */

    int ivel, imag, ietot, ieint;
    NumberOfBaryonFields = 0;

    if (UseGas) {
      FieldType[NumberOfBaryonFields++] = Density;
      FieldType[ietot=NumberOfBaryonFields++] = TotalEnergy;
      if (DualEnergyFormalism)
        FieldType[ieint=NumberOfBaryonFields++] = InternalEnergy;
      ivel = NumberOfBaryonFields;
      FieldType[NumberOfBaryonFields++] = Velocity1;
      if (GridRank > 1)
        FieldType[NumberOfBaryonFields++] = Velocity2;
      if (GridRank > 2)
        FieldType[NumberOfBaryonFields++] = Velocity3;
      imag = NumberOfBaryonFields;
      if (HydroMethod == MHD_RK)
      {
        FieldType[NumberOfBaryonFields++] = Bfield1;
        if (GridRank > 1)
          FieldType[NumberOfBaryonFields++] = Bfield2;
        if (GridRank > 2)
          FieldType[NumberOfBaryonFields++] = Bfield3;
        FieldType[NumberOfBaryonFields++] = PhiField;
      }
      if (MultiSpecies)
      {
        FieldType[DeNum    = NumberOfBaryonFields++] = ElectronDensity;
        FieldType[HINum    = NumberOfBaryonFields++] = HIDensity;
        FieldType[HIINum   = NumberOfBaryonFields++] = HIIDensity;
        FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
        FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;
        FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
        if (MultiSpecies > 1)
        {
          FieldType[HMNum    = NumberOfBaryonFields++] = HMDensity;
          FieldType[H2INum   = NumberOfBaryonFields++] = H2IDensity;
          FieldType[H2IINum  = NumberOfBaryonFields++] = H2IIDensity;
        }
        if (MultiSpecies > 2)
        {
          FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
          FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
          FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
        }
      }
    }

    if (UseMetals)
    {
      MetalNum = NumberOfBaryonFields;
      FieldType[NumberOfBaryonFields++] = Metallicity;
    }

    int phip_num;
    if (UsePoissonDivergenceCleaning){
        FieldType[phip_num=NumberOfBaryonFields++] = Phi_pField;
        FieldType[NumberOfBaryonFields++] = DebugField;
    }

    if (WritePotential)
        FieldType[NumberOfBaryonFields++] = GravPotential;

    /* Return if this doesn't concern us. */

    if (ProcessorNumber != MyProcessorNumber)
    {
        return SUCCESS;
    }

    /* Part of Parallel Root Grid IO. The initializer is called twice.
       In its first call, the grid initializer does (almost) nothing. */

    if (SetBaryonFields == 0)
      return SUCCESS;

    float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits,
      CriticalDensity = 1, BoxLength = 1, mu = Mu;

    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
         &TimeUnits, &VelocityUnits, Time);

    float HaloMassParticles[MAX_SPHERES];

    /* =====================================================================
     * S. Selg (10/2019): N-BODY REALIZATION OF A DARK MATTER HALO.
     *
     * STAGE I: Get the number of particles in order to allocate memory.
     */

    if (UseParticles == 1 && level == 0 && partitioned == 1) // DM only
    {
        int SetupLoopCount = 0;
        int part_id = 0; // global counter
        int npart = 0;   // counter inside grid

        for (SetupLoopCount = 0; SetupLoopCount <= 1; SetupLoopCount++)
        {
            if (SetupLoopCount == 1)
            {
                /* STAGE II: Delete old particles and use count
                   from previous loop to allocate new ones.  */

                if (NumberOfParticles > 0)
                    this->DeleteParticles();

                NumberOfParticles = npart;

                /* reset counters for STAGE III */

                part_id = 0;
                npart = 0;

                if (debug && (MyProcessorNumber == ROOT_PROCESSOR))
                    printf("GalaxyLiveHalo: allocating %"ISYM" particles\n", NumberOfParticles);

                this->AllocateNewParticles(NumberOfParticles);
            }

            double XLow, XHigh, YLow, YHigh, ZLow,ZHigh;

            /* variables for IC raw data */

            double rawX;    // x-coordinate
            double rawY;    // y-coordinate
            double rawZ;    // z-coordinate
            double rawVelX; // vx
            double rawVelY; // vy
            double rawVelZ; // vz
            double rawMass; // particle mass

            /* iterate halos */

            for (sphere = 0; sphere < NumberOfSpheres; sphere++)
            {
                /* grid boundaries relative to halo center in cgs units */

                XLow  = (CellLeftEdge[0][GridStartIndex[0]] - SpherePosition[sphere][0]) *
                        LengthUnits;
                XHigh = (CellLeftEdge[0][GridEndIndex[0]] + CellWidth[0][GridEndIndex[0]] - SpherePosition[sphere][0]) *
                        LengthUnits;
                YLow  = (CellLeftEdge[1][GridStartIndex[1]] - SpherePosition[sphere][1]) *
                        LengthUnits;
                YHigh = (CellLeftEdge[1][GridEndIndex[1]] + CellWidth[1][GridEndIndex[1]] - SpherePosition[sphere][1]) *
                        LengthUnits;
                ZLow  = (CellLeftEdge[2][GridStartIndex[2]] - SpherePosition[sphere][2]) *
                        LengthUnits;
                ZHigh = (CellLeftEdge[2][GridEndIndex[2]] + CellWidth[2][GridEndIndex[2]] - SpherePosition[sphere][2]) *
                        LengthUnits;

                if (debug && SetupLoopCount == 0)
                {
                  std::cout << "GalaxyLiveHalo: grid x-boundaries (kpc): " << XLow/(1e3*pc_cm) << ", " << XHigh/(1e3*pc_cm) << std::endl;
                  std::cout << "GalaxyLiveHalo: grid y-boundaries (kpc): " << YLow/(1e3*pc_cm) << ", " << YHigh/(1e3*pc_cm) << std::endl;
                  std::cout << "GalaxyLiveHalo: grid z-boundaries (kpc): " << ZLow/(1e3*pc_cm) << ", " << ZHigh/(1e3*pc_cm) << std::endl;
                }

                /* read ICs from file */

                std::ifstream particle_file;
                particle_file.open(HaloDataFile[sphere]);

                if (particle_file.is_open())
                {
                    std::string line;

                    //while (particle_file.good())
                    while(getline(particle_file, line))
                    {
                        //std::cout << part_id << line << std::endl;

                        sscanf(line.c_str(), "%lf%lf%lf%lf%lf%lf%lf",
                            &rawX,
                            &rawY,
                            &rawZ,
                            &rawVelX,
                            &rawVelY,
                            &rawVelZ,
                            &rawMass);

                        /* convert to cgs */

                        rawX *= pc_cm;
                        rawY *= pc_cm;
                        rawZ *= pc_cm;
                        rawVelX *= 1e2;
                        rawVelY *= 1e2;
                        rawVelZ *= 1e2;
                        rawMass *= SolarMass;

                        /* check if particle is inside grid */

                        bool isInGrd = (rawX >= XLow) && (rawX <= XHigh) &&
                                       (rawY >= YLow) && (rawY <= YHigh) &&
                                       (rawZ >= ZLow) && (rawZ <= ZHigh);

                        if (isInGrd)
                        {
                            //std::cout << part_id << " is in grid\n";

                            if (SetupLoopCount == 1)
                            {
                                /* STAGE III: assign particle properties */

                                ParticleMass[npart]        = rawMass / pow(CellWidth[0][GridStartIndex[0]], 3.0) /
                                                             (DensityUnits * pow(LengthUnits, 3.0));
                                ParticleNumber[npart]      = part_id;
                                ParticleType[npart]        = PARTICLE_TYPE_DARK_MATTER;
                                ParticlePosition[0][npart] = rawX / LengthUnits + SpherePosition[sphere][0];
                                ParticlePosition[1][npart] = rawY / LengthUnits + SpherePosition[sphere][1];
                                ParticlePosition[2][npart] = rawZ / LengthUnits + SpherePosition[sphere][2];
                                ParticleVelocity[0][npart] = (rawVelX + SphereVelocity[sphere][0]) / VelocityUnits;
                                ParticleVelocity[1][npart] = (rawVelY + SphereVelocity[sphere][1]) / VelocityUnits;
                                ParticleVelocity[2][npart] = (rawVelZ + SphereVelocity[sphere][2]) / VelocityUnits;
                            }
                            npart++;
                        }
                        part_id++;

                    } // END WHILE new line

                } else {
                    fprintf(stderr,"Error in GalaxyLiveHaloInitializeGrid: cannot open particle file for halo %"ISYM"\n", sphere);
                    return FAIL;
                } // END IF file is open

                particle_file.close();

            } // END FOR spheres
        } // END FOR SetupLoopCount

        if (debug && (MyProcessorNumber == ROOT_PROCESSOR))
          printf("GalaxyLiveHalo: total number of particles: %"ISYM"\n", part_id);

        if (debug)
          printf("GalaxyLiveHalo, processor: %"ISYM": number of particles on current grid: %"ISYM"\n", MyProcessorNumber, NumberOfParticles);

        // TODO: communicate between processors for parallel root grid I/O and count for individual spheres
        HaloMassParticles[sphere] = ParticleMass[0] * part_id;
        printf("GalaxyLiveHalo: DM halo mass: %"GSYM" solar masses\n", HaloMassParticles[sphere] * DensityUnits * pow(LengthUnits, 3.0) / SolarMass);

    } // END IF UseParticles

    if (UseGas)
    {
      /* ==================================================================
       * TRANSFORM DISK COORDINATES AND INTERPOLATE TO GRID
       */

      float r_max[MAX_SPHERES], r_metal[MAX_SPHERES], r_core[MAX_SPHERES];
      float cgm_profile_b[MAX_SPHERES], cgm_profile_c[MAX_SPHERES];
      float density_min[MAX_SPHERES], density_cgm_cutoff[MAX_SPHERES], tempr_cgm_cutoff[MAX_SPHERES];
      float HaloMassGas[MAX_SPHERES];
      float SoundSpeedIsoth[MAX_SPHERES];
      float SphereRotNorm[MAX_SPHERES];
      float SphereTransformMatrix[MAX_SPHERES][MAX_DIMENSION][MAX_DIMENSION];

      /* Transform disk to grid coordinates and compute basic parameters. */

      for (sphere = 0; sphere < NumberOfSpheres; sphere++)
      {
        int iter = 0;
        double a, x_max, x_vir, tempr_vir;
        double norm;
        double z_sl[3], x_sl[3];

        SphereRotNorm[sphere] = sqrt(SphereRotAxis[sphere][0]*SphereRotAxis[sphere][0] +
                                     SphereRotAxis[sphere][1]*SphereRotAxis[sphere][1] +
                                     SphereRotAxis[sphere][2]*SphereRotAxis[sphere][2]);

        z_sl[0] = SphereRotAxis[sphere][0]/SphereRotNorm[sphere];
        z_sl[1] = SphereRotAxis[sphere][1]/SphereRotNorm[sphere];
        z_sl[2] = SphereRotAxis[sphere][2]/SphereRotNorm[sphere];

        x_sl[0] = 1.0 - z_sl[0];
        x_sl[1] = -z_sl[0];
        x_sl[2] = -z_sl[1];

        norm = sqrt(x_sl[0]*x_sl[0] + x_sl[1]*x_sl[1] + x_sl[2]*x_sl[2]);

        x_sl[0] /= norm;
        x_sl[1] /= norm;
        x_sl[2] /= norm;

        SphereTransformMatrix[sphere][0][0] = x_sl[0];
        SphereTransformMatrix[sphere][2][0] = z_sl[0];
        SphereTransformMatrix[sphere][1][0] = z_sl[1]*x_sl[2] - z_sl[2]*x_sl[1];

        SphereTransformMatrix[sphere][0][1] = x_sl[1];
        SphereTransformMatrix[sphere][2][1] = z_sl[1];
        SphereTransformMatrix[sphere][1][1] = z_sl[2]*x_sl[0] - z_sl[0]*x_sl[2];

        SphereTransformMatrix[sphere][0][2] = x_sl[2];
        SphereTransformMatrix[sphere][2][2] = z_sl[2];
        SphereTransformMatrix[sphere][1][2] = z_sl[0]*x_sl[1] - z_sl[1]*x_sl[0];

        if (debug && (MyProcessorNumber == ROOT_PROCESSOR))
          std::cout << "GalaxyLiveHalo: coordinate transformation matrix:\n"
                    << SphereTransformMatrix[sphere][0][0] << ", " << SphereTransformMatrix[sphere][0][1] << ", " << SphereTransformMatrix[sphere][0][2] << "\n"
                    << SphereTransformMatrix[sphere][1][0] << ", " << SphereTransformMatrix[sphere][1][1] << ", " << SphereTransformMatrix[sphere][1][2] << "\n"
                    << SphereTransformMatrix[sphere][2][0] << ", " << SphereTransformMatrix[sphere][2][1] << ", " << SphereTransformMatrix[sphere][2][2] << std::endl;

        
        // cylindrical/spherical region in which disk data are interpolated
        r_max[sphere] = SphereRadius[sphere] * kpc_cm/LengthUnits;

        // metallicity scale length in code units
        r_metal[sphere] = SphereMetallicityScale[sphere] * kpc_cm/LengthUnits;

        if (UseCGM) {
          // CGM core radius in code units
          r_core[sphere] = CGMCoreRadius[sphere] * kpc_cm/LengthUnits;

          a = LiveHaloScaleLength[sphere]/CGMCoreRadius[sphere];
          x_max = r_max[sphere]/r_core[sphere];
          x_vir = LiveHaloVirialRadius[sphere]/CGMCoreRadius[sphere];

          // check if temperature drops below zero with some margin at sphere radius 
          if (x_max >= x_vir)
            ENZO_FAIL("Error in GalaxyLiveHalo: unphysical CGM profile, sphere radius must be smaller than virial radius.")

          // density and temperature at sphere boundary assuming pressure balance with ambient medium
          density_cgm_cutoff[sphere] = CGMCentralDensity[sphere] / (1 + x_max*x_max);
          
          // coefficients b and c of CGM temperature profile
          cgm_profile_b[sphere] = GravConst * mu * mh * LiveHaloMass[sphere] * SolarMass / (kboltz * CGMCoreRadius[sphere] * kpc_cm);          
          cgm_profile_c[sphere] = -cgm_profile_b[sphere] * (1 + a*a - (a*a - 1) * (a + x_vir) * atan(x_vir) 
                                                            + a*(a + x_vir) * (log(x_vir*x_vir + 1) - 2 * log(a + x_vir))) 
                                  / (pow((a*a + 1), 2) * (a + x_vir));
          /*             
          cgm_profile_c[sphere] = CGMCentralTemperature[sphere] -
                                  cgm_profile_b[sphere] * (1 + a*a * (1.0 - 2.0 * log(a))) / (pow((a*a + 1), 2) * a);
          */
          
          tempr_cgm_cutoff[sphere] = 
            compute_tempr_cgm(r_max[sphere]/r_core[sphere],
                              LiveHaloScaleLength[sphere]/CGMCoreRadius[sphere],
                              cgm_profile_b[sphere],
                              cgm_profile_c[sphere]);

          // check if temperature drops below zero with some margin at sphere radius 
          if (tempr_cgm_cutoff[sphere] < ShockTemperatureFloor)
            ENZO_FAIL("Error in GalaxyLiveHalo: unphysical CGM profile, temperature must be positive inside sphere.")

        } else {
          // density floor to ensure that disk pressure exceeds pressure of ambient medium
          density_min[sphere] = InitialDensity * InitialTemperature/SphereTemperature[sphere];
        }

        // isothermal speed of sound (pressure/density)
        SoundSpeedIsoth[sphere] = sqrt(kboltz * SphereTemperature[sphere] / (mu * mh)) / VelocityUnits;

        if (MyProcessorNumber == ROOT_PROCESSOR) {
          printf("GalaxyLiveHalo: disk isothermal speed of sound: %"GSYM" cm/s\n", SoundSpeedIsoth[sphere] * VelocityUnits);
          if (UseMetals)
            printf("GalaxyLiveHalo: disk metallicity scale length: %"GSYM" kpc, %"GSYM"\n", SphereMetallicityScale[sphere], r_metal[sphere]);
          printf("GalaxyLiveHalo: sphere radius: %"GSYM" kpc, %"GSYM"\n", SphereRadius[sphere], r_max[sphere]);
          if (UseCGM) {
            printf("GalaxyLiveHalo: CGM core radius: %"GSYM" kpc, %"GSYM"\n", CGMCoreRadius[sphere], r_core[sphere]);
            printf("GalaxyLiveHalo: CGM total mass: %"GSYM" solar masses\n", 
              4.0 * pi * pow(r_core[sphere] * LengthUnits, 3) * CGMCentralDensity[sphere] * DensityUnits 
              * (r_max[sphere]/r_core[sphere] - atan(r_max[sphere]/r_core[sphere])) / SolarMass);
            printf("GalaxyLiveHalo: CGM profile coefficients: %"GSYM", %"GSYM"\n", cgm_profile_b[sphere], cgm_profile_c[sphere]);
            printf("GalaxyLiveHalo: CGM central density: %"GSYM" g/cm^3\n", CGMCentralDensity[sphere] * DensityUnits);
            printf("GalaxyLiveHalo: CGM central temperature: %"GSYM" K\n", 
              compute_tempr_cgm(0,
                                LiveHaloScaleLength[sphere]/CGMCoreRadius[sphere],
                                cgm_profile_b[sphere],
                                cgm_profile_c[sphere]));
            printf("GalaxyLiveHalo: CGM virial temperature: %"GSYM" K \n", 2.0 * cgm_profile_b[sphere] / (5.0 * x_vir));
            printf("GalaxyLiveHalo: CGM density at cutoff: %"GSYM" g/cm^3, %"GSYM"\n", 
              density_cgm_cutoff[sphere] * DensityUnits, density_cgm_cutoff[sphere]);
            printf("GalaxyLiveHalo: CGM temperature at cutoff: %"GSYM" K\n", tempr_cgm_cutoff[sphere]);
          } else {
            printf("GalaxyLiveHalo: sphere density floor: %"GSYM" g/cm^3, %"GSYM"\n", density_min[sphere] * DensityUnits, density_min[sphere]);
          }
          printf("GalaxyLiveHalo: background density: %"GSYM" g/cm^3 %"GSYM"\n", InitialDensity * DensityUnits, InitialDensity);
          printf("GalaxyLiveHalo: background temperature: %"GSYM" K\n", InitialTemperature);
          printf("GalaxyLiveHalo: background magnetic field: %"GSYM" G\n",
                 sqrt(4 * pi * DensityUnits) * VelocityUnits * InitialMagnField);
        }
      } // END FOR sphere

      /* Set up the baryon field. */

      this->AllocateGrids();

      /* Loop over the grid. */

      int i, j, k;
      int n = 0;

      bool isInSphere;

      double x, y, z;
      double xpos, ypos, zpos;
      double r, rcyl, xcyl, ycyl, zcyl;
      double density, temperature, metallicity;
      double density_cgm, temperature_cgm;
      double norm, damping, vphi, Bphi;
      double Toroidal[MAX_DIMENSION], Velocity[MAX_DIMENSION], MagnField[MAX_DIMENSION];

      for (k = 0; k < GridDimension[2]; k++)
        for (j = 0; j < GridDimension[1]; j++)
          for (i = 0; i < GridDimension[0]; i++, n++)
          {
            x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
            y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
            z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

            isInSphere = false;

            //vphi = 0.0;
            density = InitialDensity;
            temperature = InitialTemperature;
            metallicity = MetallicityFloor;

            for (dim = 0; dim < MAX_DIMENSION; dim++)
            {
              Velocity[dim]  = 0.0e0;
              MagnField[dim] = 0.0e0;
              Toroidal[dim]  = 0.0e0;
            }

            norm = 1.0; // to avoid division by zero

            for (sphere = 0; sphere < NumberOfSpheres; sphere++)
            {
              // displacement from sphere center
              xpos = x - SpherePosition[sphere][0];
              ypos = y - SpherePosition[sphere][1];
              zpos = z - SpherePosition[sphere][2];

              // radial distance from center
              r = sqrt(pow(fabs(xpos), 2) + pow(fabs(ypos), 2) + pow(fabs(zpos), 2));

              if (r < r_max[sphere])
              {
                /* Terminate if spheres are not disjunct. */

                if (isInSphere)
                    ENZO_VFAIL("Error in Grid::GalaxyLiveHaloInitializeGrid: sphere %"ISYM" overlaps with other sphere", sphere);

                isInSphere = true; // raise the flag

                /* Transform to disk coordinate system (Cartesian/cylindrical). */

                xcyl = xpos*SphereTransformMatrix[sphere][0][0] +
                       ypos*SphereTransformMatrix[sphere][0][1] +
                       zpos*SphereTransformMatrix[sphere][0][2];

                ycyl = xpos*SphereTransformMatrix[sphere][1][0] +
                       ypos*SphereTransformMatrix[sphere][1][1] +
                       zpos*SphereTransformMatrix[sphere][1][2];

                zcyl = xpos*SphereTransformMatrix[sphere][2][0] +
                       ypos*SphereTransformMatrix[sphere][2][1] +
                       zpos*SphereTransformMatrix[sphere][2][2];

                // radial coordinate in central plane of the disk
                rcyl = sqrt(xcyl*xcyl + ycyl*ycyl);

                /* Interpolate disk data. */

                density = DiskTable[sphere].InterpolateEquilibriumDensityTable(rcyl*LengthUnits, zcyl*LengthUnits) / DensityUnits;
                vphi    = DiskTable[sphere].InterpolateEquilibriumVcircTable(rcyl*LengthUnits, zcyl*LengthUnits) / VelocityUnits;
                Bphi    = 0.0;
                damping = 1.0;
                
                if (UseCGM) 
                {                  
                  /* CGM model */

                  density_cgm = CGMCentralDensity[sphere] / (1 + pow(r/r_core[sphere], 2));
                  temperature_cgm = 
                    compute_tempr_cgm(r/r_core[sphere],
                                      LiveHaloScaleLength[sphere]/CGMCoreRadius[sphere],
                                      cgm_profile_b[sphere],
                                      cgm_profile_c[sphere]);
                /*
                  (1 + pow(r/r_core[sphere], 2))
                    * (HaloMassDM * SolarMass * F0(r, r_core[sphere], LiveHaloScaleLength[sphere] * kpc_cm/LengthUnits)
                       + 4 * pi * pow(CGMCoreRadius[sphere] * kpc_cm, 3) * CGMCentralDensity[sphere] * DensityUnits * F1(r, r_core[sphere]))
                    / LengthUnits;
                */

                  /* Gas pressure inside disk exceeds pressure of ambient medium and CGM. */               

                  if ((density > density_min[sphere]) && 
                      (density * SphereTemperature[sphere] > density_cgm * temperature_cgm))
                  {
                      /* Set disk temperature, metallicity, and magnetic field. */

                      temperature = SphereTemperature[sphere];
                      metallicity = max(SphereMetallicity[sphere] * exp(-rcyl/r_metal[sphere]), CGMMetallicity[sphere]);
                      Bphi = sqrt(2 * density / SphereBeta[sphere]) * SoundSpeedIsoth[sphere];
                  }
                  else
                  {
                    /* Set CGM density and temperature and exponential damping of rotation velocity if                      
                       pressure exceeds pressure of ambient medium. */                    

                    if (density_cgm * temperature_cgm > InitialDensity * InitialTemperature)
                    //if (density_cgm > InitialDensity)
                    {
                      density = density_cgm;
                      temperature = temperature_cgm;
                      metallicity = CGMMetallicity[sphere];
                    }
                    else
                    {
                      density = InitialDensity;
                    }

                    /* CGM and ambient medium are non-rotating. */

                    vphi = 0.0;
                  }
                }
                else
                {
                  /* Disk interior if gas pressure exceeds ambient pressure. */
                 
                  if (density > density_min[sphere])
                  {
                    /* Set disk temperature, metallicity, and magnetic field. */

                    temperature = SphereTemperature[sphere];
                    metallicity = max(SphereMetallicity[sphere] * exp(-rcyl/r_metal[sphere]), MetallicityFloor);
                    Bphi = sqrt(2 * density / SphereBeta[sphere]) * SoundSpeedIsoth[sphere];
                  }
                  else
                  {
                    /* Set density floor and exponential damping of rotation velocity. */
 
                    damping = exp(1 - density_min[sphere]/density);
                    density = InitialDensity;
                  }
                }

                /* Consistency check while developing. */

                //if (density * temperature < InitialDensity * InitialTemperature)
                //  ENZO_FAIL("Error in GalaxyLiveHalo: pressure below floor, debug code.")

                //if ((density * temperature <= InitialDensity * InitialTemperature) && (vphi*damping > 0.0))
                //  printf("r = %"GSYM", rcyl = %"GSYM", z = %"GSYM", vphi = %"GSYM", %"GSYM"\n", r, rcyl, z, vphi*damping, density * temperature);

                /* Define azimuthal unit vector if not at disk center. */

                if (r > 0.5*CellWidth[0][i]) 
                {
                  
                  Toroidal[0] = (SphereRotAxis[sphere][1]*zpos - SphereRotAxis[sphere][2]*ypos);
                  Toroidal[1] = (SphereRotAxis[sphere][2]*xpos - SphereRotAxis[sphere][0]*zpos);
                  Toroidal[2] = (SphereRotAxis[sphere][0]*ypos - SphereRotAxis[sphere][1]*xpos);

                  norm = sqrt(Toroidal[0]*Toroidal[0] + Toroidal[1]*Toroidal[1] + Toroidal[2]*Toroidal[2]);
                }

                /* Set rotation plus center-of-mass velocity and toroidal magnetic field. */

                for (dim = 0; dim < MAX_DIMENSION; dim++)
                {
                  Toroidal[dim] /= norm;
                  Velocity[dim]  = (vphi*Toroidal[dim] + SphereVelocity[sphere][dim] / VelocityUnits) * damping;
                  MagnField[dim] =  Bphi*Toroidal[dim];
                }
              } // END IF inside sphere
            } // END FOR sphere

            /* Set density. */

            BaryonField[iden][n] = density;

            /* If doing multi-species (HI, etc.), set these. */

            if (MultiSpecies > 0) {

              BaryonField[HIINum][n] = CosmologySimulationInitialFractionHII *
                CoolData.HydrogenFractionByMass * BaryonField[0][n];
              BaryonField[HeIINum][n] = CosmologySimulationInitialFractionHeII*
                BaryonField[0][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
              BaryonField[HeIIINum][n] = CosmologySimulationInitialFractionHeIII*
                BaryonField[0][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
              BaryonField[HeINum][n] =
                (1.0 - CoolData.HydrogenFractionByMass)*BaryonField[0][n] -
                BaryonField[HeIINum][n] - BaryonField[HeIIINum][n];

              if (MultiSpecies > 1) {
                BaryonField[HMNum][n] = CosmologySimulationInitialFractionHM*
                  BaryonField[HIINum][n]* pow(temperature,float(0.88));
                BaryonField[H2IINum][n] = CosmologySimulationInitialFractionH2II*
                  2.0*BaryonField[HIINum][n]* pow(temperature,float(1.8));
                BaryonField[H2INum][n] = CosmologySimulationInitialFractionH2I*
                  BaryonField[0][n]*CoolData.HydrogenFractionByMass*pow(301.0,5.1);
              }

              BaryonField[HINum][n] =
                CoolData.HydrogenFractionByMass*BaryonField[0][n]
                - BaryonField[HIINum][n];
              if (MultiSpecies > 1)
                BaryonField[HINum][n] -= BaryonField[HMNum][n]
                  + BaryonField[H2IINum][n] + BaryonField[H2INum][n];

              BaryonField[DeNum][n] = BaryonField[HIINum][n] +
                0.25*BaryonField[HeIINum][n] + 0.5*BaryonField[HeIIINum][n];
              if (MultiSpecies > 1)
                BaryonField[DeNum][n] += 0.5*BaryonField[H2IINum][n] -
                  BaryonField[HMNum][n];

              /* Set Deuterium species (assumed to be negligible). */

              if (MultiSpecies > 2) {
                BaryonField[DINum][n] = CoolData.DeuteriumToHydrogenRatio*
                  BaryonField[HINum][n];
                BaryonField[DIINum][n] = CoolData.DeuteriumToHydrogenRatio*
                  BaryonField[HIINum][n];
                BaryonField[HDINum][n] = CoolData.DeuteriumToHydrogenRatio*
              BaryonField[H2INum][n];
              }
            } // END IF MultiSpecies

            if (UseMetals)
                BaryonField[MetalNum][n] = metallicity*density;

            /* Set Velocities. */

            for (int dim = 0; dim < GridRank; dim++)
              BaryonField[ivel+dim][n] = Velocity[dim];

            /* Set energy (thermal and then total if necessary) and magnetic field. */

            BaryonField[ietot][n] = temperature/TemperatureUnits / ((Gamma-1.0)*mu);

            if (DualEnergyFormalism)
              BaryonField[ieint][n] = BaryonField[ietot][n];

            for (dim = 0; dim < GridRank; dim++)
              BaryonField[ietot][n] += 0.5*pow(BaryonField[ivel+dim][n], 2);

            if (HydroMethod == MHD_RK) {

              // initialize uniform field in z-direction

              BaryonField[iBx][n] = 0.0;
              BaryonField[iBy][n] = 0.0;
              BaryonField[iBz][n] = InitialMagnField;
              BaryonField[iPhi][n] = 0.0;

              // add toroidal field and update total energy

              for (int dim = 0; dim < GridRank; dim++) {
                BaryonField[imag+dim][n] += MagnField[dim];
                BaryonField[ietot][n] += 0.5*pow(BaryonField[imag+dim][n], 2)/BaryonField[0][n];
              }
            }
          } // END loop over grid
    } // END FOR UseGas

    return SUCCESS;
}

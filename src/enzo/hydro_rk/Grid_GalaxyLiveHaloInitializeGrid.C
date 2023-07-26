/***********************************************************************
/
/  INITIALIZE DISK GALAXIES WITH LIVE HALOS
/
/  written by: Kai Rodenbeck and Simon Selg
/  date:       July, 2020
/
/  PURPOSE:
/    Reads initial conditions from files and 
/    sets up one or more disk galaxies and particle halos.
/
************************************************************************/

// S. Selg (10/2019) Extra includes for file reading 
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
// ---

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
#include "HydrostaticDisk.h"

#define NTHETA 1000
#define NR 1000

/********************* PROTOTYPES *********************/

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float Time);
int CosmologyComputeExpansionFactor(float time, float *a, float *dadt);

float gasdev();

double BESSI0(double y);
double BESSI1(double y);
double BESSK0(double y);
double BESSK1(double y);

/* Set various units. */

// Used to compute Bonner-Ebert density profile
//double BE(double r);
//double q(double r);
//double Ang(double a1, double a2, double R, double r);

// Returns random velocity from Maxwellian distribution
double Maxwellian(double c_tilda, double vel_unit, double mu, double gamma);
//double ERF(double x);

//int ComputeRadialVelocity(float density, double mass, float r_init, 
//			  double VelocityUnits, double LengthUnits,
//			  float dm_velrot_corr, double radius_vr[], double Vr[], 
//			  double exterior_rho[], int Npts);


/*******************************************************/

static int CollapseTestParticleCount = 0;

static float CosmologySimulationOmegaBaryonNow       = 0.0463;
static float CosmologySimulationInitialFractionHII   = 1.2e-5;
static float CosmologySimulationInitialFractionHeII  = 1.0e-14;
static float CosmologySimulationInitialFractionHeIII = 1.0e-17;
static float CosmologySimulationInitialFractionHM    = 2.0e-9;
static float CosmologySimulationInitialFractionH2I   = 2.0e-20;
static float CosmologySimulationInitialFractionH2II  = 3.0e-14;

static int legacy = 0;

// structure for legacy model test
struct galaxy_parameters {
  int N_x, N_z;
  
  // disk scales in cgs units
  double r_s, z_s;
  double dx, dz;
  double sigma;
  double v_s;

  // halo scales in cgs units
  double halo_mass, halo_scale;
};


int grid::GalaxyLiveHaloInitializeGrid(int NumberOfSpheres,
				       EquilibriumGalaxyDisk DiskTable[MAX_SPHERES],
				       float SpherePosition[MAX_SPHERES][MAX_DIMENSION],
				       float SphereAxisRot[MAX_SPHERES][MAX_DIMENSION],
				       float SphereVelocity[MAX_SPHERES][MAX_DIMENSION],
				       float SphereRadius[MAX_SPHERES],
				       float SphereTemperature[MAX_SPHERES],
				       float SphereMetallicity[MAX_SPHERES], // dummy argument
				       float InitialTemperature,
				       float InitialDensity,
				       float InitialMagnField,
				       int UseParticles,
				       int level,
				       int SetBaryonFields,
				       int partitioned,
				       float TopGridSpacing,
				       int maxlevel,
				       int grid_traditional,
				       float grid_safety_factor)
{
	/* declarations */
	
	int dim, field, sphere, size;
	int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
		DINum, DIINum, HDINum, MetalNum;

	float ParticleMeanDensity = FLOAT_UNDEFINED;
	
	/* create fields */
	
	NumberOfBaryonFields = 0;
	FieldType[NumberOfBaryonFields++] = Density;
	FieldType[NumberOfBaryonFields++] = TotalEnergy;
	if (DualEnergyFormalism)
		FieldType[NumberOfBaryonFields++] = InternalEnergy;
	int ivel = NumberOfBaryonFields;
	FieldType[NumberOfBaryonFields++] = Velocity1;
	if (GridRank > 1) 
		FieldType[NumberOfBaryonFields++] = Velocity2;
	if (GridRank > 2)
		FieldType[NumberOfBaryonFields++] = Velocity3;
	int imag = NumberOfBaryonFields;
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
	
	// S. Selg (08/2019): Toggle output of gravitational potential
	if (WritePotential)
		FieldType[NumberOfBaryonFields++] = GravPotential;

	/* Return if this doesn't concern us. */
	
	if (ProcessorNumber != MyProcessorNumber)
	{
		return SUCCESS;
	}

	/*====================================================================
	 * S. Selg (09/2019): Part of Parallel Root Grid IO. The initializer is
	 * called twice. In its first call, the grid initializer does (almost)
	 * nothing.
	 *===================================================================*/
	if (SetBaryonFields == 0)
			      return SUCCESS;
	
	float	DensityUnits,
			LengthUnits,
			TemperatureUnits,
			TimeUnits,
			VelocityUnits,
			CriticalDensity = 1,
			BoxLength = 1,
			mu = Mu;
	
	double MyGrav=GravConst;
	
	float a, dadt, ExpansionFactor = 1;
	GetUnits(	&DensityUnits,
				&LengthUnits,
				&TemperatureUnits,
				&TimeUnits,
				&VelocityUnits,
				Time);
	if (ComovingCoordinates)
	{
		CosmologyComputeExpansionFactor(Time, &a, &dadt);
		ExpansionFactor = a/(1.0+InitialRedshift);
		CriticalDensity = 2.78e11*pow(HubbleConstantNow, 2); // in Msolar/Mpc^3
		BoxLength = ComovingBoxSize*ExpansionFactor/HubbleConstantNow;  // in Mpc
	}
	else
	{
		CriticalDensity = 2.78e11*pow(0.74,2); // in Msolar/Mpc^3 for h=0.74
		BoxLength = LengthUnits / 3.086e24;
		HubbleConstantNow = 1.0;
		OmegaMatterNow = 1.0;
		MyGrav=4*pi*GravConst*DensityUnits*TimeUnits*TimeUnits;
	}
	
	double MagnConversionFactor=sqrt(4.0*3.14159/(VelocityUnits*VelocityUnits)/DensityUnits);
	
	/* =====================================================================
	 * S. Selg (10/2019): N-BODY REALIZATION OF A DARK MATTER HALO. 
	 *
	 * STAGE I: Get the number of particles in order to allocate memory.
	 */

	if (UseParticles == 1 && level == 0 && partitioned == 1) // > 1 would not be specific 
	{		             		   // since the following 
						   // refers to DM
		int ParticleLoopCount = 0;
		int npart = 0;
		for (ParticleLoopCount = 0; ParticleLoopCount < 2; ParticleLoopCount++)
		{
			if (ParticleLoopCount == 1)
			{
				// STAGE II: Delete old particles
				if (NumberOfParticles > 0)
					this->DeleteParticles();
				NumberOfParticles = npart;
				npart = 0;
				// STAGE III: Allocate new ones
				this->AllocateNewParticles(NumberOfParticles);
			}	
		

			// GET CELL GRID FACE COORDINATES
			double grdXLow  = CellLeftEdge[0][GridStartIndex[0]];
			double grdXHigh = CellLeftEdge[0][GridEndIndex[0]] + 
				CellWidth[0][GridEndIndex[0]];
			double grdYLow  = CellLeftEdge[1][GridStartIndex[1]];
			double grdYHigh = CellLeftEdge[1][GridEndIndex[1]] +
				CellWidth[1][GridEndIndex[1]];
			double grdZLow  = CellLeftEdge[2][GridStartIndex[2]];
			double grdZHigh = CellLeftEdge[2][GridEndIndex[2]] +
				CellWidth[2][GridEndIndex[2]];
			double preCompX;	// x-coordinate
			double preCompY;        // y-coordinate
			double preCompZ;        // z-coordinate
			double preCompVelX;     // vx
			double preCompVelY;  	// vy
			double preCompVelZ;     // vz
			double preCompMass;     // particle mass
			double preCompDist;     // particle distance towards center
					// of halo
			double particle_label;    // just a number, each particle has its own

			// ITERATE HALOS
			for (sphere = 0; sphere < NumberOfSpheres; sphere++)
			{
				// WE READ THE PARTICLE ICs FROM FILE
				std::ifstream particle_file;
				if (sphere == 0)
					particle_file.open("particle_ic_a");
				else if (sphere == 1)
					particle_file.open("particle_ic_b");
				else
				{
					fprintf(stderr,"Error in GalaxyLiveHaloInitializeGrid: halo %"ISYM" cannot be initialized. Missing file!\n", sphere);
					return FAIL;
				}
				if (particle_file.is_open())
				{
					std::string line;
					getline(particle_file, line); // 1st line
			        	while (particle_file.good())
					{
						// WE'LL ASSING VARIABLES TO PARTICLE
						// DATA INFERRED FROM FILE
						// READ FILE INTO MEMORY
						sscanf(line.c_str(), "%lf%lf%lf%lf%lf%lf%lf",
							&preCompX,
							&preCompY,
							&preCompZ,
							&preCompVelX,
							&preCompVelY,
							&preCompVelZ,
						        &preCompMass);
						        //&preCompDist,
							//&particle_label);

						// convert to cgs
						preCompX *= pc_cm;
 						preCompY *= pc_cm;
 						preCompZ *= pc_cm;						
						preCompVelX *= 1e2; //1e5;
						preCompVelY *= 1e2; //1e5;
						preCompVelZ *= 1e2; //1e5;
						preCompMass *= SolarMass;

						// CHECK IF THE PARTICLE FITS ONTO THE
						// GRID
						bool isInGrd = (preCompX / LengthUnits + SpherePosition[sphere][0] >= grdXLow) &&
							(preCompX / LengthUnits + SpherePosition[sphere][0] <= grdXHigh) &&
							(preCompY / LengthUnits + SpherePosition[sphere][1] >= grdYLow) &&
							(preCompY / LengthUnits + SpherePosition[sphere][1] <= grdYHigh) &&
							(preCompZ / LengthUnits + SpherePosition[sphere][2] >= grdZLow) &&
							(preCompZ / LengthUnits + SpherePosition[sphere][2] <= grdZHigh);
						// IF IT FITS, I SITS :)
						if (isInGrd)
						{
							if (ParticleLoopCount == 1)
							{
								// STAGE IV: assign particle properties
								ParticleMass[npart]        = preCompMass / pow(CellWidth[0][GridStartIndex[0]], 3.0) / (DensityUnits * pow(LengthUnits, 3.0));
								ParticleNumber[npart]      = particle_label;
								ParticleType[npart]        = PARTICLE_TYPE_DARK_MATTER;
								ParticlePosition[0][npart] = preCompX / LengthUnits + SpherePosition[sphere][0];
								ParticlePosition[1][npart] = preCompY / LengthUnits + SpherePosition[sphere][1];
								ParticlePosition[2][npart] = preCompZ / LengthUnits + SpherePosition[sphere][2];
								ParticleVelocity[0][npart] = preCompVelX / VelocityUnits + SphereVelocity[sphere][0];
								ParticleVelocity[1][npart] = preCompVelY / VelocityUnits + SphereVelocity[sphere][1];
								ParticleVelocity[2][npart] = preCompVelZ / VelocityUnits + SphereVelocity[sphere][2];
							}
							npart++;

						}
					
						getline(particle_file, line); // advance
								      // by one
								      // line
					} // CLOSING WHILE	
				
				} else {
					fprintf(stderr, "Error in MHDGalaxyDiskInializeGrid. File for halo no. %"ISYM" not found.\n", sphere);
					return FAIL;
				}// CLOSING IF file is open
				particle_file.close();
			} // CLOSING FOR spheres
		} // CLOSING FOR ParticleCountLoop

	} // CLOSING IF DM particles
	/* Set densities */
		
	float BaryonMeanDensity, ParticleCount = 0;
	switch (UseParticles)
	{
		case 1:
			BaryonMeanDensity = CosmologySimulationOmegaBaryonNow / OmegaMatterNow;
			break;
		case 2:
			BaryonMeanDensity = 1 - CosmologySimulationOmegaBaryonNow / OmegaMatterNow;
			break;
		default:
			BaryonMeanDensity = 1.0;
	} // ENDSWITCH UseParticles
	
	
	BaryonMeanDensity = 1.0;
	
	if (ParticleMeanDensity == FLOAT_UNDEFINED)
		ParticleMeanDensity = 1.0 - BaryonMeanDensity;
	else
		BaryonMeanDensity = 1.0 - ParticleMeanDensity;
	
	/* Set up the baryon field. */
	
	/* compute size of fields */
	
	size = 1;
	for (int dim = 0; dim < GridRank; dim++)
		size *= GridDimension[dim];
	
	/* allocate fields */
		
	for (field = 0; field < NumberOfBaryonFields; field++)
		if (BaryonField[field] == NULL)
			BaryonField[field] = new float[size];
	
	/* Loop over the mesh. */
	
	float r_max[MAX_SPHERES];
	float SoundSpeed[MAX_SPHERES];
	float SphereTransformMatrix[MAX_SPHERES][MAX_DIMENSION][MAX_DIMENSION];

	double poisson_coeff = 4*pi*GravConst; // cgs units

	/* legacy method */
	hydrostatic_disk Galaxy[MAX_SPHERES];
	
	/* test disk parameters in cgs units
           - disk parameters from enzo-e/input/IsolatedGalaxy 
           - Hernquist halo: M_vir = 3e10 M_sun, r_vir = 65.57 kpc */

	galaxy_parameters test;

	if (legacy) {
	  test.N_x = 101;
	  test.N_z = 144;

	  test.r_s = 1.2856989923909257e+22;
	  test.z_s = 8.999892947045047e+21;
	  test.dx = 0.1*test.r_s;
	  test.dz = 0.1*test.z_s;

	  // disk mass = 1.3913997490703146e+41 g
	  test.sigma = 1.3396564680155456e-4;
	  test.v_s = sqrt(poisson_coeff * test.sigma * test.r_s);

	  test.halo_mass = 8.008978055391178e+43;
	  test.halo_scale = 3.2109177913790463e+22;
	}

	/* Pre-compute cloud properties before looping over mesh */
	for (sphere = 0; sphere < NumberOfSpheres; sphere++)
	{
		if(MAX_DIMENSION==3 && !(SphereAxisRot[sphere][2]==0.0 && SphereAxisRot[sphere][1]==0))
		{
			double Ang_Vel=sqrt(SphereAxisRot[sphere][0]*SphereAxisRot[sphere][0]
							+	SphereAxisRot[sphere][1]*SphereAxisRot[sphere][1]
							+	SphereAxisRot[sphere][2]*SphereAxisRot[sphere][2]);
			double z_sl[3],x_sl[3];
			
			z_sl[0]=SphereAxisRot[sphere][0]/Ang_Vel;
			z_sl[1]=SphereAxisRot[sphere][1]/Ang_Vel;
			z_sl[2]=SphereAxisRot[sphere][2]/Ang_Vel;
			
			for(int it=0;it<MAX_DIMENSION;it++)
				SphereAxisRot[sphere][it]/=Ang_Vel;
			
			
			x_sl[0]=1.0-z_sl[0];
			x_sl[1]=-z_sl[0];
			x_sl[2]=-z_sl[1];
			
			double temp_norm=sqrt(x_sl[0]*x_sl[0]+x_sl[1]*x_sl[1]+x_sl[2]*x_sl[2]);
			x_sl[0]/=temp_norm;
			x_sl[1]/=temp_norm;
			x_sl[2]/=temp_norm;
			
			SphereTransformMatrix[sphere][0][0]=x_sl[0];
			SphereTransformMatrix[sphere][2][0]=z_sl[0];
			SphereTransformMatrix[sphere][1][0]=z_sl[1]*x_sl[2]-z_sl[2]*x_sl[1];
			
			SphereTransformMatrix[sphere][0][1]=x_sl[1];
			SphereTransformMatrix[sphere][2][1]=z_sl[1];
			SphereTransformMatrix[sphere][1][1]=z_sl[2]*x_sl[0]-z_sl[0]*x_sl[2];
			
			SphereTransformMatrix[sphere][0][2]=x_sl[2];
			SphereTransformMatrix[sphere][2][2]=z_sl[2];
			SphereTransformMatrix[sphere][1][2]=z_sl[0]*x_sl[1]-z_sl[1]*x_sl[0];
		}

		r_max[sphere] = SphereRadius[sphere] * kpc_cm/LengthUnits;
		
		SoundSpeed[sphere] = sqrt(kboltz * SphereTemperature[sphere] / (mu * mh)) 
		                     / VelocityUnits;
		
		if(debug && (MyProcessorNumber == ROOT_PROCESSOR)) {
		  printf("Sphere radius: %"GSYM" kpc, %"GSYM"\n", SphereRadius[sphere], r_max[sphere]);
		  printf("Background density: %"GSYM" g/cm^3\n", InitialDensity * DensityUnits);
		  printf("Isothermal speed of sound: %"GSYM" cm/s\n", SoundSpeed[sphere] * VelocityUnits);
		  if (legacy) 
		    printf("Central disk velocity =  %"GSYM"cm/s\n", test.v_s);
		}

		/* compute disk data using legacy method */

		if (legacy) {
		  Galaxy[sphere] = hydrostatic_disk(test.N_x, test.N_z, test.dx, test.dz, test.r_s, test.z_s, 
						    test.sigma, SoundSpeed[sphere] * VelocityUnits);

		  Galaxy[sphere].set_halo_mass(test.halo_mass);
		  Galaxy[sphere].set_halo_scale(test.halo_scale);

		  double err=1.0;	
		  int counter=0;
		
		  while( err > 1.0e-6 && counter < 100 )
		  {
		    err=Galaxy[sphere].integrate();
		    counter++;
		    if(MyProcessorNumber == ROOT_PROCESSOR)
		      printf("\nLegacy model sphere %i: %ith iteration, error=%e\n", sphere, counter, err);
		  }

		  Galaxy[sphere].rotation_curve();
		}

	}// ENDFOR sphere

	int i, j, k;
	int n = 0;

	float x, y, z;
	float density, temperature;
	float Velocity[MAX_DIMENSION], MagnField[MAX_DIMENSION];

	float r, rcyl, xpos, ypos, zpos;
	float cosphi, sinphi, sintheta, vphi;
	float RotVelocity[MAX_DIMENSION];

	/*
	if((MyProcessorNumber == ROOT_PROCESSOR))
	  for (i = 0; i < GridDimension[0]; i++)
	    printf("Delta x: %"GSYM" pc\n", CellWidth[0][i]*LengthUnits/pc_cm);
	*/

	for (k = 0; k < GridDimension[2]; k++)
	  for (j = 0; j < GridDimension[1]; j++)
	    for (i = 0; i < GridDimension[0]; i++, n++) {
	      
	      x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	      y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	      z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

	      density = InitialDensity;
	      temperature = InitialTemperature;
	      for (dim = 0; dim < MAX_DIMENSION; dim++) {
		Velocity[dim] = 0;
		MagnField[dim] = 0;
	      }

	      for (sphere = 0; sphere < NumberOfSpheres; sphere++) {
	    
		// radial distance from center
		r = sqrt(pow(fabs(x-SpherePosition[sphere][0]), 2) +
			 pow(fabs(y-SpherePosition[sphere][1]), 2) +
			 pow(fabs(z-SpherePosition[sphere][2]), 2) );
		r = max(r, 0.1*CellWidth[0][0]);
	    
		if (r < r_max[sphere]) {
		  xpos = x - SpherePosition[sphere][0];
		  ypos = y - SpherePosition[sphere][1];
		  zpos = z - SpherePosition[sphere][2];
	      
                  // radial coordinate in central plane of the disk
		  rcyl = sqrt(xpos*xpos + ypos*ypos);

		  // azimuthal angle
		  cosphi = xpos/sqrt(xpos*xpos + ypos*ypos);
		  sinphi = ypos/sqrt(xpos*xpos + ypos*ypos);

		  // polar angle
		  sintheta = sqrt(xpos*xpos + ypos*ypos)/sqrt(xpos*xpos + ypos*ypos + zpos*zpos);

		  // interpolate disk data
		  if (legacy)
		  {
		    double v_g_sqr, v_p_sqr, v_dm_sqr;
		    double r_n = 0.5*rcyl*LengthUnits / (test.r_s); // normalized radial coordinate
		    
		    density = max(InitialDensity,
                                  Galaxy[sphere].intpl_rho(rcyl*LengthUnits, fabs(zpos)*LengthUnits) / DensityUnits);
		    
		    v_g_sqr  = pow(test.v_s * r_n, 2) * (BESSI0(r_n) * BESSK0(r_n) - BESSI1(r_n) * BESSK1(r_n)) 
		              / (VelocityUnits*VelocityUnits);
		    v_p_sqr  = Galaxy[sphere].intpl_v_sqr(rcyl*LengthUnits, fabs(zpos)*LengthUnits) 
		              / (VelocityUnits*VelocityUnits);
		    v_dm_sqr = Galaxy[sphere].halo_vel_sqr(rcyl*LengthUnits, 0.0)
		              / (VelocityUnits*VelocityUnits);

		    vphi = sqrt(fmax(0.0, v_g_sqr + v_p_sqr + v_dm_sqr));
		    /*
		    if (debug && fabs(ypos) <= CellWidth[1][j] && fabs(zpos) <= CellWidth[2][k])
		    {
		      printf("r_n = %"GSYM": v_g_sqr =  %"GSYM", v_p_sqr =  %"GSYM", v_dm_sqr =  %"GSYM", vphi = %"GSYM"\n", 
			     2*r_n, v_g_sqr * VelocityUnits*VelocityUnits, v_p_sqr * VelocityUnits*VelocityUnits, 
			     v_dm_sqr * VelocityUnits*VelocityUnits, vphi * VelocityUnits);
		    }
		    */
		  }
		  else 
		  {
		    density = max(InitialDensity,
  				  DiskTable[sphere].InterpolateEquilibriumDensityTable(rcyl*LengthUnits, zpos*LengthUnits) / DensityUnits);

		    vphi = DiskTable[sphere].InterpolateEquilibriumVcircTable(rcyl*LengthUnits, zpos*LengthUnits) / VelocityUnits;
		  }

		  //vphi *= (1.0 - (1.0 - exp(SmoothFactor*r/r_max[sphere]))/(1.0 - exp(SmoothFactor)));

		  RotVelocity[0] = -vphi*sinphi*sintheta;
		  RotVelocity[1] = vphi*cosphi*sintheta;
		  RotVelocity[2] = 0;

		  Velocity[0] = RotVelocity[0] + SphereVelocity[sphere][0];
		  Velocity[1] = RotVelocity[1] + SphereVelocity[sphere][1];
		  Velocity[2] = RotVelocity[2] + SphereVelocity[sphere][2];
		}
	      } // end: loop over spheres

	/* Set density. */
	
	BaryonField[0][n] = density*BaryonMeanDensity;
	
	/* If doing multi-species (HI, etc.), set these. */
	
	if (MultiSpecies > 0)
	{
		BaryonField[HIINum][n] = CosmologySimulationInitialFractionHII *
		CoolData.HydrogenFractionByMass * BaryonField[0][n] *
			sqrt(OmegaMatterNow)/
			(OmegaMatterNow*BaryonMeanDensity*HubbleConstantNow);
		
		BaryonField[HeIINum][n] = CosmologySimulationInitialFractionHeII*
			BaryonField[0][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
		BaryonField[HeIIINum][n] = CosmologySimulationInitialFractionHeIII*
			BaryonField[0][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
		BaryonField[HeINum][n] = (1.0 - CoolData.HydrogenFractionByMass)*BaryonField[0][n] -
			BaryonField[HeIINum][n] - BaryonField[HeIIINum][n];
		
		if (MultiSpecies > 1)
		{
			BaryonField[HMNum][n] = CosmologySimulationInitialFractionHM*
			BaryonField[HIINum][n]* pow(temperature,float(0.88));
			BaryonField[H2IINum][n] = CosmologySimulationInitialFractionH2II*
				2.0*BaryonField[HIINum][n]* pow(temperature,float(1.8));
			BaryonField[H2INum][n] = CosmologySimulationInitialFractionH2I*
				BaryonField[0][n]*CoolData.HydrogenFractionByMass*pow(301.0,5.1)*
				pow(OmegaMatterNow, float(1.5))/
				(OmegaMatterNow*BaryonMeanDensity)/
			HubbleConstantNow*2.0;
		}
		
		BaryonField[HINum][n] = CoolData.HydrogenFractionByMass*BaryonField[0][n]
			- BaryonField[HIINum][n];
		if (MultiSpecies > 1)
			BaryonField[HINum][n] -= BaryonField[HMNum][n]
				+ BaryonField[H2IINum][n]
				+ BaryonField[H2INum][n];
		
		BaryonField[DeNum][n] = BaryonField[HIINum][n] + 
			0.25*BaryonField[HeIINum][n] + 0.5*BaryonField[HeIIINum][n];
		if (MultiSpecies > 1)
			BaryonField[DeNum][n] += 0.5*BaryonField[H2IINum][n] - 
		BaryonField[HMNum][n];
		
		/* Set Deuterium species (assumed to be negligible). */
		
		if (MultiSpecies > 2)
		{
			BaryonField[DINum][n] = CoolData.DeuteriumToHydrogenRatio * BaryonField[HINum][n];
			BaryonField[DIINum][n] = CoolData.DeuteriumToHydrogenRatio * BaryonField[HIINum][n];
			BaryonField[HDINum][n] = CoolData.DeuteriumToHydrogenRatio * BaryonField[H2INum][n];
		}
	}
	
	/* Set Velocities. */
	
	for (int dim = 0; dim < GridRank; dim++)
	{
		BaryonField[ivel+dim][n] = Velocity[dim];
		if (HydroMethod == MHD_RK)
			BaryonField[imag+dim][n] = MagnField[dim];
	}
	if(HydroMethod == MHD_RK)
		BaryonField[imag+GridRank][n] = 0.0;//1.0e-30;
	
	/* Set energy (thermal and then total if necessary). */
	
	BaryonField[1][n] = temperature/TemperatureUnits/ ((Gamma-1.0)*mu);
	
	if (DualEnergyFormalism)
		BaryonField[2][n] = BaryonField[1][n];
	
	if (HydroMethod != Zeus_Hydro)
		for (int dim = 0; dim < GridRank; dim++)
		{
			BaryonField[1][n] += 0.5*pow(BaryonField[ivel+dim][n], 2);
			if(HydroMethod == MHD_RK)
				BaryonField[1][n] += 0.5*pow(BaryonField[imag+dim][n], 2)/BaryonField[0][n];
		}
		

		} // end loop over grid
		
	if (UseParticles && level == 0)
	{
		printf("MHDGalaxyDisk, Processor: %"ISYM": Number of particles on current grid: %"ISYM"\n", MyProcessorNumber, NumberOfParticles);
	}
	return SUCCESS;
}

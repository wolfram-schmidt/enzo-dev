/***********************************************************************
/
/  INITIALIZE DISK GALAXIES WITH LIVE HALOS
/
/  written by: Kai Rodenbeck and Simon Selg
/  date:       July, 2020
/  modified:   August, 2023
/
/  PURPOSE:
/    Sets up one or more disk galaxies and particle halos.
/
************************************************************************/

#include <stdlib.h>
#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "EquilibriumGalaxyDisk.h"

void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);

// ============================================================================
// S. Selg (11/2019): Adjustments in order to use PRGIO. Mind that there is a
// discrimination between the states of ``SetBaryonFiels``.
// ----------------------------------------------------------------------------

int GalaxyLiveHaloInitialize(FILE *fptr, FILE *Outfptr,
			  HierarchyEntry &TopGrid, TopGridData &MetaData,
			  int SetBaryonFields
		)
{
  const char *DensName = "Density";
  const char *TEName   = "TotalEnergy";
  const char *GEName   = "GasEnergy";
  const char *Vel1Name = "x-velocity";
  const char *Vel2Name = "y-velocity";
  const char *Vel3Name = "z-velocity";

  // S. Selg (08/2019): Gravitational Potential (for output)
  const char *GPotName          = "Grav_Potential";
  const char *Bfield1Name	= "Bx";
  const char *Bfield2Name	= "By";
  const char *Bfield3Name	= "Bz";
  const char *PhiName		= "Phi";

  const char *ColourName = "colour";
  const char *ElectronName = "Electron_Density";
  const char *HIName    = "HI_Density";
  const char *HIIName   = "HII_Density";
  const char *HeIName   = "HeI_Density";
  const char *HeIIName  = "HeII_Density";
  const char *HeIIIName = "HeIII_Density";
  const char *HMName    = "HM_Density";
  const char *H2IName   = "H2I_Density";
  const char *H2IIName  = "H2II_Density";
  const char *DIName    = "DI_Density";
  const char *DIIName   = "DII_Density";
  const char *HDIName   = "HDI_Density";
  const char *MetalName = "Metal_Density";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i;

  /* set default parameters */

  int DiskNumberOfSpheres = 1;
  int DiskRefineAtStart   = TRUE;
  int UseParticles    = TRUE;
  int UseGas          = TRUE;
  int UseMetals       = FALSE;

  float InitialTemperature = 1e4;  // 10.000 K
  float InitialDensity     = 1e-3; // in code units
  float InitialMagnField   = 0.0;

  float DiskRadius[MAX_SPHERES],
	    DiskTemperature[MAX_SPHERES],
        DiskBeta[MAX_SPHERES],
        DiskMetallicity[MAX_SPHERES];

  float DiskPosition[MAX_SPHERES][MAX_DIMENSION],
	    DiskRotAxis[MAX_SPHERES][MAX_DIMENSION],
        DiskVelocity[MAX_SPHERES][MAX_DIMENSION];

  EquilibriumGalaxyDisk DiskTable[MAX_SPHERES];

  for (sphere = 0; sphere < MAX_SPHERES; sphere++) {
    DiskTemperature[sphere]	= InitialTemperature;
    DiskBeta[sphere]        = huge_number; // hydrodynamical
    DiskMetallicity[sphere]	= 0;

	/* disk at rest centered in box */

    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      DiskPosition[sphere][dim] = 0.5*(DomainLeftEdge[dim] + DomainRightEdge[dim]);
 	  DiskRotAxis[sphere][dim] = 0.0;
      DiskVelocity[sphere][dim] = 0.0;
    }
    DiskRotAxis[sphere][2] = 1.0; // rotation axis in z-direction
  }

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "DiskNumberOfSpheres = %"ISYM,
				  &DiskNumberOfSpheres);
    ret += sscanf(line, "DiskRefineAtStart = %"ISYM,
				  &DiskRefineAtStart);
    ret += sscanf(line, "UseParticles = %"ISYM,
				  &UseParticles);
    ret += sscanf(line, "UseGas = %"ISYM,
				  &UseGas);
    ret += sscanf(line, "UseMetals = %"ISYM,
				  &UseMetals);
    ret += sscanf(line, "InitialTemperature = %"FSYM,
				  &InitialTemperature);
    ret += sscanf(line, "InitialDensity = %"FSYM,
				  &InitialDensity);
    ret += sscanf(line, "InitialMagnField = %"FSYM,
				  &InitialMagnField);
    if (sscanf(line, "DiskTemperature[%"ISYM"]", &sphere) > 0)
	    ret += sscanf(line, "DiskTemperature[%"ISYM"] = %"FSYM, &sphere,
					  &DiskTemperature[sphere]);
    if (sscanf(line, "DiskBeta[%"ISYM"]", &sphere) > 0)
	    ret += sscanf(line, "DiskBeta[%"FSYM"] = %"FSYM, &sphere,
					  &DiskBeta[sphere]);
    if (sscanf(line, "DiskMetallicity[%"ISYM"]", &sphere) > 0)
      	ret += sscanf(line, "DiskMetallicity[%"ISYM"] = %"FSYM, &sphere,
					  &DiskMetallicity[sphere]);
    if (sscanf(line, "DiskRadius[%"ISYM"]", &sphere) > 0)
      	ret += sscanf(line, "DiskRadius[%"ISYM"] = %"FSYM, &sphere,
					  &DiskRadius[sphere]);
    if (sscanf(line, "DiskPosition[%"ISYM"]", &sphere) > 0)
      	ret += sscanf(line, "DiskPosition[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM,
					  &sphere,
					  &DiskPosition[sphere][0],
					  &DiskPosition[sphere][1],
					  &DiskPosition[sphere][2]);
    if (sscanf(line, "DiskRotAxis[%"ISYM"]", &sphere) > 0)
      	ret += sscanf(line, "DiskRotAxis[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM,
					  &sphere,
					  &DiskRotAxis[sphere][0],
					  &DiskRotAxis[sphere][1],
					  &DiskRotAxis[sphere][2]);
    if (sscanf(line, "DiskVelocity[%"ISYM"]", &sphere) > 0)
      	ret += sscanf(line, "DiskVelocity[%"ISYM"] = %"FSYM" %"FSYM" %"FSYM,
					  &sphere,
					  &DiskVelocity[sphere][0],
					  &DiskVelocity[sphere][1],
					  &DiskVelocity[sphere][2]);

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "Disk")
	    && line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file

  if (debug && MyProcessorNumber == ROOT_PROCESSOR)
	printf("Using metals: %"ISYM"\n", UseMetals);

  for (sphere = 0; sphere < DiskNumberOfSpheres; sphere++)
    DiskTable[sphere].ReadInData("disk-hernq.dat");

  /* set up grid */

  /* =========================================================================
   * S. Selg (11/2019): Implementation of Parallel Root Grid IO (PRGIO)
   * =========================================================================
   */
  // A bit of pre-computing
  // TOP GRID SPACING
  float TopGridSpacing = float(MetaData.TopGridDims[0]);
  HierarchyEntry *CurrentGrid;
  CurrentGrid = &TopGrid;
  while (CurrentGrid != NULL) {
  	if (CurrentGrid->GridData->GalaxyLiveHaloInitializeGrid(
					DiskNumberOfSpheres,
					DiskTable,
					DiskPosition,
					DiskRotAxis,
					DiskVelocity,
					DiskRadius,
					DiskTemperature,
				    DiskBeta,
					DiskMetallicity,
					InitialTemperature,
					InitialDensity,
					InitialMagnField,
					UseParticles,
					UseGas,
					UseMetals,
					0,
					SetBaryonFields,
					1) == FAIL) {
	  ENZO_FAIL("Error in GalaxyLiveHaloInitializeGrid.");
	}
	CurrentGrid = CurrentGrid->NextGridThisLevel;
  }

  /* S. Selg (11/2019): This will be done after first initialization! (if you are
   * using prgio. */

  if (SetBaryonFields) {

	/* Convert minimum initial overdensity for refinement to mass
       (unless MinimumMass itself was actually set). */

	for (int count = 0; count < MAX_FLAGGING_METHODS; count++)
	  if (MinimumMassForRefinement[count] == FLOAT_UNDEFINED) {
		MinimumMassForRefinement[count] = MinimumOverDensityForRefinement[count];
		for (dim = 0; dim < MetaData.TopGridRank; dim++)
		  MinimumMassForRefinement[count] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
			float(MetaData.TopGridDims[dim]);
	  }

	/* (S. Selg, 10/2019) The following lines use the code from ClusterInitialize and are to
	   implement a refinement at start. If requested, refine the grids to the desired level. */

	if (DiskRefineAtStart) {
	  // Declare, initialize and fill out the LevelArray
	  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
	  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
		  LevelArray[level] = NULL;
	  AddLevel(LevelArray, &TopGrid, 0);
	  // Add levels to the maximum depth until no new levels are created,
	  // and re-initialize the level after it's created.
	  for (level = 0; level < MaximumRefinementLevel; level++)
	  {
		  if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL)
		  {
			  fprintf(stderr, "Error in RebuildHierarchy.\n");
			  return FAIL;
		  }
		  if (LevelArray[level+1] == NULL)
			  break;
		  LevelHierarchyEntry *Temp = LevelArray[level+1];
		  while (Temp != NULL)
		  {
			  if (Temp->GridData->GalaxyLiveHaloInitializeGrid(
				DiskNumberOfSpheres,
				DiskTable,
				DiskPosition,
				DiskRotAxis,
				DiskVelocity,
				DiskRadius,
				DiskTemperature,
				DiskBeta,
				DiskMetallicity,
				InitialTemperature,
				InitialDensity,
				InitialMagnField,
				UseParticles,
				UseGas,
				UseMetals,
				level,   // S. Selg (11/2019, used to be level+1)
				SetBaryonFields,
				0) == FAIL)
				{
					fprintf(stderr, "Error in GalaxyLiveHaloInitializeGrid.\n");
					return FAIL;
				}
			  Temp = Temp->NextGridThisLevel;
		  }
	  } // end: loop over levels

	  // Loop back from the bottom, restoring the consistency among levels.
	  for (level = MaximumRefinementLevel; level > 0; level --)
	  {
		  LevelHierarchyEntry *Temp = LevelArray[level];
		  while (Temp != NULL)
		  {
			  if (Temp->GridData->ProjectSolutionToParentGrid(
				*LevelArray[level-1]->GridData) == FAIL)
			  {
				  fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
				  return FAIL;
			  }
			  Temp = Temp->NextGridThisLevel;
		  }
	  }
    } // end: if (DiskRefineAtStart)

	/* set up field names and units */

	int count = 0;
	DataLabel[count++] = (char*) DensName;
	DataLabel[count++] = (char*) TEName;

	if (DualEnergyFormalism)
		DataLabel[count++] = (char*) GEName;

	DataLabel[count++] = (char*) Vel1Name;
	DataLabel[count++] = (char*) Vel2Name;
	DataLabel[count++] = (char*) Vel3Name;

	if(HydroMethod == MHD_RK)
	{
		DataLabel[count++] = (char*) Bfield1Name;
		DataLabel[count++] = (char*) Bfield2Name;
		DataLabel[count++] = (char*) Bfield3Name;
		DataLabel[count++] = (char*) PhiName;
	}
	if (MultiSpecies)
	{
		DataLabel[count++] = (char*) ElectronName;
		DataLabel[count++] = (char*) HIName;
		DataLabel[count++] = (char*) HIIName;
		DataLabel[count++] = (char*) HeIName;
		DataLabel[count++] = (char*) HeIIName;
		DataLabel[count++] = (char*) HeIIIName;
		if (MultiSpecies > 1)
		{
			DataLabel[count++] = (char*) HMName;
			DataLabel[count++] = (char*) H2IName;
			DataLabel[count++] = (char*) H2IIName;
		}
		if (MultiSpecies > 2)
		{
			DataLabel[count++] = (char*) DIName;
			DataLabel[count++] = (char*) DIIName;
			DataLabel[count++] = (char*) HDIName;
		}
	}
	if (UseMetals)
	  DataLabel[count++] = (char*) MetalName;

	// S. Selg (08/2019): toggle output of gravitational potential
	if (WritePotential)
	  DataLabel[count++] = (char*) GPotName;

	for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

	/* Write parameters to parameter output file */

	if (MyProcessorNumber == ROOT_PROCESSOR) {
	  fprintf(Outfptr, "DiskNumberOfSpheres    = %"ISYM"\n",
			  DiskNumberOfSpheres);
	  fprintf(Outfptr, "DiskRefineAtStart      = %"ISYM"\n",
			  DiskRefineAtStart);
	  fprintf(Outfptr, "UseParticles       = %"ISYM"\n",
			  UseParticles);
	  fprintf(Outfptr, "UseGas             = %"ISYM"\n",
			  UseGas);
	  fprintf(Outfptr, "UseMetals          = %"ISYM"\n",
			  UseMetals);
	  fprintf(Outfptr, "InitialTemperature = %"FSYM"\n",
			  InitialTemperature);
	  fprintf(Outfptr, "InitialDensity     = %"FSYM"\n",
			  InitialDensity);
	  fprintf(Outfptr, "InitialMagnField     = %"FSYM"\n",
			  InitialMagnField);
	  for (sphere = 0; sphere < DiskNumberOfSpheres; sphere++) {
		fprintf(Outfptr, "DiskRadius[%"ISYM"] = %"FSYM"\n", sphere,
				DiskRadius[sphere]);
		fprintf(Outfptr, "DiskTemperature[%"ISYM"] = %"FSYM"\n", sphere,
				DiskTemperature[sphere]);
		fprintf(Outfptr, "DiskBeta[%"ISYM"] = %"FSYM"\n", sphere,
				DiskBeta[sphere]);
		fprintf(Outfptr, "DiskMetallicity[%"ISYM"] = %"FSYM"\n", sphere,
				DiskMetallicity[sphere]);
		fprintf(Outfptr, "DiskPosition[%"ISYM"] = ", sphere);
		WriteListOfFloats(Outfptr, MetaData.TopGridRank,
						  DiskPosition[sphere]);
		fprintf(Outfptr, "DiskRotAxis[%"ISYM"] = ", sphere);
		WriteListOfFloats(Outfptr, MetaData.TopGridRank,
						  DiskRotAxis[sphere]);
		fprintf(Outfptr, "DiskVelocity[%"ISYM"] = ", sphere);
		WriteListOfFloats(Outfptr, MetaData.TopGridRank,
						  DiskVelocity[sphere]);
	  }
	}
  } // end: if SetBaryonFields

  return SUCCESS;

}

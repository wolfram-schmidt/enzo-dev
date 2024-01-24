/***********************************************************************
/
/  INITIALIZE a moving spherical subcluster simulation
/
/  written by: Julian Adamek
/  date:       June 2006
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine intializes a new simulation based on the parameter file.
//
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <string.h>
#include <stdio.h>
#include <math.h>

#include "ErrorExceptions.h"
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
 
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);

int CloudWindInitialize(FILE *fptr, FILE *Outfptr, 
			HierarchyEntry &TopGrid, TopGridData &MetaData, ExternalBoundary &Exterior)
{
  char *DensName = "Density";
  char *TEName = "Total_Energy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *GEName     = "Internal_Energy";
  char *MetalName   = "Metal_Density";
  char *ColourName = "colour";
  char *AveMomt1Name = "AveMomtX";
  char *AveMomt2Name = "AveMomtY";
  char *AveMomt3Name = "AveMomtZ";

  /* declarations */

  char line[MAX_LINE_LENGTH];
  int  dim, ret, level,NumberOfSubgridZones[MAX_DIMENSION], SubgridDims[MAX_DIMENSION];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION]; // CF (float)
  int ii;

  /* set default parameters */
  /* the meaning of them is explained in the parameter file */
  
  FLOAT CloudWindCutoffRadius = 0.1;      // used, in code units, for the Unbound sub-setup
  float CloudWindBeta = 0.8;
  float CloudWindCentralDensity = 2.0;
  float CloudWindExternalDensity = 1.0;
  float CloudWindVelocity = 0.0;
  float CloudWindExternalTotalEnergy = 2.0;
  float CloudWindCentralTotalEnergy = 5.0;
  float CloudWindHSETolerance = 1e-1;
  float CloudWindMetallicity = 1.0e-10;
  int   CloudWindUseMetallicityField = FALSE;
  int   CloudWindUnbound = FALSE;

  int   HereGridRank = MetaData.TopGridRank;
  FLOAT CloudWindSubgridLeft = 0.0;        // start of subgrid
  FLOAT CloudWindSubgridRight= 0.0;        // end of subgrid  
  int   CloudWindRefineAtStart = TRUE;

  /* make sure it is not meaningless */
  if (MetaData.TopGridRank == 1) {
    ENZO_VFAIL("Cannot do CloudWind in 1D.\n")
  }

  /* read input from file */

  if (debug) printf("CloudWindInitialize: reading problem-specific parameters.\n");

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "CloudWindVelocity = %"FSYM, &CloudWindVelocity);
    ret += sscanf(line, "CloudWindCutoffRadius = %"PSYM, &CloudWindCutoffRadius);
    ret += sscanf(line, "CloudWindCentralDensity = %"FSYM, &CloudWindCentralDensity);      
    ret += sscanf(line, "CloudWindExternalDensity = %"FSYM, &CloudWindExternalDensity);
    ret += sscanf(line, "CloudWindExternalTotalEnergy = %"FSYM, &CloudWindExternalTotalEnergy);
    ret += sscanf(line, "CloudWindCentralTotalEnergy = %"FSYM, &CloudWindCentralTotalEnergy);
    ret += sscanf(line, "CloudWindBeta = %"FSYM, &CloudWindBeta);
    ret += sscanf(line, "CloudWindHSETolerance = %"FSYM, &CloudWindHSETolerance);
    ret += sscanf(line, "CloudWindMetallicity = %"FSYM, &CloudWindMetallicity);
    ret += sscanf(line, "CloudWindUseMetallicityField = %"ISYM, &CloudWindUseMetallicityField);
    ret += sscanf(line, "CloudWindUnbound = %"ISYM, &CloudWindUnbound);
    ret += sscanf(line, "CloudWindSubgridLeft = %"PSYM, &CloudWindSubgridLeft);
    ret += sscanf(line, "CloudWindSubgridRight = %"PSYM, &CloudWindSubgridRight);
    ret += sscanf(line, "CloudWindRefineAtStart = %"ISYM, &CloudWindRefineAtStart);

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "CloudWind") &&
	line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);
  }

  if (debug && MyProcessorNumber == ROOT_PROCESSOR) {
    printf("CloudWindVelocity = %f\n", CloudWindVelocity);
    printf("CloudWindCutoffRadius = %f\n", CloudWindCutoffRadius);
    printf("CloudWindCentralDensity = %f\n", CloudWindCentralDensity);
    printf("CloudWindExternalDensity = %f\n", CloudWindExternalDensity);
    printf("CloudWindExternalTotalEnergy = %f\n", CloudWindExternalTotalEnergy);
    printf("CloudWindCentralTotalEnergy = %f\n", CloudWindCentralTotalEnergy);
    printf("CloudWindBeta = %f\n", CloudWindBeta);
    printf("CloudWindMetallicity = %f\n", CloudWindMetallicity);
    printf("CloudWindUseMetallicityField  = %d\n"  , CloudWindUseMetallicityField);
    printf("CloudWindHSETolerance = %f\n", CloudWindHSETolerance); 
    printf("CloudWindUnbound  = %d\n"  , CloudWindUnbound);
    printf("CloudWindSubgridLeft  = %f\n"  , CloudWindSubgridLeft);
    printf("CloudWindSubgridRight = %f\n" , CloudWindSubgridRight);
    printf("CloudWindRefineAtStart  = %d\n\n"  , CloudWindRefineAtStart);
  }

  /* set up grid */

  if (TopGrid.GridData->CloudWindInitializeGrid(CloudWindVelocity,
						CloudWindCutoffRadius,
						CloudWindCentralDensity,
						CloudWindExternalDensity,
						CloudWindExternalTotalEnergy,
						CloudWindCentralTotalEnergy,
						CloudWindBeta,
						CloudWindMetallicity,
						CloudWindUseMetallicityField,
						CloudWindHSETolerance,
						CloudWindUnbound) == FAIL) {
    ENZO_FAIL("Error in CloudWindInitializeGrid.\n")
  }

  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */

  if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
    MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
    for (int dim = 0; dim < MetaData.TopGridRank; dim++)
      MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
	float(MetaData.TopGridDims[dim]);
  }

  if (CloudWindRefineAtStart) {

    /* Declare, initialize and fill out the LevelArray. */

    LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      LevelArray[level] = NULL;
    AddLevel(LevelArray, &TopGrid, 0);

    /* Add levels to the maximum depth or until no new levels are created,
       and re-initialize the level after it is created. */

    for (level = 0; level < MaximumRefinementLevel; level++) {
      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
        fprintf(stderr, "Error in RebuildHierarchy.\n");
        return FAIL;
      }
      if (LevelArray[level+1] == NULL)
        break;
      LevelHierarchyEntry *Temp = LevelArray[level+1];
      while (Temp != NULL) {

	if (Temp->GridData->CloudWindInitializeGrid(CloudWindVelocity,
						    CloudWindCutoffRadius,
						    CloudWindCentralDensity,
						    CloudWindExternalDensity,
						    CloudWindExternalTotalEnergy,
						    CloudWindCentralTotalEnergy,
						    CloudWindBeta,
						    CloudWindMetallicity,
						    CloudWindUseMetallicityField,
						    CloudWindHSETolerance,
						    CloudWindUnbound
						    ) == FAIL) {
	  ENZO_FAIL("Error in CloudWindInitializeGrid.\n")
	}

	Temp = Temp->NextGridThisLevel;
      }
    } // end: loop over levels

    /* Loop back from the bottom, restoring the consistency amoung levels. */

    for (level = MaximumRefinementLevel; level > 0; level--) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
        if (Temp->GridData->ProjectSolutionToParentGrid(
							*LevelArray[level-1]->GridData) == FAIL) {
          fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
          return FAIL;
        }
        Temp = Temp->NextGridThisLevel;
      }
    }

  } // end: if (CloudWindRefineAtStart)


  /* here for subgrids and other unuseful stuff to test */
  //

  /* If requested, create a subgrid */

  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    NumberOfSubgridZones[dim] =
      nint((CloudWindSubgridRight - CloudWindSubgridLeft)/
           ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
            float(MetaData.TopGridDims[dim])))
      *RefineBy;

  if (NumberOfSubgridZones[0] > 0) {

    /* Error check */

    if (CloudWindRefineAtStart) {
      ENZO_VFAIL("Cannot RefineAtStart AND create subgrid.\n");
    }

    /* create a new HierarchyEntry, attach to the top grid and fill it out */

    HierarchyEntry *Subgrid    = new HierarchyEntry;
    TopGrid.NextGridNextLevel  = Subgrid;
    Subgrid->NextGridNextLevel = NULL;
    Subgrid->NextGridThisLevel = NULL;
    Subgrid->ParentGrid        = &TopGrid;

    /* compute the dimensions and left/right edges for the subgrid */
 
    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      SubgridDims[dim] = NumberOfSubgridZones[dim] + 2*NumberOfGhostZones;
      LeftEdge[dim]    = CloudWindSubgridLeft;
      RightEdge[dim]   = CloudWindSubgridRight;
    }

    /* create a new subgrid and initialize it */

    Subgrid->GridData = new grid;
    Subgrid->GridData->InheritProperties(TopGrid.GridData);
    Subgrid->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
                                   LeftEdge, RightEdge, 0);
    if (Subgrid->GridData->CloudWindInitializeGrid(CloudWindVelocity,
						   CloudWindCutoffRadius,
						   CloudWindCentralDensity,
						   CloudWindExternalDensity,
						   CloudWindExternalTotalEnergy,
						   CloudWindCentralTotalEnergy,
						   CloudWindBeta,
						   CloudWindMetallicity,
						   CloudWindUseMetallicityField,
						   CloudWindHSETolerance,
						   CloudWindUnbound) == FAIL) {
      ENZO_FAIL("Error in CloudWindInitializeGrid.\n");
    }
  }


  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = DensName;
  DataLabel[count++] = TEName;
  if (DualEnergyFormalism)
      DataLabel[count++] = GEName;
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  if (HereGridRank != 2)
      DataLabel[count++] = Vel3Name;
  if (CloudWindUseMetallicityField)
      DataLabel[count++] = MetalName;
  /* TODO  
  if (SGSModel) {
    if (ShearImproved) {
      DataLabel[count++] = AveMomt1Name;
      DataLabel[count++] = AveMomt2Name;
      DataLabel[count++] = AveMomt3Name;
    }
  } 
  */

  /* set up boundary conditions */

  Exterior.Prepare(TopGrid.GridData);

  float InflowValue[MAX_NUMBER_OF_BARYON_FIELDS], Dummy[MAX_NUMBER_OF_BARYON_FIELDS];
  int index = 0;
  InflowValue[index++] = CloudWindExternalDensity;
  InflowValue[index++] = CloudWindExternalTotalEnergy + 0.5 * CloudWindVelocity * CloudWindVelocity;
  if (DualEnergyFormalism) 
      InflowValue[index++] = CloudWindExternalTotalEnergy;
  InflowValue[index++] = CloudWindVelocity;
  InflowValue[index++] = 0.0;
  if (HereGridRank != 2)
      InflowValue[index++] = 0.0;
  if (CloudWindUseMetallicityField)
      InflowValue[index++] = 1.0e-12;
  /* TODO
  if (SGSModel) {
    if (ShearImproved) {
      InflowValue[index++] = CloudWindExternalDensity*CloudWindVelocity;
      InflowValue[index++] = 0.0;
      InflowValue[index++] = 0.0;
    }
  } 
  */

  int fieldnumber = index;
  
  for(ii = fieldnumber+1; ii < MAX_NUMBER_OF_BARYON_FIELDS; ii++) {
      InflowValue[ii] = 0.0;    
  }       
  
  if (Exterior.InitializeExternalBoundaryFace(0, inflow, outflow, InflowValue, Dummy) == FAIL) {
    fprintf(stderr, "Error in InitializeExternalBoundaryFace.\n");
    return FAIL;
  }

  if (MetaData.TopGridRank > 1)
    Exterior.InitializeExternalBoundaryFace(1, outflow, outflow, Dummy, Dummy);

  if (MetaData.TopGridRank > 2)
    Exterior.InitializeExternalBoundaryFace(2, outflow, outflow, Dummy, Dummy);

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "CloudWindVelocity = %f\n", CloudWindVelocity);
    fprintf(Outfptr, "CloudWindCutoffRadius = %f\n", CloudWindCutoffRadius);
    fprintf(Outfptr, "CloudWindCentralDensity = %f\n", CloudWindCentralDensity);
    fprintf(Outfptr, "CloudWindExternalDensity = %f\n", CloudWindExternalDensity);
    fprintf(Outfptr, "CloudWindExternalTotalEnergy = %f\n", CloudWindExternalTotalEnergy);
    fprintf(Outfptr, "CloudWindCentralTotalEnergy = %f\n", CloudWindCentralTotalEnergy);
    fprintf(Outfptr, "CloudWindBeta = %f\n", CloudWindBeta);
    fprintf(Outfptr, "CloudWindMetallicity = %f\n", CloudWindMetallicity);
    fprintf(Outfptr, "CloudWindUseMetallicityField = %f\n", CloudWindUseMetallicityField);
    fprintf(Outfptr, "CloudWindHSETolerance = %f\n", CloudWindHSETolerance); 
    fprintf(Outfptr, "CloudWindUnbound  = %f\n"  , CloudWindUnbound);
    fprintf(Outfptr, "CloudWindSubgridLeft  = %f\n"  , CloudWindSubgridLeft);
    fprintf(Outfptr, "CloudWindSubgridRight = %f\n" , CloudWindSubgridRight);
    fprintf(Outfptr, "CloudWindRefineAtStart  = %f\n\n"  , CloudWindRefineAtStart);
  }

#ifdef USE_MPI

  // BWO: this forces the synchronization of the various point source gravity
  // parameters between processors.  If this is not done, things go to pieces!

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Bcast(&PointSourceGravityConstant,1,DataType,ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(&PointSourceGravityCoreRadius,1,DataType,ROOT_PROCESSOR, MPI_COMM_WORLD);

#endif

  return SUCCESS;
}

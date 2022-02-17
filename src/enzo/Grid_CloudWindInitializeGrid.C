/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID TO A moving subcluster setup)
/
/  written by:  Julian Adamek 
/  date:        June 2006
/  modified1:   Luigi Iapichino: included a procedure for relaxing the initial 
/               model to HydroStatic Equilibrium (HSE) from Zingale et al.,
/               ApJS 143 (2002), 539. A King profile for gravity, 
/               with rho_c(DM) = 10 * rho_c(gas), is here assumed. 
/               September 2006.
/
/  modified2:   Luigi Iapichino: included support for a flat gas cloud,
/               without structure and gravity. Densities, radius, and external
/               thermal energy given as parameters, thermal energy inside is set by pressure 
/               equilibrium. September 2013.
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include <iostream>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
using namespace std;

int grid::CloudWindInitializeGrid(float CloudWindVelocity,
					 FLOAT CloudWindCutoffRadius,
					 float CloudWindCentralDensity,
					 float CloudWindExternalDensity,
					 float CloudWindExternalTotalEnergy,
					 float CloudWindCentralTotalEnergy,
					 float CloudWindBeta,
					 float CloudWindHSETolerance,
                                           int CloudWindUnbound)
{

  /* create fields */

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  int DualNum = NumberOfBaryonFields;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  int VelNum = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;
       
  /* TODO
  int SGSENum = NumberOfBaryonFields;
  if (SGSModel)                    
      FieldType[NumberOfBaryonFields++] = SGSTurbEnergyDensity;  

  int AveMomt1Num, AveMomt2Num, AveMomt3Num;	
  if (SGSModel && ShearImproved)
  {
      AveMomt1Num = NumberOfBaryonFields;
      FieldType[NumberOfBaryonFields++] = AveMomtDensity1;
      
      AveMomt2Num = NumberOfBaryonFields;
      FieldType[NumberOfBaryonFields++] = AveMomtDensity2;
      
      AveMomt3Num = NumberOfBaryonFields;
      FieldType[NumberOfBaryonFields++] = AveMomtDensity3;
  }
  */

  if (debug) cout << "NumberOfBaryonFields = " << NumberOfBaryonFields << endl;

  /* Return if this doesn't concern us. */

    if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* declarations */

  int size, dim, field, i, j, k, n;
  FLOAT r;

  float MappingUnit;
  int HSECellNumber;  
  int ii, iter;
  float *peos, phse, *dens;
  FLOAT *radius;
  float *grav, dens_zone, pres_zone, dens_map;
  float accel, drho;
  FLOAT rad, rcubed, rcore;
  bool converge;
  FILE *original_model, *hse_model;
  const int max_iter = 50000; 

/*   HSE part begins. It is executed, unless we initialize a 
     cloud without internal structure and gravity             */         

  if (CloudWindUnbound == 0) {

  HSECellNumber = 0;
  for (dim = 0; dim < GridRank; dim++) {
    if (HSECellNumber < GridDimension[dim]) {
      HSECellNumber = GridDimension[dim];
    }      
  }           
  if (HSECellNumber < 100) {
    HSECellNumber = 100;
  }

  if (MyProcessorNumber == ROOT_PROCESSOR && debug) 
      cout << "HSECellNumber = " << HSECellNumber << endl;

  peos   = new float[HSECellNumber +1];
  dens   = new float[HSECellNumber +1];
  grav   = new float[HSECellNumber +1];
  radius = new FLOAT[HSECellNumber +1];   

  /* the next is the size of a HSE cell in code units. 
     P...CoreRadius * 2 is a guess of the subcluster size, increase if necessary */ 
  MappingUnit = PointSourceGravityCoreRadius * 2.0 / ( (float) HSECellNumber);

  /* now fill the arrays with initial values. We have to assume that the filled values are
     cell averaged quantities: it is not true, but just improve it if you don't like ;-) */

  dens[0] = CloudWindCentralDensity;
  /* the next one is the pressure, as calculated from the EoS */
  peos[0] = CloudWindCentralDensity * CloudWindCentralTotalEnergy * (Gamma -1);
  rad = 0.0;

  //fprintf(stderr, "dens,pres %f %f\n", dens[0], peos[0]);

  /* Write the file with the initial model */ 

  if (MyProcessorNumber == ROOT_PROCESSOR) {
      if (debug) cout << "Writing data to original_model.dat" << endl;

      original_model = fopen("original_model.dat", "w");
      if (original_model == NULL)
	  fprintf(stderr, "Cannot open %s\n", "original_model.dat");
      
      fprintf(original_model, "Initial model\n");
      fprintf(original_model, "Radius\t\t Density\t Pressure\t\n" );
      fprintf(original_model, "%8f\t %8f\t %8f\t\n", rad, dens[0], peos[0]);
  } 

  for (ii = 1; ii <= HSECellNumber; ii++) {
      
     rad += MappingUnit;     
     dens[ii] = CloudWindCentralDensity * pow(1.0 + rad*rad/PointSourceGravityCoreRadius/
						    PointSourceGravityCoreRadius, -3.0 * CloudWindBeta/2.0);
     peos[ii] = dens[ii] * CloudWindCentralTotalEnergy * (Gamma -1);

     if (MyProcessorNumber == ROOT_PROCESSOR) 
	fprintf(original_model, "%8f\t %8f\t %8f\t\n", rad, dens[ii], peos[ii]); 

     /* Find the gravitational acceleration. Sorry for the loss of generality, but 
	a King profile is here assumed. */
     rcubed = rad * rad * rad;
     rcore = PointSourceGravityCoreRadius;
     FLOAT x = rad / rcore;
     rcubed /= (rcore*rcore*rcore * (-x / sqrt(1 + x*x)
				     + log(x + sqrt(1 + x*x))));
     
     /* the DM central density is assumed to be 10 * gas-density */             
     accel = PointSourceGravityConstant * (CloudWindCentralDensity * 10.0) / rcubed;
     
     grav[ii] = - accel * rad;

     /* iterate the HSE procedure until the tolerance is reached, but for a finite time */

     dens_zone = dens[ii];
     pres_zone = peos[ii];

     /* HSE iteration starts */

     for (iter = 0; iter < max_iter; iter++) {      
       /* Find the pressure that satisfies HSE at the interface between ii and ii-1 */
       converge = false;
       phse = peos[ii-1] + 0.5 * MappingUnit * grav[ii] * (dens_zone + dens[ii-1]);
       /* Now we have a pressure difference between  phse and pres_zone: 
	  compute the related density difference                         */ 
       drho = (phse - pres_zone) / (pres_zone/dens_zone - 0.5 * MappingUnit * grav[ii]);
       dens_zone += drho;
      
       /* Check convergence or exit */

       if (fabs(drho) < CloudWindHSETolerance * dens_zone){
	 iter = max_iter;
	 converge = true;     
       }
       
     }
     
     if (converge==false) {
       fprintf(stderr, "The HSE iteration does not converge in zone %d\n", ii);
       if (MyProcessorNumber == ROOT_PROCESSOR) fclose(original_model);
       return FAIL;
     }
     
     dens[ii] = dens_zone;
     peos[ii] = dens[ii] * CloudWindCentralTotalEnergy * (Gamma -1);         
  }
  
  if (MyProcessorNumber == ROOT_PROCESSOR) fclose(original_model);

  /* Print the results of the HSE procedure on a file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
      if (debug) cout << "Writing data to hse_model.dat" << endl;

      hse_model = fopen("hse_model.dat", "w");
      if (hse_model == NULL)
	  fprintf(stderr, "Cannot open %s\n", "hse_model.dat");
      
      fprintf(hse_model, "Model in HSE\n");
      fprintf(hse_model, "Radius\t\t Density\t Pressure\t\n" ); 
  }

  for (ii = 0; ii <= HSECellNumber; ii++) {
     
    rad = (FLOAT) ii * MappingUnit;
    radius[ii] = rad;
    if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(hse_model, "%8f\t %8f\t %8f\t\n", rad, dens[ii], peos[ii]);
    
  }

  if (MyProcessorNumber == ROOT_PROCESSOR) fclose(hse_model);

 } // if CloudWindUnbound == 0
/* HSE part ends

  /* Here is the part for filling the cells */

  if (debug) cout << "Allocating fields..." << endl;

  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  if (debug) cout << "Size = " << size << endl;

  for (field = 0; field < NumberOfBaryonFields; field++)
    BaryonField[field] = new float[size];

  if (debug) cout << "Filling cells..." << endl;  
  
  /* set density, total energy and velocity in problem dimension */
  
  for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++)
      for (i = 0; i < GridDimension[0]; i++) {
	
	/* calculate radial coordinate of cell */
	
	FLOAT xpos = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - PointSourceGravityPosition[0];
	FLOAT ypos = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - PointSourceGravityPosition[1];
	if (GridRank != 2){
	  FLOAT zpos = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - PointSourceGravityPosition[2];
	  r = sqrt((FLOAT) xpos*xpos + ypos*ypos + zpos*zpos);
	}
	else {
	  r = sqrt((FLOAT) xpos*xpos + ypos*ypos);
	}

        if (CloudWindUnbound == 0) {

	  /* Make a linear interpolation to r from the HSE model.
	     Be careful: the cell can be outside the cluster, but your initial model must arrive at central distances
	     well outside the cluster, too. In this way, if you write nonsense in dens_map, the next condition on 
	     pressures puts you on the safe side  */
	
	  ii = 0;
	  while(radius[ii] < r && ii <= HSECellNumber) {
	    ii++; 
	  }
	
	  dens_map =(ii > HSECellNumber? 0.01 : dens[ii-1] * (radius[ii] - r) / (radius[ii] - radius[ii-1]) + 
		     dens[ii] * (r - radius[ii-1]) / (radius[ii] - radius[ii-1]));
	
	  /* Fill the cells, distinguishing between cluster and background 
	     by comparing the cell pressure with the external pressure     */
	
	  if (dens_map > (CloudWindExternalTotalEnergy * CloudWindExternalDensity / 
			  CloudWindCentralTotalEnergy)) {
	    // the interface between the subcluster and the external medium is set by pressure equilibrium
	    BaryonField[0][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = dens_map;
	    BaryonField[1][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 
	      CloudWindCentralTotalEnergy;
	    BaryonField[VelNum][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 0.0;
	    BaryonField[VelNum+1][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 0.0;
	    if (GridRank != 2)
	      BaryonField[VelNum+2][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 0.0;
	    if (DualEnergyFormalism){
	      BaryonField[DualNum][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 
		CloudWindCentralTotalEnergy;
	    }	   
	  }
	  else {
	    BaryonField[0][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 
	      CloudWindExternalDensity;
	    
	    BaryonField[1][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 
	      CloudWindExternalTotalEnergy;
	  
	    BaryonField[1][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] += 
	      0.5 * CloudWindVelocity * CloudWindVelocity;

	    BaryonField[VelNum][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 
	      CloudWindVelocity;

	    BaryonField[VelNum+1][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 0.0;
	    if (GridRank != 2)
	      BaryonField[VelNum+2][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 0.0;
	    if (DualEnergyFormalism){
	      BaryonField[DualNum][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 
		CloudWindExternalTotalEnergy;
	    }
	  }
	} // end CloudWindUnbound == 0 
	else { // CloudWindUnbound == 1

	  /* Fill the cells, distinguishing between cloud and background */

	  if (r < CloudWindCutoffRadius) {    //cutoff radius
	    
	    BaryonField[0][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 
	      CloudWindCentralDensity;
	    BaryonField[1][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 
	      CloudWindExternalDensity * CloudWindExternalTotalEnergy / CloudWindCentralDensity; 
	    /* the previous relation sets the thermal energy in the cloud by pressure equilibrium */
	    BaryonField[VelNum][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 0.0;
	    BaryonField[VelNum+1][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 0.0;
	    if (GridRank != 2)
	      BaryonField[VelNum+2][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 0.0;
	    if (DualEnergyFormalism){
	      BaryonField[DualNum][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 
		CloudWindExternalDensity * CloudWindExternalTotalEnergy / CloudWindCentralDensity;
	    }	   
	  }
	  else {
	    BaryonField[0][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 
	      CloudWindExternalDensity;
	    
	    BaryonField[1][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 
	      CloudWindExternalTotalEnergy;
	    
	    BaryonField[1][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] += 
	      0.5 * CloudWindVelocity * CloudWindVelocity;
	    
	    BaryonField[VelNum][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 
	      CloudWindVelocity;
	    
	    BaryonField[VelNum+1][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 0.0;
	    if (GridRank != 2)
	      BaryonField[VelNum+2][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 0.0;
	    if (DualEnergyFormalism){
	      BaryonField[DualNum][k*GridDimension[0]*GridDimension[1] + j*GridDimension[0] + i] = 
		CloudWindExternalTotalEnergy;
	    }
	  }
	}  
      }

  /* set SGS fields */
  
  /* TODO
  if (SGSModel) { 
    if (debug) cout << "Initializing field " << SGSENum << endl;
    for (int i = 0; i < size; i++)
          BaryonField[SGSENum][i] = InitialSGSTurbEnergy;
    if (ShearImproved) {
      if (debug) cout << "Initializing fields " << AveMomt1Num << " " << AveMomt2Num << " " << AveMomt3Num << endl;
      for (int i = 0; i < size; i++) {
	BaryonField[AveMomt1Num][i] = BaryonField[0][i]*BaryonField[VelNum][i];
	BaryonField[AveMomt2Num][i] = BaryonField[0][i]*BaryonField[VelNum+1][i];
	BaryonField[AveMomt3Num][i] = BaryonField[0][i]*BaryonField[VelNum+2][i];
      }
    }
  }
  */

  if (CloudWindUnbound == 0) {
   delete  peos;
   delete  dens;
   delete  grav;
   delete  radius;  
  }

  return SUCCESS;
}

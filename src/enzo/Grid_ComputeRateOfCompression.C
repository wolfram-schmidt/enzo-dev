/***********************************************************************
/
/  GRID CLASS (COMPUTE RATE OF COMPRESSION)
/
/  written by: Wolfram Schmidt
/  date:       October 2011  
/  modified:   February 2022
/
/  PURPOSE: computes rate of compression for refinement by variability
/           also used for output if VelAnyl = 1
/
/  RETURNS: error status
/
************************************************************************/

#include <stdlib.h>
#include "preincludes.h"
#include <iostream>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "fortran.def"
#include "Grid.h"

using namespace std;

int FindField(int f, int farray[], int n);

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);


/**
 * Computes rate of compression = -D d/D t =
 *
 *     div((grad P)/rho) + 0.5(|S|^2 - omega^2) - div f + 4 pi G rho
 *
 * @return returns SUCCESS or FAIL
 */
int grid::ComputeRateOfCompression(float* RateOfCompression)
{
	if (ProcessorNumber != MyProcessorNumber)
	{
		return SUCCESS;
	}

	if(debug) cout << "[" << MyProcessorNumber << "]" << " Computing rate of compression" << endl;

	int DensNum;
	if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields)) < 0)
	{
		fprintf(stderr, "Error in grid::ComputeRateOfCompression: Cannot find density.\n");
		return FAIL;
	}

	int size = 1;
	for (int dim = 0; dim < GridRank; dim++)
	    size *= GridDimension[dim];

	int rank_sqr=GridRank*GridRank;

        for (int m=0;m<rank_sqr;++m)
	  if (JacVelWeight[m] == NULL) {
	    fprintf(stderr, "Error in grid::ComputeRateOfCompression: velocity derivative undefined.\n");
	    return FAIL;
	  }

	/// compute expansion factor and comoving cellwidth at the current problem time.

	float a = 1.0, dadt = 0.0;
  	float CoCellWidth[MAX_DIMENSION];

	if (ComovingCoordinates)
	{
		
		//
		// for Time=0 the function CosmologyComputeExpansionFactor 
		// doesn't give resonable results for a and dadt
		//
		if(Time>tiny_number)
		{
			if (CosmologyComputeExpansionFactor(Time, &a, &dadt) == FAIL)
			{
				fprintf(stderr, "Error in CosmologyComputeExpansionFactors.\n");
			}
		}

		for (int i = 0; i < GridRank; ++i)
		{
			CoCellWidth[i] = a*CellWidth[i][0];
		}
			
	}
	else
	{
		for (int i = 0; i < GridRank; ++i)
		{
			CoCellWidth[i] = CellWidth[i][0];
		}
			
	}

	/// compute divergende of pressure gradient divided by rho

	if (this->ComputeDivInvRhoGradPressure(RateOfCompression, CoCellWidth) == FAIL)
	{
		fprintf(stderr, "Error in grid->ComputeDivInvRhoGradPressure.\n");
		return FAIL;
	}

	/// compute squared norm of rate-of-strain tensor and voriticity

	float* buf1 = new float[size];
	float* buf2 = new float[size];

	if (this->ComputeRateOfStrainNormSqr(buf1, 0) == FAIL)
	{
		fprintf(stderr, "Error in grid->ComputeRateOfStrainNormSqr.\n");
		return FAIL;
	}
	if (this->ComputeVorticityNormSqr(buf2) == FAIL)
	{
	   fprintf(stderr, "Error in grid->ComputeVorticityNormSqr.\n");
	   return FAIL;
	}

	for (int n=0; n<size; ++n)
	    RateOfCompression[n] += 0.5*(buf1[n] - buf2[n]);

	/// clean up
	
	delete[] buf1;
	delete[] buf2;
	
	/// compute divergence of external acceleration field

	if (this->ComputeAccelerationFieldExternal() == FAIL)
	{
	    fprintf(stderr, "Error in grid->ComputeAccelerationFieldExternal.\n");
	    return FAIL;
	}

        if (AccelerationField[0] != NULL) {
	  
	   float* divacc = new float[size];
	  
	   if (this->ComputeDivAcceleration(divacc, CoCellWidth) == FAIL)
	   {
	      fprintf(stderr, "Error in grid->ComputeDivAcceleration.\n");
	      return FAIL;
	   }
	   
	   for (int n=0; n<size; ++n)
	       RateOfCompression[n] -= divacc[n];
	   
	   delete[] divacc;
	}
        
	/// compute compression by self-gravity

	if (SelfGravity) {
	   if (ProblemType == 70) {
	      for (int n=0; n<size; ++n)
	          RateOfCompression[n] += GravitationalConstant*(BaryonField[DensNum][n] - 1.0);
	   }  else {
	      for (int n=0; n<size; ++n)
	          RateOfCompression[n] += GravitationalConstant*BaryonField[DensNum][n]; // here some changes for the
	   }								                 // cosmological case are needed.
	}
	
	return SUCCESS;
}

/**
 * Computes divergence of the acceleration field
 *
 * @return returns SUCCESS or FAIL
 */
int grid::ComputeDivAcceleration(float* div, float* delta)
{
	if (ProcessorNumber != MyProcessorNumber)
	{
		return SUCCESS;
	}
	
	if(debug) 
	    cout << "[" << MyProcessorNumber << "]" << " Computing divergence of acceleration field " << endl;
	
	int n;

	// compute the gradient for active zones and the ghost cells next to active zones
	// this allows us to take first-order centered differences of gradients for diffusion

	int size = 1;
	int StartIndex[GridRank];
	int EndIndex[GridRank];
	for (int dim = 0; dim < GridRank; dim++) {
	    size *= GridDimension[dim];
	    StartIndex[dim] = GridStartIndex[dim]-1;
	    EndIndex[dim] = GridEndIndex[dim]+1;
	}
		
        switch (GridRank)//HB
        {
         
          case 1:  //GridRank = 1
          {
              for (int k = StartIndex[0]; k <= EndIndex[0]; ++k)
              {
                  n = k;
                  div[n] = deriv4(AccelerationField[0][n+2], AccelerationField[0][n+1],
	                          AccelerationField[0][n-1], AccelerationField[0][n-2], delta[0]);
              }
              break;
          }
          
          case 2: //GridRank = 2
          {
              int numx = GridDimension[0];
          
	      for (int j = StartIndex[1]; j <= EndIndex[1]; ++j)
                  for (int k = StartIndex[0]; k <= EndIndex[0]; ++k)
                  {
	             n = k+j*numx;
	             div[n] = deriv4(AccelerationField[0][n+2], AccelerationField[0][n+1],
	                             AccelerationField[0][n-1], AccelerationField[0][n-2], delta[0]);
 	             div[n] += deriv4(AccelerationField[1][n+2*numx], AccelerationField[1][n+1*numx],
	                              AccelerationField[1][n-1*numx], AccelerationField[1][n-2*numx], delta[1]);
	         }
            break; 
          }
          case 3: //GridRank = 3
          {
	    int numx = GridDimension[0];
	    int numy = GridDimension[1];
            int numxy = numx*numy;
            
	    for (int i = StartIndex[2]; i <= EndIndex[2]; ++i)
		for (int j = StartIndex[1]; j <= EndIndex[1]; ++j)
		    for (int k = StartIndex[0]; k <= EndIndex[0]; ++k)
		    {
			n = k+j*numx+i*numxy;
			div[n] = deriv4(AccelerationField[0][n+2], AccelerationField[0][n+1],
					AccelerationField[0][n-1], AccelerationField[0][n-2], delta[0]);
			div[n] += deriv4(AccelerationField[1][n+2*numx], AccelerationField[1][n+1*numx],
					 AccelerationField[1][n-1*numx], AccelerationField[1][n-2*numx], delta[1]);
			div[n] += deriv4(AccelerationField[2][n+2*numxy], AccelerationField[2][n+1*numxy],
					 AccelerationField[2][n-1*numxy], AccelerationField[2][n-2*numxy], delta[2]);
		    }
	    break;
	  }
	}

	return SUCCESS;
}

/**
 * Computes divergence of the acceleration field
 *
 * @return returns SUCCESS or FAIL
 */
int grid::ComputeDivInvRhoGradPressure(float* div, float* delta)
{
	if (ProcessorNumber != MyProcessorNumber)
	{
		return SUCCESS;
	}
	
	if(debug) 
	    cout << "[" << MyProcessorNumber << "]" << " Computing div (grad P)/rho " << endl;
	
	int n;

	int DensNum;
	if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields)) < 0)
	{
		fprintf(stderr, "Error in grid::ComputeRateOfCompression: Cannot find density.\n");
		return FAIL;
	}


	// compute the gradient for active zones and the ghost cells next to active zones
	// this allows us to take first-order centered differences of gradients for diffusion

	int size = 1;
	int StartIndex[GridRank];
	int EndIndex[GridRank];
	for (int dim = 0; dim < GridRank; dim++) {
	    size *= GridDimension[dim];
	    StartIndex[dim] = GridStartIndex[dim]-1;
	    EndIndex[dim] = GridEndIndex[dim]+1;
	}
		
	int result;
	float* pressure = new float[size];

	result = this->ComputePressure(Time, pressure);
        
	if (result == FAIL) {
	   fprintf(stderr, "Error in grid->ComputePressure.\n");	
	   return FAIL;
	}

        switch (GridRank)//HB
        {
         
          case 1:  //GridRank = 1
          {
	    fprintf(stderr, "Error in grid->ComputeDivInvRhoGradPressure: not implemented for 1D.\n");
	    return FAIL;

	    break;
          }
          
          case 2: //GridRank = 2
          {
	    fprintf(stderr, "Error in grid->ComputeDivInvRhoGradPressure: not implemented for 2D.\n");
	    return FAIL;

            break; 
          }
          case 3: //GridRank = 3
          {
	    int numx = GridDimension[0];
	    int numy = GridDimension[1];
            int numxy = numx*numy;

	    float* tmpx = new float[size];
	    float* tmpy = new float[size];
	    float* tmpz = new float[size];
            
	    for (int i = StartIndex[2]-1; i <= EndIndex[2]+1; ++i)
		for (int j = StartIndex[1]-1; j <= EndIndex[1]+1; ++j)
		    for (int k = StartIndex[0]-1; k <= EndIndex[0]+1; ++k)
		    {
			n = k+j*numx+i*numxy;
			/* central differences
			tmpx[n] = (pressure[n+1] - pressure[n-1])/
			          (delta[0]*(BaryonField[DensNum][n+1] + BaryonField[DensNum][n-1]));
			tmpy[n] = (pressure[n+numx] - pressure[n-numx])/
			          (delta[1]*(BaryonField[DensNum][n+numx] + BaryonField[DensNum][n-numx])); 
			tmpz[n] = (pressure[n+numxy] - pressure[n-numxy])/
			          (delta[2]*(BaryonField[DensNum][n+numxy] + BaryonField[DensNum][n-numxy]));
			*/			
                        // modification of 2nd order central difference with density averaging for adjacent cells
			tmpx[n] = (pressure[n+1] - pressure[n]  ) / (delta[0]*(BaryonField[DensNum][n+1] + BaryonField[DensNum][n]  )) +
			          (pressure[n]   - pressure[n-1]) / (delta[0]*(BaryonField[DensNum][n]   + BaryonField[DensNum][n-1]));
                        tmpy[n] = (pressure[n+numx] - pressure[n]     ) / (delta[1]*(BaryonField[DensNum][n+numx] + BaryonField[DensNum][n]     )) +
			          (pressure[n]      - pressure[n-numx]) / (delta[1]*(BaryonField[DensNum][n]      + BaryonField[DensNum][n-numx])); 
			tmpz[n] = (pressure[n+numxy] - pressure[n]      ) / (delta[1]*(BaryonField[DensNum][n+numxy] + BaryonField[DensNum][n]     )) +
			          (pressure[n]       - pressure[n-numxy]) / (delta[1]*(BaryonField[DensNum][n]       + BaryonField[DensNum][n-numxy])); 
		    }

	    for (int i = StartIndex[2]; i <= EndIndex[2]; ++i)
		for (int j = StartIndex[1]; j <= EndIndex[1]; ++j)
		    for (int k = StartIndex[0]; k <= EndIndex[0]; ++k)
		    {
			n = k+j*numx+i*numxy;
			div[n] =  deriv2(tmpx[n+1],     tmpx[n-1],     delta[0]);
			div[n] += deriv2(tmpy[n+numx],  tmpy[n-numx],  delta[1]);
			div[n] += deriv2(tmpz[n+numxy], tmpz[n-numxy], delta[2]);
		    }

	    delete[] tmpx, tmpy, tmpz;

	    break;
	  }
	}

	delete[] pressure;

	return SUCCESS;
}


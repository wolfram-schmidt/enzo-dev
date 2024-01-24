/***********************************************************************
/
/  GRID CLASS (COMPUTE JACOBIAN OF WEIGHTED VELOCITY)
/
/  written by: Wolfram Schmidt
/  date:       April 2006
/  modified:   February 2022
/
/  PURPOSE: computes Jacobian of the velocity
/           used for structural invariants (vorticity modulus,
/           norm or rate of strain, etc.) in cell flagging methods
/
/  RETURNS: error status
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

int FindField(int f, int farray[], int n);

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);


/**
 * Computes the Jacobian matrix of the velocity field
 *
 * for weighing > 0, the velocity is weighted by the mass density to the power 1/n,
 *    where n = 1 + weighing
 * if the flag fluctuation is set, the averaged velocity field is substracted provided that
 *    the global flag ShearImproved is set (currently not implemented)
 *
 * @return returns SUCCESS or FAIL
 */
  
int grid::ComputeJacobianVelocity(int weighing)
{
	if (ProcessorNumber != MyProcessorNumber)
	{
		return SUCCESS;
	}
	
	if(debug) cout << "[" << MyProcessorNumber << "]" << " Computing Jacobian of the velocity" << endl;

	int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
	if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum)==FAIL)
	{
		fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
		return FAIL;
	}
		
        /*
	int AveMomt1Num, AveMomt2Num, AveMomt3Num;

	if (ShearImproved)
	{
	   if ((AveMomt1Num = FindField(AveMomtDensity1, FieldType, NumberOfBaryonFields)) < 0) 
	   {
	       fprintf(stderr, "Cannot find average x-component of averaged momentum.\n");
	       return FAIL;
	   }
	   if ((AveMomt2Num = FindField(AveMomtDensity2, FieldType, NumberOfBaryonFields)) < 0) 
	   {
	       fprintf(stderr, "Cannot find average y-component of averaged momentum.\n");
	       return FAIL;
	   }
	   if ((AveMomt3Num = FindField(AveMomtDensity3, FieldType, NumberOfBaryonFields)) < 0) 
	   {
	       fprintf(stderr, "Cannot find average z-component of averaged momentum.\n");
	       return FAIL;
	   }
	}
	*/

	int size = 1;
	for (int dim = 0; dim < GridRank; dim++)
	    size *= GridDimension[dim];

	float* weight;

	if (weighing > 0) {

           float expn = 1.0/(1.0+weighing);

	   if (debug) cout << "[" << MyProcessorNumber << "] Computing velocities with mass weighing exponent = " << expn << endl;

	   weight = new float[size];	

	   if (weighing > 1) {
	      for (int i=0;i<size;++i)
                  weight[i] = pow(BaryonField[DensNum][i], expn);
	   } else if (weighing == 1) {
	      for (int i=0;i<size;++i)
	          weight[i] = sqrt(BaryonField[DensNum][i]);
	   }
        }

	/*
	float* inv_density;

	if (ShearImproved && fluctuation) {

	   inv_density = new float[size];

	   for (int i=0;i<size;++i)
	       inv_density[i] = 1.0/max(BaryonField[DensNum][i], tiny_number);
	}
	*/

	float* velocity[MAX_DIMENSION];

    	for (int dim = 0; dim < GridRank; dim++)
	    velocity[dim] = new float[size];

	switch(GridRank)
	{
		case 1:
		{	
		        for(int i=0;i<size;++i)
			{
			    velocity[X][i] = BaryonField[Vel1Num][i];
			}
			/*
			if (ShearImproved && fluctuation)
			    for(int i=0;i<size;++i)
			    {
				velocity[X][i] -= inv_density[i]*BaryonField[AveMomt1Num][i];
			    }
			*/
		        if (weighing > 0)
			    for(int i=0;i<size;++i)
			    {
			        velocity[X][i] *= weight[i];
			    }
			break;
		}
		case 2:
		{
			for(int i=0;i<size;++i)
			{
			    velocity[X][i] = BaryonField[Vel1Num][i];
			    velocity[Y][i] = BaryonField[Vel2Num][i];
			}
			/*
			if (ShearImproved && fluctuation)
			    for(int i=0;i<size;++i)
			    {
				velocity[X][i] -= inv_density[i]*BaryonField[AveMomt1Num][i];
				velocity[Y][i] -= inv_density[i]*BaryonField[AveMomt2Num][i];
			    }
			*/
		        if (weighing > 0)
			    for(int i=0;i<size;++i)
			    {
			        velocity[X][i] *= weight[i];
				velocity[Y][i] *= weight[i];
			    }
			break;
		}
		case 3:
		{
		        for(int i=0;i<size;++i)
			{
			    velocity[X][i] = BaryonField[Vel1Num][i];
			    velocity[Y][i] = BaryonField[Vel2Num][i];
			    velocity[Z][i] = BaryonField[Vel3Num][i];
			}
			/*
			if (ShearImproved && fluctuation)
			    for(int i=0;i<size;++i)
			    {
				velocity[X][i] -= inv_density[i]*BaryonField[AveMomt1Num][i];
				velocity[Y][i] -= inv_density[i]*BaryonField[AveMomt2Num][i];
				velocity[Z][i] -= inv_density[i]*BaryonField[AveMomt3Num][i];
			    }
			*/
		        if (weighing > 0)
			    for(int i=0;i<size;++i)
			    {
			        velocity[X][i] *= weight[i];
				velocity[Y][i] *= weight[i];
				velocity[Z][i] *= weight[i];
			    }
			break;
		}
		default:
		{
			return FAIL;
		}
	}

	if (weighing > 0) delete[] weight;
	/*
	if (ShearImproved && fluctuation) delete[] inv_density;
	*/

	// compute expansion factor and comoving cellwidth at the current problem time.

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

	//Change the start and end index so that the velocity Jacobian is also computed for
	//the first ghost cell next to the boundary of the active area. This ensures 
	//that we can compute symmetric derivatives on a 3-point stencil.
	
	int StartIndex[GridRank];
	int EndIndex[GridRank];
	
	for(int dim=0;dim<GridRank;++dim)
	{
	    StartIndex[dim]=GridStartIndex[dim];
	    EndIndex[dim]=GridEndIndex[dim];
	}
	
	int rank_sqr=GridRank*GridRank;	
		
	for (int m=0;m<rank_sqr;++m)
	{
	    if (JacVelWeight[m] == NULL) {
	       JacVelWeight[m] = new float[size];
	       for (int n=0;n<size;++n) JacVelWeight[m][n] = 0.0;
	    }
	}

	int n=0;

	switch(GridRank)
	{
		case 1:
		{
			float dx=CoCellWidth[0];

			for(int k=StartIndex[X];k<=EndIndex[X];++k)
			{
				n=k;
				
				JacVelWeight[XX][n]=deriv4(velocity[X][n+2],velocity[X][n+1],velocity[X][n-1],velocity[X][n-2],dx); //+dadt/a;
			}
			break;
		}
		case 2:
		{
			float dx=CoCellWidth[0];
			float dy=CoCellWidth[1];
			
			int numx=GridDimension[0];
			
			for(int j=StartIndex[Y];j<=EndIndex[Y];++j)
			{
				for(int k=StartIndex[X];k<=EndIndex[X];++k)
				{
					n=k+j*numx;
					
					JacVelWeight[XX][n]=deriv4(velocity[X][n+2],velocity[X][n+1],velocity[X][n-1],velocity[X][n-2],dx); //+dadt/a;
					JacVelWeight[XY][n]=deriv4(velocity[X][n+2*numx],velocity[X][n+1*numx],velocity[X][n-1*numx],velocity[X][n-2*numx],dy);
				
					JacVelWeight[YX][n]=deriv4(velocity[Y][n+2],velocity[Y][n+1],velocity[Y][n-1],velocity[Y][n-2],dx);
					JacVelWeight[YY][n]=deriv4(velocity[Y][n+2*numx],velocity[Y][n+1*numx],velocity[Y][n-1*numx],velocity[Y][n-2*numx],dy); //+dadt/a;
				}
			}
			
			break;
		}
		case 3:
		{
			float dx=CoCellWidth[0];
			float dy=CoCellWidth[1];
			float dz=CoCellWidth[2];
			
			int numx=GridDimension[0];
			int numy=GridDimension[1];
			int numxy=numx*numy;

			for(int i=StartIndex[Z]; i<=EndIndex[Z];++i)
			{
				for(int j=StartIndex[Y];j<=EndIndex[Y];++j)
				{
					for(int k=StartIndex[X];k<=EndIndex[X];++k)
					{
						n=k+j*numx+i*numxy;
							
						JacVelWeight[XX][n]=deriv4(velocity[X][n+2],velocity[X][n+1],velocity[X][n-1],velocity[X][n-2],dx); //+dadt/a;
						JacVelWeight[XY][n]=deriv4(velocity[X][n+2*numx],velocity[X][n+1*numx],velocity[X][n-1*numx],velocity[X][n-2*numx],dy);
						JacVelWeight[XZ][n]=deriv4(velocity[X][n+2*numxy],velocity[X][n+1*numxy],velocity[X][n-1*numxy],velocity[X][n-2*numxy],dz);
				
						JacVelWeight[YX][n]=deriv4(velocity[Y][n+2],velocity[Y][n+1],velocity[Y][n-1],velocity[Y][n-2],dx);
						JacVelWeight[YY][n]=deriv4(velocity[Y][n+2*numx],velocity[Y][n+1*numx],velocity[Y][n-1*numx],velocity[Y][n-2*numx],dy); //+dadt/a;
						JacVelWeight[YZ][n]=deriv4(velocity[Y][n+2*numxy],velocity[Y][n+1*numxy],velocity[Y][n-1*numxy],velocity[Y][n-2*numxy],dz);
				
						JacVelWeight[ZX][n]=deriv4(velocity[Z][n+2],velocity[Z][n+1],velocity[Z][n-1],velocity[Z][n-2],dx);
						JacVelWeight[ZY][n]=deriv4(velocity[Z][n+2*numx],velocity[Z][n+1*numx],velocity[Z][n-1*numx],velocity[Z][n-2*numx],dy);
						JacVelWeight[ZZ][n]=deriv4(velocity[Z][n+2*numxy],velocity[Z][n+1*numxy],velocity[Z][n-1*numxy],velocity[Z][n-2*numxy],dz); //+dadt/a;
					}
				}
			}
			break;
		}
		default:
		{
			return FAIL;
		}
	}

	for (int dim = 0; dim < GridRank; dim++)
	    delete[] velocity[dim];

	return SUCCESS;
}

/**
 * Computes the subgrid-scale energy in the nonlinear structural model
 * from the square of the norm of the Jacobian of the velocity field
 *
 * @return returns SUCCESS or FAIL
 */

int grid::ComputeNonLinearSGSEnergy(float* SGSEnergy)
{
	if (ProcessorNumber != MyProcessorNumber)
	{
		return SUCCESS;
	}
	
	if(debug) cout << "[" << MyProcessorNumber << "]" << " Computing nonlinear SGS energy" << endl;

	int size = 1;
	for (int dim = 0; dim < GridRank; dim++)
	    size *= GridDimension[dim];

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

	// the combined prefactor
	float CDeltaSqr = 0.5/12. * SGScoeffNLu * POW(CoCellWidth[0]*CoCellWidth[1]*CoCellWidth[2],2./3.);

	if (this->ComputeJacobianVelocityNormSqr(SGSEnergy) == FAIL)
	{
	    fprintf(stderr, "Error in grid->ComputeJacobianVelocityNormSqr.\n");
	    return FAIL;
	}
	
	for(int i=0;i<size;++i)
	{
	    SGSEnergy[i] *= CDeltaSqr;
	}

	return SUCCESS;
}

/**
 * Computes the square of the norm of the Jacobian of the velocity field
 *
 * @return returns SUCCESS or FAIL
 */

int grid::ComputeJacobianVelocityNormSqr(float* JacVelWeightNormSqr)
{
	if (ProcessorNumber != MyProcessorNumber)
	{
		return SUCCESS;
	}
	
	if(debug) cout << "[" << MyProcessorNumber << "]" << " Computing squared norm of velocity derivative" << endl;

	int size = 1;
	for (int dim = 0; dim < GridRank; dim++)
	    size *= GridDimension[dim];

	switch(GridRank)
	{
		case 1:
		{	
			for(int i=0;i<size;++i)
			{
				JacVelWeightNormSqr[i] = 2.0*JacVelWeight[XX][i]*JacVelWeight[XX][i];
			}
			
			break;
		}
		case 2:
		{
		
			for(int i=0;i<size;++i)
			{
				JacVelWeightNormSqr[i] = 2.0*(JacVelWeight[XX][i]*JacVelWeight[XX][i]+JacVelWeight[XY][i]*JacVelWeight[XY][i] +
				                        JacVelWeight[YX][i]*JacVelWeight[YX][i]+JacVelWeight[YY][i]*JacVelWeight[YY][i]);
			}
			
			break;
		}
		case 3:
		{
			for(int i=0;i<size;++i)
			{
			
				JacVelWeightNormSqr[i] =
				   2.0*(JacVelWeight[XX][i]*JacVelWeight[XX][i]+JacVelWeight[XY][i]*JacVelWeight[XY][i]+JacVelWeight[XZ][i]*JacVelWeight[XZ][i]+
				        JacVelWeight[YX][i]*JacVelWeight[YX][i]+JacVelWeight[YY][i]*JacVelWeight[YY][i]+JacVelWeight[YZ][i]*JacVelWeight[YZ][i]+
				        JacVelWeight[ZX][i]*JacVelWeight[ZX][i]+JacVelWeight[ZY][i]*JacVelWeight[ZY][i]+JacVelWeight[ZZ][i]*JacVelWeight[ZZ][i]);
			}
			
			break;
		}
		default:
		{
			return FAIL;
		}
	}
	return SUCCESS;
}

/**
 * Computes the determinate of the Jacobian matrix of the velocity field
 *
 * @return returns SUCCESS or FAIL
 */

int grid::ComputeNonLinearScalar(float* NonLinScalar)
{
	if (ProcessorNumber != MyProcessorNumber)
	{
		return SUCCESS;
	}
	
	if(debug) cout << "[" << MyProcessorNumber << "]" << " Computing nonlinear scalar" << endl;

	int size = 1;
	for (int dim = 0; dim < GridRank; dim++)
	    size *= GridDimension[dim];

	if (this->ComputeJacobianVelocityNormSqr(NonLinScalar) == FAIL)
	{
	   fprintf(stderr, "Error in grid->ComputeJacobianVelocityNormSqr\n");
	   return FAIL;
	}

	switch(GridRank)
	{
		case 1:
                {
			for(int n=0;n<size;++n)
			{
			 	NonLinScalar[n] = 0.0;
			}
			break;

                } 
		case 2:
                {
			for(int n=0;n<size;++n)
			{
				NonLinScalar[n] = 2.0*
				  ((JacVelWeight[XX][n]*JacVelWeight[XX][n]+
				    JacVelWeight[XY][n]*JacVelWeight[XY][n])*(JacVelWeight[XX][n]-JacVelWeight[YY][n]) +
				   (JacVelWeight[YX][n]*JacVelWeight[YX][n]+
				    JacVelWeight[YY][n]*JacVelWeight[YY][n])*(JacVelWeight[YY][n]-JacVelWeight[XX][n]) +
				   (JacVelWeight[XX][n]*JacVelWeight[YX][n]+
				    JacVelWeight[XY][n]*JacVelWeight[YY][n])*(JacVelWeight[XY][n]+JacVelWeight[YX][n]))/
				  max(NonLinScalar[n], tiny_number);
			}
			break;
                }
		case 3:
		{
			float fct = 1.0/3.0;

			for(int n=0;n<size;++n)
			{
				NonLinScalar[n] = 2.0*
				  ((JacVelWeight[XX][n]*JacVelWeight[XX][n]+
				    JacVelWeight[XY][n]*JacVelWeight[XY][n]+ 
				    JacVelWeight[XZ][n]*JacVelWeight[XZ][n])*fct*(2.0*JacVelWeight[XX][n]-JacVelWeight[YY][n]-JacVelWeight[ZZ][n]) +
				   (JacVelWeight[YX][n]*JacVelWeight[YX][n]+
				    JacVelWeight[YY][n]*JacVelWeight[YY][n]+ 
				    JacVelWeight[YZ][n]*JacVelWeight[YZ][n])*fct*(2.0*JacVelWeight[YY][n]-JacVelWeight[XX][n]-JacVelWeight[ZZ][n]) +
				   (JacVelWeight[ZX][n]*JacVelWeight[ZX][n]+
				    JacVelWeight[ZY][n]*JacVelWeight[ZY][n]+ 
				    JacVelWeight[ZZ][n]*JacVelWeight[ZZ][n])*fct*(2.0*JacVelWeight[ZZ][n]-JacVelWeight[XX][n]-JacVelWeight[YY][n]) +
				   (JacVelWeight[XX][n]*JacVelWeight[YX][n]+
				    JacVelWeight[XY][n]*JacVelWeight[YY][n]+ 
				    JacVelWeight[XZ][n]*JacVelWeight[YZ][n])*(JacVelWeight[XY][n]+JacVelWeight[YX][n]) +
				   (JacVelWeight[YX][n]*JacVelWeight[ZX][n]+
				    JacVelWeight[YY][n]*JacVelWeight[ZY][n]+ 
				    JacVelWeight[YZ][n]*JacVelWeight[ZZ][n])*(JacVelWeight[YZ][n]+JacVelWeight[ZY][n]) +
				   (JacVelWeight[ZX][n]*JacVelWeight[XX][n]+
				    JacVelWeight[ZY][n]*JacVelWeight[XY][n]+ 
				    JacVelWeight[ZZ][n]*JacVelWeight[XZ][n])*(JacVelWeight[ZX][n]+JacVelWeight[XZ][n]))/
				  max(NonLinScalar[n], tiny_number);
			}
			break;
		}
                default:
		{
			return FAIL;
		}
	}
	return SUCCESS;
}


/**
 * Computes the norm of the rate-of-strain scalar of the velocity field 
 * from the Jacobian of the velocity field
 * 
 * if flag tracefree is set, the norm of the trace-free part of the tensor is computed
 * 
 * @return returns SUCCESS or FAIL
 */

int grid::ComputeRateOfStrainNormSqr(float* STNormSqr, int tracefree)
{
	if (ProcessorNumber != MyProcessorNumber)
	{
		return SUCCESS;
	}
	
        if(debug) cout << "[" << MyProcessorNumber << "]" << " Computing squared norm of rate of strain" << endl;

	int size = 1;
	for (int dim = 0; dim < GridRank; dim++)
	    size *= GridDimension[dim];

	switch(GridRank)
	{
		case 1:
		{
			if (tracefree)
			{
				for(int n=0;n<size;++n)
				{
					STNormSqr[n] = 0.0;
				}
			} 
			else 
			{
				for(int n=0;n<size;++n)
				{
					STNormSqr[n] = 2.0*JacVelWeight[XX][n]*JacVelWeight[XX][n];
				}
			}
			
			break;
		}
		case 2:
		{
			if (tracefree)
			{
				for(int n=0;n<size;++n)
				{
					STNormSqr[n] = (JacVelWeight[XX][n]-JacVelWeight[YY][n])*(JacVelWeight[XX][n]-JacVelWeight[YY][n]) +
					               (JacVelWeight[XY][n]+JacVelWeight[YX][n])*(JacVelWeight[XY][n]+JacVelWeight[YX][n]);
				}
			}
			else 
			{
				for(int n=0;n<size;++n)
				{
					STNormSqr[n] = 2.0*(JacVelWeight[XX][n]*JacVelWeight[XX][n]+JacVelWeight[YY][n]*JacVelWeight[YY][n]) +
					               (JacVelWeight[XY][n]+JacVelWeight[YX][n])*(JacVelWeight[XY][n]+JacVelWeight[YX][n]);
				}
			}
			
			break;
		}
		case 3:
		{
			if (tracefree)
			{
				for(int n=0;n<size;++n)
				{
					STNormSqr[n] = 4.0*(JacVelWeight[XX][n]*JacVelWeight[XX][n] +
					                    JacVelWeight[YY][n]*JacVelWeight[YY][n] +
					                    JacVelWeight[ZZ][n]*JacVelWeight[ZZ][n] -
					                    JacVelWeight[XX][n]*JacVelWeight[YY][n] -
					                    JacVelWeight[YY][n]*JacVelWeight[ZZ][n] -
					                    JacVelWeight[ZZ][n]*JacVelWeight[XX][n])/3.0 +
					               (JacVelWeight[XY][n]+JacVelWeight[YX][n])*(JacVelWeight[XY][n]+JacVelWeight[YX][n]) +
					               (JacVelWeight[YZ][n]+JacVelWeight[ZY][n])*(JacVelWeight[YZ][n]+JacVelWeight[ZY][n]) +
					               (JacVelWeight[XZ][n]+JacVelWeight[ZX][n])*(JacVelWeight[XZ][n]+JacVelWeight[ZX][n]);
				}
			}
			else
			{
				for(int n=0;n<size;++n)
				{
					STNormSqr[n] = 2.0*(JacVelWeight[XX][n]*JacVelWeight[XX][n] +
					                    JacVelWeight[YY][n]*JacVelWeight[YY][n] +
					                    JacVelWeight[ZZ][n]*JacVelWeight[ZZ][n]) +
					               (JacVelWeight[XY][n]+JacVelWeight[YX][n])*(JacVelWeight[XY][n]+JacVelWeight[YX][n]) +
					               (JacVelWeight[YZ][n]+JacVelWeight[ZY][n])*(JacVelWeight[YZ][n]+JacVelWeight[ZY][n]) +
					               (JacVelWeight[XZ][n]+JacVelWeight[ZX][n])*(JacVelWeight[XZ][n]+JacVelWeight[ZX][n]);
				}
			}
			break;
		}
		default:
		{
			return FAIL;
		}
	}
	return SUCCESS;
}

/**
 * Computes the norm of the vorticity of the velocity field 
 * from the Jacobian of the velocity field
 *
 * @return returns SUCCESS or FAIL
 */

int grid::ComputeVorticityNormSqr(float* VortNormSqr)
{
	if (ProcessorNumber != MyProcessorNumber)
	{
		return SUCCESS;
	}
	
	if(debug) cout << "[" << MyProcessorNumber << "]" << " Computing squared norm of vorticity" << endl;

	int size = 1;
	for (int dim = 0; dim < GridRank; dim++)
	    size *= GridDimension[dim];

	switch(GridRank)
	{
		case 1:
		{
			for(int n=0;n<size;++n)
			{
				VortNormSqr[n] = 0.0;
			}
			break;
		}
		case 2:
		{
			for(int n=0;n<size;++n)
			{
				VortNormSqr[n] = (JacVelWeight[YX][n]-JacVelWeight[XY][n])*(JacVelWeight[YX][n]-JacVelWeight[XY][n]);
			}
			
			break;
		}
		case 3:
		{
			for(int n=0;n<size;++n)
			{
				VortNormSqr[n] = (JacVelWeight[ZY][n]-JacVelWeight[YZ][n])*(JacVelWeight[ZY][n]-JacVelWeight[YZ][n])+
				                 (JacVelWeight[XZ][n]-JacVelWeight[ZX][n])*(JacVelWeight[XZ][n]-JacVelWeight[ZX][n])+
				                 (JacVelWeight[YX][n]-JacVelWeight[XY][n])*(JacVelWeight[YX][n]-JacVelWeight[XY][n]);
			}
			break;
		}
		default:
		{
			return FAIL;
		}
	}
	return SUCCESS;
}

/**
 * Computes the divergence of the velocity field 
 * from the Jacobian of the velocity field
 *
 * @return returns SUCCESS or FAIL
 */

int grid::ComputeDivergence(float* Div)
{
	if (ProcessorNumber != MyProcessorNumber)
	{
		return SUCCESS;
	}
	
	if(debug) cout << "[" << MyProcessorNumber << "]" << " Computing divergence" << endl;

	int size = 1;
	for (int dim = 0; dim < GridRank; dim++)
	    size *= GridDimension[dim];

	switch(GridRank)
	{
		case 1:
		{
			for(int n=0;n<size;++n)
			{
				Div[n] = JacVelWeight[XX][n];
			}
			break;
		}
		case 2:
		{
			for(int n=0;n<size;++n)
			{
				Div[n] = JacVelWeight[XX][n]+JacVelWeight[YY][n];
			}
			
			break;
		}
		case 3:
		{
			for(int n=0;n<size;++n)
			{
				Div[n] = JacVelWeight[XX][n]+JacVelWeight[YY][n]+JacVelWeight[ZZ][n];
			}
			break;
		}
		default:
		{
			return FAIL;
		}
	}
	return SUCCESS;
}

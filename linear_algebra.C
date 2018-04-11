//=============================================================
// 
//  linear_algebra.cpp
//  
//  Function to do routine linear algebra operations.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linear_algebra.h"
#include "grid.h"
#include "Defined_Variables.h"
#include "sys_mem.h"
#include "params.h"


//=============================================================
// 
//  crossProduct()
//
//  Computes the cross product of two vectors. Assumes 3D.
//
//=============================================================
void crossProduct ( double *u, double *v, double *w )
{
  w[0] = u[1]*v[2] - u[2]*v[1];
  w[1] = u[2]*v[0] - u[0]*v[2];
  w[2] = u[0]*v[1] - u[1]*v[0];
  return;
}


//=============================================================
// 
//  dotProduct()
//
//  Computes the dot product of two vectors.
//
//=============================================================
double dotProduct ( double *u, double *v )
{
  int i;                             // Loop counter.
  double dot = 0.;                   // Dot product result.
  

  for ( i=0; i < NDIM; i++ )
    {
      dot += u[i]*v[i];
    }

  return dot;
}


//=============================================================
// 
//  dotProductGeneral()
//
//  Computes the dot product of two vectors.
//
//=============================================================
double dotProductGeneral ( double *u, double *v, int n )
{
  int i;                             // Loop counter.
  double dot = 0.;                   // Dot product result.
  

  for ( i=0; i < n; i++ )
    {
      dot += u[i]*v[i];
    }

  return dot;
}

//=============================================================
// 
//  normalizeVectorGeneral()
//
//  Normalizes a vector of arbitrary length.
//
//=============================================================
void normalizeVectorGeneral ( double *u, int n )
{
  int i;                             // Loop counter.
  double mag = 0.;                   // Vector magnitude.

  for ( i=0; i < n; i++ )
    {
      mag += (u[i]*u[i]);
    }

  mag = sqrt(mag);

  for ( i=0; i < n; i++ )
    {
      u[i] = u[i] / mag;
    }

  return;
}




//=============================================================
// 
//  linear_solve()
//
//  Solves the linear system with a Symmetric Point Iterative scheme.
//
//  GRID *grid;                              // The grid.
//  SYS_MEM *smem;                           // The problem space.
//  PARAMS p;                                // The parameters.
//
//=============================================================

double linear_solve ( GRID *grid, SYS_MEM *smem, PARAMS p )
{
  int i, j, k;                         //Loop counters.

  int d;                               //Diagonal index for each entry.
  double r[4];                         //Space to hold the residuals.

  double *rhs1 = NULL;                 //Space to hold the RHS each iteration.

  double RMS = 10.;                    //RMS error.
  double *M = NULL;                    //Space to copy the block diagonal matrix into to run LU factorization.
  double *B = NULL;                    //Space for the substitution  and solve block.
  double *X = NULL;
  double *Y = NULL;
  double *dX = NULL;


  // Pointers.
  int nn = grid->nn;
  double *RHS = smem->RHS;
  double *LHS = smem->LHS;
  int *ia = smem->ia;
  int *iau = smem->iau;
  int *ja = smem->ja;
  double *dQ = grid->dQ;

  // If the run is set to explicit, we can optimize the solution procedure.
  if ( p.method == 0 )
    {
      RMS = explicit_update( grid, smem, p );
      return RMS;
    }

  // Allocate space for the vectors.
  M = (double*)malloc(NUM_VAR*NUM_VAR*sizeof(double));
  if ( M == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'M'.\n"); exit(1); }

  B = (double*)malloc(NUM_VAR*sizeof(double));
  if ( B == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'B'.\n"); exit(1); }

  X = (double*)malloc(NUM_VAR*sizeof(double));
  if ( X == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'X'.\n"); exit(1); }

  Y = (double*)malloc(NUM_VAR*sizeof(double));
  if ( Y == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'Y'.\n"); exit(1); }

  dX = (double*)malloc(NUM_VAR*sizeof(double));
  if ( dX == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'dX'.\n"); exit(1); }
  
  rhs1 = (double*)malloc((nn+1)*NUM_VAR*sizeof(double));
  if ( rhs1 == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'rhs1'.\n"); exit(1); }

  //Now perform the loops until we converge or run out of iterations.
  for ( i=0; i < p.subits; i++ )
    {
      // Do Symmetric GS.
      if ( i % 2 == 0 )
	{
	  RMS = 0.;
	  
	  //Copy the right hand side over to a place where we can change it.
	  for ( j=1; j <= nn; j++ )
	    {
	      for ( k=0; k < NUM_VAR; k++ )
		{
		  rhs1[j*NUM_VAR+k] = RHS[j*NUM_VAR+k];
		}
	    }
	  
	  //Now loop through all the rows.
	  for ( j=1; j <= nn; j++ )
	    {
	      //Get the index in JA of the diagonal.
	      d = iau[j];

	      //Loop over the off diagonal node and do the matrix-vector multiplication and transfer
	      //the result to the RHS.
	      for ( k=ia[j]; k < ia[j+1]; k++ )
		{
		  //Skip the diagonal entry.
		  if ( ja[k] == j )  continue;
		  
		  rhs1[j*NUM_VAR+0] -= ( LHS[k*NUM_VAR*NUM_VAR+0*NUM_VAR+0]*dQ[(ja[k])*NUM_VAR+0] +
					 LHS[k*NUM_VAR*NUM_VAR+0*NUM_VAR+1]*dQ[(ja[k])*NUM_VAR+1] +
					 LHS[k*NUM_VAR*NUM_VAR+0*NUM_VAR+2]*dQ[(ja[k])*NUM_VAR+2] +
					 LHS[k*NUM_VAR*NUM_VAR+0*NUM_VAR+3]*dQ[(ja[k])*NUM_VAR+3] );

		  rhs1[j*NUM_VAR+1] -= ( LHS[k*NUM_VAR*NUM_VAR+1*NUM_VAR+0]*dQ[(ja[k])*NUM_VAR+0] +
					 LHS[k*NUM_VAR*NUM_VAR+1*NUM_VAR+1]*dQ[(ja[k])*NUM_VAR+1] +
					 LHS[k*NUM_VAR*NUM_VAR+1*NUM_VAR+2]*dQ[(ja[k])*NUM_VAR+2] +
					 LHS[k*NUM_VAR*NUM_VAR+1*NUM_VAR+3]*dQ[(ja[k])*NUM_VAR+3] );

		  rhs1[j*NUM_VAR+2] -= ( LHS[k*NUM_VAR*NUM_VAR+2*NUM_VAR+0]*dQ[(ja[k])*NUM_VAR+0] +
					 LHS[k*NUM_VAR*NUM_VAR+2*NUM_VAR+1]*dQ[(ja[k])*NUM_VAR+1] +
					 LHS[k*NUM_VAR*NUM_VAR+2*NUM_VAR+2]*dQ[(ja[k])*NUM_VAR+2] +
					 LHS[k*NUM_VAR*NUM_VAR+2*NUM_VAR+3]*dQ[(ja[k])*NUM_VAR+3] );

		  rhs1[j*NUM_VAR+3] -= ( LHS[k*NUM_VAR*NUM_VAR+3*NUM_VAR+0]*dQ[(ja[k])*NUM_VAR+0] +
					 LHS[k*NUM_VAR*NUM_VAR+3*NUM_VAR+1]*dQ[(ja[k])*NUM_VAR+1] +
					 LHS[k*NUM_VAR*NUM_VAR+3*NUM_VAR+2]*dQ[(ja[k])*NUM_VAR+2] +
					 LHS[k*NUM_VAR*NUM_VAR+3*NUM_VAR+3]*dQ[(ja[k])*NUM_VAR+3] );

		}
	      
	      // Copy the diagonal matrix out.
	      for ( k=0; k < 4; k++ )
		{
		  M[k*NUM_VAR+0] = LHS[d*NUM_VAR*NUM_VAR+k*NUM_VAR+0];
		  M[k*NUM_VAR+1] = LHS[d*NUM_VAR*NUM_VAR+k*NUM_VAR+1];
		  M[k*NUM_VAR+2] = LHS[d*NUM_VAR*NUM_VAR+k*NUM_VAR+2];
		  M[k*NUM_VAR+3] = LHS[d*NUM_VAR*NUM_VAR+k*NUM_VAR+3];
		}
	      
	      // Do the LU factorization.
	      lu(M);

	      // Copy the right hand side over.
	      B[0] = rhs1[j*NUM_VAR+0];
	      B[1] = rhs1[j*NUM_VAR+1];
	      B[2] = rhs1[j*NUM_VAR+2];
	      B[3] = rhs1[j*NUM_VAR+3];
	      
	      // Perform the backwards and forwards substitution.
	      forwardstab ( M , Y , B );
	      backstab ( M , X , Y );
	  
	      // Get the delta between current and previous solution.
	      dX[0] = X[0] - dQ[j*NUM_VAR+0];
	      dX[1] = X[1] - dQ[j*NUM_VAR+1];
	      dX[2] = X[2] - dQ[j*NUM_VAR+2];
	      dX[3] = X[3] - dQ[j*NUM_VAR+3];
	  
	      // Relax the delta and add it to the previous solution.
	      dQ[j*NUM_VAR+0] += dX[0];
	      dQ[j*NUM_VAR+1] += dX[1];
	      dQ[j*NUM_VAR+2] += dX[2];
	      dQ[j*NUM_VAR+3] += dX[3];	  
	    }
	}

      else
	{
	  RMS = 0.;
	  
	  //Copy the right hand side over to a place where we can change it.
	  for ( j=1; j <= nn; j++ )
	    {
	      for ( k=0; k < 4; k++ )
		{
		  rhs1[j*NUM_VAR+k] = RHS[j*NUM_VAR+k];
		}
	    }
	  
	  //Now loop through all the rows.
	  for ( j=nn; j >= 1; j-- )
	    {
	      //Get the index in JA of the diagonal.
	      d = iau[j];
	      
	      //Loop over the off diagonal node and do the matrix-vector multiplication and transfer
	      //the result to the RHS.
	      for ( k=ia[j]; k < ia[j+1]; k++ )
		{
		  //Skip the diagonal entry.
		  if ( ja[k] == j )  continue;

		  rhs1[j*NUM_VAR+0] -= ( LHS[k*NUM_VAR*NUM_VAR+0*NUM_VAR+0]*dQ[(ja[k])*NUM_VAR+0] +
					 LHS[k*NUM_VAR*NUM_VAR+0*NUM_VAR+1]*dQ[(ja[k])*NUM_VAR+1] +
					 LHS[k*NUM_VAR*NUM_VAR+0*NUM_VAR+2]*dQ[(ja[k])*NUM_VAR+2] +
					 LHS[k*NUM_VAR*NUM_VAR+0*NUM_VAR+3]*dQ[(ja[k])*NUM_VAR+3] );

		  rhs1[j*NUM_VAR+1] -= ( LHS[k*NUM_VAR*NUM_VAR+1*NUM_VAR+0]*dQ[(ja[k])*NUM_VAR+0] +
					 LHS[k*NUM_VAR*NUM_VAR+1*NUM_VAR+1]*dQ[(ja[k])*NUM_VAR+1] +
					 LHS[k*NUM_VAR*NUM_VAR+1*NUM_VAR+2]*dQ[(ja[k])*NUM_VAR+2] +
					 LHS[k*NUM_VAR*NUM_VAR+1*NUM_VAR+3]*dQ[(ja[k])*NUM_VAR+3] );

		  rhs1[j*NUM_VAR+2] -= ( LHS[k*NUM_VAR*NUM_VAR+2*NUM_VAR+0]*dQ[(ja[k])*NUM_VAR+0] +
					 LHS[k*NUM_VAR*NUM_VAR+2*NUM_VAR+1]*dQ[(ja[k])*NUM_VAR+1] +
					 LHS[k*NUM_VAR*NUM_VAR+2*NUM_VAR+2]*dQ[(ja[k])*NUM_VAR+2] +
					 LHS[k*NUM_VAR*NUM_VAR+2*NUM_VAR+3]*dQ[(ja[k])*NUM_VAR+3] );

		  rhs1[j*NUM_VAR+3] -= ( LHS[k*NUM_VAR*NUM_VAR+3*NUM_VAR+0]*dQ[(ja[k])*NUM_VAR+0] +
					 LHS[k*NUM_VAR*NUM_VAR+3*NUM_VAR+1]*dQ[(ja[k])*NUM_VAR+1] +
					 LHS[k*NUM_VAR*NUM_VAR+3*NUM_VAR+2]*dQ[(ja[k])*NUM_VAR+2] +
					 LHS[k*NUM_VAR*NUM_VAR+3*NUM_VAR+3]*dQ[(ja[k])*NUM_VAR+3] );

		}
	      
	      // Copy the diagonal matrix out.
	      for ( k=0; k < 4; k++ )
		{
		  M[k*NUM_VAR+0] = LHS[d*NUM_VAR*NUM_VAR+k*NUM_VAR+0];
		  M[k*NUM_VAR+1] = LHS[d*NUM_VAR*NUM_VAR+k*NUM_VAR+1];
		  M[k*NUM_VAR+2] = LHS[d*NUM_VAR*NUM_VAR+k*NUM_VAR+2];
		  M[k*NUM_VAR+3] = LHS[d*NUM_VAR*NUM_VAR+k*NUM_VAR+3];
		}
	      
	      // Do the LU factorization.
	      lu(M);

	      // Copy the right hand side over.
	      B[0] = rhs1[j*NUM_VAR+0];
	      B[1] = rhs1[j*NUM_VAR+1];
	      B[2] = rhs1[j*NUM_VAR+2];
	      B[3] = rhs1[j*NUM_VAR+3];
	      
	      // Perform the backwards and forwards substitution.
	      forwardstab ( M , Y , B );
	      backstab ( M , X , Y );
	  
	      // Get the delta between current and previous solution.
	      dX[0] = X[0] - dQ[j*NUM_VAR+0];
	      dX[1] = X[1] - dQ[j*NUM_VAR+1];
	      dX[2] = X[2] - dQ[j*NUM_VAR+2];
	      dX[3] = X[3] - dQ[j*NUM_VAR+3];
	  
	      // Relax the delta and add it to the previous solution.
	      dQ[j*NUM_VAR+0] += dX[0];
	      dQ[j*NUM_VAR+1] += dX[1];
	      dQ[j*NUM_VAR+2] += dX[2];
	      dQ[j*NUM_VAR+3] += dX[3];	  
	    }
	}

      // Compute the residual b-Ax. Start by copying the right hand side to a temp array. Then compute
      // A*x using the solution from this time step.
      for ( j=1; j <= nn; j++ )
	{
	  r[0] = RHS[j*NUM_VAR+0];
	  r[1] = RHS[j*NUM_VAR+1];
	  r[2] = RHS[j*NUM_VAR+2];
	  r[3] = RHS[j*NUM_VAR+3];

	  // Loop over the row for cv 'j' and do matrix-vector multiplication and transfer
	  // the result to the RHS.
	  for ( k=ia[j]; k < ia[j+1]; k++ )
	    {  
	      r[0] -= ( LHS[k*NUM_VAR*NUM_VAR+0*NUM_VAR+0]*dQ[(ja[k])*NUM_VAR+0] +
			LHS[k*NUM_VAR*NUM_VAR+0*NUM_VAR+1]*dQ[(ja[k])*NUM_VAR+1] +
			LHS[k*NUM_VAR*NUM_VAR+0*NUM_VAR+2]*dQ[(ja[k])*NUM_VAR+2] +
			LHS[k*NUM_VAR*NUM_VAR+0*NUM_VAR+3]*dQ[(ja[k])*NUM_VAR+3] );

	      r[1] -= ( LHS[k*NUM_VAR*NUM_VAR+1*NUM_VAR+0]*dQ[(ja[k])*NUM_VAR+0] +
			LHS[k*NUM_VAR*NUM_VAR+1*NUM_VAR+1]*dQ[(ja[k])*NUM_VAR+1] +
			LHS[k*NUM_VAR*NUM_VAR+1*NUM_VAR+2]*dQ[(ja[k])*NUM_VAR+2] +
			LHS[k*NUM_VAR*NUM_VAR+1*NUM_VAR+3]*dQ[(ja[k])*NUM_VAR+3] );

	      r[2] -= ( LHS[k*NUM_VAR*NUM_VAR+2*NUM_VAR+0]*dQ[(ja[k])*NUM_VAR+0] +
			LHS[k*NUM_VAR*NUM_VAR+2*NUM_VAR+1]*dQ[(ja[k])*NUM_VAR+1] +
			LHS[k*NUM_VAR*NUM_VAR+2*NUM_VAR+2]*dQ[(ja[k])*NUM_VAR+2] +
			LHS[k*NUM_VAR*NUM_VAR+2*NUM_VAR+3]*dQ[(ja[k])*NUM_VAR+3] );

	      r[3] -= ( LHS[k*NUM_VAR*NUM_VAR+3*NUM_VAR+0]*dQ[(ja[k])*NUM_VAR+0] +
			LHS[k*NUM_VAR*NUM_VAR+3*NUM_VAR+1]*dQ[(ja[k])*NUM_VAR+1] +
			LHS[k*NUM_VAR*NUM_VAR+3*NUM_VAR+2]*dQ[(ja[k])*NUM_VAR+2] +
			LHS[k*NUM_VAR*NUM_VAR+3*NUM_VAR+3]*dQ[(ja[k])*NUM_VAR+3] );
	      
	    }
	  
	  RMS += (r[0]*r[0] + r[1]*r[1] + r[2]*r[2] + r[3]*r[3]);
	}

      RMS /= (1.*nn);
      RMS = sqrt ( RMS );

      if ( RMS < 1.0E-15 ) break;
    }

  printf("  LINEAR_SOLVE -> Loop %d, RMS=%.15e\n",i+1,RMS);

  freenull(rhs1);
  freenull(M);
  freenull(B);
  freenull(X);
  freenull(Y);
  freenull(dX);

  return RMS;

}


//=============================================================
// 
//  explicit_update()
//
//  Update the system with explicit Euler.
//
//  GRID *grid;                              // The grid.
//  SYS_MEM *smem;                           // The problem space.
//  PARAMS p;                                // The parameters.
//
//=============================================================

double explicit_update ( GRID *grid, SYS_MEM *smem, PARAMS p )
{
  int j, k;                            //Loop counters.

  int d;                               //Diagonal index for each entry.
  double r[4];                         //Space to hold the residuals.
  double *rhs1 = NULL;                 //Space to hold the RHS each iteration.

  double RMS = 10.;                    //RMS error.
  
  double *M = NULL;                    //Space to copy the block diagonal matrix into to run LU factorization.

  double *B = NULL;                    //Space for the substitution  and solve block.
  double *X = NULL;
  double *Y = NULL;
  double *dX = NULL;


  // Pointers.
  int nn = grid->nn;
  double *RHS = smem->RHS;
  double *LHS = smem->LHS;
  int *iau = smem->iau;
  double *dQ = grid->dQ;

  // Allocate space for the vectors.
  M = (double*)malloc(NUM_VAR*NUM_VAR*sizeof(double));
  if ( M == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'M'.\n"); exit(1); }

  B = (double*)malloc(NUM_VAR*sizeof(double));
  if ( B == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'B'.\n"); exit(1); }

  X = (double*)malloc(NUM_VAR*sizeof(double));
  if ( X == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'X'.\n"); exit(1); }

  Y = (double*)malloc(NUM_VAR*sizeof(double));
  if ( Y == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'Y'.\n"); exit(1); }

  dX = (double*)malloc(NUM_VAR*sizeof(double));
  if ( dX == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'dX'.\n"); exit(1); }
  
  rhs1 = (double*)malloc((nn+1)*NUM_VAR*sizeof(double));
  if ( rhs1 == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'rhs1'.\n"); exit(1); }

  RMS = 0.;
	  
  //Copy the right hand side over to a place where we can change it.
  for ( j=1; j <= nn; j++ )
    {
      for ( k=0; k < NUM_VAR; k++ )
	{
	  rhs1[j*NUM_VAR+k] = RHS[j*NUM_VAR+k];
	}
    }
  
  //Now loop through all the rows.
  for ( j=1; j <= nn; j++ )
    {
      //Get the index in JA of the diagonal.
      d = iau[j];
  
      // Copy the diagonal matrix out.
      for ( k=0; k < 4; k++ )
	{
	  M[k*NUM_VAR+0] = LHS[d*NUM_VAR*NUM_VAR+k*NUM_VAR+0];
	  M[k*NUM_VAR+1] = LHS[d*NUM_VAR*NUM_VAR+k*NUM_VAR+1];
	  M[k*NUM_VAR+2] = LHS[d*NUM_VAR*NUM_VAR+k*NUM_VAR+2];
	  M[k*NUM_VAR+3] = LHS[d*NUM_VAR*NUM_VAR+k*NUM_VAR+3];
	}
	      
      // Do the LU factorization.
      lu(M);

      // Copy the right hand side over.
      B[0] = rhs1[j*NUM_VAR+0];
      B[1] = rhs1[j*NUM_VAR+1];
      B[2] = rhs1[j*NUM_VAR+2];
      B[3] = rhs1[j*NUM_VAR+3];
	      
      // Perform the backwards and forwards substitution.
      forwardstab ( M , Y , B );
      backstab ( M , X , Y );
	  
      // Get the delta between current and previous solution.
      dX[0] = X[0] - dQ[j*NUM_VAR+0];
      dX[1] = X[1] - dQ[j*NUM_VAR+1];
      dX[2] = X[2] - dQ[j*NUM_VAR+2];
      dX[3] = X[3] - dQ[j*NUM_VAR+3];
	  
      // Add it to the previous solution.
      dQ[j*NUM_VAR+0] += dX[0];
      dQ[j*NUM_VAR+1] += dX[1];
      dQ[j*NUM_VAR+2] += dX[2];
      dQ[j*NUM_VAR+3] += dX[3];	  
    }

  // Compute the residual b-Ax. Start by copying the right hand side to a temp array. Then compute
  // A*x using the solution from this time step.
  for ( j=1; j <= nn; j++ )
    {
      r[0] = RHS[j*NUM_VAR+0];
      r[1] = RHS[j*NUM_VAR+1];
      r[2] = RHS[j*NUM_VAR+2];
      r[3] = RHS[j*NUM_VAR+3];

      d = iau[j];

      r[0] -= ( LHS[d*NUM_VAR*NUM_VAR+0*NUM_VAR+0]*dQ[j*NUM_VAR+0] +
		LHS[d*NUM_VAR*NUM_VAR+0*NUM_VAR+1]*dQ[j*NUM_VAR+1] +
		LHS[d*NUM_VAR*NUM_VAR+0*NUM_VAR+2]*dQ[j*NUM_VAR+2] +
		LHS[d*NUM_VAR*NUM_VAR+0*NUM_VAR+3]*dQ[j*NUM_VAR+3] );
      
      r[1] -= ( LHS[d*NUM_VAR*NUM_VAR+1*NUM_VAR+0]*dQ[j*NUM_VAR+0] +
		LHS[d*NUM_VAR*NUM_VAR+1*NUM_VAR+1]*dQ[j*NUM_VAR+1] +
		LHS[d*NUM_VAR*NUM_VAR+1*NUM_VAR+2]*dQ[j*NUM_VAR+2] +
		LHS[d*NUM_VAR*NUM_VAR+1*NUM_VAR+3]*dQ[j*NUM_VAR+3] );
      
      r[2] -= ( LHS[d*NUM_VAR*NUM_VAR+2*NUM_VAR+0]*dQ[j*NUM_VAR+0] +
		LHS[d*NUM_VAR*NUM_VAR+2*NUM_VAR+1]*dQ[j*NUM_VAR+1] +
		LHS[d*NUM_VAR*NUM_VAR+2*NUM_VAR+2]*dQ[j*NUM_VAR+2] +
		LHS[d*NUM_VAR*NUM_VAR+2*NUM_VAR+3]*dQ[j*NUM_VAR+3] );
      
      r[3] -= ( LHS[d*NUM_VAR*NUM_VAR+3*NUM_VAR+0]*dQ[j*NUM_VAR+0] +
		LHS[d*NUM_VAR*NUM_VAR+3*NUM_VAR+1]*dQ[j*NUM_VAR+1] +
		LHS[d*NUM_VAR*NUM_VAR+3*NUM_VAR+2]*dQ[j*NUM_VAR+2] +
		LHS[d*NUM_VAR*NUM_VAR+3*NUM_VAR+3]*dQ[j*NUM_VAR+3] );

      RMS += (r[0]*r[0] + r[1]*r[1] + r[2]*r[2] + r[3]*r[3]);
    }
  
  RMS /= (1.*nn);
  RMS = sqrt ( RMS );

  printf("  EXPLICIT UPDATE -> RMS=%.15e\n",RMS);

  freenull(rhs1);
  freenull(M);
  freenull(B);
  freenull(X);
  freenull(Y);
  freenull(dX);

  return RMS;

}


//==============================================================================
//
//      lu()
//
//      double *A             // The matrix to be factored. Will be overwritten.
//
//------------------------------------------------------------------------------

void lu( double *A )
{
  int i,j,k;
  int n = NUM_VAR;

  for ( k=0; k < n; k++ )                   // Column Loop
    {
      // Compute the multipliers
      for ( i=k+1; i < n; i++ )
        {
          A[i*n+k] = A[i*n+k] / A[k*n+k];
        }
      
      // Apply the mulitpliers to the submatrix
      for ( j = k+1; j < n; j++ )
        {
          for ( i=k+1; i < n; i++ )
            {
              A[i*n+j] = A[i*n+j] - A[i*n+k]*A[k*n+j];
            }
        }
    }
  return;
}


//==============================================================================
//
//      backstab()
//
//      double *A         // The factored matrix.
//      double *x         // The solution.
//      double *y         // The intermediate solution from the forwards 
//                        // substituion function.
//
//------------------------------------------------------------------------------

void backstab ( double *A, double *x, double *y )
{
  int n = NUM_VAR;

  x[3] = y[3]/A[3*n+3];
  x[2] = (y[2]-A[2*n+3]*x[3])/A[2*n+2];
  x[1] = (y[1]-A[1*n+2]*x[2]-A[1*n+3]*x[3])/A[1*n+1];
  x[0] = (y[0]-A[0*n+1]*x[1]-A[0*n+2]*x[2]-A[0*n+3]*x[3])/A[0*n+0];
  return;
}

//==============================================================================
//
//      forwardstab()
//
//      double *A            // The factored matrix.
//      double *y            // The intermediate solution.
//      double *b            // The right hand side of the system.
//
//------------------------------------------------------------------------------

void forwardstab ( double *A, double *y, double *b )
{
  int n = NUM_VAR;

  y[0] = b[0];
  y[1] = b[1]-A[1*n+0]*y[0];
  y[2] = b[2]-A[2*n+0]*y[0] - A[2*n+1]*y[1];
  y[3] = b[3]-A[3*n+0]*y[0]-A[3*n+1]*y[1]-A[3*n+2]*y[2];
  return;
}

//==============================================================================
//
//      MMMult_DtG()
//
//      Matrix matrix multiplication where the left matrix is diagonal and the
//      right matrix is general.
//
//***** This function assumes linear array memory layout and that the first matrix
//      is stored as a vector!
//
//      double *A            // The diagonal matrix (vector).
//      double *B            // The general matrix.
//      double *X            // The solution matrix.
//      int stride           // The size of the vector/matrix.
//
//------------------------------------------------------------------------------

void MMMult_DtG ( double *A, double *B, double *X, int stride )
{
  int i,j;

  for ( i=0; i < stride; i++ )
    {
      for ( j=0; j < stride; j++ )
	{
	  X[i*stride + j] = A[i] * B[i*stride+j];
	}
    }

  return;
}

//==============================================================================
//
//      MMMult_GtG()
//
//      Matrix matrix multiplication both general.
//
//***** This function assumes linear array memory layout.
//
//      double *A            // The first matrix.
//      double *B            // The second matrix.
//      double *X            // The solution matrix.
//      int stride           // The size of the vector/matrix.
//
//------------------------------------------------------------------------------

void MMMult_GtG ( double *A, double *B, double *X, int stride )
{
  int i,j,k;

  for ( i=0; i < stride; i++ )
    {
      for ( j=0; j < stride; j++ )
	{
	  X[i*stride + j] = 0.;

	  for ( k=0; k < stride; k++ )
	    {
	      X[i*stride+j] += A[i*stride+k]*B[k*stride+j];
	    }
	}
    }

  return;
}

//==============================================================================
//
//      MVMult()
//
//      Matrix vector multiplication, both general.
//
//***** This function assumes linear array memory layout.
//
//      double *A            // The matrix.
//      double *X            // The vector.
//      double *B            // The solution.
//      int stride           // The size of the matrix/vector.
//
//------------------------------------------------------------------------------

void MVMult ( double *A, double *X, double *B, int stride )
{
  int i,j;

  for ( i=0; i < stride; i++ )
    {
      B[i] = 0.;
      for ( j=0; j < stride; j++ )
	{
	  B[i] += A[i*stride+j] * X[j];
	}
    }

  return;
}



//==============================================================================
//
//      QR_House()
//
//      QR factorization using Householder matrices.
//
//***** This function assumes linear array memory layout.
//
//      double *A            // The matrix.
//      double *Q            // Q is m by m.
//      int m                // The size of the matrix.
//      int n                 
//
//------------------------------------------------------------------------------

void QR_House ( double *A, double *Q, int m, int n )
{
  int i,j,k;                             // Loop counters.
  double alpha, beta, gamma, s;          // Factors.
  double *v = NULL;
  double *acol = NULL;

  v = (double*)malloc(m*sizeof(double));
  if ( v == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'v' in QR.\n"); exit(1); }
  
  acol = (double*)malloc(m*sizeof(double));
  if ( acol == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'acol' in QR.\n"); exit(1); }
  
  // Set Q to be the identity matrix.
  for ( i=0; i < m; i++ )
    {
      for ( j=0; j < m; j++ )
	{
	  Q[i*m+j] = 0.;
	}
      Q[i*m+i] = 1.;
    }

  // Do the factorization.
  for ( k=0; k < n; k++ )
    {
      // Fill in v with the current column to start.
      for ( i=0; i < k; i++ )
	{
	  v[i] = 0.;
	}

      for ( i=k; i < m; i++ )
	{
	  v[i] = A[i*n+k];
	}

      s = sign(A[k*n+k]);
      alpha = -1.*sqrt(dotProductGeneral(v,v,m)) * s;

      v[k] -= alpha;

      beta = dotProductGeneral(v,v,m);

      if ( fabs(beta) < 10.E-14 )
	continue;

      for ( j=k; j < n; j++ )
        {
          for ( i=0; i < m; i++ )
            {
              acol[i] = A[i*n+j];
            }
          
          gamma = dotProductGeneral(v,acol,m);
          
          for ( i=0; i < m; i++ )
            {
              A[i*n+j] = acol[i] - (2.*gamma/beta)*v[i];
            } 
        }

      // Apply  H to q to build up Q.
      for ( i=0; i < m; i++ )
	{
	  // Copy column over to acol.
	  for ( j=0; j < m; j++ )
	    acol[j] = Q[j*m+i];
	  
	  gamma = dotProductGeneral(v,acol,m);

	  for ( j=0; j < m; j++ )
            {
              Q[j*m+i] = acol[j] - (2.*gamma/beta)*v[j];
            }
	}
    }

  freenull(v);
  freenull(acol);

  /*
  // Debug print.

  // Transpose Q using R.
  for ( i=0; i < m; i++ )
    {
      for ( j=0; j < m; j++ )
	{
	  R[i*m+j] = Q[j*m+i];
	}
    }

  for ( i=0; i < m; i++ )
    {
      for ( j=0; j < m; j++ )
	{
	  Q[i*m+j] = R[i*m+j];
	}
    }


  printf("Printing Q:\n");
  for ( i=0; i < m; i++ )
    {
      for ( j=0; j < m; j++ )
	{
	  printf("%lf  ",Q[i*m+j]);
	}
      printf("\n");
    }
  printf("\n");

  printf("Printing R:\n");
  for ( i=0; i < m; i++ )
    {
      for ( j=0; j < n; j++ )
	{
	  printf("%lf  ",A[i*n+j]);
	}
      printf("\n");
    }
  printf("\n");

  // Now multiply the result out to check the factorization.
  double QRA[24];
  double sum = 0;
  for ( i=0; i < m; i++ )
    {
      for ( j=0; j < n; j++ )
	{
	  sum = 0.;
	  for ( k=0; k < m; k++ )
	    sum += Q[i*m+k]*A[k*n+j];
	  QRA[i*n+j] = sum;
	}
    }

  printf("Printing QR=A:\n");
  for ( i=0; i < m; i++ )
    {
      for ( j=0; j < n; j++ )
	{
	  printf("%lf  ",QRA[i*n+j]);
	}
      printf("\n");
    }
  printf("\n");
  */

  return;
}


//==============================================================================
//
//      QR_QtB()
//
//      Uses the transpose of Q from the QR factorization to orthogonalize
//      the right hand side.
//
//
//      double *Q            // The orthogonal matrix (already transposed).
//      double *B            // The right hand side.
//      double *X            // The product.
//      int m                // The size of the matrix.
//
//------------------------------------------------------------------------------

void QR_QtB ( double *Q, double *B, double *X, int m )
{
  int i,j;                             // Loop counters.

  for ( i=0; i < m; i++ )
    {
      X[i] = 0.;
    }

  for ( i=0; i < m; i++ )
    {
      for ( j=0; j < m; j++ )
	{
	  X[i] += (Q[i*m + j] * B[j]);
	}
    }

  return;
}


//==============================================================================
//
//      QR_Backsub()
//
//      Uses the QR factorization to do the backwards substitution.
//      OVERWRITES B with the solution.
//
//      double *R            // The left hand side.
//      double *B            // The orthogonalized right hand side, and solution.
//      int n                // The size of the matrix.
//
//------------------------------------------------------------------------------

void QR_Backsub ( double *R, double *B, int n )
{
  int i,j;                             // Loop counters.

  for ( i=n-1; i >= 0; i-- )
    {
      for ( j=i+1; j < n; j++ )
        {
          B[i] -= B[j]*R[i*n+j];
        }
      B[i] = B[i] / R[i*n+i];
    }

  return;
}



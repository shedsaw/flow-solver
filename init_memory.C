//=============================================================
// 
//  init_memory.C
//  
//  Functions to build the linear system.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "sys_mem.h"
#include "params.h"
#include "Defined_Variables.h"
#include "cv_calc.h"
#include "init_memory.h"


//=============================================================
// 
//  Build_CRS()
//
//  Constructs the information for the CRS arrays used in the
//  solver. This is inefficient (and redundant) but I'm going
//  ahead this way for now.
//  
//  GRID *grid                        // The grid.
//  PARAMS p                          // The parameters.
//  SYS_MEM *smem                     // The memory object.
//
//=============================================================

void Build_CRS ( GRID *grid, PARAMS p, SYS_MEM *smem )
{
  int i,j,n;                 // Loop counters.
  int sum;                   // Sum of entries.
  int ptr = 0;               // Array pointer.
  int idx1,idx2;             // Array index pointers.

  // Pointers.
  int nn = grid->nn;
  int *nsn = grid->nsn;
  int *nnsn = grid->nnsn;

  // Calculate the needed memory.
  sum = 0;
  for ( i=1; i <= nn; i++ )
    {
      sum += (nnsn[i+1] - nnsn[i] + 1 );
    }

  // This should be nnsn[nn+1] + nn!
  if ( (nnsn[nn+1]+nn) != sum )
    {
      printf("FATAL ERROR: In Build_CRS a discrepancy was found:\n");
      printf("     nnodes = %d\n",nn);
      printf("     nnsn[nn+1] = %d\n",nnsn[nn+1]);
      printf("     sum = %d\n",sum);
      exit(1);
    }

  // Allocate our memory.
  if ( smem->ia == NULL )
    smem->ia = (int*)malloc( (nn+2)*sizeof(int) );
  if ( smem->ia == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'ia'.\n"); exit(1); }

  if ( smem->iau == NULL )
    smem->iau = (int*)malloc( (nn+1)*sizeof(int) );
  if ( smem->iau == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'iau'.\n"); exit(1); }

  if ( smem->ja == NULL )
    smem->ja = (int*)malloc( sum*sizeof(int) );
  if ( smem->ja == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'ja'.\n"); exit(1); }

  // Fill in ia. Basically copy in nsn but add the node itself.
  smem->ia[0]=0;
  for ( n=1; n <= (nn+1); n++ )
    {
      smem->ia[n] = nnsn[n] + (n-1);
    }

  // Fill in iau. This will be the last entry in ia for that node.
  smem->iau[0]=0;
  for ( n=1; n <= nn; n++ )
    {
      smem->iau[n] = smem->ia[n+1] - 1;
    }

  // Fill in the ja array.
  for ( n=1; n <= nn; n++ )
    {
      idx1 = nnsn[n];
      idx2 = nnsn[n+1];

      for ( j=idx1; j < idx2; j++ )
	{
	  smem->ja[ptr] = nsn[j];
	  ptr++;
	}
      smem->ja[ptr] = n;
      ptr++;
    }

  if ( DEBUG )
    {
      printf("IA DEBUG:\n");
      for ( i=0; i <= (nn+1); i++ )
	{
	  printf("%d\n",smem->ia[i]);
	}

      printf("IAU DEBUG\n");
      for ( i=0; i <= nn; i++ )
	{
	  printf("%d\n",smem->iau[i]);
	}

      printf("JA DEBUG\n");
      for ( i=1; i <= nn; i++ )
	{
	  idx1 = smem->ia[i];
	  idx2 = smem->ia[i+1];

	  printf("%d: ",i);
	  for ( j=idx1; j < idx2; j++ )
	    {
	      printf("%d ",smem->ja[j]);
	    }
	  printf("\n");
	}
    }
  

  return;
}


//=============================================================
// 
//  build_matrix()
//
//  Constructs the left hand side matrix using CRS.
//  V/dt will appear on the diagonal and the rest of the matrix
//  will be filled in with the Jacobian if needed.
//
//  I believe that I'm essentially putting V/dtau on the diagonal
//  to increase convergence of my linear system. For unsteady flows,
//  dt (real time) is going on the RHS.
//
//  GRID *grid                        // The grid.
//  PARAMS p                          // The parameters.
//  SYS_MEM *smem                     // The memory object.  
//
//=============================================================

void build_matrix ( GRID *grid, PARAMS p, SYS_MEM *smem )
{
  int i,j,k;                         // Loop counters.
  int mdim;                          // Matrix dimension.
  double V;                          // Area of control volume.
  double cfl;                        // Current CFL to apply to iteration.
  double low_time = 1.0E10;          // The smallest time step based on stability.
  FILE *fp = NULL;

  // Pointers.
  int nn = grid->nn;
  int *ia = smem->ia;
  int *iau = smem->iau;
  
  mdim = (ia[nn+1])*NUM_VAR*NUM_VAR;

  // Allocate our memory.
  if ( smem->LHS == NULL )
    {
      smem->LHS = (double*)malloc( mdim*sizeof(double) );
      if ( smem->LHS == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'LHS'.\n"); exit(1); }
    }

  if ( smem->RHS == NULL )
    {
      smem->RHS = (double*)malloc( ((nn+1)*NUM_VAR)*sizeof(double) );
      if ( smem->RHS == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'RHS'.\n"); exit(1); }
    }

  // Clean out the memory.
  for ( i=0; i < mdim; i++ )
    {
      smem->LHS[i] = 0.;
    }

  //fp = fopen("build_matrix.dat","w");

  // Get the current eigenvalue contributions.
  cv_calc_eigenvalue_contribution( grid, p );

  // Compute the CFL number for the current iteration.
  if ( p.method == 1 )
    {
      if ( grid->citer < p.ramp )
	{
	  cfl = 1. + (1.*(grid->citer-1))/(1.*p.ramp)*(p.CFL - 1.);
	}
      else
	cfl = p.CFL;
    }
  else if ( p.method == 0 )
    {
      cfl = 1.0;
    }
  else
    {
      cfl = p.CFL;
    }
  
  // If we restart, I'm hoping to use the above logic to pick up where we left off.
  // If we restart a solution, use the given CFL number.
  //if ( p.isrestart == 1 )
  //  cfl = p.CFL;
  
  // Set the value of the cfl for the grid object.
  grid->cfl = cfl;

  // Fill in the matrix with the diagonal term ( proportional to V/dt).
  for ( j=1; j <= nn; j++ )
    {
      i = iau[j];
      V = grid->cv_area[j];

      // Currently, this is deprecated since I switched over to another formulation for the time step
      // given in the summer course notes.
      
      //rho = grid->Q[j*NUM_VAR+0];
      //u =  (grid->Q[j*NUM_VAR+1])/rho;
      //v =  (grid->Q[j*NUM_VAR+2])/rho;
      //E =  (grid->Q[j*NUM_VAR+3]);

      //P = (gamma - 1.)*E - ((gamma-1.)/2.)*rho*(u*u + v*v);

      //c = sqrt( gamma*P/rho );
      //if ( P <= 0. || rho <= 0. )
      //{
      //  printf("Either PRESSURE or DENSITY is negative in build_matrix for control volume %d\n",j);
      //  exit(0);
      //	}
      //====================================================================================================

      //grid->dt[j] = ( cfl * grid->length_scale[j] )/( sqrt(u*u + v*v) + c );  OLD METHOD

      grid->dt[j] = cfl * ( V / grid->ev_con[j] );

      if ( grid->dt[j] < low_time )
	low_time = grid->dt[j];

      if ( grid->dt[j] <= 0. )
	{
	  printf("  Negative time step for cv %d!\n",j);
	  exit(0);
	}

      if ( !isfinite( grid->dt[j] ) )
	{
	  printf("  FATAL ERROR: time step is not finite fo cv %d.\n",j);
	  fflush(stdout);
	  exit(0);
	}

      // HACK for unsteady residual.
      //if ( p.unsteady == 1 )
      //grid->dt[j] = p.mintime;

      for ( k=0; k < NUM_VAR; k++ )
	{
	  if ( p.tacc == 2 && grid->citer > (grid->rfiter+1) )
	    {
	      //smem->LHS[ i*NUM_VAR*NUM_VAR + k*NUM_VAR + k] = V/(grid->dt[j]) + 1.5 * V/p.mintime;
	      smem->LHS[ i*NUM_VAR*NUM_VAR + k*NUM_VAR + k] = 1.5 * V/p.mintime;
	    }
	  else if ( p.unsteady > 0 )
	    {
	      //smem->LHS[ i*NUM_VAR*NUM_VAR + k*NUM_VAR + k] = V/(grid->dt[j]) + V/p.mintime;
	      smem->LHS[ i*NUM_VAR*NUM_VAR + k*NUM_VAR + k] = V/p.mintime;
	    }
	  else
	    {
	      smem->LHS[ i*NUM_VAR*NUM_VAR + k*NUM_VAR + k] = V/(grid->dt[j]);
	    }
	}

      //fprintf(fp,"Node %d: \n",j-1);
      //fprintf(fp,"  Q[0] = %.15E\n",grid->Q[j*NUM_VAR+0]);
      //fprintf(fp,"  Q[1] = %.15E\n",grid->Q[j*NUM_VAR+1]);
      //fprintf(fp,"  Q[2] = %.15E\n",grid->Q[j*NUM_VAR+2]);
      //fprintf(fp,"  Q[3] = %.15E\n",grid->Q[j*NUM_VAR+3]);
      //fprintf(fp,"  V = %.15E\n",V);
      //fprintf(fp,"  p = %.15E\n",P);
      //fprintf(fp,"  c = %.15E\n",c);
      //fprintf(fp,"  length_scale = %.15E\n",grid->length_scale[j]);
      //fprintf(fp,"  dt = %.15E\n\n",dt);
      //fprintf(fp,"%.15E\n",V/dt);
    }
  //fclose(fp);

  printf("The minimum time step for stability with CFL = %.15e  is %.15e\n",cfl,low_time);
  fflush(stdout);

  return;
}

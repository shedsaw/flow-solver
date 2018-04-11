//=============================================================
// 
//  solve.C
//  
//  Functions to solve the equations for the current time step.
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
#include "residual.h"
#include "init_memory.h"
#include "gradient.h"
#include "hessian.h"
#include "compute_derivatives.h"
#include "mesh_utils.h"
#include "linear_algebra.h"
#include "reconstruct.h"
#include "curved_boundaries.h"
#include "limiter.h"
#include "solve.h"
#include "grid_io.h"
#include "cv_calc.h"
#include "maps.h"

//=============================================================
// 
//  solve()
//
//  GRID *grid;                          // The grid.
//  SYS_MEM *smem;                       // The memory for the system.
//  PARAMS p;                            // The parameters.
//  int ts;                              // The current time step.
//
//=============================================================

double solve ( GRID *grid, SYS_MEM *smem, PARAMS p, int ts )
{
  int i,j,k,m,n;                       // Loop counters.
  double rms,RMS;                      // Rms error.
  double url2;                         // Unsteady residual L2 norm.
  double dql2;                         // L2 norm for change in Q.
  int jacupdate;                       // Jacobian update flag.
  char buff[100];
  int DO_LIMITER = 0;

  dql2 = 0.;

  // build the layer of neighbors needed for the iteration level.
  if ( p.order == 3 && ts == (p.foits+p.soits + 1) )
    {
      printf("Adding first layer of neighbors........\n");
      for ( i=1; i <= grid->nn; i++ )
	{
	  Add_Layer_to_GNSN ( grid->nn, grid->num_elem, grid->c2n, grid->esn, grid->nnsn, grid->edges,
			      &(grid->gnsn), &(grid->gnnsn), grid->gnnsn1, grid->node_fill_level, i );
	}
      printf("Done\n");
      
      printf("Adding second layer of neighbors........\n");
      for ( i=1; i <= grid->nn; i++ )
	{
	  if ( grid->node_order[i] > 2 )
	    {
	      Add_Layer_to_GNSN ( grid->nn, grid->num_elem, grid->c2n, grid->esn, grid->nnsn, grid->edges,
				  &(grid->gnsn), &(grid->gnnsn), grid->gnnsn1, grid->node_fill_level, i );
	    }
	}
      printf("Done\n");

      printf("Checking to make sure the nodes all have a large enough stencil.......\n");

      for ( i=1; i <= grid->nn; i++ )
	{
	  if ( grid->node_order[i] < 3 )
	    continue;
	  
	  // Find the number of neighbors.
	  j = grid->gnnsn[ grid->gnnsn1[i+1] ] - grid->gnnsn[ ( grid->gnnsn1[i] + 1 ) ];
	  
	  if ( j < 6 )
	    {
	      for ( k=0; k < 3; k++ )
		{
		  printf("  Node %d has %d neighbors in its stencil. Adding degree %d neighbors to its stencil.......\n",i,j,
			 grid->node_fill_level[i]+1);
		  Add_Layer_to_GNSN ( grid->nn, grid->num_elem, grid->c2n, grid->esn, grid->nnsn, grid->edges,
				      &(grid->gnsn), &(grid->gnnsn), grid->gnnsn1, grid->node_fill_level, i );
	   
		  // recompute the number of neighbors.
		  j = grid->gnnsn[ grid->gnnsn1[i+1] ] - grid->gnnsn[ ( grid->gnnsn1[i] + 1 ) ];

		  if ( j >= 6 ) break;
		}

	      if ( k >= 3 )
		{
		  printf("FATAL ERROR: Node %d has attempted to add too many layers of neighbors. Halting...\n",i);
		  fflush(stdout);
		  exit(0);
		}

	      printf("Finished. Node %d now has %d neighbors in its stencil.\n",i,j);
	    }
	      
	}
    }
  else if ( p.order == 4 && ts == (p.foits+p.soits + 1) )
    {
      printf("Adding first layer of neighbors........\n");
      for ( i=1; i <= grid->nn; i++ )
	{
	  Add_Layer_to_GNSN ( grid->nn, grid->num_elem, grid->c2n, grid->esn, grid->nnsn, grid->edges,
			      &(grid->gnsn), &(grid->gnnsn), grid->gnnsn1, grid->node_fill_level, i );
	}
      printf("Done\n");

      printf("Adding second layer of neighbors........\n");
      for ( i=1; i <= grid->nn; i++ )
	{
	  if ( grid->node_order[i] > 2 )
	    {
	      Add_Layer_to_GNSN ( grid->nn, grid->num_elem, grid->c2n, grid->esn, grid->nnsn, grid->edges,
				  &(grid->gnsn), &(grid->gnnsn), grid->gnnsn1, grid->node_fill_level, i );
	    }
	}
      printf("Done\n");

      printf("Checking to make sure the nodes all have a large enough stencil.......\n");

      for ( i=1; i <= grid->nn; i++ )
	{
	  if ( grid->node_order[i] < 3 )
	    continue;
	  
	  // Find the number of neighbors.
	  j = grid->gnnsn[ grid->gnnsn1[i+1] ] - grid->gnnsn[ ( grid->gnnsn1[i] + 1 ) ];
	  
	  if ( j < 6 )
	    {
	      for ( k=0; k < 3; k++ )
		{
		  printf("  Node %d has %d neighbors in its stencil. Adding degree %d neighbors to its stencil.......\n",i,j,
			 grid->node_fill_level[i]+1);
		  Add_Layer_to_GNSN ( grid->nn, grid->num_elem, grid->c2n, grid->esn, grid->nnsn, grid->edges,
				      &(grid->gnsn), &(grid->gnnsn), grid->gnnsn1, grid->node_fill_level, i );
	   
		  // recompute the number of neighbors.
		  j = grid->gnnsn[ grid->gnnsn1[i+1] ] - grid->gnnsn[ ( grid->gnnsn1[i] + 1 ) ];

		  if ( j >= 6 ) break;
		}

	      if ( k >= 3 )
		{
		  printf("FATAL ERROR: Node %d has attempted to add too many layers of neighbors. Halting...\n",i);
		  fflush(stdout);
		  exit(0);
		}

	      printf("Finished. Node %d now has %d neighbors in its stencil.\n",i,j);
	    }
	      
	}
    }
  else if ( p.order == 4 && ts == (p.foits+p.soits+p.toits + 1) )
    {
      printf("Adding third layer of neighbors........\n");
      for ( i=1; i <= grid->nn; i++ )
	{
	  if ( grid->node_order[i] > 2 )
	    {
	      if ( grid->node_fill_level[i] == 2 )
		{
		  Add_Layer_to_GNSN ( grid->nn, grid->num_elem, grid->c2n, grid->esn, grid->nnsn, grid->edges,
				      &(grid->gnsn), &(grid->gnnsn), grid->gnnsn1, grid->node_fill_level, i );
		}
	    }
	}
      printf("Done\n");

      printf("Checking to make sure the nodes all have a large enough stencil.......\n");

      for ( i=1; i <= grid->nn; i++ )
	{
	  if ( grid->node_order[i] < 3 )
	    continue;
	  
	  // Find the number of neighbors.
	  j = grid->gnnsn[ grid->gnnsn1[i+1] ] - grid->gnnsn[ ( grid->gnnsn1[i] + 1 ) ];
	  
	  if ( j < 10 )
	    {
	      for ( k=0; k < 3; k++ )
		{
		  printf("  Node %d has %d neighbors in its stencil. Adding degree %d neighbors to its stencil.......\n",i,j,
			 grid->node_fill_level[i]+1);
		  Add_Layer_to_GNSN ( grid->nn, grid->num_elem, grid->c2n, grid->esn, grid->nnsn, grid->edges,
				      &(grid->gnsn), &(grid->gnnsn), grid->gnnsn1, grid->node_fill_level, i );
	   
		  // recompute the number of neighbors.
		  j = grid->gnnsn[ grid->gnnsn1[i+1] ] - grid->gnnsn[ ( grid->gnnsn1[i] + 1 ) ];

		  if ( j >= 10 ) break;
		}

	      if ( k >= 3 )
		{
		  printf("FATAL ERROR: Node %d has attempted to add too many layers of neighbors. Halting...\n",i);
		  fflush(stdout);
		  exit(0);
		}

	      printf("Finished. Node %d now has %d neighbors in its stencil.\n",i,j);
	    }
	      
	}
    }

  // DEBUG Code to test my flux scheme.
  if ( 0 )
    {
      for ( i=1; i <= (grid->nn + grid->nn_ghost); i++ )
	{
	  grid->Q[i*NUM_VAR+0] = 1.;
	  grid->Q[i*NUM_VAR+1] = 0.5 * grid->Q[i*NUM_VAR+0];
	  grid->Q[i*NUM_VAR+2] = 0.;
	  grid->Q[i*NUM_VAR+3] = 1.0 / ( p.gamma*(p.gamma-1.0) ) + 0.5 * ( (grid->Q[i*NUM_VAR+1]*grid->Q[i*NUM_VAR+1] + grid->Q[i*NUM_VAR+2]*grid->Q[i*NUM_VAR+2]) / grid->Q[i*NUM_VAR+0] );
	}
    }

  if ( grid->phi == NULL )
    {
      grid->phi = (double*)malloc( (grid->nn+1)*2*4*sizeof(double) );
      if ( grid->phi == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'phi'.\n"); exit(1); }
    }

  // Let's go ahead and allocate memory for the gradient and hessian and set them to zero.
  if ( grid->grad == NULL )
    {
      grid->grad = (double*)malloc((grid->nn + 1)*NDIM*NUM_VAR*sizeof(double));
      if ( grid->grad == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'grad'.\n"); exit(1); }
    }

  if ( grid->hess == NULL )
    {
      grid->hess = (double*)malloc((grid->nn + 1)*NUM_MOM*NUM_VAR*sizeof(double));
      if ( grid->hess == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'hess'.\n"); exit(1); }
    }

  if ( ts == 1 )
    {
      if ( grid->Qe == NULL )
	{
	  grid->Qe = (double*)malloc((grid->nn+grid->nn_ghost+1)*NUM_VAR*sizeof(double));

	  if ( grid->Qe ==NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'Qe'.\n"); exit(1); }
	}

      if ( grid->Source == NULL )
	{
	  grid->Source = (double*)malloc((grid->nn+grid->nn_ghost+1)*NUM_VAR*sizeof(double));

	  if ( grid->Source ==NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'Source'.\n"); exit(1); }
	}

      for ( n=1; n <= (grid->nn+grid->nn_ghost); n++ )
	{
	  grid->Qe[n*NUM_VAR+0] = 0.0;
	  grid->Qe[n*NUM_VAR+1] = 0.0;
	  grid->Qe[n*NUM_VAR+2] = 0.0;
	  grid->Qe[n*NUM_VAR+3] = 0.0;

	  grid->Source[n*NUM_VAR+0] = 0.0;
	  grid->Source[n*NUM_VAR+1] = 0.0;
	  grid->Source[n*NUM_VAR+2] = 0.0;
	  grid->Source[n*NUM_VAR+3] = 0.0;
	}
    }

  //if ( (MMS && ts == 1) && 1 ) // changed from not recon_prim
  if ( ts == 1 && ( MMS || ANNULUS || MMS_EXP ) )
    {
      MMS_Compute_Q_Exact(grid);

      sprintf(buff,"tecplot_qexact.dat");
      write_tecplot_solutionE ( buff, grid );
    }

  // Clean them out.
  for ( n=1; n <= grid->nn; n++ )
    {
      for ( i=0; i < NUM_VAR; i++ )
	{
	  for ( j=0; j < NDIM; j++ )
	    {
	      grid->grad[n*NDIM*NUM_VAR + i*NDIM + j] = 0.;
	    }
	}
    }
  
  for ( i=0; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  for ( k=0; k < NUM_MOM; k++ )
	    grid->hess[i*NUM_MOM*NUM_VAR + j*NUM_MOM + k] = 0.; 
	}
    }
  
  // We need to copy the previous solution over to the old solution vector.
  for ( i=1; i <= (grid->nn + grid->nn_ghost); i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  grid->nQ[i*NUM_VAR+j] = grid->Q[i*NUM_VAR+j];
	}
    }

  // Begin our newton iterations.
  for ( m=1; m <= p.iter; m++ )
    {
      if ( ts == p.foits )
	{
	  url2 = 0.;
	}

      // Determine whether or not to update the jacobians.
      if ( (m-1)%p.jufreq == 0 ) {  jacupdate = 1; }
      else { jacupdate = 0; }

      if ( m == 1 ) { jacupdate = 1; } // Always update on the first newton iteration.

      if ( p.method == 0 ) { jacupdate = 0; } // If explicit, never compute the jacobians.

      // Set the limiters to 1.
      for ( i=1; i <= grid->nn; i++ )
	{
	  for ( j=0; j < 8; j++ )
	    grid->phi[i*8+j] = 1.;
	}

      // ======================================================================================
      // If doing primitive variable reconstruction, I am assuming that the Q are in conserved
      // form at the beginning of the Newton iteration loop.
      // =====================================================================================
      // ===================================================================================
      // Now, if I am doing reconstruction on the primitive variables, I need to convert my
      // conserved variables to primitives before computing the gradient/hessian on them.
      // ==================================================================================
      if ( RECON_PRIM && p.ic != 3 )  // ic_3 uses primitive variables.
	{
	  for ( i=1; i <= (grid->nn + grid->nn_ghost); i++ )                     // ghost values are initially set in the ic initialization.
	    {
	      ConvertConservedToPrimitive ( p.gamma , &(grid->nQ[i*NUM_VAR]) );  // I only need to convert nQ, my current guess. Every function
	                                                                         // in the Newton loop uses nQ.
	    }
	}
      
      // Clean out the left hand side and put on V/dtau, this will also call cv_calc_eigenvalue_contribution() -> this functions assumes that the variables are in the proper
      // form (either conserved or primitive depending on the solve).
      build_matrix( grid, p, smem );

      // Check the diagonal for NaN's.
      for ( i=1; i <= grid->nn; i++ )
	{
	  j = smem->iau[i];
	  for ( k=0; k < NUM_VAR; k++ )
	    {
	      if ( !(isfinite( smem->LHS[j*NUM_VAR*NUM_VAR + k*NUM_VAR + k ])) )
		{
		  printf("Problem detected in the LHS matrix diagonal after build_matrix(). Node %d, Variable %d\n",i,k);
		}
	    }
	}

      if ( (p.order == 1 || p.order == 2) && ts > p.foits )
	Compute_Gradient_LeastSquares(grid);
      else if ( p.order == 3 )
	{
	  if ( ts <= p.foits )
	    {
	      // Do nothing.
	      ;
	    }
	  else if ( ts > p.foits && ts <= (p.foits+p.soits) )
	    Compute_Gradient_LeastSquares(grid);
	  else
	    {
	      if ( CYLINDER && CYLINDER_HACK )
		Compute_Gradient_LeastSquares(grid);

	      //Compute_Hessian(grid,p);
	      Compute_Derivatives(grid,p,3);
	    }
	}
      else if ( p.order == 4 )
	{
	  if ( ts <= p.foits )
	    {
	      // Do nothing.
	      ;
	    }
	  else if ( ts > p.foits && ts <= (p.foits+p.soits) )
	    Compute_Gradient_LeastSquares(grid);
	  else if ( ts > (p.foits+p.soits) && ts <= (p.foits+p.soits+p.toits) )
	    {
	      if ( CYLINDER && CYLINDER_HACK )
		Compute_Gradient_LeastSquares(grid);

	      //Compute_Hessian(grid,p);
	      Compute_Derivatives(grid,p,3);
	    }
	  else
	    {
	      if ( CYLINDER && CYLINDER_HACK )
		Compute_Gradient_LeastSquares(grid);

	      Compute_Derivatives(grid,p,4);
	    }
	}

      // NEED TO FIX LIMITER OR USE PRIMITIVES! HOLD OFF FOR NOW, LIMITER IS NOT USED NOW
      // ********** LIMITER FUNCTIONALITY IS FIXED!
      
      // What I think I want to try here is to apply the limiter only during the phase of the solve that reconstructs to
      // the desired accuracy as set in the parameter file.

      if ( p.order == 0 )
	{
	  DO_LIMITER = 0;
	}
      else if ( p.order == 1 || p.order == 2 )
	{
	  if ( p.limit > 0 && ts > p.foits )
	    DO_LIMITER = 1;
	}
      else if ( p.order == 3 )
	{
	  if ( p.limit > 0 && ts > (p.foits+p.soits) )
	    DO_LIMITER = 1;
	}
      else if ( p.order == 4 )
	{
	  if ( p.limit > 0 && ts > (p.foits + p.soits + p.toits) )
	    DO_LIMITER = 1;
	}

      // Get the new values for the limiters if desired.
      if ( DO_LIMITER )
	{
	  if ( p.limit == 1 )
	    compute_limiter_barth ( grid, p );

	  else if ( p.limit == 2 )
	    compute_limiter_venkatakrishnan ( grid, p );

	  else if ( p.limit == 3 )
	    compute_limiter_venkatakrishnan2 ( grid, p );

	  else if ( p.limit == 4 )
	    compute_limiter_gooch1 ( grid, p );

	  else if ( p.limit == 5 )
	    compute_limiter_gooch2 ( grid, p );
	}

      // Debug the limiters, looking for visual confirmation.
      if ( ts > (p.foits+p.soits) && 0 )
	{
	  compute_limiter_barth ( grid, p );
	  sprintf( buff, "barth_limiter.dat" );
	  write_tecplot_cv_limiter ( buff, grid );

	  compute_limiter_venkatakrishnan ( grid, p );
	  sprintf( buff, "venkatakrishnan_limiter.dat" );
	  write_tecplot_cv_limiter ( buff, grid );

	  compute_limiter_venkatakrishnan2 ( grid, p );
	  sprintf( buff, "venkatakrishnan_limiter2.dat" );
	  write_tecplot_cv_limiter ( buff, grid );

	  compute_limiter_gooch1 ( grid, p );
	  sprintf( buff, "gooch_limiter1.dat" );
	  write_tecplot_cv_limiter ( buff, grid );

	  compute_limiter_gooch2 ( grid, p );
	  sprintf( buff, "gooch_limiter2.dat" );
	  write_tecplot_cv_limiter ( buff, grid );
	}

      // Write out the limiter for debugging.
      //write_tecplot_cv_limiter ( "grid_limiter.dat", grid );

      // Compute the extrema in the field generated by the reconstruction operator.
      //Compute_Extrema ( grid, p, m );





      // MMS code that will bypass the rest of the code below it.
      if ( MMS || MMS_EXP ) // changed from not recon_prim
	{
	  if ( (ts == 1 && MMS) && ANNULUS )  // Major Hack reporting for duty!
	    {
	      grid->bbc[0] = 0;
	      grid->bbc[1] = 0;
	      grid->bbc[2] = 14;
	      grid->bbc[3] = 0;
	      grid->bbc[4] = 16;
	    }

	  // get the new source terms.
	  if ( ts == 1 && MMS )
	    MMS_Compute_Source_Terms(grid);
	  else if ( ts == 1 && MMS_EXP )
	    MMS_Compute_Source_Terms_Exp(grid);

	  if ( (ts == 1 && MMS) && ANNULUS )  // Major Hack reporting for duty!
	    {
	      grid->bbc[0] = 0;
	      grid->bbc[1] = 0;
	      grid->bbc[2] = 0;
	      grid->bbc[3] = 0;
	      grid->bbc[4] = 0;
	    }
	  
	  build_interior_residuals ( grid, smem, p, jacupdate, ts );

	  printf("Interior Residual Check.\n");
	  for ( i=1; i <= grid->nn; i++ )
	    {
	      for ( j=0; j < NUM_VAR; j++ )
		{
		  if ( !(isfinite(grid->R[i*NUM_VAR+j])) )
		    {
		      if ( i <= grid->nn )
			{
			  printf("  Found a problem at node %d in variable %d: grid->R[i*NUM_VAR+j] = %f\n",i,j,(float)grid->R[i*NUM_VAR+j]);
			  fflush(stdout);
			}
		      else
			{
			  printf("  Found a problem at ghost node %d in variable %d: grid->R[i*NUM_VAR+j] = %f\n",i,j,(float)grid->R[i*NUM_VAR+j]);
			  fflush(stdout);
			}
		    }
		}
	    }
	  
	  /*
	  // No boundary conditions. Freeze the boundary nodes at the correct value.
	  for ( i=1; i <= grid->nbedges; i++ )
	    {
	      // Grab the data from the structure.
	      n =  grid->bedges[i*5+0];

	      // Copy the variables over.
	      grid->R[n*NUM_VAR+0] = 0.;
	      grid->R[n*NUM_VAR+1] = 0.;
	      grid->R[n*NUM_VAR+2] = 0.;
	      grid->R[n*NUM_VAR+3] = 0.;
	    }
	  */

	  build_boundary_residuals ( grid, smem, p, jacupdate ); 


	  printf("Boundary Residual Check.\n");
	  for ( i=1; i <= grid->nn; i++ )
	    {
	      for ( j=0; j < NUM_VAR; j++ )
		{
		  if ( !(isfinite(grid->R[i*NUM_VAR+j])) )
		    {
		      if ( i <= grid->nn )
			{
			  printf("  Found a problem at node %d in variable %d: grid->R[i*NUM_VAR+j] = %f\n",i,j,(float)grid->R[i*NUM_VAR+j]);
			  fflush(stdout);
			}
		      else
			{
			  printf("  Found a problem at ghost node %d in variable %d: grid->R[i*NUM_VAR+j] = %f\n",i,j,(float)grid->R[i*NUM_VAR+j]);
			  fflush(stdout);
			}
		    }
		}
	    }
	  
	  // Set our initial guess to zero.
	  for ( i=0; i <= grid->nn; i++ )
	    {
	      for ( j=0; j < NUM_VAR; j++ )
		{
		  grid->dQ[i*NUM_VAR+j] = 0.;
		}
	    }
	  
	  // Build the right hand side.
	  for ( i=1; i <= grid->nn; i++ )
	    {
	      for ( j=0; j < NUM_VAR; j++ )
		{
		  smem->RHS[i*NUM_VAR+j] = -( grid->R[i*NUM_VAR+j] ) + grid->Source[i*NUM_VAR+j] ;
		}
	    }
	  
	  // Kill updates on the boundary.
	  /*
	  for ( i=1; i <= grid->nbedges; i++ )
	    {
	      // Grab the data from the structure.
	      n =  grid->bedges[i*5+0];

	      for ( j=0; j < NUM_VAR; j++ )
		smem->RHS[i*NUM_VAR+j] = 0.;
	    }
	  */

	  // Get the Unsteady residual L2 norm.
	  url2 = 0.;
	  for ( i=1; i <= grid->nn; i++ )
	    {
	      for ( j=0; j < NUM_VAR; j++ )
		{
		  url2 += ( (smem->RHS[i*NUM_VAR+j]) * (smem->RHS[i*NUM_VAR+j]) );
		}
	    }
	  
	  // Generate our solution update vector.
	  RMS = linear_solve( grid, smem, p );
	  
	  // Get the L2 norm of dQ.
	  dql2 = 0.;
	  for ( i=1; i <= grid->nn; i++ )
	    {
	      for ( j=0; j < NUM_VAR; j++ )
		{
		  dql2 += ( grid->dQ[i*NUM_VAR+j] * grid->dQ[i*NUM_VAR+j] );
		}
	    }

	  printf("  Newton Iteration %d -> Linear System RMS error = %.15e\n",m,RMS);
	  printf("                         Residual L2 Norm = %.15e    CFL = %f\n",rms,(float)grid->cfl);
	  printf("                         Unsteady Residual L2 Norm = %.15E\n",sqrt(url2));
	  printf("                         L2 Norm of dQ = %.15E\n",sqrt(dql2));
	  fflush(stdout);
	  
	  // Update the solution vector.
	  
	  
	  // dQ is the update to the conserved variables.
	  if ( RECON_PRIM && p.ic != 3 )  // ic_3 uses primitive variables.
	    {
	      for ( i=1; i <= (grid->nn + grid->nn_ghost); i++ )    
		{
		  ConvertPrimitiveToConserved ( p.gamma , &(grid->nQ[i*NUM_VAR]) );
		}
	    }
	  
	  //for ( i=0; i <= (grid->nn + grid->nn_ghost); i++ )
	  for ( i=1; i <= grid->nn; i++ )
	    {
	      for ( j=0; j < NUM_VAR; j++ )
		{
		  grid->nQ[i*NUM_VAR+j] += (grid->dQ[i*NUM_VAR+j]);		       
		}
	    }

	  continue;
	}






      // Get the resiudals. Residual code for the most part does not depend on the state of the variables. Its
      // calling functions do however (which have been adjusted).

      build_interior_residuals ( grid, smem, p, jacupdate, ts );

      if ( 0 )
	{
	  printf("After calling build interior residuals.\n");

	  for ( i=1; i <= grid->nn; i++ )
	    {
	      printf("Node %d:\n",i);
	      printf("  Var 1: %.15e\n",grid->R[i*NUM_VAR+0]);
	      printf("  Var 2: %.15e\n",grid->R[i*NUM_VAR+1]);
	      printf("  Var 3: %.15e\n",grid->R[i*NUM_VAR+2]);
	      printf("  Var 4: %.15e\n",grid->R[i*NUM_VAR+3]);
	    }
	}

      printf("Interior Residual Check.\n");
      //for ( i=1; i <= (grid->nn + grid->nn_ghost); i++ )
      for ( i=1; i <= grid->nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      if ( !(isfinite(grid->R[i*NUM_VAR+j])) )
		{
		  if ( i <= grid->nn )
		    {
		      printf("  Found a problem at node %d in variable %d: grid->R[i*NUM_VAR+j] = %f\n",i,j,(float)grid->R[i*NUM_VAR+j]);
		      fflush(stdout);
		    }
		  else
		    {
		      printf("  Found a problem at ghost node %d in variable %d: grid->R[i*NUM_VAR+j] = %f\n",i,j,(float)grid->R[i*NUM_VAR+j]);
		      fflush(stdout);
		    }
		}
	    }
	}
      
      
      // Pick your poison.
      build_boundary_residuals ( grid, smem, p, jacupdate );    //  :( Doesn't work yet. Ok, so it works for approximate Jacobians but not numerical ones.

      //fix_farfield_boundary ( grid, smem, p );

      if ( 0 )
	{
	  printf("After calling build boundary residuals.\n");

	  for ( i=1; i <= grid->nn; i++ )
	    {
	      printf("Node %d:\n",i);
	      printf("  Var 1: %.15e\n",grid->R[i*NUM_VAR+0]);
	      printf("  Var 2: %.15e\n",grid->R[i*NUM_VAR+1]);
	      printf("  Var 3: %.15e\n",grid->R[i*NUM_VAR+2]);
	      printf("  Var 4: %.15e\n",grid->R[i*NUM_VAR+3]);
	    }
	}

      printf("Boundary Residual Check.\n");
      //for ( i=1; i <= (grid->nn + grid->nn_ghost); i++ )
      for ( i=1; i <= grid->nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      if ( !(isfinite(grid->R[i*NUM_VAR+j])) )
		{
		  if ( i <= grid->nn )
		    {
		      printf("  Found a problem at node %d in variable %d: grid->R[i*NUM_VAR+j] = %f\n",i,j,(float)grid->R[i*NUM_VAR+j]);
		      fflush(stdout);
		    }
		  else
		    {
		      printf("  Found a problem at ghost node %d in variable %d: grid->R[i*NUM_VAR+j] = %f\n",i,j,(float)grid->R[i*NUM_VAR+j]);
		      fflush(stdout);
		    }
		}
	    }
	}
      

      //if ( MMS && RECON_PRIM )
      if ( 0 )
	{
	  debug_flux_integral ( grid, p );
	  return 0.;
	}
      
      //do_hard_boundary_flux_jac ( grid, smem, p );
      
      //do_soft_boundary_flux_jac ( grid, smem, p );


      // Generate the L2 norm of the residual vector.
      rms = 0.;
      for ( i=1; i <= grid->nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      rms += ( (grid->R[i*NUM_VAR+j])*(grid->R[i*NUM_VAR+j]) );
	    }
	}
      
      //rms /= (1.*grid->nn);
      //rms = sqrt ( rms );
      

      //======= DEPRECATED =============!!!!!!!!!!
      // Populate the linear system. This contructs the residual vector along with the Jacobian matrix.
      //rms = build_residuals( grid, smem, p );

      // Set our initial guess to zero.
      for ( i=0; i <= grid->nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      grid->dQ[i*NUM_VAR+j] = 0.;
	    }
	}

      // Now we construct the RHS vector.

      // ===================================================================================
      // Now, if I am doing reconstruction on the primitive variables, I need to convert
      // back to conserved variables to construct the actual linear system.
      // ==================================================================================
      if ( RECON_PRIM && p.ic != 3 )
	{
	  for ( i=1; i <= (grid->nn + grid->nn_ghost); i++ )
	    {
	      ConvertPrimitiveToConserved ( p.gamma , &(grid->nQ[i*NUM_VAR]) );
	    }
	}

      // For tacc=2, we need to start off first order for the first iteration.

      if ( ts == (grid->rfiter+1) && p.tacc == 2 )
	{
	  for ( i=1; i <= grid->nn; i++ )
	    {
	      for ( j=0; j < NUM_VAR; j++ )
		{
		  smem->RHS[i*NUM_VAR+j] = -( grid->cv_area[i] / p.mintime * ( grid->nQ[i*NUM_VAR+j] - grid->Q[i*NUM_VAR+j] )
					      + grid->R[i*NUM_VAR+j] );
		}
	    }
	}
      else if ( p.tacc == 2 && ts > 1 )  // 2nd order in time and  iteration > 1.
	{
	  for ( i=1; i <= grid->nn; i++ )
	    {
	      for ( j=0; j < NUM_VAR; j++ )
		{
		  smem->RHS[i*NUM_VAR+j] = -( grid->cv_area[i] / ( 2. *p.mintime ) *
					      ( 3.*grid->nQ[i*NUM_VAR+j] - 4.*grid->Q[i*NUM_VAR+j] + grid->pQ[i*NUM_VAR+j] )
					      + grid->R[i*NUM_VAR+j] );
		}
	    }
	}
      else if ( p.tacc == 1 )  // First order in time.
	{
	  for ( i=1; i <= grid->nn; i++ )
	    {
	      for ( j=0; j < NUM_VAR; j++ )
		{
		  smem->RHS[i*NUM_VAR+j] = -( grid->cv_area[i] / p.mintime * ( grid->nQ[i*NUM_VAR+j] - grid->Q[i*NUM_VAR+j] )
					      + grid->R[i*NUM_VAR+j] );
		}
	    }
	}
      else // Steady state.
	{
	  for ( i=1; i <= grid->nn; i++ )
	    {
	      for ( j=0; j < NUM_VAR; j++ )
		{
		  smem->RHS[i*NUM_VAR+j] = -( grid->R[i*NUM_VAR+j] );
		}
	    }
	}

      // Get the Unsteady residual L2 norm.
      url2 = 0.;
      for ( i=1; i <= grid->nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      url2 += ( (smem->RHS[i*NUM_VAR+j]) * (smem->RHS[i*NUM_VAR+j]) );
	    }
	}
      
      // Generate our solution update vector.
      RMS = linear_solve( grid, smem, p );

      // Get the L2 norm of dQ.
      dql2 = 0.;
      for ( i=1; i <= grid->nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      dql2 += ( grid->dQ[i*NUM_VAR+j] * grid->dQ[i*NUM_VAR+j] );
	    }
	}

      printf("  Newton Iteration %d -> Linear System RMS error = %.15e\n",m,RMS);
      printf("                         Residual L2 Norm = %.15e    CFL = %f\n",rms,(float)grid->cfl);
      printf("                         Unsteady Residual L2 Norm = %.15E\n",sqrt(url2));
      printf("                         L2 Norm of dQ = %.15E\n",sqrt(dql2));
      fflush(stdout);
      
      // Update the solution vector.
      //for ( i=0; i <= (grid->nn + grid->nn_ghost); i++ )
      for ( i=1; i <= grid->nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      grid->nQ[i*NUM_VAR+j] += grid->dQ[i*NUM_VAR+j];
	    }
	}
      
    }         // End of the Newton iterations.

  
  // Set pQ to be the current solution.
  for ( i=1; i <= (grid->nn + grid->nn_ghost); i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  grid->pQ[i*NUM_VAR+j] = grid->Q[i*NUM_VAR+j];
	}
    }

  // Let Q be the solution from the next time step.
  for ( i=1; i <= (grid->nn + grid->nn_ghost); i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  grid->Q[i*NUM_VAR+j] = grid->nQ[i*NUM_VAR+j];
	}
    }

  //return rms;
  return ( sqrt(dql2) );

}

//=============================================================
// 
//  Compute_Extrema()
//
//  Computes the extrema introduced in the field.
//
//  GRID *grid;                          // The grid.
//  PARAMS p;                            // The parameters.
//  int m;                               // The current Newton iteration.
//
//=============================================================

void Compute_Extrema ( GRID *grid, PARAMS p, int m )
{
  int i,j,k;                                    // Loop counters.
  int isubedge;
  int nodeL,nodeR;
  double nx,ny,len;
  double xc,yc,xmid,ymid;
  int gelem;
  int iedge;
  int real_node;
  int ghost_node;
  int b,s,bc;
  double Qmin[NUM_VAR];                        // Minimum values.
  double Qmax[NUM_VAR];                        // Maximum values.
  double qleft[4];                             // Q on the left of the face.
  double qright[4];                            // Q on the right of the face.
  double XL[2],XR[2];
  double GP[6];

  // Variables for higher order flux calculation.
  double x1,x2,x3;                       // Integration points.
  double y1,y2,y3;
  double xL,xR,yL,yR;
  double t1,t2,t3;                       // Gaussian roots.
  double xi, yi;                         // Node coordinates.
  double gp[2];                          // Gauss point.
  
  FILE *fp = NULL;

  // Open the file for writing.
  fp = fopen("extrema.dat","a");

  for ( i=0; i < NUM_VAR; i++ )
    {
      Qmin[i] = 1.0E+10;
      Qmax[i] = 1.0E-10;
    }

  // First lets get the max and min values.
  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  Qmin[j] = MIN( Qmin[j] , grid->nQ[i*NDIM+j] );
	  Qmax[j] = MAX( Qmax[j] , grid->nQ[i*NDIM+j] );
	}
    }

  // Write out some general information to start out.
  fprintf(fp,"Time step = %d , Newton Iteration = %d\n",grid->citer,m);
  fprintf(fp,"Qmax =\n  %.15E\n  %.15E\n  %.15E\n  %.15E\n\n",Qmax[0],Qmax[1],Qmax[2],Qmax[3]);
  fprintf(fp,"Qmin =\n  %.15E\n  %.15E\n  %.15E\n  %.15E\n\n",Qmin[0],Qmin[1],Qmin[2],Qmin[3]);

  // Now we can loop over the subedges and reconstruct to the Gauss points and check
  // for new extrema.

  // Loop over the subedges abd perform the flux calculation.
  for ( isubedge=1; isubedge <= grid->nsubedges; isubedge++ )
    {
      // Extract the real edge this corresponds to.
      iedge = grid->subedges[isubedge*2];

      // Get the normal vector.
      nx = grid->xn_subedges[isubedge*3 + 0];
      ny = grid->xn_subedges[isubedge*3 + 1];
      len = grid->xn_subedges[isubedge*3 + 2];

      // Get the left and right nodes.
      nodeL = grid->edges[iedge*2 + 0];
      nodeR = grid->edges[iedge*2 + 1];

      xmid = grid->xm_subedges[isubedge*NDIM+0];
      ymid = grid->xm_subedges[isubedge*NDIM+1];

      // Find the global element.
      gelem = grid->subedges[isubedge*2+1];

      // Get the element centroid.
      xc = grid->el_cent[gelem*2];
      yc = grid->el_cent[gelem*2+1];

      if ( p.order == 2  && (  (grid->citer >= (p.foits+p.soits)) || FULL_GAUSS_FLUX ) )
	{
	  // Calculate the Gauss points.
	  // For Gaussian quadrature I need to map the points in t-space
	  // onto my dual edge.

	  // I'm also pretty sure that the to be consistent with the right hand
	  // rule the dual edge is defined to be from the edge midpoint to the
	  // element centroid so that the subedge normal vector points away from
	  // the left node.

	  // To get the x,y points from t space use - x = (xL+xR)/2 + (xR-xL)/2*t
	  //                                          y = (yL+yR)/2 + (yR-yL)/2*t
	  //                                          xL,yL = edge midpoint
	  //                                          xR,yR = element centroid
	  
	  t1 = sqrt(15.0)/5.0;
	  t2 = 0.;
	  t3 = sqrt(15.0)/(-5.0);

	  x1 = (xmid+xc)*0.5 + (xc-xmid)*0.5*t1;
	  y1 = (ymid+yc)*0.5 + (yc-ymid)*0.5*t1;
      
	  x2 = (xmid+xc)*0.5 + (xc-xmid)*0.5*t2;
	  y2 = (ymid+yc)*0.5 + (yc-ymid)*0.5*t2;
      
	  x3 = (xmid+xc)*0.5 + (xc-xmid)*0.5*t3;
	  y3 = (ymid+yc)*0.5 + (yc-ymid)*0.5*t3;

	  // For each gauss point, we will reconstruct the solution.
	  
	  gp[0] = x1;  gp[1] = y1;
	  Reconstruct_Gauss ( isubedge, nodeL, nodeR, gp, qleft, qright, grid, p, grid->citer );  

	  // Check for new extrema.

	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      if ( qleft[i] < Qmin[i] )
		{
		  fprintf(fp,"New extrema found on the left at Gauss point 1 for variable %d:\n",i);
		  fprintf(fp,"subedge = %d, edge = %d, nodeL = %d, nodeR = %d\n",isubedge,iedge,nodeL,nodeR);
		  fprintf(fp," at x = %f  y = %f\n",(float)gp[0],(float)gp[1]);
		  fprintf(fp," QL = %.15e     Qmin = %.15e\n\n",qleft[i],Qmin[i]);
		}
	    }

	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      if ( qleft[i] > Qmax[i] )
		{
		  fprintf(fp,"New extrema found on the left at Gauss point 1 for variable %d:\n",i);
		  fprintf(fp,"subedge = %d, edge = %d, nodeL = %d, nodeR = %d\n",isubedge,iedge,nodeL,nodeR);
		  fprintf(fp," at x = %f  y = %f\n",(float)gp[0],(float)gp[1]);
		  fprintf(fp," QL = %.15e     Qmax = %.15e\n\n",qleft[i],Qmax[i]);
		}
	    }

	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      if ( qright[i] < Qmin[i] )
		{
		  fprintf(fp,"New extrema found on the right at Gauss point 1 for variable %d:\n",i);
		  fprintf(fp,"subedge = %d, edge = %d, nodeL = %d, nodeR = %d\n",isubedge,iedge,nodeL,nodeR);
		  fprintf(fp," at x = %f  y = %f\n",(float)gp[0],(float)gp[1]);
		  fprintf(fp," QR = %.15e     Qmin = %.15e\n\n",qright[i],Qmin[i]);
		}
	    }

	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      if ( qright[i] > Qmax[i] )
		{
		  fprintf(fp,"New extrema found on the right at Gauss point 1 for variable %d:\n",i);
		  fprintf(fp,"subedge = %d, edge = %d, nodeL = %d, nodeR = %d\n",isubedge,iedge,nodeL,nodeR);
		  fprintf(fp," at x = %f  y = %f\n",(float)gp[0],(float)gp[1]);
		  fprintf(fp," QR = %.15e     Qmax = %.15e\n\n",qright[i],Qmax[i]);
		}
	    }

	  gp[0] = x2;  gp[1] = y2;
	  Reconstruct_Gauss ( isubedge, nodeL, nodeR, gp, qleft, qright, grid, p, grid->citer );

	   for ( i=0; i < NUM_VAR; i++ )
	    {
	      if ( qleft[i] < Qmin[i] )
		{
		  fprintf(fp,"New extrema found on the left at Gauss point 2 for variable %d:\n",i);
		  fprintf(fp,"subedge = %d, edge = %d, nodeL = %d, nodeR = %d\n",isubedge,iedge,nodeL,nodeR);
		  fprintf(fp," at x = %f  y = %f\n",(float)gp[0],(float)gp[1]);
		  fprintf(fp," QL = %.15e     Qmin = %.15e\n\n",qleft[i],Qmin[i]);
		}
	    }

	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      if ( qleft[i] > Qmax[i] )
		{
		  fprintf(fp,"New extrema found on the left at Gauss point 2 for variable %d:\n",i);
		  fprintf(fp,"subedge = %d, edge = %d, nodeL = %d, nodeR = %d\n",isubedge,iedge,nodeL,nodeR);
		  fprintf(fp," at x = %f  y = %f\n",(float)gp[0],(float)gp[1]);
		  fprintf(fp," QL = %.15e     Qmax = %.15e\n\n",qleft[i],Qmax[i]);
		}
	    }

	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      if ( qright[i] < Qmin[i] )
		{
		  fprintf(fp,"New extrema found on the right at Gauss point 2 for variable %d:\n",i);
		  fprintf(fp,"subedge = %d, edge = %d, nodeL = %d, nodeR = %d\n",isubedge,iedge,nodeL,nodeR);
		  fprintf(fp," at x = %f  y = %f\n",(float)gp[0],(float)gp[1]);
		  fprintf(fp," QR = %.15e     Qmin = %.15e\n\n",qright[i],Qmin[i]);
		}
	    }

	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      if ( qright[i] > Qmax[i] )
		{
		  fprintf(fp,"New extrema found on the right at Gauss point 2 for variable %d:\n",i);
		  fprintf(fp,"subedge = %d, edge = %d, nodeL = %d, nodeR = %d\n",isubedge,iedge,nodeL,nodeR);
		  fprintf(fp," at x = %f  y = %f\n",(float)gp[0],(float)gp[1]);
		  fprintf(fp," QR = %.15e     Qmax = %.15e\n\n",qright[i],Qmax[i]);
		}
	    }

	  gp[0] = x3;  gp[1] = y3;
	  Reconstruct_Gauss ( isubedge, nodeL, nodeR, gp, qleft, qright, grid, p, grid->citer );
	  
	   for ( i=0; i < NUM_VAR; i++ )
	    {
	      if ( qleft[i] < Qmin[i] )
		{
		  fprintf(fp,"New extrema found on the left at Gauss point 3 for variable %d:\n",i);
		  fprintf(fp,"subedge = %d, edge = %d, nodeL = %d, nodeR = %d\n",isubedge,iedge,nodeL,nodeR);
		  fprintf(fp," at x = %f  y = %f\n",(float)gp[0],(float)gp[1]);
		  fprintf(fp," QL = %.15e     Qmin = %.15e\n\n",qleft[i],Qmin[i]);
		}
	    }

	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      if ( qleft[i] > Qmax[i] )
		{
		  fprintf(fp,"New extrema found on the left at Gauss point 3 for variable %d:\n",i);
		  fprintf(fp,"subedge = %d, edge = %d, nodeL = %d, nodeR = %d\n",isubedge,iedge,nodeL,nodeR);
		  fprintf(fp," at x = %f  y = %f\n",(float)gp[0],(float)gp[1]);
		  fprintf(fp," QL = %.15e     Qmax = %.15e\n\n",qleft[i],Qmax[i]);
		}
	    }

	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      if ( qright[i] < Qmin[i] )
		{
		  fprintf(fp,"New extrema found on the right at Gauss point 3 for variable %d:\n",i);
		  fprintf(fp,"subedge = %d, edge = %d, nodeL = %d, nodeR = %d\n",isubedge,iedge,nodeL,nodeR);
		  fprintf(fp," at x = %f  y = %f\n",(float)gp[0],(float)gp[1]);
		  fprintf(fp," QR = %.15e     Qmin = %.15e\n\n",qright[i],Qmin[i]);
		}
	    }

	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      if ( qright[i] > Qmax[i] )
		{
		  fprintf(fp,"New extrema found on the right at Gauss point 3 for variable %d:\n",i);
		  fprintf(fp,"subedge = %d, edge = %d, nodeL = %d, nodeR = %d\n",isubedge,iedge,nodeL,nodeR);
		  fprintf(fp," at x = %f  y = %f\n",(float)gp[0],(float)gp[1]);
		  fprintf(fp," QR = %.15e     Qmax = %.15e\n\n",qright[i],Qmax[i]);
		}
	    }

	}
      else
	{
	  // Get the Q variables and if high order, then reconstruct the solution.
	  Reconstruct( isubedge, nodeL, nodeR, qleft, qright, grid, p, grid->citer);

	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      if ( qleft[i] < Qmin[i] )
		{
		  fprintf(fp,"New extrema found on the left for variable %d:\n",i);
		  fprintf(fp,"subedge = %d, edge = %d, nodeL = %d, nodeR = %d\n",isubedge,iedge,nodeL,nodeR);
		  fprintf(fp," at x = %f  y = %f\n",(float)xmid,(float)ymid);
		  fprintf(fp," QL = %.15e     Qmin = %.15e\n\n",qleft[i],Qmin[i]);
		}
	    }

	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      if ( qleft[i] > Qmax[i] )
		{
		  fprintf(fp,"New extrema found on the left for variable %d:\n",i);
		  fprintf(fp,"subedge = %d, edge = %d, nodeL = %d, nodeR = %d\n",isubedge,iedge,nodeL,nodeR);
		  fprintf(fp," at x = %f  y = %f\n",(float)xmid,(float)ymid);
		  fprintf(fp," QL = %.15e     Qmax = %.15e\n\n",qleft[i],Qmax[i]);
		}
	    }

	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      if ( qright[i] < Qmin[i] )
		{
		  fprintf(fp,"New extrema found on the right for variable %d:\n",i);
		  fprintf(fp,"subedge = %d, edge = %d, nodeL = %d, nodeR = %d\n",isubedge,iedge,nodeL,nodeR);
		  fprintf(fp," at x = %f  y = %f\n",(float)xmid,(float)ymid);
		  fprintf(fp," QR = %.15e     Qmin = %.15e\n\n",qright[i],Qmin[i]);
		}
	    }

	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      if ( qright[i] > Qmax[i] )
		{
		  fprintf(fp,"New extrema found on the right for variable %d:\n",i);
		  fprintf(fp,"subedge = %d, edge = %d, nodeL = %d, nodeR = %d\n",isubedge,iedge,nodeL,nodeR);
		  fprintf(fp," at x = %f  y = %f\n",(float)xmid,(float)ymid);
		  fprintf(fp," QR = %.15e     Qmax = %.15e\n\n",qright[i],Qmax[i]);
		}
	    }
	}

    }

  // Do this also for the boundary edges since I extrapolate for high-order.

  // Loop over the boundary edges, do Roe scheme, and accumulate the results to the real node.
  for ( i=1; i <= grid->nbedges; i++ )
    {
      // Grab the data from the structure.
      real_node = grid->bedges[i*5+0];
      ghost_node = grid->bedges[i*5+1] + grid->nn;
      
      b = grid->bedges[i*5+3];
      bc = grid->bbc[b];
      
      // Get the normal vector information.
      nx = grid->xn_bedges[i*3+0];
      ny = grid->xn_bedges[i*3+1];
      len = grid->xn_bedges[i*3+2];
      
      if ( p.order == 2  && (  (grid->citer >= (p.foits+p.soits)) || FULL_GAUSS_FLUX ) )
	{
	  // Get the segment.
	  s = grid->bedges[i*5+4];
	  
	  // Now get the nodes attached to the edge.
	  nodeL = grid->bs[b][s][0];
	  nodeR = grid->bs[b][s][1];
	  
	  // Get the node coordinates.
	  xi = grid->x[real_node*NDIM+0];
	  yi = grid->x[real_node*NDIM+1];
	  
	  // Get the edge midpoint.
	  xmid = 0.5*( grid->x[nodeL*NDIM+0] + grid->x[nodeR*NDIM+0] );
	  ymid = 0.5*( grid->x[nodeL*NDIM+1] + grid->x[nodeR*NDIM+1] );
	  
	  if ( bc <= 10 ) // Linear boundary.
	    {
	      // For Gaussian quadrature I need to map the points in t-space
	      // onto my dual edge.
	      
	      // I'm also pretty sure that the to be consistent with the right hand
	      // rule the boundary edge is defined from node0 to node1 on the segment.
	      // The normal should thus point out of the mesh.

	      // To get the x,y points from t space use - x = (xL+xR)/2 + (xR-xL)/2*t
	      //                                          y = (yL+yR)/2 + (yR-yL)/2*t
	      // The left and right nodes depend on which node the node is. Is it on the
	      // beginning or end of the edge.
	  
	      if ( real_node == nodeL )
		{
		  xL = xi; xR = xmid;
		  yL = yi; yR = ymid;
		}
	      else
		{
		  xL = xmid; xR = xi;
		  yL = ymid; yR = yi;
		}
	      
	      t1 = sqrt(15.0)/5.0;
	      t2 = 0.;
	      t3 = sqrt(15.0)/(-5.0);

	      x1 = (xL+xR)*0.5 + (xR-xL)*0.5*t1;
	      y1 = (yL+yR)*0.5 + (yR-yL)*0.5*t1;
      
	      x2 = (xL+xR)*0.5 + (xR-xL)*0.5*t2;
	      y2 = (yL+yR)*0.5 + (yR-yL)*0.5*t2;

	      x3 = (xL+xR)*0.5 + (xR-xL)*0.5*t3;
	      y3 = (yL+yR)*0.5 + (yR-yL)*0.5*t3;

	      gp[0] = x1;  gp[1] = y1;
	      Reconstruct_Gauss_Boundary ( real_node, gp, qleft, grid, p);

	      // Check for new extrema.
	      
	      for ( j=0; j < NUM_VAR; j++ )
		{
		  if ( qleft[j] < Qmin[j] )
		    {
		      fprintf(fp,"New extrema found on bedge = %d  at Gauss point 1 for variable %d:\n",i,j);
		      fprintf(fp,"real_node = %d, ghost_node = %d, b = %d, s = %d, bc = %d\n",real_node,ghost_node,b,s,bc);
		      fprintf(fp," at x = %f  y = %f\n",(float)gp[0],(float)gp[1]);
		      fprintf(fp," Q = %.15e     Qmin = %.15e\n\n",qleft[j],Qmin[j]);
		    }
		}

	      for ( j=0; j < NUM_VAR; j++ )
		{
		  if ( qleft[j] > Qmax[j] )
		    {
		      fprintf(fp,"New extrema found on bedge = %d  at Gauss point 1 for variable %d:\n",i,j);
		      fprintf(fp,"real_node = %d, ghost_node = %d, b = %d, s = %d, bc = %d\n",real_node,ghost_node,b,s,bc);
		      fprintf(fp," at x = %f  y = %f\n",(float)gp[0],(float)gp[1]);
		      fprintf(fp," Q = %.15e     Qmax = %.15e\n\n",qleft[j],Qmax[j]);
		    }
		}
	      
	      gp[0] = x2;  gp[1] = y2;
	      Reconstruct_Gauss_Boundary ( real_node, gp, qleft, grid, p);

	      for ( j=0; j < NUM_VAR; j++ )
		{
		  if ( qleft[j] < Qmin[j] )
		    {
		      fprintf(fp,"New extrema found on bedge = %d  at Gauss point 2 for variable %d:\n",i,j);
		      fprintf(fp,"real_node = %d, ghost_node = %d, b = %d, s = %d, bc = %d\n",real_node,ghost_node,b,s,bc);
		      fprintf(fp," at x = %f  y = %f\n",(float)gp[0],(float)gp[1]);
		      fprintf(fp," Q = %.15e     Qmin = %.15e\n\n",qleft[j],Qmin[j]);
		    }
		}

	      for ( j=0; j < NUM_VAR; j++ )
		{
		  if ( qleft[j] > Qmax[j] )
		    {
		      fprintf(fp,"New extrema found on bedge = %d  at Gauss point 2 for variable %d:\n",i,j);
		      fprintf(fp,"real_node = %d, ghost_node = %d, b = %d, s = %d, bc = %d\n",real_node,ghost_node,b,s,bc);
		      fprintf(fp," at x = %f  y = %f\n",(float)gp[0],(float)gp[1]);
		      fprintf(fp," Q = %.15e     Qmax = %.15e\n\n",qleft[j],Qmax[j]);
		    }
		}

	      gp[0] = x3;  gp[1] = y3;
	      Reconstruct_Gauss_Boundary ( real_node, gp, qleft, grid, p);

	      for ( j=0; j < NUM_VAR; j++ )
		{
		  if ( qleft[j] < Qmin[j] )
		    {
		      fprintf(fp,"New extrema found on bedge = %d  at Gauss point 3 for variable %d:\n",i,j);
		      fprintf(fp,"real_node = %d, ghost_node = %d, b = %d, s = %d, bc = %d\n",real_node,ghost_node,b,s,bc);
		      fprintf(fp," at x = %f  y = %f\n",(float)gp[0],(float)gp[1]);
		      fprintf(fp," Q = %.15e     Qmin = %.15e\n\n",qleft[j],Qmin[j]);
		    }
		}

	      for ( j=0; j < NUM_VAR; j++ )
		{
		  if ( qleft[j] > Qmax[j] )
		    {
		      fprintf(fp,"New extrema found on bedge = %d  at Gauss point 3 for variable %d:\n",i,j);
		      fprintf(fp,"real_node = %d, ghost_node = %d, b = %d, s = %d, bc = %d\n",real_node,ghost_node,b,s,bc);
		      fprintf(fp," at x = %f  y = %f\n",(float)gp[0],(float)gp[1]);
		      fprintf(fp," Q = %.15e     Qmax = %.15e\n\n",qleft[j],Qmax[j]);
		    }
		}
	      
	    }
	  else   // Curved boundary.
	    {
	      if ( real_node == nodeL )
		{
		  XL[0] = grid->x[real_node*NDIM+0];
		  XL[1] = grid->x[real_node*NDIM+1];
		  XR[0] = grid->x[ghost_node*NDIM+0];
		  XR[1] = grid->x[ghost_node*NDIM+1];
		}
	      else
		{
		  XR[0] = grid->x[real_node*NDIM+0];
		  XR[1] = grid->x[real_node*NDIM+1];
		  XL[0] = grid->x[ghost_node*NDIM+0];
		  XL[1] = grid->x[ghost_node*NDIM+1];
		}

	      // Get the Gauss points.
	      curved_boundary_gauss_points ( grid, bc, XL, XR, GP );

	      curved_boundary_arclength ( grid, bc, XL, XR, &len );

	      for ( j=0; j < 3; j++ )
		{
		  // Get the normal vector.
		  curved_boundary_normal_vector ( grid, bc, &(GP[j*NDIM]), &nx, &ny );
		  
		  Reconstruct_Gauss_Boundary ( real_node, &(GP[j*NDIM]), qleft, grid, p);
		  
		  for ( k=0; k < NUM_VAR; k++ )
		    {
		      if ( qleft[k] < Qmin[k] )
			{
			  fprintf(fp,"New extrema found on bedge = %d  at Gauss point %d for variable %d:\n",i,j+1,k);
			  fprintf(fp,"real_node = %d, ghost_node = %d, b = %d, s = %d, bc = %d\n",real_node,ghost_node,b,s,bc);
			  fprintf(fp," at x = %f  y = %f\n",(float)GP[j*NDIM+0],(float)GP[j*NDIM+1]);
			  fprintf(fp," Q = %.15e     Qmin = %.15e\n\n",qleft[k],Qmin[k]);
			}
		    }
		  
		  for ( k=0; k < NUM_VAR; k++ )
		    {
		      if ( qleft[k] > Qmax[k] )
			{
			  fprintf(fp,"New extrema found on bedge = %d  at Gauss point %d for variable %d:\n",i,j+1,k);
			  fprintf(fp,"real_node = %d, ghost_node = %d, b = %d, s = %d, bc = %d\n",real_node,ghost_node,b,s,bc);
			  fprintf(fp," at x = %f  y = %f\n",(float)GP[j*NDIM+0],(float)GP[j*NDIM+1]);
			  fprintf(fp," Q = %.15e     Qmax = %.15e\n\n",qleft[k],Qmax[k]);
			}
		    }
		  
		  
		}
	      
	    }
	  
	}
    }

  // Finished checking for extrema.
  fclose(fp);
  
  return;

}



//========================================================
//
//  Short function to verify that Int_CV del*F*dA = Int_dCV F*nhat*dS
//
//========================================================

void debug_flux_integral ( GRID *grid, PARAMS p )
{
  int i=0,j=0,k=0,e=0;
  int v;                                 // Vector loop.
  int ind;                               // Component loop.
  int nodes[MAX_NUM_VERT] = {0};         // Element vertices.
  int node;                              // Node index.
  int num_vert;                          // Number of vertices for an element.
  int left_id;                           // ID of left edge.
  int right_id;                          // ID of right edge.
  double xc[NDIM];                       // Vectors.
  double xmidL[NDIM];
  double xmidR[NDIM];
  double v_left[NDIM];
  double v_right[NDIM];
  double v_middle[NDIM];
  double bary_coord[7][3];               // The barycentric coordinates.
  double wg[7];                          // The associated weights;
  double b_r = (6. - sqrt(15.))/21.;     // Precompute some the coordinates.
  double b_t = (6. + sqrt(15.))/21.;
  double b_s = (9. + 2.*sqrt(15.))/21.;
  double b_u = (9. - 2.*sqrt(15.))/21.;
  double wA = (155. - sqrt(15.))/1200.;
  double wB = (155. + sqrt(15.))/1200.;
  double norm=0., temp=0.;
  double XGP[NDIM];
  double QGP[NUM_VAR];
  double QTRI[NUM_VAR];

  double gp[7][3];                       // Gauss points for integration.
  double dA;                             // Differential area (for a triangle).

  double *del_F = NULL;

  int *node_status = NULL;

  // Curved triangle stuff.
  double ct_gp[7][2];                    // Gauss integration points in (r,s) space (xi,eta).
  double ct_nodes[7][2];                 // The x,y locations of the 6 nodes defining the p2 triangle (curved side).
  double N1,N2,N3,N4,N5,N6;              // Shape functions.
  double N1r,N1s,N2r,N2s,N3r,N3s;        // Shape function derivatives.
  double N4r,N4s,N5r,N5s,N6r,N6s;
  double xr,xs,yr,ys,jacobian;           // Transformation jacobian.
  double XI,ETA;                         // Triangle coordinates.
  int b,bc,seg,ct_flag=0;
  int ghost_node;

  //double LS;

  //LS = 1./6.;
  //LS = 1./12.;
  //LS = 1./24.;
  //LS = 1./48.;
  //LS = 1./96.;

  //printf("USING LENGTH SCALE = %.15E\n",LS);
  

  node_status = (int*)malloc((grid->nn+1)*sizeof(int));
  if ( node_status == NULL ) { printf("node_status didnt allocate.\n"); exit(0); }

  del_F = (double*)calloc((grid->nn+1)*NUM_VAR,sizeof(double));
  if ( del_F == NULL ) { printf("del_F didnt allocate.\n"); exit(0); }

  for ( i=1; i <= grid->nn; i++ )
    {
      node_status[i] = 1;  // everything is interior.
      del_F[i] = 0.;
    }

  for ( i=1; i <= grid->nbedges; i++ )
    {
      node_status[grid->bedges[i*5+0]] = 0;  // turn off boundary nodes.
    }


  // Now we calculate the integral of the flux divergence over the CV using Gaussian quadrature for triangles.
  
  // Set up the barycentric coordinates and their respective weights.
  bary_coord[0][0] = 1./3.; bary_coord[0][1] = 1./3.; bary_coord[0][2] = 1./3.;
  bary_coord[1][0] = b_r;   bary_coord[1][1] = b_r;   bary_coord[1][2] = b_s;
  bary_coord[2][0] = b_r;   bary_coord[2][1] = b_s;   bary_coord[2][2] = b_r;
  bary_coord[3][0] = b_s;   bary_coord[3][1] = b_r;   bary_coord[3][2] = b_r;
  bary_coord[4][0] = b_t;   bary_coord[4][1] = b_t;   bary_coord[4][2] = b_u;
  bary_coord[5][0] = b_t;   bary_coord[5][1] = b_u;   bary_coord[5][2] = b_t;
  bary_coord[6][0] = b_u;   bary_coord[6][1] = b_t;   bary_coord[6][2] = b_t;
      
  wg[0] = 9./40.;
  wg[1] = wA;
  wg[2] = wA;
  wg[3] = wA;
  wg[4] = wB;
  wg[5] = wB;
  wg[6] = wB;

  // Lets also map these guys onto the reference triangle for integration on the triangles with curved sides.
  for ( j=0; j < 7; j++ )
    {
      ct_gp[j][0] = (bary_coord[j][0] * 0. +
		     bary_coord[j][1] * 1. +
		     bary_coord[j][2] * 0. );
      
      ct_gp[j][1] = (bary_coord[j][0] * 0. +
		     bary_coord[j][1] * 0. +
		     bary_coord[j][2] * 1. );
      
    }  // The weights can be reused since everything is consistent in the ordering.

  // Start the process by looping over all the elements and looping through the internal median dual pieces.
  for ( e=Tri; e <= Quad; e++ )                 // Element types loop.
    {
      // Find the number of nodes and edges for this element type.
      num_vert = NumberNodesForElement(e);

      for ( j=1; j <= grid->num_elem[e]; j++ )              // Element loop.
	{
	  // Retrieve the element node ids. This is for convenience.
	  for ( k=0; k < num_vert; k++ )              // Element node extraction loop.
	    {
	      nodes[k] = grid->c2n[e][j*num_vert+k];
	    }                                         // End element node extraction loop.
	  
	  // Find the element centroid.
	  xc[0] = 0.; xc[1] = 0.;
	  for ( k=0; k < num_vert; k++ )
	    {
	      for ( i=0; i < NDIM; i++ )
		{
		  xc[i] += grid->x[NDIM*nodes[k]+i];
		}
	    }
	  xc[0] /= ((double)num_vert);
	  xc[1] /= ((double)num_vert);

	  // Loop through each node in the element.
	  for ( k=0; k < num_vert; k++ )              // Element node loop.
	    {
	      node = nodes[k];

	      // Get the two neighboring vertices.
	      left_id = (k+1)%num_vert;
	      right_id = (k+num_vert-1)%num_vert;
		  
	      left_id =  grid->c2n[e][j*num_vert+left_id];
	      right_id = grid->c2n[e][j*num_vert+right_id];
	      
	      // Get the edge midpoints.
	      for ( i=0; i < NDIM; i++ )
		{
		  xmidL[i] = grid->x[NDIM*node+i] + grid->x[NDIM*left_id+i];
		  xmidR[i] = grid->x[NDIM*node+i] + grid->x[NDIM*right_id+i];
		}

	      for ( i=0; i < NDIM; i++ )
		{
		  xmidL[i] *= 0.5;
		  xmidR[i] *= 0.5;
		}

	      // Make the vectors.
	      for ( i=0; i < NDIM; i++ )
		{
		  v_left[i]   = xmidL[i] - grid->x[NDIM*node+i];
		  v_right[i]  = xmidR[i] - grid->x[NDIM*node+i];
		  v_middle[i] = xc[i]    - grid->x[NDIM*node+i];
		}

	      ct_flag = 0;
		  
	      // If the node and its neighbor are both boundary nodes, we need to gather some information.
	      if ( grid->node_state[node] == BOUNDARY && grid->node_state[left_id] == BOUNDARY )
		{
		  // Retrieve the relevant information.
		  
		  // First find the correct bedge. It will have both nodes.
		  for ( i=1; i <= grid->nbedges; i++ )
		    {
		      if ( node == grid->bedges[i*5] )  // candidate bedge.
			{
			  b = grid->bedges[i*5+3];
			  seg = grid->bedges[i*5+4];
			  
			  if ( (  node == grid->bs[b][seg][0]   ||  node == grid->bs[b][seg][1]   ) &&
			       ( left_id == grid->bs[b][seg][0] || left_id == grid->bs[b][seg][1] )   )  // we have a winner.
			    {
			      ghost_node = grid->bedges[i*5+1] + grid->nn;
			      b = grid->bedges[i*5+3];
			      bc = grid->bbc[b];
			      seg = grid->bedges[i*5+4];
			      
			      if ( bc > 10 )
				ct_flag = 1;
			      
			      break;
			    }
			}
		    }
		}

	      if ( ct_flag )
		{
		  // Now we need to set up the six nodes used by the shape functions.
		      
		  // Node 1 is the element centroid.
		  ct_nodes[1][0] = xc[0];
		  ct_nodes[1][1] = xc[1];

		  // Node 2 is the node we are working on.
		  ct_nodes[2][0] = grid->x[NDIM*node+0];
		  ct_nodes[2][1] = grid->x[NDIM*node+1];

		  // Third node is the ghost node.
		  ct_nodes[3][0] = grid->x[NDIM*ghost_node+0];
		  ct_nodes[3][1] = grid->x[NDIM*ghost_node+1];

		  // Node 4 is the midpoint of 1,2.
		  ct_nodes[4][0] = 0.5*(ct_nodes[1][0] + ct_nodes[2][0]);
		  ct_nodes[4][1] = 0.5*(ct_nodes[1][1] + ct_nodes[2][1]);

		  // Node 5 is the midpoint of the curved boundary between node and ghost node.
		  curved_boundary_midpoint(grid, bc, &(ct_nodes[2][0]), &(ct_nodes[3][0]), &(ct_nodes[5][0]) );

		  // Node 6 is the midpoint of 1,3.
		  ct_nodes[6][0] = 0.5*(ct_nodes[1][0] + ct_nodes[3][0]);
		  ct_nodes[6][1] = 0.5*(ct_nodes[1][1] + ct_nodes[3][1]);

		  // First we need to make sure that the concavity constraint is met. Under a linear transformation
		  // to the reference element, node 5 must be greater than 0.25 in both r and s coordinates.
		  // See Flaherty's notes online or Claes Johnson's book.

		  XI = ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
			 (ct_nodes[5][0]*ct_nodes[3][1] - ct_nodes[5][1]*ct_nodes[3][0]) +
			 (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
		       ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
		         (ct_nodes[2][0]*ct_nodes[3][1] - ct_nodes[2][1]*ct_nodes[3][0]) +
		         (ct_nodes[2][0]*ct_nodes[1][1] - ct_nodes[2][1]*ct_nodes[1][0])  );
		  
		  ETA = ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			  (ct_nodes[5][0]*ct_nodes[2][1] - ct_nodes[5][1]*ct_nodes[2][0]) +
			  (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
		        ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			  (ct_nodes[3][0]*ct_nodes[2][1] - ct_nodes[3][1]*ct_nodes[2][0]) +
			  (ct_nodes[3][0]*ct_nodes[1][1] - ct_nodes[1][0]*ct_nodes[3][1])  );

		  if ( ETA < 0.25 || XI < 0.25 )  // The transformation will not be invertible.
		    {
		      printf("CRITICAL ERROR: Triangle 1 of element %d of type %d is too concave (compute cv averages: hessian.C)!\n",j,e);
		      printf("  xi = %.15e     eta = %.15e\n",XI,ETA);
		      fflush(stdout);
		      exit(1);
		    }

		  // Do this for each variable.

		  QTRI[0] = 0.; QTRI[1] = 0.; QTRI[2] = 0.; QTRI[3] = 0.;
		  
		  // Now I can loop over the Gaussian integration nodes and apply the quadratic transformation.
		  for ( i=0; i < 7; i++ )
		    {
		      XI = ct_gp[i][0];
		      ETA = ct_gp[i][1];
		      
		      // Compute basis functions.
		      N1 = 1. - 3.*XI - 3.*ETA + 4.*XI*ETA + 2.*XI*XI + 2.*ETA*ETA;
		      N2 = XI*(2.*XI - 1.);
		      N3 = ETA*(2.*ETA - 1.);
		      N4 = 4.*XI*(1. - XI - ETA);
		      N5 = 4.*XI*ETA;
		      N6 = 4.*ETA*(1. - XI - ETA);

		      // Compute their derivatives.
		      N1r = -3. + 4.*ETA + 4.*XI;
		      N1s = -3. + 4.*XI + 4.*ETA;
		      N2r = 4.*XI - 1.;
		      N2s = 0.;
		      N3r = 0.;
		      N3s = 4.*ETA - 1.;
		      N4r = 4.*(1. - 2.*XI - ETA);
		      N4s = -4.*XI;
		      N5r = 4.*ETA;
		      N5s = 4.*XI;
		      N6r = -4.*ETA;
		      N6s = 4.*(1. - 2.*ETA - XI);

		      // Compute the derivatives of x,y wrt to the transformation.
		      xr = N1r*ct_nodes[1][0] + N2r*ct_nodes[2][0] + N3r*ct_nodes[3][0] + N4r*ct_nodes[4][0]
			 + N5r*ct_nodes[5][0] + N6r*ct_nodes[6][0];
		      xs = N1s*ct_nodes[1][0] + N2s*ct_nodes[2][0] + N3s*ct_nodes[3][0] + N4s*ct_nodes[4][0]
			 + N5s*ct_nodes[5][0] + N6s*ct_nodes[6][0];
		      yr = N1r*ct_nodes[1][1] + N2r*ct_nodes[2][1] + N3r*ct_nodes[3][1] + N4r*ct_nodes[4][1]
			 + N5r*ct_nodes[5][1] + N6r*ct_nodes[6][1];
		      ys = N1s*ct_nodes[1][1] + N2s*ct_nodes[2][1] + N3s*ct_nodes[3][1] + N4s*ct_nodes[4][1]
			 + N5s*ct_nodes[5][1] + N6s*ct_nodes[6][1];

		      jacobian = xr*ys - xs*yr;

		      // Now I need to extrapolate to this point. First we get the real coordinates of the
		      // Gaussian point.
		      XGP[0] = N1*ct_nodes[1][0] + N2*ct_nodes[2][0] + N3*ct_nodes[3][0] +
			       N4*ct_nodes[4][0] + N5*ct_nodes[5][0] + N6*ct_nodes[6][0];

		      XGP[1] = N1*ct_nodes[1][1] + N2*ct_nodes[2][1] + N3*ct_nodes[3][1] +
			       N4*ct_nodes[4][1] + N5*ct_nodes[5][1] + N6*ct_nodes[6][1];

		      flux_divergence ( XGP, QGP );

		      // Accumulate.
		      for ( v=0; v < NUM_VAR; v++ )
			{
			  QTRI[v] += ( wg[i] * QGP[v] *jacobian * 0.5);
			} 
		      
		    }

		  // Add to the running total for the node.
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      del_F[node*NUM_VAR+v] += QTRI[v];
		    }
		  
		}
	      else
		{
		  // Get the area.
		  dA = 0.5*( v_left[0]*v_middle[1] - v_left[1]*v_middle[0] );

		  // Initialize the Gauss points.
		  for ( v=0; v < 7; v++ )
		    {
		      for ( ind=0; ind < NDIM; ind++ )
			{
			  gp[v][ind] = (bary_coord[v][0] * grid->x[node*NDIM+ind] +
					bary_coord[v][1] * xmidL[ind] +
					bary_coord[v][2] * xc[ind]);
			}
		      gp[v][2] = 0.;   // Always 0 since 2D vectors.
		    }

		  // Now we loop over the Gauss points and get the function value
		  // and perform the numerical integration (the area integral of the triangle).
	      
		  // Do this for each variable.

		  QTRI[0] = 0.; QTRI[1] = 0.; QTRI[2] = 0.; QTRI[3] = 0.;
		  
		  for ( i=0; i < 7; i++ )
		    {
		      // Get the value at the point.
		      XGP[0] = gp[i][0];
		      XGP[1] = gp[i][1];

		      flux_divergence ( XGP, QGP );
		      
		      for ( v=0; v < NUM_VAR; v++ )
			{
			  QTRI[v] += ( wg[i] * QGP[v] );
			} 
		    }
	      
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      QTRI[v] *= dA;
		    }
	      
		  // Accumulate the integral to the node.
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      del_F[node*NUM_VAR+v] += QTRI[v];
		    }
		}

	      // Process the second triangle.

	      ct_flag = 0;
	      
	      // If the node and its neighbor are both boundary nodes, we need to gather some information.
	      if ( grid->node_state[node] == BOUNDARY && grid->node_state[right_id] == BOUNDARY )
		{
		  // Retrieve the relevant information.
		  
		  // First find the correct bedge. It will have both nodes.
		  for ( i=1; i <= grid->nbedges; i++ )
		    {
		      if ( node == grid->bedges[i*5] )  // candidate bedge.
			{
			  b = grid->bedges[i*5+3];
			  seg = grid->bedges[i*5+4];
			  
			  if ( (   node == grid->bs[b][seg][0]   ||   node == grid->bs[b][seg][1]   ) &&
			       ( right_id == grid->bs[b][seg][0] || right_id == grid->bs[b][seg][1] )   )  // we have a winner.
			    {
			      ghost_node = grid->bedges[i*5+1] + grid->nn;
			      b = grid->bedges[i*5+3];
			      bc = grid->bbc[b];
			      seg = grid->bedges[i*5+4];
			      if ( bc > 10 )
				ct_flag = 1;
			      
			      break;
			    }
			}
		    }
		}

	      if ( ct_flag )
		{
		  // Now we need to set up the six nodes used by the shape functions.
		      
		  // Node 1 is the element centroid.
		  ct_nodes[1][0] = xc[0];
		  ct_nodes[1][1] = xc[1];

		  // Node 2 is the ghost_node.
		  ct_nodes[2][0] = grid->x[NDIM*ghost_node+0];
		  ct_nodes[2][1] = grid->x[NDIM*ghost_node+1];

		  // Third node is the node we are working on.
		  ct_nodes[3][0] = grid->x[NDIM*node+0];
		  ct_nodes[3][1] = grid->x[NDIM*node+1];

		  // Node 4 is the midpoint of 1,2.
		  ct_nodes[4][0] = 0.5*(ct_nodes[1][0] + ct_nodes[2][0]);
		  ct_nodes[4][1] = 0.5*(ct_nodes[1][1] + ct_nodes[2][1]);

		  // Node 5 is the midpoint of the curved boundary between node and ghost node.
		  curved_boundary_midpoint(grid, bc, &(ct_nodes[2][0]), &(ct_nodes[3][0]), &(ct_nodes[5][0]) );

		  // Node 6 is the midpoint of 1,3.
		  ct_nodes[6][0] = 0.5*(ct_nodes[1][0] + ct_nodes[3][0]);
		  ct_nodes[6][1] = 0.5*(ct_nodes[1][1] + ct_nodes[3][1]);

		  // First we need to make sure that the concavity constraint is met. Under a linear transformation
		  // to the reference element, node 5 must be greater than 0.25 in both r and s coordinates.
		  // See Flaherty's notes online or Claes Johnson's book.

		  XI = ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
			 (ct_nodes[5][0]*ct_nodes[3][1] - ct_nodes[5][1]*ct_nodes[3][0]) +
			 (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
		       ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
		         (ct_nodes[2][0]*ct_nodes[3][1] - ct_nodes[2][1]*ct_nodes[3][0]) +
		         (ct_nodes[2][0]*ct_nodes[1][1] - ct_nodes[2][1]*ct_nodes[1][0])  );
		  
		  ETA = ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			  (ct_nodes[5][0]*ct_nodes[2][1] - ct_nodes[5][1]*ct_nodes[2][0]) +
			  (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
		        ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			  (ct_nodes[3][0]*ct_nodes[2][1] - ct_nodes[3][1]*ct_nodes[2][0]) +
			  (ct_nodes[3][0]*ct_nodes[1][1] - ct_nodes[1][0]*ct_nodes[3][1])  );

		  if ( ETA < 0.25 || XI < 0.25 )  // The transformation will not be invertible.
		    {
		      printf("CRITICAL ERROR: Triangle 1 of element %d of type %d is too concave (compute cv averages: hessian.C)!\n",j,e);
		      printf("  xi = %.15e     eta = %.15e\n",XI,ETA);
		      fflush(stdout);
		      exit(1);
		    }

		  // Do this for each variable.

		  QTRI[0] = 0.; QTRI[1] = 0.; QTRI[2] = 0.; QTRI[3] = 0.;
		  
		  // Now I can loop over the Gaussian integration nodes and apply the quadratic transformation.
		  for ( i=0; i < 7; i++ )
		    {
		      XI = ct_gp[i][0];
		      ETA = ct_gp[i][1];
		      
		      // Compute basis functions.
		      N1 = 1. - 3.*XI - 3.*ETA + 4.*XI*ETA + 2.*XI*XI + 2.*ETA*ETA;
		      N2 = XI*(2.*XI - 1.);
		      N3 = ETA*(2.*ETA - 1.);
		      N4 = 4.*XI*(1. - XI - ETA);
		      N5 = 4.*XI*ETA;
		      N6 = 4.*ETA*(1. - XI - ETA);

		      // Compute their derivatives.
		      N1r = -3. + 4.*ETA + 4.*XI;
		      N1s = -3. + 4.*XI + 4.*ETA;
		      N2r = 4.*XI - 1.;
		      N2s = 0.;
		      N3r = 0.;
		      N3s = 4.*ETA - 1.;
		      N4r = 4.*(1. - 2.*XI - ETA);
		      N4s = -4.*XI;
		      N5r = 4.*ETA;
		      N5s = 4.*XI;
		      N6r = -4.*ETA;
		      N6s = 4.*(1. - 2.*ETA - XI);

		      // Compute the derivatives of x,y wrt to the transformation.
		      xr = N1r*ct_nodes[1][0] + N2r*ct_nodes[2][0] + N3r*ct_nodes[3][0] + N4r*ct_nodes[4][0]
			 + N5r*ct_nodes[5][0] + N6r*ct_nodes[6][0];
		      xs = N1s*ct_nodes[1][0] + N2s*ct_nodes[2][0] + N3s*ct_nodes[3][0] + N4s*ct_nodes[4][0]
			 + N5s*ct_nodes[5][0] + N6s*ct_nodes[6][0];
		      yr = N1r*ct_nodes[1][1] + N2r*ct_nodes[2][1] + N3r*ct_nodes[3][1] + N4r*ct_nodes[4][1]
			 + N5r*ct_nodes[5][1] + N6r*ct_nodes[6][1];
		      ys = N1s*ct_nodes[1][1] + N2s*ct_nodes[2][1] + N3s*ct_nodes[3][1] + N4s*ct_nodes[4][1]
			 + N5s*ct_nodes[5][1] + N6s*ct_nodes[6][1];

		      jacobian = xr*ys - xs*yr;

		      // Now I need to extrapolate to this point. First we get the real coordinates of the
		      // Gaussian point.
		      XGP[0] = N1*ct_nodes[1][0] + N2*ct_nodes[2][0] + N3*ct_nodes[3][0] +
			       N4*ct_nodes[4][0] + N5*ct_nodes[5][0] + N6*ct_nodes[6][0];

		      XGP[1] = N1*ct_nodes[1][1] + N2*ct_nodes[2][1] + N3*ct_nodes[3][1] +
			       N4*ct_nodes[4][1] + N5*ct_nodes[5][1] + N6*ct_nodes[6][1];

		      flux_divergence ( XGP, QGP );

		      // Accumulate.
		      for ( v=0; v < NUM_VAR; v++ )
			{
			  QTRI[v] += ( wg[i] * QGP[v] *jacobian * 0.5);
			} 
		      
		    }

		  // Add to the running total for the node.
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      del_F[node*NUM_VAR+v] += QTRI[v];
		    }
		  
		}
	      else
		{
		  // Get the area.
		  dA = 0.5*( v_middle[0]*v_right[1] - v_middle[1]*v_right[0] );
		  
		  // Initialize the Gauss points.
		  for ( v=0; v < 7; v++ )
		    {
		      for ( ind=0; ind < NDIM; ind++ )
			{
			  gp[v][ind] = (bary_coord[v][0] * grid->x[node*NDIM+ind] +
					bary_coord[v][1] * xmidR[ind] +
					bary_coord[v][2] * xc[ind]);
			}
		      gp[v][2] = 0.;   // Always 0 since 2D vectors.
		    }
		  
		  // Now we loop over the Gauss points and get the function value
		  // and perform the numerical integration (the area integral of the triangle).
		  
		  // Do this for each variable.
		  
		  QTRI[0] = 0.; QTRI[1] = 0.; QTRI[2] = 0.; QTRI[3] = 0.;
		  
		  for ( i=0; i < 7; i++ )
		    {
		      // Get the value at the point.
		      XGP[0] = gp[i][0];
		      XGP[1] = gp[i][1];
		  
		      flux_divergence ( XGP, QGP );
		  
		      for ( v=0; v < NUM_VAR; v++ )
			{
			  QTRI[v] += ( wg[i] * QGP[v] );
			} 
		    }
	      
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      QTRI[v] *= dA;
		    }
		  
		  // Accumulate the integral to the node.
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      del_F[node*NUM_VAR+v] += QTRI[v];
		    }
		}
	    }                                       // End element node loop.
	  
	}                                           // End element loop.
      
    }                                               // End element type loop.
  
  // Now lets add up all the contributions to the flux divergence.
  
  QTRI[0] = 0.; QTRI[1] = 0.; QTRI[2] = 0.; QTRI[3] = 0.;

  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	QTRI[j] += del_F[i*NUM_VAR+j];
    }
  
  printf("TOTAL DIVERGENCE OF THE FLUX IS\n");
  for ( i=0; i < NUM_VAR; i++ )
    {
      j = i+1;
      printf("  VARIABLE %d : %.15E\n",j,QTRI[i]);
    }
  
  // Now we can loop over and compare the difference between the two vectors.
  
  printf("Flux integral test:\n\n");

  printf("Unweighted test:\n");
  // Rho.
  printf("Density:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      norm += fabs( grid->R[i*NUM_VAR+0] - del_F[i*NUM_VAR+0] );
    }
  printf("  L1 norm of difference: %.15e\n",norm);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      temp = fabs( grid->R[i*NUM_VAR+0] - del_F[i*NUM_VAR+0] );
      norm += ( temp*temp );
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm));
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      temp = fabs( grid->R[i*NUM_VAR+0] - del_F[i*NUM_VAR+0] );
      if ( temp > norm )
	{
	  norm = temp;
	  j = i;
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm,j);

  // u
  printf("U:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      norm += fabs( grid->R[i*NUM_VAR+1] - del_F[i*NUM_VAR+1] );
    }
  printf("  L1 norm of difference: %.15e\n",norm);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      temp = fabs( grid->R[i*NUM_VAR+1] - del_F[i*NUM_VAR+1] );
      norm += ( temp*temp );
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm));
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      temp = fabs( grid->R[i*NUM_VAR+1] - del_F[i*NUM_VAR+1] );
      if ( temp > norm )
	{
	  norm = temp;
	  j = i;
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm,j);

  printf("v:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      norm += fabs( grid->R[i*NUM_VAR+2] - del_F[i*NUM_VAR+2] );
    }
  printf("  L1 norm of difference: %.15e\n",norm);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      temp = fabs( grid->R[i*NUM_VAR+2] - del_F[i*NUM_VAR+2] );
      norm += ( temp*temp );
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm));
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {      
      temp = fabs( grid->R[i*NUM_VAR+2] - del_F[i*NUM_VAR+2] );
      if ( temp > norm )
	{
	  norm = temp;
	  j = i;
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm,j);
  

  printf("Pressure:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      norm += fabs( grid->R[i*NUM_VAR+3] - del_F[i*NUM_VAR+3] );
    }
  printf("  L1 norm of difference: %.15e\n",norm);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      temp = fabs( grid->R[i*NUM_VAR+3] - del_F[i*NUM_VAR+3] );
      norm += ( temp*temp );
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm));
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      temp = fabs( grid->R[i*NUM_VAR+3] - del_F[i*NUM_VAR+3] );
      if ( temp > norm )
	{
	  norm = temp;
	  j = i;
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm,j);

  /*
  printf("Weighted by length scale:\n");
  // Rho.
  printf("Density:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	norm += fabs( grid->R[i*NUM_VAR+0] - del_F[i*NUM_VAR+0] );
    }
  printf("  L1 norm of difference: %.15e\n",norm*LS);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+0] - del_F[i*NUM_VAR+0] );
	  norm += ( temp*temp );
	}
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm)*LS);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+0] - del_F[i*NUM_VAR+0] );
	  if ( temp > norm )
	    {
	      norm = temp;
	      j = i;
	    }
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm*LS,j);

  // u
  printf("U:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	norm += fabs( grid->R[i*NUM_VAR+1] - del_F[i*NUM_VAR+1] );
    }
  printf("  L1 norm of difference: %.15e\n",norm*LS);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+1] - del_F[i*NUM_VAR+1] );
	  norm += ( temp*temp );
	}
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm)*LS);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+1] - del_F[i*NUM_VAR+1] );
	  if ( temp > norm )
	    {
	      norm = temp;
	      j = i;
	    }
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm*LS,j);

  printf("v:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	norm += fabs( grid->R[i*NUM_VAR+2] - del_F[i*NUM_VAR+2] );
    }
  printf("  L1 norm of difference: %.15e\n",norm*LS);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+2] - del_F[i*NUM_VAR+2] );
	  norm += ( temp*temp );
	}
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm)*LS);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+2] - del_F[i*NUM_VAR+2] );
	  if ( temp > norm )
	    {
	      norm = temp;
	      j = i;
	    }
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm*LS,j);
  

  printf("Pressure:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	norm += fabs( grid->R[i*NUM_VAR+3] - del_F[i*NUM_VAR+3] );
    }
  printf("  L1 norm of difference: %.15e\n",norm*LS);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+3] - del_F[i*NUM_VAR+3] );
	  norm += ( temp*temp );
	}
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm)*LS);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+3] - del_F[i*NUM_VAR+3] );
	  if ( temp > norm )
	    {
	      norm = temp;
	      j = i;
	    }
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm*LS,j);



  printf("weighted by length scale squared:\n");
  // Rho.
  printf("Density:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	norm += fabs( grid->R[i*NUM_VAR+0] - del_F[i*NUM_VAR+0] );
    }
  printf("  L1 norm of difference: %.15e\n",norm*LS*LS);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+0] - del_F[i*NUM_VAR+0] );
	  norm += ( temp*temp );
	}
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm)*LS*LS);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+0] - del_F[i*NUM_VAR+0] );
	  if ( temp > norm )
	    {
	      norm = temp;
	      j = i;
	    }
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm*LS*LS,j);

  // u
  printf("U:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	norm += fabs( grid->R[i*NUM_VAR+1] - del_F[i*NUM_VAR+1] );
    }
  printf("  L1 norm of difference: %.15e\n",norm*LS*LS);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+1] - del_F[i*NUM_VAR+1] );
	  norm += ( temp*temp );
	}
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm)*LS*LS);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+1] - del_F[i*NUM_VAR+1] );
	  if ( temp > norm )
	    {
	      norm = temp;
	      j = i;
	    }
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm*LS*LS,j);

  printf("v:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	norm += fabs( grid->R[i*NUM_VAR+2] - del_F[i*NUM_VAR+2] );
    }
  printf("  L1 norm of difference: %.15e\n",norm*LS*LS);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+2] - del_F[i*NUM_VAR+2] );
	  norm += ( temp*temp );
	}
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm)*LS*LS);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+2] - del_F[i*NUM_VAR+2] );
	  if ( temp > norm )
	    {
	      norm = temp;
	      j = i;
	    }
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm*LS*LS,j);
  

  printf("Pressure:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	norm += fabs( grid->R[i*NUM_VAR+3] - del_F[i*NUM_VAR+3] );
    }
  printf("  L1 norm of difference: %.15e\n",norm*LS*LS);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+3] - del_F[i*NUM_VAR+3] );
	  norm += ( temp*temp );
	}
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm)*LS*LS);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+3] - del_F[i*NUM_VAR+3] );
	  if ( temp > norm )
	    {
	      norm = temp;
	      j = i;
	    }
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm*LS*LS,j);



  printf("weighted by sqrt(length scale):\n");
  // Rho.
  printf("Density:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	norm += fabs( grid->R[i*NUM_VAR+0] - del_F[i*NUM_VAR+0] );
    }
  printf("  L1 norm of difference: %.15e\n",norm*sqrt(LS));
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+0] - del_F[i*NUM_VAR+0] );
	  norm += ( temp*temp );
	}
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm)*sqrt(LS));
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+0] - del_F[i*NUM_VAR+0] );
	  if ( temp > norm )
	    {
	      norm = temp;
	      j = i;
	    }
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm*sqrt(LS),j);

  // u
  printf("U:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	norm += fabs( grid->R[i*NUM_VAR+1] - del_F[i*NUM_VAR+1] );
    }
  printf("  L1 norm of difference: %.15e\n",norm*sqrt(LS));
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+1] - del_F[i*NUM_VAR+1] );
	  norm += ( temp*temp );
	}
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm)*sqrt(LS));
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+1] - del_F[i*NUM_VAR+1] );
	  if ( temp > norm )
	    {
	      norm = temp;
	      j = i;
	    }
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm*sqrt(LS),j);

  printf("v:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	norm += fabs( grid->R[i*NUM_VAR+2] - del_F[i*NUM_VAR+2] );
    }
  printf("  L1 norm of difference: %.15e\n",norm*sqrt(LS));
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+2] - del_F[i*NUM_VAR+2] );
	  norm += ( temp*temp );
	}
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm)*sqrt(LS));
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+2] - del_F[i*NUM_VAR+2] );
	  if ( temp > norm )
	    {
	      norm = temp;
	      j = i;
	    }
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm*sqrt(LS),j);
  

  printf("Pressure:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	norm += fabs( grid->R[i*NUM_VAR+3] - del_F[i*NUM_VAR+3] );
    }
  printf("  L1 norm of difference: %.15e\n",norm*sqrt(LS));
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+3] - del_F[i*NUM_VAR+3] );
	  norm += ( temp*temp );
	}
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm)*sqrt(LS));
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+3] - del_F[i*NUM_VAR+3] );
	  if ( temp > norm )
	    {
	      norm = temp;
	      j = i;
	    }
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm*sqrt(LS),j);
  */

  /*
   printf("weighted by CV area:\n");
  // Rho.
  printf("Density:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	norm += fabs( grid->R[i*NUM_VAR+0] - del_F[i*NUM_VAR+0] )*grid->cv_area[i];
    }
  printf("  L1 norm of difference: %.15e\n",norm);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+0] - del_F[i*NUM_VAR+0] )*grid->cv_area[i];
	  norm += ( temp*temp );
	}
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm));
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+0] - del_F[i*NUM_VAR+0] )*grid->cv_area[i];
	  if ( temp > norm )
	    {
	      norm = temp;
	      j = i;
	    }
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm,j);

  // u
  printf("U:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	norm += fabs( grid->R[i*NUM_VAR+1] - del_F[i*NUM_VAR+1] )*grid->cv_area[i];
    }
  printf("  L1 norm of difference: %.15e\n",norm);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+1] - del_F[i*NUM_VAR+1] )*grid->cv_area[i];
	  norm += ( temp*temp );
	}
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm));
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+1] - del_F[i*NUM_VAR+1] )*grid->cv_area[i];
	  if ( temp > norm )
	    {
	      norm = temp;
	      j = i;
	    }
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm,j);

  printf("v:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	norm += fabs( grid->R[i*NUM_VAR+2] - del_F[i*NUM_VAR+2] )*grid->cv_area[i];
    }
  printf("  L1 norm of difference: %.15e\n",norm);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+2] - del_F[i*NUM_VAR+2] )*grid->cv_area[i];
	  norm += ( temp*temp );
	}
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm));
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+2] - del_F[i*NUM_VAR+2] )*grid->cv_area[i];
	  if ( temp > norm )
	    {
	      norm = temp;
	      j = i;
	    }
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm,j);
  

  printf("Pressure:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	norm += fabs( grid->R[i*NUM_VAR+3] - del_F[i*NUM_VAR+3] )*grid->cv_area[i];
    }
  printf("  L1 norm of difference: %.15e\n",norm);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+3] - del_F[i*NUM_VAR+3] )*grid->cv_area[i];
	  norm += ( temp*temp );
	}
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm));
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( grid->R[i*NUM_VAR+3] - del_F[i*NUM_VAR+3] )*grid->cv_area[i];
	  if ( temp > norm )
	    {
	      norm = temp;
	      j = i;
	    }
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm,j);


  printf("relative error by del dot F:\n");
  // Rho.
  printf("Density:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	norm += fabs( (grid->R[i*NUM_VAR+0] - del_F[i*NUM_VAR+0]) / del_F[i*NUM_VAR+0] );
    }
  printf("  L1 norm of difference: %.15e\n",norm);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( (grid->R[i*NUM_VAR+0] - del_F[i*NUM_VAR+0])/del_F[i*NUM_VAR+0] );
	  norm += ( temp*temp );
	}
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm));
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( (grid->R[i*NUM_VAR+0] - del_F[i*NUM_VAR+0]) / del_F[i*NUM_VAR+0] );
	  if ( temp > norm )
	    {
	      norm = temp;
	      j = i;
	    }
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm,j);

  // u
  printf("U:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	norm += fabs( (grid->R[i*NUM_VAR+1] - del_F[i*NUM_VAR+1])/ del_F[i*NUM_VAR+1] );
    }
  printf("  L1 norm of difference: %.15e\n",norm);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( (grid->R[i*NUM_VAR+1] - del_F[i*NUM_VAR+1])/del_F[i*NUM_VAR+1] );
	  norm += ( temp*temp );
	}
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm));
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( (grid->R[i*NUM_VAR+1] - del_F[i*NUM_VAR+1])/del_F[i*NUM_VAR+1] );
	  if ( temp > norm )
	    {
	      norm = temp;
	      j = i;
	    }
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm,j);

  printf("v:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	norm += fabs( (grid->R[i*NUM_VAR+2] - del_F[i*NUM_VAR+2])/del_F[i*NUM_VAR+2]);
    }
  printf("  L1 norm of difference: %.15e\n",norm);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( (grid->R[i*NUM_VAR+2] - del_F[i*NUM_VAR+2])/del_F[i*NUM_VAR+2] );
	  norm += ( temp*temp );
	}
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm));
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( (grid->R[i*NUM_VAR+2] - del_F[i*NUM_VAR+2])/del_F[i*NUM_VAR+2] );
	  if ( temp > norm )
	    {
	      norm = temp;
	      j = i;
	    }
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm,j);
  

  printf("Pressure:\n");
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	norm += fabs( (grid->R[i*NUM_VAR+3] - del_F[i*NUM_VAR+3])/del_F[i*NUM_VAR+3] );
    }
  printf("  L1 norm of difference: %.15e\n",norm);
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( (grid->R[i*NUM_VAR+3] - del_F[i*NUM_VAR+3])/del_F[i*NUM_VAR+3] );
	  norm += ( temp*temp );
	}
    }
  printf("  L2 norm of difference: %.15e\n",sqrt(norm));
  
  norm = 0.;
  for ( i=1; i <= grid->nn; i++ )
    {
      if ( node_status[i] )
	{
	  temp = fabs( (grid->R[i*NUM_VAR+3] - del_F[i*NUM_VAR+3])/del_F[i*NUM_VAR+3] );
	  if ( temp > norm )
	    {
	      norm = temp;
	      j = i;
	    }
	}
    }
  printf("  L_inf norm of difference: %.15e   at node %d\n",norm,j);
  */


  freenull(node_status);
  freenull(del_F);

  return;
}


void flux_divergence ( double *x, double *q )
{
  double rho,u,v,P;
  double drdx, drdy, dudx, dudy;
  double dvdx, dvdy, dPdx, dPdy;
  double df1, df2, df3, df4;
  double gamma = 1.4;
  double pi = M_PI;

  df1 = 0.;
  df2 = 0.;
  df3 = 0.;
  df4 = 0.;

  rho = 1.0 + sin(pi*x[0])*sin(pi*x[1]);
  u = 0.5 + 0.5*sin(pi*x[0])*cos(2.*pi*x[1]);
  v = 0.5 + 0.5*cos(2.*pi*x[0])*sin(pi*x[1]);
  P = 1.0/gamma + 0.05*cos(2.*pi*x[0])*cos(2.*pi*x[1]);

  drdx = pi*cos(pi*x[0])*sin(pi*x[1]);
  drdy = pi*sin(pi*x[0])*cos(pi*x[1]);

  dudx =  0.5*pi*cos(pi*x[0])*cos(2.*pi*x[1]);
  dudy = -pi*sin(pi*x[0])*sin(2.*pi*x[1]);

  dvdx = -pi*sin(2.*pi*x[0])*sin(pi*x[1]);
  dvdy =  0.5*pi*cos(2.*pi*x[0])*cos(pi*x[1]);

  dPdx = -0.1*pi*sin(2.*pi*x[0])*cos(2.*pi*x[1]);
  dPdy = -0.1*pi*cos(2.*pi*x[0])*sin(2.*pi*x[1]);

  if ( ANNULUS )
    {
      double rhoi,r,r2,Mi,U,Ui,Ri;
      double gm1 = gamma - 1.0;
      double gm1inv = 1.0 / gm1;
      double gm1o2 = gm1 * 0.5;
      rhoi = 1.0;
      Mi = 2.0;
      Ui = 2.0;
      Ri = 2.0;

      r2 = (x[0]*x[0]) + (x[1]*x[1]);
      r = sqrt(r2);
      
      rho = pow( ( 1.0 + ((gamma - 1.0)/2.0)*Mi*Mi*( 1.0 - (Ri*Ri)/r2)) , (1.0/(gamma-1.0)) );
      rhoi = ( 1.0 + ((gamma - 1.0)/2.0)*Mi*Mi*( 1.0 - (Ri*Ri)/r2));
      U   = (Ui*Ri)/r;
      u   = (x[1]*U)/r;
      v   = (-x[0]*U)/r;
      P   = ( pow( rho , gamma) )/gamma;

      drdx = ( 2.0 * rho * gm1o2 * Mi*Mi * Ri*Ri * x[0] ) /
	     ( gm1 * r2*r2 * rhoi );

      drdy = ( 2.0 * rho * gm1o2 * Mi*Mi * Ri*Ri * x[1] ) /
	     ( gm1 * r2*r2 * rhoi );

      dudx = ( -4.0 * x[1] * Ri * x[0] ) / ( r2*r2 );
      dudy = ( 2.0*Ri )/r2 - ( 4.0*x[1]*x[1]*Ri )/(r2*r2);
      
      dvdx = ( -2.0*Ri )/r2 + ( 4.0*x[0]*x[0]*Ri )/(r2*r2);
      dvdy = ( 4.0 * x[1] * Ri * x[0] ) / ( r2*r2 );

      dPdx = ( 2.0 * ( pow( rho , gamma ) ) * gm1o2 * Mi*Mi * Ri*Ri * x[0] ) /
	     ( gm1 * r2*r2 * rhoi );

      dPdy = ( 2.0 * ( pow( rho , gamma ) ) * gm1o2 * Mi*Mi * Ri*Ri * x[1] ) /
	     ( gm1 * r2*r2 * rhoi );
    }

  
  df1 = rho*dudx + u*drdx + rho*dvdy + v*drdy;
  df2 = 2.*u*rho*dudx + u*u*drdx + dPdx + u*v*drdy + rho*u*dvdy + rho*v*dudy;
  df3 = u*v*drdx + rho*v*dudx + rho*u*dvdx + 2.*rho*v*dvdy + v*v*drdy + dPdy;
  df4 = gamma/(gamma-1.) * ( P*dudx + u*dPdx + P*dvdy + v*dPdy ) +
        0.5*( u*drdx*(u*u + v*v) + rho*dudx*(u*u + v*v) + 2.*rho*u*( u*dudx + v*dvdx) +
	      v*drdy*(u*u + v*v) + rho*dvdy*(u*u + v*v) + 2.*rho*v*( u*dudy + v*dvdy) ) ;


  q[0] = df1;
  q[1] = df2;
  q[2] = df3;
  q[3] = df4;

  return;
}


void MMS_Compute_Q_Exact ( GRID *grid )
{
  int n;                                 // Node loop counter.
  int i,j,k,e;                           // Loop counters.
  int v;                                 // Vector loop.
  int ind;                               // Component loop.
  int nodes[MAX_NUM_VERT];               // Element vertices.
  int node;                              // Node index.
  int gn;                                // Ghost node index.
  int num_vert;                          // Number of vertices for an element.
  int left_id;                           // ID of left edge.
  int right_id;                          // ID of right edge.
  double xc[NDIM];                       // Vectors.
  double xmidL[NDIM];
  double xmidR[NDIM];
  double v_left[NDIM];
  double v_right[NDIM];
  double v_middle[NDIM];
  double bary_coord[7][3];               // The barycentric coordinates.
  double wg[7];                          // The associated weights;
  double b_r = (6. - sqrt(15.))/21.;     // Precompute some the coordinates.
  double b_t = (6. + sqrt(15.))/21.;
  double b_s = (9. + 2.*sqrt(15.))/21.;
  double b_u = (9. - 2.*sqrt(15.))/21.;
  double wA = (155. - sqrt(15.))/1200.;
  double wB = (155. + sqrt(15.))/1200.;

  double XGP[NDIM];
  double QGP[NUM_VAR];
  double QTRI[NUM_VAR];

  double gp[7][3];                       // Gauss points for integration.
  double dA;                             // Differential area (for a triangle).

  // Curved triangle stuff.
  double ct_gp[7][2];                    // Gauss integration points in (r,s) space (xi,eta).
  double ct_nodes[7][2];                 // The x,y locations of the 6 nodes defining the p2 triangle (curved side).
  double N1,N2,N3,N4,N5,N6;              // Shape functions.
  double N1r,N1s,N2r,N2s,N3r,N3s;        // Shape function derivatives.
  double N4r,N4s,N5r,N5s,N6r,N6s;
  double xr,xs,yr,ys,jacobian;           // Transformation jacobian.
  double XI,ETA;                         // Triangle coordinates.
  int b,bc,seg,ct_flag=0;
  int ghost_node;

  int dummy_bcs[100];
  int reset_flag = 0;

  // This is very hacky and is for the annulus unstructured grid sequence only!
  if ( ANNULUS )
    {
      reset_flag = 0;

      
      if ( grid->bbc[2] != 14 && grid->bbc[2] != 13 )
	{
	  dummy_bcs[0] = grid->bbc[0];
	  dummy_bcs[1] = grid->bbc[1];
	  dummy_bcs[2] = grid->bbc[2];
	  dummy_bcs[3] = grid->bbc[3];
	  dummy_bcs[4] = grid->bbc[4];

	  grid->bbc[1] = 0;
	  grid->bbc[2] = 14;
	  grid->bbc[3] = 2;
	  grid->bbc[4] = 16;

	  reset_flag = 1;
	  
	  printf("RECALCULATING Control volume areas to account for curved boundaries while computing the exact solution.\n");

	  cv_calc_area_Green ( grid );

	}
      
      /*
      if ( grid->bbc[2] != 14 && grid->bbc[2] != 13 )
	{
	  dummy_bcs[0] = grid->bbc[0];
	  dummy_bcs[1] = grid->bbc[1];
	  dummy_bcs[2] = grid->bbc[2];
	  dummy_bcs[3] = grid->bbc[3];
	  dummy_bcs[4] = grid->bbc[4];
	  dummy_bcs[5] = grid->bbc[5];
	  dummy_bcs[6] = grid->bbc[6];

	  grid->bbc[1] = 2;
	  grid->bbc[2] = 13;
	  grid->bbc[3] = 0;
	  grid->bbc[4] = 2;
	  grid->bbc[5] = 15;
	  grid->bbc[6] = 0;
	  
	  reset_flag = 1;
	  
	  printf("RECALCULATING Control volume areas to account for curved boundaries while computing the exact solution.\n");

	  cv_calc_area_Green ( grid );
	}
      */

    }

  
  // Clean out memory because we are accumulating.
  for ( n=1; n <= grid->nn; n++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	grid->Qe[n*NUM_VAR+j] = 0.;
    }
  
  // Now we calculate the area-averaged value of the initial condition using Gaussian quadrature
  
  // Set up the barycentric coordinates and their respective weights.
  bary_coord[0][0] = 1./3.; bary_coord[0][1] = 1./3.; bary_coord[0][2] = 1./3.;
  bary_coord[1][0] = b_r;   bary_coord[1][1] = b_r;   bary_coord[1][2] = b_s;
  bary_coord[2][0] = b_r;   bary_coord[2][1] = b_s;   bary_coord[2][2] = b_r;
  bary_coord[3][0] = b_s;   bary_coord[3][1] = b_r;   bary_coord[3][2] = b_r;
  bary_coord[4][0] = b_t;   bary_coord[4][1] = b_t;   bary_coord[4][2] = b_u;
  bary_coord[5][0] = b_t;   bary_coord[5][1] = b_u;   bary_coord[5][2] = b_t;
  bary_coord[6][0] = b_u;   bary_coord[6][1] = b_t;   bary_coord[6][2] = b_t;
      
  wg[0] = 9./40.;
  wg[1] = wA;
  wg[2] = wA;
  wg[3] = wA;
  wg[4] = wB;
  wg[5] = wB;
  wg[6] = wB;
  
  // Lets also map these guys onto the reference triangle for integration on the triangles with curved sides.
  for ( j=0; j < 7; j++ )
    {
      ct_gp[j][0] = (bary_coord[j][0] * 0. +
		     bary_coord[j][1] * 1. +
		     bary_coord[j][2] * 0. );
      
      ct_gp[j][1] = (bary_coord[j][0] * 0. +
		     bary_coord[j][1] * 0. +
		     bary_coord[j][2] * 1. );
      
    }  // The weights can be reused since everything is consistent in the ordering.
  
  // Start the process by looping over all the elements and looping through the internal median dual pieces.
  for ( e=Tri; e <= Quad; e++ )                 // Element types loop.
    {
      // Find the number of nodes and edges for this element type.
      num_vert = NumberNodesForElement(e);

      for ( j=1; j <= grid->num_elem[e]; j++ )              // Element loop.
	{
	  // Retrieve the element node ids. This is for convenience.
	  for ( k=0; k < num_vert; k++ )              // Element node extraction loop.
	    {
	      nodes[k] = grid->c2n[e][j*num_vert+k];
	    }                                         // End element node extraction loop.
	  
	  // Find the element centroid.
	  xc[0] = 0.; xc[1] = 0.;
	  for ( k=0; k < num_vert; k++ )
	    {
	      for ( i=0; i < NDIM; i++ )
		{
		  xc[i] += grid->x[NDIM*nodes[k]+i];
		}
	    }
	  xc[0] /= ((double)num_vert);
	  xc[1] /= ((double)num_vert);

	  // Loop through each node in the element.
	  for ( k=0; k < num_vert; k++ )              // Element node loop.
	    {
	      node = nodes[k];

	      // Get the two neighboring vertices.
	      left_id = (k+1)%num_vert;
	      right_id = (k+num_vert-1)%num_vert;
		  
	      left_id = grid->c2n[e][j*num_vert+left_id];
	      right_id = grid->c2n[e][j*num_vert+right_id];
	      
	      // Get the edge midpoints.
	      for ( i=0; i < NDIM; i++ )
		{
		  xmidL[i] = grid->x[NDIM*node+i] + grid->x[NDIM*left_id+i];
		  xmidR[i] = grid->x[NDIM*node+i] + grid->x[NDIM*right_id+i];
		}

	      for ( i=0; i < NDIM; i++ )
		{
		  xmidL[i] *= 0.5;
		  xmidR[i] *= 0.5;
		}

	      // Make the vectors.
	      for ( i=0; i < NDIM; i++ )
		{
		  v_left[i]   = xmidL[i] - grid->x[NDIM*node+i];
		  v_right[i]  = xmidR[i] - grid->x[NDIM*node+i];
		  v_middle[i] = xc[i]    - grid->x[NDIM*node+i];
		}

	      ct_flag = 0;
	      
	      // If the node and its neighbor are both boundary nodes, we need to gather some information.
	      if ( grid->node_state[node] == BOUNDARY && grid->node_state[left_id] == BOUNDARY )
		{
		  // Retrieve the relevant information.
		  
		  // First find the correct bedge. It will have both nodes.
		  for ( i=1; i <= grid->nbedges; i++ )
		    {
		      if ( node == grid->bedges[i*5] )  // candidate bedge.
			{
			  b = grid->bedges[i*5+3];
			  seg = grid->bedges[i*5+4];
			  
			  if ( (  node == grid->bs[b][seg][0]   ||  node == grid->bs[b][seg][1]   ) &&
			       ( left_id == grid->bs[b][seg][0] || left_id == grid->bs[b][seg][1] )   )  // we have a winner.
			    {
			      ghost_node = grid->bedges[i*5+1] + grid->nn;
			      b = grid->bedges[i*5+3];
			      bc = grid->bbc[b];
			      seg = grid->bedges[i*5+4];
			      
			      if ( bc > 10 )
				ct_flag = 1;
			      
			      break;
			    }
			}
		    }
		}

	      
	      if ( ct_flag )
		{
		  // Now we need to set up the six nodes used by the shape functions.
		      
		  // Node 1 is the element centroid.
		  ct_nodes[1][0] = xc[0];
		  ct_nodes[1][1] = xc[1];

		  // Node 2 is the node we are working on.
		  ct_nodes[2][0] = grid->x[NDIM*node+0];
		  ct_nodes[2][1] = grid->x[NDIM*node+1]; 

		  // Third node is the ghost node.
		  //ct_nodes[3][0] = grid->x[NDIM*ghost_node+0];
		  //ct_nodes[3][1] = grid->x[NDIM*ghost_node+1];  // this might be incorrect when the edges were created from linear boundaries.

		  curved_boundary_midpoint(grid, bc, &(ct_nodes[2][0]), &(grid->x[left_id*NDIM+0]), &(ct_nodes[3][0]) );

		  // Node 4 is the midpoint of 1,2.
		  ct_nodes[4][0] = 0.5*(ct_nodes[1][0] + ct_nodes[2][0]);
		  ct_nodes[4][1] = 0.5*(ct_nodes[1][1] + ct_nodes[2][1]);

		  // Node 5 is the midpoint of the curved boundary between node and ghost node.
		  curved_boundary_midpoint(grid, bc, &(ct_nodes[2][0]), &(ct_nodes[3][0]), &(ct_nodes[5][0]) );

		  // Node 6 is the midpoint of 1,3.
		  ct_nodes[6][0] = 0.5*(ct_nodes[1][0] + ct_nodes[3][0]);
		  ct_nodes[6][1] = 0.5*(ct_nodes[1][1] + ct_nodes[3][1]);

		  // First we need to make sure that the concavity constraint is met. Under a linear transformation
		  // to the reference element, node 5 must be greater than 0.25 in both r and s coordinates.
		  // See Flaherty's notes online or Claes Johnson's book.

		  XI = ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
			 (ct_nodes[5][0]*ct_nodes[3][1] - ct_nodes[5][1]*ct_nodes[3][0]) +
			 (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
		       ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
		         (ct_nodes[2][0]*ct_nodes[3][1] - ct_nodes[2][1]*ct_nodes[3][0]) +
		         (ct_nodes[2][0]*ct_nodes[1][1] - ct_nodes[2][1]*ct_nodes[1][0])  );
		  
		  ETA = ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			  (ct_nodes[5][0]*ct_nodes[2][1] - ct_nodes[5][1]*ct_nodes[2][0]) +
			  (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
		        ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			  (ct_nodes[3][0]*ct_nodes[2][1] - ct_nodes[3][1]*ct_nodes[2][0]) +
			  (ct_nodes[3][0]*ct_nodes[1][1] - ct_nodes[1][0]*ct_nodes[3][1])  );

		  if ( ETA < 0.25 || XI < 0.25 )  // The transformation will not be invertible.
		    {
		      printf("CRITICAL ERROR: Triangle 1 of element %d of type %d is too concave (compute cv averages: hessian.C)!\n",j,e);
		      printf("  xi = %.15e     eta = %.15e\n",XI,ETA);
		      fflush(stdout);
		      exit(1);
		    }

		  // Do this for each variable.

		  QTRI[0] = 0.; QTRI[1] = 0.; QTRI[2] = 0.; QTRI[3] = 0.;
		  
		  // Now I can loop over the Gaussian integration nodes and apply the quadratic transformation.
		  for ( i=0; i < 7; i++ )
		    {
		      XI = ct_gp[i][0];
		      ETA = ct_gp[i][1];
		      
		      // Compute basis functions.
		      N1 = 1. - 3.*XI - 3.*ETA + 4.*XI*ETA + 2.*XI*XI + 2.*ETA*ETA;
		      N2 = XI*(2.*XI - 1.);
		      N3 = ETA*(2.*ETA - 1.);
		      N4 = 4.*XI*(1. - XI - ETA);
		      N5 = 4.*XI*ETA;
		      N6 = 4.*ETA*(1. - XI - ETA);

		      // Compute their derivatives.
		      N1r = -3. + 4.*ETA + 4.*XI;
		      N1s = -3. + 4.*XI + 4.*ETA;
		      N2r = 4.*XI - 1.;
		      N2s = 0.;
		      N3r = 0.;
		      N3s = 4.*ETA - 1.;
		      N4r = 4.*(1. - 2.*XI - ETA);
		      N4s = -4.*XI;
		      N5r = 4.*ETA;
		      N5s = 4.*XI;
		      N6r = -4.*ETA;
		      N6s = 4.*(1. - 2.*ETA - XI);

		      // Compute the derivatives of x,y wrt to the transformation.
		      xr = N1r*ct_nodes[1][0] + N2r*ct_nodes[2][0] + N3r*ct_nodes[3][0] + N4r*ct_nodes[4][0]
			 + N5r*ct_nodes[5][0] + N6r*ct_nodes[6][0];
		      xs = N1s*ct_nodes[1][0] + N2s*ct_nodes[2][0] + N3s*ct_nodes[3][0] + N4s*ct_nodes[4][0]
			 + N5s*ct_nodes[5][0] + N6s*ct_nodes[6][0];
		      yr = N1r*ct_nodes[1][1] + N2r*ct_nodes[2][1] + N3r*ct_nodes[3][1] + N4r*ct_nodes[4][1]
			 + N5r*ct_nodes[5][1] + N6r*ct_nodes[6][1];
		      ys = N1s*ct_nodes[1][1] + N2s*ct_nodes[2][1] + N3s*ct_nodes[3][1] + N4s*ct_nodes[4][1]
			 + N5s*ct_nodes[5][1] + N6s*ct_nodes[6][1];

		      jacobian = xr*ys - xs*yr;

		      // Now I need to extrapolate to this point. First we get the real coordinates of the
		      // Gaussian point.
		      XGP[0] = N1*ct_nodes[1][0] + N2*ct_nodes[2][0] + N3*ct_nodes[3][0] +
			       N4*ct_nodes[4][0] + N5*ct_nodes[5][0] + N6*ct_nodes[6][0];

		      XGP[1] = N1*ct_nodes[1][1] + N2*ct_nodes[2][1] + N3*ct_nodes[3][1] +
			       N4*ct_nodes[4][1] + N5*ct_nodes[5][1] + N6*ct_nodes[6][1];

		      MMS_Get_Point_Qe ( XGP, QGP );

		      // Accumulate.
		      for ( v=0; v < NUM_VAR; v++ )
			{
			  QTRI[v] += ( wg[i] * QGP[v] *jacobian * 0.5);
			} 
		      
		    }

		  // Add to the running total for the node.
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      grid->Qe[node*NUM_VAR+v] += QTRI[v];
		    }
		  
		}
	      else
		{
		  // Get the area.
		  dA = 0.5*( v_left[0]*v_middle[1] - v_left[1]*v_middle[0] );
		  
		  // Initialize the Gauss points.
		  for ( v=0; v < 7; v++ )
		    {
		      for ( ind=0; ind < NDIM; ind++ )
			{
			  gp[v][ind] = (bary_coord[v][0] * grid->x[node*NDIM+ind] +
					bary_coord[v][1] * xmidL[ind] +
					bary_coord[v][2] * xc[ind]);
			}
		      gp[v][2] = 0.;   // Always 0 since 2D vectors.
		    }

		  // Now we loop over the Gauss points and get the function value
		  // and perform the numerical integration (the area integral of the triangle).
		  
		  // Do this for each variable.
		  
		  QTRI[0] = 0.; QTRI[1] = 0.; QTRI[2] = 0.; QTRI[3] = 0.;
		  
		  for ( i=0; i < 7; i++ )
		    {
		      // Get the value at the point.
		      XGP[0] = gp[i][0];
		      XGP[1] = gp[i][1];
		      
		      MMS_Get_Point_Qe ( XGP , QGP );
		  
		      for ( v=0; v < NUM_VAR; v++ )
			{
			  QTRI[v] += ( wg[i] * QGP[v] );
			} 
		    }
	      
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      QTRI[v] *= dA;
		    }
	      
		  // Accumulate the integral to the node.
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      grid->Qe[node*NUM_VAR+v] += QTRI[v];
		    }
		}

	      // Process the second triangle.

	      ct_flag = 0;
	      
	      // If the node and its neighbor are both boundary nodes, we need to gather some information.
	      if ( grid->node_state[node] == BOUNDARY && grid->node_state[right_id] == BOUNDARY )
		{
		  // Retrieve the relevant information.
		  
		  // First find the correct bedge. It will have both nodes.
		  for ( i=1; i <= grid->nbedges; i++ )
		    {
		      if ( node == grid->bedges[i*5] )  // candidate bedge.
			{
			  b = grid->bedges[i*5+3];
			  seg = grid->bedges[i*5+4];
			  
			  if ( (   node == grid->bs[b][seg][0]   ||   node == grid->bs[b][seg][1]   ) &&
			       ( right_id == grid->bs[b][seg][0] || right_id == grid->bs[b][seg][1] )   )  // we have a winner.
			    {
			      ghost_node = grid->bedges[i*5+1] + grid->nn;
			      b = grid->bedges[i*5+3];
			      bc = grid->bbc[b];
			      seg = grid->bedges[i*5+4];
			      if ( bc > 10 )
				ct_flag = 1;
			      
			      break;
			    }
			}
		    }
		}

	      if ( ct_flag )
		{
		  // Now we need to set up the six nodes used by the shape functions.
		      
		  // Node 1 is the element centroid.
		  ct_nodes[1][0] = xc[0];
		  ct_nodes[1][1] = xc[1];

		  // Node 2 is the ghost_node.
		  //ct_nodes[2][0] = grid->x[NDIM*ghost_node+0];
		  //ct_nodes[2][1] = grid->x[NDIM*ghost_node+1];  // this might be incorrect when the edges were created from linear boundaries.

		  curved_boundary_midpoint(grid, bc, &(grid->x[right_id*NDIM+0]), &(grid->x[node*NDIM+0]), &(ct_nodes[2][0]) );

		  // Third node is the node we are working on.
		  ct_nodes[3][0] = grid->x[NDIM*node+0];
		  ct_nodes[3][1] = grid->x[NDIM*node+1];

		  // Node 4 is the midpoint of 1,2.
		  ct_nodes[4][0] = 0.5*(ct_nodes[1][0] + ct_nodes[2][0]);
		  ct_nodes[4][1] = 0.5*(ct_nodes[1][1] + ct_nodes[2][1]);

		  // Node 5 is the midpoint of the curved boundary between node and ghost node.
		  curved_boundary_midpoint(grid, bc, &(ct_nodes[2][0]), &(ct_nodes[3][0]), &(ct_nodes[5][0]) );

		  // Node 6 is the midpoint of 1,3.
		  ct_nodes[6][0] = 0.5*(ct_nodes[1][0] + ct_nodes[3][0]);
		  ct_nodes[6][1] = 0.5*(ct_nodes[1][1] + ct_nodes[3][1]);

		  // First we need to make sure that the concavity constraint is met. Under a linear transformation
		  // to the reference element, node 5 must be greater than 0.25 in both r and s coordinates.
		  // See Flaherty's notes online or Claes Johnson's book.

		  XI = ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
			 (ct_nodes[5][0]*ct_nodes[3][1] - ct_nodes[5][1]*ct_nodes[3][0]) +
			 (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
		       ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
		         (ct_nodes[2][0]*ct_nodes[3][1] - ct_nodes[2][1]*ct_nodes[3][0]) +
		         (ct_nodes[2][0]*ct_nodes[1][1] - ct_nodes[2][1]*ct_nodes[1][0])  );
		  
		  ETA = ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			  (ct_nodes[5][0]*ct_nodes[2][1] - ct_nodes[5][1]*ct_nodes[2][0]) +
			  (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
		        ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			  (ct_nodes[3][0]*ct_nodes[2][1] - ct_nodes[3][1]*ct_nodes[2][0]) +
			  (ct_nodes[3][0]*ct_nodes[1][1] - ct_nodes[1][0]*ct_nodes[3][1])  );

		  if ( ETA < 0.25 || XI < 0.25 )  // The transformation will not be invertible.
		    {
		      printf("CRITICAL ERROR: Triangle 1 of element %d of type %d is too concave (compute cv averages: hessian.C)!\n",j,e);
		      printf("  xi = %.15e     eta = %.15e\n",XI,ETA);
		      fflush(stdout);
		      exit(1);
		    }

		  // Do this for each variable.

		  QTRI[0] = 0.; QTRI[1] = 0.; QTRI[2] = 0.; QTRI[3] = 0.;
		  
		  // Now I can loop over the Gaussian integration nodes and apply the quadratic transformation.
		  for ( i=0; i < 7; i++ )
		    {
		      XI = ct_gp[i][0];
		      ETA = ct_gp[i][1];
		      
		      // Compute basis functions.
		      N1 = 1. - 3.*XI - 3.*ETA + 4.*XI*ETA + 2.*XI*XI + 2.*ETA*ETA;
		      N2 = XI*(2.*XI - 1.);
		      N3 = ETA*(2.*ETA - 1.);
		      N4 = 4.*XI*(1. - XI - ETA);
		      N5 = 4.*XI*ETA;
		      N6 = 4.*ETA*(1. - XI - ETA);

		      // Compute their derivatives.
		      N1r = -3. + 4.*ETA + 4.*XI;
		      N1s = -3. + 4.*XI + 4.*ETA;
		      N2r = 4.*XI - 1.;
		      N2s = 0.;
		      N3r = 0.;
		      N3s = 4.*ETA - 1.;
		      N4r = 4.*(1. - 2.*XI - ETA);
		      N4s = -4.*XI;
		      N5r = 4.*ETA;
		      N5s = 4.*XI;
		      N6r = -4.*ETA;
		      N6s = 4.*(1. - 2.*ETA - XI);

		      // Compute the derivatives of x,y wrt to the transformation.
		      xr = N1r*ct_nodes[1][0] + N2r*ct_nodes[2][0] + N3r*ct_nodes[3][0] + N4r*ct_nodes[4][0]
			 + N5r*ct_nodes[5][0] + N6r*ct_nodes[6][0];
		      xs = N1s*ct_nodes[1][0] + N2s*ct_nodes[2][0] + N3s*ct_nodes[3][0] + N4s*ct_nodes[4][0]
			 + N5s*ct_nodes[5][0] + N6s*ct_nodes[6][0];
		      yr = N1r*ct_nodes[1][1] + N2r*ct_nodes[2][1] + N3r*ct_nodes[3][1] + N4r*ct_nodes[4][1]
			 + N5r*ct_nodes[5][1] + N6r*ct_nodes[6][1];
		      ys = N1s*ct_nodes[1][1] + N2s*ct_nodes[2][1] + N3s*ct_nodes[3][1] + N4s*ct_nodes[4][1]
			 + N5s*ct_nodes[5][1] + N6s*ct_nodes[6][1];

		      jacobian = xr*ys - xs*yr;

		      // Now I need to extrapolate to this point. First we get the real coordinates of the
		      // Gaussian point.
		      XGP[0] = N1*ct_nodes[1][0] + N2*ct_nodes[2][0] + N3*ct_nodes[3][0] +
			       N4*ct_nodes[4][0] + N5*ct_nodes[5][0] + N6*ct_nodes[6][0];

		      XGP[1] = N1*ct_nodes[1][1] + N2*ct_nodes[2][1] + N3*ct_nodes[3][1] +
			       N4*ct_nodes[4][1] + N5*ct_nodes[5][1] + N6*ct_nodes[6][1];

		      MMS_Get_Point_Qe ( XGP, QGP );
		      
		      // Accumulate.
		      for ( v=0; v < NUM_VAR; v++ )
			{
			  QTRI[v] += ( wg[i] * QGP[v] *jacobian * 0.5);
			} 
		      
		    }

		  // Add to the running total for the node.
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      grid->Qe[node*NUM_VAR+v] += QTRI[v];
		    }
		  
		}
	      else
		{
		  // Get the area.
		  dA = 0.5*( v_middle[0]*v_right[1] - v_middle[1]*v_right[0] );
		  
		  // Initialize the Gauss points.
		  for ( v=0; v < 7; v++ )
		    {
		      for ( ind=0; ind < NDIM; ind++ )
			{
			  gp[v][ind] = (bary_coord[v][0] * grid->x[node*NDIM+ind] +
					bary_coord[v][1] * xmidR[ind] +
					bary_coord[v][2] * xc[ind]);
			}
		      gp[v][2] = 0.;   // Always 0 since 2D vectors.
		    }
	      
		  // Now we loop over the Gauss points and get the function value
		  // and perform the numerical integration (the area integral of the triangle).
	      
		  // Do this for each variable.
	      
		  QTRI[0] = 0.; QTRI[1] = 0.; QTRI[2] = 0.; QTRI[3] = 0.;
	      
		  for ( i=0; i < 7; i++ )
		    {
		      // Get the value at the point.
		      XGP[0] = gp[i][0];
		      XGP[1] = gp[i][1];
		  
		      MMS_Get_Point_Qe ( XGP , QGP );
		      
		      for ( v=0; v < NUM_VAR; v++ )
			{
			  QTRI[v] += ( wg[i] * QGP[v] );
			} 
		    }
		  
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      QTRI[v] *= dA;
		    }
		  
		  // Accumulate the integral to the node.
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      grid->Qe[node*NUM_VAR+v] += QTRI[v];
		    }
		}
	    }                                       // End element node loop.
	  
	}                                           // End element loop.
      
    }                                               // End element type loop.
  
  
  // Finished the accumulation. Now divide by the control volume area.
  for ( i=1; i<= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  grid->Qe[i*NUM_VAR+j] /= grid->cv_area[i];
	}
    }

  //if ( MMS || ANNULUS )
  if ( MMS )// if I am running a manufactured solution case, I should copy the exact solution over as my initial guess.
    {                   // The annulus case appearently needs the exact solution as an initial guess since it diverges when I give it
                        // a standard initial guess.
      for ( i=1; i <= grid->nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      grid->Q[i*NUM_VAR+j] = grid->Qe[i*NUM_VAR+j];
	    }
	}
      
      // if I am reconstructing the primitive variables, then I need to convert grid->Q to conserved since it just stored primitives.
      if ( RECON_PRIM )
	{
	  for ( i=1; i <= grid->nn; i++ )
	    {
	      ConvertPrimitiveToConserved ( 1.40 , &(grid->Q[i*NUM_VAR]) );              // Hard set gamma here. 
	    }
	}

      
      for ( i=1; i <= grid->nbedges; i++ )
	{
	  // Grab the data from the structure.
	  n =  grid->bedges[i*5+0];
	  gn = grid->bedges[i*5+1] + grid->nn;

	  // For boundary nodes, we copy over the physical node values.
	  grid->Q[gn*NUM_VAR+0] = grid->Q[n*NUM_VAR+0];
	  grid->Q[gn*NUM_VAR+1] = grid->Q[n*NUM_VAR+1];
	  grid->Q[gn*NUM_VAR+2] = grid->Q[n*NUM_VAR+2];
	  grid->Q[gn*NUM_VAR+3] = grid->Q[n*NUM_VAR+3];
	}
    }
 
  // Copy the boundary node values to its ghost node.
  for ( i=1; i <= grid->nbedges; i++ )
    {
      // Grab the data from the structure.
      n =  grid->bedges[i*5+0];
      gn = grid->bedges[i*5+1] + grid->nn;

      // Copy the variables over.
      grid->Qe[gn*NUM_VAR+0] = grid->Qe[n*NUM_VAR+0];
      grid->Qe[gn*NUM_VAR+1] = grid->Qe[n*NUM_VAR+1];
      grid->Qe[gn*NUM_VAR+2] = grid->Qe[n*NUM_VAR+2];
      grid->Qe[gn*NUM_VAR+3] = grid->Qe[n*NUM_VAR+3];
    }

  if ( ANNULUS && reset_flag )
    {
      grid->bbc[1] = dummy_bcs[1];
      grid->bbc[2] = dummy_bcs[2];
      grid->bbc[3] = dummy_bcs[3];
      grid->bbc[4] = dummy_bcs[4];
      //grid->bbc[5] = dummy_bcs[5];
      //grid->bbc[6] = dummy_bcs[6];

      printf("RECALCULATING Control volume areas with the original boundary conditions (linear).\n");

      cv_calc_area_Green ( grid );

    }


  return;
}



void MMS_Get_Point_Qe ( double *x, double *q )
{
  double pi = M_PI;

  double rho0,u0,v0,E0;
  double rho,u,v,E,P;

  //rho0 = 1.;
  //u0 = 0.5;
  //v0 = 0.5;
  //E0 = 3.;

  //rho = rho0 * ( 1. + sin(pi*x[0]) * cos(pi*x[0]) * sin(pi*x[1]) * cos(pi*x[1]) );
  //u = u0 * ( 1. + sin(2.*pi*x[0]) * cos(2.*pi*x[0]) * sin(2.*pi*x[1]) * cos(2.*pi*x[1]) );
  //v = v0 * (1. + sin(2.*pi*x[0]) * cos(2.*pi*x[0]) * sin(2.*pi*x[1]) * cos(2.*pi*x[1]) );
  //E = E0 * (1. + sin(pi*x[0])*sin(pi*x[0])*sin(pi*x[1])*sin(pi*x[1]) );

  //rho = 1.0 + sin(pi*x[0])*sin(pi*x[1]);
  //u = 0.25 + 0.25*sin(pi*x[0])*cos(2.*pi*x[1]);
  //v = 0.25 + 0.25*cos(2.*pi*x[0])*sin(pi*x[1]);
  //P = 1.0/1.4 + 0.05*cos(2.*pi*x[0])*cos(2.*pi*x[1]);

  rho = 1.0 + 0.25*sin(pi*x[0])*sin(pi*x[1]);
  u = 0.25 + 0.25*sin(pi*x[0])*cos(2.*pi*x[1]);
  v = 0.25 + 0.25*cos(2.*pi*x[0])*sin(pi*x[1]);
  P = 1.0/1.4 + 0.05*cos(2.*pi*x[0])*cos(2.*pi*x[1]);

  E = P / (1.4 - 1.) + 0.5*rho*(u*u + v*v);

  if ( MMS_EXP )
    {
      rho = exp(x[0]-1.0) * exp(x[1] - 1.0);
      u   = exp(x[0]-1.0) * exp(x[1] - 1.0);
      v   = exp(x[0]-1.0) * exp(x[1] - 1.0);
      E   = 1.5 * exp( 3.0 * ( x[0] - 1.0) ) * exp( 3.0 * ( x[1] - 1.0 ) );
      P   = ( (1.4-1.0) * 0.5 ) * exp( 3.0 * ( x[0] - 1.0) ) * exp( 3.0 * ( x[1] - 1.0 ) );
    }
  
  //if ( ANNULUS || 1 )
  if ( ANNULUS && 0 )
    {
      double rhoi,r,r2,Mi,U,Ui,Ri;
      rhoi = 1.0;
      Mi = 2.0;
      Ui = 2.0;
      Ri = 2.0;
      
      r2 = (x[0]*x[0]) + (x[1]*x[1]);
      r = sqrt(r2);
      
      rho = pow( ( 1.0 + ((1.4 - 1.0)/2.0)*Mi*Mi*( 1.0 - (Ri*Ri)/r2)) , (1.0/(1.4-1.0)) );
      U   = (Ui*Ri)/r;
      u   = (x[1]*U)/r;
      v   = (-x[0]*U)/r;
      P   = ( pow( rho , 1.4) )/1.4;
      
      E   = P / (1.4 - 1.) + 0.5*rho*(u*u + v*v);
    }

  // Store the variables.
  if ( RECON_PRIM )
    {
      q[0] = rho;
      q[1] = u;
      q[2] = v;
      q[3] = P;
    }
  else
    {
      q[0] = rho;
      q[1] = rho*u;
      q[2] = rho*v;
      q[3] = E;
    }

  return;
}


void MMS_Get_Flux_Divergence ( double *x, double *q )
{
  double pi = M_PI;
  double gamma = 1.4;
  double rho0,u0,v0,E0;
  double rho,u,v,E,P;
  double df1, df2, df3, df4;
  double drdx,drdy,dudx,dudy;
  double dvdx,dvdy,dEdx,dEdy,dPdx,dPdy;

  rho0 = 1.;
  u0 = 0.5;
  v0 = 0.5;
  E0 = 3.;

  /*
  rho = rho0 * ( 1. + sin(pi*x[0]) * cos(pi*x[0]) * sin(pi*x[1]) * cos(pi*x[1]) );
  u = u0 * ( 1. + sin(2.*pi*x[0]) * cos(2.*pi*x[0]) * sin(2.*pi*x[1]) * cos(2.*pi*x[1]) );
  v = v0 * (1. + sin(2.*pi*x[0]) * cos(2.*pi*x[0]) * sin(2.*pi*x[1]) * cos(2.*pi*x[1]) );
  E = E0 * (1. + sin(pi*x[0])*sin(pi*x[0])*sin(pi*x[1])*sin(pi*x[1]) );
  P = (gamma - 1.)*(E - 0.5*rho*(u*u + v*v));

  drdx = rho0 * pi * ( cos(pi*x[0])*cos(pi*x[0])*sin(pi*x[1])*cos(pi*x[1]) -
		       sin(pi*x[0])*sin(pi*x[0])*sin(pi*x[1])*cos(pi*x[1]) );

  drdy = rho0 * pi * ( cos(pi*x[1])*cos(pi*x[1])*sin(pi*x[0])*cos(pi*x[0]) -
		       sin(pi*x[1])*sin(pi*x[1])*sin(pi*x[0])*cos(pi*x[0]) );

  dudx = u0 * 2. * pi * ( sin(2.*pi*x[1])*cos(2.*pi*x[1]) * ( cos(2.*pi*x[0])*cos(2.*pi*x[0]) -
	 						      sin(2.*pi*x[0])*sin(2.*pi*x[0])));

  dudy = u0 * 2. * pi * ( sin(2.*pi*x[0])*cos(2.*pi*x[0]) * ( cos(2.*pi*x[1])*cos(2.*pi*x[1]) -
	 						      sin(2.*pi*x[1])*sin(2.*pi*x[1])));

  dvdx = v0 * 2. * pi * ( sin(2.*pi*x[1])*cos(2.*pi*x[1]) * ( cos(2.*pi*x[0])*cos(2.*pi*x[0]) -
	 						      sin(2.*pi*x[0])*sin(2.*pi*x[0])));

  dvdy = v0 * 2. * pi * ( sin(2.*pi*x[0])*cos(2.*pi*x[0]) * ( cos(2.*pi*x[1])*cos(2.*pi*x[1]) -
							      sin(2.*pi*x[1])*sin(2.*pi*x[1])));

  dEdx = E0 * 2. * pi * sin(pi*x[0])*cos(pi*x[0])*sin(pi*x[1])*sin(pi*x[1]);

  dEdy = E0 * 2. * pi * sin(pi*x[0])*sin(pi*x[0])*sin(pi*x[1])*cos(pi*x[1]);

  dPdx = (gamma - 1.)*( dEdx - 0.5*drdx*(u*u + v*v) - rho*(u*dudx + v*dvdx));

  dPdy = (gamma - 1.)*( dEdy - 0.5*drdy*(u*u + v*v) - rho*(u*dudy + v*dvdy));

  q[0] = rho*dudx + drdx*u + rho*dvdy + drdy*v;
  q[1] = drdx*u*u + 2.*rho*u*dudx + dPdx + drdy*u*v + rho*dudy*v + rho*u*dvdy;
  q[2] = drdx*u*v + rho*dudx*v + rho*u*dvdx + drdy*v*v + 2.*rho*v*dvdy + dPdy;
  q[3] = (E+P)*dudx + (dEdx+dPdx)*u + (E+P)*dvdy + (dEdy+dPdy)*v;
  */

  df1 = 0.;
  df2 = 0.;
  df3 = 0.;
  df4 = 0.;

  //rho = 1.0 + sin(pi*x[0])*sin(pi*x[1]);
  rho = 1.0 + 0.25*sin(pi*x[0])*sin(pi*x[1]);
  u = 0.25 + 0.25*sin(pi*x[0])*cos(2.*pi*x[1]);
  v = 0.25 + 0.25*cos(2.*pi*x[0])*sin(pi*x[1]);
  P = 1.0/gamma + 0.05*cos(2.*pi*x[0])*cos(2.*pi*x[1]);

  //drdx = pi*cos(pi*x[0])*sin(pi*x[1]);
  //drdy = pi*sin(pi*x[0])*cos(pi*x[1]);
  drdx = 0.25*pi*cos(pi*x[0])*sin(pi*x[1]);
  drdy = 0.25*pi*sin(pi*x[0])*cos(pi*x[1]);

  dudx =  0.25*pi*cos(pi*x[0])*cos(2.*pi*x[1]);
  dudy = -0.5*pi*sin(pi*x[0])*sin(2.*pi*x[1]);

  dvdx = -0.5*pi*sin(2.*pi*x[0])*sin(pi*x[1]);
  dvdy =  0.25*pi*cos(2.*pi*x[0])*cos(pi*x[1]);

  dPdx = -0.1*pi*sin(2.*pi*x[0])*cos(2.*pi*x[1]);
  dPdy = -0.1*pi*cos(2.*pi*x[0])*sin(2.*pi*x[1]);

  if ( ANNULUS && 0 )
    {
      double rhoi,r,r2,Mi,U,Ui,Ri;
      double gm1 = gamma - 1.0;
      double gm1inv = 1.0 / gm1;
      double gm1o2 = gm1 * 0.5;
      rhoi = 1.0;
      Mi = 2.0;
      Ui = 2.0;
      Ri = 2.0;

      r2 = (x[0]*x[0]) + (x[1]*x[1]);
      r = sqrt(r2);

      rho = pow( ( 1.0 + ((gamma - 1.0)/2.0)*Mi*Mi*( 1.0 - (Ri*Ri)/r2)) , (1.0/(gamma-1.0)) );
      rhoi = ( 1.0 + ((gamma - 1.0)/2.0)*Mi*Mi*( 1.0 - (Ri*Ri)/r2));
      U   = (Ui*Ri)/r;
      u   = (x[1]*U)/r;
      v   = (-x[0]*U)/r;
      P   = ( pow( rho , gamma) )/gamma;

      drdx = ( 2.0 * rho * gm1o2 * Mi*Mi * Ri*Ri * x[0] ) /
	     ( gm1 * r2*r2 * rhoi );

      drdy = ( 2.0 * rho * gm1o2 * Mi*Mi * Ri*Ri * x[1] ) /
	     ( gm1 * r2*r2 * rhoi );

      dudx = ( -4.0 * x[1] * Ri * x[0] ) / ( r2*r2 );
      dudy = ( 2.0*Ri )/r2 - ( 4.0*x[1]*x[1]*Ri )/(r2*r2);
      
      dvdx = ( -2.0*Ri )/r2 + ( 4.0*x[0]*x[0]*Ri )/(r2*r2);
      dvdy = ( 4.0 * x[1] * Ri * x[0] ) / ( r2*r2 );

      dPdx = ( 2.0 * ( pow( rho , gamma ) ) * gm1o2 * Mi*Mi * Ri*Ri * x[0] ) /
	     ( gm1 * r2*r2 * rhoi );

      dPdy = ( 2.0 * ( pow( rho , gamma ) ) * gm1o2 * Mi*Mi * Ri*Ri * x[1] ) /
	     ( gm1 * r2*r2 * rhoi );

    }

  
  df1 = rho*dudx + u*drdx + rho*dvdy + v*drdy;
  df2 = 2.*u*rho*dudx + u*u*drdx + dPdx + u*v*drdy + rho*u*dvdy + rho*v*dudy;
  df3 = u*v*drdx + rho*v*dudx + rho*u*dvdx + 2.*rho*v*dvdy + v*v*drdy + dPdy;
  df4 = gamma/(gamma-1.) * ( P*dudx + u*dPdx + P*dvdy + v*dPdy ) +
        0.5*( u*drdx*(u*u + v*v) + rho*dudx*(u*u + v*v) + 2.*rho*u*( u*dudx + v*dvdx) +
	      v*drdy*(u*u + v*v) + rho*dvdy*(u*u + v*v) + 2.*rho*v*( u*dudy + v*dvdy) ) ;


  q[0] = df1;
  q[1] = df2;
  q[2] = df3;
  q[3] = df4;



  return;
}



void MMS_Compute_Source_Terms ( GRID *grid )
{
  int n;                                 // Node loop counter.
  int i,j,k,e;                           // Loop counters.
  int v;                                 // Vector loop.
  int ind;                               // Component loop.
  int nodes[MAX_NUM_VERT];               // Element vertices.
  int node;                              // Node index.
  int num_vert;                          // Number of vertices for an element.
  int left_id;                           // ID of left edge.
  int right_id;                          // ID of right edge.
  double xc[NDIM];                       // Vectors.
  double xmidL[NDIM];
  double xmidR[NDIM];
  double v_left[NDIM];
  double v_right[NDIM];
  double v_middle[NDIM];
  double bary_coord[7][3];               // The barycentric coordinates.
  double wg[7];                          // The associated weights;
  double b_r = (6. - sqrt(15.))/21.;     // Precompute some the coordinates.
  double b_t = (6. + sqrt(15.))/21.;
  double b_s = (9. + 2.*sqrt(15.))/21.;
  double b_u = (9. - 2.*sqrt(15.))/21.;
  double wA = (155. - sqrt(15.))/1200.;
  double wB = (155. + sqrt(15.))/1200.;

  double XGP[NDIM];
  double QGP[NUM_VAR];
  double QTRI[NUM_VAR];

  double gp[7][3];                       // Gauss points for integration.
  double dA;                             // Differential area (for a triangle).

  // Curved triangle stuff.
  double ct_gp[7][2];                    // Gauss integration points in (r,s) space (xi,eta).
  double ct_nodes[7][2];                 // The x,y locations of the 6 nodes defining the p2 triangle (curved side).
  double N1,N2,N3,N4,N5,N6;              // Shape functions.
  double N1r,N1s,N2r,N2s,N3r,N3s;        // Shape function derivatives.
  double N4r,N4s,N5r,N5s,N6r,N6s;
  double xr,xs,yr,ys,jacobian;           // Transformation jacobian.
  double XI,ETA;                         // Triangle coordinates.
  int b,bc,seg,ct_flag=0;
  int ghost_node;

  
  // Clean out memory because we are accumulating.
  for ( n=1; n <= grid->nn; n++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	grid->Source[n*NUM_VAR+j] = 0.;
    }
  
  // Now we calculate the area-averaged value of the initial condition using Gaussian quadrature
  
  // Set up the barycentric coordinates and their respective weights.
  bary_coord[0][0] = 1./3.; bary_coord[0][1] = 1./3.; bary_coord[0][2] = 1./3.;
  bary_coord[1][0] = b_r;   bary_coord[1][1] = b_r;   bary_coord[1][2] = b_s;
  bary_coord[2][0] = b_r;   bary_coord[2][1] = b_s;   bary_coord[2][2] = b_r;
  bary_coord[3][0] = b_s;   bary_coord[3][1] = b_r;   bary_coord[3][2] = b_r;
  bary_coord[4][0] = b_t;   bary_coord[4][1] = b_t;   bary_coord[4][2] = b_u;
  bary_coord[5][0] = b_t;   bary_coord[5][1] = b_u;   bary_coord[5][2] = b_t;
  bary_coord[6][0] = b_u;   bary_coord[6][1] = b_t;   bary_coord[6][2] = b_t;
      
  wg[0] = 9./40.;
  wg[1] = wA;
  wg[2] = wA;
  wg[3] = wA;
  wg[4] = wB;
  wg[5] = wB;
  wg[6] = wB;
  
  // Lets also map these guys onto the reference triangle for integration on the triangles with curved sides.
  for ( j=0; j < 7; j++ )
    {
      ct_gp[j][0] = (bary_coord[j][0] * 0. +
		     bary_coord[j][1] * 1. +
		     bary_coord[j][2] * 0. );
      
      ct_gp[j][1] = (bary_coord[j][0] * 0. +
		     bary_coord[j][1] * 0. +
		     bary_coord[j][2] * 1. );
      
    }  // The weights can be reused since everything is consistent in the ordering.
  
  // Start the process by looping over all the elements and looping through the internal median dual pieces.
  for ( e=Tri; e <= Quad; e++ )                 // Element types loop.
    {
      // Find the number of nodes and edges for this element type.
      num_vert = NumberNodesForElement(e);

      for ( j=1; j <= grid->num_elem[e]; j++ )              // Element loop.
	{
	  // Retrieve the element node ids. This is for convenience.
	  for ( k=0; k < num_vert; k++ )              // Element node extraction loop.
	    {
	      nodes[k] = grid->c2n[e][j*num_vert+k];
	    }                                         // End element node extraction loop.
	  
	  // Find the element centroid.
	  xc[0] = 0.; xc[1] = 0.;
	  for ( k=0; k < num_vert; k++ )
	    {
	      for ( i=0; i < NDIM; i++ )
		{
		  xc[i] += grid->x[NDIM*nodes[k]+i];
		}
	    }
	  xc[0] /= ((double)num_vert);
	  xc[1] /= ((double)num_vert);

	  // Loop through each node in the element.
	  for ( k=0; k < num_vert; k++ )              // Element node loop.
	    {
	      node = nodes[k];

	      // Get the two neighboring vertices.
	      left_id = (k+1)%num_vert;
	      right_id = (k+num_vert-1)%num_vert;
		  
	      left_id = grid->c2n[e][j*num_vert+left_id];
	      right_id = grid->c2n[e][j*num_vert+right_id];
	      
	      // Get the edge midpoints.
	      for ( i=0; i < NDIM; i++ )
		{
		  xmidL[i] = grid->x[NDIM*node+i] + grid->x[NDIM*left_id+i];
		  xmidR[i] = grid->x[NDIM*node+i] + grid->x[NDIM*right_id+i];
		}

	      for ( i=0; i < NDIM; i++ )
		{
		  xmidL[i] *= 0.5;
		  xmidR[i] *= 0.5;
		}

	      // Make the vectors.
	      for ( i=0; i < NDIM; i++ )
		{
		  v_left[i]   = xmidL[i] - grid->x[NDIM*node+i];
		  v_right[i]  = xmidR[i] - grid->x[NDIM*node+i];
		  v_middle[i] = xc[i]    - grid->x[NDIM*node+i];
		}

	      ct_flag = 0;
	      
	      // If the node and its neighbor are both boundary nodes, we need to gather some information.
	      if ( grid->node_state[node] == BOUNDARY && grid->node_state[left_id] == BOUNDARY )
		{
		  // Retrieve the relevant information.
		  
		  // First find the correct bedge. It will have both nodes.
		  for ( i=1; i <= grid->nbedges; i++ )
		    {
		      if ( node == grid->bedges[i*5] )  // candidate bedge.
			{
			  b = grid->bedges[i*5+3];
			  seg = grid->bedges[i*5+4];
			  
			  if ( (  node == grid->bs[b][seg][0]   ||  node == grid->bs[b][seg][1]   ) &&
			       ( left_id == grid->bs[b][seg][0] || left_id == grid->bs[b][seg][1] )   )  // we have a winner.
			    {
			      ghost_node = grid->bedges[i*5+1] + grid->nn;
			      b = grid->bedges[i*5+3];
			      bc = grid->bbc[b];
			      seg = grid->bedges[i*5+4];
			      
			      if ( bc > 10 )
				ct_flag = 1;
			      
			      break;
			    }
			}
		    }
		}

	      
	      if ( ct_flag )
		{
		  // Now we need to set up the six nodes used by the shape functions.
		      
		  // Node 1 is the element centroid.
		  ct_nodes[1][0] = xc[0];
		  ct_nodes[1][1] = xc[1];

		  // Node 2 is the node we are working on.
		  ct_nodes[2][0] = grid->x[NDIM*node+0];
		  ct_nodes[2][1] = grid->x[NDIM*node+1];

		  // Third node is the ghost node.
		  //ct_nodes[3][0] = grid->x[NDIM*ghost_node+0];
		  //ct_nodes[3][1] = grid->x[NDIM*ghost_node+1];  // this might be incorrect when the edges were created from linear boundaries.

		  curved_boundary_midpoint(grid, bc, &(ct_nodes[2][0]), &(grid->x[left_id*NDIM+0]), &(ct_nodes[3][0]) );

		  // Node 4 is the midpoint of 1,2.
		  ct_nodes[4][0] = 0.5*(ct_nodes[1][0] + ct_nodes[2][0]);
		  ct_nodes[4][1] = 0.5*(ct_nodes[1][1] + ct_nodes[2][1]);

		  // Node 5 is the midpoint of the curved boundary between node and ghost node.
		  curved_boundary_midpoint(grid, bc, &(ct_nodes[2][0]), &(ct_nodes[3][0]), &(ct_nodes[5][0]) );

		  // Node 6 is the midpoint of 1,3.
		  ct_nodes[6][0] = 0.5*(ct_nodes[1][0] + ct_nodes[3][0]);
		  ct_nodes[6][1] = 0.5*(ct_nodes[1][1] + ct_nodes[3][1]);

		  // First we need to make sure that the concavity constraint is met. Under a linear transformation
		  // to the reference element, node 5 must be greater than 0.25 in both r and s coordinates.
		  // See Flaherty's notes online or Claes Johnson's book.

		  XI = ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
			 (ct_nodes[5][0]*ct_nodes[3][1] - ct_nodes[5][1]*ct_nodes[3][0]) +
			 (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
		       ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
		         (ct_nodes[2][0]*ct_nodes[3][1] - ct_nodes[2][1]*ct_nodes[3][0]) +
		         (ct_nodes[2][0]*ct_nodes[1][1] - ct_nodes[2][1]*ct_nodes[1][0])  );
		  
		  ETA = ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			  (ct_nodes[5][0]*ct_nodes[2][1] - ct_nodes[5][1]*ct_nodes[2][0]) +
			  (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
		        ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			  (ct_nodes[3][0]*ct_nodes[2][1] - ct_nodes[3][1]*ct_nodes[2][0]) +
			  (ct_nodes[3][0]*ct_nodes[1][1] - ct_nodes[1][0]*ct_nodes[3][1])  );

		  if ( ETA < 0.25 || XI < 0.25 )  // The transformation will not be invertible.
		    {
		      printf("CRITICAL ERROR: Triangle 1 of element %d of type %d is too concave (compute MMS source term: solve.C)!\n",j,e);
		      printf("  xi = %.15e     eta = %.15e\n",XI,ETA);
		      fflush(stdout);
		      exit(1);
		    }

		  // Do this for each variable.

		  QTRI[0] = 0.; QTRI[1] = 0.; QTRI[2] = 0.; QTRI[3] = 0.;
		  
		  // Now I can loop over the Gaussian integration nodes and apply the quadratic transformation.
		  for ( i=0; i < 7; i++ )
		    {
		      XI = ct_gp[i][0];
		      ETA = ct_gp[i][1];
		      
		      // Compute basis functions.
		      N1 = 1. - 3.*XI - 3.*ETA + 4.*XI*ETA + 2.*XI*XI + 2.*ETA*ETA;
		      N2 = XI*(2.*XI - 1.);
		      N3 = ETA*(2.*ETA - 1.);
		      N4 = 4.*XI*(1. - XI - ETA);
		      N5 = 4.*XI*ETA;
		      N6 = 4.*ETA*(1. - XI - ETA);

		      // Compute their derivatives.
		      N1r = -3. + 4.*ETA + 4.*XI;
		      N1s = -3. + 4.*XI + 4.*ETA;
		      N2r = 4.*XI - 1.;
		      N2s = 0.;
		      N3r = 0.;
		      N3s = 4.*ETA - 1.;
		      N4r = 4.*(1. - 2.*XI - ETA);
		      N4s = -4.*XI;
		      N5r = 4.*ETA;
		      N5s = 4.*XI;
		      N6r = -4.*ETA;
		      N6s = 4.*(1. - 2.*ETA - XI);

		      // Compute the derivatives of x,y wrt to the transformation.
		      xr = N1r*ct_nodes[1][0] + N2r*ct_nodes[2][0] + N3r*ct_nodes[3][0] + N4r*ct_nodes[4][0]
			 + N5r*ct_nodes[5][0] + N6r*ct_nodes[6][0];
		      xs = N1s*ct_nodes[1][0] + N2s*ct_nodes[2][0] + N3s*ct_nodes[3][0] + N4s*ct_nodes[4][0]
			 + N5s*ct_nodes[5][0] + N6s*ct_nodes[6][0];
		      yr = N1r*ct_nodes[1][1] + N2r*ct_nodes[2][1] + N3r*ct_nodes[3][1] + N4r*ct_nodes[4][1]
			 + N5r*ct_nodes[5][1] + N6r*ct_nodes[6][1];
		      ys = N1s*ct_nodes[1][1] + N2s*ct_nodes[2][1] + N3s*ct_nodes[3][1] + N4s*ct_nodes[4][1]
			 + N5s*ct_nodes[5][1] + N6s*ct_nodes[6][1];

		      jacobian = xr*ys - xs*yr;

		      // Now I need to extrapolate to this point. First we get the real coordinates of the
		      // Gaussian point.
		      XGP[0] = N1*ct_nodes[1][0] + N2*ct_nodes[2][0] + N3*ct_nodes[3][0] +
			       N4*ct_nodes[4][0] + N5*ct_nodes[5][0] + N6*ct_nodes[6][0];

		      XGP[1] = N1*ct_nodes[1][1] + N2*ct_nodes[2][1] + N3*ct_nodes[3][1] +
			       N4*ct_nodes[4][1] + N5*ct_nodes[5][1] + N6*ct_nodes[6][1];

		      MMS_Get_Flux_Divergence ( XGP, QGP );

		      // Accumulate.
		      for ( v=0; v < NUM_VAR; v++ )
			{
			  QTRI[v] += ( wg[i] * QGP[v] *jacobian * 0.5);
			} 
		      
		    }

		  // Add to the running total for the node.
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      grid->Source[node*NUM_VAR+v] += QTRI[v];
		    }
		  
		}
	      else
		{
		  // Get the area.
		  dA = 0.5*( v_left[0]*v_middle[1] - v_left[1]*v_middle[0] );

		  // Initialize the Gauss points.
		  for ( v=0; v < 7; v++ )
		    {
		      for ( ind=0; ind < NDIM; ind++ )
			{
			  gp[v][ind] = (bary_coord[v][0] * grid->x[node*NDIM+ind] +
					bary_coord[v][1] * xmidL[ind] +
					bary_coord[v][2] * xc[ind]);
			}
		      gp[v][2] = 0.;   // Always 0 since 2D vectors.
		    }

		  // Now we loop over the Gauss points and get the function value
		  // and perform the numerical integration (the area integral of the triangle).
	      
		  // Do this for each variable.
		  
		  QTRI[0] = 0.; QTRI[1] = 0.; QTRI[2] = 0.; QTRI[3] = 0.;
		  
		  for ( i=0; i < 7; i++ )
		    {
		      // Get the value at the point.
		      XGP[0] = gp[i][0];
		      XGP[1] = gp[i][1];
		      
		      MMS_Get_Flux_Divergence ( XGP , QGP );
		      
		      for ( v=0; v < NUM_VAR; v++ )
			{
			  QTRI[v] += ( wg[i] * QGP[v] );
			} 
		    }
	      
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      QTRI[v] *= dA;
		    }
		  
		  // Accumulate the integral to the node.
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      grid->Source[node*NUM_VAR+v] += QTRI[v];
		    }
		}
		 
	      
	      // Process the second triangle.

	      ct_flag = 0;
	      
	      // If the node and its neighbor are both boundary nodes, we need to gather some information.
	      if ( grid->node_state[node] == BOUNDARY && grid->node_state[right_id] == BOUNDARY )
		{
		  // Retrieve the relevant information.
		  
		  // First find the correct bedge. It will have both nodes.
		  for ( i=1; i <= grid->nbedges; i++ )
		    {
		      if ( node == grid->bedges[i*5] )  // candidate bedge.
			{
			  b = grid->bedges[i*5+3];
			  seg = grid->bedges[i*5+4];
			  
			  if ( (   node == grid->bs[b][seg][0]   ||   node == grid->bs[b][seg][1]   ) &&
			       ( right_id == grid->bs[b][seg][0] || right_id == grid->bs[b][seg][1] )   )  // we have a winner.
			    {
			      ghost_node = grid->bedges[i*5+1] + grid->nn;
			      b = grid->bedges[i*5+3];
			      bc = grid->bbc[b];
			      seg = grid->bedges[i*5+4];
			      if ( bc > 10 )
				ct_flag = 1;
			      
			      break;
			    }
			}
		    }
		}

	      if ( ct_flag )
		{
		  // Now we need to set up the six nodes used by the shape functions.
		      
		  // Node 1 is the element centroid.
		  ct_nodes[1][0] = xc[0];
		  ct_nodes[1][1] = xc[1];

		  // Node 2 is the ghost_node.
		  //ct_nodes[2][0] = grid->x[NDIM*ghost_node+0];
		  //ct_nodes[2][1] = grid->x[NDIM*ghost_node+1];  // this might be incorrect when the edges were created from linear boundaries.

		  curved_boundary_midpoint(grid, bc, &(grid->x[right_id*NDIM+0]), &(grid->x[node*NDIM+0]), &(ct_nodes[2][0]) );

		  // Third node is the node we are working on.
		  ct_nodes[3][0] = grid->x[NDIM*node+0];
		  ct_nodes[3][1] = grid->x[NDIM*node+1];

		  // Node 4 is the midpoint of 1,2.
		  ct_nodes[4][0] = 0.5*(ct_nodes[1][0] + ct_nodes[2][0]);
		  ct_nodes[4][1] = 0.5*(ct_nodes[1][1] + ct_nodes[2][1]);

		  // Node 5 is the midpoint of the curved boundary between node and ghost node.
		  curved_boundary_midpoint(grid, bc, &(ct_nodes[2][0]), &(ct_nodes[3][0]), &(ct_nodes[5][0]) );

		  // Node 6 is the midpoint of 1,3.
		  ct_nodes[6][0] = 0.5*(ct_nodes[1][0] + ct_nodes[3][0]);
		  ct_nodes[6][1] = 0.5*(ct_nodes[1][1] + ct_nodes[3][1]);

		  // First we need to make sure that the concavity constraint is met. Under a linear transformation
		  // to the reference element, node 5 must be greater than 0.25 in both r and s coordinates.
		  // See Flaherty's notes online or Claes Johnson's book.

		  XI = ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
			 (ct_nodes[5][0]*ct_nodes[3][1] - ct_nodes[5][1]*ct_nodes[3][0]) +
			 (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
		       ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
		         (ct_nodes[2][0]*ct_nodes[3][1] - ct_nodes[2][1]*ct_nodes[3][0]) +
		         (ct_nodes[2][0]*ct_nodes[1][1] - ct_nodes[2][1]*ct_nodes[1][0])  );
		  
		  ETA = ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			  (ct_nodes[5][0]*ct_nodes[2][1] - ct_nodes[5][1]*ct_nodes[2][0]) +
			  (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
		        ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			  (ct_nodes[3][0]*ct_nodes[2][1] - ct_nodes[3][1]*ct_nodes[2][0]) +
			  (ct_nodes[3][0]*ct_nodes[1][1] - ct_nodes[1][0]*ct_nodes[3][1])  );

		  if ( ETA < 0.25 || XI < 0.25 )  // The transformation will not be invertible.
		    {
		      printf("CRITICAL ERROR: Triangle 1 of element %d of type %d is too concave (compute MMS Source term: solve.C)!\n",j,e);
		      printf("  xi = %.15e     eta = %.15e\n",XI,ETA);
		      fflush(stdout);
		      exit(1);
		    }

		  // Do this for each variable.

		  QTRI[0] = 0.; QTRI[1] = 0.; QTRI[2] = 0.; QTRI[3] = 0.;
		  
		  // Now I can loop over the Gaussian integration nodes and apply the quadratic transformation.
		  for ( i=0; i < 7; i++ )
		    {
		      XI = ct_gp[i][0];
		      ETA = ct_gp[i][1];
		      
		      // Compute basis functions.
		      N1 = 1. - 3.*XI - 3.*ETA + 4.*XI*ETA + 2.*XI*XI + 2.*ETA*ETA;
		      N2 = XI*(2.*XI - 1.);
		      N3 = ETA*(2.*ETA - 1.);
		      N4 = 4.*XI*(1. - XI - ETA);
		      N5 = 4.*XI*ETA;
		      N6 = 4.*ETA*(1. - XI - ETA);

		      // Compute their derivatives.
		      N1r = -3. + 4.*ETA + 4.*XI;
		      N1s = -3. + 4.*XI + 4.*ETA;
		      N2r = 4.*XI - 1.;
		      N2s = 0.;
		      N3r = 0.;
		      N3s = 4.*ETA - 1.;
		      N4r = 4.*(1. - 2.*XI - ETA);
		      N4s = -4.*XI;
		      N5r = 4.*ETA;
		      N5s = 4.*XI;
		      N6r = -4.*ETA;
		      N6s = 4.*(1. - 2.*ETA - XI);

		      // Compute the derivatives of x,y wrt to the transformation.
		      xr = N1r*ct_nodes[1][0] + N2r*ct_nodes[2][0] + N3r*ct_nodes[3][0] + N4r*ct_nodes[4][0]
			 + N5r*ct_nodes[5][0] + N6r*ct_nodes[6][0];
		      xs = N1s*ct_nodes[1][0] + N2s*ct_nodes[2][0] + N3s*ct_nodes[3][0] + N4s*ct_nodes[4][0]
			 + N5s*ct_nodes[5][0] + N6s*ct_nodes[6][0];
		      yr = N1r*ct_nodes[1][1] + N2r*ct_nodes[2][1] + N3r*ct_nodes[3][1] + N4r*ct_nodes[4][1]
			 + N5r*ct_nodes[5][1] + N6r*ct_nodes[6][1];
		      ys = N1s*ct_nodes[1][1] + N2s*ct_nodes[2][1] + N3s*ct_nodes[3][1] + N4s*ct_nodes[4][1]
			 + N5s*ct_nodes[5][1] + N6s*ct_nodes[6][1];

		      jacobian = xr*ys - xs*yr;

		      // Now I need to extrapolate to this point. First we get the real coordinates of the
		      // Gaussian point.
		      XGP[0] = N1*ct_nodes[1][0] + N2*ct_nodes[2][0] + N3*ct_nodes[3][0] +
			       N4*ct_nodes[4][0] + N5*ct_nodes[5][0] + N6*ct_nodes[6][0];

		      XGP[1] = N1*ct_nodes[1][1] + N2*ct_nodes[2][1] + N3*ct_nodes[3][1] +
			       N4*ct_nodes[4][1] + N5*ct_nodes[5][1] + N6*ct_nodes[6][1];

		      MMS_Get_Flux_Divergence ( XGP, QGP );
		      
		      // Accumulate.
		      for ( v=0; v < NUM_VAR; v++ )
			{
			  QTRI[v] += ( wg[i] * QGP[v] *jacobian * 0.5);
			} 
		      
		    }
		  
		  // Add to the running total for the node.
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      grid->Source[node*NUM_VAR+v] += QTRI[v];
		    }
		  
		}
	      else
		{
		  // Get the area.
		  dA = 0.5*( v_middle[0]*v_right[1] - v_middle[1]*v_right[0] );
		  
		  // Initialize the Gauss points.
		  for ( v=0; v < 7; v++ )
		    {
		      for ( ind=0; ind < NDIM; ind++ )
			{
			  gp[v][ind] = (bary_coord[v][0] * grid->x[node*NDIM+ind] +
					bary_coord[v][1] * xmidR[ind] +
					bary_coord[v][2] * xc[ind]);
			}
		      gp[v][2] = 0.;   // Always 0 since 2D vectors.
		    }
	      
		  // Now we loop over the Gauss points and get the function value
		  // and perform the numerical integration (the area integral of the triangle).
	      
		  // Do this for each variable.
		  
		  QTRI[0] = 0.; QTRI[1] = 0.; QTRI[2] = 0.; QTRI[3] = 0.;
		  
		  for ( i=0; i < 7; i++ )
		    {
		      // Get the value at the point.
		      XGP[0] = gp[i][0];
		      XGP[1] = gp[i][1];
		      
		      MMS_Get_Flux_Divergence ( XGP , QGP );
		  
		      for ( v=0; v < NUM_VAR; v++ )
			{
			  QTRI[v] += ( wg[i] * QGP[v] );
			} 
		    }
		  
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      QTRI[v] *= dA;
		    }
		  
		  // Accumulate the integral to the node.
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      grid->Source[node*NUM_VAR+v] += QTRI[v];
		    }
		}
	    }                                       // End element node loop.
	  
	}                                           // End element loop.
      
    }                                               // End element type loop.
  
  
  return;
}


void MMS_Compute_Source_Terms_Exp ( GRID *grid )  // Special function for computing the source terms of the exponential functions as contour integrals on the cv boundaries.
{
  int i,j,k,n;                           // Loop counters.
  int B,S;                               // Boundary loop counters.
  int gelem;                             // Global element.
  int nodeL, nodeR;                      // Edge/element nodes.
  int num_vert;                          // Number of vertices for an element.
  int n0,n1;                             // Boundary segment nodes.
  double xc, yc;                         // Values at the centroid.
  double xmid, ymid;                     // Values at face midpoints.
  double nx,ny;                          // Normal vector components.
  double len;                            // Magnitude of a vector.
  double a,b,c,d;                        // Variables for the integral parameterization.
  double F1,F2,F3,F4;                    // Components of the contour integral.

  
  // Clean out memory because we are accumulating.
  for ( n=1; n <= grid->nn; n++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	grid->Source[n*NUM_VAR+j] = 0.;
    }

  // Loop over the element types and elements to get the median dual pieces.
  for ( i=Tri; i <= Quad; i++ )                                                       // Element types loop.
    {
      // Find the number of nodes and edges for this element type.
      num_vert = NumberNodesForElement(i);

      for ( j=1; j <= grid->num_elem[i]; j++ )                                        // Element loop.
	{
	  // Get the global element id.
	  gelem = LocalToGlobalID ( grid->num_elem,
				    j,
				    i );

	  // Get the centroid.
	  xc = grid->el_cent[gelem*NDIM];
          yc = grid->el_cent[gelem*NDIM+1];

	  // Loop over the element nodes (edges).
	  for ( k=0; k < num_vert; k++ )                                            // Element node loop.
	    {
	      nodeL = k;
	      nodeR = (k+1)%num_vert;
	      
	      // Get the global node numbers.
	      nodeL = grid->c2n[i][j*num_vert + nodeL];
	      nodeR = grid->c2n[i][j*num_vert + nodeR];
	      
	      // Calculate the edge midpoint.
	      xmid = 0.5*( grid->x[nodeL*NDIM+0] + grid->x[nodeR*NDIM+0] );
	      ymid = 0.5*( grid->x[nodeL*NDIM+1] + grid->x[nodeR*NDIM+1] );
	      
	      nx = yc - ymid;
	      ny = xmid - xc;

	      // Calculate the length
	      len = sqrt( nx*nx + ny*ny );

	      nx = nx / len;
	      ny = ny / len;

	      a = xmid;
	      b = xc - xmid;
	      c = ymid;
	      d = yc - ymid;

	      F1 = (-0.5) * (nx+ny) * len * (1.0 / (d+b) ) * ( exp( 2.0*a + 2.0*c - 4.0 ) -
							       exp( 2.0*a + 2.0*b + 2.0*c + 2.0*d - 4.0) );

	      F2 = (1.0/6.0) * len * (1.0 / (b+d) ) * ( exp( 3.0*a + 3.0*c - 6.0 ) * ( nx + 1.4*nx + 2.0*ny ) * ( -1.0 + exp( 3.0*b + 3.0*d ) ) );

	      F3 = (1.0/6.0) * len * (1.0 / (b+d) ) * ( exp( 3.0*a + 3.0*c - 6.0 ) * ( 2.0*nx + 1.4*ny + ny ) * ( -1.0 + exp( 3.0*b + 3.0*d ) ) );

	      F4 = (1.0/8.0) * len * (1.0 / (b+d) ) * ( exp( 4.0*a + 4.0*c - 8.0 ) * ( 2.0*nx + 1.4*nx + 2.0*ny + 1.4*ny ) * ( -1.0 + exp( 4.0*b + 4.0*d ) ) );

	      // Now we add to the left node and subtract from the right node (since the normal is flipped).
	      grid->Source[nodeL*NUM_VAR+0] += F1;
	      grid->Source[nodeL*NUM_VAR+1] += F2;
	      grid->Source[nodeL*NUM_VAR+2] += F3;
	      grid->Source[nodeL*NUM_VAR+3] += F4;
	      
	      grid->Source[nodeR*NUM_VAR+0] -= F1;
	      grid->Source[nodeR*NUM_VAR+1] -= F2;
	      grid->Source[nodeR*NUM_VAR+2] -= F3;
	      grid->Source[nodeR*NUM_VAR+3] -= F4;
	    }                                                                       // End element node loop. 
	}                                                                           // End element loop.
    }                                                                               // End element type loop.

  // Loop over the boundary edges to close of the cv's.
  for ( B=1; B <= grid->nb; B++ )                                         // Boundary loop.
    {                                             
      for ( S=1; S <= grid->nbs[B]; S++ )                                 // Segment loop.
	{                                  
	  n0 = grid->bs[B][S][0];
	  n1 = grid->bs[B][S][1];
	  
	  // the edge from n0 to n1 is consistent with RHR, so the normals will be outward pointing for the two pieces.
	  
	  xmid = (grid->x[NDIM*n0]+grid->x[NDIM*n1])/2.;
	  ymid = (grid->x[NDIM*n0+1]+grid->x[NDIM*n1+1])/2.;

	  // Get the normal and add it to the control volumes.
	  nx = grid->x[NDIM*n1+1] - ymid;
	  ny = xmid - grid->x[NDIM*n1];

	  // Its the same length and direction for both halves of the boundary segment.

	  // Calculate the length
	  len = sqrt( nx*nx + ny*ny );
	  
	  nx = nx / len;
	  ny = ny / len;

	  // First piece. We integrate from the left node (n0) to the midpoint.

	  a = grid->x[n0*NDIM+0];
	  b = xmid - a;
	  c = grid->x[n0*NDIM+1];
	  d = ymid - c;
	  
	  F1 = (-0.5) * (nx+ny) * len * (1.0 / (d+b) ) * ( exp( 2.0*a + 2.0*c - 4.0 ) -
							   exp( 2.0*a + 2.0*b + 2.0*c + 2.0*d - 4.0) );

	  F2 = (1.0/6.0) * len * (1.0 / (b+d) ) * ( exp( 3.0*a + 3.0*c - 6.0 ) * ( nx + 1.4*nx + 2.0*ny ) * ( -1.0 + exp( 3.0*b + 3.0*d ) ) );
	  
	  F3 = (1.0/6.0) * len * (1.0 / (b+d) ) * ( exp( 3.0*a + 3.0*c - 6.0 ) * ( 2.0*nx + 1.4*ny + ny ) * ( -1.0 + exp( 3.0*b + 3.0*d ) ) );
	  
	  F4 = (1.0/8.0) * len * (1.0 / (b+d) ) * ( exp( 4.0*a + 4.0*c - 8.0 ) * ( 2.0*nx + 1.4*nx + 2.0*ny + 1.4*ny ) * ( -1.0 + exp( 4.0*b + 4.0*d ) ) );

	  // Now we add to the node.
	  grid->Source[n0*NUM_VAR+0] += F1;
	  grid->Source[n0*NUM_VAR+1] += F2;
	  grid->Source[n0*NUM_VAR+2] += F3;
	  grid->Source[n0*NUM_VAR+3] += F4;


	  // Second piece. We integrate from the midpoint to the right node (n1) so as to be consistent with the RHR.

	  a = xmid;
	  b = grid->x[n1*NDIM+0] - a;
	  c = ymid;
	  d = grid->x[n1*NDIM+1] - c;

	  F1 = (-0.5) * (nx+ny) * len * (1.0 / (d+b) ) * ( exp( 2.0*a + 2.0*c - 4.0 ) -
							   exp( 2.0*a + 2.0*b + 2.0*c + 2.0*d - 4.0) );

	  F2 = (1.0/6.0) * len * (1.0 / (b+d) ) * ( exp( 3.0*a + 3.0*c - 6.0 ) * ( nx + 1.4*nx + 2.0*ny ) * ( -1.0 + exp( 3.0*b + 3.0*d ) ) );
	  
	  F3 = (1.0/6.0) * len * (1.0 / (b+d) ) * ( exp( 3.0*a + 3.0*c - 6.0 ) * ( 2.0*nx + 1.4*ny + ny ) * ( -1.0 + exp( 3.0*b + 3.0*d ) ) );
	  
	  F4 = (1.0/8.0) * len * (1.0 / (b+d) ) * ( exp( 4.0*a + 4.0*c - 8.0 ) * ( 2.0*nx + 1.4*nx + 2.0*ny + 1.4*ny ) * ( -1.0 + exp( 4.0*b + 4.0*d ) ) );

	  // Now we add to the node.
	  grid->Source[n1*NUM_VAR+0] += F1;
	  grid->Source[n1*NUM_VAR+1] += F2;
	  grid->Source[n1*NUM_VAR+2] += F3;
	  grid->Source[n1*NUM_VAR+3] += F4;
	  
	}                                                                   // End segment loop.
    }                                                                       // End boundary loop.

  // Done we have accumulated the source terms to the appropriate nodes. The source terms here represent the exact flux as computed by contour integrals. Note, here, that since we
  // have full support of the functions across the domain, QL == QR so the analytical flux calculation preserves the correct answer. Well, of course it does since any numerical flux
  // must correspond to the anayltical Euler flux when the left and right states are identical.

  return;
}

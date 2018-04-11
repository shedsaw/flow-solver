//=============================================================
// 
//  main.C
//  
//  This function acts as a driver program to run the Euler code
//  on an unstructured 2D mesh.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "grid.h"
#include "maps.h"
#include "params.h"
#include "Defined_Variables.h"
#include "sys_mem.h"
#include "grid_io.h"
#include "cv_calc.h"
#include "geometry.h"
#include "linear_algebra.h"
#include "mesh_io.h"
#include "mesh_utils.h"
#include "metrics.h"
#include "reorder.h"
#include "reconstruct.h"
#include "ic.h"
#include "gradient.h"
#include "hessian.h"
#include "liblapack_interface.h"
#include "init_memory.h"
#include "residual.h"
#include "flux.h"
#include "solve.h"
#include "curved_boundaries.h"
#include "compute_derivatives.h"
#include "limiter.h"

void compute_integral_error ( GRID *grid, PARAMS p );


int main(int argcs, char* pArgs[])
{

  //===================================================
  // Variable Declaration block.
  //===================================================

  int i, j, k;                         // Loop counters.
  int n, gn;                           // Loop counters - n for node,gn - ghost node.
  int sidx, eidx;
  int num_neighbors;

  int post_process_now;
  int terminate_now;
  int halt_foits = 0;
  int halt_soits = 0;
  int halt_toits = 0;

  int num_neighbors_needed;

  FILE *fp;                            // File pointer.
  FILE *tfp;                           // File pointer for time info.
  FILE *cpfp;                          // File pointer for core pressure (Vortex Convection only).
  FILE *command;                       // File pointer to the control file.
  FILE *bforce;                        // File pointer to the body forces.

  struct timeval start_time;           // Structure to hold the time program started running.
  struct timeval end_time;             // Structure to hold the time at the end of the current iteration.

  double dt;                           // Delta time between current iteration and when program began execution.

  char buffer[100];                    // Static string buffer.
  char buff[100];
  int buffsize = 100;
  
  GRID grid;                           // Pointer to the grid object.
  PARAMS p;                            // Parameter object.
  SYS_MEM smem;                        // System memory object.

  double *rms;                         // RMS error of the Residuals.

  double RMS;                          // RMS error returned from the linear_solve function.

  double time, cpressure;              // Core pressure at time t.

  double Bforces[2];

  double *QexC = NULL;
  double *QexP = NULL;

  //===================================================


    
  if ( argcs != 2 )
    {
      printf("Usage: ./field_solver <in.mesh>\n");
      exit(1);
    }

  read_parameters(&p);

  /*
  printf("  After reading parameters, now in main().\n");
  printf("  Number of iterations = %d\n",p.nsteps);
  printf("  Number of Newton iterations = %d\n",p.iter);
  printf("  Jacobian Update frequency = %d\n",p.jufreq);
  printf("  CFL Number = %lf\n",p.CFL);
  printf("  Number of iterations to ramp = %d\n",p.ramp);
  printf("  Gamma = %lf\n",p.gamma);
  printf("  Alpha = %lf\n",p.alpha);
  printf("  Mach Number = %lf\n",p.mach);
  printf("  Number of subiterations = %d\n",p.subits);
  printf("  Unsteady flag = %d\n",p.unsteady);
  printf("  Order of time accuracy = %d\n",p.tacc);
  printf("  Minimum time step = %lf\n",p.mintime);
  printf("  Method = %d\n",p.method);
  printf("  Order = %d\n",p.order);
  printf("  Number of first order iterations = %d\n",p.foits);
  printf("  Number of second order iterations = %d\n",p.soits);
  printf("  Number of third order iterations = %d\n,p.toits);
  printf("  Isrestart = %d\n",p.isrestart);
  printf("  Limit = %d\n",p.limit);
  printf("  Reorder = %d\n",p.reorder);
  printf("  IC = %d\n",p.ic);
  printf("  Point or CV values = %d\n",p.p_or_cv);
  printf("  Start unsteady file write at iteration = %d\n",p.start_post_pro);
  printf("  Frequency of post processing = %d\n",p.post_pro_freq);
  */
  
  // Read in the grid.
  read_mesh( pArgs[1],&grid,p);

  strcpy(grid.gridname,pArgs[1]);

  // Throw a fit if I something wrong.
  if ( CYLINDER_HACK && ( p.order<3 ) )
    {
      printf("FATAL ERROR: Mixed order reconstruction is expressly forbidden when the input order is not at least 3rd.\n");
      fflush(stdout);
      exit(0);
    }

  // Reorder the nodes.
  if ( p.isrestart == 0 && p.reorder == 1 )
    Cuthill_McKee(&grid);
      
  // Check the mesh metrics.
  metrics2D(&grid);

  // Generate the cells surrounding nodes map.
  printf("Building list of cells surrounding vertex........");
  Cells_Surrounding_Node(grid.nn,
			 grid.num_elem,
			 grid.c2n,
			 &(grid.csn),
			 &(grid.ncsn));
  printf("Done\n");

  // Generate the nodes surrounding nodes map.
  printf("Building list of nodes surrounding vertex........\n");
  Nodes_Surrounding_Node(grid.nn,
			 grid.num_elem,
			 grid.c2n,
			 grid.csn,
			 grid.ncsn,
			 &(grid.nsn),
			 &(grid.nnsn));
  printf("Done\n");


  // Generate the edges surrounding nodes map.
  printf("Building list of edges surrounding vertex........");
  Edges_Surrounding_Node(grid.nn,
			 grid.nsn,
			 grid.nnsn,
			 &(grid.nedges),
			 &(grid.edges),
			 &(grid.esn));
  printf("Done\n");
   
  // Build the cell to edge map.
  printf("Building cell to edges map........");
  Build_Cell_to_Edge_Map(grid.c2n,
			 grid.num_elem,
			 grid.nsn,
			 grid.nnsn,
			 grid.esn,
			 &(grid.c2e));
  printf("Done\n");

  // Build the element surrounding element maps.
  printf("Building elements surrounding elements map........");
  Elements_Surrounding_Element(grid.num_elem,
			       grid.c2n,
			       grid.csn,
			       grid.ncsn,
			       &(grid.ese),
			       &(grid.nese));
  printf("Done\n");


  // Generate the nodes surrounding nodes map.
  printf("Building list of nodes surrounding vertex degree 2........\n");
  Nodes_Surrounding_Node2(grid.nn,
			  grid.num_elem,
			  grid.c2n,
			  grid.esn,
			  grid.nnsn,
			  grid.edges,
			  &(grid.nsn2),
			  &(grid.nnsn2));
  printf("Done\n");

  printf("Initializing gnnsn, gnnsn1, and gnsn........\n");
  Initialize_Generic_Nodes_Surrounding_Node ( grid.nn,
					      grid.num_elem,
					      grid.c2n,
					      grid.esn,
					      grid.nnsn,
					      grid.edges,
					      &(grid.gnsn),
					      &(grid.gnnsn),
					      &(grid.gnnsn1),
					      &(grid.node_fill_level) );
  
  printf("Done\n");

  //if ( p.order == 4 || DEBUG_DERIVS )
  //if ( p.order > 2 || DEBUG_DERIVS )
  if ( 0 )
    {
      printf("Building the generic list of nodes surrounding node........\n");
      printf("Adding first layer of neighbors........\n");
      for ( i=1; i <= grid.nn; i++ )
	{
	  Add_Layer_to_GNSN ( grid.nn, grid.num_elem, grid.c2n, grid.esn, grid.nnsn, grid.edges,
			      &(grid.gnsn), &(grid.gnnsn), grid.gnnsn1, grid.node_fill_level, i );
	}
      printf("Done\n");

      printf("Adding second layer of neighbors........\n");
      for ( i=1; i <= grid.nn; i++ )
	{
	  Add_Layer_to_GNSN ( grid.nn, grid.num_elem, grid.c2n, grid.esn, grid.nnsn, grid.edges,
			      &(grid.gnsn), &(grid.gnnsn), grid.gnnsn1, grid.node_fill_level, i );
	}
      printf("Done\n");

      if ( p.order == 4 )
	{
	  printf("Adding third layer of neighbors........\n");
	  for ( i=1; i <= grid.nn; i++ )
	    {
	      Add_Layer_to_GNSN ( grid.nn, grid.num_elem, grid.c2n, grid.esn, grid.nnsn, grid.edges,
				  &(grid.gnsn), &(grid.gnnsn), grid.gnnsn1, grid.node_fill_level, i );
	    }
	  printf("Done\n");
	}

      if ( p.order == 1 || p.order == 2 )
	{
	  num_neighbors_needed = 3;
	}
      else if ( p.order == 3 )
	{
	  num_neighbors_needed = 6;
	}
      else if ( p.order == 4 )
	{
	  num_neighbors_needed = 10;
	}

      // Now check to make sure we have enough neighbors.
      printf("Checking to make sure the nodes all have a large enough stencil.......\n");
      
      for ( i=1; i <= grid.nn; i++ )
	{
	  // Find the number of neighbors.
	  num_neighbors = grid.gnnsn[ grid.gnnsn1[i+1] ] - grid.gnnsn[ ( grid.gnnsn1[i] + 1 ) ];
	  
	  //if ( num_neighbors <= NUM_MOM )  // defecient stencil.
	  if ( num_neighbors < num_neighbors_needed )
	    {
	      for ( j=0; j < 3; j++ )
		{
		  //printf("  Node %d has %d neighbors in its stencil. Adding degree %d neighbors to its stencil.......\n",i,num_neighbors,(3+j+1));
		  printf("  Node %d has %d neighbors in its stencil. Adding degree %d neighbors to its stencil.......\n",i,num_neighbors,
			 grid.node_fill_level[i]+1);
		  Add_Layer_to_GNSN ( grid.nn, grid.num_elem, grid.c2n, grid.esn, grid.nnsn, grid.edges,
				      &(grid.gnsn), &(grid.gnnsn), grid.gnnsn1, grid.node_fill_level, i );
	   
		  // recompute the number of neighbors.
		  num_neighbors = grid.gnnsn[ grid.gnnsn1[i+1] ] - grid.gnnsn[ ( grid.gnnsn1[i] + 1 ) ];

		  if ( num_neighbors >= num_neighbors_needed ) break;
		}

	      if ( j >= 3 )
		{
		  printf("FATAL ERROR: Node %d has attempted to add too many layers of neighbors. Halting...\n",i);
		  fflush(stdout);
		  exit(0);
		}

	      printf("Finished. Node %d now has %d neighbors in its stencil.\n",i,num_neighbors);
	    }
	      
	}

    }

  //if ( p.order == 4 || DEBUG )
  if ( 0 )
    {
      // Add first layer.
      printf("Adding first layer of neighbors........\n");
      for ( i=1; i <= grid.nn; i++ )
	{
	  Add_Layer_to_GNSN ( grid.nn, grid.num_elem, grid.c2n, grid.esn, grid.nnsn, grid.edges,
			      &(grid.gnsn), &(grid.gnnsn), grid.gnnsn1, grid.node_fill_level, i );
	}
      printf("Done\n");

      if ( DEBUG )
	{
	  printf("GNNSN1 DEBUG:\n");
	  for ( i=0; i <= (grid.nn+1); i++ )
	    {
	      printf("Node %d : %d\n",i,grid.gnnsn1[i]);
	    }

	  printf("GNNSN DEBUG:\n");
	  for ( i=0; i <= (grid.nn+1); i++ )
	    {
	      printf("Node %d : ",i);
	      sidx = grid.gnnsn1[i];
	      eidx = grid.gnnsn1[i+1];

	      for ( j=sidx; j < eidx; j++ )
		{
		  printf("%d ",grid.gnnsn[j]);
		}
	      printf("\n");
	    }

	  printf("GNSN DEBUG:\n");
	  for ( i=0; i <= (grid.nn+1); i++ )
	    {
	      printf("Node %d : ",i);
	      sidx = grid.gnnsn1[i];
	      eidx = grid.gnnsn1[i+1];

	      sidx = grid.gnnsn[sidx];
	      eidx = grid.gnnsn[eidx];
	      
	      for ( j=sidx; j < eidx; j++ )
		{
		  printf("%d ",grid.gnsn[j]);
		}
	      printf("\n");
	    }
	  fflush(stdout);
	}
      
      // Add second layer.
      printf("Adding second layer of neighbors........\n");
      for ( i=1; i <= grid.nn; i++ )
	{
	  Add_Layer_to_GNSN ( grid.nn, grid.num_elem, grid.c2n, grid.esn, grid.nnsn, grid.edges,
			      &(grid.gnsn), &(grid.gnnsn), grid.gnnsn1, grid.node_fill_level, i );
	}
      printf("Done\n");

      if ( DEBUG )
	{
	  printf("GNNSN1 DEBUG:\n");
	  for ( i=0; i <= (grid.nn+1); i++ )
	    {
	      printf("Node %d : %d\n",i,grid.gnnsn1[i]);
	    }

	  printf("GNNSN DEBUG:\n");
	  for ( i=0; i <= (grid.nn+1); i++ )
	    {
	      printf("Node %d : ",i);
	      sidx = grid.gnnsn1[i];
	      eidx = grid.gnnsn1[i+1];

	      for ( j=sidx; j < eidx; j++ )
		{
		  printf("%d ",grid.gnnsn[j]);
		}
	      printf("\n");
	    }

	  printf("GNSN DEBUG:\n");
	  for ( i=0; i <= (grid.nn+1); i++ )
	    {
	      printf("Node %d : ",i);
	      sidx = grid.gnnsn1[i];
	      eidx = grid.gnnsn1[i+1];

	      sidx = grid.gnnsn[sidx];
	      eidx = grid.gnnsn[eidx];
	      
	      for ( j=sidx; j < eidx; j++ )
		{
		  printf("%d ",grid.gnsn[j]);
		}
	      printf("\n");
	    }
	  fflush(stdout);
	}

      // Add third layer.
      printf("Adding third layer of neighbors........\n");
      for ( i=1; i <= grid.nn; i++ )
	{
	  Add_Layer_to_GNSN ( grid.nn, grid.num_elem, grid.c2n, grid.esn, grid.nnsn, grid.edges,
			      &(grid.gnsn), &(grid.gnnsn), grid.gnnsn1, grid.node_fill_level, i );
	}
      printf("Done\n");

      if ( DEBUG )
	{
	  printf("GNNSN1 DEBUG:\n");
	  for ( i=0; i <= (grid.nn+1); i++ )
	    {
	      printf("Node %d : %d\n",i,grid.gnnsn1[i]);
	    }

	  printf("GNNSN DEBUG:\n");
	  for ( i=0; i <= (grid.nn+1); i++ )
	    {
	      printf("Node %d : ",i);
	      sidx = grid.gnnsn1[i];
	      eidx = grid.gnnsn1[i+1];

	      for ( j=sidx; j < eidx; j++ )
		{
		  printf("%d ",grid.gnnsn[j]);
		}
	      printf("\n");
	    }

	  printf("GNSN DEBUG:\n");
	  for ( i=0; i <= (grid.nn+1); i++ )
	    {
	      printf("Node %d : ",i);
	      sidx = grid.gnnsn1[i];
	      eidx = grid.gnnsn1[i+1];

	      sidx = grid.gnnsn[sidx];
	      eidx = grid.gnnsn[eidx];
	      
	      for ( j=sidx; j < eidx; j++ )
		{
		  printf("%d ",grid.gnsn[j]);
		}
	      printf("\n");
	    }
	  fflush(stdout);
	}
    }
    
  // Get the Element Centroids.
  printf("Calculating the element centroids.......\n");
  Find_Element_Centroids(&grid);
  printf("Done.\n");

  // Build the edge/subedge data.
  printf("Building the edge/subedge data structures.......\n");
  Build_Edge_Data(&grid);
  printf("Done.\n");

  // Build the boundary edge data.
  printf("Building the boundary edge data structure.......\n");
  Build_Boundary_Edge_Data(&grid);
  printf("Done.\n");

  // Fix subedges.
  printf("Fixing subedges that touch curved boundaries.......\n");
  Fix_Subedge_Data(&grid);
  printf("Done.\n");

  // Build the Gaussian node data.
  printf("Building the Gaussian data for subedges and bedges.......\n");
  Generate_Gaussian_Data_For_Edge(&grid);
  printf("Done.\n");

  // Calculate the area of control volumes by Green's Theorem.
  cv_calc_area_Green(&grid);
  
  // Count the number of control volumes with negative displacement.
  RMS = 0.;
  for ( i=1, j=0; i <= grid.nn; i++ )
    {
      if ( grid.cv_area[i] <= 0. )
	{
	  j++;
	  printf("CV (node) %d has negative displacement - %.15e  %.15e\n",i,grid.x[2*i],grid.x[2*i+1]);
	}
      else
	{
	  RMS += grid.cv_area[i];
	}
    }
  
  printf("Detected %d CVs with negative displacement.\n",j);
  if ( j != 0 )
    printf(" Total positive area is %f\n                               or %.15e\n",(float)RMS,RMS);

  //sprintf(buff,"tecplot_area.dat");
  //write_tecplot_cv_area ( buff, &grid );

    // Check that all control volumes are closed, i.e. normal vectors sum to zero.
  printf("Checking for leaks........\n");
  cv_calc_check_closure(&grid);
  printf("Mesh is watertight.\n");
      
  // Build the length scale for the cv's.
  printf("Calculating the CV length scales.......\n");
  cv_calc_length_scale(&grid);
  printf("Done.\n");

  // Calculate the control volume moments.
  printf("Calculating the CV moments.......\n");
  cv_calc_Moments(&grid);
  printf("Done.\n");

  // Test the gradient.
  if ( DEBUG_GRAD )
    {
      Compute_Gradient_LeastSquares(&grid);
      sprintf(buff,"tecplot_gradient.dat");
      write_tecplot_derivs(buff,&grid);
    }
  if ( DEBUG_HESS )
    {
      Compute_Hessian(&grid, p);
      sprintf(buff,"tecplot_hessian.dat");
      write_tecplot_derivs(buff,&grid);
    }
  if ( DEBUG_DERIVS )
    {
      Compute_Derivatives ( &grid, p, 4 );
      sprintf(buff,"tecplot_derivs.dat");
      write_tecplot_derivs(buff,&grid);
    }

  // Initialize the problem.
  Build_CRS( &grid, p, &smem );

  if ( ANNULUS && 0 )
    {
      grid.bbc[0] = 0;
      grid.bbc[1] = 0;
      grid.bbc[2] = 14;
      grid.bbc[3] = 2;
      grid.bbc[4] = 16;
    }

  // CYLINDER HACK -
  //if ( CYLINDER && CYLINDER_HACK )
  if ( 1 ) 
    {
      for ( i=1; i <= grid.nn; i++ )
	{
	  if ( CYLINDER && CYLINDER_HACK )
	    {
	      if ( grid.x[i*NDIM+0]  > -1.0e-8 )
		{
		  grid.node_order[i] = p.order;
		}
	      else
		{
		  grid.node_order[i] = 1;
		}
	    }
	  else
	    {
	      grid.node_order[i] = p.order;
	    }
	}
      
      for ( i=1; i <= grid.nbedges; i++ )
	{
	  n =  grid.bedges[i*5+0];
	  gn = grid.bedges[i*5+1] + grid.nn;
	  
	  grid.node_order[gn] = grid.node_order[n];
	}
    }

  printf("Initializing flow field.\n");

  if ( p.isrestart )
    {
      printf("Restarting the solution. Using data from the restart file.\n");

      // Now I still need to copy the current solution over to the boundary nodes.

      printf("Copying the previous solution over to the ghost nodes.\n");

      for ( i=1; i <= grid.nbedges; i++ )
	{
	  // Grab the data from the structure.
	  n = grid.bedges[i*5+0];
	  j = grid.bedges[i*5+1] + grid.nn;
	  
	  // Copy the variables over.
	  grid.Q[j*NUM_VAR+0] = grid.Q[n*NUM_VAR+0];
	  grid.Q[j*NUM_VAR+1] = grid.Q[n*NUM_VAR+1];
	  grid.Q[j*NUM_VAR+2] = grid.Q[n*NUM_VAR+2];
	  grid.Q[j*NUM_VAR+3] = grid.Q[n*NUM_VAR+3];

	  grid.nQ[j*NUM_VAR+0] = grid.Q[n*NUM_VAR+0];
	  grid.nQ[j*NUM_VAR+1] = grid.Q[n*NUM_VAR+1];
	  grid.nQ[j*NUM_VAR+2] = grid.Q[n*NUM_VAR+2];
	  grid.nQ[j*NUM_VAR+3] = grid.Q[n*NUM_VAR+3];
	}
	  
    }
  else
    {
      if ( p.ic == 1 && p.p_or_cv == 1 )
	{
	  printf("Initializing the Vortex problem with control volume averages.\n");
	  init_cv( &grid, p );
	}
      else if ( p.ic == 2 && p.p_or_cv == 1 )
	{
	  printf("Initializing Yee's Vortex problem with control volume averages.\n");
	  init_cv( &grid, p );
	}
      else if ( p.ic == 3 && p.p_or_cv == 1 )
	init_cv( &grid, p );
      else if ( p.ic == 0 && p.p_or_cv == 1 )
	init_cv( &grid, p );
      else
	init_pv( &grid, p );
    }

  if ( ANNULUS && 0 )
    {
      grid.bbc[0] = 0;
      grid.bbc[1] = 0;
      grid.bbc[2] = 1;
      grid.bbc[3] = 2;
      grid.bbc[4] = 3;
    }

  //sprintf(buff,"tecplot_conserved_after_ic_init.dat");
  //write_tecplot_solutionC( buff, &grid );

  // Open the Iteration history files.
  sprintf(buffer,"residuals.dat");
  if ( p.isrestart )
    fp = fopen(buffer,"a");
  else
    fp = fopen(buffer,"w");

  sprintf(buffer,"time.dat");
  if ( p.isrestart )
    tfp = fopen(buffer,"a");
  else
    tfp = fopen(buffer,"w");

  // Clobber any existing extrema file.
  sprintf(buffer,"extrema.dat");
  tfp = fopen(buffer,"w");
  fclose(tfp);
  tfp = NULL;

  // Open a command file to control the solution procedure.
  sprintf(buffer,"command.dat");
  command = fopen(buffer,"w");
  fprintf(command,"#Post process now\n");
  fprintf(command,"0\n");
  fprintf(command,"#Terminate now\n");
  fprintf(command,"0\n");
  fprintf(command,"#Halt First Order iterations\n");
  fprintf(command,"0\n");
  fprintf(command,"#Halt Second Order iterations\n");
  fprintf(command,"0\n");
  fprintf(command,"#Halt Third Order iterations\n");
  fprintf(command,"0\n");
  fclose(command);

  command = NULL;

  // Declare memory to buffer the Residual RMS error from each iteration to print to a file later.
  rms = (double*)malloc( (p.nsteps+1)*sizeof(double));
  if ( rms == NULL ) { printf("MEMORY ERROR: Could not allocate 'rms'.\n"); exit(1); }

  sprintf(buff,"tecplot_conserved_ts_0.dat");
  write_tecplot_solutionC( buff, &grid );

  printf("\nStarting solution algorithm...\n\n");

  for ( i=0; i <= p.nsteps; i++ )
    {
      rms[i] = 1.E+10;
    }

  // Start the clock.
  gettimeofday(&start_time,NULL);

  if ( p.ic == 1 )                                               // Everything is in conserved so this is fine. At any rate this
    {                                                            // the first time step which is handled separately. Here grid->citer is 0.
      sprintf(buffer,"core_pressure.dat");
      if ( !(p.isrestart) )
	{
	  cpfp = fopen(buffer,"w");
	  get_core_pressure(&grid,p,&time,&cpressure);
	  fprintf(cpfp,"%.15E     %.15E\n",time,cpressure);
	  fclose(cpfp);
	}
    }

  if ( p.ic == 2 )
    {
      sprintf(buffer,"core_pressure.dat");
      if ( !(p.isrestart) )
	{
	  cpfp = fopen(buffer,"w");
	  get_core_pressure_yee(&grid,p,&time,&cpressure);
	  fprintf(cpfp,"%.15E     %.15E\n",time,cpressure);
	  fclose(cpfp);
	}
    }
  

  sprintf(buffer,"body_forces.dat");
  if ( p.isrestart )
    bforce = fopen(buffer,"a");
  else
    bforce = fopen(buffer,"w");

  // Close files for now.
  if ( fp != NULL )
    fclose(fp);
  if ( tfp != NULL )
    fclose(tfp);
  if ( bforce != NULL )
    fclose(bforce);

  fp = NULL;
  tfp = NULL;
  cpfp = NULL;
  bforce = NULL;

  i = grid.citer;

  if ( p.unsteady )
    {
      if ( i == p.start_post_pro ) // Write out the first file at the start of post processing.
	{
	  // Write out the unsteady file.
	  strcpy(buffer,pArgs[1]);
	  k = strlen(buffer);
	  
	  // the last five characters are ".mesh"
	  j = k - 5;
	  
	  if ( i < 10 )
	    sprintf(&(buffer[j]),"_0000%d_solution.dat",i);
	  else if ( i < 100 )
	    sprintf(&(buffer[j]),"_000%d_solution.dat",i);
	  else if ( i < 1000 )
	    sprintf(&(buffer[j]),"_00%d_solution.dat",i);
	  else if ( i < 10000 )
	    sprintf(&(buffer[j]),"_0%d_solution.dat",i);
	  else
	    sprintf(&(buffer[j]),"_%d_solution.dat",i);
	  
	  write_tecplot_solutionC(buffer, &grid );
	}
    }

  // Begin the main loop to solve the system.
  for ( i=(grid.citer + 1); i <= p.nsteps; i++ )
    {
      // Open files.
      sprintf(buffer,"residuals.dat");
      fp = fopen(buffer,"a");

      sprintf(buffer,"time.dat");
      tfp = fopen(buffer,"a");

      if ( p.ic == 1 || p.ic == 2 )
	{
	  sprintf(buffer,"core_pressure.dat");
	  cpfp = fopen(buffer,"a");
	}
      
      sprintf(buffer,"body_forces.dat");
      bforce = fopen(buffer,"a");

      // Compute the current solution time.
      grid.eltime = i*p.mintime;

      // Update the current time step for the grid object.
      grid.citer = i;
      grid.iter =  i;

      // Execute the solve routine.
      rms[i] = solve( &grid, &smem, p, i );

      if ( i%1 == 0 && i >= 1000 )
	{
	  //sprintf(buff,"tecplot_conserved_ts_%d.dat",i);
	  //write_tecplot_solutionC( buff, &grid );
	}
            
      fprintf(fp,"%d   %.15e\n",i,rms[i]);

      fclose(fp);
      fp = NULL;

      // Grab the time at the end of this iteration. This snippet of code was provided by Vince Betro
      // from notes given by Dr. Hyams.
      gettimeofday(&end_time,NULL);
      dt = (double)(end_time.tv_sec - start_time.tv_sec);
      dt += (double)(end_time.tv_usec - start_time.tv_usec)/1.0e+6;
      fprintf(tfp,"%.15e  %.15e\n", dt, rms[i]);

      fclose(tfp);
      tfp = NULL;
      
      //grid.nsteps++; These aren't needed since this is handled by the params file.
      //grid.iter++;

      printf("Finished Iteration %d - Residaul RMS error = %.15e\n",i,rms[i]);
      printf("Elapsed Time is = %f\n",(float)grid.eltime);
      fflush(stdout);

      // BIG OL'E HACK
      //if ( p.order == 4 && ( i > (p.foits+p.soits+p.toits) ) )
      //{
      //  if ( p.unsteady == 0 )
      //    {
	      //if ( (rms[i] < 1.0e-6) && ( (int)p.CFL <= 1000 ) )
		//p.CFL += 10.0;

	      //if ( rms[i] > rms[i-1] )
		//p.CFL -= 50.0;
      //    }
      //}
      
      if ( p.ic == 1 )
	{
	  get_core_pressure(&grid,p,&time,&cpressure);
	  fprintf(cpfp,"%.15E     %.15E\n",time,cpressure);
	  fclose(cpfp);
	  cpfp = NULL;
	}

      if ( p.ic == 2 )
	{
	  get_core_pressure_yee(&grid,p,&time,&cpressure);
	  fprintf(cpfp,"%.15E     %.15E\n",time,cpressure);
	  fclose(cpfp);
	  cpfp = NULL;
	}

      // Compute the body forces.
      Bforces[0] = 0.;  Bforces[1] = 0.;
      Compute_Body_Forces ( &grid, p, Bforces );

      // Write the body forces.
      if ( p.unsteady > 0 )
	fprintf(bforce,"%f  %.15E  %.15E\n",(float)grid.eltime,Bforces[0],Bforces[1]);
      else
	fprintf(bforce,"%f  %.15E  %.15E\n",(float)(i*0.01),Bforces[0],Bforces[1]);

      fclose(bforce);
      bforce = NULL;

      // Check the command file.
      sprintf(buffer,"command.dat");
      command = fopen(buffer,"r");

      if ( command == NULL )
	continue;

      fgets(buffer,buffsize,command);
      fgets(buffer,buffsize,command);
      sscanf(buffer,"%d",&(post_process_now));

      fgets(buffer,buffsize,command);
      fgets(buffer,buffsize,command);
      sscanf(buffer,"%d",&(terminate_now));

      fgets(buffer,buffsize,command);
      fgets(buffer,buffsize,command);
      sscanf(buffer,"%d",&(halt_foits));

      fgets(buffer,buffsize,command);
      fgets(buffer,buffsize,command);
      sscanf(buffer,"%d",&(halt_soits));

      fgets(buffer,buffsize,command);
      fgets(buffer,buffsize,command);
      sscanf(buffer,"%d",&(halt_toits));


      if ( post_process_now )
	{
	  // Write out the file.
	  strcpy(buffer,pArgs[1]);
	  k = strlen(buffer);

	  // the last five characters are ".mesh"
	  j = k - 5;

	  if ( i < 10 )
	    sprintf(&(buffer[j]),"_0000%d_solution.dat",i);
	  else if ( i < 100 )
	    sprintf(&(buffer[j]),"_000%d_solution.dat",i);
	  else if ( i < 1000 )
	    sprintf(&(buffer[j]),"_00%d_solution.dat",i);
	  else if ( i < 10000 )
	    sprintf(&(buffer[j]),"_0%d_solution.dat",i);
	  else
	    sprintf(&(buffer[j]),"_%d_solution.dat",i);

	  write_tecplot_solutionC(buffer, &grid );

	  strcpy(buffer,pArgs[1]);
	  k = strlen(buffer);

	  // the last five characters are ".mesh"
	  j = k - 5;

	  if ( i < 10 )
	    sprintf(&(buffer[j]),"_0000%d_restart.mesh",i);
	  else if ( i < 100 )
	    sprintf(&(buffer[j]),"_000%d_restart.mesh",i);
	  else if ( i < 1000 )
	    sprintf(&(buffer[j]),"_00%d_restart.mesh",i);
	  else if ( i < 10000 )
	    sprintf(&(buffer[j]),"_0%d_restart.mesh",i);
	  else
	    sprintf(&(buffer[j]),"_%d_restart.mesh",i);
	      
	  write_mesh(buffer, &grid );
	}

      if ( terminate_now )
	{
	  // Exit the time steps.
	  printf("COMMAND CONTROL: Recieved signal to terminate time steps.\n");
	  fflush(stdout);
	  break;
	}

      if ( halt_foits && (i>p.foits) )
	{
	  halt_foits = 0;
	  printf("COMMAND CONTROL: The solver is past first order iterations. No action will be taken.\n");
	  fflush(stdout);
	}

      if ( halt_soits && (i>(p.foits+p.soits)) )
	{
	  halt_soits = 0;
	  printf("COMMAND CONTROL: The solver is past second order iterations. No action will be taken.\n");
	  fflush(stdout);
	}

      if ( halt_soits && (i>(p.foits+p.soits+p.toits)) )
	{
	  halt_toits = 0;
	  printf("COMMAND CONTROL: The solver is past third order iterations. No action will be taken.\n");
	  fflush(stdout);
	}

      if ( halt_foits )
	{
	  // Skip the remaining first order iterations.
	  printf("COMMAND CONTROL: Halting further application of first order iterations.\n");
	  fflush(stdout);

	  // Reset the parameter controls.
	  p.foits = i;
	}

      if ( halt_soits )
	{
	  // Skip the remaining second order iterations.
	  printf("COMMAND CONTROL: Halting further application of second order iterations.\n");
	  fflush(stdout);

	  // Reset the parameter controls.
	  p.soits = (i - p.foits);

	  if ( p.soits < 0 )  p.soits = 0;
	}

      if ( halt_toits )
	{
	  // Skip the remaining second order iterations.
	  printf("COMMAND CONTROL: Halting further application of third order iterations.\n");
	  fflush(stdout);

	  // Reset the parameter controls.
	  p.toits = (i-p.foits-p.soits);
	  
	  if ( p.toits < 0 )  p.toits = 0;
	}
      
      fclose(command);
      command = NULL;
      sprintf(buffer,"command.dat");
      command = fopen(buffer,"w");
      fprintf(command,"#Post process now\n");
      fprintf(command,"0\n");
      fprintf(command,"#Terminate now\n");
      fprintf(command,"0\n");
      fprintf(command,"#Halt First Order iterations\n");
      fprintf(command,"0\n");
      fprintf(command,"#Halt Second Order iterations\n");
      fprintf(command,"0\n");
      fprintf(command,"#Halt Third Order iterations\n");
      fprintf(command,"0\n");
      fclose(command);
      
      command = NULL;
      
      // Now we can process unsteady files.
      if ( i >= p.start_post_pro && p.start_post_pro >= 0 )
	{
	  if ( p.post_pro_freq < 1 )
	    continue;

	  if ( (i % p.post_pro_freq) == 0 )
	    {
	      // Write out the unsteady file.
	      strcpy(buffer,pArgs[1]);
	      k = strlen(buffer);
	      
	      // the last five characters are ".mesh"
	      j = k - 5;
	      
	      if ( i < 10 )
		sprintf(&(buffer[j]),"_0000%d_solution.dat",i);
	      else if ( i < 100 )
		sprintf(&(buffer[j]),"_000%d_solution.dat",i);
	      else if ( i < 1000 )
		sprintf(&(buffer[j]),"_00%d_solution.dat",i);
	      else if ( i < 10000 )
		sprintf(&(buffer[j]),"_0%d_solution.dat",i);
	      else
		sprintf(&(buffer[j]),"_%d_solution.dat",i);
	      
	      write_tecplot_solutionC(buffer, &grid );
	      
	      // Write out a restart file.
	      strcpy(buffer,pArgs[1]);
	      k = strlen(buffer);
	      
	      // the last five characters are ".mesh"
	      j = k - 5;
	      
	      if ( i < 10 )
		sprintf(&(buffer[j]),"_0000%d_restart.mesh",i);
	      else if ( i < 100 )
		sprintf(&(buffer[j]),"_000%d_restart.mesh",i);
	      else if ( i < 1000 )
		sprintf(&(buffer[j]),"_00%d_restart.mesh",i);
	      else if ( i < 10000 )
		sprintf(&(buffer[j]),"_0%d_restart.mesh",i);
	      else
		sprintf(&(buffer[j]),"_%d_restart.mesh",i);
	      
	      write_mesh(buffer, &grid );
	      
	    }
	  else if ( i == p.start_post_pro ) // Write out the first file at the start of post processing.
	    {
	      // Write out the unsteady file.
	      strcpy(buffer,pArgs[1]);
	      k = strlen(buffer);
	      
	      // the last five characters are ".mesh"
	      j = k - 5;

	      if ( i < 10 )
		sprintf(&(buffer[j]),"_0000%d_solution.dat",i);
	      else if ( i < 100 )
		sprintf(&(buffer[j]),"_000%d_solution.dat",i);
	      else if ( i < 1000 )
		sprintf(&(buffer[j]),"_00%d_solution.dat",i);
	      else if ( i < 10000 )
		sprintf(&(buffer[j]),"_0%d_solution.dat",i);
	      else
		sprintf(&(buffer[j]),"_%d_solution.dat",i);
	      
	      write_tecplot_solutionC(buffer, &grid );

	      // Write out a restart file.
	      strcpy(buffer,pArgs[1]);
	      k = strlen(buffer);

	      // the last five characters are ".mesh"
	      j = k - 5;

	      if ( i < 10 )
		sprintf(&(buffer[j]),"_0000%d_restart.mesh",i);
	      else if ( i < 100 )
		sprintf(&(buffer[j]),"_000%d_restart.mesh",i);
	      else if ( i < 1000 )
		sprintf(&(buffer[j]),"_00%d_restart.mesh",i);
	      else if ( i < 10000 )
		sprintf(&(buffer[j]),"_0%d_restart.mesh",i);
	      else
		sprintf(&(buffer[j]),"_%d_restart.mesh",i);
	      
	      write_mesh(buffer, &grid );
	    }
	}

      //if ( rms[i] < 1.E-15 ) break;
      
    }                                                // End of iterations.
  
  if ( DEBUG ) return 0;
  
  //fclose(fp);
  //fclose(tfp);

  if ( p.ic == 1 || p.ic == 2 )
    {
      //fclose(cpfp);
    }
  
  // hack
  //if ( 0 )  -- Why did I force this out?

  // I think I had this commented out so that it wouldn't do the centerline plots for a restart. Now I'll just check that it did at least
  // one time step so that the derivatives are defined since they are needed for reconstruction calls.

  if ( p.nsteps > 0 )
    {
      if ( p.ic == 1 )
	write_vortex_pressure_centerline ( &grid, p, &grid.eltime );
      
      if ( p.ic == 2 )
	write_vortex_pressure_centerline_yee ( &grid, p, &grid.eltime );
    }

  
  sprintf(buff,"tecplot_solution_conserved.dat");
  write_tecplot_solutionC( buff, &grid );

  sprintf(buff,"tecplot_solution_primitive.dat");
  write_tecplot_solutionP( buff, &grid, p.gamma );

  write_surface_cp( &grid, p );

  // Write out the solution.
  strcpy(buffer,pArgs[1]);
  i = strlen(buffer);

  // the last five characters are ".mesh"

  j = i - 5;

  strcpy(&(buffer[j]),"_solution.mesh");

  //printf("%s\n",buffer);

  write_mesh( buffer, &grid );


  if ( NACA )
    {
      compute_integral_error ( &grid , p );
    }

  if ( CYLINDER )
    {
      double Pinf, Ptinf, Pt, vmag, clocal, mlocal;
      double ui, vi, rhoi, Ei, Pi;

      Ptinf = 1.0 / p.gamma * ( pow( ( 1 + 0.5*(p.gamma-1.)*p.mach*p.mach) , (p.gamma/(p.gamma-1.)) ) );
      Pinf = 1. / p.gamma ;

      double tA = 0.;
      double L1;
      double L2;
      double Linf;
      double error = 0.;

      L1 = 0.;
      L2 = 0.;
      Linf = 0.;

      for ( n=1; n <= grid.nn; n++ )
	{
	  tA += grid.cv_area[n];
	}

      // compute the norm of the error Q - Qe.

      for ( i=1; i <= grid.nn; i++ )
	{
	  rhoi = grid.Q[i*NUM_VAR+0];
	  ui = grid.Q[i*NUM_VAR+1]/rhoi;
	  vi = grid.Q[i*NUM_VAR+2]/rhoi;
	  Ei = grid.Q[i*NUM_VAR+3];
	  
	  Pi = (p.gamma-1.)*(Ei - 0.5*rhoi*(ui*ui + vi*vi));
	  
	  vmag = sqrt( ui*ui + vi*vi );
	  clocal = sqrt( p.gamma * Pi / rhoi );
	  mlocal = vmag / clocal;
	  
	  Pt = Pi * ( pow( ( 1. + 0.5*(p.gamma-1.)*mlocal*mlocal) , (p.gamma/(p.gamma-1.)) ) );
	  error = fabs( Pt/Ptinf - 1. );
	  
	  L1 += error;
	  L2 += error*error;
	  Linf = MAX( error , Linf );
	}
      
      L1 = L1 / ( (double)grid.nn);
      L2 = L2 / ( (double)grid.nn);
      L2 = sqrt(L2);
	  
      printf("\n\n");
      printf("Error analysis. Calculated as |Pt/Ptinf-1| / (nn)\n");
      printf("Error in total pressure:\n");
      printf("  L1 norm = %.15e\n",L1);
      printf("  L2 norm = %.15e\n",L2);
      printf("  L-infinity norm = %.15e\n\n",Linf);
      printf("\n");
      printf("%.15e\n",L1);
      printf("%.15e\n",L2);
      printf("%.15e\n",Linf);

      for ( i=1; i <= grid.nn; i++ )
	{
	  rhoi = grid.Q[i*NUM_VAR+0];
	  ui = grid.Q[i*NUM_VAR+1]/rhoi;
	  vi = grid.Q[i*NUM_VAR+2]/rhoi;
	  Ei = grid.Q[i*NUM_VAR+3];
	  
	  Pi = (p.gamma-1.)*(Ei - 0.5*rhoi*(ui*ui + vi*vi));
	  
	  vmag = sqrt( ui*ui + vi*vi );
	  clocal = sqrt( p.gamma * Pi / rhoi );
	  mlocal = vmag / clocal;
	  
	  Pt = Pi * ( pow( ( 1. + 0.5*(p.gamma-1.)*mlocal*mlocal) , (p.gamma/(p.gamma-1.)) ) );
	  error = fabs( Pt/Ptinf - 1. );
	  
	  L1 += error * grid.cv_area[i];
	  L2 += error*error * grid.cv_area[i];
	  Linf = MAX( error , Linf );
	}

      L1 = L1 / (tA);
      L2 = L2 / (tA);
      L2 = sqrt(L2);
	  
      printf("\n\n");
      printf("Error analysis. Calculated as ( |Pt/Ptinf-1|*dA / tA )\n");
      printf("Error in total pressure:\n");
      printf("  L1 norm = %.15e\n",L1);
      printf("  L2 norm = %.15e\n",L2);
      printf("  L-infinity norm = %.15e\n\n",Linf);
      printf("\n");
      printf("%.15e\n",L1);
      printf("%.15e\n",L2);
      printf("%.15e\n",Linf);

       compute_integral_error ( &grid , p );
    }

  // VORTEX CONVECTION
  if ( p.ic == 1 )  // Vortex convection problem.
    {

      if ( p.isrestart == 1 && p.nsteps == 0 )  // recompute derivatives if this is a restart.
	{
	  // for ( i=1; i <= (grid.nn + grid.nn_ghost); i++ )      // mesh_io handles this.
 	    //{
	  //for ( j=0; j < NUM_VAR; j++ )
	  //{
	  //  grid.nQ[i*NUM_VAR+j] = grid.Q[i*NUM_VAR+j];
	  //}
	  //}

	  // Set the limiters to 1.
	  for ( i=1; i <= grid.nn; i++ )
	    {
	      for ( j=0; j < 8; j++ )
		grid.phi[i*8+j] = 1.;
	    }
	  
	  if ( (p.order == 1 || p.order == 2) )
	    Compute_Gradient_LeastSquares(&grid);
	  else if ( p.order == 3 )
	    {
	      Compute_Hessian(&grid,p);
	    }
	  else
	    Compute_Derivatives(&grid,p,4);

	  // Get the new values for the limiters if desired.
	  if ( p.limit == 1 )
	    compute_limiter_barth ( &grid, p );

	  else if ( p.limit == 2 )
	    compute_limiter_venkatakrishnan ( &grid, p );

	  else if ( p.limit == 3 )
	    compute_limiter_venkatakrishnan2 ( &grid, p );

	  else if ( p.limit == 4 )
	    compute_limiter_gooch1 ( &grid, p );

	  else if ( p.limit == 5 )
	    compute_limiter_gooch2 ( &grid, p );
	}

      double tA = 0.;
      double L1[4];
      double L2[4];
      double Linf[4];
      double error = 0.;

      // Initialize the analytical solution.
      QexC = (double*)malloc((grid.nn + 1)*NUM_VAR*sizeof(double));
      if ( QexC == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'QexC' in main.C.\n"); exit(1); }

      QexP = (double*)malloc((grid.nn + 1)*NUM_VAR*sizeof(double));
      if ( QexP == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'QexP' in main.C.\n"); exit(1); }

      initialize_vortex_cv ( p, &grid, QexC );

      initialize_vortex_cv ( p, &grid, QexP );

      // Convert QexP to primitives.
      for ( i=1; i <= grid.nn; i++ )
	{
	  ConvertConservedToPrimitive ( p.gamma , &(QexP[i*NUM_VAR]) );
	}

      // Write out the exact solutions. Copy the exact solution to Q after Q is copied to pQ.
      for ( i=1; i <= grid.nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      grid.pQ[i*NUM_VAR+j] = grid.Q[i*NUM_VAR+j];
	    }
	}

      for ( i=1; i <= grid.nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      grid.Q[i*NUM_VAR+j] = QexC[i*NUM_VAR+j];
	    }
	}

      sprintf(buff,"vortex_exact_conserved.dat");
      write_tecplot_solutionC( buff, &grid );

      sprintf(buff,"vortex_exact_primitive.dat");
      write_tecplot_solutionP( buff, &grid, p.gamma );

      // Now copy the numerical solution back to grid.Q
      for ( i=1; i <= grid.nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      grid.Q[i*NUM_VAR+j] = grid.pQ[i*NUM_VAR+j];
	    }
	}

      // Get error norms.
      for ( i=0; i < NUM_VAR; i++ )
	{
	  L1[i] = 0.;
	  L2[i] = 0.;
	  Linf[i] = 0.;
	}

      for ( n=1; n <= grid.nn; n++ )
	{
	  tA += grid.cv_area[n];
	}

      // compute the norm of the error Q - Qe.

      for ( n=1; n <= grid.nn; n++ )
	{
	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      error = fabs( grid.Q[n*NUM_VAR+i] - QexC[n*NUM_VAR+i] );
	
	      L1[i] += error;
	      L2[i] += error*error;
	      Linf[i] = MAX( error , Linf[i] );
	    }
	}

      for ( i=0; i < NUM_VAR; i++ )
	{
	  L1[i] = L1[i] / ( (double)grid.nn );
	  L2[i] = L2[i] / ( (double)grid.nn );
	  L2[i] = sqrt( L2[i] );
	}
	  
      printf("\n\n");
      printf("Error analysis. Calculated as |Qn - Qe| / (nn)\n");
      printf("Density:\n");
      printf("  L1 norm = %.15e\n",L1[0]);
      printf("  L2 norm = %.15e\n",L2[0]);
      printf("  L-infinity norm = %.15e\n\n",Linf[0]);

      printf("X momentum:\n");
      printf("  L1 norm = %.15e\n",L1[1]);
      printf("  L2 norm = %.15e\n",L2[1]);
      printf("  L-infinity norm = %.15e\n\n",Linf[1]);

      printf("Y momentum:\n");
      printf("  L1 norm = %.15e\n",L1[2]);
      printf("  L2 norm = %.15e\n",L2[2]);
      printf("  L-infinity norm = %.15e\n\n",Linf[2]);

      printf("Total Energy:\n");
      printf("  L1 norm = %.15e\n",L1[3]);
      printf("  L2 norm = %.15e\n",L2[3]);
      printf("  L-infinity norm = %.15e\n\n",Linf[3]);

      printf("\n");
      printf("%.15e\n",L1[0]);
      printf("%.15e\n",L2[0]);
      printf("%.15e\n",Linf[0]);
      printf("%.15e\n",L1[1]);
      printf("%.15e\n",L2[1]);
      printf("%.15e\n",Linf[1]);
      printf("%.15e\n",L1[2]);
      printf("%.15e\n",L2[2]);
      printf("%.15e\n",Linf[2]);
      printf("%.15e\n",L1[3]);
      printf("%.15e\n",L2[3]);
      printf("%.15e\n",Linf[3]);

      for ( i=0; i < NUM_VAR; i++ )
	{
	  L1[i] = 0.;
	  L2[i] = 0.;
	  Linf[i] = 0.;
	}

      for ( n=1; n <= grid.nn; n++ )
	{
	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      error = fabs( grid.Q[n*NUM_VAR+i] - QexC[n*NUM_VAR+i] );
	
	      L1[i] += error * grid.cv_area[n];
	      L2[i] += error*error * grid.cv_area[n];
	      Linf[i] = MAX( error , Linf[i] );
	    }
	}

      for ( i=0; i < NUM_VAR; i++ )
	{
	  L1[i] = L1[i] / tA;
	  L2[i] = L2[i] / tA;
	  L2[i] = sqrt( L2[i] );
	}
	  
      printf("\n\n");
      printf("Error analysis. Calculated as ( |Qn - Qe|*dA / tA )\n");
      printf("Density:\n");
      printf("  L1 norm = %.15e\n",L1[0]);
      printf("  L2 norm = %.15e\n",L2[0]);
      printf("  L-infinity norm = %.15e\n\n",Linf[0]);

      printf("X momentum:\n");
      printf("  L1 norm = %.15e\n",L1[1]);
      printf("  L2 norm = %.15e\n",L2[1]);
      printf("  L-infinity norm = %.15e\n\n",Linf[1]);

      printf("Y momentum:\n");
      printf("  L1 norm = %.15e\n",L1[2]);
      printf("  L2 norm = %.15e\n",L2[2]);
      printf("  L-infinity norm = %.15e\n\n",Linf[2]);

      printf("Total Energy:\n");
      printf("  L1 norm = %.15e\n",L1[3]);
      printf("  L2 norm = %.15e\n",L2[3]);
      printf("  L-infinity norm = %.15e\n\n",Linf[3]);

      printf("\n");
      printf("%.15e\n",L1[0]);
      printf("%.15e\n",L2[0]);
      printf("%.15e\n",Linf[0]);
      printf("%.15e\n",L1[1]);
      printf("%.15e\n",L2[1]);
      printf("%.15e\n",Linf[1]);
      printf("%.15e\n",L1[2]);
      printf("%.15e\n",L2[2]);
      printf("%.15e\n",Linf[2]);
      printf("%.15e\n",L1[3]);
      printf("%.15e\n",L2[3]);
      printf("%.15e\n",Linf[3]);

      if ( RECON_PRIM )              // convert nQ, which is used in the reconstruction, to primitives if needed.
	{
	  for ( i=1; i <= grid.nn; i++ )
	    {
	      ConvertConservedToPrimitive ( p.gamma , &(grid.nQ[i*NUM_VAR]) );
	    }
	}
      
      compute_integral_error ( &grid , p );

      // Now, below I am expecting grid->nQ to be primitive. If set above, all is good. If not, set it now.
      if ( !RECON_PRIM )
	{
	  for ( i=1; i <= grid.nn; i++ )
	    {
	      ConvertConservedToPrimitive ( p.gamma , &(grid.nQ[i*NUM_VAR]) );
	    }
	}

      // Compute the difference between numerical and exact solutions.
      for ( i=1; i <= grid.nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      grid.Q[i*NUM_VAR+j]  = fabs( QexC[i*NUM_VAR+j] - grid.Q[i*NUM_VAR+j] );           // conserved variables.
	      grid.nQ[i*NUM_VAR+j] = fabs( QexP[i*NUM_VAR+j] - grid.nQ[i*NUM_VAR+j] );          // primitive variables.
	    }
	}

      sprintf(buff,"vortex_diff_num_exact_cons.dat");
      write_tecplot_solutionC( buff, &grid );

      // Need to copy nQ to Q. First copy Q to pQ.
      for ( i=1; i <= grid.nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      grid.pQ[i*NUM_VAR+j] = grid.Q[i*NUM_VAR+j];
	    }
	}

      for ( i=1; i <= grid.nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      grid.Q[i*NUM_VAR+j] = grid.nQ[i*NUM_VAR+j];
	    }
	}

      sprintf(buff,"vortex_diff_num_exact_prim.dat");
      write_tecplot_solutionC( buff, &grid );  // It will call them conserved names but they are primitives.

      // Now copy pQ back to Q
      for ( i=1; i <= grid.nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      grid.Q[i*NUM_VAR+j] = grid.pQ[i*NUM_VAR+j];
	    }
	}

      freenull(QexC);
      freenull(QexP);
    }
  else if ( p.ic == 2 )
    {

      if ( p.isrestart == 1 && p.nsteps == 0 )  // recompute derivatives if this is a restart.
	{
	  // for ( i=1; i <= (grid.nn + grid.nn_ghost); i++ )       // mesh_io handles this.
	  //{
	  //  for ( j=0; j < NUM_VAR; j++ )
	  //{
	  //  grid.nQ[i*NUM_VAR+j] = grid.Q[i*NUM_VAR+j];
	  //}
	  //}

	  // Set the limiters to 1.
	  for ( i=1; i <= grid.nn; i++ )
	    {
	      for ( j=0; j < 8; j++ )
		grid.phi[i*8+j] = 1.;
	    }
	  
	  if ( (p.order == 1 || p.order == 2) )
	    Compute_Gradient_LeastSquares(&grid);
	  else if ( p.order == 3 )
	    {
	      Compute_Hessian(&grid,p);
	    }
	  else
	    Compute_Derivatives(&grid,p,4);

	  // Get the new values for the limiters if desired.
	  if ( p.limit == 1 )
	    compute_limiter_barth ( &grid, p );

	  else if ( p.limit == 2 )
	    compute_limiter_venkatakrishnan ( &grid, p );

	  else if ( p.limit == 3 )
	    compute_limiter_venkatakrishnan2 ( &grid, p );

	  else if ( p.limit == 4 )
	    compute_limiter_gooch1 ( &grid, p );

	  else if ( p.limit == 5 )
	    compute_limiter_gooch2 ( &grid, p );
	}

      double tA = 0.;
      double L1[4];
      double L2[4];
      double Linf[4];
      double error = 0.;

      // Initialize the analytical solution.
      QexC = (double*)malloc((grid.nn + 1)*NUM_VAR*sizeof(double));
      if ( QexC == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'QexC' in main.C.\n"); exit(1); }

      QexP = (double*)malloc((grid.nn + 1)*NUM_VAR*sizeof(double));
      if ( QexP == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'QexP' in main.C.\n"); exit(1); }

      initialize_vortex_yee_cv ( p, &grid, QexC );
      initialize_vortex_yee_cv ( p, &grid, QexP );

      // Convert QexP to primitives.
      for ( i=1; i <= grid.nn; i++ )
	{
	  ConvertConservedToPrimitive ( p.gamma , &(QexP[i*NUM_VAR]) );
	}

      // Write out the exact solutions. Copy the exact solution to Q after Q is copied to pQ.
      for ( i=1; i <= grid.nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      grid.pQ[i*NUM_VAR+j] = grid.Q[i*NUM_VAR+j];
	    }
	}

      for ( i=1; i <= grid.nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      grid.Q[i*NUM_VAR+j] = QexC[i*NUM_VAR+j];
	    }
	}

      sprintf(buff,"vortex_yee_exact_conserved.dat");
      write_tecplot_solutionC( buff, &grid );

      sprintf(buff,"vortex_yee_exact_primitive.dat");
      write_tecplot_solutionP( buff, &grid, p.gamma );

      // Now copy the numerical solution back to grid.Q
      for ( i=1; i <= grid.nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      grid.Q[i*NUM_VAR+j] = grid.pQ[i*NUM_VAR+j];
	    }
	}

      // Get error norms.
      for ( i=0; i < NUM_VAR; i++ )
	{
	  L1[i] = 0.;
	  L2[i] = 0.;
	  Linf[i] = 0.;
	}

      for ( n=1; n <= grid.nn; n++ )
	{
	  tA += grid.cv_area[n];
	}

      // compute the norm of the error Q - Qe.
      for ( n=1; n <= grid.nn; n++ )
	{
	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      error = fabs( grid.Q[n*NUM_VAR+i] -QexC[n*NUM_VAR+i] );
	
	      L1[i] += error;
	      L2[i] += error*error;
	      Linf[i] = MAX( error , Linf[i] );
	    }
	}

      for ( i=0; i < NUM_VAR; i++ )
	{
	  L1[i] = L1[i] / ( (double)grid.nn );
	  L2[i] = L2[i] / ( (double)grid.nn );
	  L2[i] = sqrt( L2[i] );
	}
	  
      printf("\n\n");
      printf("Error analysis. Calculated as |Qn - Qe| / (nn)\n");
      printf("Density:\n");
      printf("  L1 norm = %.15e\n",L1[0]);
      printf("  L2 norm = %.15e\n",L2[0]);
      printf("  L-infinity norm = %.15e\n\n",Linf[0]);

      printf("X momentum:\n");
      printf("  L1 norm = %.15e\n",L1[1]);
      printf("  L2 norm = %.15e\n",L2[1]);
      printf("  L-infinity norm = %.15e\n\n",Linf[1]);

      printf("Y momentum:\n");
      printf("  L1 norm = %.15e\n",L1[2]);
      printf("  L2 norm = %.15e\n",L2[2]);
      printf("  L-infinity norm = %.15e\n\n",Linf[2]);

      printf("Total Energy:\n");
      printf("  L1 norm = %.15e\n",L1[3]);
      printf("  L2 norm = %.15e\n",L2[3]);
      printf("  L-infinity norm = %.15e\n\n",Linf[3]);

      printf("\n");
      printf("%.15e\n",L1[0]);
      printf("%.15e\n",L2[0]);
      printf("%.15e\n",Linf[0]);
      printf("%.15e\n",L1[1]);
      printf("%.15e\n",L2[1]);
      printf("%.15e\n",Linf[1]);
      printf("%.15e\n",L1[2]);
      printf("%.15e\n",L2[2]);
      printf("%.15e\n",Linf[2]);
      printf("%.15e\n",L1[3]);
      printf("%.15e\n",L2[3]);
      printf("%.15e\n",Linf[3]);

      for ( i=0; i < NUM_VAR; i++ )
	{
	  L1[i] = 0.;
	  L2[i] = 0.;
	  Linf[i] = 0.;
	}
      
      for ( n=1; n <= grid.nn; n++ )
	{
	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      error = fabs( grid.Q[n*NUM_VAR+i] - QexC[n*NUM_VAR+i] );
	
	      L1[i] += error * grid.cv_area[n];
	      L2[i] += error*error * grid.cv_area[n];
	      Linf[i] = MAX( error , Linf[i] );
	    }
	}

      for ( i=0; i < NUM_VAR; i++ )
	{
	  L1[i] = L1[i] / tA;
	  L2[i] = L2[i] / tA;
	  L2[i] = sqrt( L2[i] );
	}
	  
      printf("\n\n");
      printf("Error analysis. Calculated as ( |Qn - Qe|*dA / tA )\n");
      printf("Density:\n");
      printf("  L1 norm = %.15e\n",L1[0]);
      printf("  L2 norm = %.15e\n",L2[0]);
      printf("  L-infinity norm = %.15e\n\n",Linf[0]);

      printf("X momentum:\n");
      printf("  L1 norm = %.15e\n",L1[1]);
      printf("  L2 norm = %.15e\n",L2[1]);
      printf("  L-infinity norm = %.15e\n\n",Linf[1]);

      printf("Y momentum:\n");
      printf("  L1 norm = %.15e\n",L1[2]);
      printf("  L2 norm = %.15e\n",L2[2]);
      printf("  L-infinity norm = %.15e\n\n",Linf[2]);

      printf("Total Energy:\n");
      printf("  L1 norm = %.15e\n",L1[3]);
      printf("  L2 norm = %.15e\n",L2[3]);
      printf("  L-infinity norm = %.15e\n\n",Linf[3]);

      printf("\n");
      printf("%.15e\n",L1[0]);
      printf("%.15e\n",L2[0]);
      printf("%.15e\n",Linf[0]);
      printf("%.15e\n",L1[1]);
      printf("%.15e\n",L2[1]);
      printf("%.15e\n",Linf[1]);
      printf("%.15e\n",L1[2]);
      printf("%.15e\n",L2[2]);
      printf("%.15e\n",Linf[2]);
      printf("%.15e\n",L1[3]);
      printf("%.15e\n",L2[3]);
      printf("%.15e\n",Linf[3]);

      if ( RECON_PRIM )              // convert nQ, which is used in the reconstruction, to primitives if needed.
	{
	  for ( i=1; i <= grid.nn; i++ )
	    {
	      ConvertConservedToPrimitive ( p.gamma , &(grid.nQ[i*NUM_VAR]) );
	    }
	}

      compute_integral_error( &grid , p );

      // Now, below I am expecting grid->nQ to be primitive. If set above, all is good. If not, set it now.
      if ( !RECON_PRIM )
	{
	  for ( i=1; i <= grid.nn; i++ )
	    {
	      ConvertConservedToPrimitive ( p.gamma , &(grid.nQ[i*NUM_VAR]) );
	    }
	}
      
      // Compute the difference between numerical and exact solutions.
      for ( i=1; i <= grid.nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      grid.Q[i*NUM_VAR+j]  = fabs( QexC[i*NUM_VAR+j] - grid.Q[i*NUM_VAR+j] );
	      grid.nQ[i*NUM_VAR+j] = fabs( QexP[i*NUM_VAR+j] - grid.nQ[i*NUM_VAR+j] );
	    }
	}

      sprintf(buff,"vortex_yee_diff_num_exact_cons.dat");
      write_tecplot_solutionC( buff, &grid );

      // Need to copy nQ to Q. First copy Q to pQ.
      for ( i=1; i <= grid.nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      grid.pQ[i*NUM_VAR+j] = grid.Q[i*NUM_VAR+j];
	    }
	}

      for ( i=1; i <= grid.nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      grid.Q[i*NUM_VAR+j] = grid.nQ[i*NUM_VAR+j];
	    }
	}

      sprintf(buff,"vortex_yee_diff_num_exact_prim.dat");
      write_tecplot_solutionC( buff, &grid );  // It will call them conserved names but they are primitives.

      // Now copy pQ back to Q
      for ( i=1; i <= grid.nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      grid.Q[i*NUM_VAR+j] = grid.pQ[i*NUM_VAR+j];
	    }
	}

      freenull(QexC);
      freenull(QexP);
    }

  if ( ANNULUS && !MMS )
    {
      double tA = 0.;
      double L1[4];
      double L2[4];
      double Linf[4];
      double error = 0.;

      for ( i=0; i < NUM_VAR; i++ )
	{
	  L1[i] = 0.;
	  L2[i] = 0.;
	  Linf[i] = 0.;
	}

      // If I was reconstructing the primitive variables, then Qe stores the primitives so I will need to
      // convert grid->Q to primitives as well.
      if ( RECON_PRIM && p.ic != 3 )  // ic_3 uses primitive variables.
	{
	  for ( i=1; i <= grid.nn; i++ )
	    {
	      ConvertConservedToPrimitive ( p.gamma , &(grid.Q[i*NUM_VAR]) );
	      ConvertConservedToPrimitive ( p.gamma , &(grid.nQ[i*NUM_VAR]) );  // for use in reconstruction.
	    }
	}

      // compute the norm of the error Q - Qe.

      for ( n=1; n <= grid.nn; n++ )
	{
	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      error = fabs( grid.Q[n*NUM_VAR+i] - grid.Qe[n*NUM_VAR+i] );
	
	      L1[i] += error;
	      L2[i] += error*error;
	      Linf[i] = MAX( error , Linf[i] );
	    }
	}

      for ( i=0; i < NUM_VAR; i++ )
	{
	  L1[i] = L1[i] / ( (double)grid.nn );
	  L2[i] = L2[i] / ( (double)grid.nn );
	  L2[i] = sqrt( L2[i] );
	}
	  
      printf("\n\n");
      printf("Error analysis. Calculated as |Qn - Qe| / (nn)\n");
      printf("Density:\n");
      printf("  L1 norm = %.15e\n",L1[0]);
      printf("  L2 norm = %.15e\n",L2[0]);
      printf("  L-infinity norm = %.15e\n\n",Linf[0]);

      if ( RECON_PRIM )
	printf("X velocity:\n");
      else
	printf("X momentum:\n");
      printf("  L1 norm = %.15e\n",L1[1]);
      printf("  L2 norm = %.15e\n",L2[1]);
      printf("  L-infinity norm = %.15e\n\n",Linf[1]);

      if ( RECON_PRIM )
	printf("Y velocity:\n");
      else
	printf("Y momentum:\n");
      printf("  L1 norm = %.15e\n",L1[2]);
      printf("  L2 norm = %.15e\n",L2[2]);
      printf("  L-infinity norm = %.15e\n\n",Linf[2]);

      if ( RECON_PRIM)
	printf("Pressure:\n");
      else
	printf("Total Energy:\n");
      printf("  L1 norm = %.15e\n",L1[3]);
      printf("  L2 norm = %.15e\n",L2[3]);
      printf("  L-infinity norm = %.15e\n\n",Linf[3]);

      printf("\n");
      printf("%.15e\n",L1[0]);
      printf("%.15e\n",L2[0]);
      printf("%.15e\n",Linf[0]);
      printf("%.15e\n",L1[1]);
      printf("%.15e\n",L2[1]);
      printf("%.15e\n",Linf[1]);
      printf("%.15e\n",L1[2]);
      printf("%.15e\n",L2[2]);
      printf("%.15e\n",Linf[2]);
      printf("%.15e\n",L1[3]);
      printf("%.15e\n",L2[3]);
      printf("%.15e\n",Linf[3]);

      
      for ( i=0; i < NUM_VAR; i++ )
	{
	  L1[i] = 0.;
	  L2[i] = 0.;
	  Linf[i] = 0.;
	}

      for ( n=1; n <= grid.nn; n++ )
	{
	  tA += grid.cv_area[n];
	}

      // compute the norm of the error Q - Qe.

      for ( n=1; n <= grid.nn; n++ )
	{
	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      error = fabs( grid.Q[n*NUM_VAR+i] - grid.Qe[n*NUM_VAR+i] );
	
	      L1[i] += error * grid.cv_area[n];
	      L2[i] += error*error * grid.cv_area[n];
	      Linf[i] = MAX( error , Linf[i] );
	    }
	}

      for ( i=0; i < NUM_VAR; i++ )
	{
	  L1[i] = L1[i] / tA;
	  //L1[i] = L1[i] / ((double)grid.nn);
	  L2[i] = L2[i] / tA;
	  //L2[i] = L2[i] / ((double)grid.nn);
	  L2[i] = sqrt( L2[i] );
	}
	  
      printf("\n\n");
      printf("Error analysis. Calculated as ( |Qn - Qe|^p *dA / tA )^1/p\n");
      printf("Density:\n");
      printf("  L1 norm = %.15e\n",L1[0]);
      printf("  L2 norm = %.15e\n",L2[0]);
      printf("  L-infinity norm = %.15e\n\n",Linf[0]);

      if ( RECON_PRIM )
	printf("X velocity:\n");
      else
	printf("X momentum:\n");
      printf("  L1 norm = %.15e\n",L1[1]);
      printf("  L2 norm = %.15e\n",L2[1]);
      printf("  L-infinity norm = %.15e\n\n",Linf[1]);

      if ( RECON_PRIM )
	printf("Y velocity:\n");
      else
	printf("Y momentum:\n");
      printf("  L1 norm = %.15e\n",L1[2]);
      printf("  L2 norm = %.15e\n",L2[2]);
      printf("  L-infinity norm = %.15e\n\n",Linf[2]);

      if ( RECON_PRIM)
	printf("Pressure:\n");
      else
	printf("Total Energy:\n");
      printf("  L1 norm = %.15e\n",L1[3]);
      printf("  L2 norm = %.15e\n",L2[3]);
      printf("  L-infinity norm = %.15e\n\n",Linf[3]);

      printf("\n");
      printf("%.15e\n",L1[0]);
      printf("%.15e\n",L2[0]);
      printf("%.15e\n",Linf[0]);
      printf("%.15e\n",L1[1]);
      printf("%.15e\n",L2[1]);
      printf("%.15e\n",Linf[1]);
      printf("%.15e\n",L1[2]);
      printf("%.15e\n",L2[2]);
      printf("%.15e\n",Linf[2]);
      printf("%.15e\n",L1[3]);
      printf("%.15e\n",L2[3]);
      printf("%.15e\n",Linf[3]);

      compute_integral_error ( &grid, p );
      

      // Compute the difference between numerical and exact solutions.
      for ( i=1; i <= grid.nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      grid.Q[i*NUM_VAR+j]  = fabs( grid.Qe[i*NUM_VAR+j] - grid.Q[i*NUM_VAR+j] );
	    }
	}

      if ( RECON_PRIM )
	sprintf(buff,"annulus_diff_num_exact_primitive.dat");
      else
	sprintf(buff,"annulus_diff_num_exact_cons.dat");
      write_tecplot_solutionC( buff, &grid );
    }

  // changed from not recon_prim
  //if ( MMS && !RECON_PRIM )    // manufactured solution error analysis.
  //if ( MMS || MMS_EXP )
  if ( MMS || MMS_EXP || DEBUG || DEBUG_HESS || DEBUG_DERIVS )
    {
      double tA = 0.;
      double L1[4];
      double L2[4];
      double Linf[4];
      double error = 0.;

      // If I was reconstructing the primitive variables, then Qe stores the primitives so I will need to
      // convert grid->Q to primitives as well.
      if ( RECON_PRIM && p.ic != 3 )  // ic_3 uses primitive variables.
	{
	  for ( i=1; i <= grid.nn; i++ )
	    {
	      ConvertConservedToPrimitive ( p.gamma , &(grid.Q[i*NUM_VAR]) );
	      ConvertConservedToPrimitive ( p.gamma , &(grid.nQ[i*NUM_VAR]) );  // for use in reconstruction.
	    }
	}
      
      for ( i=0; i < NUM_VAR; i++ )
	{
	  L1[i] = 0.;
	  L2[i] = 0.;
	  Linf[i] = 0.;
	}

      for ( n=1; n <= grid.nn; n++ )
	{
	  tA += grid.cv_area[n];
	}

      // compute the norm of the error Q - Qe.

      for ( n=1; n <= grid.nn; n++ )
	{
	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      error = fabs( grid.Q[n*NUM_VAR+i] - grid.Qe[n*NUM_VAR+i] );
	  
	      L1[i] += error * grid.cv_area[n];
	      L2[i] += error*error * grid.cv_area[n];
	      Linf[i] = MAX( error , Linf[i] );
	    }
	}

      for ( i=0; i < NUM_VAR; i++ )
	{
	  L1[i] = L1[i] / tA;
	  L2[i] = L2[i] / tA;
	  L2[i] = sqrt(L2[i]);
	}
	  
      printf("\n\n");
      printf("Error analysis for method of manufactured solutions. Calculated as ( |Qn - Qe|^p *dA / tA )^1/p\n");
      printf("Density:\n");
      printf("  L1 norm = %.15e\n",L1[0]);
      printf("  L2 norm = %.15e\n",L2[0]);
      printf("  L-infinity norm = %.15e\n\n",Linf[0]);

      if ( RECON_PRIM )
	printf("X velocity:\n");
      else
	printf("X momentum:\n");
      printf("  L1 norm = %.15e\n",L1[1]);
      printf("  L2 norm = %.15e\n",L2[1]);
      printf("  L-infinity norm = %.15e\n\n",Linf[1]);

      if ( RECON_PRIM )
	printf("Y velocity:\n");
      else
	printf("Y momentum:\n");
      printf("  L1 norm = %.15e\n",L1[2]);
      printf("  L2 norm = %.15e\n",L2[2]);
      printf("  L-infinity norm = %.15e\n\n",Linf[2]);

      if ( RECON_PRIM)
	printf("Pressure:\n");
      else
	printf("Total Energy:\n");
      printf("  L1 norm = %.15e\n",L1[3]);
      printf("  L2 norm = %.15e\n",L2[3]);
      printf("  L-infinity norm = %.15e\n\n",Linf[3]);

      printf("\n");
      printf("%.15e\n",L1[0]);
      printf("%.15e\n",L2[0]);
      printf("%.15e\n",Linf[0]);
      printf("%.15e\n",L1[1]);
      printf("%.15e\n",L2[1]);
      printf("%.15e\n",Linf[1]);
      printf("%.15e\n",L1[2]);
      printf("%.15e\n",L2[2]);
      printf("%.15e\n",Linf[2]);
      printf("%.15e\n",L1[3]);
      printf("%.15e\n",L2[3]);
      printf("%.15e\n",Linf[3]);


      // Redo this for the interior cv's only (look at node state).
      printf("Repeating the analysis for just the interior control volumes.\n");

      for ( i=0; i < NUM_VAR; i++ )
	{
	  L1[i] = 0.;
	  L2[i] = 0.;
	  Linf[i] = 0.;
	}

      tA = 0.;

      for ( n=1; n <= grid.nn; n++ )
	{
	  if ( grid.node_state[n] == INTERIOR )
	    tA += grid.cv_area[n];
	}

      // compute the norm of the error Q - Qe.

      for ( n=1; n <= grid.nn; n++ )
	{
	  for ( i=0; i < NUM_VAR; i++ )
	    {
	      if ( grid.node_state[n] == INTERIOR )
		{
		  error = fabs( grid.Q[n*NUM_VAR+i] - grid.Qe[n*NUM_VAR+i] );
		  
		  L1[i] += error * grid.cv_area[n];
		  L2[i] += error*error * grid.cv_area[n];
		  Linf[i] = MAX( error , Linf[i] );
		}
	    }
	}
      
      for ( i=0; i < NUM_VAR; i++ )
	{
	  L1[i] = L1[i] / tA;
	  L2[i] = L2[i] / tA;
	  L2[i] = sqrt(L2[i]);
	}
      
      printf("\n\n");
      printf("Error analysis for method of manufactured solutions. Calculated as ( |Qn - Qe|^p *dA / tA )^1/p\n");
      printf("Density:\n");
      printf("  L1 norm = %.15e\n",L1[0]);
      printf("  L2 norm = %.15e\n",L2[0]);
      printf("  L-infinity norm = %.15e\n\n",Linf[0]);

      if ( RECON_PRIM )
	printf("X velocity:\n");
      else
	printf("X momentum:\n");
      printf("  L1 norm = %.15e\n",L1[1]);
      printf("  L2 norm = %.15e\n",L2[1]);
      printf("  L-infinity norm = %.15e\n\n",Linf[1]);

      if ( RECON_PRIM )
	printf("Y velocity:\n");
      else
	printf("Y momentum:\n");
      printf("  L1 norm = %.15e\n",L1[2]);
      printf("  L2 norm = %.15e\n",L2[2]);
      printf("  L-infinity norm = %.15e\n\n",Linf[2]);

      if ( RECON_PRIM)
	printf("Pressure:\n");
      else
	printf("Total Energy:\n");
      printf("  L1 norm = %.15e\n",L1[3]);
      printf("  L2 norm = %.15e\n",L2[3]);
      printf("  L-infinity norm = %.15e\n\n",Linf[3]);

      printf("\n");
      printf("%.15e\n",L1[0]);
      printf("%.15e\n",L2[0]);
      printf("%.15e\n",Linf[0]);
      printf("%.15e\n",L1[1]);
      printf("%.15e\n",L2[1]);
      printf("%.15e\n",Linf[1]);
      printf("%.15e\n",L1[2]);
      printf("%.15e\n",L2[2]);
      printf("%.15e\n",Linf[2]);
      printf("%.15e\n",L1[3]);
      printf("%.15e\n",L2[3]);
      printf("%.15e\n",Linf[3]);

      compute_integral_error ( &grid , p );

      for ( i=1; i <= grid.nn; i++ )
	{
	  for ( j=0; j < NUM_VAR; j++ )
	    {
	      grid.Q[i*NUM_VAR+j] = fabs( grid.Q[i*NUM_VAR+j] - grid.Qe[i*NUM_VAR+j] );
	    }
	}

      if ( RECON_PRIM )
	sprintf(buff,"tecplot_exact_computed_diff_prim.dat");
      else
	sprintf(buff,"tecplot_exact_computed_diff_cons.dat");

      write_tecplot_solutionC ( buff, &grid );
    }

  // Inviscid boundary check.
  if ( (ANNULUS || CYLINDER) && 0 )
    {
      int real_node, ghost_node, b, bc;
      int s,nodeL,nodeR;
      double nx,ny,len,xi,yi,xmid,ymid;
      double t1,t2,t3,w1,w2,w3;
      double x1,x2,x3,y1,y2,y3;
      double udotn[4], gp[NDIM];
      double qleft[NUM_VAR];
      double u,v;
      double XL[2],XR[2],GP[6];
      double xL,xR,yL,yR;
      double L1,L2,Linf;

      L1 = 0.;
      L2 = 0.;
      Linf = 0.;

      if (RECON_PRIM)
	{
	  for ( i=1; i <= grid.nn; i++ )
	    {
	      ConvertConservedToPrimitive ( p.gamma , &(grid.Q[i*NUM_VAR]) );
	      ConvertConservedToPrimitive ( p.gamma , &(grid.nQ[i*NUM_VAR]) );  // for use in reconstruction.
	    }
	}

      for ( i=1; i <= grid.nbedges; i++ )
	{
	  // Grab the data from the structure.
	  real_node = grid.bedges[i*5+0];
	  ghost_node = grid.bedges[i*5+1] + grid.nn;
	  
	  b = grid.bedges[i*5+3];
	  bc = grid.bbc[b];

	  if ( ANNULUS )
	    {
	      if ( bc==0 || bc==2 )
		continue;
	    }

	  if ( CYLINDER )
	    {
	      if ( bc%2 != 1 )
		continue;
	    }
	  
	  // Get the normal vector information.
	  nx = grid.xn_bedges[i*3+0];
	  ny = grid.xn_bedges[i*3+1];
	  len = grid.xn_bedges[i*3+2];

	  if ( p.order > 2 )
	    {
	      // Get the segment.
	      s = grid.bedges[i*5+4];
	  
	      // Now get the nodes attached to the edge.
	      nodeL = grid.bs[b][s][0];
	      nodeR = grid.bs[b][s][1];

	      // Get the node coordinates.
	      xi = grid.x[real_node*NDIM+0];
	      yi = grid.x[real_node*NDIM+1];

	      // Get the edge midpoint.
	      xmid = 0.5*( grid.x[nodeL*NDIM+0] + grid.x[nodeR*NDIM+0] );
	      ymid = 0.5*( grid.x[nodeL*NDIM+1] + grid.x[nodeR*NDIM+1] );

	      if ( bc <= 10 ) // Linear boundary.
		{
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
		  w1 = (5./9.)*0.5;
		  
		  t2 = 0.;
		  w2 = (8./9.)*0.5;
		  
		  t3 = sqrt(15.0)/(-5.0);
		  w3 = w1;
		  
		  x1 = (xL+xR)*0.5 + (xR-xL)*0.5*t1;
		  y1 = (yL+yR)*0.5 + (yR-yL)*0.5*t1;
		  
		  x2 = (xL+xR)*0.5 + (xR-xL)*0.5*t2;
		  y2 = (yL+yR)*0.5 + (yR-yL)*0.5*t2;
		  
		  x3 = (xL+xR)*0.5 + (xR-xL)*0.5*t3;
		  y3 = (yL+yR)*0.5 + (yR-yL)*0.5*t3;

		  gp[0] = xi;  gp[1] = yi;
		  Reconstruct_Gauss_Boundary ( real_node, gp, qleft, &grid, p);

		  if ( RECON_PRIM )
		    {
		      u = qleft[1];
		      v = qleft[2];
		    }
		  else
		    {
		      u = qleft[1] / qleft[0];
		      v = qleft[2] / qleft[0];
		    }
		      
		  udotn[0] = u*nx + v*ny;

		  L1 += fabs( udotn[0] );
		  L2 += ( udotn[0] * udotn[0] );
		  Linf = MAX( fabs(udotn[0]) , Linf );

		  gp[0] = x1;  gp[1] = y1;
		  Reconstruct_Gauss_Boundary ( real_node, gp, qleft, &grid, p);

		  if ( RECON_PRIM )
		    {
		      u = qleft[1];
		      v = qleft[2];
		    }
		  else
		    {
		      u = qleft[1] / qleft[0];
		      v = qleft[2] / qleft[0];
		    }
		      
		  udotn[1] = u*nx + v*ny;

		  L1 += fabs( udotn[1] );
		  L2 += ( udotn[1] * udotn[1] );
		  Linf = MAX( fabs(udotn[1]) , Linf );

		  gp[0] = x2;  gp[1] = y2;
		  Reconstruct_Gauss_Boundary ( real_node, gp, qleft, &grid, p);

		  if ( RECON_PRIM )
		    {
		      u = qleft[1];
		      v = qleft[2];
		    }
		  else
		    {
		      u = qleft[1] / qleft[0];
		      v = qleft[2] / qleft[0];
		    }
		      
		  udotn[2] = u*nx + v*ny;

		  L1 += fabs( udotn[2] );
		  L2 += ( udotn[2] * udotn[2] );
		  Linf = MAX( fabs(udotn[2]) , Linf );

		  gp[0] = x3;  gp[1] = y3;
		  Reconstruct_Gauss_Boundary ( real_node, gp, qleft, &grid, p);

		  if ( RECON_PRIM )
		    {
		      u = qleft[1];
		      v = qleft[2];
		    }
		  else
		    {
		      u = qleft[1] / qleft[0];
		      v = qleft[2] / qleft[0];
		    }
		      
		  udotn[3] = u*nx + v*ny;

		  L1 += fabs( udotn[3] );
		  L2 += ( udotn[3] * udotn[3] );
		  Linf = MAX( fabs(udotn[3]) , Linf );
		}
	      else
		{
		  if ( real_node == nodeL )
		    {
		      XL[0] = grid.x[real_node*NDIM+0];
		      XL[1] = grid.x[real_node*NDIM+1];
		      XR[0] = grid.x[ghost_node*NDIM+0];
		      XR[1] = grid.x[ghost_node*NDIM+1];
		    }
		  else
		    {
		      XR[0] = grid.x[real_node*NDIM+0];
		      XR[1] = grid.x[real_node*NDIM+1];
		      XL[0] = grid.x[ghost_node*NDIM+0];
		      XL[1] = grid.x[ghost_node*NDIM+1];
		    }

		  // Get the Gauss points.
		  curved_boundary_gauss_points ( &grid, bc, XL, XR, GP );
		  
		  curved_boundary_arclength ( &grid, bc, XL, XR, &len );

		  // do the node itself first.
		  gp[0] = grid.x[real_node*NDIM+0];      gp[1] = grid.x[real_node*NDIM+1];
		  curved_boundary_normal_vector ( &grid, bc, gp, &nx, &ny );
		  Reconstruct_Gauss_Boundary ( real_node, gp, qleft, &grid, p);

		  if ( RECON_PRIM )
		    {
		      u = qleft[1];
		      v = qleft[2];
		    }
		  else
		    {
		      u = qleft[1] / qleft[0];
		      v = qleft[2] / qleft[0];
		    }
		      
		  udotn[0] = u*nx + v*ny;
		  
		  L1 += fabs( udotn[0] );
		  L2 += ( udotn[0] * udotn[0] );
		  Linf = MAX( fabs(udotn[0]) , Linf );

		  for ( j=0; j < 3; j++ )
		    {
		      // Get the normal vector.
		      curved_boundary_normal_vector ( &grid, bc, &(GP[j*NDIM]), &nx, &ny );

		      Reconstruct_Gauss_Boundary ( real_node, &(GP[j*NDIM]), qleft, &grid, p);

		      if ( RECON_PRIM )
			{
			  u = qleft[1];
			  v = qleft[2];
			}
		      else
			{
			  u = qleft[1] / qleft[0];
			  v = qleft[2] / qleft[0];
			}
		      
		      udotn[j+1] = u*nx + v*ny;

		      L1 += fabs( udotn[j+1] );
		      L2 += ( udotn[j+1] * udotn[j+1] );
		      Linf = MAX( fabs(udotn[j+1]) , Linf );
		    }
		}

	      // Print out the results for the boundary edge.
	      printf("Boundary Edge %d:\n",i);
	      printf("Node %d : U dot n = %.15E\n",real_node,udotn[0]);
	      printf("Gauss Point 1: U dot n = %.15E\n",udotn[1]);
	      printf("Gauss Point 2: U dot n = %.15E\n",udotn[2]);
	      printf("Gauss Point 3: U dot n = %.15E\n",udotn[3]);
	      printf("\n");
	    }
	  else
	    {
	      // Get the node coordinates.
	      xi = grid.x[real_node*NDIM+0];
	      yi = grid.x[real_node*NDIM+1];
	      gp[0] = xi;
	      gp[1] = yi;
	      
	      //Reconstruct_Gauss_Boundary ( real_node, gp, qleft, &grid, p);
	      qleft[0] = grid.nQ[real_node*NUM_VAR+0];
	      qleft[1] = grid.nQ[real_node*NUM_VAR+1];
	      qleft[2] = grid.nQ[real_node*NUM_VAR+2];
	      qleft[3] = grid.nQ[real_node*NUM_VAR+3];
	      	      
	      if ( RECON_PRIM )
		{
		  u = qleft[1];
		  v = qleft[2];
		}
	      else
		{
		  u = qleft[1] / qleft[0];
		  v = qleft[2] / qleft[0];
		}
	      
	      udotn[0] = u*nx + v*ny;

	      L1 += fabs( udotn[0] );
	      L2 += ( udotn[0] * udotn[0] );
	      Linf = MAX( fabs(udotn[0]) , Linf );

	      // Print out the results for the boundary edge.
	      printf("Boundary Edge %d:\n",i);
	      printf("Node %d : U dot n = %.15E\n",real_node,udotn[0]);
	      printf("\n");
	    }  
	}

      // norm results.
      printf("L1 norm = %.15E\n",L1);
      printf("L2 norm = %.15E\n",sqrt(L2));
      printf("Infinity norm = %.15E\n",Linf);
    }

  if ( tfp != NULL )
    fclose(tfp);
  
  freenull(rms);

  grid.~GRID();
  smem.~SYS_MEM();

  return 0;
}


// Here I am integrating the error of the exact solution minus the reconstructed solution over the control volumes
// with respect to the true geometry (so I'm using curved boundaries regardless of the bc's that were in the input file).
void compute_integral_error ( GRID *grid, PARAMS p )
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
  double dx[NDIM];
  double Qrecon[NUM_VAR];

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

  double *error = NULL;

  double rhoi,r,r2,Mi,U,Ui,Ri;
  double rhoE,uE,vE,PE,EE,TE;       // Exact values.
  double rhoR,uR,vR,PR,ER;          // Reconstructed values.

  double Pt, Ptinf;                 // Total pressure.

  double Pinf , Tinf;
  double T;

  double tA = 0.;
  double L1[4];
  double L2[4];
  double Linf[4];

  double pi = M_PI;

  // for vortex problems
  double gamma = p.gamma;                          // Ratio of coefficients of heat.
  double gm1 = gamma - 1.0;
  double oogm1 = 1.0 / gm1;
  double mach = p.mach;                            // Mach number.
  double x0;                                       // x location of vortex at time t = 0
  double y0;                                       // y location of vortex at time t = 0
  double beta;                                     // Vortex strength.
  double fac1,fac2,fac3;                           // Terms in the computation.
  double rho0;                                     // Core density.
  double u0,v0;                                    // Core velocity.
  double P0,T0,s0;                                 // Core pressure,temperature, and entropy.

  double phi,sigma;                                // Parameter to control vortex strength.
  double upert,vpert,Tpert;                        // Isentropic vortex perturbations.
  double Mref,rhoref,Pref,Tref;
  double rhoD,PD;                                  // Dimensional quantities.
  double FC;                                       // Fixed constant.
  double angle;                                    // Angle of attack in radians.

  int num_expected_ct = 0;                         // Number of curved triangles I expect to process.
  int num_processed_ct = 0;                        // Number of curved triangles I actually processed.

  int n0, n1;
  double r1;

  int s,nodeL,nodeR,gelem,iedge;
  double xi,yi,nx,ny,len,xmid,ymid;
  double t1,t2,t3,x1,x2,x3,y1,y2,y3,w1,w2,w3,xL,xR,yL,yR;
  double XL[2],XR[2],XM[2],GP[6];
  double *u_exact = NULL;
  double *v_exact = NULL;
  double w[3];

  // for the cylinder
  double v_mag, c_local, mach_local;

  w[0] = 8./9.;
  w[1] = 5./9.;
  w[2] = 5./9.;

  Mref = 340.294065;
  Tref = 288.15;
  Pref = 101325.0;
  rhoref = 1.225;

  FC = Pref / pow(rhoref,gamma);

  rhoi = 1.0;
  Mi = 2.0;
  Ui = 2.0;
  Ri = 2.0;

  angle = p.alpha*M_PI/180.;

  // Set total pressure at inifinity. This is intended for use in the circular cylinder problem.
  u0 = p.mach * cos( angle );
  v0 = p.mach * sin( angle );
  rho0 = 1.0;

  //Ptinf = 1.0/gamma + 0.5 * rho0 * ( u0*u0 + v0*v0 );

  //Ptinf = 101325.0 + 0.5 * 1.225 * ( 0.09 * 1.4 * 101325.0 / 1.225 );

  Ptinf = 1.0 / gamma * ( pow( ( 1 + 0.5*(gamma-1.)*mach*mach) , (gamma/(gamma-1.)) ) );

  Pinf = 1. / gamma ;
  Tinf = 1.;

  // Initialize variables for the vortex problem(s).
  x0 = 5. + ( grid->citer * p.mintime * p.mach );   // Put the vortex center at the correct spot.
  y0 = 0.;
  
  sigma = 4.0;
  phi = 1.0;

  beta = 5.0;
  fac1 = beta/(2.0*M_PI);
  fac2 = -0.5*gm1/gamma*fac1*fac1;

  // Check to make sure that this function was called with something to do analysis on!
  if ( (MMS+MMS_EXP+ANNULUS+CYLINDER+p.ic) != 1 )
    {
      printf("ERROR: In calling compute_integral_error(). Compile time flags are set wrong and may cause undesired behavior.\n");
      printf("  MMS = %d\n",MMS);
      printf("  MMS_EXP = %d\n",MMS_EXP);
      printf("  ANNULUS = %d\n",ANNULUS);
      printf("  CYLINDER = %d\n",CYLINDER);
      printf("  p.ic = %d\n",p.ic);
    }

  // If RECON_PRIM was set to one, then grid->Q has already been set back to primitives. This is good since for high order reconstruction
  // the point value and derivatives will be based on primitive variables.

  // This is very hacky and is for the annulus unstructured grid sequence only!
  if ( ANNULUS )
    {
      reset_flag = 0;

      for ( i=0; i <= grid->nb; i++ )
	{
	  printf("bbc[%d] = %d\n",i,grid->bbc[i]);
	}

      // old uns sequence and the tri sequence
      
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
	  
	  printf("RECALCULATING Control volume areas to account for curved boundaries while computing the error between exact and reconstructed solutions.\n");

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
	  
	  printf("RECALCULATING Control volume areas to account for curved boundaries while computing the error between exact and reconstructed solutions.\n");

	  cv_calc_area_Green ( grid );
	}
      */

    }

  // This is also very hacky and is for the circle tri grid sequence (1,2,3) only!
  if ( CYLINDER )
    {
      reset_flag = 0;

      if ( grid->bbc[1] != 11 )
	{
	  dummy_bcs[0] = grid->bbc[0];
	  dummy_bcs[1] = grid->bbc[1];
	  dummy_bcs[2] = grid->bbc[2];
	  dummy_bcs[3] = grid->bbc[3];
	  dummy_bcs[4] = grid->bbc[4];
	  dummy_bcs[5] = grid->bbc[5];
	  dummy_bcs[6] = grid->bbc[6];
	  dummy_bcs[7] = grid->bbc[7];
	  dummy_bcs[8] = grid->bbc[8];

	  grid->bbc[1] = 11;
	  grid->bbc[2] = 11;
	  grid->bbc[3] = 22;
	  grid->bbc[4] = 22;
	  grid->bbc[5] = 11;
	  grid->bbc[6] = 11;
	  grid->bbc[7] = 22;
	  grid->bbc[8] = 22;

	  reset_flag = 1;
	  
	  printf("RECALCULATING Control volume areas to account for curved boundaries while computing the error between exact and reconstructed solutions.\n");

	  cv_calc_area_Green ( grid );

	}

    }

  // Compute the number of curved triangles I expect to process. Each segment on a boundary that is tagged as curved will contribute two curved triangles.
  for ( i=1; i <= grid->nb; i++ )
    {
      if ( grid->bbc[i] > 10 )
	{
	  num_expected_ct += ( 2 * grid->nbs[i] );
	}
    }

  // Now let's go ahead and loop over all the curved boundaries of interest (circular cylinder and the annulus) and make the points
  // are actually lying on the boundary.

  for ( i=1; i <= grid->nb; i++ )
    {
      if ( grid->bbc[i] <= 10 )
	continue;                      // skip all the linear boundaries.

      // skip everything.
      continue;

      for ( j=1; j <= grid->nbs[i]; j++ )
	{
	  n0 = grid->bs[i][j][0];
	  n1 = grid->bs[i][j][1];

	  bc = grid->bbc[i];

	  r1 = sqrt( (grid->x[n0*NDIM+0])*(grid->x[n0*NDIM+0]) + (grid->x[n0*NDIM+1])*(grid->x[n0*NDIM+1]) );
	  r2 = sqrt( (grid->x[n1*NDIM+0])*(grid->x[n1*NDIM+0]) + (grid->x[n1*NDIM+1])*(grid->x[n1*NDIM+1]) );

	  if ( bc == 11 )
	    {
	      if ( fabs( r1 - 1.0 ) > 1.0E-10 )
		{
		  printf("ERROR: Node %d is tagged with bc %d: It is at x=%.15E  y=%.15E    r=%.15E\n",n0,bc,grid->x[n0*NDIM+0],grid->x[n0*NDIM+1],r1);
		}

	      if ( fabs( r2 - 1.0 ) > 1.0E-10 )
		{
		  printf("ERROR: Node %d is tagged with bc %d: It is at x=%.15E  y=%.15E    r=%.15E\n",n1,bc,grid->x[n1*NDIM+0],grid->x[n1*NDIM+1],r2);
		}
	    }

	  if ( bc == 14 || bc == 13 )
	    {
	      if ( fabs( r1 - 2.0 ) > 1.0E-10 )
		{
		  printf("ERROR: Node %d is tagged with bc %d: It is at x=%.15E  y=%.15E    r=%.15E\n",n0,bc,grid->x[n0*NDIM+0],grid->x[n0*NDIM+1],r1);
		}

	      if ( fabs( r2 - 2.0 ) > 1.0E-10 )
		{
		  printf("ERROR: Node %d is tagged with bc %d: It is at x=%.15E  y=%.15E    r=%.15E\n",n1,bc,grid->x[n1*NDIM+0],grid->x[n1*NDIM+1],r2);
		}
	    }

	  if ( bc == 16 || bc == 15 )
	    {
	      if ( fabs( r1 - 3.0 ) > 1.0E-10 )
		{
		  printf("ERROR: Node %d is tagged with bc %d: It is at x=%.15E  y=%.15E    r=%.15E\n",n0,bc,grid->x[n0*NDIM+0],grid->x[n0*NDIM+1],r1);
		}

	      if ( fabs( r2 - 3.0 ) > 1.0E-10 )
		{
		  printf("ERROR: Node %d is tagged with bc %d: It is at x=%.15E  y=%.15E    r=%.15E\n",n1,bc,grid->x[n1*NDIM+0],grid->x[n1*NDIM+1],r2);
		}
	    }
	}
    }

  error = (double*)malloc((grid->nn + 1)*NUM_VAR*sizeof(double));
  if ( error == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'error' in compute integral error (main.C).\n"); exit(1); }
  
  // Clean out memory because we are accumulating.
  for ( n=1; n <= grid->nn; n++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	error[n*NUM_VAR+j] = 0.;
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
		      printf("CRITICAL ERROR: Triangle 1 of element %d of type %d is too concave (compute integral error: main.C)!\n",j,e);
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

		      // Here we put the quantity to be integrated.
		      
		      // First reconstruct to the Gauss node.
		      dx[0] = XGP[0] - grid->x[node*NDIM+0];
		      dx[1] = XGP[1] - grid->x[node*NDIM+1];
		      
		      generic_reconstruct ( node, Qrecon, dx, grid, p );

		      // Now compute the exact solution at the Gauss node.

		      if ( ANNULUS )
			{
			  r2 = (XGP[0])*(XGP[0]) + (XGP[1])*(XGP[1]); 
			  r = sqrt(r2);
			  
			  rhoE = pow( ( 1.0 + ((p.gamma - 1.0)/2.0)*Mi*Mi*( 1.0 - (Ri*Ri)/r2)) , (1.0/(p.gamma-1.0)) );
			  U   = (Ui*Ri)/r;
			  uE  = (XGP[1]*U)/r;
			  vE  = (-1.0 * XGP[0]*U)/r;
			  PE  = ( pow( rhoE , p.gamma) )/(p.gamma);
			  EE  = PE / (gm1) + 0.5*rhoE*(uE*uE + vE*vE);
			}

		      if ( MMS )
			{
			  rhoE = 1.0 + 0.25*sin(pi*XGP[0])*sin(pi*XGP[1]);
			  uE = 0.25 + 0.25*sin(pi*XGP[0])*cos(2.*pi*XGP[1]);
			  vE = 0.25 + 0.25*cos(2.*pi*XGP[0])*sin(pi*XGP[1]);
			  PE = 1.0/p.gamma + 0.05*cos(2.*pi*XGP[0])*cos(2.*pi*XGP[1]);
			  EE = PE / (gm1) + 0.5*rhoE*(uE*uE + vE*vE);
			}

		      if ( MMS_EXP )
			{
			  rhoE = exp(XGP[0]-1.0) * exp(XGP[1] - 1.0);
			  uE   = exp(XGP[0]-1.0) * exp(XGP[1] - 1.0);
			  vE   = exp(XGP[0]-1.0) * exp(XGP[1] - 1.0);
			  EE   = 1.5 * exp( 3.0 * ( XGP[0] - 1.0) ) * exp( 3.0 * ( XGP[1] - 1.0 ) );
			  PE   = ( (gm1) * 0.5 ) * exp( 3.0 * ( XGP[0] - 1.0) ) * exp( 3.0 * ( XGP[1] - 1.0 ) );
			}

		      if ( p.ic == 1 )
			{
			  // Initialize the variables.
			  r2 = (XGP[0] - x0)*(XGP[0] - x0) + (XGP[1] - y0)*(XGP[1] - y0);
			  
			  fac3 = exp((1.0 - r2)/2.0);
			  
			  rho0 = 1.0;
			  u0 = p.mach*1.0;          // Flow in only the x direction at reference velocity of 1.
			  v0 = 0.;
			  P0 = 1.0/gamma;
			  T0 = P0/rho0;
			  s0 = P0/pow(rho0,gamma);
			  
			  TE = T0 + fac2*fac3*fac3;
			  uE = u0 - (fac1*fac3)*(XGP[1] - y0);
			  vE = v0 + (fac1*fac3)*(XGP[0] - x0);
			  PE = pow(pow(TE, gamma)/s0,1.0/gm1);
			  rhoE = PE/TE;
			  EE = PE / (gm1) + 0.5*rhoE*(uE*uE + vE*vE);
			}

		      if ( p.ic == 2 )
			{
			  // Initialize the variables.
			  r2 = (XGP[0] - x0)*(XGP[0] - x0) + (XGP[1] - y0)*(XGP[1] - y0);

			  upert = (-sigma)/(2.*pi) * ( XGP[1] - y0 ) * exp( phi*(1.-r2) );
			  vpert =  (sigma)/(2.*pi) * ( XGP[0] - x0 ) * exp( phi*(1.-r2) );
			  Tpert = ( -1. * sigma*sigma * gm1 )/( 16. * phi * gamma * pi * pi ) * exp( 2. * phi * (1.-r2) );

			  // Compute the new values for the node.
			  TE = 1.0 + Tpert;
			  rhoE = pow( ( TE ) , oogm1 );
			  rhoD = rhoE*rhoref;
			  uE = p.mach + upert;
			  vE = 0. + vpert;
			  PD = pow( rhoD , gamma ) * FC;
			  PE = PD / (gamma * Pref );
			  rhoE = rhoD;
			  EE = PE / gm1 + 0.5*rhoE*(uE*uE + vE*vE);
			}
		      
		      // Finally the error.

		      //if ( RECON_PRIM )
		      //{
		      //  QGP[0] = fabs( rhoE - Qrecon[0] );
		      //  QGP[1] = fabs( uE -   Qrecon[1] );
		      //  QGP[2] = fabs( vE -   Qrecon[2] );
		      //  QGP[3] = fabs( PE -   Qrecon[3] );
		      //}
		      //else
		      //{
		      //  QGP[0] = fabs( rhoE -    Qrecon[0] );
		      //  QGP[1] = fabs( rhoE*uE - Qrecon[1] );
		      //  QGP[2] = fabs( rhoE*vE - Qrecon[2] );
		      //  QGP[3] = fabs( EE -      Qrecon[3] );
		      //}

		      if ( RECON_PRIM )
			{
			  QGP[0] = ( rhoE - Qrecon[0] );
			  QGP[1] = ( uE -   Qrecon[1] );
			  QGP[2] = ( vE -   Qrecon[2] );
			  QGP[3] = ( PE -   Qrecon[3] );
			}
		      else
			{
			  QGP[0] = ( rhoE -    Qrecon[0] );
			  QGP[1] = ( rhoE*uE - Qrecon[1] );
			  QGP[2] = ( rhoE*vE - Qrecon[2] );
			  QGP[3] = ( EE -      Qrecon[3] );
			}

		      // for the circular cylinder, I can do what Gooch did.
		      if ( CYLINDER )
			{
			  if ( RECON_PRIM )
			    {
			      rhoR = Qrecon[0];
			      uR = Qrecon[1];
			      vR = Qrecon[2];
			      PR = Qrecon[3];
			    }
			  else
			    {
			      rhoR = Qrecon[0];
			      uR = Qrecon[1] / rhoR;
			      vR = Qrecon[2] / rhoR;
			      ER = Qrecon[3];
			      PR = (1.4 - 1.0)*( ER - 0.5*rhoR*(uR*uR + vR*vR) );
			    }

			  // dimensionalize the results.
			  //rhoR = rhoR * 1.225;
			  //uR = uR * Mref;
			  //vR = vR * Mref;
			  //PR = PR * p.gamma * 101325.0;

			  // compute total pressure.
			  //Pt = PR + 0.5*rhoR*(uR*uR + vR*vR);

			  v_mag = sqrt( uR*uR + vR*vR );
			  c_local = sqrt( gamma * PR / rhoR );
			  mach_local = v_mag / c_local;

			  Pt = PR * ( pow( ( 1.0 + 0.5*(gamma-1.)*mach_local*mach_local) , (gamma/(gamma-1.)) ) );
			  
			  // and the error.
			  QGP[0] = ( Pt / Ptinf - 1.0);
			}

		      if ( NACA )
			{
			  if ( RECON_PRIM )
			    {
			      rhoR = Qrecon[0];
			      uR = Qrecon[1];
			      vR = Qrecon[2];
			      PR = Qrecon[3];
			    }
			  else
			    {
			      rhoR = Qrecon[0];
			      uR = Qrecon[1] / rhoR;
			      vR = Qrecon[2] / rhoR;
			      ER = Qrecon[3];
			      PR = (1.4 - 1.0)*( ER - 0.5*rhoR*(uR*uR + vR*vR) );
			    }

			  // Get the local temperature.
			  T = (gamma * PR)/rhoR;

			  // compute the change in entropy.
			  QGP[0] = ( log( T / Tinf ) - log( PR / Pinf ) );
			}

		      // Accumulate.
		      for ( v=0; v < NUM_VAR; v++ )
			{
			  QTRI[v] += ( wg[i] * QGP[v] *jacobian * 0.5);
			} 
		      
		    }

		  // Add to the running total for the node.
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      error[node*NUM_VAR+v] += QTRI[v];
		    }

		  num_processed_ct++;
		  
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
		      
		      // Here we put the quantity to be integrated.
		      
		      // First reconstruct to the Gauss node.
		      dx[0] = XGP[0] - grid->x[node*NDIM+0];
		      dx[1] = XGP[1] - grid->x[node*NDIM+1];
		      
		      generic_reconstruct ( node, Qrecon, dx, grid, p );

		      // Now compute the exact solution at the Gauss node.

		      if ( ANNULUS )
			{
			  r2 = (XGP[0])*(XGP[0]) + (XGP[1])*(XGP[1]); 
			  r = sqrt(r2);
			  
			  rhoE = pow( ( 1.0 + ((p.gamma - 1.0)/2.0)*Mi*Mi*( 1.0 - (Ri*Ri)/r2)) , (1.0/(p.gamma-1.0)) );
			  U   = (Ui*Ri)/r;
			  uE  = (XGP[1]*U)/r;
			  vE  = (-1.0 * XGP[0]*U)/r;
			  PE  = ( pow( rhoE , p.gamma) )/(p.gamma);
			  EE  = PE / (gm1) + 0.5*rhoE*(uE*uE + vE*vE);
			}

		      if ( MMS )
			{
			  rhoE = 1.0 + 0.25*sin(pi*XGP[0])*sin(pi*XGP[1]);
			  uE = 0.25 + 0.25*sin(pi*XGP[0])*cos(2.*pi*XGP[1]);
			  vE = 0.25 + 0.25*cos(2.*pi*XGP[0])*sin(pi*XGP[1]);
			  PE = 1.0/p.gamma + 0.05*cos(2.*pi*XGP[0])*cos(2.*pi*XGP[1]);
			  EE = PE / (gm1) + 0.5*rhoE*(uE*uE + vE*vE);
			}

		      if ( MMS_EXP )
			{
			  rhoE = exp(XGP[0]-1.0) * exp(XGP[1] - 1.0);
			  uE   = exp(XGP[0]-1.0) * exp(XGP[1] - 1.0);
			  vE   = exp(XGP[0]-1.0) * exp(XGP[1] - 1.0);
			  EE   = 1.5 * exp( 3.0 * ( XGP[0] - 1.0) ) * exp( 3.0 * ( XGP[1] - 1.0 ) );
			  PE   = ( (gm1) * 0.5 ) * exp( 3.0 * ( XGP[0] - 1.0) ) * exp( 3.0 * ( XGP[1] - 1.0 ) );
			}

		      if ( p.ic == 1 )
			{
			  // Initialize the variables.
			  r2 = (XGP[0] - x0)*(XGP[0] - x0) + (XGP[1] - y0)*(XGP[1] - y0);
			  
			  fac3 = exp((1.0 - r2)/2.0);
			  
			  rho0 = 1.0;
			  u0 = p.mach*1.0;          // Flow in only the x direction at reference velocity of 1.
			  v0 = 0.;
			  P0 = 1.0/gamma;
			  T0 = P0/rho0;
			  s0 = P0/pow(rho0,gamma);
			  
			  TE = T0 + fac2*fac3*fac3;
			  uE = u0 - (fac1*fac3)*(XGP[1] - y0);
			  vE = v0 + (fac1*fac3)*(XGP[0] - x0);
			  PE = pow(pow(TE, gamma)/s0,1.0/gm1);
			  rhoE = PE/TE;
			  EE = PE / (gm1) + 0.5*rhoE*(uE*uE + vE*vE);
			}

		      if ( p.ic == 2 )
			{
			  // Initialize the variables.
			  r2 = (XGP[0] - x0)*(XGP[0] - x0) + (XGP[1] - y0)*(XGP[1] - y0);

			  upert = (-sigma)/(2.*pi) * ( XGP[1] - y0 ) * exp( phi*(1.-r2) );
			  vpert =  (sigma)/(2.*pi) * ( XGP[0] - x0 ) * exp( phi*(1.-r2) );
			  Tpert = ( -1. * sigma*sigma * gm1 )/( 16. * phi * gamma * pi * pi ) * exp( 2. * phi * (1.-r2) );

			  // Compute the new values for the node.
			  TE = 1.0 + Tpert;
			  rhoE = pow( ( TE ) , oogm1 );
			  rhoD = rhoE*rhoref;
			  uE = p.mach + upert;
			  vE = 0. + vpert;
			  PD = pow( rhoD , gamma ) * FC;
			  PE = PD / (gamma * Pref );
			  rhoE = rhoD;
			  EE = PE / gm1 + 0.5*rhoE*(uE*uE + vE*vE);
			}
		      
		      // Finally the error.

		      //if ( RECON_PRIM )
		      //{
		      //  QGP[0] = fabs( rhoE - Qrecon[0] );
		      //  QGP[1] = fabs( uE -   Qrecon[1] );
		      //  QGP[2] = fabs( vE -   Qrecon[2] );
		      //  QGP[3] = fabs( PE -   Qrecon[3] );
		      //}
		      //else
		      //{
		      //  QGP[0] = fabs( rhoE -    Qrecon[0] );
		      //  QGP[1] = fabs( rhoE*uE - Qrecon[1] );
		      //  QGP[2] = fabs( rhoE*vE - Qrecon[2] );
		      //  QGP[3] = fabs( EE -      Qrecon[3] );
		      //}

		      if ( RECON_PRIM )
			{
			  QGP[0] = ( rhoE - Qrecon[0] );
			  QGP[1] = ( uE -   Qrecon[1] );
			  QGP[2] = ( vE -   Qrecon[2] );
			  QGP[3] = ( PE -   Qrecon[3] );
			}
		      else
			{
			  QGP[0] = ( rhoE -    Qrecon[0] );
			  QGP[1] = ( rhoE*uE - Qrecon[1] );
			  QGP[2] = ( rhoE*vE - Qrecon[2] );
			  QGP[3] = ( EE -      Qrecon[3] );
			}

		      // for the circular cylinder, I can do what Gooch did.
		      if ( CYLINDER )
			{
			  if ( RECON_PRIM )
			    {
			      rhoR = Qrecon[0];
			      uR = Qrecon[1];
			      vR = Qrecon[2];
			      PR = Qrecon[3];
			    }
			  else
			    {
			      rhoR = Qrecon[0];
			      uR = Qrecon[1] / rhoR;
			      vR = Qrecon[2] / rhoR;
			      ER = Qrecon[3];
			      PR = (1.4 - 1.0)*( ER - 0.5*rhoR*(uR*uR + vR*vR) );
			    }

			  // dimensionalize the results.
			  //rhoR = rhoR * 1.225;
			  //uR = uR * Mref;
			  //vR = vR * Mref;
			  //PR = PR * p.gamma * 101325.0;

			  // compute total pressure.
			  //Pt = PR + 0.5*rhoR*(uR*uR + vR*vR);

			  v_mag = sqrt( uR*uR + vR*vR );
			  c_local = sqrt( gamma * PR / rhoR );
			  mach_local = v_mag / c_local;

			  Pt = PR * ( pow( ( 1.0 + 0.5*(gamma-1.)*mach_local*mach_local) , (gamma/(gamma-1.)) ) );
			  
			  // and the error.
			  QGP[0] = ( Pt / Ptinf - 1.0);
			}

		      if ( NACA )
			{
			  if ( RECON_PRIM )
			    {
			      rhoR = Qrecon[0];
			      uR = Qrecon[1];
			      vR = Qrecon[2];
			      PR = Qrecon[3];
			    }
			  else
			    {
			      rhoR = Qrecon[0];
			      uR = Qrecon[1] / rhoR;
			      vR = Qrecon[2] / rhoR;
			      ER = Qrecon[3];
			      PR = (1.4 - 1.0)*( ER - 0.5*rhoR*(uR*uR + vR*vR) );
			    }
			  
			  // Get the local temperature.
			  T = (gamma * PR)/rhoR;
			  
			  // compute the change in entropy.
			  QGP[0] = ( log( T / Tinf ) - log( PR / Pinf ) );
			}
		  
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
		      error[node*NUM_VAR+v] += QTRI[v];
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
		      printf("CRITICAL ERROR: Triangle 1 of element %d of type %d is too concave (compute integral error: main.C)!\n",j,e);
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

		      // Here we put the quantity to be integrated.
		      
		      // First reconstruct to the Gauss node.
		      dx[0] = XGP[0] - grid->x[node*NDIM+0];
		      dx[1] = XGP[1] - grid->x[node*NDIM+1];
		      
		      generic_reconstruct ( node, Qrecon, dx, grid, p );

		      // Now compute the exact solution at the Gauss node.

		      if ( ANNULUS )
			{
			  r2 = (XGP[0])*(XGP[0]) + (XGP[1])*(XGP[1]); 
			  r = sqrt(r2);
			  
			  rhoE = pow( ( 1.0 + ((p.gamma - 1.0)/2.0)*Mi*Mi*( 1.0 - (Ri*Ri)/r2)) , (1.0/(p.gamma-1.0)) );
			  U   = (Ui*Ri)/r;
			  uE  = (XGP[1]*U)/r;
			  vE  = (-1.0 * XGP[0]*U)/r;
			  PE  = ( pow( rhoE , p.gamma) )/(p.gamma);
			  EE  = PE / (gm1) + 0.5*rhoE*(uE*uE + vE*vE);
			}

		      if ( MMS )
			{
			  rhoE = 1.0 + 0.25*sin(pi*XGP[0])*sin(pi*XGP[1]);
			  uE = 0.25 + 0.25*sin(pi*XGP[0])*cos(2.*pi*XGP[1]);
			  vE = 0.25 + 0.25*cos(2.*pi*XGP[0])*sin(pi*XGP[1]);
			  PE = 1.0/p.gamma + 0.05*cos(2.*pi*XGP[0])*cos(2.*pi*XGP[1]);
			  EE = PE / (gm1) + 0.5*rhoE*(uE*uE + vE*vE);
			}

		      if ( MMS_EXP )
			{
			  rhoE = exp(XGP[0]-1.0) * exp(XGP[1] - 1.0);
			  uE   = exp(XGP[0]-1.0) * exp(XGP[1] - 1.0);
			  vE   = exp(XGP[0]-1.0) * exp(XGP[1] - 1.0);
			  EE   = 1.5 * exp( 3.0 * ( XGP[0] - 1.0) ) * exp( 3.0 * ( XGP[1] - 1.0 ) );
			  PE   = ( (gm1) * 0.5 ) * exp( 3.0 * ( XGP[0] - 1.0) ) * exp( 3.0 * ( XGP[1] - 1.0 ) );
			}

		      if ( p.ic == 1 )
			{
			  // Initialize the variables.
			  r2 = (XGP[0] - x0)*(XGP[0] - x0) + (XGP[1] - y0)*(XGP[1] - y0);
			  
			  fac3 = exp((1.0 - r2)/2.0);
			  
			  rho0 = 1.0;
			  u0 = p.mach*1.0;          // Flow in only the x direction at reference velocity of 1.
			  v0 = 0.;
			  P0 = 1.0/gamma;
			  T0 = P0/rho0;
			  s0 = P0/pow(rho0,gamma);
			  
			  TE = T0 + fac2*fac3*fac3;
			  uE = u0 - (fac1*fac3)*(XGP[1] - y0);
			  vE = v0 + (fac1*fac3)*(XGP[0] - x0);
			  PE = pow(pow(TE, gamma)/s0,1.0/gm1);
			  rhoE = PE/TE;
			  EE = PE / (gm1) + 0.5*rhoE*(uE*uE + vE*vE);
			}

		      if ( p.ic == 2 )
			{
			  // Initialize the variables.
			  r2 = (XGP[0] - x0)*(XGP[0] - x0) + (XGP[1] - y0)*(XGP[1] - y0);

			  upert = (-sigma)/(2.*pi) * ( XGP[1] - y0 ) * exp( phi*(1.-r2) );
			  vpert =  (sigma)/(2.*pi) * ( XGP[0] - x0 ) * exp( phi*(1.-r2) );
			  Tpert = ( -1. * sigma*sigma * gm1 )/( 16. * phi * gamma * pi * pi ) * exp( 2. * phi * (1.-r2) );
			  
			  // Compute the new values for the node.
			  TE = 1.0 + Tpert;
			  rhoE = pow( ( TE ) , oogm1 );
			  rhoD = rhoE*rhoref;
			  uE = p.mach + upert;
			  vE = 0. + vpert;
			  PD = pow( rhoD , gamma ) * FC;
			  PE = PD / (gamma * Pref );
			  rhoE = rhoD;
			  EE = PE / gm1 + 0.5*rhoE*(uE*uE + vE*vE);
			}

		      // Finally the error.
		      //if ( RECON_PRIM )
		      //{
		      //  QGP[0] = fabs( rhoE - Qrecon[0] );
		      //  QGP[1] = fabs( uE -   Qrecon[1] );
		      //  QGP[2] = fabs( vE -   Qrecon[2] );
		      //  QGP[3] = fabs( PE -   Qrecon[3] );
		      //}
		      //else
		      //{
		      //  QGP[0] = fabs( rhoE -    Qrecon[0] );
		      //  QGP[1] = fabs( rhoE*uE - Qrecon[1] );
		      //  QGP[2] = fabs( rhoE*vE - Qrecon[2] );
		      //  QGP[3] = fabs( EE -      Qrecon[3] );
		      //}

		      if ( RECON_PRIM )
			{
			  QGP[0] = ( rhoE - Qrecon[0] );
			  QGP[1] = ( uE -   Qrecon[1] );
			  QGP[2] = ( vE -   Qrecon[2] );
			  QGP[3] = ( PE -   Qrecon[3] );
			}
		      else
			{
			  QGP[0] = ( rhoE -    Qrecon[0] );
			  QGP[1] = ( rhoE*uE - Qrecon[1] );
			  QGP[2] = ( rhoE*vE - Qrecon[2] );
			  QGP[3] = ( EE -      Qrecon[3] );
			}

		      // for the circular cylinder, I can do what Gooch did.
		      if ( CYLINDER )
			{
			  if ( RECON_PRIM )
			    {
			      rhoR = Qrecon[0];
			      uR = Qrecon[1];
			      vR = Qrecon[2];
			      PR = Qrecon[3];
			    }
			  else
			    {
			      rhoR = Qrecon[0];
			      uR = Qrecon[1] / rhoR;
			      vR = Qrecon[2] / rhoR;
			      ER = Qrecon[3];
			      PR = (1.4 - 1.0)*( ER - 0.5*rhoR*(uR*uR + vR*vR) );
			    }

			  // dimensionalize the results.
			  //rhoR = rhoR * 1.225;
			  //uR = uR * Mref;
			  //vR = vR * Mref;
			  //PR = PR * p.gamma * 101325.0;

			  // compute total pressure.
			  //Pt = PR + 0.5*rhoR*(uR*uR + vR*vR);

			  v_mag = sqrt( uR*uR + vR*vR );
			  c_local = sqrt( gamma * PR / rhoR );
			  mach_local = v_mag / c_local;

			  Pt = PR * ( pow( ( 1.0 + 0.5*(gamma-1.)*mach_local*mach_local) , (gamma/(gamma-1.)) ) );
			  
			  // and the error.
			  QGP[0] = ( Pt / Ptinf - 1.0);
			}

		      if ( NACA )
			{
			  if ( RECON_PRIM )
			    {
			      rhoR = Qrecon[0];
			      uR = Qrecon[1];
			      vR = Qrecon[2];
			      PR = Qrecon[3];
			    }
			  else
			    {
			      rhoR = Qrecon[0];
			      uR = Qrecon[1] / rhoR;
			      vR = Qrecon[2] / rhoR;
			      ER = Qrecon[3];
			      PR = (1.4 - 1.0)*( ER - 0.5*rhoR*(uR*uR + vR*vR) );
			    }

			  // Get the local temperature.
			  T = (gamma * PR)/rhoR;

			  // compute the change in entropy.
			  QGP[0] = ( log( T / Tinf ) - log( PR / Pinf ) );
			}
		      
		      // Accumulate.
		      for ( v=0; v < NUM_VAR; v++ )
			{
			  QTRI[v] += ( wg[i] * QGP[v] *jacobian * 0.5);
			} 
		      
		    }

		  // Add to the running total for the node.
		  for ( v=0; v < NUM_VAR; v++ )
		    {
		      error[node*NUM_VAR+v] += QTRI[v];
		    }

		  num_processed_ct++;
		  
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
			  //gp[v][ind] = (bary_coord[v][0] * grid->x[node*NDIM+ind] +
			  //	bary_coord[v][1] * xmidR[ind] +
			  //	bary_coord[v][2] * xc[ind]);

			  gp[v][ind] = (bary_coord[v][0] * grid->x[node*NDIM+ind] +
					bary_coord[v][1] * xc[ind] +
					bary_coord[v][2] * xmidR[ind]);
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
		  
		      // Here we put the quantity to be integrated.
		      
		      // First reconstruct to the Gauss node.
		      dx[0] = XGP[0] - grid->x[node*NDIM+0];
		      dx[1] = XGP[1] - grid->x[node*NDIM+1];
		      
		      generic_reconstruct ( node, Qrecon, dx, grid, p );

		      // Now compute the exact solution at the Gauss node.

		      if ( ANNULUS )
			{
			  r2 = (XGP[0])*(XGP[0]) + (XGP[1])*(XGP[1]); 
			  r = sqrt(r2);
			  
			  rhoE = pow( ( 1.0 + ((p.gamma - 1.0)/2.0)*Mi*Mi*( 1.0 - (Ri*Ri)/r2)) , (1.0/(p.gamma-1.0)) );
			  U   = (Ui*Ri)/r;
			  uE  = (XGP[1]*U)/r;
			  vE  = (-1.0 * XGP[0]*U)/r;
			  PE  = ( pow( rhoE , p.gamma) )/(p.gamma);
			  EE  = PE / (gm1) + 0.5*rhoE*(uE*uE + vE*vE);
			}

		      if ( MMS )
			{
			  rhoE = 1.0 + 0.25*sin(pi*XGP[0])*sin(pi*XGP[1]);
			  uE = 0.25 + 0.25*sin(pi*XGP[0])*cos(2.*pi*XGP[1]);
			  vE = 0.25 + 0.25*cos(2.*pi*XGP[0])*sin(pi*XGP[1]);
			  PE = 1.0/p.gamma + 0.05*cos(2.*pi*XGP[0])*cos(2.*pi*XGP[1]);
			  EE = PE / (gm1) + 0.5*rhoE*(uE*uE + vE*vE);
			}

		      if ( MMS_EXP )
			{
			  rhoE = exp(XGP[0]-1.0) * exp(XGP[1] - 1.0);
			  uE   = exp(XGP[0]-1.0) * exp(XGP[1] - 1.0);
			  vE   = exp(XGP[0]-1.0) * exp(XGP[1] - 1.0);
			  EE   = 1.5 * exp( 3.0 * ( XGP[0] - 1.0) ) * exp( 3.0 * ( XGP[1] - 1.0 ) );
			  PE   = ( (gm1) * 0.5 ) * exp( 3.0 * ( XGP[0] - 1.0) ) * exp( 3.0 * ( XGP[1] - 1.0 ) );
			}

		      if ( p.ic == 1 )
			{
			  // Initialize the variables.
			  r2 = (XGP[0] - x0)*(XGP[0] - x0) + (XGP[1] - y0)*(XGP[1] - y0);
			  
			  fac3 = exp((1.0 - r2)/2.0);
			  
			  rho0 = 1.0;
			  u0 = p.mach*1.0;          // Flow in only the x direction at reference velocity of 1.
			  v0 = 0.;
			  P0 = 1.0/gamma;
			  T0 = P0/rho0;
			  s0 = P0/pow(rho0,gamma);
			  
			  TE = T0 + fac2*fac3*fac3;
			  uE = u0 - (fac1*fac3)*(XGP[1] - y0);
			  vE = v0 + (fac1*fac3)*(XGP[0] - x0);
			  PE = pow(pow(TE, gamma)/s0,1.0/gm1);
			  rhoE = PE/TE;
			  EE = PE / (gm1) + 0.5*rhoE*(uE*uE + vE*vE);
			}

		      if ( p.ic == 2 )
			{
			  // Initialize the variables.
			  r2 = (XGP[0] - x0)*(XGP[0] - x0) + (XGP[1] - y0)*(XGP[1] - y0);

			  upert = (-sigma)/(2.*pi) * ( XGP[1] - y0 ) * exp( phi*(1.-r2) );
			  vpert =  (sigma)/(2.*pi) * ( XGP[0] - x0 ) * exp( phi*(1.-r2) );
			  Tpert = ( -1. * sigma*sigma * gm1 )/( 16. * phi * gamma * pi * pi ) * exp( 2. * phi * (1.-r2) );

			  // Compute the new values for the node.
			  TE = 1.0 + Tpert;
			  rhoE = pow( ( TE ) , oogm1 );
			  rhoD = rhoE*rhoref;
			  uE = p.mach + upert;
			  vE = 0. + vpert;
			  PD = pow( rhoD , gamma ) * FC;
			  PE = PD / (gamma * Pref );
			  rhoE = rhoD;
			  EE = PE / gm1 + 0.5*rhoE*(uE*uE + vE*vE);
			}

		      // Finally the error.
		      //if ( RECON_PRIM )
		      //{
		      //  QGP[0] = fabs( rhoE - Qrecon[0] );
		      //  QGP[1] = fabs( uE -   Qrecon[1] );
		      //  QGP[2] = fabs( vE -   Qrecon[2] );
		      //  QGP[3] = fabs( PE -   Qrecon[3] );
		      //}
		      //else
		      //{
		      //  QGP[0] = fabs( rhoE -    Qrecon[0] );
		      //  QGP[1] = fabs( rhoE*uE - Qrecon[1] );
		      //  QGP[2] = fabs( rhoE*vE - Qrecon[2] );
		      //  QGP[3] = fabs( EE -      Qrecon[3] );
		      //}

		      if ( RECON_PRIM )
			{
			  QGP[0] = ( rhoE - Qrecon[0] );
			  QGP[1] = ( uE -   Qrecon[1] );
			  QGP[2] = ( vE -   Qrecon[2] );
			  QGP[3] = ( PE -   Qrecon[3] );
			}
		      else
			{
			  QGP[0] = ( rhoE -    Qrecon[0] );
			  QGP[1] = ( rhoE*uE - Qrecon[1] );
			  QGP[2] = ( rhoE*vE - Qrecon[2] );
			  QGP[3] = ( EE -      Qrecon[3] );
			}

		      // for the circular cylinder, I can do what Gooch did.
		      if ( CYLINDER )
			{
			  if ( RECON_PRIM )
			    {
			      rhoR = Qrecon[0];
			      uR = Qrecon[1];
			      vR = Qrecon[2];
			      PR = Qrecon[3];
			    }
			  else
			    {
			      rhoR = Qrecon[0];
			      uR = Qrecon[1] / rhoR;
			      vR = Qrecon[2] / rhoR;
			      ER = Qrecon[3];
			      PR = (1.4 - 1.0)*( ER - 0.5*rhoR*(uR*uR + vR*vR) );
			    }

			  // dimensionalize the results.
			  //rhoR = rhoR * 1.225;
			  //uR = uR * Mref;
			  //vR = vR * Mref;
			  //PR = PR * p.gamma * 101325.0;

			  // compute total pressure.
			  //Pt = PR + 0.5*rhoR*(uR*uR + vR*vR);

			  v_mag = sqrt( uR*uR + vR*vR );
			  c_local = sqrt( gamma * PR / rhoR );
			  mach_local = v_mag / c_local;

			  Pt = PR * ( pow( ( 1.0 + 0.5*(gamma-1.)*mach_local*mach_local) , (gamma/(gamma-1.)) ) );
			  
			  // and the error.
			  QGP[0] = ( Pt / Ptinf - 1.0);
			}

		      if ( NACA )
			{
			  if ( RECON_PRIM )
			    {
			      rhoR = Qrecon[0];
			      uR = Qrecon[1];
			      vR = Qrecon[2];
			      PR = Qrecon[3];
			    }
			  else
			    {
			      rhoR = Qrecon[0];
			      uR = Qrecon[1] / rhoR;
			      vR = Qrecon[2] / rhoR;
			      ER = Qrecon[3];
			      PR = (1.4 - 1.0)*( ER - 0.5*rhoR*(uR*uR + vR*vR) );
			    }

			  // Get the local temperature.
			  T = (gamma * PR)/rhoR;

			  // compute the change in entropy.
			  QGP[0] = ( log( T / Tinf ) - log( PR / Pinf ) );
			}

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
		      error[node*NUM_VAR+v] += QTRI[v];
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
	  error[i*NUM_VAR+j] /= grid->cv_area[i];
	}
    }

  // Now lets get the magnitude of the error.
  for ( i=1; i<= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  error[i*NUM_VAR+j] = fabs( error[i*NUM_VAR+j] );
	}
    }

  for ( i=1; i<= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{
	  if ( !isfinite( error[i*NUM_VAR+j] ) )
	    printf(" error for node %d variable %d is not finite: %.15e\n", i,j,error[i*NUM_VAR+j] );
	}
    }

  for ( i=1; i<= grid->nn; i++ )
    {
      if ( !isfinite( grid->cv_area[i] ) )
	printf(" CV area for node %d is not finite: %.15e\n", i,grid->cv_area[i] );
    }

  if ( CYLINDER )
    {
      // Copy the error to grid->Q so I can write it out.
      
      // move grid->Q[0] to error[1].
      for (i=1; i<= grid->nn; i++ )
	{
	  error[i*NUM_VAR+1] = grid->Q[i*NUM_VAR+1];
	  grid->Q[i*NUM_VAR+0] = error[i*NUM_VAR+0];
	}

      // write it out.
      write_tecplot_solutionC ( "cylinder_error_in_total_pressure.dat", grid );

      // copy the stuff back to grid->Q.
      for (i=1; i<= grid->nn; i++ )
	{
	  grid->Q[i*NUM_VAR+0] = error[i*NUM_VAR+1];
	}
    }

  if ( NACA )
    {
      // Copy the error to grid->Q so I can write it out.
      
      // move grid->Q[0] to error[1].
      for (i=1; i<= grid->nn; i++ )
	{
	  error[i*NUM_VAR+1] = grid->Q[i*NUM_VAR+1];
	  grid->Q[i*NUM_VAR+0] = error[i*NUM_VAR+0];
	}
      
      // write it out.
      write_tecplot_solutionC ( "naca_error_in_entropy.dat", grid );
      
      // copy the stuff back to grid->Q.
      for (i=1; i<= grid->nn; i++ )
	{
	  grid->Q[i*NUM_VAR+0] = error[i*NUM_VAR+1];
	}
    }

  // Now lets compute the error analysis.

  for ( i=0; i < NUM_VAR; i++ )
    {
      L1[i] = 0.;
      L2[i] = 0.;
      Linf[i] = 0.;
    }
  
  tA = 0.0;
  for ( i=1; i <= grid->nn; i++ )
    {
      tA += grid->cv_area[i];
    }

  printf("tA (total area of the mesh) is %.15E\n",tA);

  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NUM_VAR; j++ )
	{	  
	  L1[j] += error[i*NUM_VAR+j] * grid->cv_area[i];
	  L2[j] += error[i*NUM_VAR+j] * error[i*NUM_VAR+j] * grid->cv_area[i];
	  Linf[j] = MAX( error[i*NUM_VAR+j] , Linf[j] );
	}
    }
  
  for ( i=0; i < NUM_VAR; i++ )
    {
      L1[i] = L1[i] / tA;
      L2[i] = L2[i] / tA;
      L2[i] = sqrt( L2[i] );
    }

  printf("\n\n");
  
  if ( CYLINDER ) 
    {
      // total pressure should be constant since shock-free. AIAA 2007-4194
      printf("Error for cylinder is ( Pt/Pt_inf - 1 )\n");
      printf("  L1 norm = %.15e\n",L1[0]);
      printf("  L2 norm = %.15e\n",L2[0]);
      printf("  L-infinity norm = %.15e\n\n",Linf[0]);

      printf("\n");
      printf("%.15e\n",L1[0]);
      printf("%.15e\n",L2[0]);
      printf("%.15e\n",Linf[0]);
    }
  else if ( NACA )
    {
      // entropy production.
      printf("Error for NACA 0012 is ( s_i - s_inf = ln(T_i/T_inf) - ln(P/P_inf) )\n");
      printf("  L1 norm = %.15e\n",L1[0]);
      printf("  L2 norm = %.15e\n",L2[0]);
      printf("  L-infinity norm = %.15e\n\n",Linf[0]);

      printf("\n");
      printf("%.15e\n",L1[0]);
      printf("%.15e\n",L2[0]);
      printf("%.15e\n",Linf[0]);
    }
  else
    {
      printf("Area integrated  error analysis. Calculated as ( |Error^p| *dA / tA )^1/p\n");
      printf("  where Error = 1/Area_cv INT_cv (Qexact-Qreconstructed) dA\n\n");
  
      printf("Density:\n");
      printf("  L1 norm = %.15e\n",L1[0]);
      printf("  L2 norm = %.15e\n",L2[0]);
      printf("  L-infinity norm = %.15e\n\n",Linf[0]);
      
      if ( RECON_PRIM )
	printf("X velocity:\n");
      else
	printf("X momentum:\n");
      printf("  L1 norm = %.15e\n",L1[1]);
      printf("  L2 norm = %.15e\n",L2[1]);
      printf("  L-infinity norm = %.15e\n\n",Linf[1]);
      
      if ( RECON_PRIM )
	printf("Y velocity:\n");
      else
	printf("Y momentum:\n");
      printf("  L1 norm = %.15e\n",L1[2]);
      printf("  L2 norm = %.15e\n",L2[2]);
      printf("  L-infinity norm = %.15e\n\n",Linf[2]);

      if ( RECON_PRIM)
	printf("Pressure:\n");
      else
	printf("Total Energy:\n");
      printf("  L1 norm = %.15e\n",L1[3]);
      printf("  L2 norm = %.15e\n",L2[3]);
      printf("  L-infinity norm = %.15e\n\n",Linf[3]);
      
      printf("\n");
      printf("%.15e\n",L1[0]);
      printf("%.15e\n",L2[0]);
      printf("%.15e\n",Linf[0]);
      printf("%.15e\n",L1[1]);
      printf("%.15e\n",L2[1]);
      printf("%.15e\n",Linf[1]);
      printf("%.15e\n",L1[2]);
      printf("%.15e\n",L2[2]);
      printf("%.15e\n",Linf[2]);
      printf("%.15e\n",L1[3]);
      printf("%.15e\n",L2[3]);
      printf("%.15e\n",Linf[3]);
    }

  printf("Expected to calculate %d curved triangles.\n",num_expected_ct);
  printf("Processed %d curved triangles.\n",num_processed_ct);
  
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

  if ( CYLINDER && reset_flag )
    {
      grid->bbc[1] = dummy_bcs[1];
      grid->bbc[2] = dummy_bcs[2];
      grid->bbc[3] = dummy_bcs[3];
      grid->bbc[4] = dummy_bcs[4];
      grid->bbc[5] = dummy_bcs[5];
      grid->bbc[6] = dummy_bcs[6];
      grid->bbc[7] = dummy_bcs[7];
      grid->bbc[8] = dummy_bcs[8];
      
      printf("RECALCULATING Control volume areas with the original boundary conditions (linear).\n");

      cv_calc_area_Green ( grid );

    }

  // some new code to check the velocities out of the annulus test case.
  // compute the area-averaged values of u and v using Green's Theorem since u,v have valid closed form solutions on the region.

  if ( (ANNULUS && RECON_PRIM) && 0 )
    {

      u_exact = (double*)malloc((grid->nn + 1)*sizeof(double));
      if ( u_exact == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'u_exact' in compute integral error (main.C).\n"); exit(1); }
      
      v_exact = (double*)malloc((grid->nn + 1)*sizeof(double));
      if ( v_exact == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'v_exact' in compute integral error (main.C).\n"); exit(1); }

      for ( i=0; i <= grid->nn; i++ )
	{
	  u_exact[i] = 0.;
	  v_exact[i] = 0.;
	}

      for ( i=1; i <= grid->nsubedges; i++ )
	{
	  // Find the global edge associated with the subedge.
	  iedge = grid->subedges[i*2];
	  
	  // Find the global element.
	  gelem = grid->subedges[i*2+1];
	  
	  // Now get the nodes attached to the edge.
	  nodeL = grid->edges[2*iedge+0];
	  nodeR = grid->edges[2*iedge+1];

	  // Get the node coordinates.
	  xi = grid->x[nodeL*NDIM+0];
	  yi = grid->x[nodeL*NDIM+1];

	  // Get the element centroid.
	  xc[0] = grid->el_cent[gelem*2];
	  xc[1] = grid->el_cent[gelem*2+1];

	  nx = grid->xn_subedges[i*3+0];
	  ny = grid->xn_subedges[i*3+1];
	  len = grid->xn_subedges[i*3+2];
	  
	  xmid = grid->xm_subedges[i*NDIM+0];
	  ymid = grid->xm_subedges[i*NDIM+1];
	  
	  t1 = sqrt(15.0)/5.0;
	  w1 = 5./9.;

	  t2 = 0.;
	  w2 = 8./9.;

	  t3 = sqrt(15.0)/(-5.0);
	  w3 = w1;

	  x1 = (xmid+xc[0])*0.5 + (xc[0]-xmid)*0.5*t1;
	  y1 = (ymid+xc[1])*0.5 + (xc[1]-ymid)*0.5*t1;
      
	  x2 = (xmid+xc[0])*0.5 + (xc[0]-xmid)*0.5*t2;
	  y2 = (ymid+xc[1])*0.5 + (xc[1]-ymid)*0.5*t2;
      
	  x3 = (xmid+xc[0])*0.5 + (xc[0]-xmid)*0.5*t3;
	  y3 = (ymid+xc[1])*0.5 + (xc[1]-ymid)*0.5*t3;

	  uE = w1 * ( 4. * atan( x1/y1 ) ) +
	       w2 * ( 4. * atan( x2/y2 ) ) +
	       w3 * ( 4. * atan( x3/y3 ) );

	  vE = w1 * ( -2. * log( x1*x1 + y1*y1 ) ) +
	       w2 * ( -2. * log( x2*x2 + y2*y2 ) ) +
	       w3 * ( -2. * log( x3*x3 + y3*y3 ) );

	  u_exact[nodeL] += ( uE * nx * len/2. );
	  v_exact[nodeL] += ( vE * nx * len/2. );

	  u_exact[nodeR] -= ( uE * nx * len/2. );
	  v_exact[nodeR] -= ( vE * nx * len/2. );
	}

      // Now I need to loop over and close off the boundaries.
  
      for ( i=1; i <= grid->nbedges; i++ )
	{
	  // Get the node.
	  node = grid->bedges[i*5+0];

	  // Get the boundary.
	  b = grid->bedges[i*5+3];
	  
	  // Get the segment.
	  s = grid->bedges[i*5+4];

	  bc = grid->bbc[b];

	  // Now get the nodes attached to the edge.
	  nodeL = grid->bs[b][s][0];
	  nodeR = grid->bs[b][s][1];

	  if ( bc <= 10 )
	    {
	      // Get the node coordinates.
	      xi = grid->x[node*NDIM+0];
	      yi = grid->x[node*NDIM+1];

	      // Get the edge midpoint.
	      xmid = 0.5*( grid->x[nodeL*NDIM+0] + grid->x[nodeR*NDIM+0] );
	      ymid = 0.5*( grid->x[nodeL*NDIM+1] + grid->x[nodeR*NDIM+1] );

	      // Get the length.
	      len = grid->xn_bedges[i*3+2];

	      // Retrieve the x component of the dual edge normal vector.
	      nx = grid->xn_bedges[i*3+0];
	      
	      if ( node == nodeL )
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
	      w1 = 5./9.;

	      t2 = 0.;
	      w2 = 8./9.;

	      t3 = sqrt(15.0)/(-5.0);
	      w3 = w1;
	      
	      x1 = (xL+xR)*0.5 + (xR-xL)*0.5*t1;
	      y1 = (yL+yR)*0.5 + (yR-yL)*0.5*t1;
	      
	      x2 = (xL+xR)*0.5 + (xR-xL)*0.5*t2;
	      y2 = (yL+yR)*0.5 + (yR-yL)*0.5*t2;

	      x3 = (xL+xR)*0.5 + (xR-xL)*0.5*t3;
	      y3 = (yL+yR)*0.5 + (yR-yL)*0.5*t3; 

	      uE = w1 * ( 4. * atan( x1/y1 ) ) +
		   w2 * ( 4. * atan( x2/y2 ) ) +
	           w3 * ( 4. * atan( x3/y3 ) );

	      vE = w1 * ( -2. * log( x1*x1 + y1*y1 ) ) +
	           w2 * ( -2. * log( x2*x2 + y2*y2 ) ) +
	           w3 * ( -2. * log( x3*x3 + y3*y3 ) );

	      u_exact[node] += ( uE * nx * len/2. );
	      v_exact[node] += ( vE * nx * len/2. );
	    }

	  else        // Curved boundary.
	    {
	      XL[0] = grid->x[nodeL*NDIM+0];
	      XL[1] = grid->x[nodeL*NDIM+1];
	      
	      XR[0] = grid->x[nodeR*NDIM+0];
	      XR[1] = grid->x[nodeR*NDIM+1];

	      // Get the midpoint.
	      curved_boundary_midpoint ( grid, bc, XL, XR, XM );
	      
	      if ( node == nodeL )
		{
		  XL[0] = grid->x[node*NDIM+0];
		  XL[1] = grid->x[node*NDIM+1];
		  XR[0] = XM[0];
		  XR[1] = XM[1];
		}
	      else
		{
		  XL[0] = XM[0];
		  XL[1] = XM[1];
		  XR[0] = grid->x[node*NDIM+0];
		  XR[1] = grid->x[node*NDIM+1];
		}

	      // Get the Gauss points on the curve.
	      curved_boundary_gauss_points ( grid, bc, XL, XR, GP );

	      // Get the total length of the arc segment.
	      curved_boundary_arclength ( grid, bc, XL, XR, &len );
	  
	      uE = 0.;
	      vE = 0.;

	      for ( k=0; k < 3; k++ )
		{
		  curved_boundary_normal_vector ( grid, bc, &(GP[k*NDIM]), &nx, &ny );

		  uE += ( w[k] * 4. * atan( GP[k*NDIM+0]/GP[k*NDIM+1] ) * nx * len );
		  vE += ( w[k] * -2. * log( GP[k*NDIM+0]*GP[k*NDIM+0] + GP[k*NDIM+1]*GP[k*NDIM+1] ) * nx * len );
		}
	  
	      // Adjust the values to account for the constant and the length.
	      uE = uE * 0.5;
	      vE = vE * 0.5;

	      u_exact[node] += uE;
	      v_exact[node] += vE;
	    }
	  
	}

      for ( i=1; i<= grid->nn; i++ )
	{
	  u_exact[i] /= grid->cv_area[i];
	  v_exact[i] /= grid->cv_area[i];
	}

      // now look at the error between the exact value and the computed value for the CV average.
      L1[0] = 0.;
      L2[0] = 0.;
      Linf[0] = 0.;

      L1[1] = 0.;
      L2[1] = 0.;
      Linf[1] = 0.;
      
      for ( i=1; i <= grid->nn; i++ )
	{
	  xmid = fabs( u_exact[i] - grid->nQ[i*NUM_VAR+1] );
	  ymid = fabs( v_exact[i] - grid->nQ[i*NUM_VAR+2] );
	  
	  L1[0] += ( xmid * grid->cv_area[i] );
	  L2[0] += ( xmid * xmid * grid->cv_area[i] );
	  Linf[0]  = MAX( xmid , Linf[0] );

	  L1[1] += ( ymid * grid->cv_area[i] );
	  L2[1] += ( ymid * ymid * grid->cv_area[i] );
	  Linf[1]  = MAX( ymid , Linf[1] );
	}

      L1[0] = L1[0] / tA;
      L2[0] = L2[0] / tA;
      L2[0] = sqrt(L2[0]);

      L1[1] = L1[1] / tA;
      L2[1] = L2[1] / tA;
      L2[1] = sqrt(L2[1]);

      
      printf("\n\n");
      printf("Error analysis. Contour integral of velocity components.  Calculated as ( |Qn - Qe|^p *dA / tA )^1/p\n");
      printf("X velocity:\n");
      printf("  L1 norm = %.15e\n",L1[0]);
      printf("  L2 norm = %.15e\n",L2[0]);
      printf("  L-infinity norm = %.15e\n\n",Linf[0]);

      
      printf("Y velocity:\n");
      printf("  L1 norm = %.15e\n",L1[1]);
      printf("  L2 norm = %.15e\n",L2[1]);
      printf("  L-infinity norm = %.15e\n\n",Linf[1]);

      freenull(u_exact);
      freenull(v_exact);
    }

  freenull(error);

  return;
}


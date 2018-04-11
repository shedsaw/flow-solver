//=============================================================
// 
//  cv_calc.C
//  
//  Functions to calculate various information about the control
//  volumes in the mesh.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "mesh_utils.h"
#include "Defined_Variables.h"
#include "linear_algebra.h"
#include "curved_boundaries.h"
#include "cv_calc.h"

//=============================================================
// 
//  cv_calc_area_Green()
//
//  Computes the area for control volumes using
//  Green's Theorem applied to edges created by the median dual
//  configuration.
//  
//  GRID *grid;                       // The grid.
//
//=============================================================

void cv_calc_area_Green ( GRID *grid )
{
  int i, k, b, s, n;                     // Loop counters.
  double nx;                             // X component of a normal.
  double ny;
  double xc, qc;                         // Values at the centroid.
  double xmid, ymid, qmid;               // Values at face midpoints.
  double len;                            // Length of an edge piece.
  int n0, n1;                            // Boundary segment node ids.
  int nodeL, nodeR;                      // Nodes comprising an edge.
  int bc;                                // Boundary condition on segment.
  int iedge,gelem;

  // Curved boundary stuff.
  double xL[2],xR[2],xM[2],GP[6];
  double w[3];
  double Ftemp;

  //FILE *fp = NULL;

  //fp = fopen( "naca0012_gauss_node_normals.dat","w");

  w[0] = 8./9.;
  w[1] = 5./9.;
  w[2] = w[1];

  // Check for memory.
  if ( grid->cv_area == NULL )
    {
      grid->cv_area = (double*)malloc((grid->nn + 1)*sizeof(double));
      if ( grid->cv_area == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'cv_area'.\n"); exit(1); }
    }

  // Make sure that cv_area is clean since area contribtutions are accumulated.
  for ( n=1; n <= grid->nn; n++ )
    {
      grid->cv_area[n] = 0.;
    }

  // Green's Theorem is applied with Q (the arbitrary function) as Q=x*i + 0*j so that del*Q=1 and with
  // Phi=1 so that del*Phi=0. Using a simple midpoint quadrature, this reduces the area computation to:
  // area = x * nx, where nx is the x component of the normal vector pointing out of the cv at that edge.

  // New method over the subedges to get the correct cv_area.

  for ( i=1; i <= grid->nsubedges; i++ )
    {
      // Find the global edge associated with the subedge.
      iedge = grid->subedges[i*2];
      
      // Find the global element.
      gelem = grid->subedges[i*2+1];
      
      // Now get the nodes attached to the edge.
      nodeL = grid->edges[2*iedge+0];
      nodeR = grid->edges[2*iedge+1];

      // Get the element centroid.
      xc = grid->el_cent[gelem*2];

      nx = grid->xn_subedges[i*3+0];
      ny = grid->xn_subedges[i*3+1];
      len = grid->xn_subedges[i*3+2];

      xmid = grid->xm_subedges[i*NDIM+0];
      ymid = grid->xm_subedges[i*NDIM+1];

      qc = xc;

      qmid = xmid;
      
      // Accumulate the area to the two bordering control volumes.
      grid->cv_area[nodeL] += (qc+qmid)/2. * nx*len;
      grid->cv_area[nodeR] -= (qc+qmid)/2. * nx*len;

    }

  // The above loops finished off the cv's completely surrounded by other cv's. The cv's touching
  // boundaries need to be closed off.
  for ( b=1; b <= grid->nb; b++ )
    {                                                       // Boundary loop.
      for ( s=1; s <= grid->nbs[b]; s++ )
	{                     
	  n0 = grid->bs[b][s][0];
	  n1 = grid->bs[b][s][1];
	  bc = grid->bbc[b];

	  if ( bc <= 10 ) // Linear boundaries.
	    {
	      xmid = (grid->x[NDIM*n0]+grid->x[NDIM*n1])/2.;
	      ymid = (grid->x[NDIM*n0+1]+grid->x[NDIM*n1+1])/2.;
	  
	      nx = grid->x[NDIM*n1+1] - ymid;
	      
	      // The 0.75 and 0.25 weights where suggested by Dr. Anderson so that this process would recover
	      // a linear distribution through the mesh.
	      grid->cv_area[n0] += (0.75*grid->x[NDIM*n0]+0.25*grid->x[NDIM*n1])*nx;
	      grid->cv_area[n1] += (0.75*grid->x[NDIM*n1]+0.25*grid->x[NDIM*n0])*nx;

	    }
	  else // Process curved boundary.
	    {
	      xL[0] = grid->x[NDIM*n0+0];
	      xL[1] = grid->x[NDIM*n0+1];

	      xR[0] = grid->x[NDIM*n1+0];
	      xR[1] = grid->x[NDIM*n1+1];

	      // Get the midpoint.
	      curved_boundary_midpoint ( grid, bc, xL, xR, xM );

	      // Do the first part - xL to xM
	      Ftemp = 0.;

	      // Get the Gauss points on the curve.
	      curved_boundary_gauss_points ( grid, bc, xL, xM, GP );

	      for ( k=0; k < 3; k++ )
		{
		  // Get the normal vector.
		  curved_boundary_normal_vector ( grid, bc, &(GP[k*NDIM]), &nx, &ny );

		  Ftemp += ( w[k] * GP[k*NDIM] * nx );

		  // write out the normal vector to a file to be read by gridgen.
		  //fprintf(fp,"2\n");
		  //fprintf(fp,"%.15e %.15e %.15e\n", GP[k*NDIM+0],GP[k*NDIM+1], 0.0);
		  //fprintf(fp,"%.15e %.15e %.15e\n", ( GP[k*NDIM+0] + 0.005*nx ),
		  //                            ( GP[k*NDIM+1] + 0.005*ny ), 0.0);
		}

	      // Now multiply by the arclength.
	      curved_boundary_arclength ( grid, bc, xL, xM, &len );

	      Ftemp = Ftemp * len * 0.5;
	      grid->cv_area[n0] += Ftemp;

	      // Do the same for the second piece - xM to xR
	      Ftemp = 0.;

	      // Get the Gauss points on the curve.
	      curved_boundary_gauss_points ( grid, bc, xM, xR, GP );

	      for ( k=0; k < 3; k++ )
		{
		  // Get the normal vector.
		  curved_boundary_normal_vector ( grid, bc, &(GP[k*NDIM]), &nx, &ny );

		  Ftemp += ( w[k] * GP[k*NDIM] * nx );

		  // write out the normal vector to a file to be read by gridgen.
		  //fprintf(fp,"2\n");
		  //fprintf(fp,"%.15e %.15e %.15e\n", GP[k*NDIM+0],GP[k*NDIM+1], 0.0);
		  //fprintf(fp,"%.15e %.15e %.15e\n", ( GP[k*NDIM+0] + 0.005*nx ),
		  //                                ( GP[k*NDIM+1] + 0.005*ny ), 0.0);
		}

	      // Now multiply by the arclength.
	      curved_boundary_arclength ( grid, bc, xM, xR, &len );

	      Ftemp = Ftemp * len * 0.5;
	      grid->cv_area[n1] += Ftemp;

	    }
	}
    }

  //fclose(fp);


  //Reuse nx to accumulate total area covered by all control volumes.
  nx = 0.;
  for ( n=1; n <= grid->nn; n++ )
    {
      nx += grid->cv_area[n];
    }

  printf("//============================================\n");
  printf("  Total area of mesh is < %.15e >.\n\n",nx);

  // Debug code to print to screen.
  // #if DEBUG
  #if DEBUG
  printf("CV AREA GREEN DEBUG:\n");
  for ( n=1; n <= grid->nn; n++ )
    {
      printf("%d : %.15e\n",n,grid->cv_area[n]);
    }
  #endif

  return;
}



//  NO CURVES BOUNDARY SUPPORT !!!!
//=============================================================
// 
//  cv_calc_area_Tri_Decomp()
//
//  Computes the area for control volumes by
//  decomposing the region of an element belonging to the control
//  volume into triangles and summing that area.
//  
//  GRID *grid                         // The grid.
//
//=============================================================

void cv_calc_area_Tri_Decomp ( GRID *grid )
{
  int i, j, k, n;                        // Loop counters.
  int node[MAX_NUM_VERT];                // Element vertices.
  int num_vert;                          // Number of vertices for an element.
  int num_edge;                          // Number of edges for an element.
  int left_id;                           // ID of left edge.
  int right_id;                          // ID of right edge.
  double xc[NDIM];                       // Vectors.
  double xmidL[NDIM];
  double xmidR[NDIM];
  double v_left[NDIM];
  double v_right[NDIM];
  double v_middle[NDIM];

  // Check for memory.
  if ( grid->cv_area == NULL )
    {
      grid->cv_area = (double*)malloc((grid->nn + 1)*sizeof(double));
      if ( grid->cv_area == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'cv_area'.\n"); exit(1); }
    }

  // Make sure that cv_area is clean since area contribtutions are accumulated.
  for ( n=1; n <= grid->nn; n++ )
    {
      grid->cv_area[n] = 0.;
    }


  // To calculate area, three vectors will be created for every node in an element. The vectors will be
  // appropriately crossed (left cross middle, middle cross right). The areas of the two triangles these
  // vectors represent will be accumulated to each node.

  // Start the process by looping over all the elements and looping through the internal median dual pieces.
  for ( i=Tri; i <= Quad; i++ )                 // Element types loop.
    {
      // Find the number of nodes and edges for this element type.
      num_vert = NumberNodesForElement(i);
      num_edge = NumberEdgesForElement(i);

      for ( j=1; j <= grid->num_elem[i]; j++ )              // Element loop.
	{
	  // Retrieve the element node ids. This is for convenience.
	  for ( k=0; k < num_vert; k++ )              // Element node extraction loop.
	    {
	      node[k] = grid->c2n[i][j*num_vert+k];
	    }                                         // End element node extraction loop.

	  // Find the element centroid.
	  xc[0] = 0.; xc[1] = 0.;
	  for ( k=0; k < num_vert; k++ )
	    {
	      for ( n=0; n < NDIM; n++ )
		{
		  xc[n] += grid->x[NDIM*node[k]+n];
		}
	    }
	  xc[0] /= (num_vert*1.);
	  xc[1] /= (num_vert*1.);

	  // Loop through each node in the element.
	  for ( k=0; k < num_vert; k++ )              // Element node loop.
	    {
	      // Get the two neighboring vertices.
	      left_id = (k+1)%num_vert;
	      right_id = (k+num_vert-1)%num_vert;

	      left_id = grid->c2n[i][j*num_vert+left_id];
	      right_id = grid->c2n[i][j*num_vert+right_id];
	      
	      // Get the edge midpoints.
	      for ( n=0; n < NDIM; n++ )
		{
		  xmidL[n] = grid->x[NDIM*node[k]+n] + grid->x[NDIM*left_id+n];
		  xmidR[n] = grid->x[NDIM*node[k]+n] + grid->x[NDIM*right_id+n];
		}

	       for ( n=0; n < NDIM; n++ )
		{
		  xmidL[n] *= 0.5;
		  xmidR[n] *= 0.5;
		}

	      // Make the vectors.
	      for ( n=0; n < NDIM; n++ )
		{
		  v_left[n] = xmidL[n] - grid->x[NDIM*node[k]+n];
		  v_right[n] = xmidR[n] - grid->x[NDIM*node[k]+n];
		  v_middle[n] = xc[n] - grid->x[NDIM*node[k]+n];
		}

	      // Accumulate the area to the control volume.
	      grid->cv_area[ node[k] ] += 0.5*( (v_left[0]*v_middle[1] - v_left[1]*v_middle[0]) +
						(v_middle[0]*v_right[1] - v_middle[1]*v_right[0]) );
	      
	    }                                       // End element node loop.

	}                                           // End element loop.

    }                                               // End element type loop.

      
  //Reuse xc to accumulate total area covered by all control volumes.
  xc[0] = 0.;
  for ( n=1; n <= grid->nn; n++ )
    {
      xc[0] += grid->cv_area[n];
    }
  
  printf("//============================================\n");
  printf("  Total area of mesh is < %.15e >.\n\n",xc[0]);
  
  // Debug code to print to screen.
  #if DEBUG
  printf("CV AREA TRI DECOMP DEBUG:\n");
  for ( n=1; n <= grid->nn; n++ )
    {
      printf("%d : %.15e\n",n,grid->cv_area[n]);
    }
  #endif
  
  return;
}



//=============================================================
// 
//  cv_calc_check_closure()
//
//  Loop across the control volumes and checks that the sum of
//  normal vectors is zero, indicating a closed volume (area).
//
//  GRID *grid;                       // The grid.  
//
//=============================================================

void cv_calc_check_closure( GRID *grid )
{
  int i, j, k, b, s;                     // Loop counters.
  double *v;                             // Array of vectors associated with each control volume.
  double xmid, ymid;                     // Values at face midpoints.
  double xc, yc;                         // Element centroid coordinates.
  double nx,ny;                          // Normal vector components.
  int n0, n1;                            // Boundary segment node ids.
  double mag;                            // Magnitude of the vector made from the sum around edges.
  double max_leak=0.;                    // Largest leak found in the mesh.
  int cv_leak=0;                         // Which control volume has the largest leak.
  int nodeL, nodeR;                      // Nodes comprising an edge.
  int bc;                                // Boundary condition on segment.
  int iedge, gelem;

  // Curved boundary stuff.
  double xL[2],xR[2],xM[2],GP[6];
  double w[3];
  double len;
  double Ftempx, Ftempy;

  w[0] = 8./9.;
  w[1] = 5./9.;
  w[2] = w[1];

  // Allocate space for vector sum for each control volume. Use calloc since normal vectors
  // are summed to this.
  v = (double*)calloc(NDIM*(grid->nn+1),sizeof(double));
  if ( v == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'v' in disclosure check.\n"); exit(1); }

  for ( i=1; i <= grid->nn; i++ )
    {
      for ( j=0; j < NDIM; j++ )
	v[i*NDIM+j] = 0.;
    }

  // Use the subedges to get the needed data.
  for ( i=1; i <= grid->nsubedges; i++ )
    {
      // Find the global edge associated with the subedge.
      iedge = grid->subedges[i*2];
      
      // Find the global element.
      gelem = grid->subedges[i*2+1];
      
      // Now get the nodes attached to the edge.
      nodeL = grid->edges[2*iedge+0];
      nodeR = grid->edges[2*iedge+1];

      // Get the element centroid.
      xc = grid->el_cent[gelem*2];
      yc = grid->el_cent[gelem*2+1];

      nx = grid->xn_subedges[i*3+0];
      ny = grid->xn_subedges[i*3+1];
      len = grid->xn_subedges[i*3+2];

      xmid = grid->xm_subedges[i*NDIM+0];
      ymid = grid->xm_subedges[i*NDIM+1];

      // Accumulate the normal vector to the control volumes.
      nx = nx*len;
      ny = ny*len;

      v[nodeL*NDIM+0] += nx;
      v[nodeL*NDIM+1] += ny;

      v[nodeR*NDIM+0] -= nx;
      v[nodeR*NDIM+1] -= ny;
   
    }

  /*
  // Start the process by looping over all the elements and looping through the internal median dual pieces.
  for ( i=Tri; i <= Quad; i++ )                 // Element types loop.
    {
      // Find the number of nodes and edges for this element type.
      num_vert = NumberNodesForElement(i);
      num_edge = NumberEdgesForElement(i);

      for ( j=1; j <= grid->num_elem[i]; j++ )              // Element loop.
	{
	  // Retrieve the element node ids.
	  for ( k=0; k < num_vert; k++ )              // Element node extraction loop.
	    {
	      node[k] = grid->c2n[i][j*num_vert+k];
	    }                                         // End element node extraction loop.

	  // Find the element centroid.
	  xc = 0.;  yc = 0.;
	  for ( k=0; k < num_vert; k++ )
	    {
	      xc += grid->x[NDIM*node[k]];
	      yc += grid->x[NDIM*node[k]+1];
	    }
	  xc /= (num_vert*1.);
	  yc /= (num_vert*1.);

	  // Loop through the nodes that make up that element since each represent a distinct control volume. 
	  // This process will be similar to the area calculations above with Green's Theorem. For each edge
	  // normal, it will be added to one cv and subtracted from the other.
	  for ( k=0; k < num_vert; k++ )              // Element node loop.
	    {
	      // In 2D the local edge id needed is the same as the local node id (i.e. 0,1,..).
	      // Get the edge ID.
	      edge_id = grid->c2e[i][j*num_edge+k];

	      // Get the edges nodes.
	      nodeL = grid->edges[edge_id*2];
	      nodeR = grid->edges[edge_id*2+1];

	      // Calculate the edge midpoint.
	      xmid = 0.5*( grid->x[nodeL*NDIM+0] +grid->x[nodeR*NDIM+0] );
	      ymid = 0.5*( grid->x[nodeL*NDIM+1] +grid->x[nodeR*NDIM+1] );
	      
	      // Make the normal vector.
	      nx = yc-ymid;
	      ny = xmid-xc;

	      // Accumulate the normal vector to the two bordering control volumes.
	      v[ node[k]*NDIM ] += nx;
	      v[ node[k]*NDIM + 1 ] += ny;

	      v[ (node[(k+1)%num_vert])*NDIM ] -= nx;        // The modulus logic will ensure that when
              v[ (node[(k+1)%num_vert])*NDIM + 1 ] -= ny  ;  // processing the last node, the bordering
                                                             // cv is local node 0.

	    }                                       // End element node loop.

	}                                           // End element loop.

    }                                               // End element type loop.
  */

  // The above loops finished off the cv's completely surrounded by other cv's. The cv's touching
  // boundaries need to be closed off.
  for ( b=1; b <= grid->nb; b++ )                                         // Boundary loop.
    {                                             
      for ( s=1; s <= grid->nbs[b]; s++ )                                 // Edge loop.
	{                                  
	  n0 = grid->bs[b][s][0];
	  n1 = grid->bs[b][s][1];
	  bc = grid->bbc[b];

	  if ( bc <= 10 )
	    {
	      xmid = (grid->x[NDIM*n0]+grid->x[NDIM*n1])/2.;
	      ymid = (grid->x[NDIM*n0+1]+grid->x[NDIM*n1+1])/2.;

	      // Get the normal and add it to the control volumes.
	      nx = grid->x[NDIM*n1+1] - ymid;
	      ny = xmid - grid->x[NDIM*n1];
	  
	      v[n0*NDIM] += nx;
	      v[n0*NDIM+1] += ny;
	  
	      v[n1*NDIM] += nx;
	      v[n1*NDIM+1] += ny;

	    }
	  else // Process curved boundary.
	    {
	      xL[0] = grid->x[NDIM*n0+0];
	      xL[1] = grid->x[NDIM*n0+1];

	      xR[0] = grid->x[NDIM*n1+0];
	      xR[1] = grid->x[NDIM*n1+1];

	      Ftempx = 0.;  Ftempy = 0.;

	      // Get the midpoint.
	      curved_boundary_midpoint ( grid, bc, xL, xR, xM );

	      // Get the Gauss points on the curve.
	      curved_boundary_gauss_points ( grid, bc, xL, xM, GP );

	      curved_boundary_arclength ( grid, bc, xL, xM, &len );

	      for ( k=0; k < 3; k++ )
		{
		  // Get the normal vector.
		  curved_boundary_normal_vector ( grid, bc, &(GP[k*NDIM]), &nx, &ny );

		  Ftempx += ( w[k] * nx );
		  Ftempy += ( w[k] * ny );

		  //v[n0*NDIM] += nx*len;
		  //v[n0*NDIM+1] += ny*len;
		}

	      Ftempx = Ftempx * len * 0.5;
	      Ftempy = Ftempy * len * 0.5;

	      v[n0*NDIM+0] += Ftempx;
	      v[n0*NDIM+1] += Ftempy;
	      
	      Ftempx = 0.;  Ftempy = 0.;

	      // Get the Gauss points on the curve.
	      curved_boundary_gauss_points ( grid, bc, xM, xR, GP );
	      
	      curved_boundary_arclength ( grid, bc, xM, xR, &len );

	      for ( k=0; k < 3; k++ )
		{
		  // Get the normal vector.
		  curved_boundary_normal_vector ( grid, bc, &(GP[k*NDIM]), &nx, &ny );

		  Ftempx += ( w[k] * nx );
		  Ftempy += ( w[k] * ny );

		  //v[n1*NDIM] += nx*len;
		  //v[n1*NDIM+1] += ny*len;
		}

	      Ftempx = Ftempx * len * 0.5;
	      Ftempy = Ftempy * len * 0.5;

	      v[n1*NDIM+0] += Ftempx;
	      v[n1*NDIM+1] += Ftempy;

	    }
	}                                                                   // End Edge loop.
    }                                                                       // End boundary loop.

  // Now search through the control volumes and check that there are no leaks.
  // Also, record the largest leak.
  for ( i=1; i <= grid->nn; i++ )
    {
      mag = sqrt( v[i*NDIM]*v[i*NDIM] + v[i*NDIM+1]*v[i*NDIM+1] );

      if ( mag > 1.e-14 )
	{
	  //printf("FATAL ERROR: Leak detected in control volume %d : magnitude of closure = %.15e\n",i,mag);
	  //exit(1);
	}

      if ( mag > max_leak )
	{
	  max_leak = mag;
	  cv_leak = i;
	}
    }

  printf("cv_calc_check_closure(): Magnitude of largest nonclosure = <%.15e>.\n",max_leak);
  printf("                         At control volume %d.\n",cv_leak);

  // Debug code to print to screen.
  #if DEBUG
  int n = 0;
  printf("CV CLOSURE DEBUG:\n");
  for ( n=1; n <= grid->nn; n++ )
    {
      printf("%d : %.15e   %.15e  ::  MAGNITUDE ->  %.15e\n",n,v[n*NDIM],v[n*NDIM+1],
	     sqrt( v[n*NDIM+0]*v[n*NDIM+0] + v[n*NDIM+1]*v[n*NDIM+1] ) );
    }
  #endif

  freenull(v);

  return;
}


//  NO CURVES BOUNDARY SUPPORT !!!!
//=============================================================
// 
//  cv_calc_length_scale()
//
//  Computes the length scale for the control volume using the
//  the median dual pieces.
//
//  Currently implemented:  Length scale = Area / Perimeter.
//
//  GRID *grid;                       // The grid.  
//
//=============================================================

void cv_calc_length_scale ( GRID *grid )
{
  int i, j, k, b, s, n;                  // Loop counters.
  int node[MAX_NUM_VERT];                // Element vertices.
  double xc, yc;                         // Values at the centroid.
  double xmid, ymid;                     // Values at face midpoints.
  double mag;                            // Length of an edge piece.
  int edge_id;                           // Edge ID.
  int n0, n1;                            // Boundary segment node ids.
  int num_vert;                          // Number of vertices for an element.
  int num_edge;                          // Number of edges for an element.
  int nodeL, nodeR;                      // Nodes comprising an edge.


  // Check for memory.
  if ( grid->length_scale == NULL )
    {
      grid->length_scale = (double*)malloc((grid->nn + 1)*sizeof(double));
      if ( grid->length_scale == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'length_scale'.\n"); exit(1); }
    }

  // Make sure it is clean since length contribtutions are accumulated.
  for ( n=1; n <= grid->nn; n++ )
    {
      grid->length_scale[n] = 0.;
    }
  
  // Start the process by looping over all the elements and looping through the internal median dual pieces.
  for ( i=Tri; i <= Quad; i++ )                 // Element types loop.
    {
      // Find the number of nodes and edges for this element type.
      num_vert = NumberNodesForElement(i);
      num_edge = NumberEdgesForElement(i);

      for ( j=1; j <= grid->num_elem[i]; j++ )              // Element loop.
	{
	  // Retrieve the element node ids.
	  for ( k=0; k < num_vert; k++ )              // Element node extraction loop.
	    {
	      node[k] = grid->c2n[i][j*num_vert+k];
	    }                                         // End element node extraction loop.

	  // Find the element centroid values now since does not change during the process.
	  xc = 0.;  yc = 0.;
	  for ( k=0; k < num_vert; k++ )
	    {
	      xc += grid->x[NDIM*node[k]];
	      yc += grid->x[NDIM*node[k]+1];
	    }
	  xc /= (num_vert*1.);
	  yc /= (num_vert*1.);

	  // Loop through the nodes that make up that element since each represent a distinct control volume.
	  for ( k=0; k < num_vert; k++ )              // Element node loop.
	    {
	      // In 2D the local edge id needed is the same as the local node id (i.e. 0,1,..).
	      // Get the edge ID.
	      edge_id = grid->c2e[i][j*num_edge+k];

	      // Get the edges nodes.
	      nodeL = grid->edges[edge_id*2];
	      nodeR = grid->edges[edge_id*2+1];

	      // Calculate the edge midpoint.
	      xmid = 0.5*( grid->x[nodeL*NDIM+0] +grid->x[nodeR*NDIM+0] );
	      ymid = 0.5*( grid->x[nodeL*NDIM+1] +grid->x[nodeR*NDIM+1] );

	      mag = sqrt ( (xmid-xc)*(xmid-xc) + (ymid-yc)*(ymid-yc) );

	      // Accumulate the length to the two bordering control volumes.
	      grid->length_scale[ node[k] ] += mag;
	      grid->length_scale[ node[(k+1)%num_vert] ] += mag;   // The modulus logic will ensure that when
                                                                   // processing the last node, the bordering
                                                                   // cv is local node 0.

	    }                                       // End element node loop.

	}                                           // End element loop.

    }                                               // End element type loop.

  // The above loops finished off the cv's completely surrounded by other cv's. The cv's touching
  // boundaries need to be closed off.
  for ( b=1; b <= grid->nb; b++ )
    {                                                       // Boundary loop.
      for ( s=1; s <= grid->nbs[b]; s++ )
	{                                  
	  n0 = grid->bs[b][s][0];
	  n1 = grid->bs[b][s][1];

	  xmid = (grid->x[NDIM*n0]+grid->x[NDIM*n1])/2.;
	  ymid = (grid->x[NDIM*n0+1]+grid->x[NDIM*n1+1])/2.;
	  
	  mag = sqrt ( (xmid-grid->x[NDIM*n0])*(xmid-grid->x[NDIM*n0]) +
		       (ymid-grid->x[NDIM*n0+1])*(ymid-grid->x[NDIM*n0+1]) );

	  grid->length_scale[n0] += mag;
	  grid->length_scale[n1] += mag;
	}
    }

  // Now length_scale contains only the perimeter. Process all cv's and find the length scale : Area / Perimeter.
  for ( i=1; i <= grid->nn; i++ )
    {
      grid->length_scale[i] = grid->cv_area[i] / grid->length_scale[i];
    }
  
  // Debug code to print to screen.
  #if DEBUG
  printf("LENGTH SCALE DEBUG:\n");
  for ( n=1; n <= grid->nn; n++ )
    {
      printf("%d : %.15e\n",n,grid->length_scale[n]);
    }
  #endif

  return;
}

// I MIGHT LEAVE CURVED BOUNDARY SUPPORT OUT FOR NOW SINCE THIS
// IS REALLY ONLY FOR A ROUGHT ESTIMATE.
//=============================================================
// 
//  cv_calc_eigenvalue_contribution()
//
//  Accumulates the maixmum eigenvalue time normal vector length
//  across the median dual structure of the nodes.
//
//  GRID *grid;                       // The grid.
//  PARAMS p;                         // The parameters.
//
//=============================================================

void cv_calc_eigenvalue_contribution ( GRID *grid, PARAMS p )
{
  int n;                                 // Loop counters.
  int isubedge, iedge, ibedge;           // Edge identifiers.
  int nodeL, nodeR;                      // Nodes comprising an edge.

  double nx,ny,len;                      // Normal vector.
  double rho,u,v,P,theta,c;              // Primitive / eigen values.
  double tpc, tmc;
  double max_ev;                         // Maximum eigenvalue.


  // Check for memory.
  if ( grid->ev_con == NULL )
    {
      grid->ev_con = (double*)malloc((grid->nn + 1)*sizeof(double));
      if ( grid->ev_con == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'ev_con'.\n"); exit(1); }
    }

  // Make sure it is clean since length contribtutions are accumulated.
  for ( n=1; n <= grid->nn; n++ )
    {
      grid->ev_con[n] = 0.;
    }
  
  
  // To compute the accumulation of the eigenvalues contributions, I will loop over the subedge data structure
  // for the interior pieces and then loop over the boundary edges to complete the mesh.

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

      // Calculate the primitive state for nodeL.

      if ( RECON_PRIM )
	{
	  rho = grid->nQ[nodeL*NUM_VAR + 0];
	  u =   grid->nQ[nodeL*NUM_VAR + 1];
	  v =   grid->nQ[nodeL*NUM_VAR + 2];
	  P =   grid->nQ[nodeL*NUM_VAR + 3];
	}
      else
	{
	  rho = grid->nQ[nodeL*NUM_VAR + 0];
	  u = grid->nQ[nodeL*NUM_VAR + 1] / rho;
	  v = grid->nQ[nodeL*NUM_VAR + 2] / rho;
	  P = (p.gamma - 1.0)*( grid->nQ[nodeL*NUM_VAR + 3] - 0.5*rho*(u*u + v*v) );
	}

      c = sqrt( (p.gamma*P)/rho );

      theta = u*nx + v*ny;
      tpc = theta + c;
      tmc = theta - c;

      tpc = fabs(tpc);
      tmc = fabs(tmc);
      max_ev = MAX(tpc,tmc);

      if ( !(isfinite(max_ev) ) )
	{
	  printf("Problem detected in the eigenvalue contribution of left node %d from subedge %d : theta = %.15E    c = %.15E\n",nodeL,isubedge,theta,c);
	  printf("rho = %.15E     u = %.15E     v = %.15E     P = %.15E\n",rho,u,v,P);
	  fflush(stdout);
	  exit(0);
	}

      grid->ev_con[nodeL] += ( max_ev * len );

      // Calculate the primitive state for nodeR.

      if ( RECON_PRIM )
	{
	  rho = grid->nQ[nodeR*NUM_VAR + 0];
	  u =   grid->nQ[nodeR*NUM_VAR + 1];
	  v =   grid->nQ[nodeR*NUM_VAR + 2];
	  P =   grid->nQ[nodeR*NUM_VAR + 3];
	}
      else
	{
	  rho = grid->nQ[nodeR*NUM_VAR + 0];
	  u = grid->nQ[nodeR*NUM_VAR + 1] / rho;
	  v = grid->nQ[nodeR*NUM_VAR + 2] / rho;
	  P = (p.gamma - 1.0)*( grid->nQ[nodeR*NUM_VAR + 3] - 0.5*rho*(u*u + v*v) );
	}

      c = sqrt( (p.gamma*P)/rho );

      theta = u*nx + v*ny;
      tpc = theta + c;
      tmc = theta - c;

      tpc = fabs(tpc);
      tmc = fabs(tmc);
      max_ev = MAX(tpc,tmc);

      if ( !(isfinite(max_ev) ) )
	{
	  printf("Problem detected in the eigenvalue contribution of right node %d from subedge %d : theta = %.15E    c = %.15E\n",nodeR,isubedge,theta,c);
	  printf("rho = %.15E     u = %.15E     v = %.15E     P = %.15E\n",rho,u,v,P);
	  fflush(stdout);
	  exit(0);
	}

      grid->ev_con[nodeR] += ( max_ev * len );
    }
  
  for ( ibedge=1; ibedge <= grid->nbedges; ibedge++ )
    {
      // Extract the real node on the edge.
      nodeL = grid->bedges[ibedge*5];

      // Get the normal vector.
      nx = grid->xn_bedges[ibedge*3 + 0];
      ny = grid->xn_bedges[ibedge*3 + 1];
      len = grid->xn_bedges[ibedge*3 + 2];

      // Calculate the primitive state for the node.
      if ( RECON_PRIM )
	{
	  rho = grid->nQ[nodeL*NUM_VAR + 0];
	  u =   grid->nQ[nodeL*NUM_VAR + 1];
	  v =   grid->nQ[nodeL*NUM_VAR + 2];
	  P =   grid->nQ[nodeL*NUM_VAR + 3];
	}
      else
	{
	  rho = grid->nQ[nodeL*NUM_VAR + 0];
	  u = grid->nQ[nodeL*NUM_VAR + 1] / rho;
	  v = grid->nQ[nodeL*NUM_VAR + 2] / rho;
	  P = (p.gamma - 1.0)*( grid->nQ[nodeL*NUM_VAR + 3] - 0.5*rho*(u*u + v*v) );
	}

      c = sqrt( (p.gamma*P)/rho );

      theta = u*nx + v*ny;
      tpc = theta + c;
      tmc = theta - c;

      tpc = fabs(tpc);
      tmc = fabs(tmc);
      max_ev = MAX(tpc,tmc);

      if ( !(isfinite(max_ev) ) )
	{
	  printf("Problem detected in the eigenvalue contribution of boundary node %d from bedge %d : theta = %.15E    c = %.15E\n",nodeL,ibedge,theta,c);
	  printf("rho = %.15E     u = %.15E     v = %.15E     P = %.15E\n",rho,u,v,P);
	  fflush(stdout);
	  exit(0);
	}

      grid->ev_con[nodeL] += ( max_ev * len );
    }

  // Debug code to print to screen.
  #if DEBUG
  printf("EIGENVALUE CONTRIBUTION DEBUG:\n");
  for ( n=1; n <= grid->nn; n++ )
    {
      printf("%d : %.15e\n",n,grid->ev_con[n]);
    }
  #endif

  return;
}



//=============================================================
// 
//  cv_calc_Moments()
//
//  Computes the control volume moments that appear in the reconstruction
//  scheme ( from Ollivier-Gooch ). Use Gaussian quadrature rule to
//  get high order accuracy.
//
//  GRID *grid;                       // The grid.  
//
//=============================================================

void cv_calc_Moments ( GRID *grid )
{
  int i, j, k, b, s, n;                  // Loop counters.
  double xc, yc;                         // Values at the centroid.
  double xmid, ymid;                     // Values at face midpoints.
  double len;                            // Length of an edge piece.
  double nx;                             // Normal vector component.
  double ny;
  double x1,x2,x3;                       // Integration points.
  double y1,y2,y3;
  double xL,xR,yL,yR;
  double t1,t2,t3;                       // Gaussian roots.
  double w1,w2,w3;                       // Gaussian coefficients.
  double xi, yi;                         // Node coordinates.
  double Ix,Iy,Ixx,Iyy,Ixy;              // Moment terms.
  double Ixxx,Iyyy,Ixxy,Ixyy;
  int iedge;                             // Edge ID.
  int num_vert;                          // Number of vertices for an element.
  int nodeL, nodeR;                      // Nodes comprising an edge.
  int node;                              // Generic node.
  int gelem;                             // Global element number.

  int bc;
  double XL[2],XR[2],XM[2],GP[6],w[3],XI[2];

  
  // DEBUG STUFF
  double *moments_area_integral = NULL;

  int e,v,ind;
  int nodes[MAX_NUM_VERT];               // Element vertices.
  int left_id;                           // ID of left edge.
  int right_id;                          // ID of right edge.
  double XC[NDIM];                       // Vectors.
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
  double QTRI[NUM_MOM];
  double QGP[NUM_MOM];

  double gp[7][3];                       // Gauss points for integration.
  double dA;                             // Differential area (for a triangle).

  int ct_flag;
  int ghost_node;
  int seg;
  double ct_gp[7][2];                    // Gauss integration points in (r,s) space (xi,eta).
  double ct_nodes[7][2];                 // The x,y locations of the 6 nodes defining the p2 triangle (curved side).
  double N1,N2,N3,N4,N5,N6;              // Shape functions.
  double N1r,N1s,N2r,N2s,N3r,N3s;        // Shape function derivatives.
  double N4r,N4s,N5r,N5s,N6r,N6s;
  double xr,xs,yr,ys,jacobian;           // Transformation jacobian.
  double zi,eta;                         // Triangle coordinates.

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
  /////////////////////////////////////////////////////////////////////////////

  w[0] = 8./9.;
  w[1] = 5./9.;
  w[2] = 5./9.;

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

  // Check for memory.
  if ( grid->Moments == NULL )
    {
      grid->Moments = (double*)malloc((grid->nn + 1)*NUM_MOM*sizeof(double));
      if ( grid->Moments == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'Moments'.\n"); exit(1); }
    }

  moments_area_integral = (double*)malloc((grid->nn+1)*NUM_MOM*sizeof(double));
  if ( moments_area_integral == NULL ) { printf("Could not allocate moments_area_integral.\n"); exit(1); }

  // Make sure it is clean since length contribtutions are accumulated.
  for ( n=1; n <= grid->nn; n++ )
    {
      for ( i=0; i < NUM_MOM; i++ )
	{
	  grid->Moments[n*NUM_MOM+i] = 0.;
	  moments_area_integral[n*NUM_MOM+i] = 0.;
	}
    }

  // The process is to:
  // 1. Loop over all the subedges which each represent a dual edge.
  // 2. Reconstruct the dual edge from the available information.
  // 3. Apply Gaussian quadrature to calculate the moment terms for each node (left and right).

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
      xc = grid->el_cent[gelem*2];
      yc = grid->el_cent[gelem*2+1];

      // Get the edge midpoint.
      //xmid = 0.5*( grid->x[nodeL*NDIM+0] + grid->x[nodeR*NDIM+0] );
      //ymid = 0.5*( grid->x[nodeL*NDIM+1] + grid->x[nodeR*NDIM+1] );

      nx = grid->xn_subedges[i*3+0];
      ny = grid->xn_subedges[i*3+1];
      len = grid->xn_subedges[i*3+2];

      //xmid = ny*len + xc;
      //ymid = yc - nx*len;            // should return the correct point on the curved boundary.

      xmid = grid->xm_subedges[i*NDIM+0];
      ymid = grid->xm_subedges[i*NDIM+1];

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
      w1 = 5./9.;

      t2 = 0.;
      w2 = 8./9.;

      t3 = sqrt(15.0)/(-5.0);
      w3 = w1;

      x1 = (xmid+xc)*0.5 + (xc-xmid)*0.5*t1;
      y1 = (ymid+yc)*0.5 + (yc-ymid)*0.5*t1;
      
      x2 = (xmid+xc)*0.5 + (xc-xmid)*0.5*t2;
      y2 = (ymid+yc)*0.5 + (yc-ymid)*0.5*t2;
      
      x3 = (xmid+xc)*0.5 + (xc-xmid)*0.5*t3;
      y3 = (ymid+yc)*0.5 + (yc-ymid)*0.5*t3;

      // Now I can apply the function from each of the moment pieces to the Gauss points.
    
      // 1/2 (x - x_i)^2
      Ix = w1*(x1-xi)*(x1-xi) +
	   w2*(x2-xi)*(x2-xi) +
	   w3*(x3-xi)*(x3-xi);

      Ix = 0.5 * Ix * nx * len/2.;

      // (x - x_i)(y - y_i)
      Iy = w1*(x1-xi)*(y1-yi) +
	   w2*(x2-xi)*(y2-yi) +
	   w3*(x3-xi)*(y3-yi);

      Iy = Iy * nx * len/2.;

      // 1/3 (x-x_i)^3
      Ixx = w1*(x1-xi)*(x1-xi)*(x1-xi) +
	    w2*(x2-xi)*(x2-xi)*(x2-xi) +
	    w3*(x3-xi)*(x3-xi)*(x3-xi);

      Ixx = (1./3.) * Ixx * nx * len/2.;

      // (x-x_i)(y-y_i)^2
      Iyy = w1*(x1-xi)*(y1-yi)*(y1-yi) +
	    w2*(x2-xi)*(y2-yi)*(y2-yi) +
	    w3*(x3-xi)*(y3-yi)*(y3-yi);

      Iyy = Iyy * nx * len/2.;

      // 1/2 (x-x_i)^2(y-y_i)
      Ixy = w1*(x1-xi)*(x1-xi)*(y1-yi) +
	    w2*(x2-xi)*(x2-xi)*(y2-yi) +
	    w3*(x3-xi)*(x3-xi)*(y3-yi);

      Ixy = 0.5 * Ixy * nx * len/2.;

      // (x-x_i)^3 -> 1/4 (x-x_i)^4
      Ixxx = w1*(x1-xi)*(x1-xi)*(x1-xi)*(x1-xi) +
             w2*(x2-xi)*(x2-xi)*(x2-xi)*(x2-xi) +
             w3*(x3-xi)*(x3-xi)*(x3-xi)*(x3-xi);

      Ixxx = (0.25) * Ixxx * nx * len/2.;

      // (y-y_i)^3 -> (x-x_i)*(y-y_i)^3
      Iyyy = w1*(x1-xi)*(y1-yi)*(y1-yi)*(y1-yi) +
             w2*(x2-xi)*(y2-yi)*(y2-yi)*(y2-yi) +
             w3*(x3-xi)*(y3-yi)*(y3-yi)*(y3-yi);

      Iyyy = Iyyy * nx * len/2.;

      // (x-x_i)^2*(y-y_i) -> 1/3 (x-x_i)^3(y-y_i)
      Ixxy = w1*(x1-xi)*(x1-xi)*(x1-xi)*(y1-yi) +
             w2*(x2-xi)*(x2-xi)*(x2-xi)*(y2-yi) +
             w3*(x3-xi)*(x3-xi)*(x3-xi)*(y3-yi);

      Ixxy = (1./3.) * Ixxy * nx * len/2.;

      // (x-x_i)*(y-y_i)^2 -> 1/2 (x-x_i)^2(y-y_i)^2
      Ixyy = w1*(x1-xi)*(x1-xi)*(y1-yi)*(y1-yi) +
             w2*(x2-xi)*(x2-xi)*(y2-yi)*(y2-yi) +
             w3*(x3-xi)*(x3-xi)*(y3-yi)*(y3-yi);

      Ixyy = 0.5 * Ixyy * nx * len/2.;


      // Accumulate the moments to the node.
      grid->Moments[nodeL*NUM_MOM+0] += Ix;
      grid->Moments[nodeL*NUM_MOM+1] += Iy;
      grid->Moments[nodeL*NUM_MOM+2] += Ixx;
      grid->Moments[nodeL*NUM_MOM+3] += Iyy;
      grid->Moments[nodeL*NUM_MOM+4] += Ixy;
      grid->Moments[nodeL*NUM_MOM+5] += Ixxx;
      grid->Moments[nodeL*NUM_MOM+6] += Iyyy;
      grid->Moments[nodeL*NUM_MOM+7] += Ixxy;
      grid->Moments[nodeL*NUM_MOM+8] += Ixyy;


      // Redo the calculations for the right node.
      // Get the node coordinates.
      xi = grid->x[nodeR*NDIM+0];
      yi = grid->x[nodeR*NDIM+1];

      // Now I can apply the function from each of the moment pieces to the Gauss points.
    
      // 1/2 (x - x_i)^2
      Ix = w1*(x1-xi)*(x1-xi) +
	   w2*(x2-xi)*(x2-xi) +
	   w3*(x3-xi)*(x3-xi);

      Ix = 0.5 * Ix * nx * len/2.;

      // (x - x_i)(y - y_i)
      Iy = w1*(x1-xi)*(y1-yi) +
	   w2*(x2-xi)*(y2-yi) +
	   w3*(x3-xi)*(y3-yi);

      Iy = Iy * nx * len/2.;

      // 1/3 (x-x_i)^3
      Ixx = w1*(x1-xi)*(x1-xi)*(x1-xi) +
	    w2*(x2-xi)*(x2-xi)*(x2-xi) +
	    w3*(x3-xi)*(x3-xi)*(x3-xi);

      Ixx = (1./3.) * Ixx * nx * len/2.;

      // (x-x_i)(y-y_i)^2
      Iyy = w1*(x1-xi)*(y1-yi)*(y1-yi) +
	    w2*(x2-xi)*(y2-yi)*(y2-yi) +
	    w3*(x3-xi)*(y3-yi)*(y3-yi);

      Iyy = Iyy * nx * len/2.;

      // 1/2 (x-x_i)^2(y-y_i)
      Ixy = w1*(x1-xi)*(x1-xi)*(y1-yi) +
	    w2*(x2-xi)*(x2-xi)*(y2-yi) +
	    w3*(x3-xi)*(x3-xi)*(y3-yi);

      Ixy = 0.5 * Ixy * nx * len/2.;

      // (x-x_i)^3 -> 1/4 (x-x_i)^4
      Ixxx = w1*(x1-xi)*(x1-xi)*(x1-xi)*(x1-xi) +
             w2*(x2-xi)*(x2-xi)*(x2-xi)*(x2-xi) +
             w3*(x3-xi)*(x3-xi)*(x3-xi)*(x3-xi);

      Ixxx = (0.25) * Ixxx * nx * len/2.;

      // (y-y_i)^3 -> (x-x_i)*(y-y_i)^3
      Iyyy = w1*(x1-xi)*(y1-yi)*(y1-yi)*(y1-yi) +
             w2*(x2-xi)*(y2-yi)*(y2-yi)*(y2-yi) +
             w3*(x3-xi)*(y3-yi)*(y3-yi)*(y3-yi);

      Iyyy = Iyyy * nx * len/2.;

      // (x-x_i)^2*(y-y_i) -> 1/3 (x-x_i)^3(y-y_i)
      Ixxy = w1*(x1-xi)*(x1-xi)*(x1-xi)*(y1-yi) +
             w2*(x2-xi)*(x2-xi)*(x2-xi)*(y2-yi) +
             w3*(x3-xi)*(x3-xi)*(x3-xi)*(y3-yi);

      Ixxy = (1./3.) * Ixxy * nx * len/2.;

      // (x-x_i)*(y-y_i)^2 -> 1/2 (x-x_i)^2(y-y_i)^2
      Ixyy = w1*(x1-xi)*(x1-xi)*(y1-yi)*(y1-yi) +
             w2*(x2-xi)*(x2-xi)*(y2-yi)*(y2-yi) +
             w3*(x3-xi)*(x3-xi)*(y3-yi)*(y3-yi);

      Ixyy = 0.5 * Ixyy * nx * len/2.;

      // Accumulate the moments to the nodes.
      grid->Moments[nodeR*NUM_MOM+0] -= Ix;
      grid->Moments[nodeR*NUM_MOM+1] -= Iy;
      grid->Moments[nodeR*NUM_MOM+2] -= Ixx;
      grid->Moments[nodeR*NUM_MOM+3] -= Iyy;
      grid->Moments[nodeR*NUM_MOM+4] -= Ixy;
      grid->Moments[nodeR*NUM_MOM+5] -= Ixxx;
      grid->Moments[nodeR*NUM_MOM+6] -= Iyyy;
      grid->Moments[nodeR*NUM_MOM+7] -= Ixxy;
      grid->Moments[nodeR*NUM_MOM+8] -= Ixyy;
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

	  // For Gaussian quadrature I need to map the points in t-space
	  // onto my dual edge.

	  // I'm also pretty sure that the to be consistent with the right hand
	  // rule the boundary edge is defined from node0 to node1 on the segment.
	  // The normal should thus point out of the mesh.

	  // To get the x,y points from t space use - x = (xL+xR)/2 + (xR-xL)/2*t
	  //                                          y = (yL+yR)/2 + (yR-yL)/2*t
	  // The left and right nodes depend on which node the node is. Is it on the
	  // beginning or end of the edge.

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

	  // Now I can apply the function from each of the moment pieces to the Gauss points.
	  
	  // 1/2 (x - x_i)^2
	  Ix = w1*(x1-xi)*(x1-xi) +
	       w2*(x2-xi)*(x2-xi) +
	       w3*(x3-xi)*(x3-xi);
	  
	  Ix = 0.5 * Ix * nx * len/2.;

	  // (x - x_i)(y - y_i)
	  Iy = w1*(x1-xi)*(y1-yi) +
	       w2*(x2-xi)*(y2-yi) +
	       w3*(x3-xi)*(y3-yi);

	  Iy = Iy * nx * len/2.;

	  // 1/3 (x-x_i)^3
	  Ixx = w1*(x1-xi)*(x1-xi)*(x1-xi) +
	        w2*(x2-xi)*(x2-xi)*(x2-xi) +
	        w3*(x3-xi)*(x3-xi)*(x3-xi);

	  Ixx = (1./3.) * Ixx * nx * len/2.;

	  // (x-x_i)(y-y_i)^2
	  Iyy = w1*(x1-xi)*(y1-yi)*(y1-yi) +
	        w2*(x2-xi)*(y2-yi)*(y2-yi) +
	        w3*(x3-xi)*(y3-yi)*(y3-yi);

	  Iyy = Iyy * nx * len/2.;

	  // 1/2 (x-x_i)^2(y-y_i)
	  Ixy = w1*(x1-xi)*(x1-xi)*(y1-yi) +
	        w2*(x2-xi)*(x2-xi)*(y2-yi) +
	        w3*(x3-xi)*(x3-xi)*(y3-yi);

	  Ixy = 0.5 * Ixy * nx * len/2.;

	  // (x-x_i)^3 -> 1/4 (x-x_i)^4
	  Ixxx = w1*(x1-xi)*(x1-xi)*(x1-xi)*(x1-xi) +
   	         w2*(x2-xi)*(x2-xi)*(x2-xi)*(x2-xi) +
	         w3*(x3-xi)*(x3-xi)*(x3-xi)*(x3-xi);

	  Ixxx = (0.25) * Ixxx * nx * len/2.;

	  // (y-y_i)^3 -> (x-x_i)*(y-y_i)^3
	  Iyyy = w1*(x1-xi)*(y1-yi)*(y1-yi)*(y1-yi) +
                 w2*(x2-xi)*(y2-yi)*(y2-yi)*(y2-yi) +
                 w3*(x3-xi)*(y3-yi)*(y3-yi)*(y3-yi);

	  Iyyy = Iyyy * nx * len/2.;

	  // (x-x_i)^2*(y-y_i) -> 1/3 (x-x_i)^3(y-y_i)
	  Ixxy = w1*(x1-xi)*(x1-xi)*(x1-xi)*(y1-yi) +
                 w2*(x2-xi)*(x2-xi)*(x2-xi)*(y2-yi) +
                 w3*(x3-xi)*(x3-xi)*(x3-xi)*(y3-yi);

	  Ixxy = (1./3.) * Ixxy * nx * len/2.;

	  // (x-x_i)*(y-y_i)^2 -> 1/2 (x-x_i)^2(y-y_i)^2
	  Ixyy = w1*(x1-xi)*(x1-xi)*(y1-yi)*(y1-yi) +
                 w2*(x2-xi)*(x2-xi)*(y2-yi)*(y2-yi) +
                 w3*(x3-xi)*(x3-xi)*(y3-yi)*(y3-yi);

	  Ixyy = 0.5 * Ixyy * nx * len/2.;

	  // Accumulate the moments to the node.
	  grid->Moments[node*NUM_MOM+0] += Ix;
	  grid->Moments[node*NUM_MOM+1] += Iy;
	  grid->Moments[node*NUM_MOM+2] += Ixx;
	  grid->Moments[node*NUM_MOM+3] += Iyy;
	  grid->Moments[node*NUM_MOM+4] += Ixy;
	  grid->Moments[node*NUM_MOM+5] += Ixxx;
	  grid->Moments[node*NUM_MOM+6] += Iyyy;
	  grid->Moments[node*NUM_MOM+7] += Ixxy;
	  grid->Moments[node*NUM_MOM+8] += Ixyy;
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
	  
	  // Now I can apply the function from each of the moment pieces to the Gauss points.

	  XI[0] = grid->x[node*NDIM+0];
	  XI[1] = grid->x[node*NDIM+1];

	  // Lets try to be more efficient about how we calculate them (or why I changed it from the doing each one in turn).
	  Ix = 0;
	  Iy = 0.;
	  Ixx = 0.;
	  Iyy = 0.;
	  Ixy = 0.;
	  Ixxx = 0.;
	  Iyyy = 0.;
	  Ixxy = 0.;
	  Ixyy = 0.;

	  for ( k=0; k < 3; k++ )
	    {
	      curved_boundary_normal_vector ( grid, bc, &(GP[k*NDIM]), &nx, &ny );

	      Ix +=   ( w[k] * (GP[k*NDIM+0]-XI[0]) * (GP[k*NDIM+0]-XI[0])  * nx * len );
	      Iy +=   ( w[k] * (GP[k*NDIM+0]-XI[0]) * (GP[k*NDIM+1]-XI[1])  * nx * len );
	      Ixx +=  ( w[k] * (GP[k*NDIM+0]-XI[0]) * (GP[k*NDIM+0]-XI[0]) * (GP[k*NDIM+0]-XI[0]) * nx * len );
	      Iyy +=  ( w[k] * (GP[k*NDIM+0]-XI[0]) * (GP[k*NDIM+1]-XI[1]) * (GP[k*NDIM+1]-XI[1]) * nx * len );
	      Ixy +=  ( w[k] * (GP[k*NDIM+0]-XI[0]) * (GP[k*NDIM+0]-XI[0]) * (GP[k*NDIM+1]-XI[1]) * nx * len );
	      Ixxx += ( w[k] * (GP[k*NDIM+0]-XI[0]) * (GP[k*NDIM+0]-XI[0]) * (GP[k*NDIM+0]-XI[0]) * (GP[k*NDIM+0]-XI[0]) * nx * len );
	      Iyyy += ( w[k] * (GP[k*NDIM+0]-XI[0]) * (GP[k*NDIM+1]-XI[1]) * (GP[k*NDIM+1]-XI[1]) * (GP[k*NDIM+1]-XI[1]) * nx * len );
	      Ixxy += ( w[k] * (GP[k*NDIM+0]-XI[0]) * (GP[k*NDIM+0]-XI[0]) * (GP[k*NDIM+0]-XI[0]) * (GP[k*NDIM+1]-XI[1]) * nx * len );
	      Ixyy += ( w[k] * (GP[k*NDIM+0]-XI[0]) * (GP[k*NDIM+0]-XI[0]) * (GP[k*NDIM+1]-XI[1]) * (GP[k*NDIM+1]-XI[1]) * nx * len );
	    }
	  
	  // Adjust the values to account for the constant and the length.
	  Ix = 0.5 * Ix * 0.5;
	  Iy = Iy * 0.5;
	  Ixx = (1./3.) * Ixx * 0.5;
	  Iyy = Iyy * 0.5;
	  Ixy = 0.5 * Ixy * 0.5;
	  Ixxx = (0.25) * Ixxx * 0.5;
	  Iyyy = Iyyy * 0.5;
	  Ixxy = (1./3.) * Ixxy * 0.5;
	  Ixyy = 0.5 * Ixyy * 0.5;

	  // Accumulate the moments to the node.
	  grid->Moments[node*NUM_MOM+0] += Ix;
	  grid->Moments[node*NUM_MOM+1] += Iy;
	  grid->Moments[node*NUM_MOM+2] += Ixx;
	  grid->Moments[node*NUM_MOM+3] += Iyy;
	  grid->Moments[node*NUM_MOM+4] += Ixy;
	  grid->Moments[node*NUM_MOM+5] += Ixxx;
	  grid->Moments[node*NUM_MOM+6] += Iyyy;
	  grid->Moments[node*NUM_MOM+7] += Ixxy;
	  grid->Moments[node*NUM_MOM+8] += Ixyy;
	  
	}

    }


  // Start the process by looping over all the elements and looping through the internal median dual pieces.
  for ( e=Tri; e <= Quad; e++ )                 // Element types loop.
    {
      // Find the number of nodes and edges for this element type.
      num_vert = NumberNodesForElement(e);
      //printf("  Processing element type %d\n",e);
      
      for ( j=1; j <= grid->num_elem[e]; j++ )              // Element loop.
	{
	  //printf("    Processing element number %d\n",j);
	  
	  // Retrieve the element node ids. This is for convenience.
	  for ( k=0; k < num_vert; k++ )              // Element node extraction loop.
	    {
	      nodes[k] = grid->c2n[e][j*num_vert+k];
	    }                                         // End element node extraction loop.
	      
	  // Find the element centroid.
	  XC[0] = 0.; XC[1] = 0.;
	  for ( k=0; k < num_vert; k++ )
	    {
	      for ( i=0; i < NDIM; i++ )
		{
		  XC[i] += grid->x[NDIM*nodes[k]+i];
		}
	    }
	  XC[0] /= ((double)num_vert);
	  XC[1] /= ((double)num_vert);
	  
	  // Loop through each node in the element.
	  for ( k=0; k < num_vert; k++ )              // Element node loop.
	    {
	      //printf("      Processing node %d\n",k);
	      
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
		  v_middle[i] = XC[i]    - grid->x[NDIM*node+i];
		}

	      // Process the first triangle.
	      //printf("        Processing triangle 1.\n");

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
		  ct_nodes[1][0] = XC[0];
		  ct_nodes[1][1] = XC[1];

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

		  zi = ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
			 (ct_nodes[5][0]*ct_nodes[3][1] - ct_nodes[5][1]*ct_nodes[3][0]) +
			 (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
		       ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
			 (ct_nodes[2][0]*ct_nodes[3][1] - ct_nodes[2][1]*ct_nodes[3][0]) +
			 (ct_nodes[2][0]*ct_nodes[1][1] - ct_nodes[2][1]*ct_nodes[1][0])  );

		  eta = ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			  (ct_nodes[5][0]*ct_nodes[2][1] - ct_nodes[5][1]*ct_nodes[2][0]) +
			  (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
		        ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			  (ct_nodes[3][0]*ct_nodes[2][1] - ct_nodes[3][1]*ct_nodes[2][0]) +
			  (ct_nodes[3][0]*ct_nodes[1][1] - ct_nodes[1][0]*ct_nodes[3][1])  );

		  if ( eta < 0.25 || zi < 0.25 )  // The transformation will not be invertible.
		    {
		      printf("CRITICAL ERROR: Triangle 1 of element %d of type %d is too concave (cv_calc.C)!\n",j,e);
		      printf("  zi = %.15e     eta = %.15e\n",zi,eta);
		      fflush(stdout);
		      exit(1);
		    }

		  // Do this for each variable.
		  for ( i=0; i < NUM_MOM; i++ )
		    {
		      QTRI[i] = 0.;
		    }
		  
		  // Now I can loop over the Gaussian integration nodes and apply the quadratic transformation.
		  for ( i=0; i < 7; i++ )
		    {
		      zi = ct_gp[i][0];
		      eta = ct_gp[i][1];
		      
		      // Compute basis functions.
		      N1 = 1. - 3.*zi - 3.*eta + 4.*zi*eta + 2.*zi*zi + 2.*eta*eta;
		      N2 = zi*(2.*zi - 1.);
		      N3 = eta*(2.*eta - 1.);
		      N4 = 4.*zi*(1. - zi - eta);
		      N5 = 4.*zi*eta;
		      N6 = 4.*eta*(1. - zi - eta);

		      // Compute their derivatives.
		      N1r = -3. + 4.*eta + 4.*zi;
		      N1s = -3. + 4.*zi + 4.*eta;
		      N2r = 4.*zi - 1.;
		      N2s = 0.;
		      N3r = 0.;
		      N3s = 4.*eta - 1.;
		      N4r = 4.*(1. - 2.*zi - eta);
		      N4s = -4.*zi;
		      N5r = 4.*eta;
		      N5s = 4.*zi;
		      N6r = -4.*eta;
		      N6s = 4.*(1. - 2.*eta - zi);

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

		      
		      Ix = ( XGP[0] - grid->x[NDIM*node+0]);
		      Iy = ( XGP[1] - grid->x[NDIM*node+1]);
		      Ixx = ( XGP[0] - grid->x[NDIM*node+0])*( XGP[0] - grid->x[NDIM*node+0]);
		      Iyy = ( XGP[1] - grid->x[NDIM*node+1])*( XGP[1] - grid->x[NDIM*node+1]);
		      Ixy = ( XGP[0] - grid->x[NDIM*node+0])*( XGP[1] - grid->x[NDIM*node+1]);
		      Ixxx = ( XGP[0] - grid->x[NDIM*node+0])*( XGP[0] - grid->x[NDIM*node+0])*( XGP[0] - grid->x[NDIM*node+0]);
		      Iyyy = ( XGP[1] - grid->x[NDIM*node+1])*( XGP[1] - grid->x[NDIM*node+1])*( XGP[1] - grid->x[NDIM*node+1]);
		      Ixxy = ( XGP[0] - grid->x[NDIM*node+0])*( XGP[0] - grid->x[NDIM*node+0])*( XGP[1] - grid->x[NDIM*node+1]);
		      Ixyy = ( XGP[0] - grid->x[NDIM*node+0])*( XGP[1] - grid->x[NDIM*node+1])*( XGP[1] - grid->x[NDIM*node+1]);
		      
		      QTRI[0] += ( wg[i] * Ix * jacobian * 0.5 );
		      QTRI[1] += ( wg[i] * Iy * jacobian * 0.5 );
		      QTRI[2] += ( wg[i] * Ixx * jacobian * 0.5 );
		      QTRI[3] += ( wg[i] * Iyy * jacobian * 0.5 );
		      QTRI[4] += ( wg[i] * Ixy * jacobian * 0.5 );
		      QTRI[5] += ( wg[i] * Ixxx * jacobian * 0.5 );
		      QTRI[6] += ( wg[i] * Iyyy * jacobian * 0.5 );
		      QTRI[7] += ( wg[i] * Ixxy * jacobian * 0.5 );
		      QTRI[8] += ( wg[i] * Ixyy * jacobian * 0.5 );
		    }
		  
		  // Add to the running total for the node.
		  for ( i=0; i < NUM_MOM; i++ )
		    {
		      moments_area_integral[node*NUM_MOM+i] += QTRI[i];
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
					bary_coord[v][2] * XC[ind]);
			}
		      gp[v][2] = 0.;   // Always 0 since 2D vectors.
		    }

		  // Now we loop over the Gauss points and get the function value
		  // and perform the numerical integration (the area integral of the triangle).
	      
		  for ( i=0; i < NUM_MOM; i++ )
		    QTRI[i] = 0.;
		  
		  for ( i=0; i < 7; i++ )
		    {
		      // Get the value at the point.
		      XGP[0] = gp[i][0];
		      XGP[1] = gp[i][1];
		      
		      Ix = ( XGP[0] - grid->x[NDIM*node+0]);
		      Iy = ( XGP[1] - grid->x[NDIM*node+1]);
		      Ixx = ( XGP[0] - grid->x[NDIM*node+0])*( XGP[0] - grid->x[NDIM*node+0]);
		      Iyy = ( XGP[1] - grid->x[NDIM*node+1])*( XGP[1] - grid->x[NDIM*node+1]);
		      Ixy = ( XGP[0] - grid->x[NDIM*node+0])*( XGP[1] - grid->x[NDIM*node+1]);
		      Ixxx = ( XGP[0] - grid->x[NDIM*node+0])*( XGP[0] - grid->x[NDIM*node+0])*( XGP[0] - grid->x[NDIM*node+0]);
		      Iyyy = ( XGP[1] - grid->x[NDIM*node+1])*( XGP[1] - grid->x[NDIM*node+1])*( XGP[1] - grid->x[NDIM*node+1]);
		      Ixxy = ( XGP[0] - grid->x[NDIM*node+0])*( XGP[0] - grid->x[NDIM*node+0])*( XGP[1] - grid->x[NDIM*node+1]);
		      Ixyy = ( XGP[0] - grid->x[NDIM*node+0])*( XGP[1] - grid->x[NDIM*node+1])*( XGP[1] - grid->x[NDIM*node+1]);
		      
		      QTRI[0] += ( wg[i] * Ix  );
		      QTRI[1] += ( wg[i] * Iy  );
		      QTRI[2] += ( wg[i] * Ixx );
		      QTRI[3] += ( wg[i] * Iyy );
		      QTRI[4] += ( wg[i] * Ixy );
		      QTRI[5] += ( wg[i] * Ixxx );
		      QTRI[6] += ( wg[i] * Iyyy );
		      QTRI[7] += ( wg[i] * Ixxy );
		      QTRI[8] += ( wg[i] * Ixyy );
		    }
		  
		  for ( i=0; i < NUM_MOM; i++ )
		    QTRI[i] *= dA;
		  
		  // Accumulate the integral to the node.
		  for ( i=0; i < NUM_MOM; i++ )
		    moments_area_integral[node*NUM_MOM+i] += QTRI[i];
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
		  ct_nodes[1][0] = XC[0];
		  ct_nodes[1][1] = XC[1];

		  // Node 2 is the ghost node.
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

		  zi = ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
			 (ct_nodes[5][0]*ct_nodes[3][1] - ct_nodes[5][1]*ct_nodes[3][0]) +
			 (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
		       ( (ct_nodes[1][0]*ct_nodes[3][1] - ct_nodes[1][1]*ct_nodes[3][0]) -
			 (ct_nodes[2][0]*ct_nodes[3][1] - ct_nodes[2][1]*ct_nodes[3][0]) +
			 (ct_nodes[2][0]*ct_nodes[1][1] - ct_nodes[2][1]*ct_nodes[1][0])  );

		  eta = ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			  (ct_nodes[5][0]*ct_nodes[2][1] - ct_nodes[5][1]*ct_nodes[2][0]) +
			  (ct_nodes[5][0]*ct_nodes[1][1] - ct_nodes[5][1]*ct_nodes[1][0])  ) /
		        ( (ct_nodes[1][0]*ct_nodes[2][1] - ct_nodes[1][1]*ct_nodes[2][0]) -
			  (ct_nodes[3][0]*ct_nodes[2][1] - ct_nodes[3][1]*ct_nodes[2][0]) +
			  (ct_nodes[3][0]*ct_nodes[1][1] - ct_nodes[1][0]*ct_nodes[3][1])  );

		  if ( eta < 0.25 || zi < 0.25 )  // The transformation will not be invertible.
		    {
		      printf("CRITICAL ERROR: Triangle 2 of element %d of type %d is too concave (cv_calc.C)!\n",j,e);
		      printf("  zi = %.15e     eta = %.15e\n",zi,eta);
		      fflush(stdout);
		      exit(1);
		    }

		  for ( i=0; i < NUM_MOM; i++ )
		    {
		      QTRI[i] = 0.;
		    }

		  // Now I can loop over the Gaussian integration nodes and apply the quadratic transformation.
		  for ( i=0; i < 7; i++ )
		    {
		      zi = ct_gp[i][0];
		      eta = ct_gp[i][1];
		      
		      // Compute basis functions.
		      N1 = 1. - 3.*zi - 3.*eta + 4.*zi*eta + 2.*zi*zi + 2.*eta*eta;
		      N2 = zi*(2.*zi - 1.);
		      N3 = eta*(2.*eta - 1.);
		      N4 = 4.*zi*(1. - zi - eta);
		      N5 = 4.*zi*eta;
		      N6 = 4.*eta*(1. - zi - eta);

		      // Compute their derivatives.
		      N1r = -3. + 4.*eta + 4.*zi;
		      N1s = -3. + 4.*zi + 4.*eta;
		      N2r = 4.*zi - 1.;
		      N2s = 0.;
		      N3r = 0.;
		      N3s = 4.*eta - 1.;
		      N4r = 4.*(1. - 2.*zi - eta);
		      N4s = -4.*zi;
		      N5r = 4.*eta;
		      N5s = 4.*zi;
		      N6r = -4.*eta;
		      N6s = 4.*(1. - 2.*eta - zi);

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

		      
		      Ix = ( XGP[0] - grid->x[NDIM*node+0]);
		      Iy = ( XGP[1] - grid->x[NDIM*node+1]);
		      Ixx = ( XGP[0] - grid->x[NDIM*node+0])*( XGP[0] - grid->x[NDIM*node+0]);
		      Iyy = ( XGP[1] - grid->x[NDIM*node+1])*( XGP[1] - grid->x[NDIM*node+1]);
		      Ixy = ( XGP[0] - grid->x[NDIM*node+0])*( XGP[1] - grid->x[NDIM*node+1]);
		      Ixxx = ( XGP[0] - grid->x[NDIM*node+0])*( XGP[0] - grid->x[NDIM*node+0])*( XGP[0] - grid->x[NDIM*node+0]);
		      Iyyy = ( XGP[1] - grid->x[NDIM*node+1])*( XGP[1] - grid->x[NDIM*node+1])*( XGP[1] - grid->x[NDIM*node+1]);
		      Ixxy = ( XGP[0] - grid->x[NDIM*node+0])*( XGP[0] - grid->x[NDIM*node+0])*( XGP[1] - grid->x[NDIM*node+1]);
		      Ixyy = ( XGP[0] - grid->x[NDIM*node+0])*( XGP[1] - grid->x[NDIM*node+1])*( XGP[1] - grid->x[NDIM*node+1]);
		      
		      QTRI[0] += ( wg[i] * Ix * jacobian * 0.5 );
		      QTRI[1] += ( wg[i] * Iy * jacobian * 0.5 );
		      QTRI[2] += ( wg[i] * Ixx * jacobian * 0.5 );
		      QTRI[3] += ( wg[i] * Iyy * jacobian * 0.5 );
		      QTRI[4] += ( wg[i] * Ixy * jacobian * 0.5 );
		      QTRI[5] += ( wg[i] * Ixxx * jacobian * 0.5 );
		      QTRI[6] += ( wg[i] * Iyyy * jacobian * 0.5 );
		      QTRI[7] += ( wg[i] * Ixxy * jacobian * 0.5 );
		      QTRI[8] += ( wg[i] * Ixyy * jacobian * 0.5 );
		    }
		  
		  // Add to the running total for the node.
		  for ( i=0; i < NUM_MOM; i++ )
		    {
		      moments_area_integral[node*NUM_MOM+i] += QTRI[i];
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
					bary_coord[v][2] * XC[ind]);
			}
		      gp[v][2] = 0.;   // Always 0 since 2D vectors.
		    }
	      
		  // Now we loop over the Gauss points and get the function value
		  // and perform the numerical integration (the area integral of the triangle).
		  
		  for ( i=0; i < NUM_MOM; i++ )
		    QTRI[i] = 0.;
	      
		  for ( i=0; i < 7; i++ )
		    {
		      // Get the value at the point.
		      XGP[0] = gp[i][0];
		      XGP[1] = gp[i][1];
		      
		      Ix = ( XGP[0] - grid->x[NDIM*node+0]);
		      Iy = ( XGP[1] - grid->x[NDIM*node+1]);
		      Ixx = ( XGP[0] - grid->x[NDIM*node+0])*( XGP[0] - grid->x[NDIM*node+0]);
		      Iyy = ( XGP[1] - grid->x[NDIM*node+1])*( XGP[1] - grid->x[NDIM*node+1]);
		      Ixy = ( XGP[0] - grid->x[NDIM*node+0])*( XGP[1] - grid->x[NDIM*node+1]);
		      Ixxx = ( XGP[0] - grid->x[NDIM*node+0])*( XGP[0] - grid->x[NDIM*node+0])*( XGP[0] - grid->x[NDIM*node+0]);
		      Iyyy = ( XGP[1] - grid->x[NDIM*node+1])*( XGP[1] - grid->x[NDIM*node+1])*( XGP[1] - grid->x[NDIM*node+1]);
		      Ixxy = ( XGP[0] - grid->x[NDIM*node+0])*( XGP[0] - grid->x[NDIM*node+0])*( XGP[1] - grid->x[NDIM*node+1]);
		      Ixyy = ( XGP[0] - grid->x[NDIM*node+0])*( XGP[1] - grid->x[NDIM*node+1])*( XGP[1] - grid->x[NDIM*node+1]);
		      
		      QTRI[0] += ( wg[i] * Ix  );
		      QTRI[1] += ( wg[i] * Iy  );
		      QTRI[2] += ( wg[i] * Ixx );
		      QTRI[3] += ( wg[i] * Iyy );
		      QTRI[4] += ( wg[i] * Ixy );
		      QTRI[5] += ( wg[i] * Ixxx );
		      QTRI[6] += ( wg[i] * Iyyy );
		      QTRI[7] += ( wg[i] * Ixxy );
		      QTRI[8] += ( wg[i] * Ixyy );
		    }
	      
		  for ( i=0; i < NUM_MOM; i++ )
		    QTRI[i] *= dA;
		  
		  // Accumulate the integral to the node.
		  for ( i=0; i < NUM_MOM; i++ )
		    moments_area_integral[node*NUM_MOM+i] += QTRI[i];
		}
	      
	    }                                       // End element noqde loop.
	  
	}                                           // End element loop.
      
    }                                               // End element type loop.

  // Finished the accumulation. Now divide by the control volume area.
  for ( i=1; i<= grid->nn; i++ )
    {
      for ( j=0; j < NUM_MOM; j++ )
	{
	  grid->Moments[i*NUM_MOM+j] /= grid->cv_area[i];
	  moments_area_integral[i*NUM_MOM+j] /= grid->cv_area[i];
	}
    }

  // Debug code to print to screen.
  #if 0
  printf("CV MOMENTS DEBUG:\n");
  for ( n=1; n <= grid->nn; n++ )
    {
      printf("%d :\n",n);
      printf("  Ix  = %lf\n",grid->Moments[n*NUM_MOM+0]);
      printf("  Iy  = %lf\n",grid->Moments[n*NUM_MOM+1]);
      printf("  Ixx = %lf\n",grid->Moments[n*NUM_MOM+2]);
      printf("  Iyy = %lf\n",grid->Moments[n*NUM_MOM+3]);
      printf("  Ixy = %lf\n",grid->Moments[n*NUM_MOM+4]);
      printf("  Ixxx = %lf\n",grid->Moments[n*NUM_MOM+5]);
      printf("  Iyyy = %lf\n",grid->Moments[n*NUM_MOM+6]);
      printf("  Ixxy = %lf\n",grid->Moments[n*NUM_MOM+7]);
      printf("  Ixyy = %lf\n",grid->Moments[n*NUM_MOM+8]);
    }
  #endif


  #if 0
  printf("Difference in CV MOMENTS between divergence and area integral:\n");
  for ( n=1; n <= grid->nn; n++ )
    {
      printf("%d :\n",n);
      printf("  Ix  = %.15e\n",fabs( grid->Moments[n*NUM_MOM+0] - moments_area_integral[n*NUM_MOM+0] ) );
      printf("  Iy  = %.15e\n",fabs( grid->Moments[n*NUM_MOM+1] - moments_area_integral[n*NUM_MOM+1] ) );
      printf("  Ixx = %.15e\n",fabs( grid->Moments[n*NUM_MOM+2] - moments_area_integral[n*NUM_MOM+2] ) );
      printf("  Iyy = %.15e\n",fabs( grid->Moments[n*NUM_MOM+3] - moments_area_integral[n*NUM_MOM+3] ) );
      printf("  Ixy = %.15e\n",fabs( grid->Moments[n*NUM_MOM+4] - moments_area_integral[n*NUM_MOM+4] ) );
      printf("  Ixxx = %.15e\n",fabs( grid->Moments[n*NUM_MOM+5] - moments_area_integral[n*NUM_MOM+5] ) );
      printf("  Iyyy = %.15e\n",fabs( grid->Moments[n*NUM_MOM+6] - moments_area_integral[n*NUM_MOM+6] ) );
      printf("  Ixxy = %.15e\n",fabs( grid->Moments[n*NUM_MOM+7] - moments_area_integral[n*NUM_MOM+7] ) );
      printf("  Ixyy = %.15e\n",fabs( grid->Moments[n*NUM_MOM+8] - moments_area_integral[n*NUM_MOM+8] ) );
    }
  #endif

  freenull(moments_area_integral);
  
  return;
}

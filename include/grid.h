//=============================================================
// 
//  grid.h
//  
//  A simple class that will contain all the data needed for a
//  grid object.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include "Defined_Variables.h"
 
#ifndef GRID_H
#define GRID_H

class GRID
{
  public:

  // Initialize everything to zero.
  GRID() 
    {
      nn = 0;
      nn_ghost = 0;
      nedges = 0;
      nbedges = 0;
      nsubedges = 0;
      iter = 1;
      rfiter = 0;
      nsteps = 0;
      citer = 0;
      cfl = 0.;
      eltime = 0.;
      for ( int i=0; i < NUM_ELEM_TYPE; i++ )
	num_elem[i] = 0;

      c2n = 0;
      c2e = 0;
      csn = 0;
      ncsn = 0;
      nsn = 0;
      nnsn = 0;
      nsn2 = 0;
      nnsn2 = 0;
      gnsn = 0;
      gnnsn = 0;
      gnnsn1 = 0;
      node_fill_level = 0;
      esn = 0;
      ese = 0;
      nese = 0;
      edges = 0;
      bedges = 0;
      nn_map = 0;
      node_order = 0;
      nb = 0;
      nbs = 0;
      bs = 0;
      bbc = 0;
      x = 0;
      cv_area = 0;
      nn_map = 0;
      Q = 0;
      nQ = 0;
      dQ = 0;
      pQ = 0;
      length_scale = 0;
      ev_con = 0;
      phi = 0;
      grad = 0;
      hess = 0;
      Moments = 0;
      el_cent = 0;
      subedges = 0;
      xn_edges = 0;
      xn_subedges = 0;
      xm_subedges = 0;
      xn_bedges = 0;
      xg_subedges = 0;
      xg_bedges = 0;
      xgn_bedges = 0;
      edge_weight = 0;
      point_value = 0;
      dt = 0;
      R = 0;
      Qfac = 0;
      Rfac = 0;
      Qe = 0;
      Source = 0;
      node_state = 0;
      SVD_U = 0;
      SVD_V = 0;
      SVD_S = 0;
      SVD_size = 0;
      SVDinv_order = 0;
     }
  
  ~GRID() 
    { 
      int b,s;             // Loop counters.
      
      // Free up all memory.

      if ( c2n != 0 )
	{
	  for (b=Edge; b <= Quad; b++)
	    {
	      freenull ( c2n[b] );
	    }

	  freenull ( c2n );
	}

      if ( c2e != 0 )
	{
	  for (b=Tri; b <= Quad; b++)
	    {
	      freenull ( c2e[b] );
	    }

	  freenull ( c2e );
	}

      freenull(x);
      freenull(csn);
      freenull(ncsn);
      freenull(nsn);
      freenull(nnsn);
      freenull(nsn2);
      freenull(nnsn2);
      freenull(gnsn);
      freenull(gnnsn);
      freenull(gnnsn1);
      freenull(node_fill_level);
      freenull(esn);
      freenull(ese);
      freenull(nese);
      freenull(edges);
      freenull(bedges);
      freenull(cv_area);
      freenull(nn_map);
      freenull(node_order);
      freenull(bbc);
      freenull(Q);
      freenull(nQ);
      freenull(dQ);
      freenull(pQ);
      freenull(length_scale);
      freenull(ev_con);
      freenull(dt);
      freenull(phi);
      freenull(grad);
      freenull(hess);
      freenull(Moments);
      freenull(subedges);
      freenull(xn_edges);
      freenull(xn_subedges);
      freenull(xm_subedges);
      freenull(xg_subedges);
      freenull(xg_bedges);
      freenull(xgn_bedges);
      freenull(el_cent);
      freenull(xn_bedges);
      freenull(edge_weight);
      freenull(point_value);
      freenull(R);
      freenull(Qe);
      freenull(Source);
      freenull(node_state);

      // Loop through bs array and free up memory in reverse order of allocation. This assumes that the mesh file format
      // of Dr. Karman is used.
      if ( bs != 0 )
	{
	  for ( b=1; b <= nb; b++ )
	    {
	      for ( s=1; s <= nbs[b]; s++ )
		{
		  freenull ( bs[b][s] );
		}
	      
	      freenull ( bs[b] );
	    }
      
	  freenull ( bs );
	}

      freenull(nbs);

      if ( Qfac != 0 )
	{
	  for ( b=1; b <= nn; b++ )
	    {
	      freenull(Qfac[b]);
	      freenull(Rfac[b]);
	    }
	}
      
      freenull(Qfac);
      freenull(Rfac);

      if ( SVD_U != 0 )
	{
	  for ( b=1; b <= nn; b++ )
	    {
	      freenull(SVD_U[b]);
	      freenull(SVD_V[b]);
	      freenull(SVD_S[b]);
	    }
	}
      
      freenull(SVD_U);
      freenull(SVD_V);
      freenull(SVD_S);
      freenull(SVD_size);
      
    }

  

  // Data members.

  int nn;                       // Number of nodes.
  int nn_ghost;                 // Number of ghost nodes.
  int *node_state;              // Describes the state of the node.
  int nedges;                   // Number of edges.
  int nsubedges;                // Number of subedges.
  int nbedges;                  // Number of boundary edges.
  int iter;                     // Number of newton iterations.
  int nsteps;                   // Number of time steps.
  int citer;                    // Current time step.
  int rfiter;                   // The iteration we are starting from.
  int num_elem[NUM_ELEM_TYPE];  // Number of elements of each type.
  int **c2n;                    // Cell to node array for each element type.
  int **c2e;                    // Cell to edge array for each element type.
  int *csn;                     // Cells surrounding node array.
  int *ncsn;                    // Number of cells surrounding a node - used to access the list of cells
                                //  in the csn array.
  int *nsn;                     // Nodes surrounding node array.
  int *nnsn;                    // Number of nodes surrounding a node - used to access the list of nodes
                                //  in the nsn array.
  int *nsn2;                    // Like the nsn array but includes 2nd neighbors.
  int *nnsn2;                   // Pointers into nsn2 array.
  int *gnsn;                    // Generic node surrounding node. Allow different levels of fill per node. Support arbitrary number of layers.
  int *gnnsn;                   // Pointers into gnsn array.
  int *gnnsn1;                  // Pointer array into gnnsn. Tells where the nodes start.
  int *node_fill_level;         // Desribes the number of layers of neighbors for the node.
  int *esn;                     // Edges surrounding node array. When complete, has the same structure as nsn. So, nnsn can be
                                //  be recycled to point into it also.
  int *ese;                     // Element to neighboring element connectivity map.
  int *nese;                    // Number of elements neighboring an element - pointer to the *ese structure.
  int *edges;                   // The edge structure list. Has the left and right nodes of each edge.
  int *subedges;                // The subedge structure. Has the global edge and global element associated with it.
  int *bedges;                  // The boundary edge structure. Has the associated node, the boundary its on, and the segment it belongs to.
  int *node_order;              // An array that specifies the reconstruction order to apply to the node.
  double *xn_edges;             // Holds the normal vector - nx,ny,area.
  double *xn_subedges;          // Holds the normal vector - nx,ny,area.
  double *xm_subedges;          // Holds the midpoint of the subedge.
  double *xn_bedges;            // Holds the normal vector - nx,ny,area.
  double *xg_subedges;          // Holds the coordinates of the Gaussian nodes for the subedge.
  double *xg_bedges;            // Holds the coordinates of the Gaussian nodes for the bedge.
  double *xgn_bedges;           // Holds the normal vector information for the Gaussian nodes along a boundary edge.
  int *nn_map;                  // The new numbering map that contains the new index for the nodes after (Cuthill-McKee)
                                //  reordering. Ex. - nn_map[i] -> index; 'i' is the old sequential index from when the file
                                //  was read and 'index' is the new index used in the rest of the code.
  int nb;                       // Number of boundaries.
  int *nbs;                     // Number of segments per boundary.
  int ***bs;                    // Boundary segment array - bs[boundary][segment][node].

  int *bbc;                     // Boundary boundary condition array. I know it sounds redundant. This is a hangover from my 2D code
                                // used in CFD 2.
                                //  0 - Farfield
                                //  1 - Inviscid Solid Wall

  double *x;                    // Node coordinate array.
  double *cv_area;              // Area of the control volumes.
  double *el_cent;              // Element centroids.
  double *Q;                    // The solution vector.
  double *nQ;                   // The solution vector at the next step ( in our newton iterations ).
  double *pQ;                   // The solution vector at time level n-1. For 2nd order time accuracy.
  double *dQ;                   // The solution vector update.
  double *Qe;                   // Exact CV averaged value for the manufactured solution.
  double *Source;               // Source terms for the manufactured solution.
  double *length_scale;         // The length scale for each control volume.
  double *ev_con;               // The sum over the control volume of lambda_max * len.
  double *phi;                  // The limiter for the control volumes.
  double *grad;                 // The gradient for the control volumes. This is from Green-Gauss or Least Squares.
  double *hess;                 // The hessian for the control volumes. I want to store the gradient also from my method.
  double *Moments;              // The control volume moments.
  double *edge_weight;          // The edge weights used in the least squares calculation for the gradient.
  double *point_value;          // The value of the dependent variables at the vertex.
  double cfl;                   // The current cfl number.
  double eltime;                // Elapsed time of the solution. ( Makes sense for time accurate flow )
  double *dt;                   // Time step to be applied to the nodes.
  double *R;                    // Residual vector.
  double **Qfac;                // Q^t from the QR factorization.
  double **Rfac;                // R from the QR factorization.
  double **SVD_U;               // Storage space for the Moore-Penrose pseudoinverse.
  double **SVD_V;
  double **SVD_S;
  double *SVD_size;             // Size of the pseudoinverse components.
  int SVDinv_order;             // The order at which the pseudoinverse was made.
  char gridname[120];           // The mesh input file.
};

#endif

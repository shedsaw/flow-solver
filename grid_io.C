//=============================================================
// 
//  grid_io.C
//  
//  Functions that write out the solution data and gnuplot file.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "grid_io.h"
#include "params.h"
#include "Defined_Variables.h"
#include "curved_boundaries.h"
#include "reconstruct.h"
#include "hessian.h"
#include "mesh_utils.h"

//=============================================================
// 
//  write_gnuplot()
//  
//  Writes out a gnuplot file.
//
//  char *file_name                   // The file name to write out.
//  GRID *grid                        // The grid.
//
//=============================================================

void write_gnuplot ( char *file_name, GRID *grid )
{
  int i;                                 // Loop counters.
  FILE *fp;                              // Pointer to the file.
  int n1,n2,n3,n4;                       // Node ids of elements.

  //Open the file, but kick out if the file could not be opened.
  if (( fp = fopen(file_name,"w")) == 0)
    {
      printf("\nCouldn't open <%s>\n",file_name);
      exit(0);
    }

  printf("\n//==================================================\n");
  printf("// Starting PLOT IO - WRITING GNUPLOT FILE <%s>\n",file_name);
  printf("//==================================================\n");


  // Write out triangles.
  for (i=1; i <= grid->num_elem[Tri]; i++)
    {
      n1 = grid->c2n[Tri][i*3+0];
      n2 = grid->c2n[Tri][i*3+1];
      n3 = grid->c2n[Tri][i*3+2];
      fprintf(fp,"%19.10e %19.10e 0.0\n",  grid->x[NDIM*n1], grid->x[NDIM*n1+1]);
      fprintf(fp,"%19.10e %19.10e 0.0\n",  grid->x[NDIM*n2], grid->x[NDIM*n2+1]);
      fprintf(fp,"%19.10e %19.10e 0.0\n",  grid->x[NDIM*n3], grid->x[NDIM*n3+1]);
      fprintf(fp,"%19.10e %19.10e 0.0\n\n",grid->x[NDIM*n1], grid->x[NDIM*n1+1]);
    }

  // Write out quads.
  for (i=1; i <= grid->num_elem[Quad]; i++)
    {
      n1 = grid->c2n[Quad][i*4+0];
      n2 = grid->c2n[Quad][i*4+1];
      n3 = grid->c2n[Quad][i*4+2];
      n4 = grid->c2n[Quad][i*4+3];
      fprintf(fp,"%19.10e %19.10e 0.0\n",  grid->x[NDIM*n1], grid->x[NDIM*n1+1]);
      fprintf(fp,"%19.10e %19.10e 0.0\n",  grid->x[NDIM*n2], grid->x[NDIM*n2+1]);
      fprintf(fp,"%19.10e %19.10e 0.0\n",  grid->x[NDIM*n3], grid->x[NDIM*n3+1]);
      fprintf(fp,"%19.10e %19.10e 0.0\n",  grid->x[NDIM*n4], grid->x[NDIM*n4+1]);
      fprintf(fp,"%19.10e %19.10e 0.0\n\n",grid->x[NDIM*n1], grid->x[NDIM*n1+1]);
    }

  fclose(fp);

  printf("\n//==================================================\n");
  printf("// Finished PLOT IO.\n");
  printf("//==================================================\n\n");
  
  return;
}


//=============================================================
// 
//  write_tecplot_solutionC()
//  
//  Writes out a mesh with solution in the Tecplot file format.
//  Variables written - rho,u momentum, y momentum, E.
//
//  CONSERVED VARIABLES.
//
//  char *file_name                   // The file name.
//  GRID *grid                        // The grid.
//
//=============================================================

void write_tecplot_solutionC ( char *file_name, GRID *grid )
{
  int i;                                 // Loop counter.
  int n0, n1, n2, n3;                    // Element nodes.
  FILE *fp;                              // Pointer to a file.
  double rho;                            // Fluid denstiy.
  double u;                              // X-component of momentum.
  double v;                              // Y-component of momentum.
  double E;                              // energy.

  // Open the solution file for writing. Quit if it can't be opened.
  if ((fp = fopen ( file_name ,"w")) == 0 )
    {
      printf("\nCouldn't open <%s>\n",file_name);
      exit(0);
    }

  printf("//===================================================================\n");
  printf("// Starting SOLUTION IO - CONSERVED VARS - WRITING TECPLOT SOLUTION FILE <%s>\n",file_name);
  printf("//===================================================================\n");


  // Print the standard information to the top of the file. This must be here to work properly.
  fprintf(fp,"TITLE=\"data file\"\n");
  fprintf(fp,"VARIABLES=\"X\", \"Y\", \"Rho\", \"Umom\", \"Vmom\", \"E\"\n");

  // If there are no quads in the mesh, then write out the data as triangles. If there are quads, we treat tri's as
  // degenerate quads and write out only quads.
  if ( grid->num_elem[Quad] == 0 )
    {
      fprintf(fp,"ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",grid->nn,grid->num_elem[Tri]);
  
      // Loop over all the nodes and write out the x,y locations along with the desired variables.
      for ( i=1; i <= grid->nn; i++ )
	{
	  // Grab the data
	  rho = grid->Q[i*NUM_VAR+0];
	  u = grid->Q[i*NUM_VAR+1];
	  v = grid->Q[i*NUM_VAR+2];
	  E = grid->Q[i*NUM_VAR+3];
	  fprintf(fp,"%.15e %.15e %.15e %.15e %.15e %.15e\n",grid->x[i*2], grid->x[i*2+1],
		  rho,u,v,E);
	}

      fprintf(fp,"\n");
  
      // Print out the triangle connectivity to the file.
      for ( i=1; i <= grid->num_elem[Tri]; i++ )
	{
	  n0 = grid->c2n[Tri][i*3];
	  n1 = grid->c2n[Tri][i*3+1];
	  n2 = grid->c2n[Tri][i*3+2];
	  fprintf(fp,"%d %d %d\n",n0,n1,n2);
	}
      
      fclose(fp);
    }
  else
    {
      fprintf(fp,"ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n",grid->nn,
	      (grid->num_elem[Tri])+(grid->num_elem[Quad]));
  
      // Loop over all the nodes and write out the x,y locations along with the desired variables.
      for ( i=1; i <= grid->nn; i++ )
	{
	  // Grab the data
	  rho = grid->Q[i*NUM_VAR+0];
	  u = grid->Q[i*NUM_VAR+1];
	  v = grid->Q[i*NUM_VAR+2];
	  E = grid->Q[i*NUM_VAR+3];
	  fprintf(fp,"%.15e %.15e %.15e %.15e %.15e %.15e\n",grid->x[i*2], grid->x[i*2+1],
		  rho,u,v,E);
	}
      
      fprintf(fp,"\n");
  
      // Print out element connectivity to the file. For triangles, print the last node twice.
      for ( i=1; i <= grid->num_elem[Tri]; i++ )
	{
	  n0 = grid->c2n[Tri][i*3];
	  n1 = grid->c2n[Tri][i*3+1];
	  n2 = grid->c2n[Tri][i*3+2];
	  fprintf(fp,"%d %d %d %d\n",n0,n1,n2,n2);
	}
      for ( i=1; i <= grid->num_elem[Quad]; i++ )
	{
	  n0 = grid->c2n[Quad][i*4];
	  n1 = grid->c2n[Quad][i*4+1];
	  n2 = grid->c2n[Quad][i*4+2];
	  n3 = grid->c2n[Quad][i*4+3];
	  fprintf(fp,"%d %d %d %d\n",n0,n1,n2,n3);
	}
      
      fclose(fp);
    }
  return;
}


//=============================================================
// 
//  write_tecplot_solutionP()
//  
//  Writes out a mesh with solution in the Tecplot file format.
//  Variables written - rho,u,v,P.
//
//  PRIMITIVE VARIABLES.
//
//  char *file_name                   // The file name.
//  GRID *grid                        // The grid.
//  double gamma                      // Ratio of Specific heats.
//
//=============================================================

void write_tecplot_solutionP ( char *file_name, GRID *grid, double gamma )
{
  int i;                                 // Loop counter.
  int n0, n1, n2, n3;                    // Element nodes.
  FILE *fp;                              // Pointer to a file.
  double rho;                            // Fluid denstiy.
  double u;                              // X-component of velocity.
  double v;                              // Y-component of velocity.
  double P;                              // Pressure.
  double E;                              // Energy.

  // Open the solution file for writing. Quit if it can't be opened.
  if ((fp = fopen ( file_name ,"w")) == 0 )
    {
      printf("\nCouldn't open <%s>\n",file_name);
      exit(0);
    }

  printf("//===================================================================\n");
  printf("// Starting SOLUTION IO - CONSERVED VARS - WRITING TECPLOT SOLUTION FILE <%s>\n",file_name);
  printf("//===================================================================\n");


  // Print the standard information to the top of the file. This must be here to work properly.
  fprintf(fp,"TITLE=\"data file\"\n");
  fprintf(fp,"VARIABLES=\"X\", \"Y\", \"Rho\", \"U\", \"V\", \"P\"\n");

  // If there are no quads in the mesh, then write out the data as triangles. If there are quads, we treat tri's as
  // degenerate quads and write out only quads.
  if ( grid->num_elem[Quad] == 0 )
    {
      fprintf(fp,"ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",grid->nn,grid->num_elem[Tri]);
  
      // Loop over all the nodes and write out the x,y locations along with the desired variables.
      for ( i=1; i <= grid->nn; i++ )
	{
	  // Grab the data
	  rho = grid->Q[i*NUM_VAR+0];
	  u = (grid->Q[i*NUM_VAR+1])/rho;
	  v = (grid->Q[i*NUM_VAR+2])/rho;
	  E = (grid->Q[i*NUM_VAR+3]);

	  //P = (gamma - 1.)*E - ((gamma-1.)/(2.*rho))*(u*u + v*v);
	  P = (gamma - 1.)*E - ((gamma-1.)/2.)*rho*(u*u + v*v);
	  fprintf(fp,"%.15e %.15e %.15e %.15e %.15e %.15e\n",grid->x[i*2], grid->x[i*2+1],
		  rho,u,v,P);
	}

      fprintf(fp,"\n");
  
      // Print out the triangle connectivity to the file.
      for ( i=1; i <= grid->num_elem[Tri]; i++ )
	{
	  n0 = grid->c2n[Tri][i*3];
	  n1 = grid->c2n[Tri][i*3+1];
	  n2 = grid->c2n[Tri][i*3+2];
	  fprintf(fp,"%d %d %d\n",n0,n1,n2);
	}
      
      fclose(fp);
    }
  else
    {
      fprintf(fp,"ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n",grid->nn,
	      (grid->num_elem[Tri])+(grid->num_elem[Quad]));
  
      // Loop over all the nodes and write out the x,y locations along with the desired variables.
      for ( i=1; i <= grid->nn; i++ )
	{
	  // Grab the data
	  rho = grid->Q[i*NUM_VAR+0];
	  u = (grid->Q[i*NUM_VAR+1])/rho;
	  v = (grid->Q[i*NUM_VAR+2])/rho;
	  E = (grid->Q[i*NUM_VAR+3]);

	  P = (gamma - 1.)*E - ((gamma-1.)/(2.*rho))*(u*u + v*v);
	  fprintf(fp,"%.15e %.15e %.15e %.15e %.15e %.15e\n",grid->x[i*2], grid->x[i*2+1],
		  rho,u,v,P);
	}
      
      fprintf(fp,"\n");
  
      // Print out element connectivity to the file. For triangles, print the last node twice.
      for ( i=1; i <= grid->num_elem[Tri]; i++ )
	{
	  n0 = grid->c2n[Tri][i*3];
	  n1 = grid->c2n[Tri][i*3+1];
	  n2 = grid->c2n[Tri][i*3+2];
	  fprintf(fp,"%d %d %d %d\n",n0,n1,n2,n2);
	}
      for ( i=1; i <= grid->num_elem[Quad]; i++ )
	{
	  n0 = grid->c2n[Quad][i*4];
	  n1 = grid->c2n[Quad][i*4+1];
	  n2 = grid->c2n[Quad][i*4+2];
	  n3 = grid->c2n[Quad][i*4+3];
	  fprintf(fp,"%d %d %d %d\n",n0,n1,n2,n3);
	}
      
      fclose(fp);
    }
  return;
}

//=============================================================
// 
//  write_tecplot_solutionE()
//  
//  Writes out a mesh with solution in the Tecplot file format.
//
//  Q exact (MMS) - either primitive or conserved depending
//    on the initialization.
//
//  char *file_name                   // The file name.
//  GRID *grid                        // The grid.
//
//=============================================================

void write_tecplot_solutionE ( char *file_name, GRID *grid )
{
  int i;                                 // Loop counter.
  int n0, n1, n2, n3;                    // Element nodes.
  FILE *fp;                              // Pointer to a file.
  double rho;                            // Fluid denstiy.
  double u;                              // X-component of momentum.
  double v;                              // Y-component of momentum.
  double E;                              // energy.

  // Open the solution file for writing. Quit if it can't be opened.
  if ((fp = fopen ( file_name ,"w")) == 0 )
    {
      printf("\nCouldn't open <%s>\n",file_name);
      exit(0);
    }

  printf("//===================================================================\n");
  printf("// Starting SOLUTION IO - WRITING TECPLOT SOLUTION FILE <%s>\n",file_name);
  printf("//===================================================================\n");


  // Print the standard information to the top of the file. This must be here to work properly.
  fprintf(fp,"TITLE=\"data file\"\n");
  if ( RECON_PRIM )
    fprintf(fp,"VARIABLES=\"X\", \"Y\", \"Rho\", \"U\", \"V\", \"P\"\n");
  else
    fprintf(fp,"VARIABLES=\"X\", \"Y\", \"Rho\", \"Umom\", \"Vmom\", \"E\"\n");

  // If there are no quads in the mesh, then write out the data as triangles. If there are quads, we treat tri's as
  // degenerate quads and write out only quads.
  if ( grid->num_elem[Quad] == 0 )
    {
      fprintf(fp,"ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",grid->nn,grid->num_elem[Tri]);
  
      // Loop over all the nodes and write out the x,y locations along with the desired variables.
      for ( i=1; i <= grid->nn; i++ )
	{
	  // Grab the data
	  rho = grid->Qe[i*NUM_VAR+0];
	  u = grid->Qe[i*NUM_VAR+1];
	  v = grid->Qe[i*NUM_VAR+2];
	  E = grid->Qe[i*NUM_VAR+3];
	  fprintf(fp,"%.15e %.15e %.15e %.15e %.15e %.15e\n",grid->x[i*2], grid->x[i*2+1],
		  rho,u,v,E);
	}

      fprintf(fp,"\n");
  
      // Print out the triangle connectivity to the file.
      for ( i=1; i <= grid->num_elem[Tri]; i++ )
	{
	  n0 = grid->c2n[Tri][i*3];
	  n1 = grid->c2n[Tri][i*3+1];
	  n2 = grid->c2n[Tri][i*3+2];
	  fprintf(fp,"%d %d %d\n",n0,n1,n2);
	}
      
      fclose(fp);
    }
  else
    {
      fprintf(fp,"ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n",grid->nn,
	      (grid->num_elem[Tri])+(grid->num_elem[Quad]));
  
      // Loop over all the nodes and write out the x,y locations along with the desired variables.
      for ( i=1; i <= grid->nn; i++ )
	{
	  // Grab the data
	  rho = grid->Qe[i*NUM_VAR+0];
	  u = grid->Qe[i*NUM_VAR+1];
	  v = grid->Qe[i*NUM_VAR+2];
	  E = grid->Qe[i*NUM_VAR+3];
	  fprintf(fp,"%.15e %.15e %.15e %.15e %.15e %.15e\n",grid->x[i*2], grid->x[i*2+1],
		  rho,u,v,E);
	}
      
      fprintf(fp,"\n");
  
      // Print out element connectivity to the file. For triangles, print the last node twice.
      for ( i=1; i <= grid->num_elem[Tri]; i++ )
	{
	  n0 = grid->c2n[Tri][i*3];
	  n1 = grid->c2n[Tri][i*3+1];
	  n2 = grid->c2n[Tri][i*3+2];
	  fprintf(fp,"%d %d %d %d\n",n0,n1,n2,n2);
	}
      for ( i=1; i <= grid->num_elem[Quad]; i++ )
	{
	  n0 = grid->c2n[Quad][i*4];
	  n1 = grid->c2n[Quad][i*4+1];
	  n2 = grid->c2n[Quad][i*4+2];
	  n3 = grid->c2n[Quad][i*4+3];
	  fprintf(fp,"%d %d %d %d\n",n0,n1,n2,n3);
	}
      
      fclose(fp);
    }
  return;
}


//=============================================================
// 
//  write_tecplot_cv_area()
//  
//  Writes out a 2D mesh with only control volume area as
//  the variable.
//
//  char *file_name                   // The file name.
//  GRID *grid                        // The grid.
//
//=============================================================

void write_tecplot_cv_area ( char *file_name, GRID *grid )
{
  int i;                                 // Loop counter.
  int n0, n1, n2, n3;                    // Element nodes.
  FILE *fp;                              // Pointer to a file.

  // Open the solution file for writing. Quit if it can't be opened.
  if ((fp = fopen ( file_name ,"w")) == 0 )
    {
      printf("\nCouldn't open <%s>\n",file_name);
      exit(0);
    }

  printf("//===================================================================\n");
  printf("// Starting SOLUTION IO - WRITING TECPLOT SOLUTION FILE <%s>\n",file_name);
  printf("//===================================================================\n");

  // Print the standard information to the top of the file. This must be here to work properly.
  fprintf(fp,"TITLE=\"data file\"\n");
  fprintf(fp,"VARIABLES=\"X\", \"Y\", \"CV Area\"\n");
  
  // If there are no quads in the mesh, then write out the data as triangles. If there are quads, we treat tri's as
  // degenerate quads and write out only quads.
  if ( grid->num_elem[Quad] == 0 )
    {
      fprintf(fp,"ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",grid->nn,grid->num_elem[Tri]);
  
      // Loop over all the nodes and write out the x,y locations along with the desired variables.
      for ( i=1; i <= grid->nn; i++ )
	{
	  fprintf(fp,"%.15e %.15e %.15e\n",grid->x[i*2], grid->x[i*2+1], grid->cv_area[i]);
	}

      fprintf(fp,"\n");
  
      // Print out the triangle connectivity to the file.
      for ( i=1; i <= grid->num_elem[Tri]; i++ )
	{
	  n0 = grid->c2n[Tri][i*3];
	  n1 = grid->c2n[Tri][i*3+1];
	  n2 = grid->c2n[Tri][i*3+2];
	  fprintf(fp,"%d %d %d\n",n0,n1,n2);
	}
      
      fclose(fp);
      
      return;
    }
  else
    {
      fprintf(fp,"ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n",grid->nn,
	      (grid->num_elem[Tri])+(grid->num_elem[Quad]));
  
      // Loop over all the nodes and write out the x,y locations along with the desired variables.
      for ( i=1; i <= grid->nn; i++ )
	{
	  fprintf(fp,"%.15e %.15e %.15e\n",grid->x[i*2], grid->x[i*2+1], grid->cv_area[i]);
	}
      
      fprintf(fp,"\n");
  
      // Print out element connectivity to the file. For triangles, print the last node twice.
      for ( i=1; i <= grid->num_elem[Tri]; i++ )
	{
	  n0 = grid->c2n[Tri][i*3];
	  n1 = grid->c2n[Tri][i*3+1];
	  n2 = grid->c2n[Tri][i*3+2];
	  fprintf(fp,"%d %d %d %d\n",n0,n1,n2,n2);
	}
      for ( i=1; i <= grid->num_elem[Quad]; i++ )
	{
	  n0 = grid->c2n[Quad][i*4];
	  n1 = grid->c2n[Quad][i*4+1];
	  n2 = grid->c2n[Quad][i*4+2];
	  n3 = grid->c2n[Quad][i*4+3];
	  fprintf(fp,"%d %d %d %d\n",n0,n1,n2,n3);
	}
      
      fclose(fp);
      
      return;
    }

}

//=============================================================
// 
//  write_tecplot_derivs()
//  
//  Writes out a 2D mesh with the derivatives for every vertex.
//
//  char *file_name                   // The file name.
//  GRID *grid                        // The grid.
//
//=============================================================

void write_tecplot_derivs ( char *file_name, GRID *grid )
{
  int i;                                 // Loop counter.
  int n0, n1, n2, n3;                    // Element nodes.
  double derivs[NUM_MOM];                // Analytical values.
  FILE *fp;                              // Pointer to a file.

  // Open the solution file for writing. Quit if it can't be opened.
  if ((fp = fopen ( file_name ,"w")) == 0 )
    {
      printf("\nCouldn't open <%s>\n",file_name);
      exit(0);
    }
  
  printf("//===================================================================\n");
  printf("// Starting SOLUTION IO - WRITING TECPLOT DERIVATIVE FILE <%s>\n",file_name);
  printf("//===================================================================\n");
  
  // Print the standard information to the top of the file. This must be here to work properly.
  fprintf(fp,"TITLE=\"data file\"\n");
  if ( DEBUG_DERIVS )
    {
      fprintf(fp,"VARIABLES=\"X\", \"Y\", \"Ux\", \"Uy\", \"Uxx\", \"Uyy\", \"Uxy\", \"Uxxx\", \"Uyyy\", \"Uxxy\", \"Uxyy\", \"dUx\", \"dUy\", \"dUxx\", \"dUyy\", \"dUxy\", \"dUxxx\", \"dUyyy\", \"dUxxy\", \"dUxyy\"\n");
    }
  else if ( DEBUG_HESS )
    {
      fprintf(fp,"VARIABLES=\"X\", \"Y\", \"Ux\", \"Uy\", \"Uxx\", \"Uyy\", \"Uxy\", \"dUx\", \"dUy\", \"dUxx\", \"dUyy\", \"dUxy\"\n");
    }
  else if ( DEBUG_GRAD )
    {
      fprintf(fp,"VARIABLES=\"X\", \"Y\", \"Ux\", \"Uy\", \"dUx\", \"dUy\"\n");
    }
  else
    fprintf(fp,"VARIABLES=\"X\", \"Y\"\n");
  
  // If there are no quads in the mesh, then write out the data as triangles. If there are quads, we treat tri's as
  // degenerate quads and write out only quads.
  if ( grid->num_elem[Quad] == 0 )
    {
      fprintf(fp,"ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",grid->nn,grid->num_elem[Tri]);
  
      // Loop over all the nodes and write out the x,y locations along with the desired variables.
      for ( i=1; i <= grid->nn; i++ )
	{
	  Test_Function_Derivatives(&(grid->x[i*NDIM]),derivs);

	  if ( DEBUG_DERIVS )
	    {
	      fprintf(fp,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
		      grid->x[i*2],grid->x[i*2+1],
		      grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 0],grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 1],
		      grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 2],grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 3],
		      grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 4],grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 5],
		      grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 6],grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 7],
		      grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 8],
		      fabs( derivs[0] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 0] ),
		      fabs( derivs[1] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 1] ),
		      fabs( derivs[2] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 2] ),
		      fabs( derivs[3] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 3] ),
		      fabs( derivs[4] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 4] ),
		      fabs( derivs[5] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 5] ),
		      fabs( derivs[6] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 6] ),
		      fabs( derivs[7] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 7] ),
		      fabs( derivs[8] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 8] ) );
	    }
	  else if ( DEBUG_HESS )
	    {
	      fprintf(fp,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
		      grid->x[i*2],grid->x[i*2+1],
		      grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 0],grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 1],
		      grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 2],grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 3],
		      grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 4],
		      fabs( derivs[0] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 0] ),
		      fabs( derivs[1] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 1] ),
		      fabs( derivs[2] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 2] ),
		      fabs( derivs[3] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 3] ),
		      fabs( derivs[4] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 4] ) );
	    }
	  else if ( DEBUG_GRAD )
	    {
	      fprintf(fp,"%.15e %.15e %.15e %.15e %.15e %.15e\n",grid->x[i*2],grid->x[i*2+1],
		      grid->grad[i*NDIM*NUM_VAR + 0*NDIM + 0],grid->grad[i*NDIM*NUM_VAR + 0*NDIM + 1],
		      fabs( derivs[0] - grid->grad[i*NDIM*NUM_VAR + 0*NDIM + 0] ),
		      fabs( derivs[1] - grid->grad[i*NDIM*NUM_VAR + 0*NDIM + 1] ) );
	    }
	  else
	    fprintf(fp,"%.15e %.15e\n",grid->x[i*2], grid->x[i*2+1]);
	}

      fprintf(fp,"\n");
  
      // Print out the triangle connectivity to the file.
      for ( i=1; i <= grid->num_elem[Tri]; i++ )
	{
	  n0 = grid->c2n[Tri][i*3];
	  n1 = grid->c2n[Tri][i*3+1];
	  n2 = grid->c2n[Tri][i*3+2];
	  fprintf(fp,"%d %d %d\n",n0,n1,n2);
	}
      
      fclose(fp);
      
      return;
    }
  else
    {
      fprintf(fp,"ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n",grid->nn,
	      (grid->num_elem[Tri])+(grid->num_elem[Quad]));
  
      // Loop over all the nodes and write out the x,y locations along with the desired variables.
      for ( i=1; i <= grid->nn; i++ )
	{
	  Test_Function_Derivatives(&(grid->x[i*NDIM]),derivs);

	  if ( DEBUG_DERIVS )
	    {
	      fprintf(fp,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
		      grid->x[i*2],grid->x[i*2+1],
		      grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 0],grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 1],
		      grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 2],grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 3],
		      grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 4],grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 5],
		      grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 6],grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 7],
		      grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 8],
		      fabs( derivs[0] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 0] ),
		      fabs( derivs[1] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 1] ),
		      fabs( derivs[2] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 2] ),
		      fabs( derivs[3] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 3] ),
		      fabs( derivs[4] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 4] ),
		      fabs( derivs[5] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 5] ),
		      fabs( derivs[6] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 6] ),
		      fabs( derivs[7] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 7] ),
		      fabs( derivs[8] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 8] ) );
	    }
	  else if ( DEBUG_HESS )
	    {
	      fprintf(fp,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
		      grid->x[i*2],grid->x[i*2+1],
		      grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 0],grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 1],
		      grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 2],grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 3],
		      grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 4],
		      fabs( derivs[0] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 0] ),
		      fabs( derivs[1] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 1] ),
		      fabs( derivs[2] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 2] ),
		      fabs( derivs[3] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 3] ),
		      fabs( derivs[4] - grid->hess[i*NUM_MOM*NUM_VAR + 0*NUM_MOM + 4] ) );
	    }
	  else if ( DEBUG_GRAD )
	    {
	      fprintf(fp,"%.15e %.15e %.15e %.15e %.15e %.15e\n",grid->x[i*2],grid->x[i*2+1],
		      grid->grad[i*NDIM*NUM_VAR + 0*NDIM + 0],grid->grad[i*NDIM*NUM_VAR + 0*NDIM + 1],
		      fabs( derivs[0] - grid->grad[i*NDIM*NUM_VAR + 0*NDIM + 0] ),
		      fabs( derivs[1] - grid->grad[i*NDIM*NUM_VAR + 0*NDIM + 1] ) );
	    }
	  else
	    fprintf(fp,"%.15e %.15e\n",grid->x[i*2], grid->x[i*2+1]);
	}
      
      fprintf(fp,"\n");
  
      // Print out element connectivity to the file. For triangles, print the last node twice.
      for ( i=1; i <= grid->num_elem[Tri]; i++ )
	{
	  n0 = grid->c2n[Tri][i*3];
	  n1 = grid->c2n[Tri][i*3+1];
	  n2 = grid->c2n[Tri][i*3+2];
	  fprintf(fp,"%d %d %d %d\n",n0,n1,n2,n2);
	}
      for ( i=1; i <= grid->num_elem[Quad]; i++ )
	{
	  n0 = grid->c2n[Quad][i*4];
	  n1 = grid->c2n[Quad][i*4+1];
	  n2 = grid->c2n[Quad][i*4+2];
	  n3 = grid->c2n[Quad][i*4+3];
	  fprintf(fp,"%d %d %d %d\n",n0,n1,n2,n3);
	}
      
      fclose(fp);
      
      return;
    }

}

//=============================================================
// 
//  write_tecplot_cv_limiter()
//  
//  Writes out a 2D mesh with the limiters as the variables.
//
//  char *file_name                   // The file name.
//  GRID *grid                        // The grid.
//
//=============================================================

void write_tecplot_cv_limiter ( char *file_name, GRID *grid )
{
  int i;                                 // Loop counter.
  int n0, n1, n2, n3;                    // Element nodes.
  FILE *fp;                              // Pointer to a file.

  // Open the solution file for writing. Quit if it can't be opened.
  if ((fp = fopen ( file_name ,"w")) == 0 )
    {
      printf("\nCouldn't open <%s>\n",file_name);
      exit(0);
    }

  printf("//===================================================================\n");
  printf("// Starting SOLUTION IO - WRITING TECPLOT SOLUTION FILE <%s>\n",file_name);
  printf("//===================================================================\n");

  // Print the standard information to the top of the file. This must be here to work properly.
  fprintf(fp,"TITLE=\"data file\"\n");
  //fprintf(fp,"VARIABLES=\"X\", \"Y\", \"PHI\", \"SIGMA\"\n");
  fprintf(fp,"VARIABLES=\"X\", \"Y\", \"PHI1\", \"PHI2\", \"PHI3\", \"PHI4\", \"SIGMA1\", \"SIGMA2\", \"SIGMA3\", \"SIGMA4\"\n");
  
  // If there are no quads in the mesh, then write out the data as triangles. If there are quads, we treat tri's as
  // degenerate quads and write out only quads.
  if ( grid->num_elem[Quad] == 0 )
    {
      fprintf(fp,"ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",grid->nn,grid->num_elem[Tri]);
  
      // Loop over all the nodes and write out the x,y locations along with the desired variables.
      for ( i=1; i <= grid->nn; i++ )
	{
	  //fprintf(fp,"%.15e %.15e %.15e %.15e\n",grid->x[i*2], grid->x[i*2+1], grid->phi[i*8+0],grid->phi[i*8+1]);
	  fprintf(fp,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",grid->x[i*2], grid->x[i*2+1], grid->phi[i*8+0], grid->phi[i*8+2],
		  grid->phi[i*8+4], grid->phi[i*8+6], grid->phi[i*8+1], grid->phi[i*8+3], grid->phi[i*8+5], grid->phi[i*8+7] );
	}

      fprintf(fp,"\n");
  
      // Print out the triangle connectivity to the file.
      for ( i=1; i <= grid->num_elem[Tri]; i++ )
	{
	  n0 = grid->c2n[Tri][i*3];
	  n1 = grid->c2n[Tri][i*3+1];
	  n2 = grid->c2n[Tri][i*3+2];
	  fprintf(fp,"%d %d %d\n",n0,n1,n2);
	}
      
      fclose(fp);
      
      return;
    }
  else
    {
      fprintf(fp,"ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n",grid->nn,
	      (grid->num_elem[Tri])+(grid->num_elem[Quad]));
  
      // Loop over all the nodes and write out the x,y locations along with the desired variables.
      for ( i=1; i <= grid->nn; i++ )
	{
	  //fprintf(fp,"%.15e %.15e %.15e %.15e\n",grid->x[i*2], grid->x[i*2+1], grid->phi[i*8+0],grid->phi[i*8+1]);
	  fprintf(fp,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",grid->x[i*2], grid->x[i*2+1], grid->phi[i*8+0], grid->phi[i*8+2],
		  grid->phi[i*8+4], grid->phi[i*8+6], grid->phi[i*8+1], grid->phi[i*8+3], grid->phi[i*8+5], grid->phi[i*8+7] );
	}
      
      fprintf(fp,"\n");
  
      // Print out element connectivity to the file. For triangles, print the last node twice.
      for ( i=1; i <= grid->num_elem[Tri]; i++ )
	{
	  n0 = grid->c2n[Tri][i*3];
	  n1 = grid->c2n[Tri][i*3+1];
	  n2 = grid->c2n[Tri][i*3+2];
	  fprintf(fp,"%d %d %d %d\n",n0,n1,n2,n2);
	}
      for ( i=1; i <= grid->num_elem[Quad]; i++ )
	{
	  n0 = grid->c2n[Quad][i*4];
	  n1 = grid->c2n[Quad][i*4+1];
	  n2 = grid->c2n[Quad][i*4+2];
	  n3 = grid->c2n[Quad][i*4+3];
	  fprintf(fp,"%d %d %d %d\n",n0,n1,n2,n3);
	}
      
      fclose(fp);
      
      return;
    }

}


//=============================================================
// 
//  write_surface_cp()
//  
//  Writes out a data file of the coeffecient of pressure for all
//  boundary points tagged as inviscid on a per boundary basis.
//
//  GRID *grid;                            // The grid.
//  PARAMS p;                              // The parameters.
//
//=============================================================

void write_surface_cp ( GRID *grid, PARAMS p )
{
  int b, n;                                  // Loop counters.
  int s,i,j;                                 // Loop counters.
  int bc;                                    // Boundary condition.
  int bind;                                  // Boundary index.
  int sind;                                  // Segment index.
  int n0, n1;                                // Edge nodes.
  int n0bedge, n1bedge;                      // Boundary edge indicies for segment nodes.
  int real_node,ghost_node;
  double cp;                                 // Coeffecient of pressure.
  double pinf;                               // Free-stream pressure.
  double minf;                               // Free-stream Mach number.
  double P;                                  // Pressure of CV.
  double r;                                  // Density of CV.
  double u;                                  // X-component of velocity of CV.
  double v;                                  // Y-component of velocity of CV.
  double e;                                  // Total energy of CV.
  double t1,t2,t3;                           // Gauss points and weights.
  double xmid,ymid;
  double x1,x2,x3,y1,y2,y3;                  // Gauss coordinates.
  double xL,xR,yL,yR;
  double qleft[4];
  double gp[2];
  double XL[2],XR[2],GP[6];
  char buffer[100];                          // Static string buffer.
  char sys_buff[256];                        // Buffer for calling system shell executions.
  FILE *fp;                                  // File pointer.

  int *list_of_nodes = NULL;
  int is_unique = 0;
  int nn_unique = 0;

  n0 = 0; n1 = 0; n0bedge = 0; n1bedge = 0;

  printf("//===================================================================\n");
  printf("// Starting SURFACE PRESSURE IO - WRITING DATA FILES\n");
  printf("//===================================================================\n");

  t3 = sqrt(15.0)/5.0;
  t2 = 0.;
  t1 = sqrt(15.0)/(-5.0);
  
  // Loop through all inviscid boundaries and create a separate data file. Then write out the
  // Cp values for all nodes on that boundary.

  pinf = 1./p.gamma;
	  
  minf = p.mach;
  
  for ( b=1; b <= grid->nb; b++ )
    {
      if ( (grid->bbc[b] == 0) || (grid->bbc[b]%2 == 0) )  continue;    // only process inviscid boundaries.
      
      sprintf(buffer,"surface_%d_pressure.dat",b);

      // Open the solution file for writing. Quit if it can't be opened.
      if ((fp = fopen ( buffer ,"w")) == 0 )
	{
	  printf("\nCouldn't open <%s>\n",buffer);
	  exit(0);
	}

      // set up the arrays.
      nn_unique = 0;
      
      list_of_nodes = (int*)malloc( (grid->nbs[b] * 2) * sizeof(int));
      if ( list_of_nodes == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'list_of_nodes'.\n"); exit(1); }

      if ( p.order > 6 )  // Include the gauss points on the surface if the solution is higher order. Actually, don't do this.
	{
	  // Loop over the boundary segments.
	  for ( s=1; s <= grid->nbs[b]; s++ )
	    {
	      n0 = grid->bs[b][s][0];
	      n1 = grid->bs[b][s][1];

	      bc = grid->bbc[b];

	      // Find the correct boundary edge.
	      for ( i=1; i <= grid->nbedges; i++ )
		{
		  real_node = grid->bedges[i*5+0];
		  
		  bind = grid->bedges[i*5+3];
		  sind = grid->bedges[i*5+4];

		  if ( b == bind && s == sind )
		    {
		      if ( real_node == n0 )
			n0bedge = i;
		      else
			n1bedge = i;
		    }
		}

	      // Now we have the correct boundary edges associated with each node on the boundary segment.

	      if ( bc <= 10 )  // Linear boundary.
		{
		  // Get the edge midpoint.
		  xmid = 0.5*( grid->x[n0*NDIM+0] + grid->x[n1*NDIM+0] );
		  ymid = 0.5*( grid->x[n0*NDIM+1] + grid->x[n1*NDIM+1] );

		  // Write out the pressure at the node.

		  // At this point (this function is called at the end of the solver run, so grid->Q == grid->nQ and both are in conserved form.
		  // Since the reconstruction functions use nQ let's convert it to primitives since the derivatives are in primitive form.
		  if ( RECON_PRIM )
		    {
		      ConvertConservedToPrimitive( p.gamma , &(grid->nQ[n0*NUM_VAR]) );
		      P = grid->nQ[n0*NUM_VAR+3];
		    }
		  else
		    {
		      r = grid->Q[n0*NUM_VAR+0];
		      u = grid->Q[n0*NUM_VAR+1] / r;
		      v = grid->Q[n0*NUM_VAR+2] / r;
		      e = grid->Q[n0*NUM_VAR+3];
		      P = (p.gamma - 1.)*(e - 0.5*r*(u*u + v*v));
		    }
		  
		  cp = ( P/pinf - 1. )*(2. / (p.gamma*minf*minf));
		  
		  fprintf(fp,"%f  %f\n",(float)grid->x[n0*NDIM],(float)cp);
		  
		  // Do the left node first.

		  xL = grid->x[n0*NDIM+0];
		  yL = grid->x[n0*NDIM+1];

		  xR = xmid;
		  yR = ymid;

		  x1 = (xL+xR)*0.5 + (xR-xL)*0.5*t1;
		  y1 = (yL+yR)*0.5 + (yR-yL)*0.5*t1;
		  
		  x2 = (xL+xR)*0.5 + (xR-xL)*0.5*t2;
		  y2 = (yL+yR)*0.5 + (yR-yL)*0.5*t2;
		  
		  x3 = (xL+xR)*0.5 + (xR-xL)*0.5*t3;
		  y3 = (yL+yR)*0.5 + (yR-yL)*0.5*t3;

		  gp[0] = x1;  gp[1] = y1;
		  Reconstruct_Gauss_Boundary ( n0, gp, qleft, grid, p);

		  // Get the pressure
		  if ( RECON_PRIM )
		    {
		      P = qleft[3];
		    }
		  else
		    {
		      r = qleft[0];
		      u = qleft[1] / r;
		      v = qleft[2] / r;
		      e = qleft[3];
		      P = (p.gamma - 1.)*(e - 0.5*r*(u*u + v*v));
		    }
		  
		  cp = ( P/pinf - 1. )*(2. / (p.gamma*minf*minf));
		  
		  fprintf(fp,"%f  %f\n",(float)gp[0],(float)cp);

		  gp[0] = x2;  gp[1] = y2;
		  Reconstruct_Gauss_Boundary ( n0, gp, qleft, grid, p);

		  // Get the pressure
		  if ( RECON_PRIM )
		    {
		      P = qleft[3];
		    }
		  else
		    {
		      r = qleft[0];
		      u = qleft[1] / r;
		      v = qleft[2] / r;
		      e = qleft[3];
		      P = (p.gamma - 1.)*(e - 0.5*r*(u*u + v*v));
		    }
		  
		  cp = ( P/pinf - 1. )*(2. / (p.gamma*minf*minf));
		  
		  fprintf(fp,"%f  %f\n",(float)gp[0],(float)cp);

		  gp[0] = x3;  gp[1] = y3;
		  Reconstruct_Gauss_Boundary ( n0, gp, qleft, grid, p);

		  // Get the pressure
		  if ( RECON_PRIM )
		    {
		      P = qleft[3];
		    }
		  else
		    {
		      r = qleft[0];
		      u = qleft[1] / r;
		      v = qleft[2] / r;
		      e = qleft[3];
		      P = (p.gamma - 1.)*(e - 0.5*r*(u*u + v*v));
		    }
		  
		  cp = ( P/pinf - 1. )*(2. / (p.gamma*minf*minf));
		  
		  fprintf(fp,"%f  %f\n",(float)gp[0],(float)cp);

		  // Now convert back.
		  if ( RECON_PRIM )
		    {
		      ConvertPrimitiveToConserved( p.gamma , &(grid->nQ[n0*NUM_VAR]) );
		    }

		  // Now process the right node.
		  if ( RECON_PRIM )
		    {
		      ConvertConservedToPrimitive( p.gamma , &(grid->nQ[n1*NUM_VAR]) );
		    }

		  xL = xmid;
		  yL = ymid;

		  xR = grid->x[n1*NDIM+0];
		  yR = grid->x[n1*NDIM+1];

		  x1 = (xL+xR)*0.5 + (xR-xL)*0.5*t1;
		  y1 = (yL+yR)*0.5 + (yR-yL)*0.5*t1;
		  
		  x2 = (xL+xR)*0.5 + (xR-xL)*0.5*t2;
		  y2 = (yL+yR)*0.5 + (yR-yL)*0.5*t2;
		  
		  x3 = (xL+xR)*0.5 + (xR-xL)*0.5*t3;
		  y3 = (yL+yR)*0.5 + (yR-yL)*0.5*t3;

		  gp[0] = x1;  gp[1] = y1;
		  Reconstruct_Gauss_Boundary ( n1, gp, qleft, grid, p);

		  // Get the pressure
		  if ( RECON_PRIM )
		    {
		      P = qleft[3];
		    }
		  else
		    {
		      r = qleft[0];
		      u = qleft[1] / r;
		      v = qleft[2] / r;
		      e = qleft[3];
		      P = (p.gamma - 1.)*(e - 0.5*r*(u*u + v*v));
		    }
		  
		  cp = ( P/pinf - 1. )*(2. / (p.gamma*minf*minf));
		  
		  fprintf(fp,"%f  %f\n",(float)gp[0],(float)cp);

		  gp[0] = x2;  gp[1] = y2;
		  Reconstruct_Gauss_Boundary ( n1, gp, qleft, grid, p);

		  // Get the pressure
		  if ( RECON_PRIM )
		    {
		      P = qleft[3];
		    }
		  else
		    {
		      r = qleft[0];
		      u = qleft[1] / r;
		      v = qleft[2] / r;
		      e = qleft[3];
		      P = (p.gamma - 1.)*(e - 0.5*r*(u*u + v*v));
		    }

		  cp = ( P/pinf - 1. )*(2. / (p.gamma*minf*minf));
		  
		  fprintf(fp,"%f  %f\n",(float)gp[0],(float)cp);

		  gp[0] = x3;  gp[1] = y3;
		  Reconstruct_Gauss_Boundary ( n1, gp, qleft, grid, p);

		  // Get the pressure
		  if ( RECON_PRIM )
		    {
		      P = qleft[3];
		    }
		  else
		    {
		      r = qleft[0];
		      u = qleft[1] / r;
		      v = qleft[2] / r;
		      e = qleft[3];
		      P = (p.gamma - 1.)*(e - 0.5*r*(u*u + v*v));
		    }
		  
		  cp = ( P/pinf - 1. )*(2. / (p.gamma*minf*minf));
		  
		  fprintf(fp,"%f  %f\n",(float)gp[0],(float)cp);

		  // Write out the right node data.
		  if ( RECON_PRIM )
		    {
		      P = grid->nQ[n1*NUM_VAR+3];
		    }
		  else
		    {
		      r = grid->Q[n1*NUM_VAR+0];
		      u = grid->Q[n1*NUM_VAR+1] / r;
		      v = grid->Q[n1*NUM_VAR+2] / r;
		      e = grid->Q[n1*NUM_VAR+3];
		      P = (p.gamma - 1.)*(e - 0.5*r*(u*u + v*v));
		    }
		  
		  cp = ( P/pinf - 1. )*(2. / (p.gamma*minf*minf));
		  
		  fprintf(fp,"%f  %f\n",(float)grid->x[n1*NDIM],(float)cp);

		  // and finally convert back.
		  if ( RECON_PRIM )
		    {
		      ConvertPrimitiveToConserved( p.gamma , &(grid->nQ[n1*NUM_VAR]) );
		    }
		}
	      else  // Curved boundary.
		{
		  // Do the left node first.

		  // At this point (this function is called at the end of the solver run, so grid->Q == grid->nQ and both are in conserved form.
		  // Since the reconstruction functions use nQ let's convert it to primitives since the derivatives are in primitive form.
		  if ( RECON_PRIM )
		    {
		      ConvertConservedToPrimitive( p.gamma , &(grid->nQ[n0*NUM_VAR]) );
		      P = grid->nQ[n0*NUM_VAR+3];
		    }
		  else
		    {
		      r = grid->Q[n0*NUM_VAR+0];
		      u = grid->Q[n0*NUM_VAR+1] / r;
		      v = grid->Q[n0*NUM_VAR+2] / r;
		      e = grid->Q[n0*NUM_VAR+3];
		      P = (p.gamma - 1.)*(e - 0.5*r*(u*u + v*v));
		    }
		  
		  cp = ( P/pinf - 1. )*(2. / (p.gamma*minf*minf));
		  
		  fprintf(fp,"%f  %f\n",(float)grid->x[n0*NDIM],(float)cp);

		  real_node = grid->bedges[n0bedge*5+0];
		  ghost_node = grid->bedges[n0bedge*5+1] + grid->nn;

		  XL[0] = grid->x[real_node*NDIM+0];
		  XL[1] = grid->x[real_node*NDIM+1];
		  XR[0] = grid->x[ghost_node*NDIM+0];
		  XR[1] = grid->x[ghost_node*NDIM+1];

		  // Get the Gauss points.
		  curved_boundary_gauss_points ( grid, bc, XL, XR, GP );

		  for ( j=0; j < 3; j++ )
		    {
		      Reconstruct_Gauss_Boundary ( real_node, &(GP[j*NDIM]), qleft, grid, p);
		      
		      // Get the pressure
		      if ( RECON_PRIM )
			{
			  P = qleft[3];
			}
		      else
			{
			  r = qleft[0];
			  u = qleft[1] / r;
			  v = qleft[2] / r;
			  e = qleft[3];
			  P = (p.gamma - 1.)*(e - 0.5*r*(u*u + v*v));
			}
		      cp = ( P/pinf - 1. )*(2. / (p.gamma*minf*minf));
		      
		      fprintf(fp,"%f  %f\n",(float)GP[j*NDIM],(float)cp);
		    }

		  // Now convert back.
		  if ( RECON_PRIM )
		    {
		      ConvertPrimitiveToConserved( p.gamma , &(grid->nQ[n0*NUM_VAR]) );
		    }
		  
		  // Work on the right node.
		  
		  if ( RECON_PRIM )
		    {
		      ConvertConservedToPrimitive( p.gamma , &(grid->nQ[n1*NUM_VAR]) );
		    }

		  real_node = grid->bedges[n1bedge*5+0];
		  ghost_node = grid->bedges[n1bedge*5+1] + grid->nn;

		  XR[0] = grid->x[real_node*NDIM+0];
		  XR[1] = grid->x[real_node*NDIM+1];
		  XL[0] = grid->x[ghost_node*NDIM+0];
		  XL[1] = grid->x[ghost_node*NDIM+1];

		  // Get the Gauss points.
		  curved_boundary_gauss_points ( grid, bc, XL, XR, GP );

		  for ( j=0; j < 3; j++ )
		    {
		      Reconstruct_Gauss_Boundary ( real_node, &(GP[j*NDIM]), qleft, grid, p);
		      
		      // Get the pressure
		      if ( RECON_PRIM )
			{
			  P = qleft[3];
			}
		      else
			{
			  r = qleft[0];
			  u = qleft[1] / r;
			  v = qleft[2] / r;
			  e = qleft[3];
			  P = (p.gamma - 1.)*(e - 0.5*r*(u*u + v*v));
			}
		      
		      cp = ( P/pinf - 1. )*(2. / (p.gamma*minf*minf));
		      
		      fprintf(fp,"%f  %f\n",(float)GP[j*NDIM],(float)cp);
		    }

		  // Write the data for the right node.
		  if ( RECON_PRIM )
		    {
		      P = grid->nQ[n1*NUM_VAR+3];
		    }
		  else
		    {
		      r = grid->Q[n1*NUM_VAR+0];
		      u = grid->Q[n1*NUM_VAR+1] / r;
		      v = grid->Q[n1*NUM_VAR+2] / r;
		      e = grid->Q[n1*NUM_VAR+3];
		      P = (p.gamma - 1.)*(e - 0.5*r*(u*u + v*v));
		    }
		  
		  cp = ( P/pinf - 1. )*(2. / (p.gamma*minf*minf));
		  
		  fprintf(fp,"%f  %f\n",(float)grid->x[n1*NDIM],(float)cp);
		  
		  // Now convert back.
		  if ( RECON_PRIM )
		    {
		      ConvertPrimitiveToConserved( p.gamma , &(grid->nQ[n1*NUM_VAR]) );
		    }
		}
	    }
	}
      else                               // In the below code, regardless of RECON_PRIM, grid->Q is in conserved form.
	{
	  for ( n=1; n <= grid->nbs[b]; n++ )
	    {
	      n0 = grid->bs[b][n][0];
	      n1 = grid->bs[b][n][1];

	      // Check that n0 is unique. If it is, process it.
	      is_unique = 1;

	      for ( i=0; i < nn_unique; i++ )
		{
		  if ( n0 == list_of_nodes[i] )
		    is_unique = 0;
		}

	      // if its not in the list skip it. Else, process it.
	      if ( is_unique )
		{
		  list_of_nodes[nn_unique] = n0;
		  nn_unique++;

		  r = grid->Q[n0*NUM_VAR+0];
		  u = grid->Q[n0*NUM_VAR+1] / r;
		  v = grid->Q[n0*NUM_VAR+2] / r;
		  e = grid->Q[n0*NUM_VAR+3];
		  P = (p.gamma - 1.)*(e - 0.5*r*(u*u + v*v));

		  cp = ( P/pinf - 1. )*(2. / (p.gamma*minf*minf));
	      
		  fprintf(fp,"%f  %f\n",(float)grid->x[n0*NDIM],(float)cp);
		}

	      // repeat for the other node.
	      is_unique = 1;

	      for ( i=0; i < nn_unique; i++ )
		{
		  if ( n1 == list_of_nodes[i] )
		    is_unique = 0;
		}

	      if ( is_unique )
		{
		  list_of_nodes[nn_unique] = n1;
		  nn_unique++;
	      
		  r = grid->Q[n1*NUM_VAR+0];
		  u = grid->Q[n1*NUM_VAR+1] / r;
		  v = grid->Q[n1*NUM_VAR+2] / r;
		  e = grid->Q[n1*NUM_VAR+3];
		  P = (p.gamma - 1.)*(e - 0.5*r*(u*u + v*v));

		  cp = ( P/pinf - 1. )*(2. / (p.gamma*minf*minf));
	  
		  fprintf(fp,"%f  %f\n",(float)grid->x[n1*NDIM],(float)cp);
		}
	    }

	  // done with this boundary.
	  freenull( list_of_nodes );
	  
	}
      
      fclose(fp);

      // now lets have the system sort the results in proper numeric order.
      // Warning: This step is non-portable as this will work only for properly configured UNIX type systems.

      // we have closed the file. Now lets copy it to a temp file and then sort it to the original file name.
      sprintf(sys_buff,"cp %s temp_bound_%d.dat",buffer, b);
      system(sys_buff);

      sprintf(sys_buff,"rm %s",buffer);
      system(sys_buff);

      sprintf(sys_buff,"sort -g temp_bound_%d.dat > %s",b,buffer);
      system(sys_buff);

      sprintf(sys_buff,"rm temp_bound_%d.dat",b);
      system(sys_buff);
    }
  
  return;
}

//=============================================================
// 
//  Compute_Body_Forces()
//  
//  Compute the lift and drag along surfaces that are tagged
//  as inviscid.
//
//  GRID *grid;                            // The grid.
//  PARAMS p;                              // The parameters.
//  double *Forces;                        // The body forces.
//
//=============================================================

void Compute_Body_Forces ( GRID *grid, PARAMS p, double *Forces )
{
  int i,j,b;                           // Loop counters.
  int real_node, ghost_node;           // Boundary edge nodes.
  int bc;                              // Boundary condition.
  double nx,ny,len;                    // Edge normal vector.
  double qleft[4];                     // Q on the left of the face.
  double P;                            // Pressure of CV.
  double r;                            // Density of CV.
  double u;                            // X-component of velocity of CV.
  double v;                            // Y-component of velocity of CV.
  double e;                            // Total energy of CV.
  double temp1,temp2;
  double angle;

  // Variables for higher order flux.
  double xmid, ymid;                     // Values at face midpoints.
  double x1,x2,x3;                       // Integration points.
  double y1,y2,y3;
  double xL,xR,yL,yR;
  double t1,t2,t3;                       // Gaussian roots.
  double w1,w2,w3;                       // Gaussian coefficients.
  double xi, yi;                         // Node coordinates.
  int nodeL, nodeR;                      // Nodes comprising an edge.
  int s;                                 // Boundary segment.
  double gp[2];                          // Gauss point.

  // HACK for order==2.
  PARAMS p1;

  if ( p.order == 2 )
    {
      p1.foits = p.foits;
      p1.soits = p.soits;
      p1.order = 1;
    }

  // Curved boundary stuff.
  double XL[2],XR[2],GP[6];
  double w[3];
  
  w[0] = 8./9.;
  w[1] = 5./9.;
  w[2] = w[1];

  Forces[0] = 0.;
  Forces[1] = 0.;

  angle = p.alpha * M_PI/180;  // Get in radians.
  
  // The procedure I'm going to follow is: Loop over the boundary edges, skipping
  // those that aren't marked inviscid, get the Gaussian quadrature points, get the
  // pressure at those points, and integrate. Then accumulate the total forces.

  for ( i=1; i <= grid->nbedges; i++ )
    {
      temp1 = 0.;
      temp2 = 0.;

      // Grab the data from the structure.
      real_node = grid->bedges[i*5+0];
      ghost_node = grid->bedges[i*5+1] + grid->nn;

      b = grid->bedges[i*5+3];
      bc = grid->bbc[b];

      if ( bc == 0 || (bc%2 == 0) )
	continue;
      
      // Get the normal vector information.
      nx = grid->xn_bedges[i*3+0];
      ny = grid->xn_bedges[i*3+1];
      len = grid->xn_bedges[i*3+2];

      // Flip the normal vector since I have outward pointing normals.
      nx = -nx;
      ny = -ny;
      
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

	  gp[0] = x1;  gp[1] = y1;

	  // If we reconstructed the variables as primitives, then nQ needs to be converted from conserved to primitive.

	  if ( RECON_PRIM )
	    {
	      ConvertConservedToPrimitive( p.gamma , &(grid->nQ[real_node*NUM_VAR]) ); // reconstruction uses nQ not Q.
	    }
	  
	  if ( p.order == 2 )
	    Reconstruct_Gauss_Boundary ( real_node, gp, qleft, grid, p1);
	  else
	    Reconstruct_Gauss_Boundary ( real_node, gp, qleft, grid, p);
	  
	  // Now get the pressure.
	  if ( RECON_PRIM )
	    {
	      P = qleft[3];
	    }
	  else
	    {
	      r = qleft[0];
	      u = qleft[1] / r;
	      v = qleft[2] / r;
	      e = qleft[3];
	      P = (p.gamma - 1.)*(e - 0.5*r*(u*u + v*v));
	    }

	  // Compute the body force contribution.
	  temp1 += ( w1 * ( P*ny*cos(angle) - P*nx*sin(angle) ) );
	  temp2 += ( w1 * ( P*nx*cos(angle) + P*ny*sin(angle) ) );

	  gp[0] = x2;  gp[1] = y2;
	  
	  if ( p.order == 2 )
	    Reconstruct_Gauss_Boundary ( real_node, gp, qleft, grid, p1);
	  else
	    Reconstruct_Gauss_Boundary ( real_node, gp, qleft, grid, p);

	  // Now get the pressure.
	  if ( RECON_PRIM )
	    {
	      P = qleft[3];
	    }
	  else
	    {
	      r = qleft[0];
	      u = qleft[1] / r;
	      v = qleft[2] / r;
	      e = qleft[3];
	      P = (p.gamma - 1.)*(e - 0.5*r*(u*u + v*v));
	    }

	  // Compute the body force contribution.
	  temp1 += ( w2 * ( P*ny*cos(angle) - P*nx*sin(angle) ) );
	  temp2 += ( w2 * ( P*nx*cos(angle) + P*ny*sin(angle) ) );

	  gp[0] = x3;  gp[1] = y3;
	  
	  if ( p.order == 2 )
	    Reconstruct_Gauss_Boundary ( real_node, gp, qleft, grid, p1);
	  else
	    Reconstruct_Gauss_Boundary ( real_node, gp, qleft, grid, p);

	  // Now get the pressure.
	  if ( RECON_PRIM )
	    {
	      P = qleft[3];
	    }
	  else
	    {
	      r = qleft[0];
	      u = qleft[1] / r;
	      v = qleft[2] / r;
	      e = qleft[3];
	      P = (p.gamma - 1.)*(e - 0.5*r*(u*u + v*v));
	    }
	  
	  // Compute the body force contribution.
	  temp1 += ( w3 * ( P*ny*cos(angle) - P*nx*sin(angle) ) );
	  temp2 += ( w3 * ( P*nx*cos(angle) + P*ny*sin(angle) ) );

	  // Integrate the edge.
	  temp1 = temp1 * len * 0.5;
	  temp2 = temp2 * len * 0.5;

	  // Convert the contents of nQ for real_node back.
	  if ( RECON_PRIM )
	    {
	      ConvertPrimitiveToConserved( p.gamma , &(grid->nQ[real_node*NUM_VAR]) );
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

	  if ( RECON_PRIM )
	    {
	      ConvertConservedToPrimitive( p.gamma , &(grid->nQ[real_node*NUM_VAR]) ); // reconstruction uses nQ not Q.
	    }

	  for ( j=0; j < 3; j++ )
	    {
	      // Get the normal vector.
	      curved_boundary_normal_vector ( grid, bc, &(GP[j*NDIM]), &nx, &ny );

	      // flip normals.
	      nx = -nx;
	      ny = -ny;
	      
	      if ( p.order == 2 )
		Reconstruct_Gauss_Boundary ( real_node, &(GP[j*NDIM]), qleft, grid, p1);
	      else
		Reconstruct_Gauss_Boundary ( real_node, &(GP[j*NDIM]), qleft, grid, p);
		  
	      // Now get the pressure.
	      if ( RECON_PRIM )
		{
		  P = qleft[3];
		}
	      else
		{
		  r = qleft[0];
		  u = qleft[1] / r;
		  v = qleft[2] / r;
		  e = qleft[3];
		  P = (p.gamma - 1.)*(e - 0.5*r*(u*u + v*v));
		}

	      // Compute the body force contribution.
	      temp1 += ( w[j] * ( P*ny*cos(angle) - P*nx*sin(angle) ) );
	      temp2 += ( w[j] * ( P*nx*cos(angle) + P*ny*sin(angle) ) );
	    }

	  // Integrate the edge.
	  temp1 = temp1 * len * 0.5;
	  temp2 = temp2 * len * 0.5;

	  // Convert the contents of nQ for real_node back.
	  if ( RECON_PRIM )
	    {
	      ConvertPrimitiveToConserved( p.gamma , &(grid->nQ[real_node*NUM_VAR]) );
	    }
	}

      // Done with this edge, accumulate the forces.
      Forces[0] += temp1;
      Forces[1] += temp2;

    }

  return;
}


//=============================================================
// 
//  mesh_io.C
//  
//  Functions needed for reading and writing
//  a mesh file in Dr. Steve Karman's format.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include "params.h"
#include "grid.h"
#include "mesh_io.h"
#include "Defined_Variables.h"

//=============================================================
// 
//  read_mesh()
//
//  Reads in the mesh file.
//  
//  char *file_name                   // The file to be read in.
//  GRID *grid                        // The grid.
//  PARAMS p;                         // The parameters.
//
//=============================================================

void read_mesh ( char *file_name, GRID *grid, PARAMS p )
{
  int i,n,t,b,s;                         // Loop counters.
  FILE *fp;                              // Pointer to the file.
  const int bdim = 132;                  // Default array dimension size.
  char buff[bdim];                       // Default buffer for reading in from a stream.
  int ns=0;                              // Number of boundary segments in total.
  int ecount=1;                          // Global tracker of boundary segments.
  double cfl;                            // CFL number.

  //Open the file, but kick out if the file could not be opened.
  if (( fp = fopen(file_name,"r")) == 0)
    {
      printf("\nCouldn't open <%s>\n",file_name);
      exit(0);
    }

  printf("//==================================================\n");
  printf("// Starting MESH IO - READING MESH <%s>\n",file_name);
  printf("//==================================================\n");

  //Read number of nodes from grid.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(grid->nn));
  printf("  Number of points = %d",grid->nn);
  fflush(stdout);

  //Allocate space for the coordinates and for the control volume displacements as well as other data
  //members that depend on number of nodes.

  grid->node_state = (int*)malloc( (grid->nn+1)*sizeof(int));
  if ( grid->node_state == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'node_state'.\n"); exit(1); }

  grid->node_order = (int*)malloc( (grid->nn+1)*sizeof(int));
  if ( grid->node_order == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'node_order'.\n"); exit(1); }

  grid->x = (double*)malloc( (grid->nn+1)*NDIM * sizeof(double) );
  if ( grid->x == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'x'.\n"); exit(1); }

  grid->Q = (double*)malloc( (grid->nn+1)*NUM_VAR * sizeof(double) );
  if ( grid->Q == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'Q'.\n"); exit(1); }

  grid->nQ = (double*)malloc( (grid->nn+1)*NUM_VAR * sizeof(double) );
  if ( grid->nQ == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'nQ'.\n"); exit(1); }

  grid->pQ = (double*)malloc( (grid->nn+1)*NUM_VAR * sizeof(double) );
  if ( grid->nQ == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'pQ'.\n"); exit(1); }

  grid->dQ = (double*)malloc( (grid->nn+1)*NUM_VAR * sizeof(double) );
  if ( grid->dQ == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'dQ'.\n"); exit(1); }

  grid->dt = (double*)malloc( (grid->nn+1)*NUM_VAR * sizeof(double) );
  if ( grid->dt == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'dt'.\n"); exit(1); }

  grid->R = (double*)malloc( (grid->nn+1)*NUM_VAR * sizeof(double) );
  if ( grid->R == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'R'.\n"); exit(1); }

  // Allocate space for the element types in c2n.
  grid->c2n = (int**)malloc(NUM_ELEM_TYPE*sizeof(int*));
  if ( grid->c2n == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'c2n'.\n"); exit(1); }
  
  // Initialize the pointers.
  for ( i=Edge; i <= Quad; i++ )
    {
      grid->c2n[i] = NULL;
    }

  //Read in coordinates.
  for (n=1; n <= grid->nn; n++)
  {
    fgets(buff,bdim,fp);
    sscanf(buff,"%lf %lf",&(grid->x[n*NDIM]),&(grid->x[n*NDIM+1]));
  }

  //Read in number of blocks.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff, "%d",&i);
  printf("\n  Number of blocks = %d",i);
  fflush(stdout);

  if ( i != 1 )
    {
      printf("FATAL ERROR: Multiple blocks in the mesh format is not supported!\n");
      exit(1);
    }

  //Read in number of triangles.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(grid->num_elem[Tri]));
  printf("\n  Number of triangles = %d",grid->num_elem[Tri]);
  fflush(stdout);

  //Allocate space for triangle connectivity.
  grid->c2n[Tri] = (int*)malloc(((grid->num_elem[Tri]+1)*3)*sizeof(int));
  if ( grid->c2n[Tri] == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'c2n[tri]'.\n"); exit(1); }

  //Read in triangle connectivity retaining the Fortran indexing.
  for (t=1; t <= grid->num_elem[Tri]; t++)
    {
      fgets(buff,bdim,fp);
      sscanf(buff,"%d %d %d",&(grid->c2n[Tri][t*3+0]),&(grid->c2n[Tri][t*3+1]),&(grid->c2n[Tri][t*3+2]));
    }

  //Read in number of quads.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(grid->num_elem[Quad]));
  printf("\n  Number of quads = %d",grid->num_elem[Quad]);
  fflush(stdout);
  
  //Allocate space for quad connectivity.
  grid->c2n[Quad] = (int*)malloc(((grid->num_elem[Quad]+1)*4)*sizeof(int));
  if ( grid->c2n[Quad] == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'c2n[quad]'.\n"); exit(1); }

  //Read in quad connectivity retaining the Fortran indexing.
  for (t=1; t <= grid->num_elem[Quad]; t++)
    {
      fgets(buff,bdim,fp);
      sscanf(buff,"%d %d %d %d",&(grid->c2n[Quad][t*4+0]),&(grid->c2n[Quad][t*4+1]),
	                        &(grid->c2n[Quad][t*4+2]),&(grid->c2n[Quad][t*4+3]));
    }

  //Read in number of boundaries.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(grid->nb));
  printf("\n  Number of boundaries = %d",grid->nb);

  //Allocate space for segments per boundary and boundary storage.
  grid->nbs = (int*)malloc((grid->nb+1)*sizeof(int));
  if ( grid->nbs == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'nbs'.\n"); exit(1); }

  grid->bs = (int***)malloc((grid->nb+1)*sizeof(int**));
  if ( grid->bs == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'bs'.\n"); exit(1); }

  //Loop through for each boundary and store the segments that make up that boundary.
  for (b=1; b <= grid->nb; b++)
  {
    fgets(buff,bdim,fp);
    fgets(buff,bdim,fp);
    sscanf(buff,"%d",&(grid->nbs[b]));
  
    printf("\n  Boundary %d has %d segment(s).",b,grid->nbs[b]);

    grid->bs[b] = (int**)malloc((grid->nbs[b]+1)*sizeof(int*));
    if ( grid->bs[b] == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'bs[%d]'.\n",b); exit(1); }
    
    //Loop through the segments and store the node ids retaining Fortran indexing.
    for (i=1; i <= grid->nbs[b]; i++)
    {
      grid->bs[b][i] = (int*)malloc(2*sizeof(int));
      if ( grid->bs[b][i] == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'bs[%d][%d]'.\n",b,i); exit(1); }

      fgets(buff,bdim,fp);
      sscanf(buff,"%d %d",&(grid->bs[b][i][0]), &(grid->bs[b][i][1]));
    }
  }

  // Copy the ***bs structure to c2n[Edges]. First count up the total number of segments.
  for ( b=1; b <= grid->nb; b++ )
    {
      ns += grid->nbs[b];
    }

  grid->num_elem[Edge] = ns;

  // Allocate space for the boundary edges. The information stored is left node, right node, and boundary ID.
  grid->c2n[Edge] = (int*)malloc((ns+1)*2*sizeof(int));
  if ( grid->c2n[Edge] == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'c2n[Edge]'.\n"); exit(1); }

  // Now copy the information over.
  for ( b=1; b <= grid->nb; b++ )
    {
      for ( s=1; s <= grid->nbs[b]; s++ )
	{
	  grid->c2n[Edge][ecount*2+0] = grid->bs[b][s][0];
	  grid->c2n[Edge][ecount*2+1] = grid->bs[b][s][1];

	  // Increment the number of Edges.
	  ecount++;
	}
    }

  // Check to see if its a restart.
  if ( p.isrestart == 1 )
    {
      // Read in the solution.

      // Get the number of constants.
      fgets(buff,bdim,fp);
      fgets(buff,bdim,fp);  // Right now I'm not going to do anything with this. So it better be zero.

      fgets(buff,bdim,fp); // iteration
      
      // This should now be two. Current iteration and CFL number.
      sscanf(buff,"%d",&i);
      printf("\nThe last iteration finished is %d\n",i);
      grid->citer = i;
      grid->rfiter = i;

      fgets(buff,bdim,fp); // CFL
      sscanf(buff,"%lf",&cfl);
      printf("The last CFL number applied is %.15E\n",cfl);
      grid->cfl = cfl;

      // Get the number of variables for each node. Should be NUM_VAR!
      fgets(buff,bdim,fp);
      fgets(buff,bdim,fp);
      sscanf(buff,"%d",&i);

      if ( i != NUM_VAR )
	{
	  printf("CRITICAL ERROR: In read_mesh() while reading in the restart solution, the number of variables is not NUM_VAR. n = %d.\n",i);
	  fflush(stdout);
	}

      for ( i=1; i <= grid->nn; i++ )
	{
	  fgets(buff,bdim,fp);
	  sscanf(buff,"%lf %lf %lf %lf",&(grid->Q[i*NUM_VAR+0]),&(grid->Q[i*NUM_VAR+1]),&(grid->Q[i*NUM_VAR+2]),&(grid->Q[i*NUM_VAR+3]));
	}

      // Copy the restart values to nQ as well.
      for ( i=1; i <= grid->nn; i++ )
	{
	  grid->nQ[i*NUM_VAR+0] = grid->Q[i*NUM_VAR+0];
	  grid->nQ[i*NUM_VAR+1] = grid->Q[i*NUM_VAR+1];
	  grid->nQ[i*NUM_VAR+2] = grid->Q[i*NUM_VAR+2];
	  grid->nQ[i*NUM_VAR+3] = grid->Q[i*NUM_VAR+3];
	}

      if ( grid->phi == NULL )
	{
	  grid->phi = (double*)malloc( (grid->nn+1)*2*4*sizeof(double) );
	  if ( grid->phi == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'phi'.\n"); exit(1); }
	}

      for ( i=1; i <= grid->nn; i++ )
	{
	  for ( n=0; n < 8; n++ )
	    grid->phi[i*8+n] = 1.;
	}
      
    }

  fclose(fp);
  
  // Now read in the boundary conditions from the bc file.
  if (( fp = fopen("bc.dat","r")) == 0)
    {
      printf("\nCouldn't open <%s>\n","bc.dat");
      exit(0);
    }

  grid->bbc = (int*)malloc((grid->nb+1)*sizeof(int));
  if ( grid->bbc == NULL ) { printf("MEMORY ERROR: COULD NOT ALLOCATE 'bbc'.\n"); exit(1); }
  
  for ( b=1; b <= grid->nb; b++ )
    {
      fgets(buff,bdim,fp);
      sscanf(buff,"%d",&(grid->bbc[b]));
    }

  fclose(fp);
  
  printf("\n//==================================================\n");
  printf("// Finished MESH IO.\n");
  printf("//==================================================\n\n");

  return;
}


//=============================================================
// 
//  write_mesh()
//  
//  Writes out a mesh file including the solution.
//
//  char *file_name                   // The file name to write out.
//  GRID *grid                        // The grid.
//
//=============================================================

void write_mesh ( char *file_name, GRID *grid )
{
  int i,n,b;                             // Loop counters.
  FILE *fp;                              // Pointer to the file.

  //Open the file, but kick out if the file could not be opened.
  if (( fp = fopen(file_name,"w")) == 0)
    {
      printf("\nCouldn't open <%s>\n",file_name);
      exit(0);
    }

  printf("//==================================================\n");
  printf("// Starting MESH IO - WRITING MESH <%s>\n",file_name);
  printf("//==================================================\n");

  //Write the new physical mesh file.
  fprintf(fp,"#Number of grid points\n");
  fprintf(fp,"%d\n",grid->nn);
  for (n=1; n <= grid->nn; n++)
    {
      fprintf(fp,"%.15e %.15e\n",grid->x[n*NDIM],grid->x[n*NDIM+1]);
    }
  
  fprintf(fp,"#Number of blocks\n");
  fprintf(fp,"%d\n",1);
	  
  fprintf(fp,"#Number of triangles\n");
  fprintf(fp,"%d\n",grid->num_elem[Tri]);

  //Write triangle connectivity.
  for (i=1; i <= grid->num_elem[Tri]; i++)
    {
      fprintf(fp,"%d %d %d\n",grid->c2n[Tri][i*3+0],grid->c2n[Tri][i*3+1],grid->c2n[Tri][i*3+2]);
    }
  
  fprintf(fp,"#Number of quads\n");
  fprintf(fp,"%d\n",grid->num_elem[Quad]);

  //Write quad connectivity.
  for (i=1; i <= grid->num_elem[Quad]; i++)
    {
      fprintf(fp,"%d %d %d %d\n",grid->c2n[Quad][i*4+0],grid->c2n[Quad][i*4+1],
	                         grid->c2n[Quad][i*4+2],grid->c2n[Quad][i*4+3]);
    }

  fprintf(fp,"#Number of boundaries\n");
  fprintf(fp,"%d\n",grid->nb);
  
  //Write boundary information.
  for (b=1; b <= grid->nb; b++)
    {
      fprintf(fp,"#Number of edges for boundary %d\n",b);
      fprintf(fp,"%d\n",grid->nbs[b]);
      
      //Write nodes making up this segment.
      for (i=1; i <= grid->nbs[b]; i++)
	{
	  fprintf(fp,"%d %d\n",(grid->bs[b][i][0]),(grid->bs[b][i][1]));
	}
    }

  fprintf(fp,"#Number of Constants\n");
  fprintf(fp,"%d\n",2);
  fprintf(fp,"%d\n",grid->citer);
  fprintf(fp,"%.15E\n",grid->cfl);

  fprintf(fp,"#Number of Variables\n");
  fprintf(fp,"%d\n",NUM_VAR);
  for (n=1; n <= grid->nn; n++)
    {
      fprintf(fp,"%.15e %.15e %.15e %.15e\n",grid->Q[n*NUM_VAR+0],grid->Q[n*NUM_VAR+1],grid->Q[n*NUM_VAR+2],grid->Q[n*NUM_VAR+3]);
    }

  fclose(fp);
  
  return;
}

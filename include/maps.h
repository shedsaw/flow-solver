//=============================================================
// 
//  maps.h
//  
//  Function prototypes for map building functions.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include "Defined_Variables.h"
#include "mesh_utils.h"

#ifndef MAPS_H
#define MAPS_H

void Cells_Surrounding_Node (const int nnodes,
			     int *num_elem,
			     int **c2n,
			     int **csn,
			     int **ncsn);

void Nodes_Surrounding_Node ( const int nnodes,
			      int *num_elem,
			      int **c2n,
			      int *csn,
			      int *ncsn,
			      int **nsn,
			      int **nnsn);

void Edges_Surrounding_Node ( int nnodes,
			      int *nsn,
			      int *nnsn,
			      int *nedges,
			      int **edges,
			      int **esn);

void Build_Cell_to_Edge_Map ( int **c2n,
			      int *num_elem,
			      int *nsn,
			      int *nnsn,
			      int *esn,
			      int ***c2e);

void Elements_Surrounding_Element ( int *num_elem,
				    int **c2n,
				    int *csn,
				    int *ncsn,
				    int **ese,
				    int **nese);

void Nodes_Surrounding_Node2 ( const int nnodes,
			       int *num_elem,
			       int **c2n,
			       int *esn,
			       int *nnsn,
			       int *edges,
			       int **nsn2,
			       int **nnsn2);

void Initialize_Generic_Nodes_Surrounding_Node ( const int nnodes,
						 int *num_elem,
						 int **c2n,
						 int *esn,
						 int *nnsn,
						 int *edges,
						 int **gnsn,
						 int **gnnsn,
						 int **gnnsn1,
						 int **node_fill_level );

void Add_Layer_to_GNSN ( const int nnodes,
			 int *num_elem,
			 int **c2n,
			 int *esn,
			 int *nnsn,
			 int *edges,
			 int **gnsn,
			 int **gnnsn,
			 int *gnnsn1,
			 int *node_fill_level,
			 int the_node );



#endif

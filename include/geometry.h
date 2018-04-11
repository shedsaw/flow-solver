//=============================================================
// 
//  geometry.h
//  
//  Function prototypes for geometry building functions.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include "Defined_Variables.h"
#include "grid.h"

#ifndef GEOMETRY_H
#define GEOMETRY_H


void Find_Element_Centroids ( GRID *grid );

void Build_Edge_Data ( GRID *grid );

void Build_Boundary_Edge_Data ( GRID *grid );

void Fix_Subedge_Data ( GRID *grid );

void Generate_Gaussian_Data_For_Edge ( GRID *grid );

#endif

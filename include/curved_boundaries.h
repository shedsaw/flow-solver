//=============================================================
// 
//  curved_boundaries.h
//  
//  Function prototypes for calculations on curved boundaries.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include "Defined_Variables.h"
#include "mesh_utils.h"
#include "grid.h"
#include "params.h"

#ifndef CURVED_BOUNDARIES_H
#define CURVED_BOUNDARIES_H

double naca_function ( double x, int bc );

double naca_arclen_function ( double x, void *params );

void curved_boundary_midpoint ( GRID *grid, int bc, double *xL, double *xR, double *xM );

void curved_boundary_gauss_points ( GRID *grid, int bc, double *xL, double *xR, double *GP );

void curved_boundary_normal_vector ( GRID *grid, int bc, double *x, double *nx, double *ny );

void curved_boundary_full_normal_vector ( GRID *grid, int bc, double *x, double *nx, double *ny );

void curved_boundary_arclength ( GRID *grid, int bc, double *xL, double *xR, double *arclength );

#endif

//=============================================================
// 
//  ic.h
//  
//  Function prototypes for the initial condition functions.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include "Defined_Variables.h"
#include "grid.h"
#include "params.h"

#ifndef IC_H
#define IC_H

void init_pv ( GRID *grid, PARAMS p );

void init_cv ( GRID *grid, PARAMS p );

void ic_0 ( PARAMS p, double *x, double *q );

void ic_1 ( PARAMS p, double *x, double *q );

void ic_2 ( PARAMS p, double *x, double *q );

void ic_3 ( PARAMS p, double *x, double *q );

void initialize_vortex_pv ( PARAMS p, GRID *grid, double *Q );

void initialize_vortex_yee_pv ( PARAMS p, GRID *grid, double *Q );

void initialize_vortex_cv ( PARAMS p, GRID *grid, double *Q );

void initialize_vortex_yee_cv ( PARAMS p, GRID *grid, double *Q );

void get_core_pressure ( GRID *grid, PARAMS p, double *time, double *core_pressure );

void get_core_pressure_yee ( GRID *grid, PARAMS p, double *time, double *core_pressure );

void write_vortex_pressure_centerline ( GRID *grid, PARAMS p, double *time );

void write_vortex_pressure_centerline_yee ( GRID *grid, PARAMS p, double *time );

#endif

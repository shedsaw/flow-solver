//=============================================================
// 
//  grid_io.h
//  
//  Prototype all functions needed for writing solution data.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include "Defined_Variables.h"
#include "grid.h"
#include "params.h"

#ifndef GRID_IO_H
#define GRID_IO_H

void write_gnuplot ( char *file_name, GRID *grid );

void write_tecplot_solutionC ( char *file_name, GRID *grid );

void write_tecplot_solutionP ( char *file_name, GRID *grid, double gamma );

void write_tecplot_solutionE ( char *file_name, GRID *grid );

void write_tecplot_cv_area ( char *file_name, GRID *grid );

void write_tecplot_cv_limiter ( char *file_name, GRID *grid );

void write_tecplot_derivs ( char *file_name, GRID *grid );

void write_surface_cp ( GRID *grid, PARAMS p );

void Compute_Body_Forces ( GRID *grid, PARAMS p, double *Forces );

#endif

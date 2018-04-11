//=============================================================
// 
//  cv_calc.h
//  
//  Function prototypes for calculations on control volumes.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include "Defined_Variables.h"
#include "mesh_utils.h"
#include "grid.h"

#ifndef CV_CALC_H
#define CV_CALC_H

void cv_calc_area_Green ( GRID *grid );

void cv_calc_area_Tri_Decomp ( GRID *grid );

void cv_calc_check_closure( GRID *grid );

void cv_calc_length_scale ( GRID *grid );

void cv_calc_eigenvalue_contribution ( GRID *grid, PARAMS p );

void cv_calc_Moments ( GRID *grid );

#endif

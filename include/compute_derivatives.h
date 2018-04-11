//=============================================================
// 
//  compute_derivatives.h
//  
//  Function prototypes to generate the first/second/third 
//  derivatives of the conserved/primitive variables.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include "Defined_Variables.h"

#ifndef COMPUTE_DERIVATIVES_H
#define COMPUTE_DERIVATIVES_H

void Compute_Derivatives ( GRID *grid, PARAMS params, int order_desired );

#endif

//=============================================================
// 
//  gradient.h
//  
//  Function prototypes to generate the gradient field of the
//  conserved variables.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include "Defined_Variables.h"

#ifndef GRADIENT_H
#define GRADIENT_H

void Compute_Gradient_Green( GRID *grid );

void Compute_Gradient_LeastSquares( GRID *grid );

void Compute_Edge_Weights( GRID *grid );

#endif

//=============================================================
// 
//  hessian.h
//  
//  Function prototypes to generate the first/second derivatives
//  of the conserved variables.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include "Defined_Variables.h"

#ifndef HESSIAN_H
#define HESSIAN_H

void Compute_Hessian ( GRID *grid , PARAMS params);

double Test_Function ( double *x );

double Test_Function_Divergent ( double *x );

void Test_Function_Derivatives ( double *x, double *derivs );

void Compute_CV_Averages ( GRID *grid, double *cv_avg );

void Get_Annulus_Derivatives ( double *X, double *derivs );

void Annulus_Test_Function ( double *x, double *Q );

void Test_Annulus_Derivatives ( GRID *grid , PARAMS p );

#endif


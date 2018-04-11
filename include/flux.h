//=============================================================
// 
//  flux.h
//  
//  Function prototypes for computing the flux.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include "Defined_Variables.h"
#include "grid.h"
#include "params.h"
#include "sys_mem.h"

#ifndef FLUX_H
#define FLUX_H

void Incompressible_Euler_Flux ( double *Q, double *Eulerflux, double nx, double ny, double len, double gamma );

int Roe_flux ( double nx, double ny, double len, double gamma,
	       double qleft[4], double qright[4], double flux[4]);

int Roe_flux_centered ( double nx, double ny, double len, double gamma,
			double qleft[4], double qright[4], double flux[4]);

int Van_Leer_flux ( double nx, double ny, double len, double gamma,
		    double qleft[4], double qright[4],
		    double flux[4], double dfdqp[][4], double dfdqm[][4] );

#endif

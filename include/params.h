//=============================================================
// 
//  params.h
//  
//  A simple class that will contain the data from the parameter
//  file.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include "Defined_Variables.h"

#ifndef PARAMS_H
#define PARAMS_H

class PARAMS
{
 public:

  PARAMS()
    {
      nsteps = 0;
      iter = 0;
      jufreq = 0;
      CFL = 1.;
      ramp = 0;
      gamma = 1.4;
      alpha = 0.;
      mach = 0.;
      subits = 20;
      unsteady = 0;
      tacc = 0;
      mintime = 0.;
      method = 0;
      order = 0;
      foits = 0;
      soits = 0;
      toits = 0;
      isrestart = 0;
      limit = 0;
      reorder = 0;
      ic = 0;
      p_or_cv = 0;
      start_post_pro = 0;
      post_pro_freq = 0;
    }

  int nsteps;            // Number of timesteps.
  int iter;              // Number of newton iterations.
  int jufreq;            // Jacobian update frequency.
  double CFL;            // Desired CFL number.
  int ramp;              // Number of iterations to ramp the CFL over.
  double gamma;          // Ratio of specific heats.
  double alpha;          // Angle of attack.
  double mach;           // Mach number.
  int subits;            // Number of iterations of the linear solver.
  int unsteady;          // Indicates steady or unsteady flow.
  int tacc;              // Level of time accuracy. 0-steady, 1- 1st, 2 - 2nd
  double mintime;        // Time step for unsteady flow.
  int method;            // Indicates the time marching scheme.
  int order;             // Indicates the order of accuracy for the spatial scheme.
  int foits;             // Number of first order iterations.
  int soits;             // Number of second order iterations.
  int toits;             // Number of third order iterations.
  int isrestart;         // Indicate restart.
  int limit;             // Indicates the limiting scheme.
  int reorder;           // Indicate vertex reordering.
  int ic;                // The initial conditions to use.
  int p_or_cv;           // Indicates point value or cv averages for the ic.
  int start_post_pro;    // Which iteration to start post processing.
  int post_pro_freq;     // Frequency of post processing.
  
};


void read_parameters(PARAMS *p);


#endif

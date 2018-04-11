//=============================================================
// 
//  params.C
//  
//  Functions to handle the parameter file.
//
//  Written by - Shane Sawyer
//
//=============================================================

#include <stdio.h>
#include <stdlib.h>
#include "Defined_Variables.h"
#include "params.h"


//=============================================================
// 
//  read_parameters()
//
//  Reads in the parameter file.
//  
//  PARAMS p;                        // The parameter object.
//
//=============================================================

void read_parameters(PARAMS *p)
{
  FILE *fp;                              // File pointer.
  const int bdim = 132;                  // Default array dimension size.
  char buff[bdim];                       // Default buffer for reading in from a stream.

  //Open the file, but kick out if the file could not be opened.
  if (( fp = fopen("euler.params","r")) == 0)
    {
      printf("\nCouldn't open <%s>\n","euler.params");
      exit(0);
    }

  printf("//==================================================\n");
  printf("// Reading in the parameter file.\n");
  printf("//==================================================\n");

  //Read number of time steps.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(p->nsteps));
  printf("  Number of time steps = %d\n",p->nsteps);
  fflush(stdout);

  //Read number of newton iterations.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(p->iter));
  printf("  Number of newton iterations = %d\n",p->iter);
  fflush(stdout);

  //Read in the jacobian update frequency.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(p->jufreq));
  printf("  Jacobian Update Frequency = %d\n",p->jufreq);
  fflush(stdout);

  //Read in CFL number.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%lf",&(p->CFL));
  printf("  CFL Number = %f\n",(float)p->CFL);
  fflush(stdout);

  //Read number of iterations to ramp CFL over.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(p->ramp));
  printf("  Number of iterations to ramp = %d\n",p->ramp);
  fflush(stdout);

  //Read in gamma.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%lf",&(p->gamma));
  printf("  Gamma = %f\n",(float)p->gamma);
  fflush(stdout);

  //Read in angle of attack.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%lf",&(p->alpha));
  printf("  Alpha = %f\n",(float)p->alpha);
  fflush(stdout);

  //Read Mach number.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%lf",&(p->mach));
  printf("  Mach Number = %f\n",(float)p->mach);
  fflush(stdout);

  //Read number of sub-iterations.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(p->subits));
  printf("  Number of subiterations = %d\n",p->subits);
  fflush(stdout);

  //Read in steady/unsteady.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(p->unsteady));
  if ( p->unsteady == 0 )
    {
      printf("  Flow is steady-state.\n");
    }
  else
    {
      printf("  Flow is unsteady.\n");
    }
  fflush(stdout);

  //Read in the order of time.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(p->tacc));
  printf("  Time Accuracy order = %d\n",p->tacc);
  fflush(stdout);

  //Read the minimum time step (for unsteady flow).
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%lf",&(p->mintime));
  printf("  Time step = %f\n",(float)p->mintime);
  fflush(stdout);
  
  //Read in method.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(p->method));
  printf("  Method = %d\n",p->method);
  fflush(stdout);

  //Read in order.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(p->order));
  printf("  Order = %d\n",p->order);
  fflush(stdout);

  //Read number of first order iterations.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(p->foits));
  printf("  Number of first order iterations = %d\n",p->foits);
  fflush(stdout);

  //Read number of second order iterations.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(p->soits));
  printf("  Number of second order iterations = %d\n",p->soits);
  fflush(stdout);

  //Read number of third order iterations.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(p->toits));
  printf("  Number of third order iterations = %d\n",p->toits);
  fflush(stdout);

  //Read in isrestart.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(p->isrestart));
  printf("  Isrestart = %d\n",p->isrestart);
  fflush(stdout);

  //Read in limit.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(p->limit));
  printf("  Limit = %d\n",p->limit);
  fflush(stdout);

  //Read in reorder.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(p->reorder));
  printf("  Reorder = %d\n",p->reorder);
  fflush(stdout);

  //Read in initial conditions.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(p->ic));
  printf("  IC = %d\n",p->ic);
  fflush(stdout);

  //Read in p_or_cv.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(p->p_or_cv));
  printf("  Point or CV values = %d\n",p->p_or_cv);
  fflush(stdout);

  //Read in when to start writing out unsteady files.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(p->start_post_pro));
  printf("  Start unsteady file write at iteration = %d\n",p->start_post_pro);
  fflush(stdout);

  //Read in frequency of solution writes.
  fgets(buff,bdim,fp);
  fgets(buff,bdim,fp);
  sscanf(buff,"%d",&(p->post_pro_freq));
  printf("  Frequency of post processing = %d\n",p->post_pro_freq);
  fflush(stdout);

  fclose(fp);
  
  printf("\n//==================================================\n");
  printf("// Finished reading parameters.\n");
  printf("// Caution: Some values may be reset if they conflict!\n");
  printf("//==================================================\n\n");

  /*
  printf("  Number of iterations = %d\n",p->nsteps);
  printf("  Number of Newton iterations = %d\n",p->iter);
  printf("  Jacobian Update frequency = %d\n",p->jufreq);
  printf("  CFL Number = %f\n",(float)p->CFL);
  printf("  Number of iterations to ramp = %d\n",p->ramp);
  printf("  Gamma = %f\n",(float)p->gamma);
  printf("  Alpha = %f\n",(float)p->alpha);
  printf("  Mach Number = %f\n",(float)p->mach);
  printf("  Number of subiterations = %d\n",p->subits);
  printf("  Unsteady flag = %d\n",p->unsteady);
  printf("  Order of time accuracy = %d\n",p->tacc);
  printf("  Minimum time step = %f\n",(float)p->mintime);
  printf("  Method = %d\n",p->method);
  printf("  Order = %d\n",p->order);
  printf("  Number of first order iterations = %d\n",p->foits);
  printf("  Number of second order iterations = %d\n",p->soits);
  printf("  Number of third order iterations = %d\n",p->toits);
  printf("  Isrestart = %d\n",p->isrestart);
  printf("  Limit = %d\n",p->limit);
  printf("  Reorder = %d\n",p->reorder);
  printf("  IC = %d\n",p->ic);
  printf("  Point or CV values = %d\n",p->p_or_cv);
  printf("  Start unsteady file write at iteration = %d\n",p->start_post_pro);
  printf("  Frequency of post processing = %d\n",p->post_pro_freq);
  */

  return;
}

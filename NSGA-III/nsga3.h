#ifndef _NSGA3_
#define _NSGA3_

#include "dominance.h"
#include "common.h"
#include "param.h"
#include "plot.h"

typedef struct {
  EMO_NDSort nd;
  EMO_List lst1, lst2; // temporary lists
  double *norm;        // normalized objective functions
  double *a;           // Intercepts of the hyper-plane
  double *one;         // Vector of ones
  double *axis;        // coordinate axes
  int wsize;           // number of weight vectors
  double *W;           // weight vectors
  double *waxis;       // weight vectors that are closest to the axes
  int *niche;          // niche count
  double *xtrm;        // matrix of extreme point
  double *inv;         // Inverse matrix
  double *min;         // ideal point
  int *pi;             // closest reference point for each individual
  double *dist;        // distance to the reference point
  int *filter;         // selected individuals
  double *cv;          // constraint value of individuals
  char *wfile;
} EMO_NSGA3;

void EMO_NSGA3_alloc(EMO_NSGA3 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_NSGA3_free(EMO_NSGA3 *alg);
void EMO_NSGA3_run(EMO_NSGA3 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);

#endif


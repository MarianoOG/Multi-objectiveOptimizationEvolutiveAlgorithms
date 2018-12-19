/* Author: Elmer Jesús Morales Orozco
 *         Universidad de la Sierra Sur,
 *         Miahuatlán, Oaxaca.
 *
 * Contribution: Mariano Orozco García
 *               UPIITA-IPN
 *               Mexico City
 *
 * August 2015
 */

#ifndef _MOEAD_H
#define _MOEAD_H

#include "utility.h"
#include "common.h"
#include "param.h"

typedef struct {
  int niche;       /* neighborhood size */
  int wsize;       /* number of weight vectors = population size */
  double *W;       /* weight vectors */
  double *ideal;   /* ideal point */
  double *diff;    /* temporary vector for subtraction */
  EMO_Utility utl; /* Utility function: tchebycheff, pbi, etc. */
  EMO_List *lnn;   /* array of lists for storing neighbors */
  double **sort;   /* temporary array for sorting */
  double *cv;      /* constraint violation */
  double s1;       /* parameter s1 */
  double s2;       /* parameter s2 */
  int ssize;
  char *wfile;
  char *utility;
} EMO_MOEAD;

void EMO_MOEAD_alloc(EMO_MOEAD *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_MOEAD_free(EMO_MOEAD *alg);
void EMO_MOEAD_run(EMO_MOEAD *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);

#endif


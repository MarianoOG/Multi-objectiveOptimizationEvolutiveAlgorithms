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

#include <stdlib.h>

#include "moead.h"
#include "numeric.h"
#include "vector.h"
#include "niche.h"
#include "evop.h"
#include "io.h"

#define _MAXCHAR 2000

/* Load specific parameters for the algorithm */
void EMO_MOEAD_load_param(EMO_MOEAD *alg, EMO_Param *param, int nvar) {

  if(!EMO_Parser_get_double(param->parser, &param->Pc, "pc")) {
    printf("Error, pc is not defined in the configuration file.\n");
    exit(1);
  }

  if(!EMO_Parser_get_double(param->parser, &param->Pm, "pm")) {
    printf("Error, pm is not defined in the configuration file.\n");
    exit(1);
  }

  if(!EMO_Parser_get_double(param->parser, &param->Nc, "nc")) {
    printf("Error, nc is not defined in the configuration file.\n");
    exit(1);
  }

  if(!EMO_Parser_get_double(param->parser, &param->Nm, "nm")) {
    printf("Error, nm is not defined in the configuration file.\n");
    exit(1);
  }

  param->Pm = (param->Pm == -1)? 1.0 / (double) nvar : param->Pm ;
  EMO_Debug_printf(param->dbg, "pm updated %f", param->Pm);

  if(!EMO_Parser_get_int(param->parser, &alg->niche, "moead_niche")) {
    printf("Error, moead_niche is not defined in the configuration file.\n");
    exit(1);
  }

  if(!EMO_Parser_get_char(param->parser, alg->wfile, "wfile")) {
    printf("Error, wfile is not defined in the configuration file.\n");
    exit(1);
  }

  if(!EMO_Parser_get_char(param->parser, alg->utility, "moead_utility")) {
    printf("Error, moead_utility is not defined in the configuration file.\n");
    exit(1);
  }

  // <MOG>
  if(!EMO_Parser_get_double(param->parser, &alg->s1, "moead_s1")) {
    printf("Error, moead_s1 is not defined in the configuration file.\n");
    exit(1);
  }

  if(!EMO_Parser_get_double(param->parser, &alg->s2, "moead_s2")) {
    printf("Error, moead_s2 is not defined in the configuration file.\n");
    exit(1);
  }
  // </MOG>

}

void EMO_MOEAD_alloc(EMO_MOEAD *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int i;

  printf("MOEA/D\n");

  if((alg->wfile = (char *) malloc(sizeof(char) * _MAXCHAR)) == NULL) {
    printf("Error, not enough memory in MOEA/D.\n");
    exit(1);
  }

  if((alg->utility = (char *) malloc(sizeof(char) * _MAXCHAR)) == NULL) {
    printf("Error, not enough memory in MOEA/D.\n");
    exit(1);
  }

  EMO_MOEAD_load_param(alg, param, mop->nvar);

  alg->wsize = 0;
  alg->W = EMO_File_read(NULL, &alg->wsize, &mop->nobj, alg->wfile, 0);

  if(alg->niche < 0 || alg->niche > alg->wsize) {
    printf("Error, niche must be less or equal than the number of weight vectors.\n");
    exit(1);
  }

  #ifdef EMO_MPI
  EMO_Population_alloc(pop, mop, alg->wsize, alg->wsize);
  #else
  EMO_Population_alloc(pop, mop, alg->wsize, 1);
  #endif

  if((alg->ideal = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in moead.\n");
    exit(1);
  }

  if((alg->diff = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in moead.\n");
    exit(1);
  }

  EMO_Utility_alloc(&alg->utl, param, mop->nobj, alg->utility);

  if((alg->lnn = (EMO_List *) malloc(sizeof(EMO_List) * alg->wsize)) == NULL) {
    printf("Error, not enogh memory in moead\n");
    exit(1);
  }

  for(i = 0; i < alg->wsize; i++)
    EMO_List_alloc(&alg->lnn[i], alg->niche);

  if((alg->sort = (double **) malloc(sizeof(double *) * alg->wsize)) == NULL) {
    printf("Error, not enough memory in moead.\n");
    exit(1);
  }

  for(i = 0; i < alg->wsize; i++) {
    if((alg->sort[i] = (double *) malloc(sizeof(double) * 2)) == NULL) {
      printf("Error, not enough memory in moead.\n");
      exit(1);
    }
  }

  alg->ssize = alg->wsize;

  // <MOG>
  if((alg->cv = (double *) malloc(sizeof(double *) * pop->size)) == NULL) {
    printf("Error, not enough memory in moead.\n");
    exit(1);
  }
  // </MOG>

  EMO_knn(alg->lnn, alg->sort, alg->W, alg->wsize, mop->nobj, alg->niche, 1);
}

void EMO_MOEAD_free(EMO_MOEAD *alg) {
  int i;

  free(alg->wfile);
  free(alg->utility);
  free(alg->W);
  free(alg->ideal);
  free(alg->diff);
  EMO_Utility_free(&alg->utl);

  for(i = alg->ssize - 1; i > -1; i--) {
    EMO_List_free(&alg->lnn[i]);
    free(alg->sort[i]);
  }

  free(alg->lnn);
  free(alg->sort);
}

void EMO_MOEAD_run(EMO_MOEAD *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int i, j, k, l, x, y, p1, p2;
  double v1, v2, tau, vmin, vmax;

  EMO_minBound(alg->ideal, pop->obj, NULL, pop->mu, mop->nobj);

/*  for(i = 0; i < mop->nobj; i++) {
    printf("ideal[] -> %f  ", alg->ideal[i]);
  } */

  x = pop->mu * mop->nvar;
  y = pop->mu * mop->nobj;

/*  printf("Population\n");
  for(i = 0; i < pop->mu; i++) {
    for(j = 0; j < mop->nobj; j++) {
      printf("%f ", pop->obj[i*mop->nobj + j]);
    }
  }*/

  /*printf("W\n");
  for(i = 0; i < pop->mu; i++) {
    printf("%d ->  ", i);
    for(j = 0; j < mop->nobj; j++) {
      printf("%f ", alg->W[i*mop->nobj + j]);
    }
    printf("\n");
    EMO_List_print(&alg->lnn[i], stdout, "");
  }*/

  // </MOG>
  if(mop->ncon != 0) {
    for(i = 0; i < pop->mu; i++) {
      alg->cv[i] = 0.0;

      if(pop->vio[i] != 0) {
        for(j = 0; j < mop->ncon; j++) {
          if(pop->con[i*mop->ncon + j] < 0.0)
            alg->cv[i] += -pop->con[i*mop->ncon + j];
        }
      }

      // printf("%d: cv: %f - vio: %d\n", i, alg->cv[i], pop->vio[i]);

      if(i==0){
        vmin = alg->cv[i];
        vmax = alg->cv[i];
      }
      else {
        vmin = min(alg->cv[i],vmin);
        vmax = max(alg->cv[i],vmax);
      }
    }
    // printf("vmin: %f, vmax: %f\n", vmin, vmax);
  }
  // </MOG>

  while(!EMO_Stop_end(param->stop)) {

    for (i = 0; i < pop->mu; i++) {
      p1 = EMO_Rand_int1(param->rand, 0, alg->niche - 1);
      p2 = EMO_Rand_int1(param->rand, 0, alg->niche - 1);
      // printf("p1: %d\n", p1);
      // printf("p2: %d\n\n", p2);
      EMO_List_get(&alg->lnn[i], &p1, p1);
      EMO_List_get(&alg->lnn[i], &p2, p2);
      // printf("P1: %d\n", p1);
      // printf("P2: %d\n\n", p2);

      p1 *= mop->nvar;
      p2 *= mop->nvar;

      EMO_crossSBX(pop->var+x, pop->vdummy, pop->var + p1, pop->var + p2, param->rand, mop, param->Pc, param->Nc);
      EMO_mutatePolynom(pop->var+x, param->rand, mop, param->Pm, param->Nm);
      //mop->f(mop, pop->obj + y, pop->var + x);
      EMO_Population_evaluate(pop, mop, pop->mu, 1);

      //PRINT EMO_vprint(stdout, pop->obj + y,mop->nobj, "new child");

      for(j=0; j< mop->nobj; j++) {
        if(pop->obj[y + j] < alg->ideal[j])
          alg->ideal[j] = pop->obj[y + j];
      }

      //PRINT EMO_vprint(stdout, alg->ideal, mop->nobj, "ideal");

      // EMO_List_print(&alg->lnn[i], NULL, "lnn");

      // <MOG>
      if(mop->ncon != 0) {
        alg->cv[pop->mu] = 0.0;

        if(pop->vio[pop->mu] != 0) {
          for(j = 0; j < mop->ncon; j++) {
            if(pop->con[pop->mu*mop->ncon + j] < 0.0)
              alg->cv[pop->mu] += -pop->con[pop->mu*mop->ncon + j];
          }
        }

        // printf("cv[pop->mu]: %f - vio[pop->mu]: %d\n", alg->cv[pop->mu], pop->vio[pop->mu]);

        vmin = min(alg->cv[pop->mu],vmin);
        vmax = max(alg->cv[pop->mu],vmax);

        // printf("vmin: %f, vmax: %f\n", vmin, vmax);

        tau = vmin + 0.3 * (vmax - vmin);
      }
      // </MOG>

      for(j = 0; j < alg->niche; j++) {
        // EMO_List_print(&alg->lnn[i], NULL, "lnn");
        EMO_List_get(&alg->lnn[i], &l, j);
        k = l * mop->nobj;
        // printf("%d, Vecino[%d]  %d\n", i, j, l);

        //EMO_vprint(stdout, &alg->W[k], mop->nobj, "w");

        EMO_vdiff(alg->diff, pop->obj+y , alg->ideal, mop->nobj);
        //v1 = alg->utl.uf(&alg->utl, &alg->W[k], pop->obj+y);
        v1 = alg->utl.uf(&alg->utl, &alg->W[k], alg->diff);

        EMO_vdiff(alg->diff, pop->obj+k , alg->ideal, mop->nobj);
        //v2 = alg->utl.uf(&alg->utl, &alg->W[k], pop->obj+k);
        v2 = alg->utl.uf(&alg->utl, &alg->W[k], alg->diff);

        // <MOG>
        if(mop->ncon != 0) {
          if(alg->cv[pop->mu] < tau)
            v1 += alg->s1 * alg->cv[pop->mu] * alg->cv[pop->mu];
          else
            v1 += alg->s1 * tau * tau + alg->s2 * (alg->cv[pop->mu] - tau);

          if(alg->cv[l] < tau)
            v2 += alg->s1 * alg->cv[l] * alg->cv[l];
          else
            v2 += alg->s1 * tau * tau + alg->s2 * (alg->cv[l] - tau);
        }
        // </MOG>

        //PRINT EMO_vprint(stdout, alg->W + k, mop->nobj, "w");
        //PRINT EMO_vprint(stdout, pop->obj+k, mop->nobj, "obj+k");
        //PRINT EMO_vprint(stdout, pop->obj+y, mop->nobj, "obj+y");
        //PRINT printf("v1 %f, v2 %f\n", v1, v2);

        //EMO_vprint(stdout, pop->obj + k, mop->nobj, "vecino");

        //printf("fit %d,%d: h %f vs l %f \n",i, j, v1, v2);

        if(v1 < v2) {
          if(mop->ncon == 0)
            EMO_Population_copy(pop, NULL, NULL, mop, l, pop->mu);
          else {
            EMO_Population_copy(pop, NULL, alg->cv, mop, l, pop->mu);

            vmin = alg->cv[0];
            vmax = alg->cv[0];
            for(j = 1; j < pop->mu; j++) {
              vmin = min(alg->cv[j],vmin);
              vmax = max(alg->cv[j],vmax);
            }
            // printf("----cambio----\n");
            // printf("vmin: %f, vmax: %f\n", vmin, vmax);
          }

        }
      }
    }

    EMO_Plot_run(param->plot, pop->obj, pop->mu, mop->feval, 0);
  }
}

#undef _MAXCHAR

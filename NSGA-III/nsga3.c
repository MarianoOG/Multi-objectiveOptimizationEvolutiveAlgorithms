#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include "nsga3.h"
#include "random.h"
#include "evop.h"
#include "utility.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "io.h"

#define _MAXCHAR 2000

/* Load specific parameters for the algorithm */
void NSGA3_load_param(EMO_NSGA3 *alg, EMO_Param *param, int nvar) {
  if(!EMO_Parser_get_int(param->parser, &param->mu, "psize")) {
    printf("Error, psize is not defined in the configuration file.\n");
    exit(1);
  }

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

  if(!EMO_Parser_get_char(param->parser, alg->wfile, "wfile")) {
    printf("Error, wfile is not defined in the configuration file.\n");
    exit(1);
  }
}


void vclose(double *x, double *vset, int *filter, int size, double *y, int n, double (*dist)(double *, double *, int)) {
  int i, j, k, l, m, flag;
  double vmin, d;

  m = 0;

  for(i = 0; i < n; i++) {
    vmin = DBL_MAX;

    l = i * n;

    for(j = 0; j < size; j++) {

      //RHG, primer frente if(filter == NULL || filter[j] == 0) {

      if(filter == NULL || filter[j]) {
        d = dist(vset + j*n, y + l, n);

        if(d < vmin) {
          flag = 0;

          /* Check for duplicated elements */
          for(k = 0; k < i; k++) {
            if(EMO_vdist(vset + j*n, x + k*n, n) < 1e-4) {
              flag = 1;
              break;
            }
          }

          if(!flag) {
            vmin = d;
            m = j;
          }
        }
      }
    }

    if(vmin != DBL_MAX) {
      memcpy(x + l, vset + m*n, sizeof(double) * n);

      /* Non negative entries */
      for(k = 0; k < n; k++) {
        if(x[l+k] < 0)
          x[l+k] = 0;
      }
    }
  }
}

// Selecciona puntos extremos de vset en x
// filter es rank
void vclose2(double *x, double *tmp, double *vset, int *filter, int size, double *y, int n) {
  int i, j, k, l, m = 0, flag; //o
  double vmin, d;
  //EMO_Utility utl;  // dtlz1, los primeros objetivos se hacen cero, al multiplicar por asf
  //EMO_UtilityFunction utlf;
  //EMO_mallocUtility(&utl, &utlf, "asf", n);

//printf("vclose\n");

  for(i = 0; i < n; i++) {
    vmin = DBL_MAX;

    l = i * n;
    //EMO_vprint(stdout, y+l, n, "axis");

    for(j = 0; j < size; j++) {

      if(filter == NULL || filter[j] == 0) {   // RHG version anterior solo aplicaba para los del primer frente
      //if(filter == NULL || filter[j]) {  // RHG version anterior solo aplicaba para los del primer frente

        //d = utlf(&utl, y + l, vset + j * n);
        //d = EMO_vdot(vset + j * n, y + l, n) / EMO_vnorm(vset + j*n, 2.0, n);
        //printf("cos theta %f\n", d);
        EMO_vorth(tmp, vset + j * n, y + l, n);
        d = EMO_vnorm(tmp, 2.0, n);

        if(d < vmin) {

          //printf("min %f\n", d);
          //EMO_vprint(stdout, vset+j*n, n, "vec");

          flag = 0;

          /* Check for duplicated or similar elements */
          for(k = 0; k < i; k++) {
            if(EMO_vdist(vset + j * n, x + k * n, n) < 1e-4) { //== 0  //< 1e-4
              flag = 1;
              break;
            }

            /* Check for duplicated components of the previous extreme points */
            /*for(o = 0; o < n; o++) {
              if(fabs(x[k*n + o]) > 1e-6 && x[k*n + o] == vset[j*n + o]) {
                flag = 1;
                break;
              }
            }*/
          }

/*          for(k = 0; k < n; k++) { dañino, no usar
            if(fabs(vset[j*n + k]) == 0.0) {
              flag = 1;
              break;
            }
          }*/


          /*for(k = 1; k < n; k++) {
            if(fabs(vset[j*n + k]) > 1e-6 && vset[j*n + k - 1] == vset[j*n + k]) {
              flag = 1;
              break;
            }
          }*/

          // Near to zero vector
          if(EMO_vnorm(vset + j * n, 2.0, n) < 1e-4) {
            flag = 1;
            break;
          }

          // Diagonal zero
          /*if(fabs(vset[j*n + i]) < 1e-4) {
            flag = 1;
            break;
          }*/

          if(!flag) {
            vmin = d;
            m = j;
            //printf("final vmin %f\n", d);
          }
        }
      }
    }

    //printf("%d: vmin %f, m %d\n", i, vmin, m);
    //if(vmin != DBL_MAX) EMO_vprint(stdout, vset + m*n, n,"v");

    if(vmin == DBL_MAX) {
      //if(gen < 1) {
        memset(x + l, 0, sizeof(double) * n);
        x[l + i] = 1.0;
        //EMO_mprint(stdout, "ncd", x + l, n, 1);
        printf("no hubo candidatos %d.\n", i);
      //}
    }
    else {
      memcpy(x + l, vset + m*n, sizeof(double) * n);
      //printf("vmin %f\n", vmin);
      /* Non negative entries */
/*      for(k = 0; k < n; k++) {
        if(x[l+k] < 0) {
          x[l+k] = 0;
          printf("negativo\n");
        }
      }*/
    }
  }

  //EMO_mprint(stdout, "xtrm", x, n, n);
  //printf(".----------------------\n");

  //EMO_freeUtility(&utl);
}

/* Normalization */
void normalize_nsga3(EMO_NSGA3 *alg, EMO_Population *pop, int nobj) {
  int i, j, k;

  /* Calculate the ideal point */
  EMO_minBound(alg->min, pop->obj, alg->filter, pop->size, nobj);
  //EMO_vprint(stdout, alg->min, nobj, "minbound");

  /* Translate objectives */
  for(i = 0; i < pop->size; i++) {
    if(alg->filter[i]) {
      for(j = 0; j < nobj; j++) {
        k = i * nobj + j;
        alg->norm[k] = pop->obj[k] - alg->min[j];
      }
    }
  }

  /* Find extreme points */ // waxis, axis ojo
    ///EMO_findmaxminBound(knn->max, knn->min, data, filter, size, knn->dim);

  vclose2(alg->xtrm, alg->a, alg->norm, alg->nd.rank, pop->size, alg->waxis, nobj);  // mejor distribucion, estable

  //for(i = 0; i < nobj; i++)
  //  EMO_vprint(stdout, alg->xtrm + i*nobj, nobj, "xtrm");

  if(EMO_minverse(alg->inv, alg->xtrm, nobj, 1, 0, 0) == 0) {

    //for(i = 0; i < nobj; i++)
    //  EMO_vprint(stdout, alg->inv + i*nobj, nobj, "inv");

    EMO_matmul(alg->a, alg->inv, alg->one, nobj, nobj, 1);
    //EMO_vprint(stdout, alg->a, nobj, "a");

    for(i = 0; i < pop->size; i++) {
      if(alg->filter[i]) {
        for(j = 0; j < nobj; j++)
          alg->norm[i*nobj + j] *= fabs(alg->a[j]);
      }
    }
  }
}

// Associate each individual w.r.t. a reference point, it also stores the distance to it.
// RHG pop se usa una vez
void associate_nsga3(EMO_NSGA3 *alg, EMO_Population *pop, int nobj) {
  int i, j, idx = 0;
  double d, vmin;

  for(i = 0; i < pop->size; i++) {
    vmin = DBL_MAX;

    if(!alg->filter[i]) continue;

    for(j = 0; j < alg->wsize; j++) {
      EMO_vorth(alg->a, alg->norm + i * nobj, alg->W + j * nobj, nobj);
      d = EMO_vnorm(alg->a, 2.0, nobj);

      if(d < vmin) {
        vmin = d;
        idx = j;
      }
    }
    alg->pi[i] = idx;
    alg->dist[i] = vmin;

    //printf("associate %d: %d, %f\n", i, alg->pi[i], alg->dist[i]);

  }
}

/* Selects n individuals */ //RHG, pop se usa una vez
void niching(EMO_NSGA3 *alg, EMO_Param *param, EMO_Population *pop, int nobj, int n, int lastfront) {
  int i, j, k, w, idx = 0;
  double vmin;

  //printf("por seleccionar %d, last front %d\n", n, lastfront);

  memset(alg->niche, 0, sizeof(int) * alg->wsize);

  for(i = 0; i < pop->size; i++)
    //if(alg.filter[i])
    if(alg->nd.rank[i] < lastfront)
      alg->niche[alg->pi[i]]++;


  // DEBUG printf("------------------------------------------------------ faltan %d\n", missing);
  //n = end_front - ini_front;
  //initList(&lst1, wsize);
  //initList(&lst2, n);
  EMO_List_clear(&alg->lst1);

  i = 0;

  while(i < n) {

    /* Look for the reference points that are isolated by the first fronts */
    if(alg->lst1.size == 0) {
      vmin = EMO_min(NULL, alg->niche, NULL, alg->wsize);

      for(j = 0; j < alg->wsize; j++)
        if(alg->niche[j] == vmin)
          EMO_List_queue(&alg->lst1, j);
    }

    /* Selects a random reference point */
    EMO_List_get(&alg->lst1, &w, EMO_Rand_int1(param->rand, 0, alg->lst1.size-1));

    // DEBUG printf("Selecciona punto de referencia: %d\n", w);

    EMO_List_clear(&alg->lst2);

    /* Look for the individuals in the last front that are near to the isolated reference points */
    for(j = 0; j < alg->nd.front[lastfront].size; j++) {
      EMO_List_get(&alg->nd.front[lastfront], &k, j);
      if(alg->filter[k] == 0 && alg->pi[k] == w)
        EMO_List_queue(&alg->lst2, k);
    }

    // DEBUG printf("Ind: ");
    // DEBUG lprint(alg.lst2);

    if(alg->lst2.size == 0) {  /* If there is no individual, the reference point is removed */
      alg->niche[w] = INT_MAX;
    }
    else {
      if(alg->niche[w] == 0) {   /* There is no individual in the first fronts that are associated */
        vmin = DBL_MAX;         /* with the reference point */

        for(j = 0; j < alg->lst2.size; j++) {
          EMO_List_get(&alg->lst2, &k, j);

          // DEBUG printf("ind %d, w %f\n", f, ind[f].dweight);

//          EMO_min(alg.dist[f],

          if(alg->dist[k] < vmin) {  /* Select the individual with the shortest distance */
            vmin = alg->dist[k];
            idx = k;
          }
        }
      }
      else {
        EMO_List_get(&alg->lst2, &idx, EMO_Rand_int1(param->rand, 0, alg->lst2.size-1));  /* Select a random individual */
        // DEBUG printf("rnd indf %d\n", idx);
      }

      // DEBUG printf("%d final %d\n", i, idx);
      alg->filter[idx] = 1;
      alg->niche[w]++;
      i++;
    }

    if(!EMO_List_remove(&alg->lst1, w)) {
      printf("Error al eliminar en lista.\n");
      exit(1);
    }
  }

  EMO_List_clear(&alg->lst1);
  EMO_List_clear(&alg->lst2);
}

// < MOG >

int EMO_NSGA3_calculate_cv(double *cv, double *con, int ncon, int *vio, int size) {
  int i, j, n;
  double *g;

  n = 0;

  for(i = 0; i < size; i++) {

    cv[i] = 0.0;

    if(vio[i] == 0) {
      n++;
    }
    else {
      g = con + i*ncon;
      for(j = 0; j < ncon; j++) {
        if(g[j] < 0.0)
          cv[i] += -g[j]; // Bracket operator
      }
    }

    // printf("%i: cv: %f - vio: %d\n", i, cv[i], vio[i]);
  }

  return n;
}

int EMO_Dominance_constraint3(double *x, double *y, int n, double cvx, double cvy, EMO_Dominance r) {

  if(cvx == 0.0 && cvy == 0.0)
    return r(x,y,n);
  else if (cvx == 0.0 && cvy != 0.0)
    return 1;
  else if (cvx != 0.0 && cvy == 0.0)
    return -1;
  else {
    if (cvx < cvy)
      return 1;
    else if (cvx > cvy)
      return -1;
    else
      return 2; // No constraint-value dominance
  }

  return 0;
}

/* Routine to perform non-dominated sorting of NSGA-III for handling constraints */
void EMO_NDSort_run3(EMO_NDSort *nd, double *obj, int nobj, double *cv, int *filter, int size) {
  int i, j, v, f;

  EMO_List_clear(&nd->front[0]);

  if(filter == NULL) {
    for(i = 0; i < size; i++) {
      nd->n[i] = 0;

      for(j = 0; j < size; j++) {
        if(i == j) continue;

        v = EMO_Dominance_constraint3(obj+i*nobj, obj+j*nobj, nobj,
                                     cv[i], cv[j], EMO_Dominance_strict);

        if(v == 1) EMO_List_queue(&nd->S[i] , j);
        if(v == -1) nd->n[i]++;
      }

      if(nd->n[i] == 0) {
        nd->rank[i] = 0;
        EMO_List_queue(&nd->front[0], i);
      }
      //pop->ind[i].d = pop->ind[i].n;
    }
  }
  else {
    for(i = 0; i < size; i++) {

      nd->n[i] = 0;

      if(!filter[i]) continue;

      for(j = 0; j < size; j++) {
        if(i == j || !filter[i]) continue;

        v = EMO_Dominance_constraint3(obj+i*nobj, obj+j*nobj, nobj,
                                     cv[i], cv[j], EMO_Dominance_strict);

        if(v == 1) EMO_List_queue(&nd->S[i] , j);
        if(v == -1) nd->n[i]++;
      }

      if(nd->n[i] == 0) {
        nd->rank[i] = 0;
        EMO_List_queue(&nd->front[0], i);
      }
      //pop->ind[i].d = pop->ind[i].n;
    }
  }

  f = 0;

  while(nd->front[f].size != 0) {
    EMO_List_clear(&nd->front[f+1]);

    for(i = 0; i < nd->front[f].size; i++) {
      EMO_List_get(&nd->front[f], &v, i);

      while(nd->S[v].size != 0) {
        EMO_List_dequeue(&nd->S[v], &j);
        nd->n[j]--;

        if(nd->n[j] == 0) {
          nd->rank[j] = f + 1;
          EMO_List_queue(&nd->front[f+1], j);
        }
      }
    }
    f++;
    //lprint(front[f]);
  }

  nd->nfront = f;
}

int tournament_selection(EMO_Rand *rnd, int mu, int *vio, double *cv) {
  int p1, p2;

  p1 = EMO_Rand_int1(rnd, 0, mu-1);
  p2 = EMO_Rand_int1(rnd, 0, mu-1);

  // printf("tournament: %d %d\n", p1, p2);
  // printf("vio: %d %d\n", vio[p1], vio[p2]);
  // printf("cv: %f %f\n", cv[p1], cv[p2]);

  if (vio[p1] == 0 && vio[p2] != 0)
    return p1;
  else if (vio[p1] != 0 && vio[p2] == 0)
    return p2;
  else if (vio[p1] != 0 && vio[p2] != 0) {
    if (cv[p1] > cv[p2])
      return p2;
    else if (cv[p1] < cv[p2])
      return p1;
    else {
      if (EMO_Rand_flip(rnd, 0.5))
        return p1;
      else
        return p2;
    }
  }
  else { // if (vio[p1] == 0 && vio[p2] == 0) {
    if (EMO_Rand_flip(rnd, 0.5))
      return p1;
    else
      return p2;
  }

  return 0;
}

// < MOG />

void EMO_NSGA3_alloc(EMO_NSGA3 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int i, n;

  printf("NSGA3\n");

  if((alg->wfile = (char *) malloc(sizeof(char) * _MAXCHAR)) == NULL) {
    printf("Error, not enough memory in NSGA3.\n");
    exit(1);
  }

  NSGA3_load_param(alg, param, mop->nvar);

  EMO_Population_alloc(pop, mop, param->mu, param->mu);
  EMO_NDSort_alloc(&alg->nd, pop->size);
  EMO_List_alloc(&alg->lst1, pop->mu);
  EMO_List_alloc(&alg->lst2, pop->mu);

  if((alg->norm = (double *) calloc(sizeof(double), pop->size * mop->nobj)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  if((alg->a = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  if((alg->one = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  for(i = 0; i < mop->nobj; i++)
    alg->one[i] = 1.0;

  n = mop->nobj * mop->nobj;

  if((alg->axis = (double *) malloc(sizeof(double) * n)) == NULL) {  // identity
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  EMO_vaxes(alg->axis, mop->nobj);

  alg->wsize = 0;
  alg->W = EMO_File_read(NULL, &alg->wsize, &mop->nobj, alg->wfile, 0);

  if((alg->waxis = (double *) malloc(sizeof(double) * n)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  vclose(alg->waxis, alg->W, NULL, alg->wsize, alg->axis, mop->nobj, EMO_vdist);

  //EMO_mprint(stdout, "waxis", alg->waxis, mop->nobj, mop->nobj);
  //exit(1);

  if((alg->niche = (int *) malloc(sizeof(int) * alg->wsize)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  if((alg->xtrm = (double *) malloc(sizeof(double) * n)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  if((alg->inv = (double *) malloc(sizeof(double) * n)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  if((alg->min = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  if((alg->pi = (int *) calloc(sizeof(int), pop->size)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  if((alg->dist = (double *) calloc(sizeof(double), pop->size)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  if((alg->filter = (int *) calloc(sizeof(int), pop->size)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  if((alg->cv = (double *) calloc(sizeof(double), pop->size)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }
}

void EMO_NSGA3_free(EMO_NSGA3 *alg) {
  free(alg->wfile);
  EMO_NDSort_free(&alg->nd);
  EMO_List_free(&alg->lst1);
  EMO_List_free(&alg->lst2);
  free(alg->norm);
  free(alg->a);
  free(alg->one);
  free(alg->axis);
  free(alg->W);
  free(alg->waxis);
  free(alg->niche);
  free(alg->xtrm);
  free(alg->inv);
  free(alg->min);
  free(alg->pi);
  free(alg->dist);
  free(alg->filter);
}

void EMO_NSGA3_run(EMO_NSGA3 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int  p1, p2, v1, v2, i, j, nfeasible, cont, lastfront, lfcont, *shuffle;  //o1, o2,

  if(mop->ncon != 0)
    EMO_NSGA3_calculate_cv(alg->cv, pop->con, mop->ncon, pop->vio, pop->mu);

  while(!EMO_Stop_end(param->stop)) {

    /* Offspring generation */
    for(i = 0; i < pop->lambda; i+=2) {

      /* Select two parents  */
      if(mop->ncon == 0) {
        p1 = EMO_Rand_int1(param->rand, 0, pop->mu-1);
        while ((p2 = EMO_Rand_int1(param->rand, 0, pop->mu-1)) == p1);
      }
      else {
        p1 = tournament_selection(param->rand, pop->mu, pop->vio, alg->cv);
        while ((p2 = tournament_selection(param->rand, pop->mu, pop->vio, alg->cv)) == p1);
      }

      // printf("%d: parents %d %d\n", i, p1, p2);

      p1 *= mop->nvar;
      p2 *= mop->nvar;

      v1 = pop->mu + i;  // o1 =
      v2 = pop->mu + i + 1;  // o2 =
      v1 *= mop->nvar;
      v2 *= mop->nvar;
      //o1 *= mop->nobj;
      //o2 *= mop->nobj;

      /* Generate an offspring by variation operators */
      EMO_crossSBX(pop->var+v1, pop->var+v2, pop->var+p1, pop->var+p2, param->rand, mop, param->Pc, param->Nc);
      EMO_mutatePolynom(pop->var+v1, param->rand, mop, param->Pm, param->Nm);
      EMO_mutatePolynom(pop->var+v2, param->rand, mop, param->Pm, param->Nm);
      //mop->f(mop, pop->obj+o1, pop->var+v1);
      //mop->f(mop, pop->obj+o2, pop->var+v2);
      EMO_Population_evaluate(pop, mop, pop->mu + i, 2);
    }

    /* Reduce population */
    if(mop->ncon == 0) {
      EMO_NDSort_run(&alg->nd, pop->obj, mop->nobj, NULL, pop->size);  /* Non-dominated sorting algorithm */

      //printf("size %d, front %d\n", pop->size, alg->nd.nfront);

      j = cont = 0;

      do {
        cont += alg->nd.front[j++].size;
      } while(cont < pop->mu);

      lastfront = j - 1;

      //printf("cont %d, lastfront %d\n", cont, lastfront);

      /* Select solutions from the first fronts */
      memset(alg->filter, 0, sizeof(int) * pop->size);

      for(i = 0; i < pop->size; i++) {
        if(alg->nd.rank[i] <= lastfront) {
          alg->filter[i] = 1;
        }
      }

      if(cont != pop->lambda) {
        normalize_nsga3(alg, pop, mop->nobj);
        associate_nsga3(alg, pop, mop->nobj);

        // Unselect last front
        for(i = 0; i < alg->nd.front[lastfront].size; i++) {
          EMO_List_get(&alg->nd.front[lastfront], &j, i);
          alg->filter[j] = 0;
        }

        niching(alg, param, pop, mop->nobj, pop->mu - cont + alg->nd.front[lastfront].size, lastfront);
      }

    }
    else {
      nfeasible = EMO_NSGA3_calculate_cv(alg->cv, pop->con, mop->ncon, pop->vio, pop->size);
      EMO_NDSort_run3(&alg->nd, pop->obj, mop->nobj, alg->cv, NULL, pop->size);

      // Printing feasibles
      printf("Feasible: %d/%d, nfront: %d\n", nfeasible, pop->size, alg->nd.nfront);
      /* for(i = 0; i < pop->size; i++) {
        printf("i: %d,\tvio: %d,\tcv: %f,\trank: %d\n",
          i, pop->vio[i], alg->cv[i], alg->nd.rank[i]);
      }*/

      if (nfeasible >= pop->mu) {

        j = cont = 0;

        do {
          cont += alg->nd.front[j++].size;
        } while(cont < pop->mu);

        lastfront = j - 1;

        //printf("cont %d, lastfront %d\n", cont, lastfront);

        /* Select solutions from the first fronts */
        memset(alg->filter, 0, sizeof(int) * pop->size);

        for(i = 0; i < pop->size; i++) {
          if(alg->nd.rank[i] <= lastfront) {
            alg->filter[i] = 1;
          }
        }

        if(cont != pop->lambda) {
          normalize_nsga3(alg, pop, mop->nobj);
          associate_nsga3(alg, pop, mop->nobj);

          // EMO_vprint(stdout, alg->min, mop->nobj, "min: ");
          // EMO_vprint(stdout, alg->xtrm, mop->nobj, "xtrm: ");

          // Unselect last front
          for(i = 0; i < alg->nd.front[lastfront].size; i++) {
            EMO_List_get(&alg->nd.front[lastfront], &j, i);
            alg->filter[j] = 0;
          }

          niching(alg, param, pop, mop->nobj, pop->mu - cont + alg->nd.front[lastfront].size, lastfront);
        }

      } else {

        j = cont = 0;

        do {
          cont += alg->nd.front[j++].size;
        } while(cont < pop->mu);

        lastfront = j - 1;

        /* Select solutions from the first fronts */
        memset(alg->filter, 0, sizeof(int) * pop->size);

        if((shuffle = (int *) malloc(sizeof(int) * alg->nd.front[lastfront].size)) == NULL) {
          printf("Error, not enough memory in NSGA3.\n");
          exit(1);
        }

        for(i = 0; i < alg->nd.front[lastfront].size; i++) {
          if(cont > pop->mu) {
            shuffle[i] = 0;
            cont--;
          }
          else
            shuffle[i] = 1;
        }

        /* printf("shuffle: ");
        for(i = 0; i < alg->nd.front[lastfront].size; i++) {
          printf("%d, ", shuffle[i]);
        }
        printf("\n"); */

        EMO_Rand_shuffle(param->rand, shuffle, alg->nd.front[lastfront].size);

        /* printf("shuffle: ");
        for(i = 0; i < alg->nd.front[lastfront].size; i++) {
          printf("%d, ", shuffle[i]);
        }
        printf("\n"); */

        lfcont = 0;

        for(i = 0; i < pop->size; i++) {
          if(alg->nd.rank[i] <= lastfront)
            alg->filter[i] = 1;
          if(alg->nd.rank[i] == lastfront) {
            alg->filter[i] = shuffle[lfcont];
            lfcont++;
          }
        }

        // printf("lfcont: %d\n", lfcont);
        // printf("%d == %d\n", cont, pop->mu);
      }
    }

    /* Reduce population */
    EMO_Population_survive(pop, NULL, NULL, mop, &alg->lst1, &alg->lst2, alg->filter);
    EMO_Plot_run(param->plot, pop->obj, pop->mu, mop->feval, 0);

  }
}

#undef _MAXCHAR
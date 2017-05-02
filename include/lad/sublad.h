//
// Created by M. Engler on 07/12/16.
//

#ifndef SUBLAD_H
#define SUBLAD_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/time.h>
#include <sys/resource.h>

#ifndef __cplusplus

int* matchedWithV; // matchedWithV[matchedWithU[u]]=u
int* nbPred; // nbPred[i] = nb of predecessors of the ith vertex of V in the DAG
int** pred; // pred[i][j] = jth predecessor the ith vertex of V in the DAG
int* nbSucc; // nbSucc[i] = nb of successors of the ith vertex of U in the DAG
int** succ; // succ[i][j] = jth successor of the ith vertex of U in the DAG
int* listV;
int* listU;
int* listDV;
int* listDU;

#endif

// define boolean type as char
#define c_true 1
#define c_false 0
#define c_bool char

typedef struct{
  c_bool isDirected; // false iff for each edge (i,j), there exists an edge (j,i)
  c_bool isLabelled; // true if labels are associetd with vertices and edges
  c_bool* isLoop; // isLoop[i] = true if there is a loop on vertex i
  int nbVertices; // Number of vertices
  int* vertexLabel; // if isLabelled then vertexLabel[i] = label associated with vertex i
  int* nbAdj;    // nbAdj[i] = number of vertices j such that (i,j) or (j,i) is an edge
  int* nbPred;   // nbPred[i] = number of vertices j such that (j,i) is an edge and (i,j) is not an edge
  int* nbSucc;   // nbSucc[i] = number of vertices j such that (i,j) is an edge and (j,i) is not an edge
  int** adj;     // forall j in [0..nbAdj[i]-1], adj[i][j] = jth vertex adjacent to i
  char** edgeDirection;	// if both (i,j) and (j,i) are edges then edgeDirection[i][j] = 3
  // else if (i,j) is an edge then edgeDirection[i][j] = 1
  // else if (j,i) is an edge then edgeDirection[i][j] = 2
  // else (neither (i,j) nor (j,i) is an edge) edgeDirection[i][j] = 0
  int** edgeLabel; // if isLabelled then edgeLabel[i][j] = label associated with edge (i,j)
} Tgraph;

typedef struct{
  int *nbVal;    // nbVal[u] = number of values in D[u]
  int *firstVal; // firstVal[u] = pos in val of the first value of D[u]
  int *val;      // val[firstVal[u]..firstVal[u]+nbVal[u]-1] = values of D[u]
  int **posInVal;
  // If v in D[u] then firstVal[u] <= posInVal[u][v] < firstVal[u]+nbVal[u]
  //                   and val[posInVal[u][v]] = v
  // otherwise posInVal[u][v] >= firstVal[u]+nbVal[u]
  int **firstMatch; // firstMatch[u][v] = pos in match of the first vertex of the covering matching of G_(u,v)
  int *matching; // matching[firstMatch[u][v]..firstMatch[u][v]+nbAdj[u]-1] = covering matching of G_(u,v)
  int nextOutToFilter; // position in toFilter of the next pattern node whose domain should be filtered (-1 if no domain to filter)
  int lastInToFilter; // position in toFilter of the last pattern node whose domain should be filtered
  int *toFilter;  // contain all pattern nodes whose domain should be filtered
  c_bool *markedToFilter;    // markedToFilter[u]=c_true if u is in toFilter; c_false otherwise
  int* globalMatchingP; // globalMatchingP[u] = node of Gt matched to u in globalAllDiff(Np)
  int* globalMatchingT; // globalMatchingT[v] = node of Gp matched to v in globalAllDiff(Np) or -1 if v is not matched
} Tdomain;

#ifdef __cplusplus
extern "C" {
#endif

extern int sub_iso(Tgraph *Gp, Tgraph *Gt, c_bool *iso, int *map, int verbose, int timeLimit, c_bool induced, c_bool firstSol);

#ifdef __cplusplus
}
#endif

#endif //SUBLAD_H

//
// Created by M. Engler on 06/12/16.
//

#include <sublad.h>
#include "subgraph_isomorphism.h"

Tgraph* mogli::translate_graph(const Molecule &mol) {
  // reads data in fileName and create the corresponding graph

  int i, j;
  Tgraph* graph = (Tgraph*)malloc(sizeof(Tgraph));

  graph->nbVertices = mol.get_atom_count();

  graph->isLabelled = c_true;
  graph->vertexLabel = (int*)calloc(graph->nbVertices,sizeof(int));
  memset(graph->vertexLabel,0,graph->nbVertices*sizeof(int));
  graph->edgeLabel = (int**)calloc(graph->nbVertices,sizeof(int*));

  graph->isLoop = (c_bool*)calloc(graph->nbVertices,sizeof(c_bool));
  graph->nbAdj = (int*)calloc(graph->nbVertices,sizeof(int));
  graph->nbPred = (int*)calloc(graph->nbVertices,sizeof(int));
  graph->nbSucc = (int*)calloc(graph->nbVertices,sizeof(int));
  graph->edgeDirection = (char**)malloc(graph->nbVertices*sizeof(char*));
  graph->adj = (int**)malloc(graph->nbVertices*sizeof(int*));
  for (i=0; i<graph->nbVertices; i++){
    graph->isLoop[i] = c_false;
    graph->adj[i] = (int*)malloc(graph->nbVertices*sizeof(int));
    graph->edgeDirection[i] = (char*)malloc(graph->nbVertices*sizeof(char));
    memset(graph->edgeDirection[i],0,graph->nbVertices*sizeof(char));
    memset(graph->nbAdj,0,graph->nbVertices*sizeof(int));
    memset(graph->nbPred,0,graph->nbVertices*sizeof(int));

    graph->edgeLabel[i] = (int*)malloc(graph->nbVertices*sizeof(int));
    memset(graph->edgeLabel[i],0,graph->nbVertices*sizeof(int));

  }

  Graph::NodeMap<int> node_ids(mol.get_graph());
  i = 0;
  for (NodeIt u = mol.get_node_iter(); u != lemon::INVALID; ++u, ++i) {
    node_ids[u] = i;
  }

  for (NodeIt u = mol.get_node_iter(); u != lemon::INVALID; ++u) {
    const int _u = node_ids[u];
    graph->vertexLabel[_u] = mol.get_color(u);
    int degree = 0;
    graph->isLoop[_u] = c_false;
    for (IncEdgeIt e = mol.get_inc_edge_iter(u); e != lemon::INVALID; ++e) {
      Node v = mol.get_opposite_node(u, e);
      const int _v = node_ids[v];
      graph->edgeLabel[_u][_v] = 0;
      graph->nbPred[_v]++;
      assert(graph->nbAdj[_u] < graph->nbVertices);
      graph->adj[_u][graph->nbAdj[_u]++] = _v;
      graph->edgeDirection[_u][_v] = 3;
      ++degree;
    }
    graph->nbSucc[_u] = degree;
  }

  graph->isDirected = c_false;
  for (i=0; i<graph->nbVertices; i++){
    for (j=0; j<graph->nbAdj[i]; j++){
      graph->isDirected = graph->isDirected || (graph->edgeDirection[i][graph->adj[i][j]]!=3);
    }
  }
  return graph;
}

void mogli::free_graph(Tgraph *graph) {

  for (int i=0; i<graph->nbVertices; i++){
    free(graph->adj[i]);
    free(graph->edgeDirection[i]);
    free(graph->edgeLabel[i]);
  }

  free(graph->vertexLabel);
  free(graph->edgeLabel);
  free(graph->isLoop);
  free(graph->nbAdj);
  free(graph->nbPred);
  free(graph->nbSucc);
  free(graph->edgeDirection);
  free(graph->adj);

  free(graph);
}

bool mogli::are_subgraph_isomorphic(const Molecule &mol_small, const Molecule &mol_large, int map[]) {
  Tgraph* gp = translate_graph(mol_small);
  Tgraph* gt = translate_graph(mol_large);

  assert(!gp->isDirected);
  assert(!gt->isDirected);

  c_bool iso = c_false;
  sub_iso(gp, gt, &iso, map, 0, 60, c_true, c_true);

  free_graph(gp);
  free_graph(gt);

  return iso == c_true;
}

bool mogli::are_subgraph_isomorphic(Tgraph* graph_small, Tgraph* graph_large, int map[]) {
  c_bool iso = c_false;
  sub_iso(graph_small, graph_large, &iso, map, 0, 60, c_true, c_true);
  return iso == c_true;
}
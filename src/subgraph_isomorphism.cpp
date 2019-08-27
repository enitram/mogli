////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    mogli - molecular graph library                                                                                 //
//                                                                                                                    //
//    Copyright (C) 2016-2019  Martin S. Engler                                                                       //
//                                                                                                                    //
//    This program is free software: you can redistribute it and/or modify                                            //
//    it under the terms of the GNU Lesser General Public License as published                                        //
//    by the Free Software Foundation, either version 3 of the License, or                                            //
//    (at your option) any later version.                                                                             //
//                                                                                                                    //
//    This program is distributed in the hope that it will be useful,                                                 //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of                                                  //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                                    //
//    GNU General Public License for more details.                                                                    //
//                                                                                                                    //
//    You should have received a copy of the GNU Lesser General Public License                                        //
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.                                          //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "subgraph_isomorphism.h"

#include <assert.h>
#include <malloc.h>


Tgraph* mogli::translate_graph(const Molecule &mol, IntVector &node_ids) {
  // reads data in fileName and create the corresponding graph

  int i, j;
  auto* graph = (Tgraph*)malloc(sizeof(Tgraph));

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

  IntToIntMap reverse_node_ids;
  i = 0;
  for (NodeIt u = mol.get_node_iter(); u != lemon::INVALID; ++u, ++i) {
    reverse_node_ids[mol.get_id(u)] = i;
    node_ids.push_back(mol.get_id(u));
  }
  assert(node_ids.size() == mol.get_atom_count());

  auto table = mol.get_perdiodic_table();
  for (NodeIt u = mol.get_node_iter(); u != lemon::INVALID; ++u) {
    const int _u = reverse_node_ids[mol.get_id(u)];
    graph->vertexLabel[_u] = table.get_equivalency_class(mol.get_color(u));
    int degree = 0;
    graph->isLoop[_u] = c_false;
    for (IncEdgeIt e = mol.get_inc_edge_iter(u); e != lemon::INVALID; ++e) {
      Node v = mol.get_opposite_node(u, e);
      const int _v = reverse_node_ids[mol.get_id(v)];
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

void mogli::translate_maps(const IntVector &node_ids_small, const IntVector &node_ids_large,
                           const int in_iso_map[], IntToIntMap &out_iso_map) {
  for (int i = 0; i < node_ids_small.size(); ++i) {
    int _to = in_iso_map[i];
    if (_to > -1) {
      assert(_to < node_ids_large.size());
      int to = node_ids_large.at(_to);
      int from = node_ids_small.at(i);
      out_iso_map[from] = to;
    }
  }
}

bool mogli::are_subgraph_isomorphic(const Molecule &mol_small, const Molecule &mol_large, IntToIntMap isomorphism_map) {
  IntVector node_ids_small, node_ids_large;
  Tgraph* gp = translate_graph(mol_small, node_ids_small);
  Tgraph* gt = translate_graph(mol_large, node_ids_large);

  assert(gp->isDirected == c_false);
  assert(gt->isDirected == c_false);

  assert(gp->nbVertices > 0);
  assert(gt->nbVertices > 0);

  c_bool iso = c_false;
  int n = mol_large.get_atom_count();
  int map[n];
  std::fill_n(map, n, -1);
  sub_iso(gp, gt, &iso, map, 0, 60, c_true, c_true);

  free_graph(gp);
  free_graph(gt);

  if (iso == c_true) {
    translate_maps(node_ids_small, node_ids_large, map, isomorphism_map);
    return true;
  } else {
    return false;
  }
}

bool mogli::are_subgraph_isomorphic(Tgraph* graph_small, Tgraph* graph_large, int isomorphism_map[]) {
  assert(graph_large->isDirected == c_false);
  assert(graph_small->isDirected == c_false);

  assert(graph_large->nbVertices > 0);
  assert(graph_small->nbVertices > 0);

  c_bool iso = c_false;
  sub_iso(graph_small, graph_large, &iso, isomorphism_map, 0, 60, c_true, c_true);
  return iso == c_true;
}
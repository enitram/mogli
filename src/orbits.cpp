//
// Created by M. Engler on 02/12/16.
//

#include "orbits.h"

using namespace mogli;

void Orbits::init(const Molecule &mol) {

  Node v = mol.get_node_iter();
  assert(v != lemon::INVALID);

  NodeToBoolMap visited(mol.get_graph(), false);
  ShortToNodeVectorMap colorMap;
  ShortSet colorSet;

  dfs(v, mol, visited, colorSet, colorMap);

}

void Orbits::dfs(const Node& current, const Molecule& mol, NodeToBoolMap& visited, ShortSet& colorSet,
                 ShortToNodeVectorMap& colorMap) {
  visited[current] = true;
  unsigned short color = mol.get_color(current);
  colorSet.insert(color);
  if (colorMap.find(color) == colorMap.end()) {
    NodeVector vector;
    colorMap.insert(std::pair<int, NodeVector>(color, vector));
  }
  colorMap.at(color).push_back(current);

  for (IncEdgeIt e = mol.get_inc_edge_iter(current); e != lemon::INVALID; ++e) {
    Node w = mol.get_opposite_node(current, e);
    if (!visited[w]) {
      dfs(w, mol, visited, colorSet, colorMap);
    }
  }

}

void Orbits::orbitsNauty(const Molecule &mol, const ShortSet &colorSet, const ShortToNodeVectorMap &colorMap,
                    const unsigned int atom_count) {
  DYNALLSTAT(int,lab,lab_sz);
  DYNALLSTAT(int,ptn,ptn_sz);
  DYNALLSTAT(int,orbits,orbits_sz);
  DYNALLSTAT(graph,ng,g_sz);
  DYNALLSTAT(graph,cg,cg_sz);

  static DEFAULTOPTIONS_GRAPH(options);
  statsblk stats;

  options.getcanon = FALSE;
  options.defaultptn = FALSE;

  int m = SETWORDSNEEDED(atom_count);
  nauty_check(WORDSIZE,m,atom_count,NAUTYVERSIONID);

  DYNALLOC1(int,lab,lab_sz,atom_count,"malloc");
  DYNALLOC1(int,ptn,ptn_sz,atom_count,"malloc");
  DYNALLOC1(int,orbits,orbits_sz,atom_count,"malloc");

  DYNALLOC2(graph,ng,g_sz,atom_count,m,"malloc");
  EMPTYGRAPH(ng,m,atom_count);

  int i = 0;
  NodeToIntMap nodes(mol.get_graph());
  NodeVector first_order;

  for (ShortSet::iterator it=colorSet.begin(), end = colorSet.end(); it != end; ++it) {
    NodeVector vector = colorMap.at(*it);
    for (NodeVector::iterator it2 = vector.begin(), end2 = vector.end(); it2 != end2; ++it2) {
      first_order.push_back(*it2);
      nodes[*it2] = i;
      lab[i] = i;
      ptn[i] = 1;
      ++i;
    }
    ptn[i-1] = 0;
  }

  for(EdgeIt e = mol.get_edge_iter(); e!=lemon::INVALID; ++e) {
    int u = nodes[mol.get_u(e)];
    int v = nodes[mol.get_v(e)];
    ADDONEEDGE(ng,u,v,m);
  }

  DYNALLOC2(graph,cg,cg_sz,atom_count,m,"malloc");
  densenauty(ng,lab,ptn,orbits,&options,&stats,m,atom_count,cg);

  for (i = 0; i < atom_count; ++i) {
    _orbits[first_order[i]] = orbits[i];
  }

}

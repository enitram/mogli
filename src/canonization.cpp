//
// Created by martin on 10/20/16.
//

#include "canonization.h"

void mogli::Canonization::init(const Molecule &mol) {

  const Graph& g = mol.getGraph();

  NodeToBoolMap visited(g, false);
  ShortToNodeVectorMap colorMap;
  NodeIt v(g);
  bool is_tree = true;

  assert(v != lemon::INVALID);

  dfs(v, v, mol, g, visited, colorMap, is_tree);

  if (!is_tree) {
    //canonTree(g, colorMap);
    canonNauty(g, colorMap, mol.getAtomCount());
  } else {
    canonNauty(g, colorMap, mol.getAtomCount());
  }
}

void mogli::Canonization::dfs(const Node& current, const Node& last,
                              const Molecule& mol, const Graph& g, NodeToBoolMap& visited,
                              ShortToNodeVectorMap& colorMap, bool& is_tree) {
  visited[current] = true;
  unsigned short color = mol.getColor(current);
  _colors.insert(color);
  if (colorMap.find(color) == colorMap.end()) {
    NodeVector vector;
    colorMap.insert(std::pair<int, NodeVector>(color, vector));
  }
  colorMap.at(color).push_back(current);

  for (IncEdgeIt e(g, current); e != lemon::INVALID; ++e) {
    Node w = g.oppositeNode(current, e);
    if (!visited[w]) {
      dfs(w, current, mol, g, visited, colorMap, is_tree);
    } else if (w != last) {
      is_tree = false;
    }
  }

}

void mogli::Canonization::canonNauty(const Graph& g, const ShortToNodeVectorMap &colorMap, const unsigned int atom_count) {

  DYNALLSTAT(int,lab,lab_sz);
  DYNALLSTAT(int,ptn,ptn_sz);
  DYNALLSTAT(int,orbits,orbits_sz);
  DYNALLSTAT(graph,ng,g_sz);
  DYNALLSTAT(graph,cg,cg_sz);

  static DEFAULTOPTIONS_GRAPH(options);
  statsblk stats;

  options.getcanon = TRUE;
  options.defaultptn = FALSE;

  int m = SETWORDSNEEDED(atom_count);
  nauty_check(WORDSIZE,m,atom_count,NAUTYVERSIONID);

  DYNALLOC1(int,lab,lab_sz,atom_count,"malloc");
  DYNALLOC1(int,ptn,ptn_sz,atom_count,"malloc");
  DYNALLOC1(int,orbits,orbits_sz,atom_count,"malloc");

  DYNALLOC2(graph,ng,g_sz,atom_count,m,"malloc");
  EMPTYGRAPH(ng,m,atom_count);

  int i = 0;
  for (ShortSet::iterator it=_colors.begin(), end = _colors.end(); it != end; ++it) {
    NodeVector vector = colorMap.at(*it);
    for (NodeVector::iterator it2 = vector.begin(), end2 = vector.end(); it2 != end2; ++it2) {
      lab[i] = g.id(*it2);
      ptn[i] = 1;
      ++i;
    }
    ptn[i-1] = 0;
  }

  for(Graph::EdgeIt e(g); e!=lemon::INVALID; ++e) {
    int u = g.id(g.u(e));
    int v = g.id(g.v(e));
    ADDONEEDGE(ng,u,v,m);
  }

  DYNALLOC2(graph,cg,cg_sz,atom_count,m,"malloc");
  densenauty(ng,lab,ptn,orbits,&options,&stats,m,atom_count,cg);

//  std::cout << "lab" << std::endl;
//  for (i = 0; i < atom_count; ++i) {
//    std::cout << lab[i] << std::endl;
//  }

  // TODO in which order does nauty return the canonical graphs?
  // TODO if our tree alg has a different order, could we accidentally return the same canonization?
//  std::cout << "cg" << std::endl;
//  for (i = 0; i < m*atom_count; ++i) {
//    std::cout << std::bitset<WORDSIZE>(static_cast<unsigned long>(cg[i])) << std::endl;
//  }

  for (i = 0; i < m*atom_count; ++i) {
    _canonization.push_back(static_cast<unsigned long>(cg[i]));
  }

}

void mogli::Canonization::canonTree(const Graph &g, const ShortToNodeVectorMap &colorMap) {
  // TODO write own canonization method
}

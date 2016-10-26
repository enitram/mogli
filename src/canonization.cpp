//
// Created by martin on 10/20/16.
//

#include "canonization.h"

void mogli::Canonization::init(const Molecule &mol) {

  Node v = mol.get_node_iter();
  assert(v != lemon::INVALID);

  NodeToBoolMap visited(mol.get_graph(), false);
  ShortToNodeVectorMap colorMap;
  ShortSet colorSet;
  bool is_tree = true;

  dfs(v, v, mol, visited, colorSet, colorMap, is_tree);

  for (ShortSet::iterator it=colorSet.begin(), end = colorSet.end(); it != end; ++it) {
    _colors.push_back(*it);
  }

  if (!is_tree) {
    //canonTree(g, colorMap);
    canonNauty(mol, colorMap, mol.get_atom_count());
  } else {
    canonNauty(mol, colorMap, mol.get_atom_count());
  }
}

void mogli::Canonization::dfs(const Node& current, const Node& last,
                              const Molecule& mol, NodeToBoolMap& visited,
                              ShortSet& colorSet, ShortToNodeVectorMap& colorMap,
                              bool& is_tree) {
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
      dfs(w, current, mol, visited, colorSet, colorMap, is_tree);
    } else if (w != last) {
      is_tree = false;
    }
  }

}

void mogli::Canonization::canonNauty(const Molecule& mol, const ShortToNodeVectorMap &colorMap, const unsigned int atom_count) {

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
  for (ShortVector::iterator it=_colors.begin(), end = _colors.end(); it != end; ++it) {
    NodeVector vector = colorMap.at(*it);
    for (NodeVector::iterator it2 = vector.begin(), end2 = vector.end(); it2 != end2; ++it2) {
      lab[i] = mol.get_id(*it2);
      ptn[i] = 1;
      ++i;
    }
    ptn[i-1] = 0;
  }

  for(EdgeIt e = mol.get_edge_iter(); e!=lemon::INVALID; ++e) {
    int u = mol.get_id(mol.get_u(e));
    int v = mol.get_id(mol.get_v(e));
    ADDONEEDGE(ng,u,v,m);
  }

  DYNALLOC2(graph,cg,cg_sz,atom_count,m,"malloc");
  densenauty(ng,lab,ptn,orbits,&options,&stats,m,atom_count,cg);

//  std::cout << "lab" << std::endl;
//  for (i = 0; i < atom_count; ++i) {
//    std::cout << lab[i] << std::endl;
//  }

  // TODO test with nodes > WORDSIZE
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

void mogli::Canonization::canonTree(const Molecule &mol, const ShortToNodeVectorMap &colorMap) {
  // TODO write own canonization method
}

//
// Created by M. Engler on 10/20/16.
//

#include <bitset>
#include "../include/canonization.h"

void mogli::Canonization::init(const Molecule &mol) {

  Node v = mol.get_node_iter();
  assert(v != lemon::INVALID);

  NodeToBoolMap visited(mol.get_graph(), false);
  ShortToNodeVectorMap colorMap;
  ShortSet colorSet;
  bool is_tree = true;

  dfs(v, v, mol, visited, colorSet, colorMap, is_tree);

  if (!is_tree) {
    //canonTree(g, colorMap);
    canonNauty(mol, colorSet, colorMap, mol.get_atom_count());
  } else {
    canonNauty(mol, colorSet, colorMap, mol.get_atom_count());
  }
}

void mogli::Canonization::init(const mogli::Molecule &mol, const mogli::Canonization::NodeToBoolMap &filter, const Node& root) {
  const FilterNodes subgraph(mol.get_graph(), filter);

  Node v = FilteredNodeIt(subgraph);
  assert(v != lemon::INVALID);

  FilteredNodeToBoolMap visited(subgraph, false);
  ShortToNodeVectorMap colorMap;
  ShortSet colorSet;
  bool is_tree = true;
  unsigned int node_count = 0;

  dfs(v, v, mol, subgraph, visited, colorSet, colorMap, is_tree, node_count);

  if (is_tree) {
    //canonTree(g, colorMap);
    canonNauty(mol, subgraph, colorSet, colorMap, root, node_count);
  } else {
    canonNauty(mol, subgraph, colorSet, colorMap, root, node_count);
  }
}

void mogli::Canonization::dfs(const Node& current, const Node& last, const Molecule& mol,
                              NodeToBoolMap& visited, ShortSet& colorSet, ShortToNodeVectorMap& colorMap,
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

void mogli::Canonization::dfs(const Node& current, const Node& last, const Molecule& mol,
                              const FilterNodes& subgraph, NodeToBoolMap& visited, ShortSet& colorSet,
                              ShortToNodeVectorMap& colorMap, bool& is_tree, unsigned int& node_count) {
  ++node_count;
  visited[current] = true;
  unsigned short color = mol.get_color(current);
  colorSet.insert(color);
  if (colorMap.find(color) == colorMap.end()) {
    NodeVector vector;
    colorMap.insert(std::pair<int, NodeVector>(color, vector));
  }
  colorMap.at(color).push_back(current);

  for (FilteredIncEdgeIt e = FilteredIncEdgeIt(subgraph, current); e != lemon::INVALID; ++e) {
    Node w = subgraph.oppositeNode(current, e);
    if (!visited[w]) {
      dfs(w, current, mol, subgraph, visited, colorSet, colorMap, is_tree, node_count);
    } else if (w != last) {
      is_tree = false;
    }
  }

}

void mogli::Canonization::canonNauty(const Molecule& mol,
                                     const ShortSet &colorSet,
                                     const ShortToNodeVectorMap &colorMap,
                                     const unsigned int atom_count) {

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
  NodeToIntMap nodes(mol.get_graph());
  NodeVector first_order;

  for (ShortSet::iterator it=colorSet.begin(), end = colorSet.end(); it != end; ++it) {
    NodeVector vector = colorMap.at(*it);
    for (NodeVector::iterator it2 = vector.begin(), end2 = vector.end(); it2 != end2; ++it2) {
      _colors.push_back(*it);
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
    _node_order.push_back(mol.get_id(first_order[lab[i]]));
  }

  for (i = 0; i < m*atom_count; ++i) {
    _canonization.push_back(static_cast<unsigned long>(cg[i]));
  }

}

void mogli::Canonization::canonNauty(const Molecule& mol,
                                     const FilterNodes &subgraph,
                                     const ShortSet &colorSet,
                                     const ShortToNodeVectorMap &colorMap,
                                     const Node& root,
                                     const unsigned int atom_count) {
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


  NodeToIntMap nodes(mol.get_graph());
  NodeVector first_order;

  _colors.push_back(mol.get_color(root));
  first_order.push_back(root);
  nodes[root] = 0;
  lab[0] = 0;
  ptn[0] = 0;

  int i = 1;
  for (ShortSet::iterator it=colorSet.begin(), end = colorSet.end(); it != end; ++it) {
    NodeVector vector = colorMap.at(*it);
    for (NodeVector::iterator it2 = vector.begin(), end2 = vector.end(); it2 != end2; ++it2) {
      if (*it2 == root)
        continue;
      _colors.push_back(*it);
      first_order.push_back(*it2);
      nodes[*it2] = i;
      lab[i] = i;
      ptn[i] = 1;
      ++i;
    }
    ptn[i-1] = 0;
  }

  for(FilteredEdgeIt e = FilteredEdgeIt(subgraph); e!=lemon::INVALID; ++e) {
    int u = nodes[subgraph.u(e)];
    int v = nodes[subgraph.v(e)];
    ADDONEEDGE(ng,u,v,m);
  }

  DYNALLOC2(graph,cg,cg_sz,atom_count,m,"malloc");
  densenauty(ng,lab,ptn,orbits,&options,&stats,m,atom_count,cg);

  for (i = 0; i < atom_count; ++i) {
    _node_order.push_back(mol.get_id(first_order[lab[i]]));
  }

  for (i = 0; i < m*atom_count; ++i) {
    _canonization.push_back(static_cast<unsigned long>(cg[i]));
  }
}

void mogli::Canonization::canonTree(const Molecule &mol, const ShortToNodeVectorMap &colorMap) {
  // TODO write tree canonization method
  // TODO in which order does nauty return the canonical graphs, if our tree alg has a different order, could we accidentally return the same canonization?
  // TODO compare tree, nauty, saucy & bliss
}

void mogli::Canonization::canonTree(const FilterNodes &subgraph, const ShortToNodeVectorMap &colorMap) {
 
}


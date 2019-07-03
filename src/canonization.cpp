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

#include "../include/canonization.h"

// public functions

const bool mogli::Canonization::is_isomorphic(const Canonization &other) const {
  const ShortVector& colors2 = other.get_colors();

  if (_colors.size() != colors2.size())
    return false;

  const LongVector& canonization2 = other.get_canonization();

  if (_canonization.size() != canonization2.size())
    return false;

  for (auto i1 = _colors.begin(), i2 = colors2.begin(),
           ie1 = _colors.end(), ie2 = colors2.end();
       i1 != ie1 && i2 != ie2; ++i1, ++i2) {
    if (*i1 != *i2)
      return false;
  }

  for (auto i1 = _canonization.begin(), i2 = canonization2.begin(),
           ie1 = _canonization.end(), ie2 = canonization2.end();
       i1 != ie1 && i2 != ie2; ++i1, ++i2) {
    if (*i1 != *i2)
      return false;
  }

  return true;
}

// protected functions

void mogli::Canonization::init(const Molecule &mol) {

  Node v = mol.get_node_iter();
  assert(v != lemon::INVALID);

  NodeToBoolMap visited(mol.get_graph(), false);
  ShortToNodeVectorMap colorMap;
  ShortSet colorSet;

  dfs(v, v, mol, visited, colorSet, colorMap);
  canonNauty(mol, colorSet, colorMap, mol.get_atom_count());
}

void mogli::Canonization::init(const mogli::Molecule &mol, const mogli::Canonization::NodeToBoolMap &filter, const Node& root) {
  const FilterNodes subgraph(mol.get_graph(), filter);

  Node v = FilteredNodeIt(subgraph);
  assert(v != lemon::INVALID);

  FilteredNodeToBoolMap visited(subgraph, false);
  ShortToNodeVectorMap colorMap;
  ShortSet colorSet;
  unsigned int node_count = 0;

  dfs(v, v, mol, subgraph, visited, colorSet, colorMap, node_count);
  canonNauty(mol, subgraph, colorSet, colorMap, root, node_count);
}

void mogli::Canonization::dfs(const Node& current, const Node& last, const Molecule& mol,
                              NodeToBoolMap& visited, ShortSet& colorSet, ShortToNodeVectorMap& colorMap) {
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
      dfs(w, current, mol, visited, colorSet, colorMap);
    }
  }

}

void mogli::Canonization::dfs(const Node& current, const Node& last, const Molecule& mol,
                              const FilterNodes& subgraph, NodeToBoolMap& visited, ShortSet& colorSet,
                              ShortToNodeVectorMap& colorMap, unsigned int& node_count) {
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
      dfs(w, current, mol, subgraph, visited, colorSet, colorMap, node_count);
    }
  }

}

void mogli::Canonization::canonNauty(const Molecule& mol,
                                     const ShortSet &colorSet,
                                     const ShortToNodeVectorMap &colorMap,
                                     unsigned int atom_count) {

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

  for (const auto & it : colorSet) {
    NodeVector vector = colorMap.at(it);
    for (const auto & it2 : vector) {
      _colors.push_back(it);
      first_order.push_back(it2);
      nodes[it2] = i;
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
                                     unsigned int atom_count) {
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
  for (const auto & it : colorSet) {
    NodeVector vector = colorMap.at(it);
    for (const auto & it2 : vector) {
      if (it2 == root)
        continue;
      _colors.push_back(it);
      first_order.push_back(it2);
      nodes[it2] = i;
      lab[i] = i;
      ptn[i] = 1;
      ++i;
    }
    ptn[i-1] = 0;
  }

  for(auto e = FilteredEdgeIt(subgraph); e!=lemon::INVALID; ++e) {
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

/*
 * product.h
 *
 *  Created on: 21-jan-2014
 *      Author: M. El-Kebir
 *
 *  Modified by M. Engler on 26/10/16.
 */

#ifndef MOGLI_PRODUCT_H
#define MOGLI_PRODUCT_H

#include "molecule.h"
#include "canonization.h"
#include "isomorphism.h"
#include <set>
#include <vector>

namespace mogli {

  class Product {

  private:
    typedef std::vector<NodeVector> NodeMatrix;
    typedef typename Graph::template EdgeMap<bool> EdgeToBoolMap;
    typedef typename Graph::template NodeMap<int> NodeToIntMap;

    typedef typename Graph::template NodeMap<Node> NodeToNodeMap;
    typedef std::multiset<int> IntSet;
    typedef typename Graph::template NodeMap<IntSet> NodeToIntSetMap;
    typedef typename Graph::template NodeMap<NodeVector> NodeToNodeVectorMap;

    typedef typename Graph::template NodeMap<Canonization> NodeToCanonizationMap;


    // old typdefs
    typedef typename Graph::template NodeMap<IntSet> IntSetNodeMap;
    typedef typename Graph::template NodeMap<NodeVector> NodeVectorMap;

    // end old typedefs

    Molecule& _mol1;
    Molecule& _mol2;
    const int _shell;
    Graph _g;
    NodeToNodeMap _mol1ToG;
    NodeToNodeMap _mol2ToG;
    NodeToNodeMap _gToMol1;
    NodeToNodeMap _gToMol2;
    EdgeToBoolMap _connectivityEdge;
    NodeToNodeVectorMap _g1ToDeg1Neighbors;
    NodeToNodeVectorMap _g2ToDeg1Neighbors;

    int _numNodes;
    int _numEdges;

  public:
    Product(Molecule& mol1, Molecule& mol2, int shell, int gen)
      : _mol1(mol1)
      , _mol2(mol2)
      , _shell(shell)
      , _g()
      , _mol1ToG(mol1.get_graph(), lemon::INVALID)
      , _mol2ToG(mol2.get_graph(), lemon::INVALID)
      , _gToMol1(_g)
      , _gToMol2(_g)
      , _connectivityEdge(_g)
      , _g1ToDeg1Neighbors(mol1.get_graph())
      , _g2ToDeg1Neighbors(mol2.get_graph())
      , _numNodes(0)
      , _numEdges(0) {
      if (gen == 0)
        generate();
      else
        generate_old();
    }

    const Graph& get_graph() const {
      return _g;
    }

    const std::string to_string() const {
      std::stringstream buffer;
      buffer << "nodes:[";
      NodeIt v = NodeIt(_g);
      if (v != lemon::INVALID) {
        buffer << _g.id(v) << ":" << _mol1.get_label(_gToMol1[v]) << "x" << _mol2.get_label(_gToMol2[v]);
        ++v;
        for (; v != lemon::INVALID; ++v) {
          buffer << ", " << _g.id(v) << ":" << _mol1.get_label(_gToMol1[v]) << "x" << _mol2.get_label(_gToMol2[v]);
        }
      }
      buffer << "],edges:[";
      EdgeIt e = EdgeIt(_g);
      if (e != lemon::INVALID) {
        std::string u1 = _mol1.get_label(_gToMol1[_g.u(e)]);
        std::string u2 = _mol2.get_label(_gToMol2[_g.u(e)]);
        std::string v1 = _mol1.get_label(_gToMol1[_g.v(e)]);
        std::string v2 = _mol2.get_label(_gToMol2[_g.v(e)]);
        buffer << "(" << u1 << "x" << u2 << "," << v1 << "x" << v2 << ")";
        ++e;
        for (; e != lemon::INVALID; ++e) {
          std::string u1 = _mol1.get_label(_gToMol1[_g.u(e)]);
          std::string u2 = _mol2.get_label(_gToMol2[_g.u(e)]);
          std::string v1 = _mol1.get_label(_gToMol1[_g.v(e)]);
          std::string v2 = _mol2.get_label(_gToMol2[_g.v(e)]);
          buffer << ", (" << u1 << "x" << u2 << "," << v1 << "x" << v2 << ")";
        }
      }
      buffer << "]";
      std::string str = buffer.str();
      return str;
    }
    
  private:

    void generate0();

    void generate();

    // old functions

    void generate_old();

    void determineDegrees(const Graph& g, NodeToIntMap& deg);

    void generate(const Molecule& mol,
             const NodeToIntMap& deg,
             IntSetNodeMap& intSet,
             IntSetNodeMap& degSet);

    void dfs(const NodeToIntMap& deg,
             const Node v, const int depth,
             const Molecule& mol,
             NodeToBoolMap& visited,
             IntSet& s, IntSet& ds);

    void generateDeg1NeighborSet(const Graph& g,
                            const NodeToIntMap& deg,
                            NodeVectorMap& deg1NeighborMap);

    //end old functions

    void generate_subgraph_canonization(Molecule &mol, const Node &v, NodeToCanonizationMap &map);

    void dfs(Molecule &mol, const Node &v, int depth, NodeToBoolMap &visited, NodeToBoolMap &filter);

  };

  inline std::ostream& operator<<(std::ostream& outputStream, const Product& product) {
    outputStream << product.to_string();
    return outputStream;
  }

} // namespace mogli

#endif // MOGLI_PRODUCT_H

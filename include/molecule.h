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

#ifndef MOGLI_MOLECULE_H
#define MOGLI_MOLECULE_H

#include <istream>
#include <ostream>

#include <lemon/adaptors.h>
#include "periodictable.h"
#include "types.h"


namespace mogli {

  typedef typename Graph::template NodeMap<int> NodeToIntMap;
  typedef typename Graph::template EdgeMap<bool> EdgeToBoolMap;
  typedef typename Graph::template EdgeMap<IntSet> EdgeToIntSetMap;
  typedef lemon::SubGraph<Graph, NodeToBoolMap, EdgeToBoolMap> SubGraph;

  /**
   * LGF-formatter for reading and writing molecular graphs.
   */
  class LGFIOConfig {

  public:

    /**
     * Initialize an LGF formatter.
     *
     * @param[in] id_property       Name of the atom ID column.
     * @param[in] color_property    Name of the element number column.
     */
    LGFIOConfig(std::string id_property, std::string color_property) :
        _id_property(std::move(id_property)),
        _color_property(std::move(color_property)),
        _bool_node_props(),
        _int_node_props(),
        _double_node_props(),
        _string_node_props() {}

  private:

    std::string _id_property;
    std::string _color_property;
    StringSet _bool_node_props;
    StringSet _int_node_props;
    StringSet _double_node_props;
    StringSet _string_node_props;

    LGFIOConfig() :
      _id_property(std::string("label")),
      _color_property(std::string("atomType")),
      _bool_node_props(),
      _int_node_props(),
      _double_node_props(),
      _string_node_props() {
      add_double_node_prop("partial_charge");
      add_string_node_prop("label2");
      add_double_node_prop("coordX");
      add_double_node_prop("coordY");
      add_double_node_prop("coordZ");
      add_int_node_prop("initColor");
    }

  public:

    /**
     * Add a bool atom property column.
     *
     * @param[in] property  Column name.
     * @return              Updated LGF formatter.
     */
    LGFIOConfig& add_bool_node_prop(const std::string& property) {
      _bool_node_props.insert(property);
      return *this;
    }

    /**
     * Add an int atom property column.
     *
     * @param[in] property  Column name.
     * @return              Updated LGF formatter.
     */
    LGFIOConfig& add_int_node_prop(const std::string& property) {
      _int_node_props.insert(property);
      return *this;
    }

    /**
     * Add a double atom property column
     *
     * @param[in] property  Column name.
     * @return              Updated LGF formatter.
     */
    LGFIOConfig& add_double_node_prop(const std::string& property) {
      _double_node_props.insert(property);
      return *this;
    }

    /**
     * Add a string atom property column.
     *
     * @param[in] property  Column name.
     * @return              Updated LGF formatter.
     */
    LGFIOConfig& add_string_node_prop(const std::string& property) {
      _string_node_props.insert(property);
      return *this;
    }

    /**
     * Returns the name of the atom ID column.
     *
     * @return  Column name.
     */
    const std::string &get_id_property() const {
      return _id_property;
    }

    /**
     * Returns the name of the element number column.
     *
     * @return  Column name.
     */
    const std::string &get_color_property() const {
      return _color_property;
    }

    /**
     * Returns the names of all bool atom property columns.
     *
     * @return  Column names.
     */
    const StringSet &get_bool_node_props() const {
      return _bool_node_props;
    }

    /**
     * Returns the names of all int atom property columns.
     *
     * @return  Column names.
     */
    const StringSet &get_int_node_props() const {
      return _int_node_props;
    }

    /**
     * Returns the names of all double atom property columns.
     *
     * @return  Column names.
     */
    const StringSet &get_double_node_props() const {
      return _double_node_props;
    }

    /**
     * Returns the names of all string atom property columns.
     *
     * @return  Column names.
     */
    const StringSet &get_string_node_props() const {
      return _string_node_props;
    }

    /**
     * @brief Returns the default LGF formatter.
     *
     * Default formatting is: label | atomType | partial_charge | label2 | coordX | coordY | coordZ | initColor
     *
     * @return  Default LGF formatter.
     */
    static LGFIOConfig &get_default() {
      static LGFIOConfig instance;
      return instance;
    }

  };

  struct sort_tuple {
    bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
      if (left.first == right.first) {
        return left.second < right.second;
      }
      return left.first < right.first;
    }
  };

  /**
   * Molecular graph.
   */
  class Molecule {

  public:

  protected:

    typedef typename Graph::template NodeMap<Any> NodeToAnyMap;

    typedef typename UniquePtrMap<std::string, NodeToAnyMap>::type StringToAnyTypeMapMap;

    typedef typename Graph::template NodeMap<unsigned short> NodeToUShortMap;
    typedef typename Graph::template NodeMap<double> NodeToDoubleMap;
    typedef typename Graph::template NodeMap<std::string> NodeToStringMap;

    typedef typename UniquePtrMap<int, Node>::type IntToNodeMap;

    Graph _g;

    unsigned int _atom_count;

    NodeToUShortMap _colors;
    NodeToIntMap _node_to_id;
    IntToNodeMap _id_to_node;
    int _max_uid;

    StringToAnyTypeMapMap _properties;

    PeriodicTable &_perdiodic_table;

    short _is_connected;

  public:

    /**
     * Empty constructor.
     */
    Molecule() : _g(),
                 _colors(_g),
                 _node_to_id(_g),
                 _max_uid(0),
                 _atom_count(0),
                 _is_connected(-1),
                 _properties(),
                 _perdiodic_table(PeriodicTable::get_default()) {}

    /**
     * Constructor with custom periodic table.
     *
     * @param[in] periodic_table    Periodic table.
     */
    explicit Molecule(PeriodicTable &periodic_table) : _g(),
                                             _colors(_g),
                                             _node_to_id(_g),
                                             _max_uid(0),
                                             _atom_count(0),
                                             _is_connected(-1),
                                             _properties(),
                                             _perdiodic_table(periodic_table) {}

    PeriodicTable &get_perdiodic_table() const {
      return _perdiodic_table;
    }

    /**
     * Returns the names of all atom properties.
     *
     * @param[out] properties   Property names.
     */
    void get_properties(StringVector& properties) const {
      for (const auto & it : _properties) {
        properties.push_back(it.first);
      }

    }

    /**
     * Returns atom property with this name.
     *
     * @param[in] node      Atom.
     * @param[in] property  Property name.
     * @return              Property value.
     */
    const Any get_property(Node node, const std::string& property) const {
      return _properties.at(property)->operator[](node);
    }

    /**
     * Set atom property.
     *
     * @param[in] node      Atom.
     * @param[in] property  Property name.
     * @param[in] value     Property value.
     */
    void set_property(Node node, const std::string& property, Any value) {
      if (_properties.count(property) == 0) {
        _properties[property] = std::make_unique<NodeToAnyMap>(_g);
      }
      _properties[property]->operator[](node) = std::move(value);
    }

    // add atoms & edges

    /**
     * Add new atom.
     *
     * @param[in] element   Element type.
     * @return              Atom.
     */
    const Node add_atom(std::string element) {
      return add_atom(_perdiodic_table.get_number(std::move(element)));
    }

    /**
     * Add new atom.
     *
     * @param[in] id        Atom ID.
     * @param[in] element   Element type.
     * @return              Atom.
     */
    const Node add_atom(int id, std::string element) {
      return add_atom(id, _perdiodic_table.get_number(std::move(element)));
    }

    /**
     * Add new atom.
     *
     * @param[in] color Element number.
     * @return          Atom.
     */
    const Node add_atom(unsigned short color) {
      return add_atom(_max_uid+1, color);
    }

    /**
     * Add new atom.
     *
     * @param[in] id    Atom ID.
     * @param[in] color Element type.
     * @return          Atom.
     */
    const Node add_atom(int id, unsigned short color);

    /**
     * Add bond.
     *
     * @param[in] u Atom.
     * @param[in] v Atom.
     */
    const Edge add_edge(const Node &u, const Node &v) {
      _is_connected = -1;
      return _g.addEdge(u,v);
    }

    // getters & iterators

    /**
     * Returns the number of atoms.
     *
     * @return  Number of atoms.
     */
    const unsigned int get_atom_count() const {
      return _atom_count;
    }

    const Graph &get_graph() const {
      return _g;
    }

    /**
     * Returns an iterator over all atoms.
     *
     * @return  Atom iterator.
     */
    const NodeIt get_node_iter() const {
      return NodeIt(_g);
    }

    /**
     * Returns an iterator over all bonds.
     *
     * @return  Bond iterator.
     */
    const EdgeIt get_edge_iter() const {
      return EdgeIt(_g);
    }

    /**
     * Returns an iterator over all incident bonds.
     *
     * @return  Bond iterator.
     */
    const IncEdgeIt get_inc_edge_iter(const Node &node) const {
      return IncEdgeIt(_g, node);
    }

    /**
     * Returns the atom on the opposite side of the bond.
     *
     * @param[in] node  Atom.
     * @param[in] edge  Bond.
     * @return          Opposite atom.
     */
    const Node get_opposite_node(const Node &node, const Edge &edge) const {
      return _g.oppositeNode(node, edge);
    }

    /**
     * Returns the first atom of this bond.
     *
     * @param[in] edge  Bond.
     * @return          First atom.
     */
    const Node get_u(const Edge &edge) const {
      return _g.u(edge);
    }

    /**
     * Returns the second atom of this bond.
     *
     * @param[in] edge  Bond.
     * @return          Second atom.
     */
    const Node get_v(const Edge &edge) const {
      return _g.v(edge);
    }

    /**
     * Test if atom with this ID exists.
     *
     * @param[in] id    Atom ID.
     * @return          True, if atom with this ID exists, false otherwise.
     */
    const bool has_node_with_id(int id) const {
      return _id_to_node.find(id) != _id_to_node.end();
    }

    /**
     * Returns atom with this ID.
     *
     * @param[in] id    Atom ID.
     * @return          Atom.
     */
    const Node get_node_by_id(int id) const {
      return *_id_to_node.at(id);
    }

    /**
     * Returns the ID of this atom.
     *
     * @param[in] node  Atom.
     * @return          Atom ID.
     */
    const int get_id(const Node &node) const {
      return _node_to_id[node];
    }

    /**
     * Returns the element number of this atom.
     *
     * @param[in] node  Atom.
     * @return          Element number.
     */
    const unsigned short get_color(const Node &node) const {
      return _colors[node];
    }

    /**
     * Returns the element type of this atom.
     *
     * @param[in] node  Atom.
     * @return          Element type.
     */
    const std::string get_element(const Node &node) const {
      return _perdiodic_table.get_element(get_color(node));
    }

    const std::string get_color_name(const Node &node) const {
      return _perdiodic_table.get_color(get_color(node));
    }

    // I/O

    /**
     * Export GML.
     *
     * @param[in]  label Graph name.
     * @param[out] out   Output stream.
     * @param[in]  raw   Write raw strings.
     */
    void write_gml_stream(const std::string& label, std::ostream &out, bool raw = false);

    /**
     * Export GML.
     *
     * @param[in] label Graph name.
     * @return          Molecular graph in default LGF format.
     */
    std::string write_gml(const std::string& label);

    /**
     * Export default formatted LGF.
     *
     * @param[out] out  Output stream.
     */
    void write_lgf_stream(std::ostream &out);

    /**
     * Export LGF in given format.
     *
     * @param[out] out  Output stream.
     * @param[in] config    LGF formatter.
     */
    void write_lgf_stream(std::ostream &out, const LGFIOConfig& config);

    /**
     * Export default formatted LGF.
     *
     * @return Molecular graph in default LGF format.
     */
    std::string write_lgf();

    /**
     * Export LGF in given format.
     *
     * @param[in] config    LGF formatter.
     * @return              Molecular graph in given LGF format.
     */
    std::string write_lgf(const LGFIOConfig &config);

    /**
     * Import default formatted LGF.
     *
     * @param[in] in    Input stream.
     */
    void read_lgf_stream(std::istream &in);

    /**
     * Import LGF in given format.
     *
     * @param[in] in        Input stream.
     * @param[in] config    LGF formatter.
     */
    void read_lgf_stream(std::istream &in, const LGFIOConfig &config);

    /**
     * Import default formatted LGF.
     *
     * @param[in] in    Molecular graph in LGF format.
     */
    void read_lgf(const std::string &in);

    /**
     * Import LGF in given format.
     *
     * @param[in] in        Molecular graph in LGF format.
     * @param[in] config    LGF formatter.
     */
    void read_lgf(const std::string &in, const LGFIOConfig &config);

    // connected?

    /**
     * Test if molecular graph is connected.
     *
     * @return  True, if molecular graph is connected, false otherwise.
     */
    const bool is_connected() {
      if (_is_connected < 0) {
        _is_connected = is_connected0();
      }
      return _is_connected == 1;
    }

    /**
     * Returns the connected components of the molecular graph.
     *
     * @param[out] components   Connected components.
     */
    void get_connected_components(SharedPtrVector<Molecule>::type &components);

    /**
     * @brief Balanced split of the molecular graph.
     *
     * Tries to split the molecule into overlapping smaller components of roughly equal size (less or equal max_size),
     *
     * @param[in] max_size      Maximal size of the resulting components.
     * @param[in] shell         Shell size (Overlap at the split-regions).
     * @param[out] components   Resulting components.
     */
    int split(int max_size, int shell, SharedPtrVector<Molecule>::type &components);

    /**
     * Test if this molecular graph is isomorphic to another molecular graph.
     *
     * @param[in] other Other molecular graph.
     * @return          True, if they are isomorphic, false otherwise.
     */
    const bool is_isomorphic(Molecule &other) const;

    /**
     * Export molecular graph to dot (graphviz) format.
     *
     * @return  Graph in dot format.
     */
    virtual const std::string print_dot() const;

    /**
     * Export molecular graph to dot (graphviz) format with selected atom properties.
     *
     * @param[in] properties    Atom properties to print.
     * @return                  Graph in dot format.
     */
    const std::string print_dot(const StringVector &properties) const;

    /**
     * Export molecular graph to dot (graphviz) format.
     *
     * @param[out] out  Output stream.
     */
    virtual const void print_dot(std::ostream& out) const;

    /**
     * Export molecular graph to dot (graphviz) format with selected atom properties.
     *
     * @param[out] out          Output stream.
     * @param[in]  properties   Atom properties to print.
     */
    const void print_dot(std::ostream &out, const StringVector &properties) const;

  protected:

    void balanced_cut(const std::shared_ptr<NodeToBoolMap>& subgraph_map,
                      int max_size, int shell,
                      SharedPtrVector<NodeToBoolMap>::type &partitions);

    int get_bctree(const std::shared_ptr<SubGraph>& subgraph,
                   Graph& bctree,
                   NodeToIntMap& bctree_node_weights, NodeToIntMap& bctree_2_bicon_comp,
                   SharedPtrVector<std::set<Node>>::type & biconnected_nodes,
                   EdgeToIntSetMap& bctree_shell_ids_u, EdgeToIntSetMap& bctree_shell_ids_v, int shell);

    void get_subgraph(Graph& bctree,
                      NodeToIntMap& bctree_components,
                      NodeToIntMap& bctree_2_bicon_comp,
                      SharedPtrVector<std::set<Node>>::type &biconnected_nodes,
                      IntSet& shell_ids,
                      int component,
                      const std::shared_ptr<NodeToBoolMap>& subgraph_map);

    void find_shell_nodes(const std::shared_ptr<SubGraph>& subgraph,
                          NodeToIntMap& bctree_2_bicon_comp,
                          SharedPtrVector<std::set<Node>>::type &biconnected_nodes,
                          Graph& bctree,
                          int shell,
                          EdgeToIntSetMap& bctree_shell_ids_u,
                          EdgeToIntSetMap& bctree_shell_ids_v);

    void bfs_find_shell(const std::shared_ptr<SubGraph>& subgraph,
                        const Node &v, int shell,
                        SubGraph::NodeMap<IntSet>& shell_ids);

    Edge get_cut_edge(Graph& bctree, NodeToIntMap &node_weights,
                      EdgeToIntSetMap& shell_ids_u, EdgeToIntSetMap& shell_ids_v, int N);

    virtual short is_connected0();

    void dfs(const Node& current, NodeToBoolMap& visited);

  };

}

#endif //MOGLI_MOLECULE_H

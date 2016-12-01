# mogli
molecular graph library


## molecules

A simple example of how to work with the molecule class:

    // create a molecule and register its properties
    Molecule mol2;
    mol2.add_string_property("label").add_string_property("label2").add_double_property("partial_charge");
    
    // read molecule from lgf
    std::ifstream ifs("./data/min_1.lgf", std::ifstream::in);
    mol2.read_lgf_stream(ifs);
    ifs.close();

    // iterate over all nodes and print the properties
    for (NodeIt n = mol2.get_node_iter(); n != lemon::INVALID; ++n) {
        std::cout << mol2.get_iacm_element(n) << " "
                  << mol2.get_string_property(n, "label") << " "
                  << mol2.get_string_property(n, "label2") << " "
                  << mol2.get_double_property(n, "partial_charge") << std::endl;
    }

String properties are supposed to be unique, so that you can find a node by its string property. However, the molecule class does not check if a property is actually unique. If the string property is not unique, the method will return the last node, where the property was set.

    // find the node labelled "H1"
    Node u = mol.get_node_by_string_property("label", "H1");
    
Check if the molecular graph is connected:

    // is the graph connected?
    std::cout << mol.is_connected() << std::endl; 

Iterating over all neighbors of a node:

    // iterate neighbors w of node v
    for (IncEdgeIt e = mol.get_inc_edge_iter(v); e != lemon::INVALID; ++e) {
            Node w = mol.get_opposite_node(v, e);
    }
    
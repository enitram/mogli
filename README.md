# mogli
The molecular graph library

## Dependencies

- LEMON 1.3
- Boost 1.60.0

## Compiling

First, get Boost (including `boost-python`) 1.60.0 or newer using your package manager or download and manually install from here: http://www.boost.org/users/download/
If you have an older version of Boost, change line 16 in the mogli `CMakeLists.txt`

    find_package( Boost 1.60.0 REQUIRED COMPONENTS python)
     
to your version (Use older versions at your own risk. In theory they should work fine, but have not been tested).
 
Second, LEMON 1.3 needs to be installed:

    wget http://lemon.cs.elte.hu/pub/sources/lemon-1.3.tar.gz
    tar xvzf lemon-1.3.tar.gz
    cd lemon-1.3
    cmake -DCMAKE_INSTALL_PREFIX=~/lemon
    make install

Note: On Mac OS 10.9, comment out the following two lines and add the code below at line 159 in the LEMON `CMakeLists.txt` before doing make install.

    #ADD_SUBDIRECTORY(demo) 
    #ADD_SUBDIRECTORY(tools)

    if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
        set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++ " )
    endif()

You can remove the LEMON sources now, i.e., `rm -rf lemon-1.3`.

Now you can compile mogli:
    
    cd /path/to/mogli
    git submodule init
    git submodule update
    mkdir build
    cd build
    cmake ..
    make

In case auto-detection of LEMON fails, do

    cmake -DLIBLEMON_ROOT=/path/to/lemon ..

## Molecules

A simple example of how to work with the molecule class:

    // create a molecule and register its properties
    Molecule mol;
    
    // read molecule from lgf
    std::ifstream ifs("./data/min_1.lgf", std::ifstream::in);
    mol.read_lgf_stream(ifs);
    ifs.close();

    // iterate over all nodes and print the properties
    for (NodeIt n = mol.get_node_iter(); n != lemon::INVALID; ++n) {
        std::cout << mol.get_element(n) << " "
                  << mol.get_id(n) << " "
                  << boost::any_cast<std::string>(mol.get_property(n, "label2")) << " "
                  << boost::any_cast<double>(mol.get_property(n, "partial_charge")) << std::endl;
    }
    
Check if the molecular graph is connected:

    // is the graph connected?
    std::cout << mol.is_connected() << std::endl; 

Iterating over all neighbors of a node:

    // iterate neighbors w of node v
    for (IncEdgeIt e = mol.get_inc_edge_iter(v); e != lemon::INVALID; ++e) {
            Node w = mol.get_opposite_node(v, e);
    }
    

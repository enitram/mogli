# mogli
The molecular graph library.

The focus of mogli is finding [maximal common fragments](https://doi.org/10.7287/peerj.preprints.3250v1). Maximal common
fragments are common subgraphs of two molecular graphs, where not just atom pairs are matched but also their
neighborhood. The size of neighborhood is determined by the shell size *k* - the maximal distance (number of bonds) from
the central atom of the neighborhood. The core atoms of a fragment are the atoms more than *k* bonds away from the
fragment border, all other atoms are shell atoms.
    
mogli provides both a C++ and a Python API and is using the following awesome third party libraries:

* [lemon graph library](https://lemon.cs.elte.hu/)
* [nauty graph isomorphism solver](http://users.cecs.anu.edu.au/~bdm/nauty/)
* [LAD subgraph isomorphism solver](https://perso.liris.cnrs.fr/christine.solnon/LAD.html)
* [boost::dynamic_bitset](https://www.boost.org/doc/libs/1_71_0/libs/dynamic_bitset/dynamic_bitset.html)
* [msgpack-c](https://github.com/msgpack/msgpack-c)
* [pybind11](https://github.com/pybind/pybind11)
* [Catch2](https://github.com/catchorg/Catch2)

## Installation

### C++

Build shared library `libcmogli.so`:

```bash
mkdir build
cd build
cmake .. -DBUILD_PYTHON=OFF
cmake --build cmogli
```

### Python

Installation is easily done with:

```bash
python setup.py install
```

## Tests

mogli's unit tests cover all basic functions. Usually it is not necessary to run the tests, unless you changed the code.

Compile the unit tests with:

```bash
mkdir build_test
cd build_test
cmake .. -DBUILD_TESTS=ON
cmake --build unit_tests
```

Then, run all unit tests:

```bash
ctest -V
```

Running only the C++ tests:

```bash
ctest -V -L cpp
```

or the Python tests:

```bash
ctest -V -L python
```

## Usage

### The LGF file format

mogli imports and exports molecules from and to the lemon graph format (LGF). The format is table-based and
a standard mogli LGF file looks like this:

```
@nodes
partial_charge label label2 atomType coordX coordY coordZ initColor
-0.048         0     C1     12       -0.765 -0.000 0.000  1
0.016          1     H3     20       -1.164 -0.813 0.619  2
0.016          2     H4     20       -1.164 -0.129 -1.013 3
0.016          3     H5     20       -1.164 0.942  0.394  4
@edges
    label
0 1 0
0 2 1
0 3 2
```

The atomType are the IACM atom types used by the [ATB](https://dx.doi.org/10.1021/ct200196m) database.

### C++

#### Molecules

A simple example of how to work with the molecule class:

```cpp
#include <molecule.h>
using namespace mogli;

...

// create an empty molecule
Molecule mol;

// read molecule data from standard lgf
std::ifstream ifs("./data/min_1.lgf", std::ifstream::in);
mol.read_lgf_stream(ifs);
ifs.close();

// iterate over all atoms and print the properties
for (NodeIt n = mol.get_node_iter(); n != lemon::INVALID; ++n) {
    std::cout << mol.get_element(n) << " "
              << mol.get_id(n) << " "
              << std::get<std::string>(mol.get_property(n, "label2")) << " "
              << std::get<double>(mol.get_property(n, "partial_charge")) << std::endl;
}
```
    
Check if the molecular graph is connected:

```cpp
// is the graph connected?
std::cout << mol.is_connected() << std::endl; 
```

Iterating over all neighbors of an atom:

```cpp
// iterate neighbors w of node v
for (IncEdgeIt e = mol.get_inc_edge_iter(v); e != lemon::INVALID; ++e) {
        Node w = mol.get_opposite_node(v, e);
}
```

#### Custom LGF files

Importing custom LGF files is also supported. For example, this LGF file has an atom ID column "id",
an atom type column "element" and bool and int atom property columns.

```
@nodes
id element bool int
0  0       0    1
1  0       1    2
...
```

The custom LGF file can be imported with:

```cpp
LGFIOConfig config("id", "element");
config.add_bool_node_prop("bool").add_int_node_prop("int")

std::ifstream ifs("./data/min_2.lgf", std::ifstream::in);
mol.read_lgf_stream(ifs, config);
ifs.close();
```

#### Periodic table of elements

By default, mogli uses the IACM atom types used by the [ATB](https://dx.doi.org/10.1021/ct200196m) database. The
atom types are managed by the PeriodicTable class. To create a custom periodic table, it is possible to create a new
periodic table from scratch, change the default periodic table or copy and then change the default periodic table:

```cpp
#include <periodictable.h>
...

// copy the default IACM periodic table
PeriodicTable table(PeriodicTable::get_default());

// add a new element "Foo"
table.add_uncolored(132, "Foo");

// make "H" and "HC" equivalent
table.make_equivalent({20, 21});

// create an empty molecule with this custom periodic table
Molecule mol(table);
```

#### Maximal common fragments

Maximal common fragments are common subgraphs of two molecular graphs, where not just atom pairs are matched but also
their neighborhood. The size of neighborhood is determined by the shell size.

Enumerating all maximal common fragments:

```cpp
#include <mcf.h>
...

FragmentVector fragments;
MatchVector matches1, matches2;
int shell = 1;
int timeout = 10;
// t: False, if timeout occured, true otherwise. 
bool t = maximal_common_fragments(mol1,      // Matched against mol2 to find maximal common fragments.
                                  mol2,      // Matched against mol1 to find maximal common fragments.
                                  fragments, // Vector of resulting fragments.
                                  matches1,  // Resulting mappings from fragment atom IDs to atom IDs in mol1.
                                  matches2,  // Resulting mappings from fragment atom IDs to atom IDs in mol2.
                                  shell,     // Shell size.
                                  timeout);  // Timeout in seconds.
```

To only retrieve the largest common fragments, set the maximum parameter of `maximal_common_fragments` to `true`.

### Python

#### Molecules

A simple example of how to work with the molecule class:

```python
from mogli import Molecule

# create an empty molecule
mol = Molecule()

# read molecule data from standard lgf
with open('./data/min_1.lgf') as f:
    mol.read_lgf(f.read())

# iterate over all atoms and print the properties
for n in mol.get_node_iter():
    print(mol.get_element(n),
          mol.get_id(n),
          mol.get_property(n, 'label2'),
          mol.get_property(n, 'partial_charge'))
```
    
Check if the molecular graph is connected:

```python
# is the graph connected?
print(mol.is_connected()) 
```

Iterating over all neighbors of an atom:

```python
# iterate neighbors w of atom v
for e in mol.get_inc_edge_iter(v):
    w = mol.get_opposite_node(v, e)
```

#### Custom LGF files

Importing custom LGF files is also supported. For example, this LGF file has an atom ID column "id",
an atom type column "element" and bool and int atom property columns.

```
@nodes
id element bool int
0  0       0    1
1  0       1    2
...
```

The custom LGF file can be imported with:

```python
from mogli import LGFIOConfig
...

config = LGFIOConfig('id', 'element')
config.add_bool_node_prop('bool').add_int_node_prop('int')

with open('./data/min_2.lgf') as f:
    mol.read_lgf(f.read(), config)
```

#### Periodic table of elements

By default, mogli uses the IACM atom types used by the [ATB](https://dx.doi.org/10.1021/ct200196m) database. The
atom types are managed by the PeriodicTable class. To create a custom periodic table, it is possible to create a new
periodic table from scratch, change the default periodic table or copy and then change the default periodic table:

```python
from mogli import PeriodicTable
...

# copy the default IACM periodic table
table = PeriodicTable(PeriodicTable.get_default())

# add a new element "Foo"
table.add_uncolored(132, 'Foo')

# make "H" and "HC" equivalent
table.make_equivalent(20, 21)

# create an empty molecule with this custom periodic table
mol = Molecule(table)
```

#### Maximal common fragments

Maximal common fragments are common subgraphs of two molecular graphs, where not just atom pairs are matched but also
their neighborhood. The size of neighborhood is determined by the shell size.

Enumerating all maximal common fragments:

```python
from mogli import maximal_common_fragments
...

shell, timeout = 1, 10
# t:         False, if timeout occured, true otherwise.
# fragments: Vector of resulting fragments.
# matches1:  Resulting mappings from fragment atom IDs to atom IDs in mol1.
# matches2:  Resulting mappings from fragment atom IDs to atom IDs in mol2.
t, fragments, matches1, matches2 = maximal_common_fragments(mol1,    # Matched against mol2 to find maximal common fragments.
                                                            mol2,    # Matched against mol1 to find maximal common fragments.
                                                            shell,   # Shell size.
                                                            timeout) # Timeout in seconds.
```

To only retrieve the largest common fragments, set the maximum parameter of `maximal_common_fragments` to `true`.
    
## Known issues

If you get this error:

```bash
fatal error: Python.h: No such file or directory
```
    
You need to install python-devel (Fedora) or python-dev (Ubuntu).

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

#include <chrono>
#include <dirent.h>
#include <fstream>

#include <lemon/arg_parser.h>
#include <lemon/connectivity.h>
#include "bronkerbosch.h"
#include "fragment.h"
#include "mcf.h"


using namespace mogli;

std::tuple<float, float, bool, int> run(Product& product, int min_core, int max_core, FragmentVector& fragments) {
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  BronKerbosch bk(product, min_core, max_core, false);
  bk.run(-1);

  int degeneracy = static_cast<int>(bk.computeDegeneracy());

  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<float> bk_duration = t2 - t1;

  NodeVectorVector cliques = bk.getMaxCliques();
  bool c = cliques.size() == 1 && cliques[0].size() == lemon::countNodes(product.get_graph());
  for (NodeVectorVector::const_iterator it = cliques.begin(), end = cliques.end(); it != end; ++it) {
    IntToIntMap g_to_mol1, g_to_mol2;
    auto fragment = std::make_shared<Fragment>(product, *it, g_to_mol1, g_to_mol2);

    if (fragment->get_core_atom_count() > 1) {
      fragments.push_back(fragment);
    }
  }

  std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<float> gen_duration = t3 - t2;

  return std::make_tuple(bk_duration.count(), gen_duration.count(), c, degeneracy);
}

void print_int_list(std::vector<int> list) {
  std::cout << "\"[";
  if (list.size() > 0) {
    std::cout << list[0];
    for (int i = 1; i < list.size(); ++i) {
      std::cout << ", " << list[i];
    }
  }
  std::cout << "]\"\t";
}

int main(int argc, char** argv) {

  lemon::ArgParser ap(argc, (const char *const *) argv);

  int gen_int = 1;
  int shell = 3;
  int min_core = 3;
  int max_core = std::numeric_limits<int>::max();
  int molid1 = 0;
  int molid2 = 0;

  std::string in_path = "";
  std::string out_path = "../data/runtime";

  ap.refOption("b", "Product-graph generation type\n"
          "     0 - No opt\n"
          "     1 - Uncon rule\n"
          "     2 - Deg-1 rule\n"
          "     3 - Uncon & Deg-1 rule", gen_int, true)
      .refOption("x", "Molid 1", molid1, true)
      .refOption("y", "Molid 2", molid2, true)
      .refOption("i", "Input directory", in_path)
      .refOption("m", "Minimal core size", min_core)
      .refOption("s", "Shell size", shell);
  ap.parse();
  Product::GenerationType gen_type = static_cast<Product::GenerationType >(gen_int);

  Molecule mol1, mol2;
  std::stringstream fname1, fname2;
  fname1 << in_path << "/" << molid1 << ".lgf";
  fname2 << in_path << "/" << molid2 << ".lgf";

  std::ifstream ifs(fname1.str(), std::ifstream::in);
  mol1.read_lgf_stream(ifs);
  ifs.close();
  std::ifstream ifs2(fname2.str(), std::ifstream::in);
  mol2.read_lgf_stream(ifs2);
  ifs2.close();

  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  Product product(mol1, mol2, shell, gen_type, min_core);

  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<float> pg_duration = t2 - t1;

  std::vector<int> components, complete, bc_components, bc_complete;

  int pg_size = lemon::countNodes(product.get_graph());
  if (pg_size == 0)
    return 0;

  std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();

  float bk_time = 0;
  float gen_time = 0;
  FragmentVector fragments;
  IntVector degeneracy;
  if (product.get_components() == 1) {

    const Graph& g = product.get_graph();
    int V = lemon::countNodes(g);
    components.push_back(V);

    if (product.is_complete()) {
      complete.push_back(V);
    }

    Graph::EdgeMap<int> bc(g);
    int bcc = lemon::biNodeConnectedComponents(g, bc);
    if (bcc == 0) {
      bc_components.push_back(V);
      if (product.is_complete()) {
        bc_complete.push_back(V);
      }
    } else {
      for (int k = 0; k < bcc; ++k) {
        std::set<int> ids;
        int BCE = 0;
        for (EdgeIt e(g); e != lemon::INVALID; ++e) {
          if (bc[e] == k) {
            ids.insert(g.id(g.u(e)));
            ids.insert(g.id(g.v(e)));
            ++BCE;
          }
        }

        int BCV = ids.size();
        bc_components.push_back(BCV);
        if (BCE == ((BCV * (BCV - 1)) / 2)) {
          bc_complete.push_back(BCV);
        }
      }
    }

    std::tuple<float, float, bool, int> timing = run(product, min_core, max_core, fragments);
    bk_time = std::get<0>(timing);
    gen_time = std::get<1>(timing);
    degeneracy.push_back(std::get<3>(timing));
    assert(product.is_complete() == std::get<2>(timing) && [](int gen_int, int molid1, int molid2){
      std::cout << std::to_string(gen_int) << ": " << std::to_string(molid1) << "x" << std::to_string(molid2) << std::endl;
      return true;
    });
  } else {
    for (int c = 0; c < product.get_components(); ++c) {
      Product component(product, c);

      const Graph& g = component.get_graph();
      int V = lemon::countNodes(g);
      components.push_back(V);

      if (component.is_complete()) {
        complete.push_back(V);
      }

      Graph::EdgeMap<int> bc(g);
      int bcc = lemon::biNodeConnectedComponents(g, bc);
      if (bcc == 0) {
        bc_components.push_back(V);
        if (component.is_complete()) {
          bc_complete.push_back(V);
        }
      } else {
        for (int k = 0; k < bcc; ++k) {
          std::set<int> ids;
          int BCE = 0;
          for (EdgeIt e(g); e != lemon::INVALID; ++e) {
            if (bc[e] == k) {
              ids.insert(g.id(g.u(e)));
              ids.insert(g.id(g.v(e)));
              ++BCE;
            }
          }

          int BCV = ids.size();
          bc_components.push_back(BCV);
          if (BCE == ((BCV * (BCV - 1)) / 2)) {
            bc_complete.push_back(BCV);
          }
        }
      }

      std::tuple<float, float, bool, int> timing = run(component, min_core, max_core, fragments);
      bk_time += std::get<0>(timing);
      gen_time += std::get<1>(timing);
      degeneracy.push_back(std::get<3>(timing));
      assert(component.is_complete() == std::get<2>(timing) && [](int gen_int, int molid1, int molid2){
        std::cout << std::to_string(gen_int) << ": " << std::to_string(molid1) << "x" << std::to_string(molid2) << std::endl;
        return true;
      });
    }
  }

  if (fragments.empty())
    return 0;

  float pg_time = pg_duration.count();
  float runtime = pg_time + bk_time + gen_time;

  std::cout << pg_size << "\t" << fragments.size() << "\t"
            << pg_time << "\t" << bk_time << "\t"
            << gen_time << "\t" << runtime << "\t";
  print_int_list(degeneracy);
  print_int_list(components);
  print_int_list(complete);
  print_int_list(bc_components);
  print_int_list(bc_complete);
  std::cout << std::endl;

  return 0;

}
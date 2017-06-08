//
// Created by M. Engler on 26/10/16.
//

#include "../include/bronkerbosch.h"
#include <chrono>
#include <dirent.h>
#include <lemon/arg_parser.h>
#include <boost/algorithm/string/replace.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include "../include/fragment.h"
#include "../include/mcf.h"


using namespace mogli;

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

  Product product(mol1, mol2, shell, gen_type, min_core, max_core);

  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<float> pg_duration = t2 - t1;

  unsigned long pg_size = 0;
  for (NodeIt v(product.get_graph()); v != lemon::INVALID; ++v) {
    ++pg_size;
  }
  if (pg_size == 0)
    return 0;

  std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();

  BronKerbosch bk(product, min_core, max_core, false);
  bk.run();

  std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<float> bk_duration = t4 - t3;

  std::chrono::high_resolution_clock::time_point t5 = std::chrono::high_resolution_clock::now();

  NodeVectorVector cliques = bk.getMaxCliques();
  FragmentVector fragments;
  for (NodeVectorVector::const_iterator it = cliques.begin(), end = cliques.end(); it != end; ++it) {
    IntToIntMap g_to_mol1, g_to_mol2;
    boost::shared_ptr<Fragment> fragment = boost::make_shared<Fragment>(product, *it, g_to_mol1, g_to_mol2);

    if (fragment->get_core_atom_count() > 1) {
      fragments.push_back(fragment);
    }
  }

  std::chrono::high_resolution_clock::time_point t6 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<float> gen_duration = t6 - t5;

  if (fragments.size() == 0)
    return 0;


  std::cout << pg_size << "\t" << pg_duration.count() << "\t"
            << bk_duration.count() << "\t" << gen_duration.count() << std::endl;

  return 0;

}
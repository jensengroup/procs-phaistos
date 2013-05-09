// test_camshift.cpp --- Test of camshift energy class
// Copyright (C) 2011 Anders Steen Christensen
//
// This file is part of Phaistos
//
// Phaistos is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Phaistos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with Phaistos.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <string.h>

#include <boost/type_traits/is_base_of.hpp>

#include "protein/chain_fb.h"
#include "protein/definitions.h"
#include "protein/iterators/pair_iterator_chaintree.h"

#include "energy/energy.h"
#include "energy/procs_full.h"
#include "energy/observable.h"
#include "energy/observable_collection.h"

//using namespace definitions;
using namespace phaistos;

void test_procs(std::string pdb_filename, std::string bcs_filename) {

     // Create chain from PDB filename
     ChainFB chain(pdb_filename,definitions::ALL_ATOMS);
     // Add atoms missing in the pdb structure
     chain.add_atoms(definitions::ALL_PHYSICAL_ATOMS);

     // Create Energy class
     Energy<ChainFB> energy(&chain);

     // Add term
     //Term_procs_full::Settings settings_procs;
     Term_procs_full<ChainFB>::Settings settings_procs;
     settings_procs.weight = 1;
     settings_procs.bcsFilename = bcs_filename;
     energy.add_term(new Term_procs_full<ChainFB>(&chain, settings_procs));

     energy.evaluate();
     std::cout <<"\n\n"<< energy << "\n";
}


int main(int argc, char *argv[]) {


     if (argc < 3) {
          std::cerr << "Too few arguments -- confused and bailing out!" << std::endl;
          std::cerr << "Usage: ./test_procs myfile.bcs myfile.pdb" << std::endl;
          exit(1);
     }

     if (argc > 3) {
          std::cerr << "Too many arguments -- confused and bailing out!" << std::endl;
          std::cerr << "Usage: ./test_procs myfile.bcs myfile.pdb" << std::endl;
          exit(1);
     }


     std::string pdb_filename = argv[2];
     std::string bcs_filename = argv[1];

     test_procs(pdb_filename, bcs_filename);
}


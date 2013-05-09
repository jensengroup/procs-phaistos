// procs_full.h -- Backbone amide proton chemical shift energy function.
// Copyright (C) 2010 by Anders Christensen.
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

/*
	Short description:
	=====================================================

	This energy function compares experimental NMR chemical
	shifts for the backbone amide hydrogen to a chemical shift
	calculated using a fast QM based formula. The calculation
	is based on a backbone term, four bonds in the local 
	hydrogen bonding network and a ring-current term.

	The method is based upon the following papers:

	Parker, Houk and Jensen, 
	J. AM. CHEM. SOC. 2006, 128, 9863-9872

	with the backbone term replaced by that of 
	Eszter Czinki and Attila G. Csaszar,
	J. of Mol. Struc. 675 (2004) 107–116

	Bonds to non-amides are treated using the method of
	Michael Barfield, J. AM. CHEM. SOC. 2002, 124, 4158-4168	

	The backbone term is based on the point-dipol formulation 
	due to Pople but with new accurate parameters.
	(not published yet).

	NOTES: You must supply a chemical shift file containing 
	the appropriate chemical shifts. Also, this energy function
	needs ALL_PHYSICAL_ATOMS in order to find correct bonds.


*/


#ifndef PROCS_FULL
#define PROCS_FULL

#include "energy/energy_term.h"
//#include "utils/sidechain_enums.h"
#include "procs_utils.h"

#include <iostream>
#include <vector>
#include <cmath>
//#include <readbcs.h>

namespace phaistos {

//const bool display_output = true;
const bool display_output = false;
const double Stdev_exp = 0.3;


template <typename CHAIN_TYPE>

//class Term_procs_full: public EnergyTerm<CHAIN_TYPE>  {
class Term_procs_full: public EnergyTermCommon<Term_procs_full<CHAIN_TYPE>, CHAIN_TYPE> {

     // For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<Term_procs_full<CHAIN_TYPE>,CHAIN_TYPE> EnergyTermCommon;     


private:

   double current_energy;
public:
   std::vector<Atom_shift> Shift_list;
    const class Settings: public EnergyTerm<CHAIN_TYPE>::Settings {
    public:

         std::string bcsFilename;
         Settings(std::string bcsFilename="")
              : bcsFilename(bcsFilename) {}

         friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
              o << "bcs-filename:" << settings.bcsFilename << "\n";
              o << static_cast<const typename EnergyTerm<CHAIN_TYPE>::Settings>(settings);
              return o;
         }                    
    } settings;

// Constructor
     Term_procs_full(CHAIN_TYPE *chain, const Settings &settings=Settings(),
                     RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "procs_full", settings, random_number_engine), 
            settings(settings) {

// This is a FullAtom energy: ChainFB CHAIN_TYPE required
      if ( ! boost::is_base_of<ChainFB,CHAIN_TYPE>::value ) {
         std::cerr << "This is a FullAtom energy: ChainFB CHAIN_TYPE required\n";
         assert(false);
      }

// Readin .bcs f
	Get_bcsfile(settings.bcsFilename, Shift_list);
      }

//copy constructor
     Term_procs_full(const Term_procs_full &other, 
                     RandomNumberEngine *random_number_engine,
                     int thread_index, CHAIN_TYPE *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            Shift_list(other.Shift_list),
//            it_settings(other.it_settings),
            settings(other.settings) {
     } 

double evaluate(MoveInfo *moveInfo=NULL) {

     // Import protein definitions (such as residue names)
     using namespace definitions;

   int debug_flag = this->settings.debug;
   if (debug_flag > 0)     std::cout << " ----------------==== BEGIN EXECUTING PROCS-FULL ====----------------" << std::endl;
   if (debug_flag > 10) std::cout << "BEGIN PROCS INIT" << std::endl;
//Make vector of aromatic rings
      std::vector<Class_Ring> Rings;
      std::vector<Data_Point> Results_list;
      InitRings(this->chain, Rings);
      int i1;
      int ishift = 0;

//Init correct atoms, so we don't need to iterate over residues more than once.
//Positive hydrogens
      std::vector<Atom*> AmideBB_H;
      std::vector<Atom*> AmideSC_H;
      std::vector<Atom*> AlcoholSC_H;
      std::vector<Atom*> Charged_H;
      std::vector<Atom*> OtherSC_H;
//Negative 
      std::vector<Atom*> CarbonylBB_O;
      std::vector<Atom*> AmideSC_O;
      std::vector<Atom*> CarboxylSC_O;
      std::vector<Atom*> AlcoholSC_O;
//Backbone O and H, and terminals

      for(ResidueIterator<CHAIN_TYPE> res1(*(this->chain)); !(res1).end(); ++res1) {
//         std::cout << "Loading residue:" << (*res1).residue_type << (*res1).index << " ..." << std::endl;
         if (res1->has_atom(H)) {
            if ((*res1).index > 0) {
               AmideBB_H.push_back((*res1)[H]);
            } else {
               Charged_H.push_back((*res1)[H]);
            }
         }
         if (res1->has_atom(H1)) {
            Charged_H.push_back((*res1)[H1]);
         }
         if (res1->has_atom(H2)) {
            Charged_H.push_back((*res1)[H2]);
         }
         if (res1->has_atom(H2)) {
            Charged_H.push_back((*res1)[H3]);
         }
         if (!(res1->has_atom(OXT))) {
            CarbonylBB_O.push_back((*res1)[O]);
         } else if (res1->has_atom(OXT)) {
            CarboxylSC_O.push_back((*res1)[O]);
            CarboxylSC_O.push_back((*res1)[OXT]);
         }
//Alcohols
         if ((*res1).residue_type == TYR) {
            AlcoholSC_H.push_back((*res1)[HH]);
            AlcoholSC_O.push_back((*res1)[OH]);
         } else if ((*res1).residue_type == THR) {
            AlcoholSC_H.push_back((*res1)[HG1]);
            AlcoholSC_O.push_back((*res1)[OG1]);
         } else if ((*res1).residue_type == SER) {
            AlcoholSC_H.push_back((*res1)[HG]);
            AlcoholSC_O.push_back((*res1)[OG]);
         }
//Amides
         else if ((*res1).residue_type == ASN) {
            AmideSC_H.push_back((*res1)[HD21]);
            AmideSC_H.push_back((*res1)[HD22]);
            AmideSC_O.push_back((*res1)[OD1]);
         } else if ((*res1).residue_type == GLN) {
            AmideSC_H.push_back((*res1)[HE21]);
            AmideSC_H.push_back((*res1)[HE22]);
            AmideSC_O.push_back((*res1)[OE1]);
//Carboxylic acids
         } else if ((*res1).residue_type == GLU) {
            CarboxylSC_O.push_back((*res1)[OE1]);
            CarboxylSC_O.push_back((*res1)[OE2]);
         } else if ((*res1).residue_type == ASP) {
            CarboxylSC_O.push_back((*res1)[OD1]);
            CarboxylSC_O.push_back((*res1)[OD2]);
//Positive charged
         } else if ((*res1).residue_type == HIS) {
            if ((res1->has_atom(HD1))&&(res1->has_atom(HE2))) {
               Charged_H.push_back((*res1)[HD1]);
               Charged_H.push_back((*res1)[HE2]);
            }
         } else if ((*res1).residue_type == LYS) {
            Charged_H.push_back((*res1)[HZ1]);
            Charged_H.push_back((*res1)[HZ2]);
            Charged_H.push_back((*res1)[HZ3]);
         } else if ((*res1).residue_type == ARG) {
            Charged_H.push_back((*res1)[HH11]);
            Charged_H.push_back((*res1)[HH21]);
            Charged_H.push_back((*res1)[HH12]);
            Charged_H.push_back((*res1)[HH22]);
//            OtherSC_H.push_back((*res1)[HE]);
//Other positive
//         } else if ((*res1).residue_type == TRP) {
//            OtherSC_H.push_back((*res1)[HE1]);
//            OtherSC_H.push_back((*res1)[HD1]);
         }
      }
//End of list making
//Do the main loop
      if (debug_flag > 10) std::cout << "DONE PROCS INIT" << std::endl;
      for(ResidueIterator<CHAIN_TYPE> res1(*(this->chain)); !res1.end(); ++res1) {
         if (debug_flag > 10) std::cout << "BEGIN BACKBONE" << std::endl;
         double dB = 0;
         double dRC = 0;
         double dP = 0;
         double dPe = 0;
         double dS = 0;
         double dSe = 0;
         int Primary_bond_type = 0; // 0 = No bond, 1 = Backbone oxygen, 2 = carboxylic acid, 3 = alcohol, 4 = amide sc
         i1=(*res1).index;
         if (display_output) std::cout << (*res1).residue_type << i1+1;
//Skip all residues that don't have an amide proton (and explicitly Proline), also the N-term residue (#1) is skipped since it doesn't hold an amide proton.
         if (((*res1).residue_type != PRO)&&(i1!=0)&&(res1->has_atom(H))) {
            int i2;
            Vector_3D H1_pos;
            H1_pos=(*res1)[H]->position;
// Backbone term         
            double phi = res1->get_phi();
            double psi = res1->get_psi();
            dB = dBackbone_Czinki(phi, psi);
//            std::cout << ":" << std::endl << "Backbone: " << dB << std::endl;
// Sum over ring current terms
            for (std::vector<Class_Ring>::iterator ring1 = Rings.begin(); ring1 < Rings.end(); ring1++) {
               Vector_3D R_Vector;
               R_Vector = H1_pos - ring1->XYZ;
               double R_Theta = calc_angle(H1_pos, ring1->XYZ, ring1->Direction);
//                  dRC = dRC + (1.0 - 3.0*std::pow(cos(R_Theta), 2))/std::pow(R_Vector.norm(), 3)*(ring1->Intensity)*22.34;
               dRC = dRC + (1.0 - 3.0*std::pow(cos(R_Theta), 2))/std::pow(R_Vector.norm(), 3)*(ring1->Intensity)*27.4689;
            }
            if (debug_flag > 10) std::cout << "END BACKBONE" << std::endl;
            if (debug_flag > 10) std::cout << "BEGIN FIND HB1" << std::endl;
//               std::cout << "Ring current:     " << dRC << std::endl;
//               std::cout << "Primary:          " << O2_HB->residue->residue_type << O2_HB->residue->index+1 << "    " << dP << std::endl;
//Find closest primary bonding residue (if any). Store the hydrogen bonding partner in the *O2_HB pointer.
//At this stage, we only take the closest hydrogen bonding partner, in case there are more than one.
            Atom *O2_HB;
//WORKAROUND OR FIX WARNING:
            O2_HB=(*res1)[O];
//            int Primary_bond_type = 0; // 0 = No bond, 1 = Backbone oxygen, 2 = carboxylic acid, 3 = alcohol, 4 = amide sc
            float closest_bond = HB_cutoff;
//Loop over backbone amides
//	    bool found_O2_HB = false;
            for (std::vector<Atom*>::iterator Atom_O = CarbonylBB_O.begin(); Atom_O<CarbonylBB_O.end(); ++Atom_O) {
//               std::cout << (((*Atom_O)->position)- H1_pos).norm() << std::endl;
               if( (((*Atom_O)->position)- H1_pos).norm() < closest_bond) {
                  i2 = (*Atom_O)->residue->index;
                  if ((i1 != i2) && (i1-1 != i2)) {
                     O2_HB = (*Atom_O);
                     closest_bond = (((*Atom_O)->position)- H1_pos).norm();
                     Primary_bond_type = 1;
                  }
               }
            }
//Loop over caarboxylic sidechains and c-term
            for (std::vector<Atom*>::iterator Atom_O = CarboxylSC_O.begin(); Atom_O<CarboxylSC_O.end(); ++Atom_O) {
               if( (((*Atom_O)->position)- H1_pos).norm() < closest_bond) {
                  O2_HB = (*Atom_O);
                  closest_bond = (((*Atom_O)->position)- H1_pos).norm();
                  Primary_bond_type = 2;
               }
            }
//Loop over alcohols
            for (std::vector<Atom*>::iterator Atom_O = AlcoholSC_O.begin(); Atom_O<AlcoholSC_O.end(); ++Atom_O) {
               if( (((*Atom_O)->position)- H1_pos).norm() < closest_bond) {
                  O2_HB = (*Atom_O);
                  closest_bond = (((*Atom_O)->position)- H1_pos).norm();
                  Primary_bond_type = 3;
               }
            }
//Loop over sidechain amides
            for (std::vector<Atom*>::iterator Atom_O = AmideSC_O.begin(); Atom_O<AmideSC_O.end(); ++Atom_O) {
               if( (((*Atom_O)->position)- H1_pos).norm() < closest_bond) {
                  O2_HB = (*Atom_O);
                  closest_bond = (((*Atom_O)->position)- H1_pos).norm();
                  Primary_bond_type = 4;
               }
            }

//Primary bond terms,
            if (debug_flag > 10) std::cout << "END FIND HB1" << std::endl;
//Case 0: No bond:
            if (Primary_bond_type == 0) {
               dP = 2.04;
//Case 1: Backbone amide
            } else if (Primary_bond_type == 1) {
               if (debug_flag > 100) std::cout << "BEGIN PRIMARY BONDTYPE 1" << std::endl;
               Vector_3D C2_pos, N3_pos, O2_pos;
               C2_pos=(*(O2_HB->residue))[C]->position;
               N3_pos=(*(O2_HB->residue->get_neighbour(+1)))[N]->position;
               O2_pos=O2_HB->position;
               double rOH = (H1_pos-O2_pos).norm();
               double rho = calc_dihedral(H1_pos, O2_pos, C2_pos, N3_pos);
               double theta = calc_angle(H1_pos, O2_pos, C2_pos);
               dP = d1HB_Barfield(rho, theta, rOH);
// Primary extended bond term this function updates the dPe term - only available for backbone-backbone hydrogen bonds, thus far
               Update_Primary_extended(CarboxylSC_O, CarbonylBB_O, AmideSC_O, AlcoholSC_O, O2_HB, dPe);
//Case 2: Carboxylic acid
               if (debug_flag > 10) std::cout << "END PRIMARY BONDTYPE 1" << std::endl;
            } else if (Primary_bond_type == 2) {
               if (debug_flag > 10) std::cout << "BEGIN PRIMARY BONDTYPE 2" << std::endl;
               Vector_3D O2_pos, C2_pos, C3_pos;
               O2_pos=O2_HB->position;
//WORKAROUND OR FIX WARNING:
               C2_pos=O2_HB->position;
               C3_pos=O2_HB->position;
               double rOH = (H1_pos-O2_pos).norm();
               if (debug_flag > 100) std::cout << "GET BONDTYPE 2 PARAMS" << std::endl;
               if ( ((*O2_HB).atom_type == OE1)||((*O2_HB).atom_type == OE2) ) {
                  C2_pos=(*(O2_HB->residue))[CD]->position;
                  C3_pos=(*(O2_HB->residue))[CG]->position;
               } else if ( ((*O2_HB).atom_type == OD1)||((*O2_HB).atom_type == OD2) ) {
                  C2_pos=(*(O2_HB->residue))[CG]->position;
                  C3_pos=(*(O2_HB->residue))[CB]->position;
               } else if ( ((*O2_HB).atom_type == O)||((*O2_HB).atom_type == OXT) ) {
                  C2_pos=(*(O2_HB->residue))[C]->position;
                  C3_pos=(*(O2_HB->residue))[CA]->position;
               }
               if (debug_flag > 100) std::cout << "GOT BONDTYPE 2 PARAMS, NOW CALCULATING dP" << std::endl;
               double rho = calc_dihedral(H1_pos, O2_pos, C2_pos, C3_pos);
               double theta = calc_angle(H1_pos, O2_pos, C2_pos);               
               dP = d1HB_Acetate(rho, theta, rOH);
//Case 3: Alcohol
               if (debug_flag > 10) std::cout << "END PRIMARY BONDTYPE 2" << std::endl;
            } else if (Primary_bond_type == 3) {
               if (debug_flag > 10) std::cout << "BEGIN PRIMARY BONDTYPE 3" << std::endl;
               Vector_3D O2_pos=Vector_3D(0.0, 0.0, 0.0), C2_pos=Vector_3D(0.0, 0.0, 0.0), H3_pos=Vector_3D(0.0, 0.0, 0.0);
               O2_pos=O2_HB->position;
               double rOH = (H1_pos-O2_pos).norm();
               if ((*O2_HB).atom_type == OH) {
                  C2_pos=(*(O2_HB->residue))[CZ]->position;
                  H3_pos=(*(O2_HB->residue))[HH]->position;
               } else if ((*O2_HB).atom_type == OG) {
                  C2_pos=(*(O2_HB->residue))[CB]->position;
                  H3_pos=(*(O2_HB->residue))[HG]->position;
               } else if ((*O2_HB).atom_type == OG1) {
                  C2_pos=(*(O2_HB->residue))[CB]->position;
                  H3_pos=(*(O2_HB->residue))[HG1]->position;
               }
               double rho = calc_dihedral(H1_pos, O2_pos, C2_pos, H3_pos);
               double theta = calc_angle(H1_pos, O2_pos, C2_pos);
               dP = d1HB_Methanol(rho, theta, rOH);
//Case 4: Amide side chain (GLN or ASN)
               if (debug_flag > 10) std::cout << "END PRIMARY BONDTYPE 3" << std::endl;
            } else if (Primary_bond_type == 4) {
               if (debug_flag > 10) std::cout << "BEGIN PRIMARY BONDTYPE 4" << std::endl;
               Vector_3D O2_pos, C2_pos, N3_pos;
               if( ((*O2_HB).residue)->residue_type == GLN ) {
                  O2_pos=O2_HB->position;
                  C2_pos=(*(O2_HB->residue))[CD]->position;
                  N3_pos=(*(O2_HB->residue))[NE2]->position;
               } else { 
                  O2_pos=O2_HB->position;
                  C2_pos=(*(O2_HB->residue))[CG]->position;
                  N3_pos=(*(O2_HB->residue))[ND2]->position;
               }
               double rOH = (H1_pos-O2_pos).norm();
               double rho = calc_dihedral(H1_pos, O2_pos, C2_pos, N3_pos);
               double theta = calc_angle(H1_pos, O2_pos, C2_pos);
               dP = d1HB_Barfield(rho, theta, rOH);
               if (debug_flag > 10) std::cout << "END PRIMARY BONDTYPE 4" << std::endl;


            }


            if (debug_flag > 10) std::cout << "BEGIN SECONDARY" << std::endl; //DEBUG -----------------------------
//DEBUG NOTATION
//            dP = 2.04;
//DEBUG NOTATION
// Secondary bond term
            Residue *ResidueBefore = (*res1).get_neighbour(-1);
            Atom *O0=(*ResidueBefore)[O];
            Update_Secondary(CarboxylSC_O, CarbonylBB_O, AmideSC_O, AlcoholSC_O, AmideBB_H, AmideSC_H, AlcoholSC_H, Charged_H, OtherSC_H, O0, dS, dSe);
            if (debug_flag > 10) std::cout << "END SECONDARY" << std::endl; //DEBUG -----------------------------
	 
         }
         if (debug_flag > 10)  std::cout << "BEGIN DATA GET" << std::endl; //DEBUG -----------------------------
         double dFinal = 0.0;
         if (dB > 0) {

// CURRENT VERSION OF THE ENERGY EXPRESSION. USING A QM-FITTED LINEAR CORRECTING TO THE NEW BACKBONE TERM
// AND A SLIGHTLY HIGHER B REFERENCE VALUE FOR THE RING CURRENT TERM (SEE UPCOMING PAPERS OR WIKI)

            dFinal = (dB+0.445)*0.878 + dRC + dP + dPe + dS + dSe;

// TEST-CASES FOR VARIOUS CORRECTIONS TO THE TOTAL ENERGY TERM
//            dFinal = dB + dRC + dP + dPe + dS + dSe;
//            dFinal = dB + dRC + dP + dPe + dS + dSe - log(dP + 1.0)*(dB - 5.36)/5.28 - 0.076;
	    
         }

         if (display_output) std::cout << "   Calc: " << dFinal;
         double dExp;
         int resExp;

         dExp = Shift_list[ishift].value;
         resExp = Shift_list[ishift].res;
         ishift = ishift +1;
// Check if the residue number of the .pdb file and .bcs file match. If not, let's cry for help.
         if (resExp != (i1+1)) {
            if (display_output) std::cout << "CRITICAL ERROR, MISMATCH IN DATA FILES!" << std::endl;
            break;
         }
         if (display_output) std::cout << "   Exp: " << dExp << "(" << resExp << ")  " << Primary_bond_type << "   "<<  dP << std::endl;

         if ((dExp > 0) && (dFinal > 0)){
            Data_Point Result;
            Result.Exp = dExp;
            Result.Calc = dFinal;
            //Estimations of uncertainties for the different bins of bond-types. Amide proton bonds to:
            if (Primary_bond_type == 0) Result.Stdev = 2.0; //Solvent exposed 2.0
            if (Primary_bond_type == 1) Result.Stdev = 0.3; //Backbone amide 0.3
            if (Primary_bond_type == 2) Result.Stdev = 0.8; //Side chain carboxylic acid  0.8
            if (Primary_bond_type == 3) Result.Stdev = 0.8; //Side chain alcohol 0.8
            if (Primary_bond_type == 4) Result.Stdev = 0.6; //Side chain amide 0.6

            Results_list.push_back(Result);
//            std::cout << dRC << "   " << dB << "   " << (dP + dPe + dS + dSe) << "   " << dFinal << "  " << dExp << std::endl; 
         }
//         std::cout << "Total:                          " << dFinal << std::endl;
//         std::cout << ".........................................................." << std::endl;
      }

      double Energy = 0;

      for (std::vector<Data_Point>::iterator this_point=Results_list.begin(); this_point <Results_list.end(); this_point++){
           double this_energy = 0.5*pow(this_point->Stdev,-2)*(this_point->Exp - this_point->Calc)*(this_point->Exp - this_point->Calc)+ 0.5*std::log(2 * M_PI * pow(this_point->Stdev,2));
         Energy = Energy + this_energy;
      }

      if (display_output) {

           double RMSD = 0;
           double Exp_mean = 0;
           double Calc_mean = 0;
           double Sum_err_sqr= 0;
     
     
           for (std::vector<Data_Point>::iterator this_point=Results_list.begin(); this_point <Results_list.end(); this_point++){
              Exp_mean = Exp_mean + this_point->Exp;
              Calc_mean = Calc_mean + this_point->Calc;
              RMSD = RMSD + (this_point->Exp - this_point->Calc)*(this_point->Exp - this_point->Calc);
              Sum_err_sqr = Sum_err_sqr + (this_point->Exp - this_point->Calc)*(this_point->Exp - this_point->Calc);
           }
     
           Exp_mean = Exp_mean/Results_list.size();
           Calc_mean = Calc_mean/Results_list.size();
           double Co_var = 0;
           double Exp_stdev = 0;
           double Calc_stdev = 0;
           for (std::vector<Data_Point>::iterator this_point=Results_list.begin(); this_point <Results_list.end(); this_point++){
              Co_var = Co_var + (this_point->Exp - Exp_mean)*(this_point->Calc - Calc_mean);
              Exp_stdev = Exp_stdev + (this_point->Exp - Exp_mean)*(this_point->Exp - Exp_mean);
              Calc_stdev = Calc_stdev + (this_point->Calc - Calc_mean)*(this_point->Calc - Calc_mean);
           }
           double Pearson = Co_var/pow((Exp_stdev*Calc_stdev), 0.5);
           RMSD = pow((RMSD/Results_list.size()), 0.5);
           if (debug_flag > 10) std::cout << "END DATA GET" << std::endl; //DEBUG -----------------------------
           if (display_output) std::cout << "Correlation:        r = " << Pearson << std::endl;
           if (display_output) std::cout << "Chemical shift:  RMSD = " << RMSD << " ppm " << std::endl;
           if (display_output) std::cout << "Energy:             E = " << Energy << std::endl;
           if (debug_flag > 0) std::cout << "----------------==== END EXECUTING PROCS-FULL ====----------------" << std::endl;

      }

      return Energy;
   }
};


}

#endif

namespace module_procs {

// Module: energy term initialization
struct EnergyInitialization {


     // Constructor - general case: do nothing
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     EnergyInitialization(const Options &options, CHAIN_TYPE *chain, DBN_TYPE *dbn,
                                       Energy<CHAIN_TYPE> *energy, std::vector<RandomNumberEngine *> *random_number_generators, 
                                       std::string prefix="") {
     }

     // Constructor - template specific case
     template <typename DBN_TYPE>
     EnergyInitialization(const Options &options, ChainFB *chain, DBN_TYPE *dbn,
                                       Energy<ChainFB> *energy, std::vector<RandomNumberEngine *> *random_number_generators, 
                                       std::string prefix="") {

          Options::OptionValue option;

          // ProCS
          option = options[prefix+"-procs"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef typename Term_procs_full<ChainFB>::Settings Settings;

               // Add energy term
               energy->add_term(new Term_procs_full<ChainFB>(chain,
                                                               options.get_settings<Settings>(option,i)));
          }
     }

};


}

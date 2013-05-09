namespace module_procs {

// Module: energy term initialization
struct EnergyOptions {

     // Constructor - general case: do nothing
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   CHAIN_TYPE *chain,
                   DBN_TYPE *dbn) {                    
     }

     // Constructor - ChainFB specific case
     template <typename DBN_TYPE>
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   ChainFB *chain,
                   DBN_TYPE *dbn) {

          // Import namespace for make_vector
          using namespace boost::fusion;


          // ProCS
          for (int counter = occurrences[prefix+"-procs"]; counter > 0; counter--) {

               typedef typename Term_procs_full<ChainFB>::Settings Settings;
               boost::shared_ptr<Settings> settings(new Settings());

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "procs (" + prefix + ")",
                         prefix+"-procs", settings,
                         make_vector(
                              make_vector(std::string("bcs-filename"),
                                          std::string("path to file containing ???"),
                                          &settings->bcsFilename)
                              )),
                    super_group, counter==1);
          }
    }
};

}

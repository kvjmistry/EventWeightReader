////////////////////////////////////////////////////////////////////////
// Class:       EventWeightReader
// Plugin Type: analyzer (art v2_05_01)
// File:        EventWeightReader_module.cc
//
// Generated at Wed Aug 22 05:28:06 2018 by Krishan Mistry using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

// Default art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
// #include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

// LARSOFT includes
#include "lardataobj/Simulation/SimPhotons.h" 
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h" 
#include "larsim/Simulation/LArG4Parameters.h"                          
// #include "larcore/Geometry/Geometry.h" 
#include "lardataobj/RecoBase/Hit.h" 
#include "nutools/ParticleNavigation/ParticleList.h" 
#include "nutools/ParticleNavigation/EmEveIdCalculator.h" 
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" 
#include "nusimdata/SimulationBase/MCTruth.h"

#include "uboone/EventWeight/MCEventWeight.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// ROOT includes
#include "TH1.h" 
#include "TH2.h" 
#include "TH3.h" 
#include "TTree.h" 
#include "TDatabasePDG.h" 
#include "TParticlePDG.h" 
#include "TCanvas.h" 
#include "TVectorT.h" 
#include "TMatrixT.h" 
#include "TMatrixDUtils.h" 
#include "TGraph.h" 
#include "TF1.h" 
#include "TMath.h"


class EventWeightReader;


class EventWeightReader : public art::EDAnalyzer {
public:
  explicit EventWeightReader(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  EventWeightReader(EventWeightReader const &) = delete;
  EventWeightReader(EventWeightReader &&) = delete;
  EventWeightReader & operator = (EventWeightReader const &) = delete;
  EventWeightReader & operator = (EventWeightReader &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void AddWeights(std::vector<double> N, int Iterations, int Universes);
  

private:

  // Declare member data here.
  // TTree* DataTree;
  // TTree* KeepTree;
  std::map<std::string, std::vector<double> > weights; // Map (Model Name, Weight Vector in each universe)
  int run, subrun, evt;

  // TH1D *hInterestingWeights;
  // TH1D *hDiscardWeights;

  std::vector<double> WeightList;           // List of all weights in one vector. 

  std::vector<std::string> KeepProcess;     // Genie processes that are non-zero
  std::vector<std::string> DiscardProcess;  // Processes that have no effect. 

  int Iterations{0};                        // Number of times looped event to see how quickly module is running

  const int Universes{10};                  // The number of universes simulated

  std::ofstream Genie_Weights_file_wNames;  // file with genie processes and names

  // Vectors for cross section calculation. 
  std::vector<double> N_gen, N_sig, MC_x_sec;
  std::vector<int> N_gen_evt, N_sig_evt; // vectors of event numbers for signal, selected and gen

  // Flux
  const double flux{4.19844e+10};

  // Num Targets
  const double targets{3.50191e+31};


};

// A function that loops over all the parameter weights and universes and re-weights the desired events. 
EventWeightReader::AddWeights(std::vector<double> N, int Iterations, int Universes){

  auto GenieEW_Handle = e.getValidHandle<std::vector<evwgh::MCEventWeight>>("mcweight"); // Request the mcweight data product

  if(GenieEW_Handle.isValid()) {
      std::cout << "[Analyze] GenieEW_Handle is valid" << std::endl; 
  }

  std::vector<evwgh::MCEventWeight> const& GenieEWvec(*GenieEW_Handle);

  for (evwgh::MCEventWeight const& GenieEW: GenieEWvec) { // Loop over weight handles.
    
    weights = GenieEW.fWeight; // Grab the weights 
  }

  // Initialise the size of the counter if it is the first event loop. 
  if (Iteration == 1) { N.resize( weights.size() * Universes ) }; // Resize to number of parameters * Universes. 
  
  // Loop over the weights and print their name and values. 
  for (auto const& it : weights) {
    GenieNames = it.first; 
    //std::cout << "\n" << GenieNames << std::endl;
    
    // Loop over each universe
    for (unsigned int i = 0; i < it.second.size(); i++){ 
      
      //std::cout << it.second[i] << "\t";
      WeightList.push_back(it.second[i]); // Add weights to a vector

      N[it] += it.second[i]; // Add weight to vector of counters.

      // Fill histograms and fill vectors which have non 1 and 1 values 
      if (it.second[i] == 1.0){
        
        // hDiscardWeights->Fill(GenieNames,1); // Reject
        
        // look to see if already found string in vector
        if (std::find(DiscardProcess.begin(), DiscardProcess.end(), it.first) != DiscardProcess.end()){}
        else { DiscardProcess.push_back(it.first); }
      
      } // End discard condition
      
      else {
        
        // hInterestingWeights->Fill(GenieNames,1); // Pass
        
        if (std::find(KeepProcess.begin(), KeepProcess.end(), it.first) != KeepProcess.end()){} // look to see if already found string in vector
        else { KeepProcess.push_back(it.first); }
           
      } // end Keep condition

    } // loop over each universe
    
    // std::cout << std::endl;
    
  } // END loop over weights 


}



EventWeightReader::EventWeightReader(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void EventWeightReader::beginJob()
{
  // Implementation of optional member function here.
  // Access ART's TFileService, which will handle histograms/trees/etc.
	// art::ServiceHandle<art::TFileService> tfs;

  // hInterestingWeights = 		tfs->make<TH1D>("hInterestingWeights", "Weights to Keep",	20,  0., 10.);
  // hDiscardWeights     = 		tfs->make<TH1D>("hDiscardWeights",     "Weights to Discard",		20  ,0., 10.);

  // Create the TTree and add relavent branches
	// DataTree = tfs->make<TTree>("EventTree","EventTree"); 
  // KeepTree = tfs->make<TTree>("KeepTree","KeepTree"); 

  // Event Information
	// DataTree->Branch("run", &run);
	// DataTree->Branch("subrun",&subrun);
  // DataTree->Branch("event",&evt);
  // DataTree->Branch("weights",&weights);
  // KeepTree->Branch("keep",&KeepProcess);
  // KeepTree->Branch("discard",&DiscardProcess);

}
void EventWeightReader::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  // Determine event ID, run and subrun 
  run = e.id().run();
  subrun = e.id().subRun();
  evt = e.id().event();
  
  TString GenieNames; // Temp string for displaying genie names

  std::cout << "++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "Iteration\t" << Iterations<< std::endl;
  std::cout << "++++++++++++++++++++++++++++++++++" << std::endl;
  Iterations++;


  // Load in the file containing the event number for Signal_Generated, N_gen or Signal_Selected N_sig
  
  // ++++++++++++++++++++ N_gen +++++++++++++++++++++++++++++
  std::ifstream fN_gen;
	fN_gen.open("filename_events.txt");
	
  if (!fN_gen.good()) { // Check if the file opened correctly
		cerr << "Error: Gen file could not be opened" << endl;
    exit(1);
	}

    double temp{0};

    if (fN_gen.is_open()) {
      while ( !fN_gen.eof()) {
        fN_gen >> temp;
        N_gen_evt.push_back(temp);
    	}
	    fN_gen.close();
  	}
    // ++++++++++++++++++++ N_sig +++++++++++++++++++++++++++++
    std::ifstream fN_sig;
  	fN_sig.open("filename_events.txt");
	
    if (!fN_sig.good()) { // Check if the file opened correctly
		cerr << "Error: sig file could not be opened" << endl;
    exit(1);
	  }

    double temp{0};

    if (fN_sig.is_open()) {
      while ( !fN_sig.eof()) {
        fN_sig >> temp;
        N_sig_evt.push_back(temp);
    	}
	    fN_sig.close();
  	}

  
  // Choose wheather to populate which varible depending on the event number
  if ( std::find( N_gen_evt.begin(), N_gen_evt.end(), evt) != N_gen_evt.end()){
    
    // Found event in the N_gen vector and so add weights for this event
    AddWeights(N_gen, Iterations, Universes);
  }
  else if ( std::find( N_sig_evt.begin(), N_sig_evt.end(), evt) != N_sig_evt.end()){
    
    // Found event in the N_sig vector nd so add weights for this event
    AddWeights(N_sig, Iterations, Universes);

  }

// DataTree->Fill();

}

void EventWeightReader::endJob()
{
  // Implementation of optional member function here.
  // KeepTree->Fill();

  // Print out which pricesses to keep and discard
  std::cout << "\n=======================" << std::endl;
  std::cout << "Keep genie Processes\n" << std::endl;
  std::cout << "=======================\n" << std::endl;
  for (unsigned int i = 0; i< KeepProcess.size(); i++){
    std::cout << KeepProcess[i] << std::endl;
  }

  std::cout << "\n=======================" << std::endl;
  std::cout << "Discard genie Processes\n" << std::endl;
  std::cout << "=======================\n" << std::endl;
  for (unsigned int i = 0; i< DiscardProcess.size(); i++){
    std::cout << DiscardProcess[i] << std::endl;
  }


  // Open a file with the new x section values in
  std::ofstream MC_weighted_xsec_file;
  MC_weighted_xsec_file.open("MC_weighted_xsec_file.txt");
  
  MC_x_sec.resize(N_gen.size()); // Resize

  // Calculate the New Cross section. 
  for (unsigned int i{0}; i < N_gen.size(); i++){

    MC_x_sec[i] =  N_gen[i] / ( flux * targets);
    std::cout << "New X-sections [10^-39 cm^2]\t" << MC_x_sec[i] << std::endl;

    MC_weighted_xsec_file << MC_x_sec[i] << "\n"; // Put the new cross sections into a text file. 
  }


}

DEFINE_ART_MODULE(EventWeightReader)

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

#include <iostream>
#include <fstream>


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
  void AddWeights(std::vector<double> N, int Iterations, int Universes, art::Event const & e);
  void ReadEvents(const char *filename, std::vector<int> N_evt );

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
  std::vector<double> N_gen, N_sig, N_bkd, N_sel ,MC_x_sec;
  std::vector<int> N_gen_evt, N_sig_evt, N_bkd_evt, N_sel_evt; // Vectors of event numbers for generated, signal, background, selected

  // Flux
  const double flux{4.19844e+10};

  // Num Targets
  const double targets{3.50191e+31};

};

// A function that loops over all the parameter weights and universes and re-weights the desired events. 
void EventWeightReader::AddWeights(std::vector<double> N, int Iterations, int Universes, art::Event const & e){
  TString GenieNames; // Temp string for displaying genie names

  auto GenieEW_Handle = e.getValidHandle<std::vector<evwgh::MCEventWeight>>("mcweight"); // Request the mcweight data product

  if(GenieEW_Handle.isValid()) {
      std::cout << "[Analyze] GenieEW_Handle is valid" << std::endl; 
  }

  std::vector<evwgh::MCEventWeight> const& GenieEWvec(*GenieEW_Handle);

  for (evwgh::MCEventWeight const& GenieEW: GenieEWvec) { // Loop over weight handles.
    
    weights = GenieEW.fWeight; // Grab the weights 
  }

  // Initialise the size of the counter if it is the first event loop. 
  if (Iterations == 1) { N.resize( weights.size() * Universes ); } // Resize to number of parameters * Universes. 
  
  int loop_counter{0};
  // Loop over the weights and print their name and values. 
  for (auto const& it : weights) {
    GenieNames = it.first; 
    //std::cout << "\n" << GenieNames << std::endl;
    
    // Loop over each universe
    for (unsigned int i = 0; i < it.second.size(); i++){ 
      
      //std::cout << it.second[i] << "\t";
      WeightList.push_back(it.second[i]); // Add weights to a vector

      N[ loop_counter + i ] += it.second[i]; // Add weight to vector of counters.

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
    loop_counter+=Universes;
  } // END loop over weights 


}

// Function that reads in the event numbers from a text file and adds those event numbers to a vector. 
void EventWeightReader::ReadEvents(const char *filename, std::vector<int> N_evt ){
	
  std::ifstream fileIN; 

  fileIN.open(filename); // Open the file
	
  if (!fileIN.good()) { // Check if the file opened correctly
		std::cerr << "Error: file:\t" << filename <<"\tcould not be opened" << std::endl;
    exit(1);
	}

    double temp{0}; // Use a temp var to get the values and push back

    if (fileIN.is_open()) { 
      
      while ( !fileIN.eof()) { // loop over lines in file
        
        fileIN >> temp;        // Add number to temp var
        N_evt.push_back(temp);
    	}
	    
      fileIN.close();
  	}

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

  std::cout << "++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "Iteration\t" << Iterations<< std::endl;
  std::cout << "++++++++++++++++++++++++++++++++++" << std::endl;
  Iterations++;


  // Load in the file containing the event number for Signal_Generated, N_gen or Signal_Selected N_sig
  
  // ++++++++++++++++++++ N_gen +++++++++++++++++++++++++++++
  
  ReadEvents("Gen_events.txt", N_gen_evt ); 

  // ++++++++++++++++++++ N_sig +++++++++++++++++++++++++++++
  
  ReadEvents("Sig_events.txt", N_sig_evt ); 

  // ++++++++++++++++++++ N_sel +++++++++++++++++++++++++++++
  
  ReadEvents("Sel_events.txt", N_sel_evt ); 

  // ++++++++++++++++++++ N_bkd +++++++++++++++++++++++++++++
  
  ReadEvents("Bkd_events.txt", N_bkd_evt ); 

  
  // Choose wheather to populate which varible depending on the event number
  
  // ++++++++++++++++++++ N_gen +++++++++++++++++++++++++++++
  if ( std::find( N_gen_evt.begin(), N_gen_evt.end(), evt) != N_gen_evt.end()){
    
    // Found event in the N_gen vector and so add weights for this event
    AddWeights(N_gen, Iterations, Universes, e);
  }
  // ++++++++++++++++++++ N_sig +++++++++++++++++++++++++++++
  else if ( std::find( N_sig_evt.begin(), N_sig_evt.end(), evt) != N_sig_evt.end()){
    
    // Found event in the N_sig vector and so add weights for this event
    AddWeights(N_sig, Iterations, Universes, e);
  }
  // ++++++++++++++++++++ N_sel +++++++++++++++++++++++++++++
  else if ( std::find( N_sel_evt.begin(), N_sel_evt.end(), evt) != N_sel_evt.end()){
    
    // Found event in the N_sel vector and so add weights for this event
    AddWeights(N_sel, Iterations, Universes, e);
  }
  // ++++++++++++++++++++ N_bkd +++++++++++++++++++++++++++++
  else if ( std::find( N_bkd_evt.begin(), N_bkd_evt.end(), evt) != N_bkd_evt.end()){
    
    // Found event in the N_bkd vector and so add weights for this event
    AddWeights(N_bkd, Iterations, Universes, e);
  }
  else {
    std::cout << "\n+++++++++++++++++++++++" << std::endl;
    std::cout << "Event not found..." << std::endl;
    std::cout << "\n+++++++++++++++++++++++" << std::endl;
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
    double efficiency = N_sig[i] / N_gen[i]; 
    std::cout << "+++++++\n efficiency\t" << efficiency << std::endl;

    MC_x_sec[i] =  (N_sel[i] - N_bkd[i]) / ( efficiency * flux * targets);
    std::cout << "New X-sections [10^-39 cm^2]\t" << MC_x_sec[i]/1e-39 << std::endl;

    MC_weighted_xsec_file << MC_x_sec[i] << "\n"; // Put the new cross sections into a text file. 
  }


}

DEFINE_ART_MODULE(EventWeightReader)

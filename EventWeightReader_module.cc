////////////////////////////////////////////////////////////////////////
// Class:       EventWeightReader
// Plugin Type: analyzer (art v2_05_01)
// File:        EventWeightReader_module.cc
//
// Generated at Wed Aug 22 05:28:06 2018 by Krishan Mistry using cetskelgen
// from cetlib version v1_21_00.
//
// This LArsoft module will read in a set of weights from a root file input
// and apply these weights on an event by event basis. It will then calculate
// a new value of the cross section and put the new values into a txt
// file which can be used to calculate the systematic error.
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
#include "TString.h"

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
  void AddWeights(std::vector<double> &N, int &Iterations, int Universes, art::Event const & e);
  void ReadEvents(const char *filename, std::vector<int> &N_evt );

private:

  // Declare member data here.
  int run, subrun, evt;
  std::map<std::string, std::vector<double> > weights; // Map (Model Name, Weight Vector in each universe)
  
  std::vector<double> WeightList;           // List of all weights in one vector. 

  std::vector<std::string> KeepProcess;     // Genie processes that are non-zero
  std::vector<std::string> DiscardProcess;  // Processes that have no effect. 

  int Iterations{1};                        // Number of times looped event to see how quickly module is running

  const int Universes{100};                   // The number of universes simulated

  double total_in{0}; // Total number of matched events 
  double tot_gen{0};  // counter for the total number of gen events read in
  double tot_sel{0};  // counter for the total number of sel events read in
  double tot_sig{0};  // counter for the total number of sig events read in
  double tot_bkg{0};  // counter for the total number of bkg events read in
  double tot_filt{0};  // counter for the total number of filtered events read in
  double tot_unmatched{0};  // counter for the total number of filtered events read in

  // Vectors for cross section calculation. 
  std::vector<double> N_gen, N_sig, N_bkg, N_sel ,MC_x_sec, Data_x_sec, Efficiency; // Vectors of num of events with new weights for generated, signal, background, selected
  std::vector<int> N_gen_evt, N_sig_evt, N_bkg_evt, N_sel_evt, N_filt;          // Vectors of event numbers from infile for generated, signal, background, selected, filtinlist
  std::vector<std::string> GenieNames;
  // Flux
  const double flux_mc{};
  const double flux_data{};

  // Num Targets
  const double targets_mc{};
  const double targets_data{};
  

  // DATA
  const double intime_cosmics_bkg{};              // Number of intime cosmics for background
  const double num_selected_data{};              // The number of selected events in data
  const double intime_cosmic_scale_factor{}; // Scale factor to apply to the intime cosimic background
  const double mc_scale_factor{};               // Scale factor to apply to the mc background

  // DEBUG
  bool DEBUG{true};

  // Histograms
  TGraph* gModelXsec;             		              // Graph of the Model vs cross section

  // TTree
  TTree *DataTree;

};

// A function that loops over all the parameter weights and universes and re-weights the desired events. 
void EventWeightReader::AddWeights(std::vector<double> &N, int &Iterations, int Universes, art::Event const & e){

  auto GenieEW_Handle = e.getValidHandle<std::vector<evwgh::MCEventWeight>>("mcweight"); // Request the mcweight data product

  if(GenieEW_Handle.isValid()) {
      std::cout << "[Analyze] GenieEW_Handle is valid" << std::endl; 
  }

  std::vector<evwgh::MCEventWeight> const& GenieEWvec(*GenieEW_Handle);

  for (evwgh::MCEventWeight const& GenieEW: GenieEWvec) { // Loop over weight handles.
    
    weights = GenieEW.fWeight; // Grab the weights 
  }

  // Initialise the size of the counter if it is the first event loop. 
  if (N.size() == 0 ) { N.resize( weights.size() * Universes ); } // Resize to number of parameters * Universes. 

  int loop_counter = 0;

  // Loop over the parameter models . 
  for (auto const& it : weights) {
    
    // Loop over each universe for a parameter
    for (unsigned int i = 0; i < it.second.size(); i++){ 
      GenieNames.push_back(it.first); 
      //std::cout <<"Universe:\t " << GenieNames <<"\tindex:\t" << loop_counter + i  << std::endl;
      
      //std::cout << it.second[i] << "\t";
      WeightList.push_back(it.second[i]); // Add weights to a vector

      N[ loop_counter + i ] += it.second[i]; // Add weight to vector of counters.
      // std::cout << "Weight\t"<< it.second[i] <<std::endl;

      // Fill vectors which have non 1 and 1 values to see what processes we can discard. 
      if (it.second[i] == 1.0){
        
        // look to see if already found string in vector
        if (std::find(DiscardProcess.begin(), DiscardProcess.end(), it.first) != DiscardProcess.end()){}
        else { DiscardProcess.push_back(it.first); }
      
      } // End discard condition
      
      else {
        
        if (std::find(KeepProcess.begin(), KeepProcess.end(), it.first) != KeepProcess.end()){} // look to see if already found string in vector
        else { KeepProcess.push_back(it.first); }
           
      } // end Keep condition

    } // loop over each universe
    
    // std::cout << std::endl;
    loop_counter += Universes;

  } // END loop over parameters
}

// Function that reads in the event numbers from a text file and adds those event numbers to a vector. 
void EventWeightReader::ReadEvents(const char *filename, std::vector<int> &N_evt ){
	
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
	art::ServiceHandle<art::TFileService> tfs;

  // Energy Calibration scatter plot
	gModelXsec = tfs->makeAndRegister<TGraph>("gModelXsec", "Universe Cross sections;Model; Re-calculated Cross Section [cm^2]");
	gModelXsec->SetMarkerStyle(kFullDotMedium);
	gModelXsec->SetMarkerSize(6);
	gModelXsec->SetLineWidth(0);

  // Create the TTree and add relavent branches
  DataTree = tfs->make<TTree>("XSectionTree","XSectionTree");
  DataTree->Branch("Sig",         &N_sig);
  DataTree->Branch("Bkg",         &N_bkg);
  DataTree->Branch("Sel",         &N_sel);
  DataTree->Branch("Gen",         &N_gen);
  DataTree->Branch("eff",         &Efficiency);
  DataTree->Branch("MCXSec",      &MC_x_sec);
  DataTree->Branch("DataXSec",    &Data_x_sec);
  DataTree->Branch("GenieNames",  &GenieNames);

  // Load in the file containing the event number for Signal_Generated, N_gen or Signal_Selected N_sig
  if (DEBUG) std::cout << "\nNow reading in event files!" << std::endl;

  // ++++++++++++++++++++ N_gen +++++++++++++++++++++++++++++
  
  ReadEvents("Gen_events.txt", N_gen_evt ); 

  // ++++++++++++++++++++ N_sig +++++++++++++++++++++++++++++
  
  ReadEvents("Sig_events.txt", N_sig_evt ); 

  // ++++++++++++++++++++ N_sel +++++++++++++++++++++++++++++
  
  ReadEvents("Sel_events.txt", N_sel_evt ); 

  // ++++++++++++++++++++ N_bkg +++++++++++++++++++++++++++++
  
  ReadEvents("Bkg_events.txt", N_bkg_evt ); 

  // ++++++++++++++++++++ Filtered list total +++++++++++++++++++++++++++++
  
  ReadEvents("FilteredList.txt", N_filt ); 


  if (DEBUG) std::cout << "Size of Generated vector:  \t"<< N_gen_evt.size() - 1 << std::endl;
  if (DEBUG) std::cout << "Size of signal vector:     \t"<< N_sig_evt.size() - 1 << std::endl;
  if (DEBUG) std::cout << "Size of selected vector:   \t"<< N_sel_evt.size() - 1 << std::endl;
  if (DEBUG) std::cout << "Size of background vector: \t"<< N_bkg_evt.size() - 1 << std::endl;
  if (DEBUG) std::cout << "Size of input Filtered List: \t"<< N_filt.size() - 1 << std::endl;

  if (DEBUG) std::cout << "\nFinished reading in event files!\n" << std::endl;

  std::cout << "++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "WARNING!! Remember to change number of universes to correct amount otherwise this will segfault :("<<  std::endl;
  std::cout << "++++++++++++++++++++++++++++++++++" << std::endl;


}
void EventWeightReader::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  // Determine event ID, run and subrun 
  run =     e.id().run();
  subrun =  e.id().subRun();
  evt =     e.id().event();

  std::cout << "++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "Event processed:\t" << Iterations<< std::endl;
  std::cout << "++++++++++++++++++++++++++++++++++" << std::endl;
  Iterations++;

  // Choose wheather to populate which varible depending on the event number

  // ++++++++++++++++++++ N_gen +++++++++++++++++++++++++++++
  if ( std::find( N_gen_evt.begin(), N_gen_evt.end(), evt) != N_gen_evt.end()){
    
    if (DEBUG) std::cout << "Matched a Generated Event\t" << std::endl;
    tot_gen++;

    // Found event in the N_gen vector and so add weights for this event
    AddWeights(N_gen, Iterations, Universes, e);
  }
  // ++++++++++++++++++++ N_sel +++++++++++++++++++++++++++++
  if ( std::find( N_sel_evt.begin(), N_sel_evt.end(), evt) != N_sel_evt.end()){
    
    if (DEBUG) std::cout << "Matched a Selected Event\t" << std::endl;

    tot_sel++;

    // Found event in the N_sel vector and so add weights for this event
    AddWeights(N_sel, Iterations, Universes, e);
  }
  // ++++++++++++++++++++ N_sig +++++++++++++++++++++++++++++
  if ( std::find( N_sig_evt.begin(), N_sig_evt.end(), evt) != N_sig_evt.end()){
    
    if (DEBUG) std::cout << "Matched a Signal Selected Event\t" << std::endl;

    tot_sig++;

    // Found event in the N_sig vector and so add weights for this event
    AddWeights(N_sig, Iterations, Universes, e);
  }
  // ++++++++++++++++++++ N_bkg +++++++++++++++++++++++++++++
  if ( std::find( N_bkg_evt.begin(), N_bkg_evt.end(), evt) != N_bkg_evt.end()){
    
    if (DEBUG) std::cout << "Matched a Background Selected Event\t" << std::endl;

    tot_bkg++;

    // Found event in the N_bkg vector and so add weights for this event
    AddWeights(N_bkg, Iterations, Universes, e);
  }

  // ++++++++++++++++++++ Filtered list check +++++++++++++++++++++++++++++
  if ( std::find( N_filt.begin(), N_filt.end(), evt) != N_filt.end()){
    
    if (DEBUG) std::cout << "Matched a Filtered Event\t" << std::endl;

    tot_filt++;
  }
  else{
    if (DEBUG) std::cout << "Unmatched Event\t" << std::endl;
    tot_unmatched++;
  }
  

total_in = tot_gen + tot_bkg;
  
}

void EventWeightReader::endJob()
{

  // Implementation of optional member function here.
  std::cout << "\nBeginning END JOB..." << std::endl;
  
  // // Print out which pricesses to keep and discard
  // std::cout << "\n=======================" << std::endl;
  // std::cout << "Keep genie Processes\n" << std::endl;
  // std::cout << "=======================\n" << std::endl;
  // for (unsigned int i = 0; i< KeepProcess.size(); i++){
  //   std::cout << KeepProcess[i] << std::endl;
  // }

  // std::cout << "\n=======================" << std::endl;
  // std::cout << "Discard genie Processes\n" << std::endl;
  // std::cout << "=======================\n" << std::endl;
  // for (unsigned int i = 0; i< DiscardProcess.size(); i++){
  //   std::cout << DiscardProcess[i] << std::endl;
  // }

  std::cout << "\n=======================" << std::endl;
  std::cout << "Read in Events\n" << std::endl;
  std::cout << "=======================\n" << std::endl;
  std::cout << "\nN_gen:\t" << tot_gen << std::endl;
  std::cout << "\nN_sel:\t" << tot_sel << std::endl;
  std::cout << "\nN_sig:\t" << tot_sig << std::endl;
  std::cout << "\nN_bkg:\t" << tot_bkg << std::endl;
  std::cout << "\nTotal:\t" << total_in << std::endl;
  std::cout << "\nTotal Filtered Events:\t" << tot_filt << std::endl;
  std::cout << "\nTotal Unmatched Events:\t" << tot_unmatched << std::endl;

  std::cout << "\n=======================" << std::endl;
  std::cout << "CV Cross sections\n" << std::endl;
  std::cout << "=======================\n" << std::endl;
  std::cout << "\nMC:\t" << (tot_sel + 2 - tot_bkg) / ( (tot_sig + 2 ) / (tot_gen + 4) * flux_mc * targets_mc) << std::endl;
  std::cout << "\nData:\t" << ( num_selected_data - ((tot_bkg) * mc_scale_factor + intime_cosmics_bkg * intime_cosmic_scale_factor )) / ( ((tot_sig + 2)/(tot_gen + 4 )) * flux_data * targets_data) << std::endl;

  std::cout << "\n=======================" << std::endl;
  std::cout << "Now Re-Calculating cross sections" << std::endl;
  std::cout << "=======================\n" << std::endl;

  // Open a file with the new x section values in
  std::ofstream MC_weighted_xsec_file;
  MC_weighted_xsec_file.open("MC_weighted_xsec_file.txt");

  std::ofstream Data_weighted_xsec_file;
  Data_weighted_xsec_file.open("Data_weighted_xsec_file.txt");

  MC_x_sec.resize(N_gen.size()); // Resize
  Data_x_sec.resize(N_gen.size());
  Efficiency.resize(N_gen.size());

  // Calculate the new Cross section. 
  for (unsigned int i{0}; i < N_gen.size(); i++){

    // Recalculations due to not weighting non MC genie stuff and other bugs
    N_gen[i] = N_gen[i] + 4  ;  // The plus 4 is for the tpc obj of size zero bug
    N_sel[i] = tot_sel  + 2; ;  // ANDY F: do not reweight the selected events in MC so it is like data
    N_sig[i] = N_sig[i] + 2;    // Missing two events from somewhere
    N_bkg[i] = N_bkg[i];

    Efficiency[i] = N_sig[i] / N_gen[i];  // 0.0884133 CV efficiency
    if (DEBUG) std::cout << "\n+++++++\nEfficiency\t" << Efficiency[i] << std::endl;

    std::cout << "\nN_gen:\t" << N_gen[i] << std::endl;
    std::cout << "N_sel:\t" << N_sel[i] << std::endl;
    std::cout << "N_sig:\t" << N_sig[i] << std::endl;
    std::cout << "N_bkg:\t" << N_bkg[i] << "\n"<< std::endl;

    MC_x_sec[i] =  (N_sel[i] - N_bkg[i]) / ( Efficiency[i] * flux_mc * targets_mc);
    if (DEBUG) std::cout << "New MC X-section [10^-39 cm^2]\t\t" << MC_x_sec[i]/1e-39 << std::endl;

    MC_weighted_xsec_file << MC_x_sec[i] << "\n"; // Put the new MC cross sections into a text file. 
    // std::cout << "MC_x_sec[i]:\t" << MC_x_sec[i] << std::endl;

    double num_bkg_data = N_bkg[i] * mc_scale_factor + intime_cosmics_bkg * intime_cosmic_scale_factor; // scale the number of background events

    std::cout << "Num Bkg Data:\t" << num_bkg_data << "\tintime_cosmics_bkg * intime_cosmic_scale_factor\t" << intime_cosmics_bkg * intime_cosmic_scale_factor << std::endl;
    
    Data_x_sec[i] = (num_selected_data - num_bkg_data) / ( Efficiency[i] * flux_data * targets_data  );

    Data_weighted_xsec_file << Data_x_sec[i] << "\n"; // Put the new Data cross sections into a text file. 
    if (DEBUG) std::cout << "New Data X-section [10^-39 cm^2]\t" << Data_x_sec[i]/1e-39 << std::endl;

    // Now fill the model vs cross section histogram and ttree.
    gModelXsec->SetPoint(i, i, MC_x_sec[i]);
    
  }
  DataTree->Fill();

}

DEFINE_ART_MODULE(EventWeightReader)
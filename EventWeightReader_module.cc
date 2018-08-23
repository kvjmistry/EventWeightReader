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

private:

  // Declare member data here.
  TTree* DataTree;
  TTree* KeepTree;
  std::map<std::string, std::vector<double> > weights;
  int run, subrun, evt;

  TH1D *hInterestingWeights;
  TH1D *hDiscardWeights;

  std::vector<std::string> KeepProcess;
  std::vector<std::string> DiscardProcess;


};


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

  hInterestingWeights = 		tfs->make<TH1D>("hInterestingWeights", "Weights to Keep",	20,  0., 10.);
  hDiscardWeights     = 		tfs->make<TH1D>("hDiscardWeights",     "Weights to Discard",		20  ,0., 10.);

  // Create the TTree and add relavent branches
	DataTree = tfs->make<TTree>("EventTree","EventTree"); 
  KeepTree = tfs->make<TTree>("KeepTree","KeepTree"); 

  // Event Information
	DataTree->Branch("run", &run);
	DataTree->Branch("subrun",&subrun);
  DataTree->Branch("event",&evt);
  DataTree->Branch("weights",&weights);
  KeepTree->Branch("keep",&KeepProcess);
  KeepTree->Branch("discard",&DiscardProcess);

}
void EventWeightReader::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  // Determine event ID
  run = e.id().run();
  subrun = e.id().subRun();
  evt = e.id().event();
  
  TString GenieNames;

  auto GenieEW_Handle = e.getValidHandle<std::vector<evwgh::MCEventWeight>>("mcweight"); // Request the mcweight data product

  std::vector<evwgh::MCEventWeight> const& GenieEWvec(*GenieEW_Handle);  

  for (evwgh::MCEventWeight const& GenieEW: GenieEWvec) {
    // Grab the weights 
    weights = GenieEW.fWeight;
  }
  
  // Loop over te weights and print their name and values. 
  for (auto const& it : weights) {
    GenieNames = it.first; 
    std::cout << "\n" << GenieNames << std::endl;
    
    for (unsigned int i = 0; i < it.second.size(); i++){
      std::cout << it.second[i] << "\t";

      // Fill histograms whcih have non 1 and 1 values 
      if (it.second[i] == 1){
        hDiscardWeights->Fill(GenieNames,1); // Reject
        
        // look to see if already found string in vector
        if (std::find(DiscardProcess.begin(), DiscardProcess.end(), it.first) != DiscardProcess.end()){}
        else {
          DiscardProcess.push_back(it.first);
        }
      }

      else {
        hInterestingWeights->Fill(GenieNames,1); // Pass
        
        // look to see if already found string in vector
        if (std::find(KeepProcess.begin(), KeepProcess.end(), it.first) != KeepProcess.end()){}
        else {
          KeepProcess.push_back(it.first);
          }
         
      
      }// end else 

    } // end loop over weights. 
    std::cout << std::endl;
    
  }

  



DataTree->Fill();

}



void EventWeightReader::endJob()
{
  // Implementation of optional member function here.
  //KeepTree->Fill();

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
}

DEFINE_ART_MODULE(EventWeightReader)

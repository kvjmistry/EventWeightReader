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
#include "TFile.h"

#include <iostream>
#include <fstream>

#include "Event_List.h"

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
	void AddWeights(std::vector<Event_List> &N, art::Event const & e);
	void ReadEvents(const char *filename, std::vector<int> &N_evt );

private:

	// Declare member data here.
	int run, subrun, evt;
	std::map<std::string, std::vector<double> > weights; // Map (Model Name, Weight Vector in each universe)

	int Iterations{1};          // Number of times looped event to see how quickly module is running

	double total_in{0}; 		// Total number of matched events 
	double tot_gen{0};  		// counter for the total number of gen events read in
	double tot_sel{0};  		// counter for the total number of sel events read in
	double tot_sig{0};  		// counter for the total number of sig events read in
	double tot_bkg{0};  		// counter for the total number of bkg events read in

	// Vectors for cross section calculation. 
	std::vector<int> N_gen_evt, N_sig_evt, N_bkg_evt, N_sel_evt, N_filt;          // Vectors of event numbers from infile for generated, signal, background, selected, filtinlist
	std::vector<Event_List> N_gen, N_sig, N_bkg, Data_x_sec, Efficiency;   // Event lists of events with new weights for generated, signal, background, selected

	Event_List temp_gen, temp_sig, temp_bkg, temp_xsec, temp_eff;       // Temporary objects for filling the ttree with

	// Flux
	const double flux_mc{4.19844e+10};
	const double flux_data{5.45797e+9};

	// Num Targets
	const double targets_mc{3.50191e+31};
	const double targets_data{3.4723e+31};
	
	// DATA
	const double intime_cosmics_bkg{81};              // Number of intime cosmics for background
	const double num_selected_data{203};              // The number of selected events in data
	const double intime_cosmic_scale_factor{0.97918}; // Scale factor to apply to the intime cosimic background
	const double mc_scale_factor{0.13};               // Scale factor to apply to the mc background

	// DEBUG
	bool DEBUG{true};

	// TTree
	TTree *DataTree;
	// TTree *TestTree;

	// labels
	std::vector<std::string> labels_genie {
		"genie_qema",
		"genie_ncelAxial",       "genie_ncelEta",
		"genie_ccresAxial",      "genie_ccresVector",
		"genie_ncresAxial",      "genie_ncresVector",
		"genie_cohMA",           "genie_cohR0",
		"genie_NonResRvp1pi",    "genie_NonResRvbarp1pi",
		"genie_NonResRvp2pi",    "genie_NonResRvbarp2pi",
		"genie_ResDecayGamma" ,  "genie_ResDecayTheta",
		"genie_NC",
		"genie_DISAth",          "genie_DISBth",          "genie_DISCv1u",       "genie_DISCv2u",
		// "genie_AGKYxF",       "genie_AGKYpT",          // removed due to just giving back weights of 1
		"genie_FormZone",
		"genie_FermiGasModelKf",
		"genie_IntraNukeNmfp",   "genie_IntraNukeNcex",   "genie_IntraNukeNel",
		"genie_IntraNukeNinel",  "genie_IntraNukeNabs",   "genie_IntraNukeNpi",
		"genie_IntraNukePImfp",  "genie_IntraNukePIcex",  "genie_IntraNukePIel",
		"genie_IntraNukePIinel", "genie_IntraNukePIabs",  "genie_IntraNukePIpi"
	};
	
	std::vector<std::string> labels_model { "model_q0q3_ccmec", "model_q0q3_ccqe" };
	std::vector<std::string> labels_reinteractions { "reinteractions_proton", "reinteractions_piplus", "reinteractions_piminus" };
	

};

// A function that loops over all the parameter weights and universes and re-weights the desired events. 
void EventWeightReader::AddWeights(std::vector<Event_List> &N, art::Event const & e){

	auto GenieEW_Handle = e.getValidHandle<std::vector<evwgh::MCEventWeight>>("mcweight"); // Request the mcweight data product

	if (GenieEW_Handle.isValid())  std::cout << "[Analyze] GenieEW_Handle is valid" << std::endl; 
	
	std::vector<evwgh::MCEventWeight> const& GenieEWvec(*GenieEW_Handle);

	//std::cout << "Size of genueEWvec: " <<  GenieEWvec.size() << "  If this size is not equal to 1 then we have a problem!!!"<< std::endl;

	for (evwgh::MCEventWeight const& GenieEW: GenieEWvec) weights = GenieEW.fWeight; // Grab the weights 

	// Loop over the labels
	for (unsigned int j=0; j < N.size(); j++){
		
		// Loop over the labels in the event weight object
		for (auto const& it : weights) {
				
			if (it.first.find(N.at(j).label.c_str()) != std::string::npos) { // match up the correct parameter

				// Insert some debug info to see if what we are doing here is correct
				// std::cout << it.first << "  " << N.at(j).label << std::endl;

				if (N.at(j).N_reweight.size() == 0) N.at(j).N_reweight.resize(it.second.size()); // resize to the number of universes if we havent already done so
			
				// Loop over each universe for a parameter and weight the event
				for (unsigned int i = 0; i < it.second.size(); i++) N.at(j).N_reweight.at(i) += it.second.at(i) ; // Weight the event
			}
		} // END loop over parameters
	}
}

// Function that reads in the event numbers from a text file and adds those event numbers to a vector. 
void EventWeightReader::ReadEvents(const char *filename, std::vector<int> &N_evt ){
	
	std::ifstream fileIN; 

	fileIN.open(filename); // Open the file
	
	if (!fileIN.good()) { // Check if the file opened correctly
		std::cerr << "Error: file:\t" << filename <<"\tcould not be opened" << std::endl;
		exit(1);
	}

	double temp, temp_run, temp_sr; // Use a temp var to get the values and push back

	if (fileIN.is_open()) { 
		
		while ( fileIN >> temp_run >> temp_sr >> temp) { // loop over lines in file
			
			N_evt.push_back(temp);
		}
		
		fileIN.close();
	}

}

// Same as Read events function, but now split the selected events into their new catagories. 
void ReadEventList(const char *filename, std::vector<int> &N_sig_evt, std::vector<int> &N_bkg_evt, std::vector<int> &N_sel_evt ){

	// event number, classifier type, mc_nu_id
	std::vector<int>         N_evtnum;
	std::vector<std::string> class_type;
	std::vector<int>         mc_nu_id;

	std::ifstream fileIN;

	fileIN.open(filename); // Open the file

	if (!fileIN.good()) {  // Check if the file opened correctly
			std::cerr << "Error: file:\t" << filename <<"\tcould not be opened" << std::endl;
			exit(1);
	}

	int temp_evtnum,  temp_mc_nu_id;
	std::string temp_class_type;

	if (fileIN.is_open()) {

		// loop over lines in file
		while ( fileIN >> temp_evtnum >> temp_class_type >> temp_mc_nu_id) {

			N_evtnum.push_back(temp_evtnum);
			class_type.push_back(temp_class_type);
			mc_nu_id.push_back(temp_mc_nu_id);
		}

		fileIN.close();
	}

	// Now got the info, recatagorise and get the relavent events 
	for (unsigned int i = 0; i < N_evtnum.size(); i++){
		// Push back selected events
		N_sel_evt.push_back(N_evtnum[i]);
		
		// Signal
		if ((class_type[i].compare(0,6,"nue_cc") == 0 && class_type[i] != "nue_cc_out_fv" && class_type[i] != "nue_cc_mixed") || (class_type[i].compare(0,10,"nue_bar_cc") == 0  && class_type[i] != "nue_bar_mixed") ) {
			N_sig_evt.push_back(N_evtnum[i]);
		}
		// Dirt
		else if (class_type[i] == "Dirt"){
		}
		// Background
		else {
			N_bkg_evt.push_back(N_evtnum[i]);
		}
	}

}

EventWeightReader::EventWeightReader(fhicl::ParameterSet const & p) : EDAnalyzer(p) {}

void EventWeightReader::beginJob() {
	// Implementation of optional member function here.
	// Access ART's TFileService, which will handle histograms/trees/etc.
	art::ServiceHandle<art::TFileService> tfs;

	// Create the TTree and add relavent branches
	DataTree = tfs->make<TTree>("XSectionTree","XSectionTree");
	DataTree->Branch("Sig",         &temp_sig);
	DataTree->Branch("Bkg",         &temp_bkg);
	DataTree->Branch("Gen",         &temp_gen);
	DataTree->Branch("eff",         &temp_eff);
	DataTree->Branch("xsec",        &temp_xsec);
	// Load in the file containing the event number for Signal_Generated, N_gen or Signal_Selected N_sig
	if (DEBUG) std::cout << "\nNow reading in event files!" << std::endl;

	// ++++++++++++++++++++ N_gen +++++++++++++++++++++++++++++
	
	ReadEvents("generated_events.txt", N_gen_evt ); 

	ReadEventList("selected_events.txt", N_sig_evt, N_bkg_evt, N_sel_evt ); 

	if (DEBUG) std::cout << "Size of Generated vector:  \t"<< N_gen_evt.size() << std::endl;
	if (DEBUG) std::cout << "Size of signal vector:     \t"<< N_sig_evt.size() << std::endl;
	if (DEBUG) std::cout << "Size of selected vector:   \t"<< N_sel_evt.size() << std::endl;
	if (DEBUG) std::cout << "Size of background vector: \t"<< N_bkg_evt.size() << std::endl;

	if (DEBUG) std::cout << "\nFinished reading in event files!\n" << std::endl;

	std::cout << "++++++++++++++++++++++++++++++++++" << std::endl;
	std::cout << "Now creating the event list vector..."<<  std::endl;

	for (unsigned int i = 0; i < labels_genie.size(); i++) {

		N_gen.push_back(Event_List(labels_genie[i],"genie"));
		N_sig.push_back(Event_List(labels_genie[i],"genie"));
		N_bkg.push_back(Event_List(labels_genie[i],"genie"));
		Data_x_sec.push_back(Event_List(labels_genie[i],"genie"));
		Efficiency.push_back(Event_List(labels_genie[i],"genie"));
	}

	for (unsigned int i = 0; i < labels_model.size(); i++) {

		N_gen.push_back(Event_List(labels_model[i],"model"));
		N_sig.push_back(Event_List(labels_model[i],"model"));
		N_bkg.push_back(Event_List(labels_model[i],"model"));
		Data_x_sec.push_back(Event_List(labels_model[i],"model"));
		Efficiency.push_back(Event_List(labels_model[i],"model"));
	}

	for (unsigned int i = 0; i < labels_reinteractions.size(); i++) {

		N_gen.push_back(Event_List(labels_reinteractions[i],"reinteractions"));
		N_sig.push_back(Event_List(labels_reinteractions[i],"reinteractions"));
		N_bkg.push_back(Event_List(labels_reinteractions[i],"reinteractions"));
		Data_x_sec.push_back(Event_List(labels_reinteractions[i],"reinteractions"));
		Efficiency.push_back(Event_List(labels_reinteractions[i],"reinteractions"));
	}
	std::cout << "Done creating the event list vector!"<<  std::endl;
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
	std::cout << "Event processed:\t" << Iterations << std::endl;
	std::cout << "++++++++++++++++++++++++++++++++++" << std::endl;
	Iterations++;

	// Choose wheather to populate which varible depending on the event number

	// ++++++++++++++++++++ N_gen +++++++++++++++++++++++++++++
	if ( std::find( N_gen_evt.begin(), N_gen_evt.end(), evt) != N_gen_evt.end()){
		
		if (DEBUG) std::cout << "Matched a Generated Event\t" << std::endl;
		tot_gen++;

		// Found event in the N_gen vector and so add weights for this event
		AddWeights(N_gen, e);
	}
	// ++++++++++++++++++++ N_sel +++++++++++++++++++++++++++++
	if ( std::find( N_sel_evt.begin(), N_sel_evt.end(), evt) != N_sel_evt.end()){
		
		if (DEBUG) std::cout << "Matched a Selected Event\t" << std::endl;

		tot_sel++;

		// Found event in the N_sel vector and so add weights for this event
		// AddWeights(N_sel, e); // for now we only care about data x sec so forget about this
	}
	// ++++++++++++++++++++ N_sig +++++++++++++++++++++++++++++
	if ( std::find( N_sig_evt.begin(), N_sig_evt.end(), evt) != N_sig_evt.end()){
		
		if (DEBUG) std::cout << "Matched a Signal Selected Event\t" << std::endl;

		tot_sig++;

		// Found event in the N_sig vector and so add weights for this event
		AddWeights(N_sig, e);
	}
	// ++++++++++++++++++++ N_bkg +++++++++++++++++++++++++++++
	if ( std::find( N_bkg_evt.begin(), N_bkg_evt.end(), evt) != N_bkg_evt.end()){
		
		if (DEBUG) std::cout << "Matched a Background Selected Event\t" << std::endl;

		tot_bkg++;

		// Found event in the N_bkg vector and so add weights for this event
		AddWeights(N_bkg, e);
	}

	total_in = tot_gen + tot_bkg;
	
}

void EventWeightReader::endJob() {

	// Implementation of optional member function here.
	std::cout << "\nBeginning END JOB..." << std::endl;
	
	std::cout << "\n=======================" << std::endl;
	std::cout << "Read in Events\n" << std::endl;
	std::cout << "=======================\n" << std::endl;
	std::cout << "\nN_gen:\t" << tot_gen << std::endl;
	std::cout << "\nN_sel:\t" << tot_sel << std::endl;
	std::cout << "\nN_sig:\t" << tot_sig << std::endl;
	std::cout << "\nN_bkg:\t" << tot_bkg << std::endl;
	std::cout << "\nTotal:\t" << total_in << std::endl;

	// std::cout << "\n=======================" << std::endl;
	// std::cout << "CV Cross sections\n" << std::endl;
	// std::cout << "=======================\n" << std::endl;
	// std::cout << "\nMC:\t" << (tot_sel + 2 - tot_bkg) / ( (tot_sig + 2 ) / (tot_gen + 4) * flux_mc * targets_mc) << std::endl;
	// std::cout << "\nData:\t" << ( num_selected_data - ((tot_bkg) * mc_scale_factor + intime_cosmics_bkg * intime_cosmic_scale_factor )) / ( ((tot_sig + 2)/(tot_gen + 4 )) * flux_data * targets_data) << std::endl;

	// std::cout << "\n=======================" << std::endl;
	// std::cout << "Now Re-Calculating cross sections" << std::endl;
	// std::cout << "=======================\n" << std::endl;

	// // Open a file with the new x section values in
	// std::ofstream MC_weighted_xsec_file;
	// MC_weighted_xsec_file.open("MC_weighted_xsec_file.txt");

	// std::ofstream Data_weighted_xsec_file;
	// Data_weighted_xsec_file.open("Data_weighted_xsec_file.txt");

	// MC_x_sec.resize(N_gen.size()); // Resize
	// Data_x_sec.resize(N_gen.size());
	// Efficiency.resize(N_gen.size());

	// // Calculate the new Cross section. 
	// for (unsigned int i{0}; i < N_gen.size(); i++){

	// 	// Recalculations due to not weighting non MC genie stuff and other bugs
	// 	N_gen[i] = N_gen[i] + 4  ;  // The plus 4 is for the tpc obj of size zero bug
	// 	N_sel[i] = tot_sel  + 2; ;  // ANDY F: do not reweight the selected events in MC so it is like data
	// 	N_sig[i] = N_sig[i] + 2;    // Missing two events from somewhere
	// 	N_bkg[i] = N_bkg[i];

	// 	Efficiency[i] = N_sig[i] / N_gen[i];  // 0.0884133 CV efficiency
	// 	if (DEBUG) std::cout << "\n+++++++\nEfficiency\t" << Efficiency[i] << std::endl;

	// 	std::cout << "\nN_gen:\t" << N_gen[i] << std::endl;
	// 	std::cout << "N_sel:\t" << N_sel[i] << std::endl;
	// 	std::cout << "N_sig:\t" << N_sig[i] << std::endl;
	// 	std::cout << "N_bkg:\t" << N_bkg[i] << "\n"<< std::endl;

	// 	MC_x_sec[i] =  (N_sel[i] - N_bkg[i]) / ( Efficiency[i] * flux_mc * targets_mc);
	// 	if (DEBUG) std::cout << "New MC X-section [10^-39 cm^2]\t\t" << MC_x_sec[i]/1e-39 << std::endl;

	// 	MC_weighted_xsec_file << MC_x_sec[i] << "\n"; // Put the new MC cross sections into a text file. 
	// 	// std::cout << "MC_x_sec[i]:\t" << MC_x_sec[i] << std::endl;

	// 	double num_bkg_data = N_bkg[i] * mc_scale_factor + intime_cosmics_bkg * intime_cosmic_scale_factor; // scale the number of background events

	// 	std::cout << "Num Bkg Data:\t" << num_bkg_data << "\tintime_cosmics_bkg * intime_cosmic_scale_factor\t" << intime_cosmics_bkg * intime_cosmic_scale_factor << std::endl;
		
	// 	Data_x_sec[i] = (num_selected_data - num_bkg_data) / ( Efficiency[i] * flux_data * targets_data  );

	// 	Data_weighted_xsec_file << Data_x_sec[i] << "\n"; // Put the new Data cross sections into a text file. 
	// 	if (DEBUG) std::cout << "New Data X-section [10^-39 cm^2]\t" << Data_x_sec[i]/1e-39 << std::endl;
		
	// }

	for (unsigned int i=0; i < N_gen.size(); i++){
		temp_gen  = N_gen[i];
		temp_sig  = N_sig[i];
		temp_bkg  = N_bkg[i];
		temp_eff  = Efficiency[i];
		temp_xsec = Data_x_sec[i];
		DataTree->Fill();
	}
	
}

DEFINE_ART_MODULE(EventWeightReader)

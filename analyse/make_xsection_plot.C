// This program will read in the calculated cross section values after running the eventweightreader module.
// It will then do some plotting.

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TGraphErrors.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TExec.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLine.h"

#include <iostream>
#include <fstream>
#include <string>
#include "../Event_List.h"
#include "functions.h"
#include <iomanip>

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void make_xsection_plot(){

	TFile* f;
	TTree* xsectree;

	// First read in the ttree and loop over the entries
	bool boolfile  = GetFile(f , "../NuMIEventWeight.root"); if (boolfile == false) gSystem->Exit(0); // File with the weights
	bool booltree  = GetTree(f, xsectree ,"microboonewvtof/XSectionTree"); // Sig, Bkg, Gen, eff, xsec
	
	Sig  = new Event_List;
	Bkg  = new Event_List;
	Gen  = new Event_List;
	eff  = new Event_List;
	xsec = new Event_List;

	xsectree -> SetBranchAddress("Sig",  &Sig);
	xsectree -> SetBranchAddress("Bkg",  &Bkg);
	xsectree -> SetBranchAddress("Gen",  &Gen);
	xsectree -> SetBranchAddress("eff",  &eff);
	xsectree -> SetBranchAddress("xsec", &xsec);

	const int tree_total_entries = xsectree->GetEntries();
	std::cout << "Total Labels: " << tree_total_entries << std::endl;

	// CV -- should remove the hardcoded values
	double CV_Sig  = 642;
	double CV_Bkg  = 356;
	double CV_Gen  = 7103;
	double CV_eff  = 0.0903843;
	double CV_xsec = 4.67e-39;

	double Err_Sig,  Err_Sig_p1sig,  Err_Sig_m1sig;
	double Err_Bkg,  Err_Bkg_p1sig,  Err_Bkg_m1sig;
	double Err_Gen,  Err_Gen_p1sig,  Err_Gen_m1sig;
	double Err_eff,  Err_eff_p1sig,  Err_eff_m1sig;
	double Err_xsec, Err_xsec_p1sig, Err_xsec_m1sig;

	int n_uni_genie, n_uni_model, n_uni_reinteractions;

	int num_genie = 0;

	// ----------------------
	//		Event loop
	// ----------------------
	std::cout << "Starting loop over the reweighters..." << std::endl;
	for (int label = 0; label < tree_total_entries; label++){
			xsectree->GetEntry(label);

			// Erase the genie and model part of the label to tidy up
			if ( Gen -> reweighter == "genie" || Gen -> reweighter == "model") {
				Sig  -> label.erase(0,6);
				Bkg  -> label.erase(0,6);
				Gen  -> label.erase(0,6);
				eff  -> label.erase(0,6);
				xsec -> label.erase(0,6);
			}
			
			std::cout << "-----------------------------"<< std::endl;
			std::cout << "Reweighter: " << Gen -> reweighter        << std::endl;
			std::cout << "Label: "      << Gen -> label             << std::endl;
			std::cout << "Mode: "       << Gen -> mode              << std::endl;
			std::cout << "Num of Uni: " << Gen -> N_reweight.size() << std::endl;

			if ( Gen -> N_reweight.size() == 0) continue; // check to see if will break the plots though 

			// Get the number of universes
			if (Gen -> reweighter == "genie" && Gen -> mode == "multisim") n_uni_genie =  Gen -> N_reweight.size();
			if (Gen -> reweighter == "model" && Gen -> mode == "multisim") n_uni_model =  Gen -> N_reweight.size();
			if (Gen -> reweighter == "reinteractions" && Gen -> mode == "multisim") n_uni_reinteractions =  Gen -> N_reweight.size();
			
			// multisim
			if (Gen -> mode == "multisim"){
				// Get the standard deviation of the universe 
				Err_Sig  = STD(Sig  -> N_reweight,  CV_Sig);
				Err_Bkg  = STD(Bkg  -> N_reweight,  CV_Bkg);
				Err_Gen  = STD(Gen  -> N_reweight,  CV_Gen);
				Err_eff  = STD(eff  -> N_reweight,  CV_eff);
				Err_xsec = STD(xsec -> N_reweight, CV_xsec);
				
				Print_Error("Sig" , CV_Sig,  "     Err", Err_Sig);
				Print_Error("Bkg" , CV_Bkg,  "     Err", Err_Bkg);
				Print_Error("Gen" , CV_Gen,  "    Err", Err_Gen);
				Print_Error("Eff" , CV_eff,  "  Err", Err_eff);
				Print_Error("xsec", CV_xsec, "Err", Err_xsec);
			}
			// Unisim
			else {
				// Get the +/- 1 sigma uncertainties
				Unisim_err(Sig  -> N_reweight,  CV_Sig,  Err_Sig_p1sig,  Err_Sig_m1sig);
				Unisim_err(Bkg  -> N_reweight,  CV_Bkg,  Err_Bkg_p1sig,  Err_Bkg_m1sig);
				Unisim_err(Gen  -> N_reweight,  CV_Gen,  Err_Gen_p1sig,  Err_Gen_m1sig);
				Unisim_err(eff  -> N_reweight,  CV_eff,  Err_eff_p1sig,  Err_eff_m1sig);
				Unisim_err(xsec -> N_reweight,  CV_xsec, Err_xsec_p1sig, Err_xsec_m1sig);
				
				Print_Error("Sig" , CV_Sig,  Err_Sig_p1sig,  Err_Sig_m1sig);
				Print_Error("Bkg" , CV_Bkg,  Err_Bkg_p1sig,  Err_Bkg_m1sig);
				Print_Error("Gen" , CV_Gen,  Err_Gen_p1sig,  Err_Gen_m1sig);
				Print_Error("Eff" , CV_eff,  Err_eff_p1sig,  Err_eff_m1sig);
				Print_Error("xsec", CV_xsec, Err_xsec_p1sig, Err_xsec_m1sig);
			}

			std::cout << "-----------------------------\n"<< std::endl;

			// GENIE
			if ( Gen -> reweighter == "genie"){
				num_genie+=2; // Use this to get the correct bin index in the interaction params to set the bin error
				
				// Multisim
				if (Gen -> mode == "multisim"){	
					hCV_genie_Sig->Fill(Sig -> label.c_str(), 0);
					hCV_genie_Sig->SetBinError(label, 100 * Err_Sig / CV_Sig );

					hCV_genie_Bkg->Fill(Bkg -> label.c_str(), 0);
					hCV_genie_Bkg->SetBinError(label, 100 * Err_Bkg / CV_Bkg );

					hCV_genie_Gen->Fill(Gen -> label.c_str(), 0);
					hCV_genie_Gen->SetBinError(label, 100 * Err_Gen / CV_Gen );

					hCV_genie_eff->Fill(eff -> label.c_str(), 0);
					hCV_genie_eff->SetBinError(label, 100 * Err_eff / CV_eff );

					hCV_genie_xsec->Fill(xsec -> label.c_str(), 0);
					hCV_genie_xsec->SetBinError(label, 100 * Err_xsec / CV_xsec );
				}
				// Unisim
				else {
					hCV_genie_Sig_p1sig->Fill(Sig -> label.c_str(), 100 * Err_Sig_p1sig / CV_Sig);
					hCV_genie_Sig_m1sig->Fill(Sig -> label.c_str(), 100 * Err_Sig_m1sig / CV_Sig);

					hCV_genie_Bkg_p1sig->Fill(Bkg -> label.c_str(), 100 * Err_Bkg_p1sig / CV_Bkg);
					hCV_genie_Bkg_m1sig->Fill(Bkg -> label.c_str(), 100 * Err_Bkg_m1sig / CV_Bkg);

					hCV_genie_Gen_p1sig->Fill(Gen -> label.c_str(), 100 * Err_Gen_p1sig / CV_Gen);
					hCV_genie_Gen_m1sig->Fill(Gen -> label.c_str(), 100 * Err_Gen_m1sig / CV_Gen);

					hCV_genie_eff_p1sig->Fill(eff -> label.c_str(), 100 * Err_eff_p1sig / CV_eff);
					hCV_genie_eff_m1sig->Fill(eff -> label.c_str(), 100 * Err_eff_m1sig / CV_eff);

					hCV_genie_xsec_p1sig->Fill(xsec -> label.c_str(), 100 * Err_xsec_p1sig / CV_xsec);
					hCV_genie_xsec_m1sig->Fill(xsec -> label.c_str(), 100 * Err_xsec_m1sig / CV_xsec);

				}

			}
			// Interaction
			else {
				// multisim
				hCV_interaction_Sig->Fill(Sig -> label.c_str(), 0);
				hCV_interaction_Sig->SetBinError(label - num_genie + 1, 100 * Err_Sig / CV_Sig );

				hCV_interaction_Bkg->Fill(Bkg -> label.c_str(), 0);
				hCV_interaction_Bkg->SetBinError(label - num_genie + 1, 100 * Err_Bkg / CV_Bkg );

				hCV_interaction_Gen->Fill(Gen -> label.c_str(), 0);
				hCV_interaction_Gen->SetBinError(label - num_genie + 1, 100 * Err_Gen / CV_Gen );

				hCV_interaction_eff->Fill(eff -> label.c_str(), 0);
				hCV_interaction_eff->SetBinError(label - num_genie + 1, 100 * Err_eff / CV_eff );

				hCV_interaction_xsec->Fill(xsec -> label.c_str(), 0);
				hCV_interaction_xsec->SetBinError(label - num_genie + 1, 100 * Err_xsec / CV_xsec );

			}

	} // End loop over labels

	gStyle->SetOptStat(0);
	gSystem->Exec("if [ ! -d \"plots\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir plots; fi");

	// Plot the genie histograms
	hist_options(c_genie_Sig,  hCV_genie_Sig,  hCV_genie_Sig_p1sig,  hCV_genie_Sig_m1sig,  15, n_uni_genie,  "plots/genie_sig.pdf"  );  // genie_sig
	hist_options(c_genie_Bkg,  hCV_genie_Bkg,  hCV_genie_Bkg_p1sig,  hCV_genie_Bkg_m1sig,  15, n_uni_genie,  "plots/genie_bkg.pdf"  );  // genie_bkg
	hist_options(c_genie_Gen,  hCV_genie_Gen,  hCV_genie_Gen_p1sig,  hCV_genie_Gen_m1sig,  15, n_uni_genie,  "plots/genie_gen.pdf"  );  // genie_gen
	hist_options(c_genie_eff,  hCV_genie_eff,  hCV_genie_eff_p1sig,  hCV_genie_eff_m1sig,  15, n_uni_genie,  "plots/genie_eff.pdf"  );  // genie_eff
	hist_options(c_genie_xsec, hCV_genie_xsec, hCV_genie_xsec_p1sig, hCV_genie_xsec_m1sig, 15, n_uni_genie, "plots/genie_xsec.pdf" );  // genie_xsec

	// Interactions
	hist_options(c_interaction_Sig,  hCV_interaction_Sig,  5, n_uni_model, n_uni_reinteractions, "plots/interaction_sig.pdf" ); // Interaction_sig
	hist_options(c_interaction_Bkg,  hCV_interaction_Bkg,  5, n_uni_model, n_uni_reinteractions, "plots/interaction_bkg.pdf" ); // Interaction_bkg
	hist_options(c_interaction_Gen,  hCV_interaction_Gen,  5, n_uni_model, n_uni_reinteractions, "plots/interaction_gen.pdf" ); // Interaction_gen
	hist_options(c_interaction_eff,  hCV_interaction_eff,  5, n_uni_model, n_uni_reinteractions, "plots/interaction_eff.pdf" ); // Interaction_eff
	hist_options(c_interaction_xsec, hCV_interaction_xsec, 5, n_uni_model, n_uni_reinteractions, "plots/interaction_xsec.pdf"); // Interaction_xsec
	
}
// END MAIN
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


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

	TH1D *hCV_genie_Sig  = new TH1D("hCV_genie_Sig", "Genie Sig;; Percentage Error [%]",   34, 0, 34);
	TH1D *hCV_genie_Bkg  = new TH1D("hCV_genie_Bkg", "Genie Bkg;; Percentage Error [%]",   34, 0, 34);
	TH1D *hCV_genie_Gen  = new TH1D("hCV_genie_Gen", "Genie Gen;; Percentage Error [%]",   34, 0, 34);
	TH1D *hCV_genie_eff  = new TH1D("hCV_genie_eff", "Genie Eff;; Percentage Error [%]",   34, 0, 34);
	TH1D *hCV_genie_xsec = new TH1D("hCV_genie_xsec","Genie X Sec;; Percentage Error [%]", 34, 0, 34);

	TCanvas * c_genie_Sig  = new TCanvas();
	TCanvas * c_genie_Bkg  = new TCanvas();
	TCanvas * c_genie_Gen  = new TCanvas();
	TCanvas * c_genie_eff  = new TCanvas();
	TCanvas * c_genie_xsec = new TCanvas();
	
	TH1D *hCV_other_Sig  = new TH1D("hCV_other_Sig","Other Sig;; Percentage Error [%]",    5, 0, 5);
	TH1D *hCV_other_Bkg  = new TH1D("hCV_other_Bkg","Other Bkg;; Percentage Error [%]",    5, 0, 5);
	TH1D *hCV_other_Gen  = new TH1D("hCV_other_Gen","Other Gen;; Percentage Error [%]",    5, 0, 5);
	TH1D *hCV_other_eff  = new TH1D("hCV_other_eff","Other Eff;; Percentage Error [%]",    5, 0, 5);
	TH1D *hCV_other_xsec = new TH1D("hCV_other_xsec","Other X Sec;; Percentage Error [%]", 5, 0, 5);

	TCanvas * c_other_Sig  = new TCanvas();
	TCanvas * c_other_Bkg  = new TCanvas();
	TCanvas * c_other_Gen  = new TCanvas();
	TCanvas * c_other_eff  = new TCanvas();
	TCanvas * c_other_xsec = new TCanvas();

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

	// CV -- will want to remove the hardcoded values
	double CV_Sig  = 622;
	double CV_Bkg  = 427;
	double CV_Gen  = 7130;
	double CV_eff  = 0.0872;
	double CV_xsec = 4.5e-39;


	int num_genie = 0;

	// ----------------------
	//		Event loop
	// ----------------------
	std::cout << "Starting loop over the labels..." << std::endl;
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
			std::cout << "Num of Uni: " << Gen -> N_reweight.size() << std::endl;
			

			// Get the standard deviation of the universe 
			double Err_Sig  = STD(Sig  -> N_reweight,  CV_Sig);
			double Err_Bkg  = STD(Bkg  -> N_reweight,  CV_Bkg);
			double Err_Gen  = STD(Gen  -> N_reweight,  CV_Gen);
			double Err_eff  = STD(eff  -> N_reweight,  CV_eff);
			double Err_xsec = STD(xsec -> N_reweight, CV_xsec);

			Print_Error("Sig" , CV_Sig,  "     Err", Err_Sig);
			Print_Error("Bkg" , CV_Bkg,  "     Err", Err_Bkg);
			Print_Error("Gen" , CV_Gen,  "    Err", Err_Gen);
			Print_Error("Eff" , CV_eff,  "  Err", Err_eff);
			Print_Error("xsec", CV_xsec, "Err", Err_xsec);
			std::cout << "-----------------------------\n"<< std::endl;

			if ( Gen -> reweighter == "genie"){
				num_genie++;

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
			else {
				hCV_other_Sig->Fill(Sig -> label.c_str(), 0);
				hCV_other_Sig->SetBinError(label - num_genie + 1, 100 * Err_Sig / CV_Sig );

				hCV_other_Bkg->Fill(Bkg -> label.c_str(), 0);
				hCV_other_Bkg->SetBinError(label - num_genie + 1, 100 * Err_Bkg / CV_Bkg );

				hCV_other_Gen->Fill(Gen -> label.c_str(), 0);
				hCV_other_Gen->SetBinError(label - num_genie + 1, 100 * Err_Gen / CV_Gen );

				hCV_other_eff->Fill(eff -> label.c_str(), 0);
				hCV_other_eff->SetBinError(label - num_genie + 1, 100 * Err_eff / CV_eff );

				hCV_other_xsec->Fill(xsec -> label.c_str(), 0);
				hCV_other_xsec->SetBinError(label - num_genie + 1, 100 * Err_xsec / CV_xsec );

			}

	} // End loop over labels

	gStyle->SetOptStat(0);

	c_genie_Sig  ->cd();
	gPad->SetBottomMargin(0.2);
	hCV_genie_Sig->LabelsOption("v");
	hCV_genie_Sig->Draw("PE1");

	c_genie_Bkg  ->cd();
	gPad->SetBottomMargin(0.2);
	hCV_genie_Bkg->LabelsOption("v");
	hCV_genie_Bkg->Draw("PE1");

	c_genie_Gen  ->cd();
	gPad->SetBottomMargin(0.2);
	hCV_genie_Gen->LabelsOption("v");
	hCV_genie_Gen->Draw("PE1");

	c_genie_eff  ->cd();
	gPad->SetBottomMargin(0.2);
	hCV_genie_eff->LabelsOption("v");
	hCV_genie_eff->Draw("PE1");

	c_genie_xsec ->cd();
	gPad->SetBottomMargin(0.2);
	hCV_genie_xsec->LabelsOption("v");
	hCV_genie_xsec->Draw("PE1");


	c_other_Sig  ->cd();
	gPad->SetBottomMargin(0.225);
	hCV_other_Sig->LabelsOption("v");
	hCV_other_Sig->Draw("PE1");

	c_other_Bkg  ->cd();
	gPad->SetBottomMargin(0.225);
	hCV_other_Bkg->LabelsOption("v");
	hCV_other_Bkg->Draw("PE1");

	c_other_Gen  ->cd();
	gPad->SetBottomMargin(0.225);
	hCV_other_Gen->LabelsOption("v");
	hCV_other_Gen->Draw("PE1");

	c_other_eff  ->cd();
	gPad->SetBottomMargin(0.225);
	hCV_other_eff->LabelsOption("v");
	hCV_other_eff->Draw("PE1");

	c_other_xsec ->cd();
	gPad->SetBottomMargin(0.225);
	hCV_other_xsec->LabelsOption("v");
	hCV_other_xsec->Draw("PE1");

}
// END MAIN
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
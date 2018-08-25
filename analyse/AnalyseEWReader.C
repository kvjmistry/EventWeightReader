// Program that loads in the ttree created from the NueKinematics Module and makes some plots as does the Larsoft Module. 
// WARNING: make sure a plots folder exits in the current directory otherwise the histogram will not get saved and get an error. 
// Execute over the NueKinematics.root file which has been run over the larsoft module NueKinematics_module.cc
// This file can be run by simply typing root AnalyseNueKinematics.C

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

#include <iostream>
#include <fstream>
#include <string>

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                                       Function Definitions
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                                       Main Function
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AnalyseEWReader() {



    // ++++++++++++++++++++++++++++++++++s++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                  Initialize
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TFile *MyFile = new TFile("plots/Plots_EWReader.root","RECREATE");
    if ( MyFile->IsOpen() ) printf("File opened successfully\n");

    // VAriables
    double mc_xsec_weight{0};
    std::vector<double> WeightList; // List of all weights in one vector. 

    //const double full_error = sqrt( pow(stat_err, 2) + pow(sys_err, 2) );
	//const double full_error = sqrt( pow(0, 2) + pow(sys_err, 2) );
    const double full_error{1.5e-39};
    const double data_xsec{5.34e-39};
    const double mc_xsec{4.83e-39};

    // Trees and Histograms
    TTree *Tree = new TTree("Tree", "Tree");
    TBranch *B_mc_xsec_weight = Tree->Branch("mc_xsec_weight", &mc_xsec_weight, "mc_xsec_weight/D");
    

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                  Load in weights
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	
	std::ifstream weightfile;
	weightfile.open("Genie_weights.txt");
	
    if (!weightfile.good()) {
		// Print error message and exit
		cerr << "Error: weight file could not be opened" << endl;
		gSystem->Exit(1);
	}

    double temp{0};

    if (weightfile.is_open()) {
        while ( !weightfile.eof()) {
            weightfile >> temp;
            
            WeightList.push_back(temp);
    	}
    
	    weightfile.close();
  	}


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                  Loop over weights and reweigh x section MC
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	

	double mc_xec_err;

	for (unsigned int i = 0; i <WeightList.size(); i++ ){ // WeightList.size()
		//std::cout << WeightList[i] << std::endl;
        mc_xsec_weight = mc_xsec * WeightList[i]; 

        //std::cout << mc_xsec_weight/1.0e-39 << "   "<< WeightList[i] <<std::endl;

		//TLine * Tline = new TLine(0.9, mc_xsec_weight, 1.1, mc_xsec_weight);
        //Tline->SetLineColor(kGray);
		//Tline->Draw("same");
        //Tree->Fill();

        mc_xec_err+=  (mc_xsec_weight - mc_xsec)*(mc_xsec_weight - mc_xsec);  

  	}

    // Calculate the standard deviation of the weights for the errorbar
    mc_xec_err = std::sqrt(mc_xec_err/WeightList.size());
    std::cout << "***********************" << std::endl;
    std::cout << "Error on mc xsec\t" << mc_xec_err << std::endl; 
    std::cout << "***********************" << std::endl;

    


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                  Make TGraph for MC x section
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
    TCanvas * c_integrated = new TCanvas();
	c_integrated->cd();
    
    double x_MC[2] = {0.9, 1.1};
	double y_MC[2] = {4.83e-39,4.83e-39};
	double ex_MC[2] = {0.0,0.0};
	double ey_MC[2] = {1.0e-39,1.0e-39} ;

   
    TGraphErrors * MCxsec_graph = new TGraphErrors(2,x_MC ,y_MC , ex_MC, ey_MC);
    MCxsec_graph->SetTitle("Integrated Nue + Nue-bar Xsec");
    MCxsec_graph->GetYaxis()->SetTitle("#sigma [cm^2]");
    MCxsec_graph->GetXaxis()->SetLimits(0.9, 1.1);//5 GeV
	MCxsec_graph->SetMinimum(0); //3.6e-39
    MCxsec_graph->SetMaximum(10e-39); // 7.2e-39
    MCxsec_graph->SetFillColor(kCyan-10);
    MCxsec_graph->Draw("a3");

	TLine * line = new TLine(0.9, mc_xsec, 1.1, mc_xsec);// CV cross section
	line->SetLineColor(kCyan+1);
    line->SetLineWidth(2);
	line->Draw("same");

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                 Create TGraph for x sec DATA
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	

    const int n = 1;
	double x[n] = {1.0};
	double y[n] = {data_xsec};
	double ex[n] = {0.0};
	double ey[n] = {full_error};
	
	TGraphErrors * xsec_graph = new TGraphErrors(n, x, y, ex, ey);
	
	xsec_graph->SetMarkerStyle(3);
    
	xsec_graph->Draw("*, same");
    
	TLegend * leg = new TLegend();
	leg->AddEntry(xsec_graph, "Data Cross Section");
	leg->AddEntry(line, "GENIE");
	leg->Draw("SAME");

	c_integrated->Print("plots/integrated_xsec.eps");


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                     END
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    gSystem->Exit(1); // Quit ROOT

} // END MAIN FUNCTION


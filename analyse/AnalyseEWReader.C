// This program will read in the calculated cross section values after running the eventweightreader module.
// It will then calculate the standard deviation of the cross sections calculated and do some plotting.

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
// Function that reads in the event numbers from a text file and adds those event numbers to a vector. 
void EventWeightReader::ReadEvents(const char *filename, std::vector<int> MC_Xsec ){
	
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
        MC_Xsec.push_back(temp);
    	}
	    
      fileIN.close();
  	}

}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                                       Main Function
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AnalyseEWReader() {
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                  Initialize
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TFile *MyFile = new TFile("plots/Plots_EWReader.root","RECREATE");
    if ( MyFile->IsOpen() ) printf("File opened successfully\n");

    // Variables

    const double full_error{1.5e-39};
    const double data_xsec {5.34e-39}; // Data Cross section
    const double mc_xsec_CV{4.83e-39}; // central value for the x section

    // Trees and Histograms
    // TTree *Tree = new TTree("Tree", "Tree");
    // TBranch *B_mc_xsec= Tree->Branch("mc_xsec_values", &mc_xsec, "mc_xsec_weight/D");
    

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                  Load in the re-weighted cross sections
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	std::vector<double> mc_xsec;

    ReadEvents("MC_weighted_xsec_file.txt", mc_xsec );

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                  Loop over cross sectons and calculate the Standard deviation
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	

	double mc_xsec_err; // systematic error from all the cross sections

	for (unsigned int i = 0; i < mc_xsec.size(); i++ ){ 

		//TLine * Tline = new TLine(0.9, mc_xsec[i], 1.1, mc_xsec[i]);
        //Tline->SetLineColor(kGray);
		//Tline->Draw("same");
        //Tree->Fill();

        mc_xsec_err +=  (mc_xsec - mc_xsec_CV) * (mc_xsec- mc_xsec_CV);  

  	}

    // Calculate the standard deviation of the weights for the errorbar
    mc_xsec_err = std::sqrt( mc_xsec_err / mc_xsec.size() );
    std::cout << "***********************" << std::endl;
    std::cout << "Calculated GENIE Sys error on mc xsec\t" << mc_xsec_err << std::endl; 
    std::cout << "***********************" << std::endl;
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                  Make TGraph for MC x section
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
    TCanvas * c_integrated = new TCanvas();
	c_integrated->cd();
    
    double x_MC[2] = {0.9, 1.1};
	double y_MC[2] = {mc_xsec_CV,mc_xsec_CV};
	double ex_MC[2] = {0.0,0.0};
	double ey_MC[2] = {mc_xec_err,mc_xec_err} ;

    TGraphErrors * MCxsec_graph = new TGraphErrors(2,x_MC ,y_MC , ex_MC, ey_MC);
    MCxsec_graph->SetTitle("Integrated Nue + Nue-bar Xsec");
    MCxsec_graph->GetYaxis()->SetTitle("#sigma [cm^2]");
    MCxsec_graph->GetXaxis()->SetLimits(0.9, 1.1);
	MCxsec_graph->SetMinimum(0);        // was: 3.6e-39
    MCxsec_graph->SetMaximum(10e-39);   // was: 7.2e-39
    MCxsec_graph->SetFillColor(kCyan-10);
    MCxsec_graph->Draw("a3");

	TLine * line = new TLine(0.9, mc_xsec, 1.1, mc_xsec); // Plot the CV cross section
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
    
    // Legend doesnt seem to work :(
	TLegend * leg = new TLegend();
    leg->AddEntry(line, "GENIE","f");
	leg->AddEntry(xsec_graph, "Data Cross Section", "f");
	leg->Draw("SAME");

	c_integrated->Print("plots/integrated_xsec.eps");


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                     END
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    gSystem->Exit(1); // Quit ROOT

} // END MAIN FUNCTION


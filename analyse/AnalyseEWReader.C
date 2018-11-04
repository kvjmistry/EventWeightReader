// This program will read in the calculated cross section values after running the eventweightreader module.
// It will then calculate the systematic error of the cross sections calculated and do some plotting.

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

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                                       Function Definitions
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function that reads in the event numbers from a text file and adds those event numbers to a vector. 
void ReadEvents(const char *filename, std::vector<double> &MC_Xsec ){
	
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

// Function that calculates the standard deviation // obsolete here
 double STD_Calculator(std::vector<double> Reweighted_vec, double CV  ){
    double Err{0};

    for (unsigned int i = 0; i < Reweighted_vec.size(); i++ ){ 

        Err +=  (Reweighted_vec[i] - CV) * (Reweighted_vec[i]- CV);  
    }

    return (std::sqrt( Err / Reweighted_vec.size() ) );
}

// Function that adds errors in quaderature. 
double Quadrature(std::vector<double> Reweighted_vec, double CV, std::vector<std::string> &GenieNames ){
    double Err{0};
    
    std::vector<double> diff ;           // Difference 
    std::vector<double> maxdiff ;        // The biggest difference out of pm 1 sigma
    std::vector<std::string> Model ;     // Model names slimmed down to correct s

    for (int i = 0 ; i < Reweighted_vec.size(); i++ ) diff.push_back(std::abs(Reweighted_vec[i] - CV)); // Fill diff vector

    // Loop over and get the max error and add to an array, also fill model vec if it is not empty.
    for (int i = 0 ; i < Reweighted_vec.size()/2; i++ ){
        int j = i + 1;
        if (diff[ ( j*2 ) - 2 ] > diff[ j*2 -1 ]) {             // +1 sig > -1 Sig
            maxdiff.push_back(  diff[ (j*2) - 2 ]);
        }
        else {
            maxdiff.push_back( diff[ j*2 -1 ]);       // -1 Sig > + 1 Sig

        }

        Model.push_back(GenieNames[j*2 -1]);

        // Erase the genie name in the string
        Model[i].erase(0, 6);
        Model[i].erase(Model[i].size()-6);

        // std::cout  << diff[ ( j*2 ) - 2 ] << "  " << diff[ j*2 -1 ] << "  " << maxdiff[i] << std::endl;
        
        // Add errors in quadrature, but only for ones which are wanted
        if (  Model[i] == "NC" || Model[i] == "qevec" || Model[i] == "ResDecayEta" || Model[i] == "FermiGasModelSf" || Model[i] == "all") {
            continue;
        }

        
        Err+= maxdiff[i] * maxdiff[i] ; // square and add to error calcualtion


    }

    return ( std::sqrt(Err) );

}

void plotError(std::vector<double> Reweighted_vec, double CV, std::vector<std::string> &GenieNames,const char * print_name,const char *title, const char * varname){
    std::vector<double> diff ;           // Difference 
    std::vector<double> pSig ;           // plus 1 sigma
    std::vector<double> mSig ;           // minus 1 sigma
    std::vector<std::string> Model ;     // Model names slimmed down to correct size. 

    double pErr{0}; //  +1  Sigma percentage error
    double mErr{0}; //  -1  Sigma percentage error


    for (int i = 0 ; i < Reweighted_vec.size(); i++ ) diff.push_back( Reweighted_vec[i] - CV ); // Fill abs of diff vector


    // Loop over and get the max error and add to an array, also fill model vec.
    for (int i = 0 ; i < Reweighted_vec.size()/2; i++ ){
        int j = i + 1;  
        pSig.push_back(  diff[ (j*2) - 2 ]  );  // +1 sig 
        mSig.push_back(  diff[ j*2 -1 ]     );  // -1 Sig
        Model.push_back( GenieNames[  (j*2) - 2  ] );    // Slimmed model

        // Erase the genie name in the string
        Model[i].erase(0, 6);
        Model[i].erase(Model[i].size()-6);

        // Add errors in quadrature, but only for ones which are wanted
        if (  Model[i] == "NC" || Model[i] == "qevec" || Model[i] == "ResDecayEta" || Model[i] == "FermiGasModelSf" || Model[i] == "all") {
            continue;
        }
        // std::cout << Model[i] << std::endl;
        pErr+=  ( (100 *  diff[ (j*2) - 2 ]) / CV ) * ( (100 *  diff[ (j*2) - 2 ]) / CV ) ;
        mErr+=  ( (100 *  diff[ j*2 -1 ]) / CV ) * ( (100 *  diff[ j*2 -1 ]) / CV ); 

        //std::cout  << pSig[i] << " " << mSig[i]<< std::endl;
    }

    

    // Create histogram
    int num_bins = Model.size() - 2 ;
    TH1D *hPlus = new TH1D("hPlus",title, num_bins, 0, num_bins);
    TH1D *hMinus = new TH1D("hMinus",title, num_bins, 0, num_bins);
    TH1D *hCV = new TH1D("hCV",title, num_bins, 0, num_bins);



    // Fill histogram choosing to only fill the x-labels with a sizeable error 
    std::string blank = " ";
    for (int i = 0 ; i < pSig.size(); i++ ) {
        if ( Model[i] == "NC" || Model[i] == "qevec" || Model[i] == "ResDecayEta" || Model[i] == "FermiGasModelSf" || Model[i] == "all") continue;

        if ( std::abs((100 * pSig[i]) / CV ) >= 0.0){
            hPlus->Fill(Model[i].c_str(), (100 * pSig[i]) / CV   );
            hMinus->Fill(Model[i].c_str(), (100 * mSig[i]) / CV );
            hCV->Fill(Model[i].c_str(), (100* CV) / CV - 100 );

        }
       /*  
        else {
            hPlus->Fill(blank.c_str(), (100 * pSig[i]) / CV   );
            hMinus->Fill(blank.c_str(), (100 * mSig[i]) / CV );
            hCV->Fill(blank.c_str(), (100* CV) / CV - 100 );
            

        }
        blank = blank + " ";  
         */
    }

    // Fill a space between the models and total
    hPlus->Fill(" ", 0);
    hMinus->Fill(" ", 0);
    hCV->Fill(" ", 0);
    
    hPlus->Fill("Total Error", std::sqrt(pErr));
    hMinus->Fill("Total Error", -1 * std::sqrt(mErr));
    hCV->Fill("Total Error", 0);


    // Create canvas and choose histogram width 
    TCanvas *c_Percent = new TCanvas();
	c_Percent->cd();
    gStyle->SetOptStat(0);
    gPad->SetBottomMargin(0.32);

    int x_range_u{num_bins};
    float x_range_d{0};
    float y_range_u{15};
    float y_range_d{-15};

    hPlus->SetLineColor(kGreen+2);
    hPlus->SetAxisRange(x_range_d,x_range_u, "X");
    hPlus->SetAxisRange(y_range_d,y_range_u, "Y");
    hPlus->LabelsOption("v");
    
    hMinus->SetLineColor(kRed+2);
    hMinus->SetAxisRange(x_range_d,x_range_u, "X");
    hMinus->SetAxisRange(y_range_d,y_range_u, "Y");
    hMinus->LabelsOption("v");

    hCV->SetLineColor(kBlack);
    hCV->SetAxisRange(x_range_d,x_range_u, "X");
    hCV->SetAxisRange(y_range_d,y_range_u, "Y");
    hCV->LabelsOption("v");
    
    // Create the Legend
    TLegend *legend = new TLegend(0.15,0.75,0.3,0.85);
    legend->AddEntry(hPlus,"+1 #sigma","l");
    legend->AddEntry(hMinus,"-1 #sigma","l");
    
    // Draw and print the histogram
    hCV ->Draw("hist");
    hPlus ->Draw("hist, same");
    hMinus ->Draw("hist, same");
    legend->Draw();
    
    c_Percent->Print(print_name);
    c_Percent->Close();

    std::cout << "-----------------------" << std::endl;
    std::cout << "Calculated GENIE Sys error on  " <<  varname  << ":\t"
    << CV << "\t+"<< std::sqrt(pErr)*CV /100 <<"\t-"<< std::sqrt(mErr)*CV /100  << std::setprecision(3) <<std::endl; 
    std::cout << "-----------------------" << std::endl;

    delete hPlus;
    delete hMinus;
    delete hCV;

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
    const double full_error{};
    
    // Central Values
    const double data_xsec_CV {}; // Data Cross section 4.84338e-39 (from coltons output, seems to be unreproducable)
    const double mc_xsec_CV{};    // central value for the x section
    const double Sel_CV{};
    const double Gen_CV{};
    const double Sig_CV{};
    const double Bkg_CV{};
    const double Eff_CV{Sig_CV/Gen_CV};
    
    // Errors
    double mc_xsec_err;   // systematic error from all the cross sections in MC
    double data_xsec_err; // systematic error from all the cross sections in Data
    double Sel_err;   // systematic error of the Sel events
    double Sig_err;   // systematic error of the Sig events
    double Gen_err;   // systematic error of the Gen events
    double Bkg_err;   // systematic error of the Bkg events
    double eff_err;   // systematic error of the efficiency

    // Other Variables
    std::vector<double> mc_xsec;
    std::vector<double> data_xsec;

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                  Read In TTree and fill variables for use later
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	std::vector<double> N_gen, N_sig, N_bkg, N_sel, Efficiency;
    std::vector<std::string> GenieNames;
    
    //Open the NuMIEventWeight_module_out File
	TFile *fileIN = TFile::Open("NuMIEventWeight_module_out.root");
	if (fileIN == 0) {
	  // if we cannot open the file, print an error message and return immediatly
	  printf("Error: cannot open NuMIEventWeight_module_out.root!\n");
	  return;
	}
    
	// Create the tree reader and its data containers
    TTreeReader myReader("microboonewvtof/XSectionTree", fileIN);

    // Vectors contain the various variables recalculated for each universe. 
	TTreeReaderValue<std::vector<double>> Sig(myReader, "Sig");
    TTreeReaderValue<std::vector<double>> Bkg(myReader, "Bkg");
    TTreeReaderValue<std::vector<double>> Sel(myReader, "Sel");
    TTreeReaderValue<std::vector<double>> Gen(myReader, "Gen");
    TTreeReaderValue<std::vector<double>> eff(myReader, "eff");
    TTreeReaderValue<std::vector<std::string>> GenieNamesRV(myReader, "GenieNames");
    TTreeReaderValue<std::vector<double>> MCXSec(myReader, "MCXSec");
    TTreeReaderValue<std::vector<double>> Data_x_sec(myReader, "DataXSec");
    
    // Create vectors which we can use to plot and calculate with
    while (myReader.Next()) {
        for (int k = 0; k < (*Gen).size(); k++){
            N_gen.push_back( (*Gen)[k] );
            N_sig.push_back((*Sig)[k]);
            N_sel.push_back((*Sel)[k]);
            N_bkg.push_back((*Bkg)[k]);
            Efficiency.push_back((*eff)[k]);
            mc_xsec.push_back( (*MCXSec)[k]);
            data_xsec.push_back((*Data_x_sec)[k]);
            GenieNames.push_back((*GenieNamesRV)[k]);
        }
    }

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                  Do the same again but for the Multisim cross sections
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	std::vector<double> N_gen_msim, N_sig_msim, N_bkg_msim, Efficiency_msim, MC_xsec_msim, Data_xsec_msim; // msim for multisim 
    
    //Open the NuMIEventWeight_module_out File
	TFile *fileIN_msim = TFile::Open("NuMIEventWeight_module_out_multisim.root");
	if (fileIN_msim == 0) {
	  // if we cannot open the file, print an error message and return immediatly
	  printf("Error: cannot open NuMIEventWeight_module_out_multisim.root!\n");
	  return;
	}
    
	// Create the tree reader and its data containers
    TTreeReader myReader_msim("microboonewvtof/XSectionTree", fileIN_msim);

    // Vectors contain the various variables recalculated for each universe. 
	TTreeReaderValue<std::vector<double>> Sig_msim(myReader_msim, "Sig");
    TTreeReaderValue<std::vector<double>> Bkg_msim(myReader_msim, "Bkg");
    TTreeReaderValue<std::vector<double>> Gen_msim(myReader_msim, "Gen");
    TTreeReaderValue<std::vector<double>> eff_msim(myReader_msim, "eff");
    TTreeReaderValue<std::vector<double>> MCXSec_msim(myReader_msim, "MCXSec");
    TTreeReaderValue<std::vector<double>> Data_x_sec_msim(myReader_msim, "DataXSec");
    
    // Create vectors which we can use to plot and calculate with
    while (myReader_msim.Next()) {
        for (int k = 0; k < (*Gen_msim).size(); k++){
            N_gen_msim.push_back( (*Gen_msim)[k] );
            N_sig_msim.push_back((*Sig_msim)[k]);
            N_bkg_msim.push_back((*Bkg_msim)[k]);
            Efficiency_msim.push_back((*eff_msim)[k]);
            MC_xsec_msim.push_back((*MCXSec_msim)[k]);
            Data_xsec_msim.push_back((*Data_x_sec_msim)[k]);
        }
    }

    // Now calculate the standard deviations of each of the quantities and then 

    double msim_gen, msim_sig, msim_bkg,msim_eff, msim_mcxsec, msim_dataxsec; 
    
    msim_gen        = STD_Calculator(N_gen_msim, Gen_CV);
    msim_bkg        = STD_Calculator(N_bkg_msim, Bkg_CV);
    msim_sig        = STD_Calculator(N_sig_msim, Sig_CV);
    msim_mcxsec     = STD_Calculator(MC_xsec_msim, mc_xsec_CV);
    msim_dataxsec   = STD_Calculator(Data_xsec_msim, data_xsec_CV);
    msim_eff        = STD_Calculator(Efficiency_msim, Eff_CV);
   

    std::cout << "-------MULTISIM------" << std::setprecision(3) << std::endl;
    std::cout << "Calculated GENIE Sys error on Gen:\t\t" << msim_gen << " -->\t\t" <<100 * msim_gen / Gen_CV << std::setprecision(3) <<" %" <<std::endl; 
    std::cout << "Calculated GENIE Sys error on Sig:\t\t" << msim_sig << " -->\t" <<100 * msim_sig / Sig_CV<< std::setprecision(3) <<" %" <<std::endl; 
    std::cout << "Calculated GENIE Sys error on Bkg:\t\t" << msim_bkg << " -->\t" << 100 * msim_bkg / Bkg_CV<< std::setprecision(3) <<" %" <<std::endl; 
    std::cout << "Calculated GENIE Sys error on Eff:\t\t" << msim_eff << " -->\t" <<100 * msim_eff / Eff_CV<< std::setprecision(3) <<" %" <<std::endl; 
    std::cout << "Calculated GENIE Sys error on MC XSEC:\t\t" << msim_mcxsec<< " -->\t" <<100 * msim_mcxsec / mc_xsec_CV<< std::setprecision(3) <<" %" <<std::endl; 
    std::cout << "Calculated GENIE Sys error on DATA XSEC:\t" << msim_dataxsec << " -->\t" <<100 * msim_dataxsec / data_xsec_CV<< std::setprecision(3) <<" %" <<std::endl; 
    std::cout << "-----------------------" << std::endl;


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                  Loop over cross sectons and calculate the cross section error.
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	


    // Calculate the standard deviation of the weights for the errorbar
    mc_xsec_err   = Quadrature(mc_xsec, mc_xsec_CV, GenieNames);
    data_xsec_err = Quadrature(data_xsec, data_xsec_CV, GenieNames);
    
    std::cout << "*********UNISIM********" << std::endl;
    std::cout << "Calculated GENIE Sys error on  * MC *  xsec\t" << mc_xsec_err   << "\t-->\t"<< 100 * mc_xsec_err/mc_xsec_CV <<" %" <<std::endl; 
    std::cout << "Calculated GENIE Sys error on * DATA * xsec\t" << data_xsec_err << "\t-->\t"<< 100 * data_xsec_err/data_xsec_CV <<" %" <<std::endl; 
    std::cout << "***********************" << std::endl;
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                  Make TGraph for MC x section
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
    TCanvas * c_integrated = new TCanvas();
	c_integrated->cd();
    
    double x_MC[2] = {0.9, 1.1};
	double y_MC[2] = {mc_xsec_CV,mc_xsec_CV};
	double ex_MC[2] = {0.0,0.0};
	double ey_MC[2] = {mc_xsec_err,mc_xsec_err} ;

    TGraphErrors *MCxsec_graph = new TGraphErrors(2,x_MC ,y_MC , ex_MC, ey_MC);
    MCxsec_graph->SetTitle("Integrated Nue + Nue-bar Xsec");
    MCxsec_graph->GetYaxis()->SetTitle("#sigma [cm^2]");
    MCxsec_graph->GetXaxis()->SetLimits(0.9, 1.1);
	MCxsec_graph->SetMinimum(0);        // was: 
    MCxsec_graph->SetMaximum(10e-39);   // was: 
    MCxsec_graph->SetFillColor(kCyan-10);
    MCxsec_graph->Draw("a3");

	TLine *line = new TLine(0.9, mc_xsec_CV, 1.1, mc_xsec_CV); // Plot the CV cross section
	line->SetLineColor(kCyan+1);
    line->SetLineWidth(2);
	line->Draw("same");

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                 Create TGraph for x sec DATA
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	

    const int n = 1;
	double x[n] = {1.0};
	double y[n] = {data_xsec_CV};
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

	c_integrated->Print("plots/integrated_xsec.pdf");
    c_integrated->Close();


 
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                    Plot efficiency for each of the models
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // put for loop here to hide comments away ;)
    for (int i=0; i < 1; i++){
        // TFile hfile("hfile.root","RECREATE");

        // // create histograms
        // TH1D *h[GenieNames.size()];
        // char name[20];
        // char title[100];
        // for (int i=0; i < GenieNames.size(); i++) {
        //     std::string fullpath;
        //     std::string index;
        //     index = std::to_string(i);
        //     fullpath = GenieNames[i] + index; 
        //     sprintf(name,fullpath.c_str());
        //     sprintf(title,fullpath.c_str(), i);
        //     h[i] = new TH1D(name,title,1,0,1);
        // }

        // // fill histograms
        // for (int i = 0; i < GenieNames.size() ;i++) {
        //     h[i]->Fill(GenieNames[i].c_str(), Efficiency[i]);
        //     h[i]->SetOption("hist");
        // }

        // TH1D *h_CV[GenieNames.size()];
        // char name1[20];
        // char title1[100];
        // for (int i=0; i < GenieNames.size(); i++) {
        //     std::string fullpath;
        //     std::string index;
        //     index = std::to_string(i);
        //     fullpath = GenieNames[i] + index + "CV";
        //     sprintf(name1,fullpath.c_str());
        //     sprintf(title1,fullpath.c_str());
            
        //     h_CV[i] = new TH1D(name1,title1,1,0,1);
        // }

        // // fill histograms
        // for (int i = 0; i < GenieNames.size() ;i++) {
        //     h_CV[i]->Fill(GenieNames[i].c_str(), Eff_CV);
        //     h_CV[i]->SetOption("hist");
        //     h_CV[i]->SetAxisRange(0.08,0.1, "Y");
        //     h_CV[i]->SetLineColor(kBlack);
        // }

        // TCanvas *c_eff[38];
        // for (int i = 0; i < 38 ;i++) {
        //     c_eff[i] = new TCanvas();
        //     c_eff[i]->cd();
        //     h_CV[i]->Draw();
        //     h[(i*2)-1]->Draw("same");
        //     h[i*2]->Draw("same");
        //     std::string namex;
        //     std::string path = "plots/eff_plots/h";
        //     std::string extension = ".eps"
        //     index = std::to_string(i);
        //     namex =  path + index + extension;
        //     c_integrated->Print(name);
        //     c_eff[i]->Close();
        // }

        // // save histograms
        // hfile.Write();
    }

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                 Plot the difference from CV in a Histogram for easy viewing for each model
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    plotError(mc_xsec, mc_xsec_CV, GenieNames, "plots/MC_X_Sec.pdf","Genie Model Errors on MC Cross Section ; ; Percentage Change from CV [%] ","MC X Sec" );
    plotError(data_xsec, data_xsec_CV, GenieNames, "plots/Data_X_Sec.pdf","Genie Model Errors on Data Cross Section ; ; Percentage Change from CV [%] ","Data X Sec" );
    plotError(N_gen, Gen_CV, GenieNames, "plots/Generated_reweight.pdf","Genie Model Errors on Generated Events ; ; Percentage Change from CV [%] ","Gen" );
    plotError(N_sig, Sig_CV, GenieNames, "plots/Signal_reweight.pdf","Genie Model Errors on Sel Signal Events ; ; Percentage Change from CV [%] ","Sig" );
    plotError(N_sel, Sel_CV, GenieNames, "plots/Selected_reweight.pdf","Genie Model Errors on Selected Events ; ; Percentage Change from CV [%] ","Sel" );
    plotError(N_bkg, Bkg_CV, GenieNames, "plots/Background_reweight.pdf","Genie Model Errors on Background Events ; ; Percentage Change from CV [%] ","Bkg" );
    plotError(Efficiency, Eff_CV, GenieNames, "plots/Efficiency_reweight.pdf","Genie Model Errors on Efficiency ; ; Percentage Change from CV [%] ","Eff" );
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                             END
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    gSystem->Exit(1); // Quit ROOT

} // END MAIN FUNCTION


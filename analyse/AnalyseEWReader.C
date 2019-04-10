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

// Use this function to append a branch onto a ttree
void treeAddBranch(std::vector<std::vector<double>> new_v_v, const char * name) {
	std::vector<double> new_v;

	TFile *f = new TFile("plots/KrishTree.root", "UPDATE");
	// TFile f("KrishTree.root", "update");


	TTree *DataTree = (TTree*)f->Get("DataTree");

	TBranch *newBranch = DataTree->Branch(name, &new_v); 

	// read the number of entries in the tree, should be equal to the number of parameters
	Long64_t nentries = DataTree->GetEntries(); 

	// Loop over the vector of vector and 
	for (Long64_t i = 0; i < nentries; i++) {
		new_v = new_v_v[i];
		newBranch->Fill();
		new_v.clear();
	}

	DataTree->Write("", TObject::kOverwrite); // save only the new version of the tree
	f->Close();
}

// ++++++++UNISIM FUNCTIONS++++++++
// Function that calculates the standard deviation
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
		if ( 
				Model[i] == "NC" || Model[i] == "qevec" || Model[i] == "ResDecayEta" || Model[i] == "FermiGasModelSf" || Model[i] == "all" || Model[i] == "AGKYpT" || Model[i] == "AGKYxF" ||
				Model[i] == "DISAth" || Model[i] == "DISBth" || Model[i] == "DISCv1u" || Model[i] == "DISCv2u"
				// Model[i] == "FermiGasModelKf" || Model[i] == "FormZone"
				// Model[i] == "qema" || Model[i] == "ccresAxial"
				// Model[i] == "ncresVector" || Model[i] == "ncresAxial" 
				// Model[i] == "ccresVector" || Model[i] == "ccresAxial" ||
				// Model[i] == "IntraNukeNmfp" || Model[i] == "IntraNukeNcex" || Model[i] == "IntraNukeNel" || Model[i] == "IntraNukeNinel" ||
				// Model[i] == "IntraNukeNabs" || Model[i] == "IntraNukeNpi" || Model[i] == "IntraNukePImfp" || Model[i] == "IntraNukePIcex" || 
				// Model[i] == "IntraNukePIel" || Model[i] == "IntraNukePIinel" || Model[i] == "IntraNukePIabs" || Model[i] == "IntraNukePIpi"
				// Model[i] == "NonResRvp1pi" || Model[i] == "NonResRvbarp1pi" || Model[i] == "NonResRvp2pi" || Model[i] == "NonResRvbarp2pi" 
		   ) 
		{

			continue;
		}
		// else continue;
		// std::cout << Model[i] <<std::endl;


		Err+= maxdiff[i] * maxdiff[i] ; // square and add to error calcualtion


	}

	return ( std::sqrt(Err) );

}

// Function that plots the pm1sigma variations on the same plot 
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
		if ( 
				//Model[i] == "DISAth" || Model[i] == "DISBth" || Model[i] == "DISCv1u" || Model[i] == "DISCv2u" || 
				Model[i] == "NC" || Model[i] == "qevec" || Model[i] == "ResDecayEta" || Model[i] == "FermiGasModelSf" || Model[i] == "all") {
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
		if ( 
				//Model[i] == "DISAth" || Model[i] == "DISBth" || Model[i] == "DISCv1u" || Model[i] == "DISCv2u" || 
				Model[i] == "NC" || Model[i] == "qevec" || Model[i] == "ResDecayEta" || Model[i] == "FermiGasModelSf" || Model[i] == "all") continue;

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

void FillIndex(std::vector<double> Reweighted_vec, int N_Univ, int &index, TTree *DataTree){
	// Loop by the number of parameters times
	for (int i = 0 ; i < Reweighted_vec.size()/N_Univ; i++ ){
		DataTree->Fill();
		index ++;
	}

}

// +++++++MULTISUM FUNCTIONS+++++++
// Use this function to calculate the standard deviation for multiple parameters not one and return a vector 
std::vector<double> STD_Calculator(std::vector<double> Reweighted_vec, double CV, int N_Univ, const char * name ){
	std::vector<double> Std_Vec;

	std::vector<std::vector<double>> temp_vec; // Create a vector of vector to loop

	// Loop by the number of parameters times
	for (int i = 0 ; i < Reweighted_vec.size()/N_Univ; i++ ){
		int j = i + 1;  

		std::vector<double> temp; // temp vector for std calculation

		// Do the standard deviation calculation for each parameter
		for (int k=0; k<N_Univ; k++){ temp.push_back( Reweighted_vec[  j * N_Univ - N_Univ + k    ] );}

		temp_vec.push_back(temp); // add the vector to the loop

		// Std_Vec.push_back( STD_Calculator( temp, CV )); // Return Values
		Std_Vec.push_back( 100 * STD_Calculator( temp, CV ) / CV); // Return Percentages

		temp.clear(); // make sure temp is empty

	}

	treeAddBranch(temp_vec, name );
	temp_vec.clear();

	return Std_Vec;
}

// Sum in quadrature the std of the multisims 
double Quadrature(std::vector<double> errors){

	double Err{0};
	for (int i = 0; i < errors.size() ;i++){
		Err+= errors[i] * errors[i] ;

	}

	return (std::sqrt(Err));
}

// Overloaded function to plot the multisim errors on top for comparison. 
void plotError(std::vector<double> Reweighted_vec, double CV, std::vector<std::string> &GenieNames,const char * print_name,const char *title, const char * varname, std::vector<double> errors){
	std::vector<double> diff ;           // Difference 
	std::vector<double> pSig ;           // plus 1 sigma
	std::vector<double> mSig ;           // minus 1 sigma
	std::vector<std::string> Model ;     // Model names slimmed down to correct size. 

	std::vector<std::string> Model_msim ; // Smaller not including the unsimulated parameters

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
		if ( 
				//Model[i] == "DISAth" || Model[i] == "DISBth" || Model[i] == "DISCv1u" || Model[i] == "DISCv2u" || 
				Model[i] == "NC" || Model[i] == "qevec" || Model[i] == "ResDecayEta" || Model[i] == "FermiGasModelSf" || Model[i] == "all" || Model[i] == "AGKYpT" || Model[i] == "AGKYxF") {
			continue;
		}
		// std::cout << Model[i] << std::endl;
		pErr+=  ( (100 *  diff[ (j*2) - 2 ]) / CV ) * ( (100 *  diff[ (j*2) - 2 ]) / CV ) ;
		mErr+=  ( (100 *  diff[ j*2 -1 ]) / CV ) * ( (100 *  diff[ j*2 -1 ]) / CV ); 
		Model_msim.push_back(Model[i]);

		//std::cout  << pSig[i] << " " << mSig[i]<< std::endl;
	}

	// Create histogram
	int num_bins = Model.size() - 4 ;
	TH1D *hPlus = new TH1D("hPlus",title, num_bins, 0, num_bins);
	TH1D *hMinus = new TH1D("hMinus",title, num_bins, 0, num_bins);
	TH1D *hCV = new TH1D("hCV",title, num_bins, 0, num_bins);
	TH1D *hErr = new TH1D("hErr",title, num_bins, 0, num_bins);

	std::cout << "Removing the first two elements of the errors array (which should be AGKY models, delete if no longer needed)"<< std::endl;
	errors.erase (errors.begin(),errors.begin()+2);

	// Set the bin errors of the multisim standard deviation to compare with pm1sigma
	for (int i = 0 ; i < Model_msim.size() ; i++ ) {
		hErr->Fill(Model_msim[i].c_str(), 0 );
		hErr->SetBinError(i+1, errors[i]); // i +1 otherwise the bin errors wont align

	}

	// Calculate the total error for the multisims too
	double tot_err_msim = Quadrature(errors);

	hErr->Fill("  ", 0);
	hErr->Fill("Total Error", 0);
	hErr->SetBinError(num_bins, tot_err_msim);

	// Fill histogram choosing to only fill the x-labels with a sizeable error 
	// std::string blank = " ";
	for (int i = 0 ; i < pSig.size(); i++ ) {
		if ( 
				//Model[i] == "DISAth" || Model[i] == "DISBth" || Model[i] == "DISCv1u" || Model[i] == "DISCv2u" || 
				Model[i] == "NC" || Model[i] == "qevec" || Model[i] == "ResDecayEta" || Model[i] == "FermiGasModelSf" || Model[i] == "all" || Model[i] == "AGKYpT" || Model[i] == "AGKYxF") continue;

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



	// for (int i = 0; i < errors.size(); i++) hErr->SetBinError(i, errors[i]) ;

	// Create canvas and choose histogram width 
	TCanvas *c_Percent = new TCanvas();
	c_Percent->cd();
	gStyle->SetOptStat(0);
	gPad->SetBottomMargin(0.2);

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

	hErr->SetAxisRange(x_range_d,x_range_u, "X");
	hErr->SetAxisRange(y_range_d,y_range_u, "Y");
	hErr->LabelsOption("v");
	hErr->SetMarkerStyle(20);
	hErr->SetMarkerSize(0.1);

	// Create the Legend
	TLegend *legend = new TLegend(0.15,0.75,0.3,0.85);
	legend->AddEntry(hPlus,"+1 #sigma","l");
	legend->AddEntry(hMinus,"-1 #sigma","l");
	legend->AddEntry(hErr,"250 univ multisim STD","e");

	// Draw and print the histogram
	hCV ->Draw("hist");
	hPlus ->Draw("hist, same");
	hMinus ->Draw("hist, same");
	hErr->Draw("PE1, same");

	legend->Draw();

	c_Percent->Print(print_name);
	c_Percent->Close();

	std::cout << "-----------------------" << std::setprecision(3) << std::endl;
	std::cout << "Calculated GENIE Sys error on  " <<  varname  << ":\t"
		<< CV << "\t+"<< std::sqrt(pErr)*CV /100 <<"\t-"<< std::sqrt(mErr)*CV /100 << "\ttot err msim:\t" << tot_err_msim<<  " %"<< std::setprecision(3) <<std::endl; 
	std::cout << "-----------------------" << std::endl;

	delete hPlus;
	delete hMinus;
	delete hCV;
	delete hErr;

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

	// Central Values
	const double data_xsec_CV {4.89936e-39}; // Data Cross section 4.84338e-39 (from coltons output, seems to be unreproducable)
	const double mc_xsec_CV{4.83114e-39};    // central value for the x section
	const double Sel_CV{948.};
	const double Gen_CV{7103.};
	const double Sig_CV{619.};
	const double Bkg_CV{329.};
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
	int Universes{250};
	std::cout << "WARNING the number of universes is:\t" << Universes << "\t This value is hardcoded and should be checked" << std::endl;

	// Multisim branches
	std::vector<double> msim_gen, msim_sig, msim_bkg,msim_eff, msim_mcxsec, msim_dataxsec; 




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

	// Open the NuMIEventWeight_module_out File, choose one of these
	// NuMIEventWeight_module_out.root
	// NuMIEventWeight_module_out_multisim.root // 250 universes with "all" always no DIS
	// NuMIEventWeight_module_out_multisim_100univ.root
	// NuMIEventWeight_module_out_multisim_250univ_all_noFSI.root
	// NuMIEventWeight_module_out_multisim_individual.root
	// NuMIEventWeight_module_out_multisim_withDIS.root
	// NuMIEventWeight_module_out_multisim_250univ_all_noNonRes.root
	// NuMIEventWeight_module_out_multisim_250univ_all_onlyFSI.root
	// NuMIEventWeight_module_out_multisim_250univ_all_only_FGKf_FormZ.root
	// NuMIEventWeight_module_out_multisim_250univ_all_only_QEMA_CCRes.root
	// NuMIEventWeight_module_out_multisim_250univ_all_only_NonRes.root
	// NuMIEventWeight_module_out_multisim_250univ_all_only_CCRes.root
	// NuMIEventWeight_module_out_multisim_250univ_all_noFSI_CCRes.root
	// NuMIEventWeight_module_out_multisim_250univ_all_noFSI_CCRes_NCRes.root
	// NuMIEventWeight_module_out_multisim_250univ_all_only_NonRes.root
	TFile *fileIN_msim = TFile::Open("NuMIEventWeight_module_out_multisim_individual.root");
	if (fileIN_msim == 0) {
		// if we cannot open the file, print an error message and return immediatly
		printf("Error: cannot open NuMIEventWeight_module_out_multisim_individual.root!\n");
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

	// TTree definition
	TFile *Krishfile = new TFile("plots/KrishTree.root","RECREATE");
	TTree *DataTree = new TTree("DataTree", "DataTree");

	// Index and size the ttree for reading an writing, these label the genie names in order of their display on the plot
	int index{1};
	TBranch *bintex = DataTree->Branch("index", &index);
	FillIndex(N_gen_msim, Universes, index, DataTree);

	// Write the indexes in place
	DataTree->Write(); 
	Krishfile->Close();

	std::cout << "Size of the multisim gen vector:\t" << N_gen_msim.size() << std::endl;

	// Now calculate the standard deviations of each of the quantities ++MULTSIMS++
	msim_gen        = STD_Calculator(N_gen_msim, Gen_CV, Universes, "Gen");
	msim_bkg        = STD_Calculator(N_bkg_msim, Bkg_CV, Universes, "Bkg");
	msim_sig        = STD_Calculator(N_sig_msim, Sig_CV, Universes, "Sig");
	msim_eff        = STD_Calculator(Efficiency_msim, Eff_CV, Universes, "Eff");
	msim_mcxsec     = STD_Calculator(MC_xsec_msim, mc_xsec_CV, Universes, "MCxsec");
	msim_dataxsec   = STD_Calculator(Data_xsec_msim, data_xsec_CV, Universes, "DataXSec");

	// Print results for genie all multisim mode only
	if (msim_mcxsec.size() == 1){
		std::cout << "-------MULTISIM GENIE ALL------" << std::setprecision(3) << std::endl;
		std::cout << "Number of Universes Simulated:\t" << MC_xsec_msim.size()  << std::endl;
		std::cout << "                                  \t\t" << "Value" << "\t\t" <<"Percentage"<< std::setprecision(3) <<" %" <<std::endl; 
		std::cout << "Calculated GENIE Sys error on Gen:\t\t" << msim_gen[0] << " -->\t\t" <<100 * msim_gen[0] / Gen_CV << std::setprecision(3) <<" %" <<std::endl; 
		std::cout << "Calculated GENIE Sys error on Sig:\t\t" << msim_sig[0] << " -->\t" <<100 * msim_sig[0] / Sig_CV<< std::setprecision(3) <<" %" <<std::endl; 
		std::cout << "Calculated GENIE Sys error on Bkg:\t\t" << msim_bkg[0] << " -->\t" << 100 * msim_bkg[0] / Bkg_CV<< std::setprecision(3) <<" %" <<std::endl; 
		std::cout << "Calculated GENIE Sys error on Eff:\t\t" << msim_eff[0] << " -->\t" <<100 * msim_eff[0] / Eff_CV<< std::setprecision(3) <<" %" <<std::endl; 
		std::cout << "Calculated GENIE Sys error on MC XSEC:\t\t" << msim_mcxsec[0]<< " -->\t" <<100 * msim_mcxsec[0] / mc_xsec_CV<< std::setprecision(3) <<" %" <<std::endl; 
		std::cout << "Calculated GENIE Sys error on DATA XSEC:\t" << msim_dataxsec[0] << " -->\t" <<100 * msim_dataxsec[0] / data_xsec_CV<< std::setprecision(3) <<" %" <<std::endl; 
		std::cout << "-----------------------" << std::endl;
	}


	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//                                                  Loop over cross sectons and calculate the cross section error.
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	

	// Calculate the standard deviation of the weights for the errorbar
	mc_xsec_err   = Quadrature(mc_xsec, mc_xsec_CV, GenieNames);
	data_xsec_err = Quadrature(data_xsec, data_xsec_CV, GenieNames);

	std::cout << "*********UNISIM********" << std::endl;
	std::cout << "Calculated GENIE Sys error on  * MC *  xsec\t" << std::setprecision(5) << mc_xsec_err   << "\t-->\t"<< 100 * mc_xsec_err/mc_xsec_CV <<" %" <<std::endl; 
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
	MCxsec_graph->SetMinimum(0);        // was: 3.6e-39
	MCxsec_graph->SetMaximum(10e-39);   // was: 7.2e-39
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
	//                                 Plot the difference from CV in a Histogram for easy viewing for each model
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if (msim_mcxsec.size() == 1){ // If only using unisim
		plotError(mc_xsec, mc_xsec_CV, GenieNames, "plots/MC_X_Sec.pdf","MC Cross Section ; ; Percentage Change from CV [%] ","MC X Sec" );
		plotError(data_xsec, data_xsec_CV, GenieNames, "plots/Data_X_Sec.pdf","Data Cross Section ; ; Percentage Change from CV [%] ","Data X Sec" );
		plotError(N_gen, Gen_CV, GenieNames, "plots/Generated_reweight.pdf","Generated Events ; ; Percentage Change from CV [%] ","Gen" );
		plotError(N_sig, Sig_CV, GenieNames, "plots/Signal_reweight.pdf","Sel Signal Events ; ; Percentage Change from CV [%] ","Sig" );
		plotError(N_sel, Sel_CV, GenieNames, "plots/Selected_reweight.pdf","Selected Events ; ; Percentage Change from CV [%] ","Sel" );
		plotError(N_bkg, Bkg_CV, GenieNames, "plots/Background_reweight.pdf","Background Events ; ; Percentage Change from CV [%] ","Bkg" );
		plotError(Efficiency, Eff_CV, GenieNames, "plots/Efficiency_reweight.pdf","Efficiency ; ; Percentage Change from CV [%] ","Eff" );
	}
	else{ // if you want to plot the multisim variations on top of each parameter
		plotError(mc_xsec, mc_xsec_CV, GenieNames, "plots/MC_X_Sec.pdf","MC Cross Section ; ; Percentage Change from CV [%] ","MC X Sec", msim_mcxsec );
		plotError(data_xsec, data_xsec_CV, GenieNames, "plots/Data_X_Sec.pdf","Data Cross Section ; ; Percentage Change from CV [%] ","Data X Sec", msim_dataxsec  );
		plotError(N_gen, Gen_CV, GenieNames, "plots/Generated_reweight.pdf","Generated Events ; ; Percentage Change from CV [%] ","Gen", msim_gen  );
		plotError(N_sig, Sig_CV, GenieNames, "plots/Signal_reweight.pdf","Sel Signal Events ; ; Percentage Change from CV [%] ","Sig", msim_sig  );
		plotError(N_sel, Sel_CV, GenieNames, "plots/Selected_reweight.pdf","Selected Events ; ; Percentage Change from CV [%] ","Sel"  );
		plotError(N_bkg, Bkg_CV, GenieNames, "plots/Background_reweight.pdf","Background Events ; ; Percentage Change from CV [%] ","Bkg", msim_bkg  );
		plotError(Efficiency, Eff_CV, GenieNames, "plots/Efficiency_reweight.pdf","Efficiency ; ; Percentage Change from CV [%] ","Eff", msim_eff  );    

	}

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//                                                             END
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	gSystem->Exit(1); // Quit ROOT

} // END MAIN FUNCTION


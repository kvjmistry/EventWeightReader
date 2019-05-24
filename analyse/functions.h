
// Histogram and canvas definitions
TH1D *hCV_genie_Sig  = new TH1D("hCV_genie_Sig", "Genie Sig;; Percentage Uncertainty [%]",   38, 0, -1);
TH1D *hCV_genie_Bkg  = new TH1D("hCV_genie_Bkg", "Genie Bkg;; Percentage Uncertainty [%]",   38, 0, -1);
TH1D *hCV_genie_Gen  = new TH1D("hCV_genie_Gen", "Genie Gen;; Percentage Uncertainty [%]",   38, 0, -1);
TH1D *hCV_genie_eff  = new TH1D("hCV_genie_eff", "Genie Eff;; Percentage Uncertainty [%]",   38, 0, -1);
TH1D *hCV_genie_xsec = new TH1D("hCV_genie_xsec","Genie X Sec;; Percentage Uncertainty [%]", 38, 0, -1);

TH1D *hCV_genie_Sig_p1sig  = new TH1D("hCV_genie_Sig_p1sig", "Genie Sig;; Percentage Uncertainty [%]",   38, 0, -1);
TH1D *hCV_genie_Bkg_p1sig  = new TH1D("hCV_genie_Bkg_p1sig", "Genie Bkg;; Percentage Uncertainty [%]",   38, 0, -1);
TH1D *hCV_genie_Gen_p1sig  = new TH1D("hCV_genie_Gen_p1sig", "Genie Gen;; Percentage Uncertainty [%]",   38, 0, -1);
TH1D *hCV_genie_eff_p1sig  = new TH1D("hCV_genie_eff_p1sig", "Genie Eff;; Percentage Uncertainty [%]",   38, 0, -1);
TH1D *hCV_genie_xsec_p1sig = new TH1D("hCV_genie_xsec_p1sig","Genie X Sec;; Percentage Uncertainty [%]", 38, 0, -1);

TH1D *hCV_genie_Sig_m1sig  = new TH1D("hCV_genie_Sig_m1sig", "Genie Sig;; Percentage Uncertainty [%]",   38, 0, -1);
TH1D *hCV_genie_Bkg_m1sig  = new TH1D("hCV_genie_Bkg_m1sig", "Genie Bkg;; Percentage Uncertainty [%]",   38, 0, -1);
TH1D *hCV_genie_Gen_m1sig  = new TH1D("hCV_genie_Gen_m1sig", "Genie Gen;; Percentage Uncertainty [%]",   38, 0, -1);
TH1D *hCV_genie_eff_m1sig  = new TH1D("hCV_genie_eff_m1sig", "Genie Eff;; Percentage Uncertainty [%]",   38, 0, -1);
TH1D *hCV_genie_xsec_m1sig = new TH1D("hCV_genie_xsec_m1sig","Genie X Sec;; Percentage Uncertainty [%]", 38, 0, -1);

TCanvas * c_genie_Sig  = new TCanvas();
TCanvas * c_genie_Bkg  = new TCanvas();
TCanvas * c_genie_Gen  = new TCanvas();
TCanvas * c_genie_eff  = new TCanvas();
TCanvas * c_genie_xsec = new TCanvas();

TH1D *hCV_interaction_Sig  = new TH1D("hCV_interaction_Sig","Interactions Sig;; Percentage Uncertainty [%]",    6, 0, 6);
TH1D *hCV_interaction_Bkg  = new TH1D("hCV_interaction_Bkg","Interactions Bkg;; Percentage Uncertainty [%]",    6, 0, 6);
TH1D *hCV_interaction_Gen  = new TH1D("hCV_interaction_Gen","Interactions Gen;; Percentage Uncertainty [%]",    6, 0, 6);
TH1D *hCV_interaction_eff  = new TH1D("hCV_interaction_eff","Interactions Eff;; Percentage Uncertainty [%]",    6, 0, 6);
TH1D *hCV_interaction_xsec = new TH1D("hCV_interaction_xsec","Interactions X Sec;; Percentage Uncertainty [%]", 6, 0, 6);

TCanvas * c_interaction_Sig  = new TCanvas();
TCanvas * c_interaction_Bkg  = new TCanvas();
TCanvas * c_interaction_Gen  = new TCanvas();
TCanvas * c_interaction_eff  = new TCanvas();
TCanvas * c_interaction_xsec = new TCanvas();

TH1D *h_genie_all_univ       = new TH1D("genie_all_univ", "Genie X Sec -- all;Universe;Percentage Uncertainty [%]", 10, 0, 1000); 
TCanvas * c_genie_all_univ   = new TCanvas();


// ------------------------------------------------------------------------------------------------------------
bool GetTree(TFile* f, TTree* &T, TString string){
        T = (TTree*)(f->Get(string));
        if (T == NULL) {
                std::cout << "\nfailed to get:\t" << string << "\tThis tree might not exist in the file\n" << std::endl;
                return false;
        }
        else {
                return true;
        }
}
// ------------------------------------------------------------------------------------------------------------
bool GetFile(TFile* &f , TString string){
        f = TFile::Open(string);

        if (f == NULL) {
                std::cout << "failed to get:\t" << string << "\tThis file might not exist in the file" << std::endl;
                return false;
        }
        else {
                return true;
        }
}
// ------------------------------------------------------------------------------------------------------------
double STD (std::vector<double> N, double CV  ){
	double Err{0};

	for (unsigned int i = 0; i < N.size(); i++ ) Err +=  (N[i] - CV) * (N[i]- CV);  
	
	return (std::sqrt( Err / N.size() ) );
}
// ------------------------------------------------------------------------------------------------------------
void Unisim_err(std::vector<double> N, double CV, double &p1sig, double &m1sig  ){

        // + 1 sigma difference 
        p1sig =  N.at(0) - CV;

        // - 1 sigma difference 
        m1sig =  N.at(1) - CV;
        
	return;
}
// ------------------------------------------------------------------------------------------------------------
// Multisim
void Print_Error(std::string CV_str, double CV, std::string Err_str, double Err){
    std::cout << CV_str << ": " << CV << " "<< Err_str << ": " << 100 * Err / CV << " %"<< std::endl;
}
// Unisim
// ------------------------------------------------------------------------------------------------------------
void Print_Error(std::string CV_str, double CV, double p1sig, double m1sig){
    std::cout << CV_str << ": " << CV << " "<< "+1 sig"<< ": " << 100 * p1sig / CV << " %" << " "<< "-1 sig"<< ": " << 100 * m1sig / CV << " %"<< std::endl;
}
// ------------------------------------------------------------------------------------------------------------
void hist_options(TCanvas* c, TH1D* msim, TH1D* p1sig, TH1D* m1sig, double axis_range, int n_uni, const char *print_name ){
        c ->cd();
	gPad->SetBottomMargin(0.2);

        TLegend *legend = new TLegend(0.15, 0.75, 0.3, 0.85);
	legend->AddEntry(p1sig ,"+1 #sigma","l");
	legend->AddEntry(m1sig, "-1 #sigma","l");
        legend->AddEntry(msim,  Form("%i Univ Multisim STD", n_uni),"e");
        legend->SetTextFont(62);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
	
        p1sig->LabelsOption("av");
	p1sig->GetYaxis()->SetRangeUser(-1 * axis_range, axis_range);
	p1sig->SetLineColor(kGreen+1);
        p1sig->SetLineWidth(2);
        p1sig->SetLabelFont(62);
	p1sig->Draw("hist,same");
        
        m1sig->LabelsOption("av");
        m1sig->GetYaxis()->SetRangeUser(-1 * axis_range, axis_range);
	m1sig->SetLineColor(kRed+1);
        m1sig->SetLineWidth(2);
        m1sig->SetLabelFont(62);
	m1sig->Draw("hist,same");

        msim->LabelsOption("av");
        msim->GetYaxis()->SetRangeUser(-1 * axis_range, axis_range);
        msim->Draw("PE1,same");

        legend->Draw();

        c->Print(print_name);
}
// ------------------------------------------------------------------------------------------------------------
void hist_options(TCanvas* c, TH1D* msim, double axis_range, int n_uni_model, int n_uni_reint, const char *print_name){

        c ->cd();
	gPad->SetBottomMargin(0.25);

        TLegend *legend = new TLegend(0.55, 0.75, 0.85, 0.85);
        legend->AddEntry(msim,  Form("CCMEC/CCQE: %i Univ Multisim STD", n_uni_model),"e");
        legend->AddEntry(msim,  Form("Reinteraction: %i Univ Multisim STD", n_uni_reint),"e");
        legend->SetTextFont(62);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);


        msim->LabelsOption("v");
        // msim->GetYaxis()->SetRangeUser(-1 * axis_range, axis_range);
        msim->SetLabelFont(62);
        msim->Draw("PE1,same");

        legend->Draw();

        c->Print(print_name);


}
// ------------------------------------------------------------------------------------------------------------
void make_genie_universe_plot(TCanvas* c, TH1D* h, const char *print_name ){
        c->cd();
        std::vector<int> universe = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900};
	std::vector<double> uncertainty = { 5.68, 4.98, 4.70, 4.73, 4.60, 4.63, 4.63, 4.62, 4.68, 4.62};

        // Fill the histogram
	for (unsigned i =0; i < universe.size(); i++)
		h->Fill(universe[i], uncertainty[i] );
	

	h->SetLineWidth(2);
	h->SetLineColor(kBlack);

	gPad->SetBottomMargin(0.1);
	h->Draw("hist");

        c->Print(print_name);


}

// Histogram and canvas definitions
TH1D *hCV_genie_Sig  = new TH1D("hCV_genie_Sig", "Genie Sig;; Percentage Uncertainty [%]",   36, 0, -1);
TH1D *hCV_genie_Bkg  = new TH1D("hCV_genie_Bkg", "Genie Bkg;; Percentage Uncertainty [%]",   36, 0, -1);
TH1D *hCV_genie_Gen  = new TH1D("hCV_genie_Gen", "Genie Gen;; Percentage Uncertainty [%]",   36, 0, -1);
TH1D *hCV_genie_eff  = new TH1D("hCV_genie_eff", "Genie Eff;; Percentage Uncertainty [%]",   36, 0, -1);
TH1D *hCV_genie_xsec = new TH1D("hCV_genie_xsec","Genie X Sec;; Percentage Uncertainty [%]", 36, 0, -1);

TH1D *hCV_genie_Sig_p1sig  = new TH1D("hCV_genie_Sig_p1sig", "Genie Sig;; Percentage Uncertainty [%]",   36, 0, -1);
TH1D *hCV_genie_Bkg_p1sig  = new TH1D("hCV_genie_Bkg_p1sig", "Genie Bkg;; Percentage Uncertainty [%]",   36, 0, -1);
TH1D *hCV_genie_Gen_p1sig  = new TH1D("hCV_genie_Gen_p1sig", "Genie Gen;; Percentage Uncertainty [%]",   36, 0, -1);
TH1D *hCV_genie_eff_p1sig  = new TH1D("hCV_genie_eff_p1sig", "Genie Eff;; Percentage Uncertainty [%]",   36, 0, -1);
TH1D *hCV_genie_xsec_p1sig = new TH1D("hCV_genie_xsec_p1sig","Genie X Sec;; Percentage Uncertainty [%]", 36, 0, -1);

TH1D *hCV_genie_Sig_m1sig  = new TH1D("hCV_genie_Sig_m1sig", "Genie Sig;; Percentage Uncertainty [%]",   36, 0, -1);
TH1D *hCV_genie_Bkg_m1sig  = new TH1D("hCV_genie_Bkg_m1sig", "Genie Bkg;; Percentage Uncertainty [%]",   36, 0, -1);
TH1D *hCV_genie_Gen_m1sig  = new TH1D("hCV_genie_Gen_m1sig", "Genie Gen;; Percentage Uncertainty [%]",   36, 0, -1);
TH1D *hCV_genie_eff_m1sig  = new TH1D("hCV_genie_eff_m1sig", "Genie Eff;; Percentage Uncertainty [%]",   36, 0, -1);
TH1D *hCV_genie_xsec_m1sig = new TH1D("hCV_genie_xsec_m1sig","Genie X Sec;; Percentage Uncertainty [%]", 36, 0, -1);

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

TGraph   *g_genie_all_univ   = new TGraph(); 
TCanvas * c_genie_all_univ   = new TCanvas();

// Float_t bins[] = { 0, 1, 1.01, 2, 3, 3.01, 4, 5, 5.01, 6 };
TH1D *h_uncertainties_1       = new TH1D("Uncertainties_1",";; Percentage Uncertainty [%]", 3, 0, 3);
TH1D *h_uncertainties_2       = new TH1D("Uncertainties_2",";; Percentage Uncertainty [%]", 3, 0, 3);
TCanvas * c_uncertainties     = new TCanvas();


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
double Quadrature(std::vector<double> vec){

	double err{0.};
	for (unsigned i=0; i < vec.size(); i++){
		err+= vec[i] * vec[i];
	}
	
	// std::cout << std::sqrt(err) << std::endl;
	return (std::sqrt(err));

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

	// c->Print(print_name);
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

	// c->Print(print_name);


}
// ------------------------------------------------------------------------------------------------------------
void make_genie_universe_plot(TCanvas* c, TGraph* g, const char *print_name ){
	c->cd();
	std::vector<int> universe = { 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
	std::vector<double> uncertainty = { 5.68, 4.98, 4.70, 4.73, 4.60, 4.63, 4.63, 4.62, 4.68, 4.62};

	// Fill the histogram
	for (size_t i = 0; i < universe.size(); i++)
		g->SetPoint(i, universe.at(i), uncertainty.at(i) );
	
	
	g->SetTitle("Genie X Sec -- all;Universe;Percentage Uncertainty [%]");
	g->SetMarkerStyle(20);
	g->SetMarkerSize(1);
	g->SetMarkerColor(kBlue);
	g->SetLineWidth(0);
	g->Draw("");
	

	// c->Print(print_name);

}
// ------------------------------------------------------------------------------------------------------------
double GetMax(double m1sig, double p1sig){
	
	if (m1sig > p1sig) return m1sig;
	else return p1sig; 
	
}
// ------------------------------------------------------------------------------------------------------------
// Function to make the final uncertainty plot
void Make_uncertainty_plot(TCanvas* c, TH1D* h1, TH1D* h2, const char *print_name, std::pair<std::string, double> genie_all,
				std::pair<std::string, double> genie_indiv, std::pair<std::string, double> ccmec,
				std::pair<std::string, double> ccqe, std::pair<std::string, double> reinteractions_all,
				std::pair<std::string, double> reinteractions_indiv, int n_uni_model, int n_uni_reint, int n_uni_genie  ){

	c->cd();

	TPaveText *pt_gall = new TPaveText(.165,-0.4,.445,-0.0);
	pt_gall->AddText("genie all");
	pt_gall->SetBorderSize(0);
	pt_gall->SetFillColor(0);
	pt_gall->SetFillStyle(0);

	TPaveText *pt_gindiv = new TPaveText(0.5,-0.4,0.85,-0.0);
	pt_gindiv->AddText("genie indiv");
	pt_gindiv->SetBorderSize(0);
	pt_gindiv->SetFillColor(0);
	pt_gindiv->SetFillStyle(0);

	TPaveText *pt_ccmec = new TPaveText(1.09,-0.4,1.47,-0.0);
	pt_ccmec->AddText("ccmec");
	pt_ccmec->SetBorderSize(0);
	pt_ccmec->SetFillColor(0);
	pt_ccmec->SetFillStyle(0);

	TPaveText *pt_ccqe = new TPaveText(1.49,-0.4,1.89,-0.0);
	pt_ccqe->AddText("ccqe");
	pt_ccqe->SetBorderSize(0);
	pt_ccqe->SetFillColor(0);
	pt_ccqe->SetFillStyle(0);

	TPaveText *pt_reint_all = new TPaveText(2.11,-0.4,2.49,-0.0);
	pt_reint_all->AddText("reinteraction all");
	pt_reint_all->SetBorderSize(0);
	pt_reint_all->SetFillColor(0);
	pt_reint_all->SetFillStyle(0);

	TPaveText *pt_reint_indiv = new TPaveText(2.50,-0.4,2.90, -0.0);
	pt_reint_indiv->AddText("reinteraction indiv");
	pt_reint_indiv->SetBorderSize(0);
	pt_reint_indiv->SetFillColor(0);
	pt_reint_indiv->SetFillStyle(0);

	TLegend *legend = new TLegend(0.60, 0.65, 0.90, 0.85);
	legend->AddEntry(h1,  Form("Genie: %i Univ Multisim STD", n_uni_genie),"");
	legend->AddEntry(h1,    Form("CCMEC/CCQE: %i Univ Multisim STD", n_uni_model),"");
	legend->AddEntry(h1,    Form("Reinteraction: %i Univ Multisim STD", n_uni_reint),"");
	legend->SetTextFont(62);
	legend->SetBorderSize(0);
	legend->SetFillStyle(0);

	h1->GetXaxis()->SetLabelSize(0);
	h2->GetXaxis()->SetLabelSize(0);

	// gPad->SetBottomMargin(0.25);
	gStyle->SetHistMinimumZero();

	h1->SetFillColor(kGreen+1);
	h1->SetBarWidth(0.4);
	h1->SetBarOffset(0.1);
	h1->SetStats(0);
	h1->SetMaximum(8);
	h1->GetXaxis()->SetLabelOffset(0.03);
	
	h2->SetFillColor(kGreen-9);
	h2->SetBarWidth(0.4);
	h2->SetBarOffset(0.5);
	h2->SetStats(0);

	h1->SetBinContent(1, genie_all.second);
	
	// h1->GetXaxis()->SetBinLabel(1, genie_all.first.c_str());
	
	h2->SetBinContent(1, genie_indiv.second);
	// h2->GetXaxis()->SetBinLabel(1, genie_indiv.first.c_str());
	
	h1->SetBinContent(2, ccmec.second);
	// h1->GetXaxis()->SetBinLabel(2, ccmec.first.c_str());
	
	h2->SetBinContent(2, ccqe.second);
	// h2->GetXaxis()->SetBinLabel(2, ccqe.first.c_str());
	
	h1->SetBinContent(3 , reinteractions_all.second);
	// h1->GetXaxis()->SetBinLabel(3, reinteractions_all.first.c_str());
	
	h2->SetBinContent(3 ,reinteractions_indiv.second);
	// h2->GetXaxis()->SetBinLabel(3, reinteractions_indiv.first.c_str());
	

	h1->Draw("b");
	h2->Draw("b same");
	pt_gall->Draw();
	pt_gindiv->Draw();
	pt_ccmec->Draw();
	pt_ccqe->Draw();
	pt_reint_all->Draw();
	pt_reint_indiv->Draw();

	legend->Draw();

	c->Print(print_name);
}
// ------------------------------------------------------------------------------------------------------------
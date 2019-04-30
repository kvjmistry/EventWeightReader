
// Histogram and canvas definitions
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
void Print_Error(std::string CV_str, double CV, std::string Err_str, double Err){
    std::cout << CV_str << ": " << CV << " "<< Err_str << ": " << 100 * Err / CV << " %"<< std::endl;
}
// ------------------------------------------------------------------------------------------------------------
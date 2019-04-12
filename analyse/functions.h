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
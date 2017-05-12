void bmethods(string filename = "TTBarProcessorLeft.root")
{
	TCanvas * c1 = new TCanvas("c1", "The 3d view",0,0,500,500);
	TFile *_file0 = TFile::Open(filename.c_str());
	string cut = "MCPDG == 5 && MCMass > 180";
	string cuts = cut+" &&";
	int totaltag = Stats->Draw("B1momentum",cut.c_str());
	int singletag  = Stats->Draw("B1momentum",(cuts+"B1VtxTag != 0 && B2VtxTag == 0 && B1KaonTag == 0 && B2KaonTag ==0").c_str());
	    singletag += Stats->Draw("B1momentum",(cuts+"B1VtxTag == 0 && B2VtxTag != 0 && B1KaonTag == 0 && B2KaonTag ==0").c_str());
	    singletag += Stats->Draw("B1momentum",(cuts+"B1VtxTag == 0 && B2VtxTag == 0 && B1KaonTag != 0 && B2KaonTag ==0").c_str());
	    singletag += Stats->Draw("B1momentum",(cuts+"B1VtxTag == 0 && B2VtxTag == 0 && B1KaonTag == 0 && B2KaonTag !=0").c_str());
	int doubletag  = Stats->Draw("B1momentum",(cuts+"B1VtxTag != 0 && B2VtxTag != 0 && B1KaonTag == 0 && B2KaonTag ==0").c_str());
	    doubletag += Stats->Draw("B1momentum",(cuts+"B1VtxTag == 0 && B2VtxTag != 0 && B1KaonTag != 0 && B2KaonTag ==0").c_str());
	    doubletag += Stats->Draw("B1momentum",(cuts+"B1VtxTag == 0 && B2VtxTag == 0 && B1KaonTag != 0 && B2KaonTag !=0").c_str());
	    doubletag += Stats->Draw("B1momentum",(cuts+"B1VtxTag != 0 && B2VtxTag == 0 && B1KaonTag == 0 && B2KaonTag !=0").c_str());
	    doubletag += Stats->Draw("B1momentum",(cuts+"B1VtxTag != 0 && B2VtxTag == 0 && B1KaonTag != 0 && B2KaonTag ==0").c_str());
	    doubletag += Stats->Draw("B1momentum",(cuts+"B1VtxTag == 0 && B2VtxTag != 0 && B1KaonTag == 0 && B2KaonTag !=0").c_str());
	int tripletag  = Stats->Draw("B1momentum",(cuts+"B1VtxTag != 0 && B2VtxTag != 0 && B1KaonTag != 0 && B2KaonTag ==0").c_str());
	    tripletag += Stats->Draw("B1momentum",(cuts+"B1VtxTag != 0 && B2VtxTag != 0 && B1KaonTag == 0 && B2KaonTag !=0").c_str());
	    tripletag += Stats->Draw("B1momentum",(cuts+"B1VtxTag != 0 && B2VtxTag == 0 && B1KaonTag != 0 && B2KaonTag !=0").c_str());
	    tripletag += Stats->Draw("B1momentum",(cuts+"B1VtxTag == 0 && B2VtxTag != 0 && B1KaonTag != 0 && B2KaonTag !=0").c_str());
	int quadtag  = Stats->Draw("B1momentum",(cuts+"B1VtxTag != 0 && B2VtxTag != 0 && B1KaonTag != 0 && B2KaonTag !=0").c_str());
	int zerotag  = Stats->Draw("B1momentum",(cuts+"B1VtxTag == 0 && B2VtxTag == 0 && B1KaonTag == 0 && B2KaonTag ==0").c_str());
	std::cout << "Total events:\t" << totaltag << endl;
	std::cout << "Zero   tag:\t" <<   zerotag << "(" <<   (float)zerotag/totaltag*100 << "%)" << endl;
	std::cout << "Single tag:\t" << singletag << "(" << (float)singletag/totaltag*100 << "%)" << endl;
	std::cout << "Double tag:\t" << doubletag << "(" << (float)doubletag/totaltag*100 << "%)" << endl;
	std::cout << "Triple tag:\t" << tripletag << "(" << (float)tripletag/totaltag*100 << "%)" << endl;
	std::cout << "Quadro tag:\t" <<   quadtag << "(" <<   (float)quadtag/totaltag*100 << "%)" << endl;
	std::cout << "Check:     \t" << singletag+doubletag+tripletag+quadtag+zerotag << endl;
}


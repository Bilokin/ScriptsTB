#include <iomanip>
string BKG_PATH = "/exp/flc/bilokin/Training/bbbar-bkg/notag/";
string SIG_PATH = "/exp/flc/bilokin/Training/bbbar-250GeV/recoverytest/";

void beff(string filename = "TTBarProcessorLeftSignalExcludedFree.root")
{
	cout.precision(2);
	std::setw(4);
	std::setprecision(5);
	cout.setf(ios::fixed);
	//TCanvas * c1 = new TCanvas("c1", "The 3d view",0,0,500,500);
	TCanvas*c1 = new TCanvas("c1", "Data-MC",0,0,500,500);
	float signalLumiLeft = 13.53;
	int nBkgLeft = 8;
	string genCuts = "MCPDG == 5 && MCMass > 180";
	string zreturnCuts = "InvMass > 180 && maxPhotonEnergy < 40 ";
	string bCuts = "B1btag > 0.8 && B2btag > 0.3";
	string mCuts = "B1mass + B2mass < 120";
	string sCuts = "Sphericity < 0.2 ";
	//string bkgFileLeft[5] = {"ZZhadronicLeft.root","ZZsemileptonicLeft.root","WWhadronicLeft.root", "WWsemileptonicLeft.root", "HZhadronicLeft.root"};
	string bkgFileLeft[8] = {"ZZhadronicLeft.root","WWhadronicLeft.root", "HZhadronicLeft.root","ZZsemileptonicLeft.root","ZReturnLeft3.root", "WWsemileptonicLeft.root","ZZWWhadronicLeft.root", "CCLeft.root        "};
	cout << "File \t\t\t Initial \t Btag \t\t\t InvMass \t Mass \t\t Sphericity \n";
	for (unsigned int i = 0; i < nBkgLeft; i++) 
	{

		float initialn = getEvents(bkgFileLeft[i]);
		float btagn = getEvents(bkgFileLeft[i], bCuts.c_str());
		float zreturnn = getEvents(bkgFileLeft[i], (bCuts + "&&" +zreturnCuts).c_str());
		float mn = getEvents(bkgFileLeft[i], (bCuts + "&&" +zreturnCuts + "&&" +mCuts).c_str());
		float sn = getEvents(bkgFileLeft[i], (bCuts + "&&" +zreturnCuts + "&&" +mCuts+ "&&" +sCuts).c_str());
		cout << "" << bkgFileLeft[i] << "\t" << (int)initialn;
		cout << "\t\t" << (int)btagn << setw(2) << " (" << btagn/initialn*100 << "%)" << setw(6) ;
		cout << "\t" << (int)zreturnn << " (" << zreturnn/initialn*100 << "%)";
		cout << "\t" << (int)mn << " (" << mn/initialn*100 << "%)" ;
		cout << "\t" << (int)sn << " (" << sn/initialn*100 << "%)" << endl;

	}
	cout << "----------------------------------------------------------------------------------\n";
	float initialn = getEvents(filename,"",1);
	float btagn = getEvents(filename, bCuts.c_str(),1);
	float zreturnn = getEvents(filename, (bCuts + "&&" +zreturnCuts).c_str(),1);
	float mn = getEvents(filename, (bCuts + "&&" +zreturnCuts + "&&" +mCuts).c_str(),1);
	float sn = getEvents(filename, (bCuts + "&&" +zreturnCuts + "&&" +mCuts+ "&&" +sCuts).c_str(),1);
	cout << "" << "Signal            " << "\t" << (int)initialn;
	cout << "\t\t" << (int)btagn << setw(2) << " (" << btagn/initialn*100 << "%)" << setw(6) ;
	cout << "\t" << (int)zreturnn << " (" << zreturnn/initialn*100 << "%)";
	cout << "\t" << (int)mn << " (" << mn/initialn*100 << "%)" ;
	cout << "\t" << (int)sn << " (" << sn/initialn*100 << "%)" << endl;
	
}
float getEvents(string filename, string cuts = "", bool signal = false)
{
	string fullPath = BKG_PATH + filename;
	if (signal) 
	{
		fullPath = SIG_PATH + filename;
	}
	TFile * file = TFile::Open(fullPath.c_str());
	file->cd();
	int nReco = Stats->Draw("MCMass",cuts.c_str(),"same");
	file->Close();
	return nReco;
}

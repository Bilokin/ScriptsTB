
string BKG_PATH = "/exp/flc/bilokin/Training/bbbar-bkg/notag/";

void bcrossection(string filename = "TTBarProcessorLeft.root")
{
	//TCanvas * c1 = new TCanvas("c1", "The 3d view",0,0,500,500);
	//float signalLumiLeft = 13.53;
	float signalLumiLeft = 250.;
	const int nBkgLeft = 5;
	string bkgFileLeft[5] = {"ZZhadronicLeft.root","ZZsemileptonicLeft.root","WWhadronicLeft.root", "WWsemileptonicLeft.root", "HZhadronicLeft.root"};
	//float bkgLumiLeft[5] = {250., 250., 72.3*270.5/1074.5, 102.2*156./1919., 1000.};
	float bkgLumiLeft[5] = {250., 250., 72.3, 102.2, 1000.};

	float signalLumiRight = 250.0;
	//float signalLumiRight = 20.0;
	const int nBkgRight = 3;
	string bkgFileRight[3] = {"ZZhadronicRight.root","ZZsemileptonicRight.root", "HZhadronicRight.root"};
	float bkgLumiRight[3] = {250., 250., 1000.};

	TFile * file = TFile::Open(filename.c_str());
	string genCuts = "MCPDG == 5 && MCMass > 180";
	//string cuts = "InvMass > 180 && maxPhotonEnergy < 40";
	//string cuts = "B1btag > 0.85 && B2btag > 0.85 && InvMass > 180 && maxPhotonEnergy < 40 && B1mass + B2mass < 140 && Sphericity < 0.15 ";
	string cuts = "B1btag > 0.85 && B2btag > 0.85 && InvMass > 200 && maxPhotonEnergy < 40 && B1mass + B2mass < 100 && Sphericity < 0.1 ";
	//string cuts = "InvMass > 180 && bbbarAngle > 2.6 && maxPhotonEnergy < 40 && B1mass + B2mass < 150 && Sphericity < 0.35 && B1Y < 1 && B2mass < 70";
	//string cuts = "B2mass > 5 && InvMass >30+150* abs(B2costheta) && InvMass > 150 &&  maxPhotonEnergy < 40 && B2chargeBalance < 0.98 && B2chargeBalance > 0.05 && bbbarAngle > 2.6 && B1mass + B2mass < 150 && Sphericity < 0.3 && B1Y < 0.8";
	//string cuts = " InvMass > 150 &&  maxPhotonEnergy < 40 && B2chargeBalance < 0.98 && B2chargeBalance > 0.05 && bbbarAngle > 2.6 && InvMass >30+150* abs(B2costheta) ";
	//string cuts = "InvMass > 150 &&  maxPhotonEnergy < 40 && B2chargeBalance < 0.98 && B2chargeBalance > 0.05 && bbbarAngle > 2.6 && InvMass >30+150* abs(B2costheta) ";
	string cuts_true = genCuts + " && " + cuts;
	
	int nMCB = GenTree->Draw("qMCcostheta[0]",genCuts.c_str());
	float nRecoTrue = Stats->Draw("B1costheta",cuts_true.c_str(),"same");
	float nRecoTrueAfterBcuts = Stats->Draw("B1costheta", genCuts.c_str(),"same");
	float nReco = Stats->Draw("B1costheta", cuts.c_str(),"same");
	
	cout << "Gen: " << nMCB << endl;
	cout << "After btag: " << nRecoTrueAfterBcuts << " (" << nRecoTrueAfterBcuts/nMCB*100 << "%)" << endl;
	//cout << "Reco efficiency: " << nReco << " (" << nReco/nMCB *100 << "%)" << endl;
	int nBkgReco =  0;//getBkgNevents("ZZhadronicLeft.root", signalLumi / 250.);
	for (unsigned int i = 0; i < nBkgRight; i++) 
	{
		nBkgReco += getBkgNevents(bkgFileRight[i], cuts,signalLumiRight / bkgLumiRight[i]);
	}//*/
	/*for (unsigned int i = 0; i < nBkgLeft; i++) 
	{
		nBkgReco += getBkgNevents(bkgFileLeft[i], cuts,signalLumiLeft / bkgLumiLeft[i]);
	}*/
	cout << "Bkg: " << nBkgReco << " (" << nBkgReco/(nReco + nBkgReco) *100 << "%)" << endl;
	
	cout << "Reco purity: " << nReco << " (" << nRecoTrue/(nReco+nBkgReco) *100 << "%)" << endl;
	cout << "Reco_true efficiency: " << nRecoTrue << " (" << nRecoTrue/nMCB*100 << "%)" << endl;

	//file->Close();
}
int getBkgNevents(string filename = "ZZhadronicLeft.root",string cuts = "InvMass > 180 && maxPhotonEnergy < 40", float ratio = 1.)
{
	string fullPath = BKG_PATH + filename;
	TFile * file = TFile::Open(fullPath.c_str());
	file->cd();
	int nReco = Stats->Draw("B1costheta",cuts.c_str(),"same");
	Summary->Draw("nEvents >> nEvt");
	nReco *= ratio;
	cout << "\tFile:" << filename << "\tnEvents:\t" << nReco << " (" << nEvt->GetMean() << ")" << endl;
	file->Close();
	return nReco;
}


string BKG_PATH = "/exp/flc/bilokin/Training/bbbar-bkg/exclude/";
string SIGNAL_PATH = "/exp/flc/bilokin/Training/bbbar-250GeV/recoverytest/";
void bsigma(string filenamel = "TTBarProcessorLeftTight.root", string filenamer = "TTBarProcessorRightTight.root")
{
	float signalLumiLeft = 13.53;
	float signalLumiRight = 20.0;
	float ZZLumiLeft = 250.09;
	float ZZLumiRight = 250.14;
	float HZLumiLeft = 1000.0;
	float HZLumiRight = 1000.0;
	float ZZWWLumiLeft = 86.77;
	float ZZWWLumiRight = 251.56;
	TFile * filel = TFile::Open((SIGNAL_PATH+filenamel).c_str());
	string genCutsCM = "MCPDG == 5 && MCMass > 200";
	string genCutsCMC = "MCPDG == 4 && MCMass > 200";
	string genCutsZ = "MCPDG == 5 && MCMass < 120";
	string genCutsZC = "MCPDG == 4 && MCMass < 120";
	float bbbarcml = GenTree->Draw("MCMass",genCutsCM.c_str());
	float ccbarcml = GenTree->Draw("MCMass",genCutsCMC.c_str());
	float bbbarzl = GenTree->Draw("MCMass",genCutsZ.c_str());
	cout << bbbarcml << endl;
	TFile * filer = TFile::Open((SIGNAL_PATH+filenamer).c_str());
	filer->cd();
	float bbbarcmr = GenTree->Draw("MCMass",genCutsCM.c_str());
	float ccbarcmr = GenTree->Draw("MCMass",genCutsCMC.c_str());
	float bbbarzr = GenTree->Draw("MCMass",genCutsZ.c_str());
	cout << bbbarcmr << endl;
	cout << "Process |\tUnpol\teLpR\teRpL\n";
	Print("bbbar  ", bbbarcml / signalLumiLeft, bbbarcmr / signalLumiRight);
	Print("ccbar  ", ccbarcml / signalLumiLeft, ccbarcmr / signalLumiRight);
	Print("g bbbar", bbbarzl / signalLumiLeft, bbbarzr / signalLumiRight);
	TFile * filezzl = TFile::Open((BKG_PATH+"ZZhadronicLeft.root").c_str());
	filezzl->cd();
	float zzcml = GenTree->Draw("MCMass","");
	TFile * filezzsr = TFile::Open((BKG_PATH+"ZZhadronicRight.root").c_str());
	filezzsr->cd();
	float zzcmr = GenTree->Draw("MCMass","");
	Print("ZZ had ", zzcml / ZZLumiLeft, zzcmr / ZZLumiRight);
	TFile * filezzsl = TFile::Open((BKG_PATH+"ZZsemileptonicLeft.root").c_str());
	filezzsl->cd();
	float zzcmsl = GenTree->Draw("MCMass","");
	TFile * filezzsr = TFile::Open((BKG_PATH+"ZZsemileptonicRight.root").c_str());
	filezzsr->cd();
	float zzcmsr = GenTree->Draw("MCMass","");
	Print("ZZ sl    ", zzcmsl / ZZLumiLeft, zzcmsr / ZZLumiRight);
	TFile * filehzl = TFile::Open((BKG_PATH+"HZhadronicLeft.root").c_str());
	filezzl->cd();
	float hzcml = GenTree->Draw("MCMass","");
	TFile * filehzsr = TFile::Open((BKG_PATH+"HZhadronicRight.root").c_str());
	filehzsr->cd();
	float hzcmr = GenTree->Draw("MCMass","");
	Print("HZ had ", hzcml / HZLumiLeft, hzcmr / HZLumiRight);
	TFile * filezwl = TFile::Open((BKG_PATH+"ZZWWhadronicLeft.root").c_str());
	filezwl->cd();
	float zwcml = GenTree->Draw("MCMass","");
	float zwcmr = 56562;
	Print("ZZWW had", zwcml / ZZWWLumiLeft, zwcmr / ZZWWLumiRight);
}
void Print(string name, float left, float right)
{
	cout << name << ":\t" <<  (left+right)/4. << "\t" << left << "\t" << right << endl;
	
}

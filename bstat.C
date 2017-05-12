void bstat(string filename = "TTBarProcessorLeft.root")
{
	TFile * file = TFile::Open(filename.c_str());
	TH1F * hist = new TH1F("tmp","",10000,0,10000);
	int totalNumber = 0;
	int generatedUsedNumber = 0;
	int afterBtagNumber = 0;
	int afterNlNumber = 0;
	int chargedBNumber = 0;
	int KNumber = 0;
	int afterKinematicsNumber = 0;
	int afterMassNumber = 0;
	int afterGenMassNumber = 0;
	Summary->SetBranchAddress("nEvents",&totalNumber);
	Summary->SetBranchAddress("nGenUsed",&generatedUsedNumber);
	Summary->SetBranchAddress("nAfterBtagCuts",&afterBtagNumber);
	Summary->SetBranchAddress("nAfterGenMassCuts",&afterGenMassNumber);
	Summary->SetBranchAddress("nAfterLeptonCuts",&afterNlNumber);
	Summary->SetBranchAddress("nChargedB",&chargedBNumber);
	Summary->SetBranchAddress("nKaons",&KNumber);
	Summary->SetBranchAddress("nAfterKinematicCuts",&afterKinematicsNumber);
	Summary->SetBranchAddress("nAfterMassCuts",&afterMassNumber);
	Summary->GetEvent(0);
	//int ccbarbkg = 
	cout << "---------------------\n";
	cout << "---------------------\n";
	printStat("Total gen events", totalNumber, totalNumber);
	printStat("Total bbar events", generatedUsedNumber, totalNumber);
	totalNumber = generatedUsedNumber;
	printStat("After mass cuts", afterGenMassNumber, totalNumber);
	cout << "---------------------\n";
	//printStat("After mass cuts", afterMassNumber, totalNumber);
	printStat("After btag cuts", afterBtagNumber, totalNumber);
	printStat("Number of Kaons", KNumber, totalNumber);
	printStat("Number of charged B", chargedBNumber, totalNumber);
	printStat("After kinematics", afterKinematicsNumber, totalNumber);
	//printStat("Normalized to leptons", afterMassNumber*0.60 * totalNumber / afterBtagNumber, totalNumber);

	cout << "---------------------\n";
	cout << "---------------------\n";
	file->Close();
}
void printStat(string name, int number, int totalnumber)
{
	cout << name << ":\t\t" << number << "\t" << (float)number / totalnumber  * 100 << "%\n";
}

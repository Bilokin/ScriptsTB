void comparator(string observable, string cuts)
{
	string signalPath = "/exp/flc/bilokin/Training/bbbar-250GeV/recoverytest/TTBarProcessorLeftTightTest.root";
	string bkgPath = "/exp/flc/bilokin/Training/bbbar-bkg/ZZhadronicLeftTest.root";
	//string bkgPath = "/exp/flc/bilokin/Training/bbbar-bkg/loose/WWsemileptonicLeft.root";
	string signalCuts = "InvMass > 180 && maxPhotonEnergy < 40 ";
	string finalCuts =signalCuts + " &&" + cuts;
	string bkgToDraw = observable + ">> bkg";
	string signalToDraw = observable + ">> signal";
	TFile * fileSignal = TFile::Open(signalPath.c_str());
	fileSignal->cd();
	float nSignalCuts =  Stats->Draw(signalToDraw.c_str(), signalCuts.c_str());
	float nSignalAfterNewCuts = Stats->Draw(signalToDraw.c_str(), finalCuts.c_str());
	signal->GetXaxis()->SetTitle(observable.c_str());
	makePretty(signal, kBlue);
	TFile * fileBkg = TFile::Open(bkgPath.c_str());
	fileBkg->cd();
	int nBkgAfterNewCuts = Stats->Draw(bkgToDraw.c_str(), finalCuts.c_str(), "same");
	makePretty(bkg,kRed);
	cout << "Purity: \t" << nSignalAfterNewCuts / (nBkgAfterNewCuts + nSignalAfterNewCuts) * 100 << "%\n";
	cout << "Efficiency:\t" << nSignalAfterNewCuts / (nSignalCuts) * 100 << "%\n";
}
void makePretty(TH1 * vtxTotal, int color, int line = 0)
{
	vtxTotal->SetLineWidth(3);
	vtxTotal->SetLineStyle(line);
	vtxTotal->SetLineColor(color);
	vtxTotal->SetMinimum(0);
	vtxTotal->SetStats(0);
}


string BKG_PATH = "/exp/flc/bilokin/Training/bbbar-bkg/notag/";
void btag(string filename = "TTBarProcessorLeft.root")
{
	const int npoints = 20;
	float purity[npoints];
	float efficiency[npoints];
	float btags[npoints];
	bool left = 1;
	TCanvas *c2 = new TCanvas("c2","Exponential fit for number of pads in clusters",200,10,500,500);

	float signalLumiLeft = 13.53;
	const int nBkgLeft = 5;
	string bkgFileLeft[5] = {"ZZhadronicLeft.root","ZZsemileptonicLeft.root","WWhadronicLeft.root", "WWsemileptonicLeft.root", "HZhadronicLeft.root"};
	//float bkgLumiLeft[5] = {250., 250., 72.3*270.5/1074.5, 102.2*156./1919., 1000.};
	float bkgLumiLeft[5] = {250., 250., 72.3, 102.2, 1000.};

	float signalLumiRight = 20.0;
	const int nBkgRight = 3;
	string bkgFileRight[3] = {"ZZhadronicRight.root","ZZsemileptonicRight.root", "HZhadronicRight.root"};
	float bkgLumiRight[3] = {250., 250., 1000.};

	TFile * file = TFile::Open(filename.c_str());
	string genCuts = "MCPDG == 5 && MCMass > 180";
	//string cuts = "InvMass > 0";
	//string cuts = "InvMass > 180 && maxPhotonEnergy < 40";
	string cuts = "InvMass > 180 && maxPhotonEnergy < 40 && B1mass + B2mass < 140 && Sphericity < 0.15";
	//string cuts = "InvMass > 200 && maxPhotonEnergy < 40 && B1mass + B2mass < 100 && Sphericity < 0.1 ";
	//string cuts = "InvMass > 180 && bbbarAngle > 2.6 && maxPhotonEnergy < 40 && B1mass + B2mass < 150 && Sphericity < 0.35 && B1Y < 1 && B2mass < 70";
	//string cuts = "B2mass > 5 && InvMass >30+150* abs(B2costheta) && InvMass > 150 &&  maxPhotonEnergy < 40 && B2chargeBalance < 0.98 && B2chargeBalance > 0.05 && bbbarAngle > 2.6 && B1mass + B2mass < 150 && Sphericity < 0.3 && B1Y < 0.8";
	//string cuts = " InvMass > 150 &&  maxPhotonEnergy < 40 && B2chargeBalance < 0.98 && B2chargeBalance > 0.05 && bbbarAngle > 2.6 && InvMass >30+150* abs(B2costheta) ";
	//string cuts = "InvMass > 150 &&  maxPhotonEnergy < 40 && B2chargeBalance < 0.98 && B2chargeBalance > 0.05 && bbbarAngle > 2.6 && InvMass >30+150* abs(B2costheta) ";
	//string cuts_true = genCuts + " && " + cuts;
	string btag1 = "&& B1btag > ";
	string btag2 = "&& B2btag > ";
	float btagcut = 0.0;
	int nMCB = GenTree->Draw("qMCcostheta[0]",genCuts.c_str());
	//int nMCB = Stats->Draw("qMCcostheta[0]",(genCuts+"&&"+cuts).c_str());
	file->Close();
	for (unsigned int i = 0; i < npoints; i++) 
	{
		float nBkg = 0.;
		string cutTrue = cuts + btag2 + intToStr(btagcut)+ btag1 + intToStr(btagcut) + "&&" + genCuts ;
		string cutReco = cuts + btag2 + intToStr(btagcut)+ btag1 + intToStr(btagcut);
		TFile * fileS = TFile::Open(filename.c_str());
		fileS->cd();
		float nRecoTrueAfterBcuts = getSignalNevents(filename, cutTrue);//Stats->Draw("B1costheta", cutTrue.c_str(),"same");
		float nReco = getSignalNevents(filename, cutReco);//Stats->Draw("B1costheta", cutReco.c_str(),"same");
		if (left) 
		{
			for (unsigned int j = 0; j < nBkgLeft; j++) 
			{
				nBkg+= getBkgNevents(bkgFileLeft[j], cutReco, signalLumiLeft / bkgLumiLeft[j]);
			}
		}
		else 
		{
			for (unsigned int j = 0; j < nBkgRight; j++) 
			{
				nBkg+= getBkgNevents(bkgFileRight[j], cutReco, signalLumiRight / bkgLumiRight[j]);
			}
		}
		fileS->Close();
		cout <<"Btag " << btagcut << ": True events " << nRecoTrueAfterBcuts << " of " << nReco << "+" << nBkg << endl;
		purity[i] = nRecoTrueAfterBcuts / (nReco + nBkg);
		efficiency[i] = nRecoTrueAfterBcuts / nMCB;
		cout << "P = " << purity[i] << "; E = " << efficiency[i]  << endl;
		btags[i] = btagcut;
		btagcut += 1./npoints;
	}
	TH1F* hist = new TH1F("radiusvsenergy",";b-tag cut value;Purity, Efficiency",npoints,0,1);
	hist->SetMaximum(1);
	hist->SetMinimum(0);
	hist->SetStats(0);

	hist->Draw();
	TGraph * pe =    getGraph(npoints,purity,efficiency, kRed);
	TGraph * btagp = getGraph(npoints,btags, purity, kBlue);
	TGraph * btage = getGraph(npoints,btags, efficiency, kGreen);
	btagp->Draw("lp");
	btage->Draw("lpsame");
	//pe->Draw("lpsame");
	TLegend *legendMean = new TLegend(0.24,0.18,0.57,0.36,NULL,"brNDC");
	legendMean->SetFillColor(kWhite);
	legendMean->SetBorderSize(0);
	legendMean->AddEntry(btagp,"btag - purity","pl");
	legendMean->AddEntry(btage,"btag - efficiency","pl");
	//legendMean->AddEntry(pe,"purity - efficiency","lp");
	legendMean->Draw();//*/
}
TGraph * getGraph(int npoints, float * x, float *y, int color = kBlue)
{
	TGraph * gr = new TGraph(npoints, x, y);
	//TGraph * gr = new TGraph(npoints,btags,efficiency);
	gr->SetFillStyle(3001);
	gr->SetLineColor(color);
	gr->SetLineWidth(1);
	gr->SetMarkerColor(color);
	gr->SetMarkerStyle(20);
	return gr;
}
int getSignalNevents(string filename = "ZZhadronicLeft.root",string cuts = "InvMass > 180 && maxPhotonEnergy < 40")
{
	TFile * fileS = TFile::Open(filename.c_str());
	fileS->cd();
	float nRecoTrueAfterBcuts = Stats->Draw("B1costheta", cuts.c_str(),"same");
	fileS->Close();
	return nRecoTrueAfterBcuts;
}
int getBkgNevents(string filename = "ZZhadronicLeft.root",string cuts = "InvMass > 180 && maxPhotonEnergy < 40", float ratio = 1.)
{
	string fullPath = BKG_PATH + filename;
	TFile * fileB = TFile::Open(fullPath.c_str());
	fileB->cd();
	int nReco = Stats->Draw("B1costheta",cuts.c_str(),"same");
	Summary->Draw("nEvents >> nEvt");
	nReco *= ratio;
	cout << "\tFile:" << filename << "\tnEvents:\t" << nReco << " (" << nEvt->GetMean() << ")" << endl;
	fileB->Close();
	return nReco;
}
string intToStr(float a)
{
	stringstream ss;
	ss << a;
	return str = ss.str();
}

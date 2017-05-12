
#include <unistd.h>
#include <iostream>
#include "/users/flc/bilokin/Processors/Macros/top/qp.C"
#define MAXV 8
string BKG_PATH = "/exp/flc/bilokin/Training/bbbar-bkg/factor/";

void b2asymmetry(string filename = "TTBarProcessorLeft.root", TCanvas * c1 = NULL, bool drawCorrection = true)
{
	TFile * file = TFile::Open(filename.c_str());
	float signalLumi = 13.53;
	float signalLumiR = 20.;
	int bin_e = 30;
	int max_e = 1;
	if (!c1) 
	{
		c1 = new TCanvas("c1", "Data-MC",0,0,500,500);
	}
	bool helicity = (filename == "TTBarProcessorLeft.root" || filename == "TTBarProcessorLeftTight.root" || filename == "TTBarProcessorLeftSmall.root" || filename == "TTBarProcessorLeftTight-3-4-5-6.root" || filename == "TTBarProcessorLeftFree.root" || filename == "TTBarProcessorLeftTight-3-4new.root") ;

	//c1->Divide(2,1);
	TH1F * cosReco = new TH1F("cosReco", "E(Ntracks)", bin_e,-1.0,max_e);
	//TH1F * correctedP = new TH1F("correctedP", "E(Ntracks)", bin_e,-1.0,max_e);
	TH1F * cosCCbar = new TH1F("cosCCbar", "E(Ntracks)", bin_e,-1.0,max_e);
	TH1F * cosRecoBkgMass = new TH1F("cosRecoBkgMass", "E(Ntracks)", bin_e,-1.0,max_e);
	TH1F * cosRecoBkgZZ = new TH1F("cosRecoBkgZZ", "E(Ntracks)", bin_e,-1.0,max_e);
	cosReco->Sumw2();
	cosRecoBkgMass->Sumw2();
	TH1F * cosGen = new TH1F("cosGen", ";cos#theta_{b}", bin_e,-1.0,max_e);

	cosGen->Sumw2();
	//correctedP->Draw();
	TTree * normaltree = Stats;
	//cosReco->SetLineColor(kBlue);
	cosReco->SetLineWidth(3);
	makePretty(cosGen, kGreen+1, 2);
	makePretty(cosReco, kBlack);
	makePretty(cosCCbar, kOrange, 2);
	makePretty(cosRecoBkgMass, kGray, 2);
	int forward = GenTree->Draw("qMCcostheta >> cosGen","qMCcostheta > 0 && MCMass > 200 && MCPDG == 5");
	int backward = GenTree->Draw("qMCcostheta >> +cosGen","qMCcostheta < 0 && MCMass > 200&& MCPDG == 5" );
	//int forward = Stats->Draw("qMCcostheta >> cosGen","qMCcostheta > 0 && InvMass > 200 && maxPhotonEnergy < 40 && MCPDG == 5");
	//int backward = Stats->Draw("qMCcostheta >> +cosGen","qMCcostheta < 0 && InvMass > 200 && maxPhotonEnergy < 40 && MCPDG == 5" );
	
//	string cuts = "&& MCMass > 200 && methodUsed";
	//string cuts = "&& InvMass > 180 && maxPhotonEnergy < 40 && B1mass < 60 + 20*B1costheta &&  B2mass < 60 + 20*B2costheta && methodUsed";
	//string cuts = "&& InvMass > 180 && maxPhotonEnergy < 40 && B1mass + B2mass < 60 && methodUsed";
	string cuts1 = "&&InvMass > 180 && maxPhotonEnergy < 40 && B1mass + B2mass < 130 && methodRefused == 0  && methodUsed == 1";
	string cuts2 = "&&InvMass > 180 && maxPhotonEnergy < 40 && B1mass + B2mass < 130 && methodRefused == 0  && methodUsed > 1";
	vector<string> cuts;
	cuts.push_back(cuts1);
	cuts.push_back(cuts2);
	TH1F * correctedP =  new TH1F("cosRecoBkgZZ", "E(Ntracks)", bin_e,-1.0,max_e);
	//string cuts = "&& InvMass > 180 && maxPhotonEnergy < 40 && B1mass + B2mass < 150 && Sphericity < 0.25 && methodUsed > 1"; //eRpL
	{
		int i = 0;
		file->cd();
		TH1F * cosReco1 = new TH1F("cosReco1", "E(Ntracks)", bin_e,-1.0,max_e);
		Stats->Project("cosReco1","qCostheta1",("qCostheta > -1.0"+cuts[i]).c_str());
		Stats->Project("+cosReco1","qCostheta2",("qCostheta > -1.0"+ cuts[i]).c_str());
		cout <<"HERE!\n";
		//int recoforward = normaltree->Draw("-Top1costheta*UsedBTVCM >> +cosReco", "(qCostheta > 0 && W1mass > 0 && UsedBTVCM !=0 && (MCBOscillation > -1 && MCBBarOscillation > -1 ))"); //UsedBTVCM
		//int recobackward = normaltree->Draw("-Top1costheta*UsedBTVCM >> +cosReco", "(qCostheta < 0 && qCostheta > -1.0 && W1mass > 0 && UsedBTVCM !=0 && (MCBOscillation > -1 && MCBBarOscillation > -1 ))");
		string zcuts = "MCMass < 180 " + cuts[i];
		string cccuts = "MCPDG == 4" + cuts[i];
		normaltree->Draw("qCostheta>> cosRecoBkgMass", zcuts.c_str());
		normaltree->Draw("qCostheta>> cosCCbar", cccuts.c_str());
		TH1F * cosRecoBkgZZ1 = new TH1F("cosRecoBkgZZ1", "E(Ntracks)", bin_e,-1.0,max_e);
		if (helicity) 
		{
			addBkg(cosRecoBkgZZ1,"ZZhadronicLeft",cuts[i],bin_e, signalLumi / 250.); 
			addBkg(cosRecoBkgZZ1,"ZZsemileptonicLeft",cuts[i],bin_e, signalLumi / 250.);
			addBkg(cosRecoBkgZZ1,"WWsemileptonicLeft",cuts[i],bin_e, signalLumi / 100.);
			addBkg(cosRecoBkgZZ1,"WWhadronicLeft",cuts[i],bin_e, signalLumi / (72.));
			addBkg(cosRecoBkgZZ1,"HZhadronicLeft",cuts[i],bin_e, signalLumi / 1000.);
		}
		else 
		{
			addBkg(cosRecoBkgZZ1,"ZZhadronicRight",cuts[i],bin_e, signalLumiR / 250.);
			addBkg(cosRecoBkgZZ1,"ZZsemileptonicRight",cuts[i],bin_e, signalLumiR / 250.);
			addBkg(cosRecoBkgZZ1,"HZhadronicRight",cuts[i],bin_e, signalLumiR / 1000.);
		}
		cosReco->Add(cosRecoBkgZZ1);//*/
		cosRecoBkgZZ->Add(cosRecoBkgZZ1);
		file->cd();
		if (drawCorrection) 
		{
			float pglobal = getPGlobal(bin_e);
			correctedP->Add(computeCorrectionF(cosReco1, pglobal, bin_e));
		}
		cosReco->Add(cosReco1);
		cout <<"END!\n";
	}
	/*{
		int i = 1;
		file->cd();
		TH1F * cosReco2 = new TH1F("cosReco2", "E(Ntracks)", bin_e,-1.0,max_e);
		Stats->Project("cosReco2","qCostheta1",("qCostheta > -1.0"+cuts[i]).c_str());
		Stats->Project("+cosReco2","qCostheta2",("qCostheta > -1.0"+ cuts[i]).c_str());
		cout <<"HERE!\n";
		//int recoforward = normaltree->Draw("-Top1costheta*UsedBTVCM >> +cosReco", "(qCostheta > 0 && W1mass > 0 && UsedBTVCM !=0 && (MCBOscillation > -1 && MCBBarOscillation > -1 ))"); //UsedBTVCM
		//int recobackward = normaltree->Draw("-Top1costheta*UsedBTVCM >> +cosReco", "(qCostheta < 0 && qCostheta > -1.0 && W1mass > 0 && UsedBTVCM !=0 && (MCBOscillation > -1 && MCBBarOscillation > -1 ))");
		string zcuts = "MCMass < 180 " + cuts[i];
		string cccuts = "MCPDG == 4" + cuts[i];
		normaltree->Draw("qCostheta>> +cosRecoBkgMass", zcuts.c_str());
		normaltree->Draw("qCostheta>> +cosCCbar", cccuts.c_str());
		TH1F * cosRecoBkgZZ2 = new TH1F("cosRecoBkgZZ2", "E(Ntracks)", bin_e,-1.0,max_e);
		if (helicity) 
		{
			addBkg(cosRecoBkgZZ2,"ZZhadronicLeft",cuts[i],bin_e, signalLumi / 250.); 
			addBkg(cosRecoBkgZZ2,"ZZsemileptonicLeft",cuts[i],bin_e, signalLumi / 250.);
			addBkg(cosRecoBkgZZ2,"WWsemileptonicLeft",cuts[i],bin_e, signalLumi / 100.);
			addBkg(cosRecoBkgZZ2,"WWhadronicLeft",cuts[i],bin_e, signalLumi / (72.));
			addBkg(cosRecoBkgZZ2,"HZhadronicLeft",cuts[i],bin_e, signalLumi / 1000.);
		}
		else 
		{
			addBkg(cosRecoBkgZZ2,"ZZhadronicRight",cuts[i],bin_e, signalLumiR / 250.);
			addBkg(cosRecoBkgZZ2,"ZZsemileptonicRight",cuts[i],bin_e, signalLumiR / 250.);
			addBkg(cosRecoBkgZZ2,"HZhadronicRight",cuts[i],bin_e, signalLumiR / 1000.);
		}
		cosReco2->Add(cosRecoBkgZZ);///
		cosRecoBkgZZ->Add(cosRecoBkgZZ2);
		file->cd();
		if (drawCorrection) 
		{
			float pglobal = getPGlobal(bin_e, 3);
			correctedP->Add(computeCorrectionF(cosReco2, pglobal, bin_e));
		}
		cosReco->Add(cosReco2);
	}*/
	makePretty(cosRecoBkgZZ,kRed,2);
	makePretty(correctedP, kBlue);
	float fitRange = 0.95;
	cosGen->SetStats(0);
	if (filename == "TTBarProcessorRight.root") 
	{
		cosGen->SetTitle("e_{R}p_{L}");
	}
	if (filename == "TTBarProcessorLeft.root") 
	{
		cosGen->SetTitle("e_{L}p_{R}");
	}
	TF1 * fgen = new TF1("fgen","pol2",-1,1);
	TF1 * freco = new TF1("freco","pol2",-fitRange, fitRange);
	fgen->SetLineColor(kGreen);
	fgen->SetLineStyle(3);
	freco->SetLineStyle(3);
		//cosGen->SetTitle("e_{L}p_{R} minivector original");
	int mid = (float)bin_e / 2;
	//cosGen->Scale(1./ cosGen->GetEntries());
	cosGen->Scale(1./ (cosGen->GetBinContent(mid) + cosGen->GetBinContent(mid+1)) * (cosReco->GetBinContent(mid) + cosReco->GetBinContent(mid+1)));
	//cosReco->Scale(1./ cosReco->GetEntries());
	cosGen->Fit("fgen","Q");
	cosReco->Fit("freco", "QR");
	cosGen->SetMinimum(0);
	cosGen->Draw("he");
	fgen->Draw("same");
	//cosGen->Draw();
	cosGen->SetMinimum(0);
	cosReco->Draw("samee");
	if (drawCorrection) 
	{
		correctedP->Draw("samee");
		correctedP->Fit("freco", "QR");
	}
	cosRecoBkgMass->Draw("same");
	cosRecoBkgZZ->Draw("same");
	cosCCbar->Draw("same");
	TLegend *legendMean2 = new TLegend(0.15,0.6,0.5,0.85,NULL,"brNDC");
        legendMean2->SetFillColor(kWhite);
        legendMean2->SetBorderSize(0);
        legendMean2->AddEntry(cosGen,"Generated","f");
        legendMean2->AddEntry(cosReco,"Reconstructed","f");
	if (drawCorrection) 
	{
		legendMean2->AddEntry(correctedP,"Corrected","f");
	}
        legendMean2->AddEntry(cosRecoBkgMass,"Z return background","f");
        legendMean2->AddEntry(cosCCbar,"c#bar{c} background","f");
        legendMean2->AddEntry(cosRecoBkgZZ,"ZZ ZH WW background","f");
	legendMean2->Draw();
	float afbgen = (float)(forward - backward) / (float) (forward + backward);
	//float afbreco = (float)(recoforward - recobackward) / (float) (recoforward + recobackward);
	cout << "--------------------------------------------------------------\n";
	cout << "--------------------------------------------------------------\n";
	//cout << "Integral: " << ( cosReco->Integral(11,20) - cosReco->Integral(1,10) ) / cosReco->Integral()<< endl;
	std::cout << "Afb gen: " << afbgen << " N: " << forward + backward <<  "\n";
	//std::cout << "Afb reco: " << afbreco << " N: " << recoforward + recobackward << "(" << afbreco / afbgen *100 << "%)"  << "\n";
	std::cout << "Chi2: " << cosReco->Chi2Test(cosGen,"UUNORMCHI2/NDF") << " N: " << cosReco->GetEntries() << "\n";
	cout << "--------------------------------------------------------------\n";
	//cout << "Integral: " << fgen->Integral(-1,0) << " " << fgen->Integral(0,1) << endl;
	float afbgenf = (fgen->Integral(0,1) - fgen->Integral(-1,0)) / (fgen->Integral(0,1) + fgen->Integral(-1,0));
	float afbrecof = (freco->Integral(0,1) - freco->Integral(-1,0)) / (freco->Integral(0,1) + freco->Integral(-1,0));

	cout << "Afb gen functional: " << afbgenf << endl;
	cout << "Afb reco functional: " << afbrecof << "(" << afbrecof / afbgenf *100 << "%)"   << endl;
	cout << "--------------------------------------------------------------\n";
	std::cout << "Intergral: " << cosReco->Integral() << "\n";
	cout << "--------------------------------------------------------------\n";
	getAfb(freco, 1);
	//file->Close();
}
void makePretty(TH1 * vtxTotal, int color, int line = 0)
{
	vtxTotal->SetLineWidth(3);
	vtxTotal->SetLineStyle(line);
	vtxTotal->SetLineColor(color);
	vtxTotal->SetMinimum(0);
	vtxTotal->SetStats(0);
}
void getAfb(TF1* freco, float range)
{
	float a = freco->GetParameter(2);
	float b = freco->GetParameter(1);
	float c = freco->GetParameter(0);

	float da = freco->GetParError(2);
	float db = freco->GetParError(1);
	float dc = freco->GetParError(0);
	float I = 2*range* (pow(range,2) * a / 3 + c);
	float Afb = b*range/(2./3 * a * range + 2 * c);

	//float Afb = b*range/freco->Integral(-range,range);
	float deltaI = sqrt(4./9*pow(range,6) * da*da + 4*pow(range,2) * dc*dc);
	float deltaAfb = range / I * sqrt(db*db + b*b / I/I * deltaI * deltaI);
	cout << "I\t= " << I << " +- " << deltaI << endl;
	cout << "Diff\t= " << freco->Integral(-range,range) -  I << " +- " << deltaI << endl;
	cout << "Afb\t= " << Afb << " +- " << deltaAfb << "(" << deltaAfb / Afb * 100 << "%)" << endl;
}
void addBkg(TH1 * vtxTotal, string name, string cuts, int bin_e, float ratio = 1.0)
{
	string root = ".root";
	string fullcuts = "qCostheta > -1" + cuts;
	string totalpath = BKG_PATH + name + root;
	TFile * fileBkg = TFile::Open(totalpath.c_str());
	fileBkg->cd();
	TH1F * cosRecoBkgZZtmp = new TH1F(name.c_str(), "E(Ntracks)", bin_e,-1.0,1.0);
	Stats->Draw(("qCostheta>> "+name).c_str(), fullcuts.c_str());
	vtxTotal->Add(cosRecoBkgZZtmp, ratio);
}

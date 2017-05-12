#include <unistd.h>
#include <iostream>
#define MAXV 8
void cutsbbar(string filename = "TTBarProcessorLeft.root")
{
	TFile * file = TFile::Open(filename.c_str(),"READ");
	int bin_e = 100;
	int max_e = 260;
	int min_e = 50;
	TH1F * mcmass = new TH1F("mcmass", "Generated mass;m_{b#bar{b}}, GeV", bin_e,min_e,max_e);
	TH1F * mcmasscut = new TH1F("mcmasscut", ";m_{b#bar{b}}, GeV", bin_e,min_e,max_e);
	
	TTree * normaltree = Stats;
	mcmass->SetLineWidth(3);
	mcmasscut->SetLineWidth(3);
	mcmasscut->SetLineColor(kGreen);
	mcmass->SetLineColor(kBlue+1);
	
	int initcmpole = GenTree->Draw("MCMass","MCMass>0 && MCMass>200 && MCPDG == 5");
	GenTree->Draw("MCMass >> mcmass","MCMass>0  && MCPDG == 5" );
	//string cut500 = " && InvMass > 300 && bbbarAngle > 3 && bbbarP < 50";
	//string cut250 = " && InvMass > 200 && bbbarAngle > 3 && bbbarP < 20 && maxPhotonEnergy < 40";
	string cut250 = " && InvMass > 150 &&  maxPhotonEnergy < 40 && B2chargeBalance < 0.98 && B2chargeBalance > 0.05 && bbbarAngle > 2.6 && InvMass >30+150* abs(B2costheta) ";
	//string cut250 = " && InvMass > 180 &&  maxPhotonEnergy < 40";
	int cmpole = normaltree->Draw("MCMass >> mcmasscut",string("MCMass>200 && MCPDG == 5"+cut250).c_str() );
	int zpole = normaltree->Draw("MCMass >> +mcmasscut",string("MCMass>0 &&MCMass < 200"+cut250).c_str() );
	int total = normaltree->Draw("MCMass",string("InvMass >0 "+cut250).c_str() );
	std::cout << "Purity: " << (float)cmpole / total * 100 << "% (" <<total - cmpole-zpole << ")\n";

	std::cout << "Efficiency: " << (float)cmpole / (initcmpole) * 100 << "%\n";
	mcmass->Draw();
	mcmasscut->Draw("same");
}

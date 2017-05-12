
#include "/users/flc/bilokin/Processors/Macros/top/basymmetry.C"
#include "/users/flc/bilokin/Processors/Macros/top/btag.C"
void berrors(string filename = "TTBarProcessorLeft.root")
{
	
	TCanvas * c1 = new TCanvas("c1","",500,500);
	int bin_e = 20;
	int max_e = 1;
	TFile * file = TFile::Open(filename.c_str());
	const int points = 21;
	float errors[points];
	float events[points];
	float fitRange = 0.9;
	string cuts = "methodUsed > 0 && InvMass > 180 && maxPhotonEnergy < 40 && B1mass + B2mass < 130 && Sphericity < 0.25 && methodRefused == 0"; 
	for (unsigned int i = 0; i < points; i++) 
	{
		int restricted = 8000 * (i+1);
		TH1F * cosReco = new TH1F("cosReco", "E(Ntracks)", bin_e,-1.0,max_e);
		int reco  = Stats->Draw("qCostheta>>  cosReco", cuts.c_str(),"",restricted);
		cout << "NEvents: " << reco << endl;
		cosReco->Draw();
		TF1 * freco = new TF1("freco","pol2",-fitRange, fitRange);
		cosReco->Fit("freco", "QR");
		float error = 0.0;
		cout << "Total:" << restricted << endl;
		float afb = getAfb(freco, 1, error);
		errors[i] = error;
		events[i] = reco;
	}
	TGraph * gr = getGraph(points,events,errors,kBlue);
	TF1 * fg = new TF1("fg","[0]/sqrt(x)",0, 20000);
	gr->Fit("fg");
	
	cout << fg->Eval(20000) << endl;
	gr->Draw();
}

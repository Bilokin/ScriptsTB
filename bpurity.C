
void bpurity(string filename = "TTBarProcessorLeft.root")
{
	TFile *_file0 = TFile::Open(filename.c_str());
	int bin_e = 7;
	int max_e = bin_e;
	TH1F * all = new TH1F("all", "E(Ntracks)", bin_e,0,max_e);
	TH1F * good = new TH1F("good", "E(Ntracks)", bin_e,0,max_e);
	TH1F * acc = new TH1F("acc", "E(Ntracks)", bin_e,0,max_e);
	TCanvas * c1 = new TCanvas("c1", "The 3d view",0,0,500,500);
	//string cut = "MCPDG == 5 && MCMass > 200";
	string cut = "InvMass > 180 && maxPhotonEnergy < 40";
	Stats->Draw("methodTaken >>all",(cut + "&& methodUsed == 1 &&  methodRefused == 0").c_str());
	Stats->Draw("methodTaken >>acc",(cut + "&& methodUsed == 1 &&  methodRefused == 0").c_str());
	Stats->Draw("methodTaken >> good",(cut + "&& methodUsed == 1  &&  methodRefused == 0 && methodCorrect == 1").c_str(), "same");
	Stats->Draw("methodSameCharge >> +all",(cut + "&& methodRefused == 1 &&  methodUsed == 0").c_str(), "same");
	good->Divide(all);
	good->SetMaximum(1);
	good->Draw();
	//all->Draw();
	float P = 0.00;
	int total = 0;
	for (unsigned int i = 0; i < 6; i++) 
	{
		int bin = i + 2;
		float p_i = good->GetBinContent(bin);
		int N_ai = acc->GetBinContent(bin);
		cout << "p_" << i << " = " << sqrt(p_i ) << "; N_ai = " << N_ai<< endl;
		P += sqrt(p_i )* N_ai;
		total += N_ai;
	}
	P /= total;
	cout << "AVERAGE P: " << P << endl;
}

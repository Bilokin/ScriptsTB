void _gen(string filename = "TTBarProcessorLeft.root")
{
	TFile * file = TFile::Open(filename.c_str());
	float max_e = 250;
	float min_e = 50;
	int bin_e = 150;

	TH1F * genMass = new TH1F("genMass", ";m_{bb}", bin_e, min_e,max_e);
	GenTree->Draw("MCMass >> genMass","MCMass > 50 && MCPDG == 5","he");
	TF1 * breit = new TF1("breit","[0]/( pow(x*x - [1]*[1],2) + [1]*[1]*[2]*[2] )*exp(x*[3])+[4]/(x-[5])/(x-[5])",min_e,max_e);
	breit->SetParLimits(0,9.,4.5e+10);
	breit->SetParLimits(1,89.,94.5);
	breit->SetParLimits(2,2.,2.8);
	//breit->SetParLimits(3,2.e-2,2.8e-2);
	breit->SetParLimits(5,240.,260);
	TF1 * wigner=new TF1("wigner","[0]/( pow(x*x - [1]*[1],2) + [1]*[1]*[2]*[2] )*exp(x*[3])",min_e,max_e);
	//wigner->SetParLimits(1,200.,260.5);
	//wigner->SetParLimits(2,-10.,-1.5e-2);
	cout << "w: " << breit->GetParameter(1) << endl;
	breit->SetLineColor(kBlue);
	wigner->SetLineColor(kGreen);
	genMass->Fit("breit","R");
	wigner->SetParameters(breit->GetParameter(0), breit->GetParameter(1), breit->GetParameter(2), breit->GetParameter(3));

	//genMass->Fit("wigner","Rsame");
	//TF1 * total = new TF1("breitwigner","[0]/( pow(x*x - [1]*[1],2) + [1]*[1]*[2]*[2] ) + [3]/( pow(x*x - [4]*[4],2) + [4]*[4]*[5]*[5] )",min_e,max_e);
	/*TF1 * total = new TF1("breitwigner","breit+wigner",min_e,max_e);
	total->SetParLimits(0,9.,4.5e+10);
	total->SetParLimits(1,89.,94.5);
	total->SetParLimits(2,2.49,3.5);
	total->SetParLimits(3,9.,14.5e+10);
	total->SetParLimits(4,200.,260.5);
	total->SetParLimits(5,4.,4.5);*/
	//genMass->Fit("breitwigner");
	//breit->Draw("same");
	wigner->Draw("same");
	//total->Draw("same");
}

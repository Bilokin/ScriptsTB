void mcorrection()
{
	TCanvas * c1 = new TCanvas("c1", "Data-MC",0,0,500,500);
	TFile * file = TFile::Open("/exp/flc/bilokin/Training/bbbar-250GeV/recoverytest/TTBarProcessorLeftExcludedNew.root");
	Stats->Draw("B1mass:abs(B1costheta)","");
	float x4[1] = {91.18};
	float sx4[1] = {0.1};
	float dcos = 0.02;
	int npoints = 1./dcos;
	TH1F * costmp = new TH1F("costmp", "", npoints, 0,1);
	TH1F * costmp2 = new TH1F("costmp2", "", npoints, 0,1);
	TH1F * costmp3 = new TH1F("costmp3", "; |cos#theta|;<m^{jet}_{1}+m^{jet}_{2}> [GeV] ", npoints, 0,1);
	TH1F * masstmp = new TH1F("masstmp", "", npoints, 0,250);
	TH1F * masstmp2 = new TH1F("masstmp2", "", npoints, 0,200);
	TH1F * masstmp3 = new TH1F("masstmp3", ";|cos#theta|; Sum of jet masses [GeV]", npoints, 0,200);
	TH1F * costmp4 = new TH1F("costmp4", "", npoints, 0,1);
	TH1F * sphtmp3 = new TH1F("sphtmp3", "", npoints, 0,4);
	string cutbase = "abs(B1costheta)";
	string cutbase2 = "abs(B2costheta)";
	string cutbase3 = "abs(ThrustCos)";
	string base = " InvMass > 180 && methodUsed && ";
	for (unsigned int i = 0; i < npoints; i++) 
	{
		float low = i * dcos;
		float high = (i+1) * dcos;
		//Stats->Draw("B1Jetmomentum >> masstmp",(cutbase + "> " + toStr(low) + " && " + cutbase + " < " + toStr(high) ).c_str());
		/*Stats->Draw("InvMass >> masstmp",(base + cutbase + "> " + toStr(low) + " && " + cutbase + " < " + toStr(high) ).c_str());
		costmp->SetBinContent(i+1,  masstmp->GetMean());
		costmp->SetBinError(i+1,  masstmp->GetMeanError());
		Stats->Draw("B2mass >> masstmp2",(base+ cutbase2 + "> " + toStr(low) + " && " + cutbase2 + " < " + toStr(high) ).c_str());
		cout << "Range " << low << " to " << high << " mass: " << masstmp->GetMean() << endl;
		costmp2->SetBinContent(i+1,  masstmp2->GetMean());
		costmp2->SetBinError(i+1,  masstmp2->GetMeanError());*/
		//Stats->Draw("B1mass +B2mass - 1.5* abs(ThrustCos) + 2.71 * abs(ThrustCos)>> masstmp3",(cutbase3 + "> " + toStr(low) + " && " + cutbase3 + " < " + toStr(high) ).c_str());
		/*Stats->Draw("B1mass +B2mass >> masstmp3",( base + cutbase3 + "> " + toStr(low) + " && " + cutbase3 + " < " + toStr(high) ).c_str());
		costmp3->SetBinContent(i+1,  masstmp3->GetMean());pow(abs(ThrustCos),3)
		costmp3->SetBinError(i+1,  masstmp3->GetMeanError());*/
		//Stats->Draw("Sphericity+7.136e-02*abs(ThrustCos)-6.43e-01*pow(abs(ThrustCos),2)+1.7628e+00*pow(abs(ThrustCos),3)-1.9431e+00*pow(abs(ThrustCos),4)+7.7519e-01*pow(abs(ThrustCos),5) >> sphtmp3",( base + cutbase3 + "> " + toStr(low) + " && " + cutbase3 + " < " + toStr(high) ).c_str());
		Stats->Draw("bbbarAngle >> sphtmp3",( base + cutbase + "> " + toStr(low) + " && " + cutbase + " < " + toStr(high) ).c_str());
		costmp4->SetBinContent(i+1,  sphtmp3->GetMean());
		costmp4->SetBinError(i+1,  sphtmp3->GetMeanError());
	}
	costmp3->SetMinimum(0);
	costmp3->SetLineWidth(3);
	costmp3->SetLineStyle(1);
	costmp3->SetLineColor(kBlue);
	costmp3->SetFillColor(kBlue);
	costmp3->SetFillStyle(3004);
	costmp3->SetStats(0);

	//costmp->Draw();
	//costmp3->Fit("pol2");
	//costmp2->Draw("same");
	costmp3->Draw("he");
	TCanvas * c2 = new TCanvas("c2", "Data-MC",0,500,500,500);
	costmp4->Draw();
	TCanvas * c3 = new TCanvas("c3", "Data-MC",500,500,500,500);
	costmp->Draw();

}
string toStr(float pi)
{
	stringstream stream;
	stream << pi;
	string s = stream.str();
	return s;
}


void topgamma(string filename = "TTBarProcessorLeft.root")
{
	TCanvas * c1 = new TCanvas("c1", "Data-MC",0,0,500,500);
	TFile * file = TFile::Open(filename.c_str(),"READ");
	int bin_e = 100;
	int max_e = 3;
	float min_e = 0.8;
	TH2F * top = new TH2F("top", ";Hadronic #gamma_{t};Leptonic #gamma_{t}", bin_e,min_e,max_e, bin_e,min_e,max_e);
	Stats->Draw("Top2gamma:Top1gamma >> top", "Top2gamma < 3");
	TF1 * cuts = new TF1("cuts","2.5-x",min_e, max_e);
	cuts->SetLineStyle(2);
	cuts->SetLineWidth(3);
	//cuts->SetLineColor(kGray+1);
	top->SetStats(0);
	top->Draw("colz");
	cuts->Draw("same");
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.12);
}

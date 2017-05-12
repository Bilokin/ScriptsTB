#include "../vertex/ratioplot.C"
void chargeTVCM(string filenameLeft = "TTBarProcessorLeft.root", string filenameMC = "MCTest.root")
{
	
	int bin_e = 50;
	float max_e = 5;
	float min_e = 0;
	int min_n = 0;
	int max_n = 13;
	TFile * fileL = TFile::Open(filenameLeft.c_str());

	//TCanvas *c1 = new TCanvas("c1", "Data-MC",0,0,1000,500);
	//c1->Divide(2,1);
	TH1F * correct = new TH1F("correct", "Reco;d_{B}, mm", bin_e,min_e,max_e);
	correct->Sumw2();
	TH1F * incorrect = new TH1F("incorrect", "Reco;d_{B}, mm", bin_e,min_e,max_e);
	incorrect->Sumw2();
	//c1->cd(1);
	//Stats->Draw("Top1bdistance >> correct1","Top1bTVCM != 0","");Stats->Draw("Top2bdistance >> correct1","Top2bTVCM != 0","");
	TH1F * allRecotracks = new TH1F("allRecotracks", "Reco;N_{tracks}", max_n,min_n,max_n);
	TH1F * badRecotracks = new TH1F("badRecotracks", "Reco;N_{tracks}", max_n,min_n,max_n);
	Stats->Draw("Top1bdistance>> correct","Top1bTVCM != 0","");// && Top1bcharge == 0
	Stats->Draw("Top2bdistance>> +correct","Top2bTVCM != 0","");
	Stats->Draw("Top1bntracks >> allRecotracks","Top1bTVCM != 0","");// && Top1bcharge == 0
	Stats->Draw("Top2bntracks >> +allRecotracks","Top2bTVCM != 0","");
	makePretty(correct,kBlue+1);
	makePretty(incorrect,kBlack);
	makePretty(allRecotracks, kBlue+1);
	makePretty(badRecotracks, kBlack);
	//correct->Draw("h");
	Stats->Draw("Top1bdistance>> incorrect","Top1bTVCM < 0","same");
	Stats->Draw("Top2bdistance>> +incorrect","Top2bTVCM < 0","same");
	Stats->Draw("Top1bntracks >> badRecotracks","Top1bTVCM < 0","");// && Top1bcharge == 0
	Stats->Draw("Top2bntracks >> +badRecotracks","Top2bTVCM < 0","");
	cout << "Total entries: " << correct->GetEntries() << " incorrect: " << incorrect->GetEntries() << "(" << incorrect->GetEntries() / correct->GetEntries()*100 << "%)" << endl;
	TFile * filemc = TFile::Open(filenameMC.c_str());
	TCanvas *c4 = new TCanvas("c4", "Data-MC",0,0,500,700);
	TCanvas *c3 = new TCanvas("c3", "Data-MC",0,0,500,700);
	TCanvas *c2 = new TCanvas("c2", "Data-MC",0,0,500,700);
	TCanvas *c1 = new TCanvas("c1", "Data-MC",0,0,500,700);
	TH1F * allGentracks = new TH1F("allGentracks", "Gen;N_{tracks}", max_n,min_n,max_n);
	TH1F * badGentracks = new TH1F("badGentracks", "Gen;N_{tracks}", max_n,min_n,max_n);
	TH1F * correctMC = new TH1F("correctMC", "Gen;d_{B}, mm", bin_e,min_e,max_e);
	TH1F * incorrectMC = new TH1F("incorrectMC", "Gen;d_{B}, mm", bin_e,min_e,max_e);
	Stats->Draw("btotalnumber >> allGentracks","btotalnumber > 0  && boscillation != 0");
	Stats->Draw("bbartotalnumber >> +allGentracks","bbartotalnumber > 0  && bbaroscillation != 0");
	Stats->Draw("bbartotalnumber >> badGentracks","bbartotalnumber > 0 && bbaroscillation == -1");
	Stats->Draw("btotalnumber >> +badGentracks","btotalnumber > 0 && boscillation == -1");

	Stats->Draw("bIPdistance + bdistance >> correctMC","boscillation != 0 && bcharge == 0");
	Stats->Draw("bIPdistance + bdistance >> incorrectMC","boscillation == -1 && bcharge == 0");
	Stats->Draw("bbarIPdistance + bbardistance >> +correctMC","bbaroscillation != 0 && bcharge == 0");
	Stats->Draw("bbarIPdistance + bbardistance >> +incorrectMC","bbaroscillation == -1 && bcharge == 0");
	correctMC->SetLineWidth(3);
	correctMC->SetLineColor(kBlue+1);
	correctMC->SetMinimum(0.01);
	incorrectMC->SetLineWidth(3);
	incorrectMC->SetLineColor(kBlack);
	makePretty(allGentracks,kBlue+1);
	makePretty(badGentracks,kBlack);
	//correctMC->Draw("h");
	//incorrectMC->Draw("hsame");
	//allGentracks->Draw("h");
	//badGentracks->Draw("hsame");
	ratioplot(incorrect, correct, c4);
	ratioplot(badGentracks, allGentracks,c1);
	ratioplot(badRecotracks, allRecotracks, c2);
	ratioplot(incorrectMC, correctMC, c3);
	cout << "Total MC entries: " << allGentracks->GetEntries() << " incorrect: " << badGentracks->GetEntries() << "(" << badGentracks->GetEntries() / allGentracks->GetEntries()*100 << "%)" << endl;
}
void makePretty(TH1 * vtxTotal, int color)
{
	vtxTotal->SetLineWidth(3);
	vtxTotal->SetLineColor(color);
	vtxTotal->SetMinimum(0);
	vtxTotal->SetStats(0);
}

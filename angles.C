

#include <unistd.h>
#include <iostream>
#define MAXV 8
void angles(string filenameLeft = "TTBarProcessorLeft.root", string filenameRight = "TTBarProcessorRight.root")
{
	int bin_e = 150;
	float max_e = 3.14;
	float min_e = 0;
	TH1F * topangleL = new TH1F("topangleL1", "E(Ntracks)", bin_e,min_e,max_e);
	TH1F * topangleR = new TH1F("topangleR1", ";m_{W|t}, GeV", bin_e,min_e,max_e);
	
	TFile * fileL = TFile::Open(filenameLeft.c_str(),"READ");
	//fileL->cd();
	
	Stats->Project("topangleL1","cos(MCTopBangle)","MCTopBangle > 0");
	topangleL1->Rebin();
	topangleL1->SetLineWidth(3);
	topangleL1->SetLineColor(kBlue+1);
	topangleL1->DrawNormalized("");
	
	TFile * fileR = TFile::Open(filenameRight.c_str(),"READ");
	fileR->cd();
	Stats->Project("topangleR1","cos(MCTopBangle)","MCTopBangle > 0","same");
	topangleR1->SetLineWidth(3);
	topangleR1->SetLineColor(kGreen);
	topangleR1->DrawNormalized("same");
	cout << "eLpR: " << topangleL1->GetMean() << " eRpL: " << topangleR1->GetMean() << endl;

}

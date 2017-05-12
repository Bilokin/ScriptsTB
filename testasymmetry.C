
#include <unistd.h>
#include <iostream>
#define MAXV 8
void testasymmetry(string filename = "TTBarProcessor.root", TCanvas * c1 = NULL)
{
	TFile * file = TFile::Open(filename.c_str());
	int bin_e = 20;
	int max_e = 1;
	if (!c1) 
	{
		c1 = new TCanvas("c1", "Data-MC",0,0,600,600);
	}
	//c1->Divide(2,1);
	TH1F * cosReco = new TH1F("cosReco", "E(Ntracks)", bin_e,-1.0,max_e);
	//cosReco->Sumw2();
	TH1F * cosGen = new TH1F("cosGen", ";cos(#theta_{top})", bin_e,-1.0,max_e);

	TTree * normaltree = Stats;
	//cosReco->SetLineColor(kBlue);
	cosReco->SetLineWidth(3);
	cosGen->SetLineWidth(3);
	cosGen->SetLineStyle(2);
	cosGen->SetLineColor(kGreen+1);
	int forward = normaltree->Draw("qMCcostheta >> cosGen","qMCcostheta > 0");
	int backward = normaltree->Draw("qMCcostheta >> +cosGen","qMCcostheta < 0");
	
	/*int forward = normaltree->Draw("cosMCTop >> cosGen","mctag == 1 && cosMCTop > 0");
	int backward = normaltree->Draw("cosMCTop >> +cosGen","mctag == 1 && cosMCTop < 0");
	forward += normaltree->Draw("-cosMCTopBar >> +cosGen","mctag == 1 && cosMCTopBar < 0");
	backward += normaltree->Draw("-cosMCTopBar >> +cosGen","mctag == 1 && cosMCTopBar > 0");
	
	int recoforward = normaltree->Draw("-cosTop1 >> cosReco","qTop1 > 0 && chiHad < 15 && btagTop1 > 0.950 && cosTop1 < 0");
	int recobackward = normaltree->Draw("-cosTop1 >> +cosReco","qTop1 > 0 && chiHad < 15 && btagTop1 > 0.950 && cosTop1 > 0");
	recoforward += normaltree->Draw("cosTop1 >> +cosReco","qTop1 < 0 && chiHad < 15 && btagTop1 > 0.950 && cosTop1 > 0");
	recobackward += normaltree->Draw("cosTop1 >> +cosReco","qTop1 < 0 && chiHad < 15 && btagTop1 > 0.950 && cosTop1 < 0");
	*/
	int recoforward = normaltree->Draw("Top1costheta >> cosReco", "(log(Top1bmomentum)/Top1bntracks/2)*(Top1costheta > 0 && Top1bmomentum > 1.0 && Top1bntracks > 2)");
	int recobackward = normaltree->Draw("Top1costheta >> +cosReco", "(log(Top1bmomentum)/Top1bntracks/2)*(Top1costheta < 0 && Top1bmomentum > 1.0 && Top1bntracks > 2)");
	cosGen->SetMinimum(0);
	cosGen->SetStats(0);
	if (filename == "TTBarProcessorRight.root") 
	{
		cosGen->SetTitle("e_{R}p_{L}");
	}
	if (filename == "TTBarProcessorLeft.root") 
	{
		cosGen->SetTitle("e_{L}p_{R}");
	}
	//cosGen->DrawNormalized();
	//cosGen->Draw();
	//cosGen->SetMinimum(0);
	//cosReco->DrawNormalized("samee");
	cosReco->Draw("e");
	TLegend *legendMean2 = new TLegend(0.15,0.75,0.6,0.85,NULL,"brNDC");
        legendMean2->SetFillColor(kWhite);
        legendMean2->SetBorderSize(0);
        legendMean2->AddEntry(cosGen,"Generated","f");
        legendMean2->AddEntry(cosReco,"Reconstructed","f");
	legendMean2->Draw();
cout << "Integral: " << ( cosReco->Integral(10,20) - cosReco->Integral(0,9) ) / cosReco->Integral()<< endl;
	std::cout << "Afb gen: " << (float)(forward - backward) / (float) (forward + backward) << " N: " << forward + backward <<  "\n";
	std::cout << "Afb reco: " << (float)(recoforward - recobackward) / (float) (recoforward + recobackward) << " N: " << recoforward + recobackward  << "\n";
	//file->Close();
}
	

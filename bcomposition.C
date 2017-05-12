#include <sstream> 
void bcomposition(string filename = "TTBarProcessorLeft.root")
{
	string names[6] = {"VTX+VTX","KAON+KAON","VTX1+KAON1","VTX2+KAON2","VTX1+KAON2","VTX2+KAON1"};
	TFile *_file0 = TFile::Open(filename.c_str());
	TCanvas * c1 = new TCanvas("c1", "Data-MC",0,0,500,500);
	c1->Divide(3,2);
	int bin_e = 7;
	int max_e = bin_e;
	c1->cd(1);
	TH1F * all = new TH1F("all", "E(Ntracks)", bin_e,0,max_e);
	TH1F * good = new TH1F("good", "E(Ntracks)", bin_e,0,max_e);
	TH1F * acc = new TH1F("acc", "E(Ntracks)", bin_e,0,max_e);
	TH1F * bad = new TH1F("bad", "E(Ntracks)", bin_e,0,max_e);
	string cut = "InvMass > 180 && maxPhotonEnergy < 40";
	Stats->Draw("methodTaken >>all",(cut + "&& methodUsed == 1 &&  methodRefused == 0").c_str());
	Stats->Draw("methodTaken >>acc",(cut + "&& methodUsed == 1 &&  methodRefused == 0").c_str());
	Stats->Draw("methodTaken >> good",(cut + "&& methodUsed == 1  &&  methodRefused == 0 && methodCorrect == 1").c_str(), "same");
	Stats->Draw("methodSameCharge >> +all",(cut + "&& methodRefused == 1 &&  methodUsed == 0").c_str(), "same");
	Stats->Draw("methodSameCharge >> bad",(cut + "&& methodRefused == 1 &&  methodUsed == 0").c_str(), "same");
	//good->Divide(all);
	int bin_c = 20;
	int max_c = 1;
	TH1F * cor = new TH1F("cor", "E(Ntracks)", bin_c,-1,max_c);
	GenTree->Draw("qMCcostheta >> cosGen","qMCcostheta > -1 && MCMass > 200 && MCPDG == 5");
	TF1 * fgen = new TF1("fgen","pol2",-1,1);
	cosGen->Fit("fgen","Q");
	float afbg = getAfb(fgen, 1);
	fgen->SetLineColor(kGreen);
	fgen->SetLineStyle(3);
	string taken  = "methodTaken == ";
	string nottaken  = "methodSameCharge == ";
	string used  = "&& methodUsed == 1 &&  methodRefused == 0 &&";
	string refused  = "&& methodUsed == 0 &&  methodRefused == 1 &&";
	int mid = (float)bin_c / 2;
	float fitrange = 0.95;
	cout << "Afbg: " << afbg << endl;

	for (unsigned int i = 0; i < 6; i++) 
	{
		int bin = i + 2;
		c1->cd(i+1);
		float N_ai = acc->GetBinContent(bin);
		float N_ri = bad->GetBinContent(bin);
		//float p_i = 0.5 + sqrt((N_ai - N_ri)/(N_ri + N_ai))/2;
		float p_i = sqrt(good->GetBinContent(bin) / all->GetBinContent(bin));
		cout << "\t" << names[i]<<endl; string name = names[i];
		cout << "\tN_ai = " << N_ai << "; N_ri = " << N_ri << ";\tp_i = " << p_i << endl;
		TH1F * cosReco = new TH1F("cosReco", ";cos#theta_{b}", bin_c,-1.0,max_c);
		makePretty(cosReco, kBlack);
		TH1F * cosGen = new TH1F("cosGen", (name+";cos#theta_{b}").c_str(), bin_c,-1.0,max_c);
		TH1F * cosRef = new TH1F("cosRef", (name+";cos#theta_{b}").c_str(), bin_c,-1.0,max_c);
		makePretty(cosGen, kGreen+1, 2);
		makePretty(cosRef, kRed+1, 2);
		GenTree->Draw("qMCcostheta >> cosGen","qMCcostheta > -1 && MCMass > 200 && MCPDG == 5");
		Stats->Draw("qCostheta1>> cosReco", (cut + used+taken + ToString(i+1)).c_str(), "same");
		Stats->Draw("qCostheta2>> +cosReco", (cut + used+taken + ToString(i+1)).c_str(),"same");
		Stats->Draw("B1costheta>> cosRef", (cut + refused+nottaken + ToString(i+1)).c_str(), "same");
		Stats->Draw("B2costheta>> +cosRef", (cut + refused+nottaken + ToString(i+1)).c_str(),"same");
		cosGen->Scale(1.0/ (cosGen->GetBinContent(mid) + cosGen->GetBinContent(mid+1)) * (cosReco->GetBinContent(mid) + cosReco->GetBinContent(mid+1)));
		TH1F * corrected = computeCorrection(cosReco, p_i,bin_c, cosRef);
		makePretty(corrected, kBlue+1);
		TF1 * freco_i = new TF1("freco_i","pol2",-fitrange, fitrange);
		freco_i->SetLineStyle(3);
		corrected->Fit("freco_i", "QR");
		cor->Fit("freco_i", "QR");
		float afbr_i = getAfb(freco_i, fitrange);

		cout << "\tResult: " << afbr_i / afbg *100 << "%\n"; 
		cosGen->Draw();
		cosReco->Draw("samee");
		corrected->Draw("samee");
		cosRef->Draw("same");
		cor->Add(corrected);
	}
	TCanvas * c2 = new TCanvas("c2", "Data-MC",0,0,500,500);
	c2->cd();
	TH1F * cosGen1 = new TH1F("cosGen1", ";cos#theta_{b}", bin_c,-1.0,max_c);
	GenTree->Draw("qMCcostheta >> cosGen1","qMCcostheta > -1 && MCMass > 200 && MCPDG == 5");
	makePretty(cosGen1, kGreen+1, 2);
	makePretty(cor, kBlue+1);
	cosGen1->Scale(1.05*(float)cor->Integral() / cosGen1->Integral());
	TF1 * frecof = new TF1("frecof","pol2", -fitrange, fitrange);
	frecof->SetLineStyle(3);
	cor->Fit("frecof", "QR");
	float afbr = getAfb(frecof, 0.95);
	cout << "RESULT: " << afbr/afbg * 100 << "%\n";
	cosGen1->Draw();
	cor->Draw("same");
}
TH1F * computeCorrection(TH1F * right, float p, int bin_e = 20, TH1F * ref = NULL)
{
	
	string hname = "correctedP";
	TH1F * corrected = new TH1F(hname.c_str(),";cos",bin_e,-1, 1);
	int niterations = (float) bin_e/2;
	int total = 0;
	for (unsigned int i = 1; i < niterations+1; i++) 
	{
		int j = bin_e + 1 - i;
		float naminus_i = right->GetBinContent(i);
		float naplus_i = right->GetBinContent(j);
		float p_i = p;
		float q_i = 1. - p_i;
		float nplus_i = Nplus(p_i,q_i, naplus_i, naminus_i) ;
		float nminus_i = Nplus(q_i,p_i, naplus_i, naminus_i);
		if (ref) 
		{
			cout << "\tN_total: " << naminus_i + naplus_i + ref->GetBinContent(j) + ref->GetBinContent(i) << " vs " << nplus_i + nminus_i << endl;
			float diff = nminus_i + nplus_i - (naminus_i + naplus_i + ref->GetBinContent(j) + ref->GetBinContent(i));
			nplus_i = nplus_i -  diff / 2;
			nminus_i = nminus_i -  diff / 2;
		}
		total += nplus_i + nminus_i;
		//float n2plus_i = p_i*p_i * nplus_i  + q_i*q_i * nminus_i;
		//float n2minus_i =  p_i*p_i*nminus_i  + q_i*q_i * nplus_i;
		float n2plus_i = p_i*p_i * nplus_i ;// + q_i*q_i * nminus_i;
		float n2minus_i = p_i*p_i * nminus_i;//  + q_i*q_i * nplus_i;
		float error_minus = sigma(naminus_i, naplus_i, p_i,q_i, 1, 1);
		float error_plus = sigma(naplus_i, naminus_i, p_i,q_i, 1, 1);
		corrected->SetBinContent(i, n2minus_i);
		corrected->SetBinContent(j, n2plus_i);
		corrected->SetBinError(i,error_minus);
		corrected->SetBinError(j, error_plus);
//		cout << "Corr: " << n2minus_i << endl;
	}
	return corrected;
}
void makePretty(TH1 * vtxTotal, int color, int line = 0)
{
	vtxTotal->SetLineWidth(3);
	vtxTotal->SetLineStyle(line);
	vtxTotal->SetLineColor(color);
	vtxTotal->SetMinimum(0);
	vtxTotal->SetStats(0);
}
string ToString(int i)
{
	ostringstream s;
	s << i;
	return s.str();
}
float Nplus(float p, float q, float Naplus, float Naminus)
{
	return (Naplus * p*p - Naminus * q*q)/(p*p*p*p - q*q*q*q);
}
float sigma(float naplus, float naminus, float p, float q, float nr, float n)
{
	float alpha = abs(1/(p*p*p*p-q*q*q*q));
	//float sigmap = sigmaP(nr,n);//sqrt(dpdNa*dpdNa*n  + dpdNr * dpdNr * nr);
	return alpha * sqrt( naplus * p*p*p*p + naminus * q*q*q*q );//+ sigmap * sigmap * pow(2*p*naplus + 2 *q* naminus ,2));
}
float getAfb(TF1* freco, float range)
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
	//cout << "I\t= " << I << " +- " << deltaI << endl;
	//cout << "Diff\t= " << freco->Integral(-range,range) -  I << " +- " << deltaI << endl;
	//cout << "Afb\t= " << Afb << " +- " << deltaAfb << "(" << deltaAfb / Afb * 100 << "%)" << endl;
	return Afb;
}

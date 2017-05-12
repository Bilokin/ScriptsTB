#include <sstream>
void qp(string filename = "TTBarProcessorLeft.root")
{
	TCanvas * c1 = new TCanvas("c1","",1000,500);
	c1->Divide(2,1);
	c1->cd(1);
	TFile * file = TFile::Open(filename.c_str());
	int bin_e = 60;
	int max_e = 1.;
	int middle = (float)bin_e/2;
	TH1F * right = new TH1F("right",";cos",bin_e,-max_e, max_e); 
	TH1F * rightg = new TH1F("rightg",";cos",bin_e,-max_e, max_e); 
	TH1F * wrong = new TH1F("wrong",";cos",bin_e,-max_e, max_e);
	TH1F * gen = new TH1F("gen",";cos#theta_{b}",bin_e,-max_e, max_e);
	makePrettyQ(gen,kGreen);
	gen->SetStats(0);
	TH1F * corrected = new TH1F("corrected",";cos#theta",bin_e,-max_e, max_e);
	TH1F * pglobalHist = new TH1F("pglobal",";cos#theta;p",bin_e,-max_e, max_e);
	TH1F * polarAngle = new TH1F("polarAngle",";cos",bin_e,-max_e, max_e); 
	polarAngle->SetStats(0);
	makePrettyQ(polarAngle,kBlack);
	float Nrefused = Stats->Draw("methodSameCharge","methodSameCharge > 0 && MCMass > 200");
	float Naccepted = Stats->Draw("methodTaken","methodTaken > 0 && MCMass > 200");
	std::cout << "Naccepted: " << Naccepted <<" Nrefused: " << Nrefused << endl;
	float N = Nrefused + Naccepted;
	float p = 0.5 + sqrt(1-2*Nrefused/N)/2;
	float q = 1.-p;
	//float pvalue = getPGlobal();
	Stats->Draw("qCostheta1>>right","(methodTaken > 0)&& MCMass > 200");
	Stats->Draw("qCostheta2>>+right","(methodTaken > 0) && MCMass > 200");
	//Stats->Draw("qCostheta1>>right","(methodTaken != 3 && methodTaken != 4)&& MCMass > 200");
	//Stats->Draw("qCostheta2>>+right","(methodTaken != 3 && methodTaken != 4) && MCMass > 200");*/
	Stats->Draw("qCostheta1>>polarAngle","MCMass > 200 && methodUsed ");
	Stats->Draw("qCostheta2>>+polarAngle","MCMass > 200 && methodUsed ");
	GenTree->Draw("qMCcostheta>>gen","MCMass > 200 && MCPDG == 5");
	Stats->Draw("B1costheta >> wrong","(methodSameCharge >0) && MCMass > 200");
	Stats->Draw("B2costheta >> +wrong","(methodSameCharge >0) && MCMass > 200");
	//Stats->Draw("B1costheta >> wrong","(methodSameCharge != 3 && methodSameCharge !=4) && MCMass > 200");
	//Stats->Draw("B2costheta >> +wrong","(methodSameCharge != 3 && methodSameCharge != 4) && MCMass > 200");

/*	Stats->Draw("B1costheta >> +wrong","(methodSameCharge == 3 || methodSameCharge ==4) && MCMass > 200");
	Stats->Draw("B2costheta >> +wrong","(methodSameCharge == 3 || methodSameCharge ==4) && MCMass > 200");
*/
	std::cout << "Global: p = "<< p << "; q = " << q << endl;
	makePrettyQ(right,kBlack);
	makePrettyQ(rightg,kBlack);
	makePrettyQ(wrong,kRed);
	gen->SetLineStyle(2);
	std::cout << "-------------------------\n";
	corrected = computeCorrection(right,wrong,bin_e, "", pglobalHist);
	makePrettyQ(corrected,kBlue);
	gen->Scale(1./(gen->GetBinContent(middle) + gen->GetBinContent(middle+1)) * (right->GetBinContent(middle) + right->GetBinContent(middle+1)));
	//corrected->Scale(1./corrected->Integral() * right->GetEntries());
	gen->Draw("");
	right->Draw("samee");
	corrected->Draw("samee");
	wrong->Draw("samee");
	//*/

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//return;
/*	std::cout << "-------------------------\n";
	TH1F * correctedg = new TH1F("correctedg",";cos",bin_e,-max_e, max_e);
	TH1F * refusedg = new TH1F("refusedg",";cos",bin_e,-max_e, max_e);
	makePretty(correctedg,kBlue);
	makePretty(refusedg,kRed);
	string cutssame = "MCMass > 200 && methodSameCharge == ";
	string cutsused = "MCMass > 200 && methodTaken == ";
	string cutsall = "MCMass > 200 && methodTaken == ";
	for (unsigned int i = 1; i < 7; i++) 
	{
		string n = str(i);
		if(i != 1) continue;
		float Nrefused_i = Stats->Draw("methodSameCharge",string(cutssame+n).c_str());
		float Naccepted_i = Stats->Draw("methodTaken",string(cutsused+n).c_str());
		float p_i = 0.5 + sqrt(1-2*Nrefused_i/(Nrefused_i + Naccepted_i))/2;
		float q_i = 1. - p_i;
		std::cout << "Naccepted_"<<n<<": " << Naccepted_i <<" Nrefused_"<<n<<": " << Nrefused_i << endl;
		std::cout << "Local: p_"<<n<<" = "<< p_i << "; q_"<<n<<" = " << q_i << endl << endl;
		TH1F * rightn = new TH1F(string("right"+n).c_str(),";cos",bin_e,-max_e, max_e); 
		TH1F * wrongn = new TH1F(string("wrong"+n).c_str(),";cos",bin_e,-max_e, max_e);
		TH1F * alln = new TH1F(string("all"+n).c_str(),";cos",bin_e,-max_e, max_e);
		makePretty(rightn,kBlack);
		makePretty(wrongn,kRed);
		Stats->Draw(string("qCostheta1>>right"+n).c_str(),string(cutsused+n).c_str());
		
		Stats->Draw(string("qCostheta2>>+right"+n).c_str(),string(cutsused+n).c_str());
		Stats->Draw(string("B1costheta >> wrong"+n).c_str(),string(cutssame+n).c_str(),"same");
		Stats->Draw(string("B2costheta >> +wrong"+n).c_str(),string(cutssame+n).c_str(),"same");
		Stats->Draw(string("B1costheta >> all"+n).c_str(),string(cutsall+n).c_str(),"same");
		Stats->Draw(string("B2costheta >> +all"+n).c_str(),string(cutsall+n).c_str(),"same");
		//wrongn->Divide(alln);
		//cout << "W: " << wrongn->Integral() <<endl;
		//wrongn->Scale(Nrefused_i/wrongn->Integral() );
		TH1F * correctedl = computeCorrection(rightn,wrongn,bin_e,n);
		std::cout << "Adding " << correctedl->GetEntries() << " events\n";
		rightg->Add(rightn);
		correctedg->Add(correctedl);
		refusedg->Add(wrongn);		
	}
	std::cout << "-------------------------\n";
	std::cout << "********************\nCHECK:\n";
	TH1F * corrected2 = computeCorrection(rightg,refusedg,bin_e,"check");
	makePretty(corrected2,kBlue);
	cout << "Integral: " << rightg->GetEntries() << endl;
	std::cout << "********************\n";
	c1->cd(2);
	//gen->Scale(1./gen->GetEntries() * rightg->GetEntries());
	rightg->Draw("e");
	gen->Draw("same");
	//correctedg->Scale(1./(correctedg->GetBinContent(middle) + correctedg->GetBinContent(middle+1)) * (rightg->GetBinContent(middle) + rightg->GetBinContent(middle+1)));
	//corrected2->Scale(1./(corrected2->GetBinContent(middle) + corrected2->GetBinContent(middle+1)) * (rightg->GetBinContent(middle) + rightg->GetBinContent(middle+1)));
	//correctedg->Scale(1./correctedg->Integral() * rightg->GetEntries());
	//refusedg->Scale(1./refusedg->Integral() * right->GetEntries()/4.);
	//correctedg->Draw("samee");///
	correctedg->Draw("samee");///
	refusedg->Draw("samee");///
	float result = getAfb(corrected, bin_e);
	std::cout << "CORRECTED AFB: " << result << " (" << result / getAfb(gen, bin_e) * 100 << "%)"<< "\n";*/
	c1->cd(2);
	TH1F * correctedg = computeCorrectionP(polarAngle, pglobal, bin_e);
	//TH1F * correctedg = computeCorrectionF(polarAngle, pvalue, bin_e);
	makePrettyQ(correctedg,kBlue);
	gen->Scale(1./(gen->GetBinContent(middle) + gen->GetBinContent(middle+1)) * (polarAngle->GetBinContent(middle) + polarAngle->GetBinContent(middle+1)));
	gen->Draw("");
	polarAngle->Draw("samee");
	gen->Draw("same");
	correctedg->Draw("samee");
	float result = getAfb(correctedg, bin_e);
	std::cout << "Histogram Afb: " << result << " (" << result / getAfb(gen, bin_e) * 100 << "%)"<< "\n";
	TF1 * freco = new TF1("freco","pol2",-0.9,0.9);
	freco->SetLineStyle(3);
	freco->SetLineColor(kBlue);
	correctedg->Fit("freco", "QR");
	float afbrecof = (freco->Integral(0,1) - freco->Integral(-1,0)) / (freco->Integral(0,1) + freco->Integral(-1,0));
	std::cout << "CORRECTED AFB: " << afbrecof << " (" << afbrecof / getAfb(gen, bin_e) * 100 << "%)"<< "\n";
	TLegend *legendMean2 = new TLegend(0.15,0.65,0.6,0.85,NULL,"brNDC");
        legendMean2->SetFillColor(kWhite);
        legendMean2->SetBorderSize(0);
        legendMean2->AddEntry(gen,"Generated","f");
        legendMean2->AddEntry(right,"Reconstructed","f");
        legendMean2->AddEntry(correctedg,"Corrected","f");
        //legendMean2->AddEntry(refusedg,"Refused","f");
	legendMean2->Draw();
	TCanvas * c2 = new TCanvas("c2","",500,500);
	c2->cd();
	pglobal->Draw("e");
	 equation();

}
template<class T>
string str(T t)
{
	ostringstream convert;
	convert << t;
	return convert.str();
}
float getAfb(TH1F * right, int bin_e = 20)
{
	int niterations = (float) bin_e/2;
	float naminus = 0.0;
	float naplus = 0.0;

	for (unsigned int i = 1; i < niterations+1; i++) 
	{
		int j = bin_e + 1 - i;
		float naminus_i = right->GetBinContent(i);
		float naplus_i = right->GetBinContent(j);
		naminus += naminus_i;
		naplus += naplus_i;

	}
	float afbmeasured =  (naplus-naminus)/(naplus+naminus);
	return afbmeasured;
}
float Nplus(float p, float q, float Naplus, float Naminus)
{
	return (Naplus * p*p - Naminus * q*q)/(p*p*p*p - q*q*q*q);
}
float Nminus(float p, float q, float Naplus, float Naminus)
{
	return (Naplus * q*q - Naminus * p*p)/(q*q*q*q - p*p*p*p );
}
void makePrettyQ(TH1 * vtxTotal, int color)
{
	vtxTotal->SetLineWidth(3);
	vtxTotal->SetLineColor(color);
	vtxTotal->SetMinimum(0);
	vtxTotal->SetStats(0);
}
float getPGlobal(int bin_e = 20, float & error = 0)
{
	TH1F * myp = new TH1F("correctedP",";cos",bin_e,-1, 1);	
	//if (nmethods == 1) 
	{
		//Stats->Draw("methodTaken>>goodEvents","InvMass > 180 && maxPhotonEnergy < 40  && B1mass + B2mass < 130 && methodUsed == 1 && methodRefused == 0");
		//Stats->Draw("methodSameCharge>>badEvents","InvMass > 180 && maxPhotonEnergy < 40  && B1mass + B2mass < 130 && methodRefused == 1 && methodUsed == 0");
		Stats->Draw("methodTaken>>goodEvents","InvMass > 180 && maxPhotonEnergy < 40  && B1mass + B2mass < 120&& Sphericity < 0.2  && methodUsed == 1 && methodRefused == 0");
		Stats->Draw("methodSameCharge>>badEvents","InvMass > 180 && maxPhotonEnergy < 40  && B1mass + B2mass < 120 && Sphericity < 0.2 && methodRefused == 1 && methodUsed == 0");
		//Stats->Draw("methodTaken>>goodEvents","InvMass > 180 && maxPhotonEnergy < 40  && B1mass + B2mass < 130 ");
		//Stats->Draw("methodSameCharge>>badEvents","InvMass > 180 && maxPhotonEnergy < 40  && B1mass + B2mass < 130");
	}
	int b1 = badEvents->GetBinContent(4);
	int b2 = badEvents->GetBinContent(5);
	//badEvents->SetBinContent(4,  b1);
	//badEvents->SetBinContent(5,  b2);
	float averageP = 0.;
	float averagedP = 0.;
	float nevents = 0;
	float p_VTX = 0.;
	float p_KAON = 0.;
	float pvalues[6];
	for (unsigned int i = 2; i < 4; i++) 
	{
		float Naccepted_i = goodEvents->GetBinContent(i);
		float Nrefused_i = badEvents->GetBinContent(i);
		std::cout << "Refused: Nr = " << Nrefused_i <<endl;
		std::cout << "Accepted: Na = " << Naccepted_i <<endl;
		float p_i = 0.5 + sqrt((Naccepted_i - Nrefused_i)/(Nrefused_i + Naccepted_i))/2;
		/*pvalues[i - 2] = p_i;
		if (i > 3) 
		{
			if (i == 4 || i == 5) 
			{
				//continue;
			}
			float p_VTX = pvalues[0];
			float p_KAON = pvalues[1];
			float Nr = (p_VTX * (1 - p_KAON) + p_KAON * (1 - p_VTX)) * (Naccepted_i + Nrefused_i);
			//float Na = Nrefused_i + Naccepted_i - Nr; //(p_VTX * p_KAON + p_KAON * p_VTX) * (Naccepted_i + Nrefused_i);
			float Na = (p_VTX * p_KAON + (1-p_KAON) *(1- p_VTX)) * (Naccepted_i + Nrefused_i);
			//std::cout << "Predicted: Nr = " << Nr <<endl;
			//std::cout << "Predicted: Na = " << Na <<endl;
			//float pk = (Naccepted_i * p_VTX - Nrefused_i) / (Naccepted_i + Nrefused_i) / (2*p_VTX -1);
			//cout << "p_k : " << pk << endl;
			//cout << "Refused real for " << i-1 << ": " << Nrefused_i << " but should be: " << Nr << "\n";
			//cout << "Accepted real for " << i-1 << ": " << Naccepted_i << " but should be: " << Na << "\n";
			cout << "p_" << i-1 << ": " << p_i; 
			//p_i = 0.5 + sqrt((Na -  Nr)/( Nr + Na))/2;
			cout <<" but should be: " << p_i << "\n";
			//continue;
		}*/
		float dp = sqrt(Naccepted_i * Nrefused_i * Nrefused_i + Nrefused_i * Naccepted_i * Naccepted_i) / (2*p_i-1) / pow(Nrefused_i + Naccepted_i,2);
		averagedP += dp * Naccepted_i;
		cout << "P_"<<i-1<<": " << p_i << " dP_i: " << dp << endl;
		averageP += p_i * Naccepted_i;
		nevents += Naccepted_i;
	}
	averageP /= nevents;
	error = averagedP / nevents;
	cout << "AVERAGE P: " << averageP << endl;
	return averageP;
}
TH1F * computeCorrectionF(TH1F * right, float p, float dp, int bin_e = 20, string name = "", TH1F * ref = NULL)
{
	string hname = "correctedP" + name;
	TH1F * corrected = new TH1F(hname.c_str(),";cos",bin_e,-1, 1);
	int niterations = (float) bin_e/2;
	right->Draw();
	int nref = 0;
	for (unsigned int i = 1; i < niterations+1; i++) 
	{
		int j = bin_e + 1 - i;
		float naminus_i = right->GetBinContent(i);
		float naplus_i = right->GetBinContent(j);
		float p_i = p;
		float q_i = 1. - p_i;
		float nplus_i = Nplus(p_i,q_i, naplus_i, naminus_i) ;
		float nminus_i = Nminus(p_i,q_i, naplus_i, naminus_i);
		/*if (ref) 
		{
			nref +=  ref->GetBinContent(j) + ref->GetBinContent(i);
			cout << "\tN_total: " << naminus_i + naplus_i + ref->GetBinContent(j) + ref->GetBinContent(i) << " vs " << nplus_i + nminus_i << endl;
			float diff = nminus_i + nplus_i - (naminus_i + naplus_i + ref->GetBinContent(j) + ref->GetBinContent(i));
			nplus_i = nplus_i -  diff / 2;
			nminus_i = nminus_i -  diff / 2;
		}*/
		float n2plus_i = p_i*p_i * nplus_i ;// + q_i*q_i * nminus_i;
		float n2minus_i = p_i*p_i * nminus_i;//  + q_i*q_i * nplus_i;
		//float error_minus = sigma(naminus_i, naplus_i, p_i,q_i, 1, 1, dp);
		//float error_plus = sigma(naplus_i, naminus_i, p_i,q_i, 1, 1, dp);
		float error_minus = sqrt(n2minus_i + 4*nminus_i*nminus_i*p*p*dp*dp);
		float error_plus =  sqrt(n2plus_i + 4*nplus_i*nplus_i*p*p*dp*dp);
		corrected->SetBinContent(i, n2minus_i);
		corrected->SetBinContent(j, n2plus_i);
		corrected->SetBinError(i,error_minus);
		corrected->SetBinError(j, error_plus);
//		cout << "Corr: " << n2minus_i << endl;
	}
	cout << "\tN_ref: " << nref << endl;
	return corrected;
}
TH1F * computeCorrectionP(TH1F * right, TH1F * p, int bin_e = 20, string name = "")
{
	string hname = "correctedP" + name;
	TH1F * corrected = new TH1F(hname.c_str(),";cos",bin_e,-1, 1);
	int niterations = (float) bin_e/2;
	for (unsigned int i = 1; i < niterations+1; i++) 
	{
		int j = bin_e + 1 - i;
		float naminus_i = right->GetBinContent(i);
		float naplus_i = right->GetBinContent(j);
		float p_i = p->GetBinContent(i);
		float q_i = 1-p_i;
		float sp_i = p->GetBinError(i);
		float q_i = 1. - p_i;
		float nplus_i = Nplus(p_i,q_i, naplus_i, naminus_i) ;
		float nminus_i = Nminus(p_i,q_i, naplus_i, naminus_i);
		float n2plus_i = p_i*p_i * nplus_i ;
		float n2minus_i =  p_i*p_i*nminus_i ;
		//float error_minus = sqrt( nminus_i + (pow(q_i*q_i/p_i/p_i, 2) +pow(2*(2/p_i/p_i - 1)*q_i/p_i, 2) *sp_i )* nplus_i );
		//float error_plus = sqrt( n2minus_i );
		//float error_minus = sqrt( 4*p_i*p_i * sp_i * sp_i * nminus_i * nminus_i + p_i*p_i * p_i*p_i * nminus_i);
		//float error_plus = sqrt( 4*p_i*p_i * sp_i * sp_i * nplus_i * nplus_i + p_i*p_i * p_i*p_i * nplus_i);
		float error_minus = sqrt( naminus_i +  4*q_i*q_i * sp_i * sp_i * nplus_i * nplus_i);
		float error_plus = sqrt( naplus_i +  4*q_i*q_i * sp_i * sp_i * nminus_i * nminus_i);

		//float error_minus = sqrt( nminus_i * nminus_i * 4*p_i*p_i * sp_i * sp_i + pow(p_i*p_i *sigma(naminus_i, naplus_i, p_i,q_i, 1, 1),2)) ;
		//float error_plus = sqrt( nplus_i * nplus_i * 4*p_i*p_i * sp_i * sp_i + pow(p_i*p_i *sigma(naplus_i, naminus_i, p_i,q_i, 1, 1),2) );
		corrected->SetBinContent(i, n2minus_i);
		corrected->SetBinContent(j, n2plus_i);
		corrected->SetBinError(i,error_minus);
		corrected->SetBinError(j, error_plus);
	}
	return corrected;
}
TH1F * computeCorrection(TH1F * right, TH1F * wrong, int bin_e = 20, string n = "",  TH1F * pglobal = NULL)
{
	if (right->GetEntries() <  wrong->GetEntries()) 
	{
		std::cout << "ERROR " << right->GetEntries() << " " << wrong->GetEntries() << "!\n";
		return right;
	}
	TH1F * corrected = new TH1F(string("corrected"+n).c_str(),";cos",bin_e,-1, 1);
	//float ntminus = 0.0;
	//float ntplus = 0.0;
	float naminus = 0.0;
	float naplus = 0.0;
	float nplus = 0.0;
	float nminus = 0.0;
	int niterations = (float) bin_e/2;
	for (unsigned int i = 1; i < niterations+1; i++) 
	{
		int j = bin_e + 1 - i;
		float Naccepted_i =( right->GetBinContent(i) + right->GetBinContent(j));
		float Nrefused_i = (wrong->GetBinContent(i) + wrong->GetBinContent(j));
		//float p_i = 0.5 + sqrt(1-2*Nrefused_i/(Nrefused_i + Naccepted_i))/2;
		float p_i = 0.5 + sqrt((Naccepted_i - Nrefused_i)/(Nrefused_i + Naccepted_i))/2;
		float q_i = 1. - p_i;
		std::cout << "Local: p_i = "<< p_i << "; q_i = " << q_i << endl;
		float naminus_i = right->GetBinContent(i);
		float naplus_i = right->GetBinContent(j);
		float nplus_i = Nplus(p_i,q_i, naplus_i, naminus_i) ;
		float nminus_i = Nminus(p_i,q_i, naplus_i, naminus_i);
		float n2plus_i = p_i*p_i * nplus_i ;//+ q_i*q_i * nminus_i;
		float n2minus_i =  p_i*p_i*nminus_i ;//+  q_i*q_i *nplus_i;

		float error_minus = sigma(naminus_i, naplus_i, p_i,q_i, Nrefused_i, Naccepted_i) ;
		float error_plus =  sigma(naplus_i, naminus_i, p_i,q_i, Nrefused_i, Naccepted_i) ;
		std::cout << "Measured: Na+ = " << naplus_i << "; Na- = " << naminus_i <<endl;
		std::cout << "Refused: Nr = " << Nrefused_i <<endl;
		std::cout << "True: N+ = " << nplus_i << "; N- = " << nminus_i <<endl;
		std::cout << "Corrected: N+ = " << n2plus_i << "; N- = " << n2minus_i << "; sN- = " << error_minus <<endl << endl;
		corrected->SetBinContent(i, n2minus_i);
		corrected->SetBinError(i,error_minus);

		corrected->SetBinContent(j, n2plus_i);
		corrected->SetBinError(j, error_plus);
		if (pglobal) 
		{
			pglobal->SetBinError(i, sigmaP(Nrefused_i, Naccepted_i));
			pglobal->SetBinContent(i, p_i);
			pglobal->SetBinError(j, sigmaP(Nrefused_i, Naccepted_i));
			pglobal->SetBinContent(j, p_i);
		}
		naminus += naminus_i;
		naplus += naplus_i;
		nplus += n2plus_i;
		nminus += n2minus_i;
	}
	std::cout << "-------------------------\n";
	float afbmeasured =  (naplus-naminus)/(naplus+naminus);
	float afbcorrected =  (nplus-nminus)/(nplus+nminus);
	return corrected;
}
float sigmaP(float nr, float n)
{
	float dpdNa = (1-2*nr/n < .0001)? 0: nr / (2* n*n *sqrt(1-2*nr/n));
	float dpdNr = (1-2*nr/n < .0001)? 0: 1. / (2* n *sqrt(1-2*nr/n));
	float sigmap = sqrt(dpdNa*dpdNa*n  + dpdNr * dpdNr * nr);
	return sigmap;
}
float sigma(float naplus, float naminus, float p, float q, float nr, float n, float dp = 0.)
{
	float alpha = abs(1/(p*p*p*p-q*q*q*q));
	float sigmap = dp;//sigmaP(nr,n);//sqrt(dpdNa*dpdNa*n  + dpdNr * dpdNr * nr);
	return alpha * sqrt( naplus * p*p*p*p + naminus * q*q*q*q + sigmap * sigmap * pow(2*p*naplus + 2 *q* naminus ,2));
	return sqrt(naplus); //+ 2*naplus*naplus*p*p*dp*dp);
}
float equation()
{
	cout << "**********************************\n";
	float N0pm = 0.493;
	float Npm = 0.193;
	float N00 = 0.313;
	float Na = 0.223;
	float a = (N00 - N0pm/2)/2;
	float b = N0pm/4 -N00;
	float c = Npm +N00/2 - Na;

	float d = b*b - 4 * a * c;
	float p1  = 1./2/a * (-b - sqrt(d));
	float p2 = 1./2/a * (-b + sqrt(d));
	cout << "p_v1 = " << p1 << "; p_v2 = " << p2 << endl;
	cout << "0 = " << a*p1*p1 + b* p1+c << endl;
	float Nr = N0pm / 2 * (1 - p1) *p1 + (1-p1)*(1-p1)*N00 / 2+ Npm * p1 * (1-p1)*(1-p1)/2;
	cout << "Nr = " << Nr << endl;
	return p1;
}

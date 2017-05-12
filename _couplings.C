#include "TMath.h"
void _couplings()
{
	float pi = TMath::Pi();
	float range = 280;
	TCanvas * c1 = new TCanvas("c1", "Data-MC",0,0,500,500);
	float s = 91.2;
	float s2 = 250;
	float Mz = 91.1876;
	float Gz = 2.5;
	float sqrtrb = sqrt(0.9868);
	float sqrtrl = sqrt(0.998);
	float kb = 1.0065;
	float kl = 1.001;
	float alpha = 1./128;
	float sw = 0.2329;
	float swe = 0.2315;
	float cw = 1 - sw;
	float cwe = 1 - swe;
	float factor = 1;
	float All = (-0.5+swe)*(-0.5+sw/3)/sqrt(sw*cw)/sqrt(swe*cwe); //(-0.5+sw)*(-0.5+sw/3)/sw/cw; // sqrt(sw*cw)/sqrt(swe*cwe)
	float Alr = (-0.5+swe)*(sw/3)/sqrt(sw*cw)/sqrt(swe*cwe); //(-0.5+sw)*(sw/3*factor)/cw/sw;
	float Arl = (swe)*(-0.5+sw/3)/sqrt(sw*cw)/sqrt(swe*cwe); //(sw)*(-0.5+sw/3)/sw/cw;
	float Arr = (swe)*(sw/3)/sqrt(sw*cw)/sqrt(swe*cwe); //(sw)*(sw/3*factor)/sw/cw;
	
	//float All = sqrtrb * sqrtrl *(-0.5+kl*sw)*(-0.5+kb*sw/3)/sqrt(kl * sw*(1-kl*sw))/sqrt(kb * sw*(1-kb*sw)); //(-0.5+sw)*(-0.5+sw/3)/sw/cw; // sqrt(sw*cw)/sqrt(swe*cwe)
	//float Alr = sqrtrb * sqrtrl *(-0.5+kl*sw)*(kb*sw/3)/sqrt(kl * sw*(1-kl*sw))/sqrt(kb * sw*(1-kb*sw)); //(-0.5+sw)*(sw/3*factor)/cw/sw;
	//float Arl = sqrtrb *sqrtrl *  (kl*sw)*(-0.5+kb*sw/3)/sqrt(kl * sw*(1-kl*sw))/sqrt(kb * sw*(1-kb*sw)); //(sw)*(-0.5+sw/3)/sw/cw;
	//float Arr = sqrtrb * sqrtrl * (kl*sw)*(kb*sw/3)/sqrt(kl * sw*(1-kl*sw))/sqrt(kb * sw*(1-kb*sw)); //(sw)*(sw/3*factor)/sw/cw;
	
	float Afbl = 3./4*(Q(All, s2) - Q(Alr, s2))/(Q(All, s2) + Q(Alr, s2));
	float Afbr = 3./4*(Q(Arr, s2) - Q(Arl, s2))/(Q(Arl, s2) + Q(Arr, s2));
	float Afb = 3./4*(Q(Arr, s) + Q(All, s) - Q(Arl, s) -  Q(Alr, s))/(Q(Arl, s) + Q(Alr, s) + Q(All, s) + Q(Arr, s));

	//float xsection = 3*4*TMath::Pi()/3*alpha * alpha / s / s * (Q(Arl, s) + Q(Alr, s) + Q(All, s) + Q(Arr, s)) ;
	float xsection = 3 / s / s * (Q(Arl, s) + Q(Alr, s) + Q(All, s) + Q(Arr, s)) *  pi * alpha * alpha /3;
	float xsectionl = 3 / s2 / s2 * ( Q(Alr, s2) + Q(All, s2)) * 4 * pi * alpha * alpha /3;
	float xsectionr = 3 / s2 / s2 * ( Q(Arl, s2) + Q(Arr, s2)) * 4 * pi * alpha * alpha /3;
	//float dx = xsection*0.0035;
	float dx = xsection*0.002;
	float min = 20;
	TF1 *fAfb2 = new TF1("fa2","3./4*(Q([0], x) - Q([1], x))/(Q([0], x) + Q([1], x))",min,range);
	fAfb2->SetLineColor(kGreen);
	fAfb2->SetParameters(All,Alr);
	fAfb2->SetLineStyle(5);
	fAfb2->SetTitle(";#sqrt{s} [GeV];A_{FB}");
	fAfb2->Draw();
	TF1 *fAfb = new TF1("fa1","3./4*(Q([0], x) - Q([1], x))/(Q([0], x) + Q([1], x))",min,range);
	fAfb->SetParameters(Arr,Arl);
	fAfb->SetLineColor(kBlue);
	fAfb->SetLineStyle(2);
	fAfb->Draw("same");
	TF1 *fAfb3 = new TF1("fa3","3./4*(Q([0], x) + Q([1], x)  - Q([2], x) -  Q([3], x))/(Q([0], x) + Q([1], x) +  Q([2], x) + Q([3], x))",min,range);
	fAfb3->SetParameters(All,Arr, Alr, Arl);
	fAfb3->Draw("same");
	int n = 2;
	float x[2] = {250, 250};
	float y[2] = {0.699, 0.28};
	float sx[2] = {0.1, 0.1};
	float sy[2] = {0.0028, 0.01};
	float x2[2] = {91.26, 190.7};
	float y2[2] = {0.098, 0.51};
	float sy2[2] = {0.0017, 0.058};
	float x3[4] = {29, 35, 44, 58};
	float y3[4] = {-0.052, -0.214, -0.46, -0.58};
	float sy3[4] = {0.081, 0.05, 0.147, 0.078};
	float sx3[4] = {0.1, 0.1, 0.1, 0.1};
	TGraphErrors * results = new  TGraphErrors(n,x,y,sx,sy);
	results->SetLineColor(kRed);
	results->SetLineWidth(3);
	results->SetMarkerStyle(21);
	results->SetMarkerColor(kRed);
	TGraphErrors * results2 = new  TGraphErrors(n,x2,y2,sx,sy2);
	results2->SetLineWidth(3);
	results2->SetMarkerStyle(20);
	results2->SetMarkerColor(kMagenta);
	results2->SetLineColor(kMagenta);
	TGraphErrors * results3 = new  TGraphErrors(4,x3,y3,sx3,sy3);
	results3->SetLineWidth(3);
	results3->SetMarkerStyle(22);
	results3->SetLineColor(kGray+1);
	results3->SetMarkerColor(kGray+1);
	results->Draw("samep");
	results2->Draw("samep");
	results3->Draw("samep");
	
	TLegend *legendMean2 = new TLegend(0.65,0.2,0.85,0.4,NULL,"brNDC");
	legendMean2->SetFillColor(kWhite);
	legendMean2->SetBorderSize(0);
	legendMean2->AddEntry(results,"ILC 250 GeV","p");
	legendMean2->AddEntry(results2,"LEP","p");
	legendMean2->AddEntry(results3,"Others","p");
	legendMean2->AddEntry(fAfb3,"e^{-}e^{+}","l");
	legendMean2->AddEntry(fAfb2,"e_{L}^{-}e_{R}^{+}","l");
	legendMean2->AddEntry(fAfb,"e_{R}^{-}e_{L}^{+}","l");
	legendMean2->Draw();

	TCanvas * c2 = new TCanvas("c2", "Data-MC",500,0,500,500);
	TF1 *fsigmal = new TF1("fs","3 / x / x * (Q([0], x) + Q([1], x)) * 4 * 3.14 * [2] * [2] /3 * 1e15 * 0.3894e-3 ",min,range);
	fsigmal->SetParameters(Alr,All,alpha);
	fsigmal->SetLineColor(kGreen);
	fsigmal->SetLineStyle(5);
	TF1 *fsigmar = new TF1("fs1","3 / x / x * (Q([0], x) + Q([1], x)) * 4 * 3.14 * [2] * [2] /3 * 1e15 * 0.3894e-3 ",min,range);
	fsigmar->SetParameters(Arr,Arl,alpha);
	fsigmar->SetLineColor(kBlue);
	fsigmar->SetLineStyle(2);

	TF1 *fsigma2 = new TF1("fs2","3 / x / x * (Q([0], x) + Q([1], x) + Q([2], x) + Q([3], x)) *  3.14 * [4] * [4] /3 * 1e15 * 0.3894e-3 ",min,range);
	fsigma2->SetParameters(Arr,Arl, Alr,All,alpha);
	//fsigma2->SetLineStyle(2);
	fsigma2->SetMinimum(1e3);
	fsigma2->SetTitle(";#sqrt{s} [GeV];#sigma_{I} [fb]");
	float x4[1] = {91.18};
	float sx4[1] = {0.1};
	float y4[1] = {fsigma2->Eval(91.18)};
	float sy4[1] = {y4[1] * 0.003};
	float x5[2] = {250,250};
	float sx2[2] = {0.1, 0.1};
	float y5[2] = {fsigmal->Eval(x5[0]),fsigmar->Eval(x5[0])};
	float sy5[2] = {fsigmal->Eval(x5[0])*0.005,fsigmar->Eval(x5[0])*0.01};
	TGraphErrors * results4 = new  TGraphErrors(1,x4,y4,sx4,sy4);
	results4->SetLineWidth(3);
	results4->SetMarkerStyle(20);
	results4->SetLineColor(kMagenta);
	results4->SetMarkerColor(kMagenta);
	TGraphErrors * results5 = new  TGraphErrors(2,x5,y5,sx2,sy5);
	results5->SetLineColor(kRed);
	results5->SetLineWidth(3);
	results5->SetMarkerStyle(21);
	results5->SetMarkerColor(kRed);
	
	/*TF1 * polar = new TF1("fa","1/[2]/[2]*((Q([0], [2]) - Q([1], [2]))*2*x + (1+x*x)*(Q([0], [2]) + Q([1], [2])) )",-1,1);
	polar->SetLineColor(kBlue);
	polar->SetLineStyle(2);
	polar->SetParameters(Arr,Arl,s2);
	TF1 * polarStd = new TF1("fb","1/[2]/[2]*((Q([0], [2]) - Q([1], [2]))*2*x + (1+x*x)*(Q([0], [2]) + Q([1], [2])))",-1,1);
	polarStd->SetTitle("Differential cross section;cos#theta;#frac{d#sigma}{dcos#theta}");
	polarStd->SetMinimum(0);
	polarStd->SetParameters((sw)*(sw/3)/sw/cw, Arl ,s2);
	polarStd->SetLineColor(kBlue+1);*/
	//polarStd->Draw();
	//polar->Draw("same");
	gPad->SetLogy();
	fsigma2->Draw();
	fsigmal->Draw("same");
	fsigmar->Draw("same");
	results4->Draw("samep");
	results5->Draw("samep");
	TLegend *legendMean2 = new TLegend(0.65,0.52,0.85,0.7,NULL,"brNDC");
	legendMean2->SetFillColor(kWhite);
	legendMean2->SetBorderSize(0);
	legendMean2->AddEntry(results,"ILC 250 GeV","p");
	legendMean2->AddEntry(results2,"LEP","p");
	legendMean2->AddEntry(fAfb3,"e^{-}e^{+}","l");
	legendMean2->AddEntry(fAfb2,"e_{L}^{-}e_{R}^{+}","l");
	legendMean2->AddEntry(fAfb,"e_{R}^{-}e_{L}^{+}","l");
	legendMean2->Draw();
	cout << "Afb = " << fAfb3->Eval(s)
	     << "; dAfb = " << (fAfb3->Eval(s) - y2[0]) / sy2[0]
	     << "; Afbl = "   << Afbl 
	     << "; Afbr = " << Afbr
	     << endl;
	cout << "s = " << xsection * 1e15 * 0.3894e-3
	     << " +- " << dx * 1e15 * 0.3894e-3
	     << "; sl = " << xsectionl * 1e15 * 0.3894e-3
	     << "; sr = " << xsectionr * 1e15 * 0.3894e-3
	     << endl;

	float xmin = xsection - dx;
	float xmax = xsection + dx; 
	
	float amin = y2[0] - sy2[0];
	float amax = y2[0] + sy2[0];
	LEPConstrains(xmin, xmax, amin, amax);
	float ddx = 0.002;
	float dda = 0.01;
	float xminl = xsectionl *(1-ddx);
	float xmaxl = xsectionl *(1+ddx); 

	float aminl = y[0] - sy[0];
	float amaxl = y[0] + sy[0];
	cout << "Afbl low = " << aminl << "; Afbl high = " << amaxl << "\n";
	float xminr = xsectionr *(1-ddx*2.5);
	float xmaxr = xsectionr *(1+ddx*2.5); 

	float aminr = y[1] - sy[1];
	float amaxr = y[1] + sy[1];
	cout << "Afbr low = " << aminr << "; Afbr high = " << amaxr << "\n";
	ILCConstrains(xminl, xmaxl, xminr, xmaxr, aminl, amaxl, aminr, amaxr);//*/
}
void LEPConstrains(float xmin, float xmax, float amin, float amax)
{
	float pi = TMath::Pi();
	float alpha = 1./128;
	float s = 91.26;
	float sw = 0.2319;
	float cw = 1 - sw;
	TCanvas * c3 = new TCanvas("c3", "Data-MC",0,500,500,500);
	int tries = 1000;
	float high = 1;
	float highl = 0.2;
	TH2F * good = new TH2F("g00d", "; factor; factor", tries, -highl, highl, tries, -high, high);
	TH2F * goodA = new TH2F("g00dA", "; factor; factor", tries, -highl, highl, tries, -high, high);
	TH2F * bad = new TH2F("b0d", "; #delta g_{L}^{Z}/g_{L}^{Z}; #delta g_{R}^{Z}/g_{R}^{Z}", tries, -highl, highl, tries, -high, high);
	bad->SetStats(0);
	cout << "g_L = " << (-0.5+sw/3) / (4*sqrt(sw*cw)) << "; g_R = " << (sw/3) / (4*sqrt(sw*cw)) << "\n";
	float factorrmin = 10;
	float factorrmax = 0;
	for (unsigned int i = 0; i < tries; i++) 
	{
		for (unsigned int j = 0; j < tries; j++) 
		{
			float factorl = -highl + i * 2*highl/tries;
			float factorr = -high + j * 2*high/tries;
			float All2 = (-0.5+sw)*(-0.5+sw/3)/sw/cw * (1+factorl);

			float Arl2 = (sw)*(-0.5+sw/3)/sw/cw * (1+factorl);
			float Alr2 =  (-0.5+sw)*(sw/3*(1+factorr))/sw/cw;

			float Arr2 = (sw)*(sw/3*(1+factorr))/sw/cw;
			float xsection = 3 / s / s * (Q(Arl2, s) + Q(Alr2, s) + Q(All2, s) + Q(Arr2, s)) *  pi * alpha * alpha /3;
			float Afb = 3./4*(Q(Arr2, s) + Q(All2, s) - Q(Arl2, s) -  Q(Alr2, s))/(Q(Arl2, s) + Q(Alr2, s) + Q(All2, s) + Q(Arr2, s));
			int goodTries = 0;
			if (xsection > xmin && xsection < xmax) 
			{
				good->Fill(factorl, factorr);
				goodTries++;
				//cout << "Factorr = " << factorr ;
				//cout << " is ok\n";
			}
			if (Afb > amin && Afb < amax) 
			{
				//cout << " is ok\n";
				goodTries++;

				goodA->Fill(factorl, factorr);
			}
			if (goodTries == 2 && factorl > -0.5 && factorr > -0.5) 
			{
				if (factorr < factorrmin) 
				{
					factorrmin = factorr;
				}
				if (factorr > factorrmax) 
				{
					factorrmax = factorr;
				}
				
			}
		}
		
	}
	makePretty(good, kBlue);
	makePretty(goodA, kRed);
	
	cout << "Min: " << factorrmin << " Max: " << factorrmax << " diff: " << factorrmax - factorrmin << "\n";
	bad->Draw("p");
	goodA->Draw("samep");
	good->Draw("samep");
	gPad->SetLeftMargin(0.14);
	gPad->SetBottomMargin(0.14);
	bad->GetYaxis()->SetTitleOffset(1.3);
	TLegend *legendMean = new TLegend(0.65,0.2,0.85,0.3,NULL,"brNDC");
	legendMean->SetFillColor(kWhite);
	legendMean->SetBorderSize(0);
	legendMean->AddEntry(good,"#sigma","f");
	legendMean->AddEntry(goodA,"A_{FB}","f");
	legendMean->Draw();
}
void ILCConstrains(float xminl, float xmaxl, float xminr, float xmaxr, float aminl, float amaxl, float aminr, float amaxr)
{
	float pi = TMath::Pi();
	float alpha = 1./128;
	float s = 250;
	float sw = 0.2319;
	float cw = 1 - sw;
	TCanvas * c4 = new TCanvas("c4", "Data-MC",500,500,500,500);
	int tries = 1000;
	float high = 1;

	float highl = 0.2;
	TH2F * goodXl = new TH2F("g00dXl", "; factor; factor", tries, -highl, highl, tries, -high, high);
	TH2F * goodXr = new TH2F("g00dXr", "; factor; factor", tries, -highl, highl, tries, -high, high);
	TH2F * goodAl = new TH2F("g00dAl", "; factor; factor", tries, -highl, highl, tries, -high, high);
	TH2F * goodAr = new TH2F("g00dAr", "; factor; factor", tries, -highl, highl, tries, -high, high);
	TH2F * bad = new TH2F("b0d", "; #delta g_{L}^{Z}/g_{L}^{Z};  #delta g_{R}^{Z}/g_{R}^{Z}", tries, -highl, highl, tries, -high, high);
	bad->SetStats(0);
	cout << "g_L = " << (-0.5+sw/3) / (4*sqrt(sw*cw)) << "; g_R = " << (sw/3) / (4*sqrt(sw*cw)) << "\n";
	float factorrmin = 10;
	float factorrmax = 0;
	float g_l = (-0.5+sw/3)/sqrt(sw*cw);
	float g_r = (sw/3)/sqrt(sw*cw);
	for (unsigned int i = 0; i < tries; i++) 
	{
		for (unsigned int j = 0; j < tries; j++) 
		{
			float factorl = -highl + i * 2*highl/tries;
			float factorr = -high + j * 2*high/tries;
			/*float F_a = 1./4/sqrt(sw)/sqrt(cw) * (1+factora);
			float F_v = (-1 + 4./3 * sw)/4/sqrt(sw)/sqrt(cw) * (1+factorv);

			float All2 = (-0.5+sw)*(F_v - F_a)/sqrt(sw)/sqrt(cw) ;
			float Arl2 = (sw)*(F_v - F_a)/sqrt(sw)/sqrt(cw) ;
			float Alr2 =  (-0.5+sw)*(F_v + F_a)/sqrt(sw)/sqrt(cw);
			float Arr2 = (sw)*(F_v + F_a)/sqrt(sw)/sqrt(cw);

			float factorl = (g_l - (F_v - F_a)) / g_l;
			float factorr = (g_r - (F_v + F_a)) / g_r;*/

			float All2 = (-0.5+sw)*(-0.5+sw/3)/sw/cw * (1+factorl);
			float Arl2 = (sw)*(-0.5+sw/3)/sw/cw * (1+factorl);
			float Alr2 =  (-0.5+sw)*(sw/3*(1+factorr))/sw/cw;
			float Arr2 = (sw)*(sw/3*(1+factorr))/sw/cw;
			float xsectionl = 3 / s / s * ( Q(Alr2, s) + Q(All2, s)) * 4 * pi * alpha * alpha /3;
			float xsectionr = 3 / s / s * ( Q(Arl2, s) + Q(Arr2, s)) * 4 * pi * alpha * alpha /3;
			float Afbl = 3./4*(Q(All2, s) - Q(Alr2, s))/(Q(All2, s) + Q(Alr2, s));
			float Afbr = 3./4*(Q(Arr2, s) - Q(Arl2, s))/(Q(Arl2, s) + Q(Arr2, s));
			int goodEntries = 0;
			if (xsectionl > xminl && xsectionl < xmaxl) 
			{
				goodXl->Fill(factorl, factorr);
				//goodEntries++;
			}
			if (Afbl > aminl && Afbl < amaxl) 
			{
				goodAl->Fill(factorl, factorr);
				goodEntries++;
			}
			if (xsectionr > xminr && xsectionr < xmaxr) 
			{
				goodXr->Fill(factorl, factorr);
				//goodEntries++;
			}
			if (Afbr > aminr && Afbr < amaxr) 
			{
				goodAr->Fill(factorl, factorr);
				goodEntries++;
			}
			if (goodEntries == 2) 
			{
				if (factorr < factorrmin) 
				{
					factorrmin = factorr;
				}
				if (factorr > factorrmax) 
				{
					factorrmax = factorr;
				}
			}
		}
		
	}
	makePretty(goodXl, kBlue+1);
	makePretty(goodXr, kBlue-1);
	makePretty(goodAr, kRed-1);
	makePretty(goodAl, kRed+1);
	cout << "Min: " << factorrmin << " Max: " << factorrmax << " diff: " << factorrmax - factorrmin << "\n";
	bad->Draw("p");
	goodAl->Draw("samep");
	goodAr->Draw("samep");
	goodXr->Draw("samep");
	goodXl->Draw("samep");
	gPad->SetLeftMargin(0.14);
	gPad->SetBottomMargin(0.14);
	bad->GetYaxis()->SetTitleOffset(1.3);
	TLegend *legendMean = new TLegend(0.65,0.2,0.85,0.4,NULL,"brNDC");
	legendMean->SetFillColor(kWhite);
	legendMean->SetBorderSize(0);
	legendMean->AddEntry(goodXl,"#sigma^{L}","f");
	legendMean->AddEntry(goodXr,"#sigma^{R}","f");
	legendMean->AddEntry(goodAl,"A^{L}_{FB}","f");
	legendMean->AddEntry(goodAr,"A^{R}_{FB}","f");
	legendMean->Draw();
}
void makePretty(TH2F * good, int color)
{
	good->SetMarkerColor(color);
	good->SetFillColor(color);
	good->SetMarkerSize(0.5);
	good->SetMarkerStyle(20);
}

float Q(float A, float sqrts = 250)
{
	float s = sqrts*sqrts;
	float Mz = 91.1876;
	float Gz = 2.5;
	return 1./9 + (2./3 * A * s * (s - Mz*Mz) + A*A*s*s)/( (s - Mz*Mz)*(s - Mz*Mz) + Gz*Gz*Mz*Mz );
}

#include <math.h>
#include "/users/flc/bilokin/Processors/Macros/top/basymmetry.C"
#include <iostream>
#include <iomanip>
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"
#include "TMath.h"
#include "TFile.h"
#include "TCut.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLorentzVector.h"
using namespace std;
void _formfactors()
{
	
	TCanvas * c1 = new TCanvas("c1", "Data-MC",0,0,500,500);
	float sqrts = 250;
	float s = sqrts * sqrts;
	float pi = TMath::Pi();
	float Mz = 91.1876;
	float Gz = 2.5;
	float alpha = 1./128;
	//float sw = 0.2329;
	float sw = 0.23116;
	float cw = 1 - sw;
	float A = 4*pi*alpha*alpha/3/s;

	float Qf = -1./3; float If = -0.5; float mf = 4.18;
	//float Qf = 2./3; float If = 0.5; float mf = 175;
	float gamma = sqrts / 2 / mf;
	float beta =  2 * sqrt( s/4 - mf*mf) / sqrts;
	float F1Vg = Qf;
	float F1Ag = 0;

	float F1Vz = (If - 2*Qf * sw)/sqrt(4*cw*sw) ;
	float F1Az = - If / sqrt(4*cw*sw);
	
	float FeL = (-0.5 + sw)/sqrt(cw*sw);
	float FeR =  sw/sqrt(cw*sw);

	float RePropagator = s*(s-Mz*Mz) / ((s-Mz*Mz)*(s-Mz*Mz) + Gz*Mz*Gz*Mz);
	float ImPropagator = - s*Gz*Mz / ((s-Mz*Mz)*(s-Mz*Mz) + Gz*Mz*Gz*Mz);

	float mF1V2L = pow(- F1Vg + FeL * F1Vz * RePropagator,2)  + pow(FeL * F1Vz * ImPropagator,2);
	float mF1A2L = pow(FeL * F1Az * RePropagator,2)  + pow(FeL * F1Az * ImPropagator,2);

	float mF1VmF1AL = -F1Vg * FeL * F1Az * RePropagator +  FeL * F1Az  *  FeL * F1Vz *(RePropagator * RePropagator + ImPropagator * ImPropagator);
	
	float Al = mF1V2L + beta*beta*mF1A2L;
	float Bl = -4*beta*mF1VmF1AL;
	float Cl = mF1V2L / gamma / gamma;

	float mF1V2R = pow(- F1Vg + FeR * F1Vz * RePropagator,2)  + pow(FeR * F1Vz * ImPropagator,2);
	float mF1A2R = pow(FeR * F1Az * RePropagator,2)  + pow(FeR * F1Az * ImPropagator,2);

	float mF1VmF1AR = F1Vg * FeR * F1Az * RePropagator -  FeR * F1Az  *  FeR * F1Vz *(RePropagator * RePropagator + ImPropagator * ImPropagator);

	float Ar = mF1V2R +beta*beta* mF1A2R;
	float Br = -4*beta*mF1VmF1AR;
	float Cr = mF1V2R / gamma / gamma;
	
	float N = A*3*3/4*beta*  1e15 * 0.3894e-3 ;
	
	cout <<"----------------FF----------------\n";
	//cout << "RePropagator: " << RePropagator << endl;
	//cout << "ImPropagator: " << ImPropagator << endl;
	cout << "Al: " << Al*N << " Bl: " << Bl*N << " Cl: " << Cl*N << endl;
	cout << "Ar: " << Ar*N << " Br: " << Br*N << " Cr: " << Cr*N << endl;

	TF1 * dsigmaL = new TF1("dsigmaL","([0] *(1+x*x) + [1]*x + [2]*(1-x*x))*[3]",-1,1);
	//TF1 * dsigmaL = new TF1("dsigmaL","dsigma(x,[0],1)",-1,1);
	dsigmaL->SetParameters(Al,Bl,Cl,N);
	//dsigmaL->SetParameter(0,250);
	TF1 * dsigmaR = new TF1("dsigmaR","([0] *(1+x*x) + [1]*x + [2]*(1-x*x))*[3]",-1,1);
	//TF1 * dsigmaR = new TF1("dsigmaR","dsigma(x,[0],0)",-1,1);
	//dsigmaR->SetParameter(0,250);
	dsigmaR->SetParameters(Ar,Br,Cr,N);
	dsigmaL->Draw("");
	dsigmaR->Draw("same");
	float afbl = (dsigmaL->Integral(0,1) - dsigmaL->Integral(-1,0)) / (dsigmaL->Integral(-1,1));
	float afbr = (dsigmaR->Integral(0,1) - dsigmaR->Integral(-1,0)) / (dsigmaR->Integral(-1,1));
	cout << "AfbL: " << afbl << " Afbr: " << afbr << endl;
	cout << "sigmaTL: " << dsigmaL->Integral(-1,1) << " sigmaTr: " << dsigmaR->Integral(-1,1) << endl;
	cout << "Events 250 fb-1 left:  " <<  dsigmaL->Integral(-1,1)*250 << " right: " <<  dsigmaR->Integral(-1,1)*250  << endl;
	float lumiLeft = 500 * 0.675;
	float lumiRight = 500 * 0.225;
	float Pe = 0.8;
	float Pp = 0.3;
	float effectiveSigmaL = 0.25*( (1+Pe*Pp)*(dsigmaL->Integral(-1,1)+dsigmaR->Integral(-1,1)) + (-Pe-Pp)*(dsigmaR->Integral(-1,1)-dsigmaL->Integral(-1,1)) );
	float effectiveSigmaR = 0.25*( (1+Pe*Pp)*(dsigmaL->Integral(-1,1)+dsigmaR->Integral(-1,1)) + (Pe+Pp)*(dsigmaR->Integral(-1,1)-dsigmaL->Integral(-1,1)) );
	cout << "Effective sigmaTL: " << effectiveSigmaL << " sigmaTr: " <<effectiveSigmaR << endl;
	cout << "Events H20 1 run left:  " <<  effectiveSigmaL*lumiLeft << " right: " <<  effectiveSigmaR*lumiRight  << endl;
	cout << "Error factor left: " << sqrt(dsigmaL->Integral(-1,1)*250) / sqrt(effectiveSigmaL*lumiLeft)<< " right: " <<  sqrt(dsigmaR->Integral(-1,1)*250) / sqrt(effectiveSigmaR*lumiRight) << endl;
	float selectionLeft = 0.4;
	float selectionRight = 0.33;
	float lumiError = 0.001;
	cout << "Error sigma left: "<< sqrt(1/(effectiveSigmaL*lumiLeft*selectionLeft) +lumiError*lumiError) <<  " right: " <<  sqrt(1/(effectiveSigmaR*lumiRight*selectionRight) +lumiError*lumiError) << endl;
	cout << "F1VL: " << - F1Vg + FeL * F1Vz * (RePropagator) << " F1AL: " << FeL * F1Az * (RePropagator) << endl;
	cout << "F1VR: " << - F1Vg + FeR * F1Vz * RePropagator << " F1AR: " << FeR * F1Az * RePropagator << endl;
	cout <<"----------------ISR---------------\n";
	TCanvas * c2 = new TCanvas("c2", "Data-MC",0,500,1200,400);
	c2->Divide(3,1);
	c2->cd(1);
	float bisr = 2*alpha/pi*(log(s/0.5e-3/0.5e-3)-1);
	cout <<"Beta: " << N << endl;
	//TF1 * isr = new TF1("isr","[0]*pow((250-x)/125,[0]-1)",0,sqrts);
	//isr->SetParameter(0,bisr);
	//isr->Draw();
	TF1 * isr2 = new TF1("isr2","exp([0]/pow((250-x),[0]))",0,sqrts);
	isr2->SetParameters(2.e-01, 3.);
	TFile * file = TFile::Open("/exp/flc/bilokin/Training/bbbar-250GeV/recoverytest/TTBarProcessorRightSignalAll.root");
	TTree * GenTree =(TTree *) file->Get("GenTree");
	float nbins = 10;
	TH1F *  h1 = new TH1F("h1", "histo from a gaussian", nbins, 180, sqrts);
	int ntotal = GenTree->Draw("MCMass >> h1","MCMass > 180 && MCMass < 250");
	h1->Draw();
	TH1F * total_isr = new TH1F("total","",100,-1,1);
	makePretty(total_isr,kRed);
	TH1F * total = new TH1F("total","",50,-1,1);
	makePretty(total,kGreen);
	bool left = 0;
	//float Nevents = 194747*2;//40487;//125 * dsigmaR->Integral(-1,1);
	//float Nevents = 203248;//40487;//125 * dsigmaR->Integral(-1,1);
	float Nevents = 208825;//40487;//125 * dsigmaR->Integral(-1,1);
	//float Nevents = 40487;//125 * dsigmaR->Integral(-1,1);
	float Nevents = 43306;//125 * dsigmaR->Integral(-1,1);
	//float Nevents = 5670.22 * 100 ;//125 * dsigmaR->Integral(-1,1);
	TF1 * dsigmaR2 = new TF1("dsigmaR2","mnm(x,[0],[1], 0,0)",-1,1);
	dsigmaR2->SetParameter(0,sqrts);
	dsigmaR2->SetParameter(1,left);
	//dsigmaR2->Draw();
	//dsigmaL2->SetParameter(0,sqrts);
	total->FillRandom("dsigmaR2",Nevents);
	float oldA = dsigmaR2->Integral(-1,1) ;
	cout << "Nominal sigmaT: " << oldA << endl;
	//TF1 * dsigmaR = new TF1("dsigmaR","mnm(x,250,1,[0])",-1,1);
	TF1 * fit = new TF1("fit","(1+x*x)*[0]+[1]*x+(1-x*x)*[2]",-1,1);
	//TF1 * fit = new TF1("fit","(1+x*x)*[0]+[1]*x",-1,1);
	//dsigmaR->SetParameter(0,sqrts);
	total->Fit("fit");
	float newA = 0;
	cout << "P: " << mnm(0.5, 250, 1) << endl;
	float error = 0;
	getAfb(fit, 1, error, left);
	total->Draw();//*/
	float Al = 0;
	float Bl = 0;
	float Cl = 0;
	float dx = dsigma(0,250,1,Al,Bl,Cl);
	cout << "A: " << Al <<" B: " << Bl << " C: " << Cl << endl;
	float Ar = 0;
	float Br = 0;
	float Cr = 0;
	float dx2 = dsigma(0,250,0,Ar,Br,Cr);
	cout << "A: " << Ar <<" B: " << Br << " C: " << Cr << endl;
	//ILC();
	//cout << "Nominal Afb left: " << AFB(sqrts, 1) << " right: " << AFB(sqrts, 0) << endl;
	//TF1 * afbrf = new TF1("Afbrf","AFB(x,0)",10,260);
	//TF1 * afblf = new TF1("Afblf","AFB(x,1)",10,260);
	//afbrf->Draw();
	//afblf->Draw("same");
	/*for (unsigned int i = 0; i < nbins; i++) 
	{
		int nevents = h1->GetBinContent(i+1);
		int sqrtsprime = h1->GetBinCenter(i+1);
		dsigmaL2->SetParameter(0,sqrtsprime);
		total_isr->FillRandom("dsigmaL2",nevents);
		float Aprime = dsigmaL2->Integral(-1,1); 
		cout << "E: " << sqrtsprime << " N: " << nevents  << " A: " << Aprime << endl;
		newA += Aprime * nevents / ntotal;
	}
	cout << "ISR A: " << newA << "(" << newA/oldA*100 << "%)" << endl;
	total->Draw();
	total_isr->Draw("same");
	cout << "No ISR events: " << total->GetEntries() << " ISR events: " <<  total_isr->GetEntries()  << endl;
	TF1 * fit = new TF1("fit","(1+x*x)*[0]+[1]*x+(1-x*x)*[2]",-1,1);
	total->Fit("fit");
	total_isr->Fit("fit");//*/

}
void ILC()
{
	float Al = 0;
	float Bl = 0;
	float Cl = 0;
	float dx = dsigma(0,250,1,Al,Bl,Cl);
	cout << "A: " << Al <<" B: " << Bl << " C: " << Cl << endl;
	float Ar = 0;
	float Br = 0;
	float Cr = 0;
	float dx2 = dsigma(0,250,0,Ar,Br,Cr);
	cout << "A: " << Ar <<" B: " << Br << " C: " << Cr << endl;
	float factorl = 1.12;
	float dAl = 1.71e-03*factorl;
	float dBl = 3.72e-03*factorl;
	float dCl = 1.e-03*factorl;
	float factorr = 1.75;
	float dAr = 1.01e-03*factorr;
	float dBr = 2.44e-03*factorr;
	float dCr = 0.9e-03*factorr;
	int nsteps = 10;

	float precision = 0.005;
	
	float p1Vz = 0.02;
	float p2Vz = 0.003;
	
	float p1Vg = 0.016;
	float p2Vg = 0.002;
	
	float p1Az = 0.01;
	float p1Ag = 0.004;
	
	float highy = precision;
	float highx = precision*2.5;
	TH2F * F1VgF2Vg = new TH2F("F1VgF2Vg", "; F_{1V}^{#gamma};F_{2V}^{#gamma}", 1000, -p1Vg, p1Vg, 1000, -p2Vg, p2Vg);
	TH2F * F1VzF2Vz = new TH2F("F1VzF2Vz", "; F_{1V}^{z};F_{2V)^{z}", 1000, -p1Vz, p1Vz, 1000, -p2Vz, p2Vz);
	TH2F * F1AgF1Az = new TH2F("F1AgF1Az", "; F_{1A}^{#gamma};F_{1A}^{z}", 1000, -p1Ag,p1Ag, 1000, -p1Az, p1Az);
	TH1F * check = new TH1F("check", "; factor; factor", 100, -highx, highx);
	TH2F * check2 = new TH2F("check2", ";#delta F_{1V}^{#gamma}/F_{1V}^{#gamma};#delta F_{2V}^{#gamma}/F_{2V}^{#gamma}", 100, -p1Vg, p1Vg, 100, -p2Vg, p2Vg);
	TH2F * check3 = new TH2F("check3", ";#delta F_{1V}^{z}/F_{1V}^{z};#delta F_{2V}^{z}/F_{2V}^{z}", 100, -p1Vz, p1Vz, 100, -p2Vz,p2Vz);
	TH2F * check4 = new TH2F("check4", ";#delta F_{1A}^{#gamma}/F_{1A}^{#gamma} ;#delta F_{1A}^{z}/F_{1A}^{z}", 100, -p1Ag, p1Ag, 100, -p1Az, p1Az);
	//vector<float*> * allowed = new vector<float*>();
	float index = 0;
	//while (index < nsteps) 
	for (unsigned int i = 0; i < nsteps; i++) 
	{
		int kick = i%10000;
		float al = 0; float bl = 0; float cl = 0;
		float ar = 0; float br = 0; float cr = 0;
		TRandom3 *eG1 = new TRandom3();
		float toStart = i;///eG1->Uniform(1000);
		eG1->SetSeed(toStart);
		float dF1Vg = getUni(eG1, -0.013, 0.015);//eG1->Uniform(p1Vg*2)-p1Vg;
		eG1->SetSeed(abs(dF1Vg)/dF1Vg+kick);
		float dF1Ag = getUni(eG1, -0.011, 0.01);//eG1->Uniform(p1Ag*2)-p1Ag;
		eG1->SetSeed(abs(dF1Ag)/dF1Ag+i+kick);
		float dF1Vz = getUni(eG1, -0.015, 0.02);//eG1->Uniform(p1Vz*2)-p1Vz;
		eG1->SetSeed(abs(dF1Vz)/dF1Vz+i*2+kick);
		float dF1Az = getUni(eG1, -0.011, 0.01);//eG1->Uniform(p1Az*2)-p1Az;
		eG1->SetSeed(abs(dF1Az)/dF1Az+2*kick);
		float dF2Vg =  getUni(eG1, -0.0012, 0.002);//G1->Uniform(p2Vg*2)-p2Vg;
		eG1->SetSeed(abs(dF2Vg)/dF2Vg);
		float dF2Vz = getUni(eG1, -0.002, 0.003);//eG1->Uniform(p2Vz*2)-p2Vz;
		//cout << "\tdF1Vg: " << dF1Vg << " dF1Ag: " << dF1Ag << " dF1Vz: " << dF1Vz << " dF1Az: " << dF1Az << " dF2Vg: " << dF2Vg << " dF2Vz: " << dF2Vz << endl;
		dsigma(0.5,250,1, al, bl, cl, dF1Vg, dF1Ag, dF1Vz, dF1Az, dF2Vg, dF2Vz);
		dsigma(0.5,250,0, ar, br, cr, dF1Vg, dF1Ag, dF1Vz, dF1Az, dF2Vg, dF2Vz);
		//check->Fill(dF1Vz);
		//check2->Fill(dF1Vg, dF2Vg);
		//check3->Fill(dF1Vz, dF2Vz);
		//check4->Fill(dF1Ag, dF1Az);
		if (al > Al - dAl && al < Al + dAl && 
		    bl > Bl - dBl && bl < Bl + dBl && 
		    cl > Cl - dCl && cl < Cl + dCl && 
		    ar > Ar - dAr && ar < Ar + dAr && 
		    br > Br - dBr && br < Br + dBr && 
		    cr > Cr - dCr && cr < Cr + dCr) 
		{
			//cout << "\t\tOk!\n";
			//float ok[6] = {dF1Vg, dF1Ag, dF1Vz, dF1Az, dF2Vg, dF2Vz};
			F1VgF2Vg->Fill(dF1Vg,dF2Vg);
			F1VzF2Vz->Fill(dF1Vz,dF2Vz);
			F1AgF1Az->Fill(dF1Ag, dF1Az);
			index++;
			//allowed->push_back(ok);
		}
		//cout << "\tA: " << al <<" B: " << bl << " C: " << cl << endl;
		//cout << "\tA: " << ar <<" B: " << br << " C: " << cr << endl;
		delete eG1;
	}
	cout << "Nevents: " << index << endl;
	makePretty(F1VgF2Vg,kBlue);
	makePretty(F1VzF2Vz,kBlue);
	makePretty(F1AgF1Az,kBlue);
	makePretty(check2,kGray);
	makePretty(check3,kGray);
	makePretty(check4,kGray);
	//check->Draw("e");
	check2->Draw();
	F1VgF2Vg->Draw("same");
	gPad->SetLeftMargin(0.17);
	gPad->SetBottomMargin(0.14);
	gPad->SetRightMargin(0.05);
	check2->GetYaxis()->SetTitleSize(.048);
	check2->GetXaxis()->SetTitleSize(.048);
	check2->GetYaxis()->SetTitleOffset(1.8);
	c2->cd(2);
	check3->Draw();
	F1VzF2Vz->Draw("same");
	gPad->SetLeftMargin(0.17);
	gPad->SetBottomMargin(0.14);
	gPad->SetRightMargin(0.05);
	check3->GetYaxis()->SetTitleSize(.048);
	check3->GetXaxis()->SetTitleSize(.048);
	check3->GetYaxis()->SetTitleOffset(1.8);
	c2->cd(3);
	check4->Draw();
	F1AgF1Az->Draw("same");
	gPad->SetLeftMargin(0.17);
	gPad->SetBottomMargin(0.14);
	gPad->SetRightMargin(0.05);
	check4->GetYaxis()->SetTitleSize(.048);
	check4->GetXaxis()->SetTitleSize(.048);
	check4->GetYaxis()->SetTitleOffset(1.8);
}
float mnm(float x, float sqrts = 250., bool left = 1, float dF1Vg = 0., float dF1Ag = 0., float dF1Vz = 0., float dF1Az = 0., float dF2Vg = 0., float dF2Vz = 0.)
{
	float s = sqrts * sqrts;
	float pi = TMath::Pi();
	float Mz = 91.1876;
	float Gz = 2.5;
	float alpha = 1./128;
	float sw = 0.23116;
	float cw = 1 - sw;
	float As = 4*pi*alpha*alpha/3/s;

	float Qf = -1./3; float If = -0.5; float mf = 4.18;
	//float Qf = 2./3; float If = 0.5; float mf = 175;
	float gamma = sqrts / 2 / mf;
	float beta =  2 * sqrt( s/4 - mf*mf) / sqrts;
	float F1Vg = Qf*(1 + dF1Vg);
	float F1Ag = 0 + dF1Ag;
	float F2Vg = 0 + dF2Vg;

	float F2Vz = 0 + dF2Vz;
	float F1Vz = (If - 2*Qf * sw)/sqrt(4*cw*sw)*(1  + dF1Vz);
	float F1Az = - If / sqrt(4*cw*sw)*(1 + dF1Az);

	float FeL = (-0.5 + sw)/sqrt(cw*sw);
	float FeR =  sw/sqrt(cw*sw);
	float RePropagator = s*(s-Mz*Mz) / ((s-Mz*Mz)*(s-Mz*Mz) + Gz*Mz*Gz*Mz);
	float ImPropagator = - s*Gz*Mz / ((s-Mz*Mz)*(s-Mz*Mz) + Gz*Mz*Gz*Mz);

	float N = As*3*3/4*beta*  1e15 * 0.3894e-3 ;
	if (left) 
	{
		float mF1V2L = pow(- (F1Vg + F2Vg) + FeL * (F1Vz + F2Vz) * RePropagator,2)  + pow(FeL * (F1Vz + F2Vz) * ImPropagator,2);
		float mF1A2L = pow(- F1Ag + FeL * F1Az * RePropagator,2)  + pow(FeL * F1Az * ImPropagator,2);

		float mF1VmF1AL = F1Ag * F1Vg - F1Ag * FeL * F1Vz*RePropagator - F1Vg * FeL * F1Az * RePropagator +  FeL * F1Az  *  FeL * F1Vz *(RePropagator * RePropagator + ImPropagator * ImPropagator);
		float mF1V2Lg = pow(- (F1Vg/gamma + F2Vg*gamma) + FeL * (F1Vz/gamma + F2Vz*gamma) * RePropagator,2)  + pow(FeL * (F1Vz/gamma + F2Vz*gamma) * ImPropagator,2);
		
		float A = mF1V2L + beta*beta*mF1A2L;
		float B = -4*beta*mF1VmF1AL;
		float C = mF1V2Lg;

		return N*((1+x*x) * A + B*x + C*(1-x*x));
	}
	float mF1V2R = pow(- (F1Vg + F2Vg) + FeR * (F1Vz + F2Vz) * RePropagator,2)  + pow(FeR * (F1Vz + F2Vz) * ImPropagator,2);
	float mF1A2R = pow(- F1Ag + FeR * F1Az * RePropagator,2)  + pow(FeR * F1Az * ImPropagator,2);

	float mF1VmF1AR =-F1Ag * F1Vg + F1Ag * FeR * F1Vz*RePropagator + F1Vg * FeR * F1Az * RePropagator -  FeR * F1Az  *  FeR * F1Vz *(RePropagator * RePropagator + ImPropagator * ImPropagator);
	float mF1V2Rg = pow(- (F1Vg/gamma + F2Vg*gamma) + FeR * (F1Vz/gamma + F2Vz*gamma) * RePropagator,2)  + pow(FeR * (F1Vz/gamma + F2Vz*gamma) * ImPropagator,2);

	float A = mF1V2R +beta*beta* mF1A2R;
	float B = -4*beta*mF1VmF1AR;
	float C = mF1V2Rg;
	
	return N*((1+x*x) * A + B*x + C*(1-x*x));///*/
}
float dsigma(float x, float sqrts, bool left, float &A, float &B, float &C,  float dF1Vg = 0, float dF1Ag = 0, float dF1Vz = 0, float dF1Az = 0, float dF2Vg = 0, float dF2Vz = 0)
{
	float s = sqrts * sqrts;
	float pi = TMath::Pi();
	float Mz = 91.1876;
	float Gz = 2.5;
	float alpha = 1./128;
	float sw = 0.23116;
	float cw = 1 - sw;
	float As = 4*pi*alpha*alpha/3/s;

	float Qf = -1./3; float If = -0.5; float mf = 4.18;
	//float Qf = 2./3; float If = 0.5; float mf = 175;
	float gamma = sqrts / 2 / mf;
	float beta =  2 * sqrt( s/4 - mf*mf) / sqrts;
	float F1Vg = Qf*(1 + dF1Vg);
	float F1Ag = 0 + dF1Ag;
	float F2Vg = 0 + dF2Vg;

	float F2Vz = 0 + dF2Vz;
	float F1Vz = (If - 2*Qf * sw)/sqrt(4*cw*sw)*(1  + dF1Vz);
	float F1Az = - If / sqrt(4*cw*sw)*(1 + dF1Az);

	float FeL = (-0.5 + sw)/sqrt(cw*sw);
	float FeR =  sw/sqrt(cw*sw);
	float RePropagator = s*(s-Mz*Mz) / ((s-Mz*Mz)*(s-Mz*Mz) + Gz*Mz*Gz*Mz);
	float ImPropagator = - s*Gz*Mz / ((s-Mz*Mz)*(s-Mz*Mz) + Gz*Mz*Gz*Mz);

	float N = As*3*3/4*beta*  1e15 * 0.3894e-3 ;
	if (left) 
	{
		float mF1V2L = pow(- (F1Vg + F2Vg) + FeL * (F1Vz + F2Vz) * RePropagator,2)  + pow(FeL * (F1Vz + F2Vz) * ImPropagator,2);
		float mF1A2L = pow(- F1Ag + FeL * F1Az * RePropagator,2)  + pow(FeL * F1Az * ImPropagator,2);

		float mF1VmF1AL = F1Ag * F1Vg - F1Ag * FeL * F1Vz*RePropagator - F1Vg * FeL * F1Az * RePropagator +  FeL * F1Az  *  FeL * F1Vz *(RePropagator * RePropagator + ImPropagator * ImPropagator);
		float mF1V2Lg = pow(- (F1Vg/gamma + F2Vg*gamma) + FeL * (F1Vz/gamma + F2Vz*gamma) * RePropagator,2)  + pow(FeL * (F1Vz/gamma + F2Vz*gamma) * ImPropagator,2);
		
		A = mF1V2L + beta*beta*mF1A2L;
		B = -4*beta*mF1VmF1AL;
		C = mF1V2Lg;

		return N*((1+x*x) * A + B*x + C*(1-x*x));
	}
	float mF1V2R = pow(- (F1Vg + F2Vg) + FeR * (F1Vz + F2Vz) * RePropagator,2)  + pow(FeR * (F1Vz + F2Vz) * ImPropagator,2);
	float mF1A2R = pow(- F1Ag + FeR * F1Az * RePropagator,2)  + pow(FeR * F1Az * ImPropagator,2);

	float mF1VmF1AR =-F1Ag * F1Vg + F1Ag * FeR * F1Vz*RePropagator + F1Vg * FeR * F1Az * RePropagator -  FeR * F1Az  *  FeR * F1Vz *(RePropagator * RePropagator + ImPropagator * ImPropagator);
	float mF1V2Rg = pow(- (F1Vg/gamma + F2Vg*gamma) + FeR * (F1Vz/gamma + F2Vz*gamma) * RePropagator,2)  + pow(FeR * (F1Vz/gamma + F2Vz*gamma) * ImPropagator,2);

	A = mF1V2R +beta*beta* mF1A2R;
	B = -4*beta*mF1VmF1AR;
	C = mF1V2Rg;
	
	return N*((1+x*x) * A + B*x + C*(1-x*x));
}
void makePretty2(TH1 * vtxTotal, int color, int line = 0)
{
	vtxTotal->SetLineWidth(3);
	vtxTotal->SetLineStyle(line);
	vtxTotal->SetLineColor(color);
	vtxTotal->SetMinimum(0);
	vtxTotal->SetStats(0);
}
void makePretty2(TH2F * good, int color)
{
	good->SetStats(0);
	good->SetMarkerColor(color);
	good->SetFillColor(color);
	good->SetMarkerSize(0.5);
	good->SetMarkerStyle(20);
}
float getUni(TRandom * gen, float x1, float x2)
{
	float random = gen->Uniform(fabs(x2-x1));
	//cout << random + x1 << endl;
	return random + x1;
}


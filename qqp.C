#include <sstream>
void qqp(string filename = "TTBarProcessorLeft.root")
{
	TCanvas * c1 = new TCanvas("c1","",1000,500);
	c1->Divide(2,1);
	c1->cd(1);
	TFile * file = TFile::Open(filename.c_str());
	int bin_e = 20;
	int max_e = 1;
	TH1F * right = new TH1F("right",";cos",bin_e,-max_e, max_e); 
	TH1F * wrong = new TH1F("wrong",";cos",bin_e,-max_e, max_e);
	TH1F * gen = new TH1F("gen",";cos",bin_e,-max_e, max_e);
	TH1F * corrected = new TH1F("corrected",";cos",bin_e,-max_e, max_e);
	float Nrefused = Stats->Draw("methodSameCharge","methodSameCharge > 0 && MCMass > 200");
	float Naccepted = Stats->Draw("methodUsed","methodUsed > 0 && MCMass > 200");
	std::cout << "Naccepted: " << Naccepted <<" Nrefused: " << Nrefused << endl;
	float N = Nrefused + Naccepted;
	float p = 0.5 + sqrt(1-2*Nrefused/N)/2;
	float q = 1.-p;
	Stats->Draw("qCostheta>>right","methodUsed > 0 && MCMass > 200");
	GenTree->Draw("qMCcostheta>>gen","MCMass > 200");
	Stats->Draw("B1costheta >> wrong","methodSameCharge > 0 && MCMass > 200");
	Stats->Draw("B2costheta >> +wrong","methodSameCharge > 0 && MCMass > 200");
	std::cout << "Global: p = "<< p << "; q = " << q << endl;
	makePretty(right,kBlack);
	makePretty(wrong,kRed);
	makePretty(gen,kGreen);
	std::cout << "-------------------------\n";
	corrected = computeCorrection(right,wrong,bin_e);
	makePretty(corrected,kBlue);
	gen->Scale(1./gen->GetEntries() * right->GetEntries());
	corrected->Scale(1./corrected->Integral() * right->GetEntries());
	gen->Draw("");
	right->Draw("samee");
	corrected->Draw("samee");
	wrong->Draw("samee");
	c1->cd(2);
	std::cout << "-------------------------\n";
	TH1F * correctedg = new TH1F("correctedg",";cos",bin_e,-max_e, max_e);
	makePretty(correctedg,kBlue);
	string cutssame = "MCMass > 200 && methodSameCharge == ";
	string cutsused = "MCMass > 200 && methodUsed == ";
	for (unsigned int i = 3; i < 5; i++) 
	{
		string n = str(i);

		float Nrefused_i = Stats->Draw("methodSameCharge",string(cutssame+n).c_str());
		float Naccepted_i = Stats->Draw("methodUsed",string(cutsused+n).c_str());
		float p_i = getP(Naccepted_i, Nrefused_i);//0.5 + sqrt(1-2*Nrefused_i/(Nrefused_i + Naccepted_i))/2;
		float q_i = 1. - p_i;
		float qp_i = Nrefused_i / (2*p_i * (Nrefused_i + Naccepted_i)) + p_i - 1;  
		std::cout << "Naccepted_"<<n<<": " << Naccepted_i <<" Nrefused_"<<n<<": " << Nrefused_i << endl;
		std::cout << "Local: p_"<<n<<" = "<< p_i << "; q_"<<n<<" = " << q_i << "; q' = " << qp_i << endl << endl;
		TH1F * rightn = new TH1F(string("right"+n).c_str(),";cos",bin_e,-max_e, max_e); 
		TH1F * wrongn = new TH1F(string("wrong"+n).c_str(),";cos",bin_e,-max_e, max_e);
		makePretty(rightn,kBlack);
		makePretty(wrongn,kRed);
		Stats->Draw(string("qCostheta>>right"+n).c_str(),string(cutsused+n).c_str());
		Stats->Draw(string("B1costheta >> wrong"+n).c_str(),string(cutssame+n).c_str(),"same");
		Stats->Draw(string("B2costheta >> +wrong"+n).c_str(),string(cutssame+n).c_str(),"same");
		TH1F * correctedl = computeCorrection(rightn,wrongn,bin_e,n);
		std::cout << "Adding " << correctedl->GetEntries() << " events\n";
		correctedg->Add(correctedl);

	}
	std::cout << "-------------------------\n";
	c1->cd(2);
	gen->Draw("");
	right->Draw("samee");
	correctedg->Scale(1./correctedg->Integral() * right->GetEntries());
	correctedg->Draw("samee");//*/
	//std::cout << "CORRECTED AFB: " << getAfb(correctedg, bin_e) << "\n";
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
float getP(float Na, float Nr)
{
	float N = Na+Nr;
	return sqrt(1./2/N)*sqrt(Na + sqrt(Na*Na-Nr*Nr));
}
float Nplus(float p, float q, float Naplus, float Naminus)
{
	return (Naplus * p*p - Naminus * q*q)/(p*p*p*p - q*q*q*q);
}
float Nminus(float p, float q, float Naplus, float Naminus)
{
	return (Naplus * q*q - Naminus * p*p)/(q*q*q*q - p*p*p*p );
}
void makePretty(TH1 * vtxTotal, int color)
{
	vtxTotal->SetLineWidth(3);
	vtxTotal->SetLineColor(color);
	vtxTotal->SetMinimum(0);
	vtxTotal->SetStats(0);
}
TH1F * computeCorrection(TH1F * right, TH1F * wrong, int bin_e = 20, string n = "")
{
	if (right->GetEntries() <  wrong->GetEntries()) 
	{
		std::cout << "ERROR " << right->GetEntries() << " " << wrong->GetEntries() << "!\n";
		return right;
	}
	TH1F * corrected = new TH1F(string("corrected"+n).c_str(),";cos",bin_e,-1, 1);
	float ntminus = 0.0;
	float ntplus = 0.0;
	float naminus = 0.0;
	float naplus = 0.0;
	float nplus = 0.0;
	float nminus = 0.0;
	int niterations = (float) bin_e/2;
	for (unsigned int i = 1; i < niterations+1; i++) 
	{
		int j = bin_e + 1 - i;
		float Naccepted_i = right->GetBinContent(i) + right->GetBinContent(j);
		float Nrefused_i = wrong->GetBinContent(i) + wrong->GetBinContent(j);
		float p_i = getP(Naccepted_i, Nrefused_i);//0.5 + sqrt(1-2*Nrefused_i/(Nrefused_i + Naccepted_i))/2;
		float q_i = 1. - p_i;
		float qp_i = Nrefused_i / (2*p_i * (Nrefused_i + Naccepted_i)) + p_i - 1;  
		std::cout << "Local: p_i = "<< p_i << "; q_i = " << q_i << "; q' = " << qp_i << endl;
		float ntminus_i = gen->GetBinContent(i);
		float ntplus_i = gen->GetBinContent(j);
		float naminus_i = right->GetBinContent(i);
		float naplus_i = right->GetBinContent(j);
		float nplus_i = Nplus(p_i,q_i, naplus_i, naminus_i);
		float nminus_i = abs(Nminus(p_i,q_i, naplus_i, naminus_i));
		float error_minus = sigma(nminus_i, nplus_i, p_i,q_i, Nrefused_i, Naccepted_i);
		std::cout << "Measured: Na+ = " << naplus_i << "; Na- = " << naminus_i <<endl;
		std::cout << "Corrected: N+ = " << nplus_i << "; N- = " << nminus_i << "; sN- = " << error_minus <<endl;
		corrected->SetBinContent(i, nminus_i);
		corrected->SetBinError(i,error_minus);

		corrected->SetBinContent(j, nplus_i);
		corrected->SetBinError(j,sigma(nplus_i, nminus_i, p_i,q_i, Nrefused_i, Naccepted_i));
		naminus += naminus_i;
		naplus += naplus_i;
		nplus += nplus_i;
		nminus += nminus_i;
		ntminus += ntminus_i;
		ntplus += ntplus_i;
	}
	std::cout << "-------------------------\n";
	float afbgen = (ntplus-ntminus)/(ntplus+ntminus);
	float afbmeasured =  (naplus-naminus)/(naplus+naminus);
	float afbcorrected =  (nplus-nminus)/(nplus+nminus);
	std::cout << "Gen Afb = " << afbgen << endl;
	std::cout << "Measured Afb = " << afbmeasured << " (" << afbmeasured/afbgen*100 <<")%" << "; Afb* = " << afbcorrected << " (" << afbcorrected/afbgen*100 <<")%" << endl;
	return corrected;
}
float sigma(float naplus, float naminus, float p, float q, float nr, float n)
{
	float alpha = abs(1/(p*p*p*p-q*q*q*q));
	//float n = naplus + naminus;
	float dpdNa = nr / (2* n*n *sqrt(1-2*nr/n));
	//std::cout << "N: " << n << " Nr: " << nr << " dp: " << dpdNa <<  endl;
	return sqrt(alpha * alpha* ( naplus * p*p*p*p + naminus * q*q*q*q ) + naplus);
	//return alpha * sqrt( naplus * p*p*p*p + naminus * q*q*q*q );
	//return alpha*sqrt( naplus*pow(p*p+ 2*p*naplus*dpdNa - 2*q* naminus*dpdNa,2) + naminus * pow(2*p*naplus*dpdNa -q*q - 2*q* naminus*dpdNa,2) );
	//return alpha*sqrt( naplus*pow(p*p+ 2*p*naplus*dpdNa,2) + naminus * pow( -q*q - 2*q* naminus*dpdNa,2) );
}

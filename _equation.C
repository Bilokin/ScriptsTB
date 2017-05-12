float _equation()
{
	float Bpm = 0.45;
	float B0 = 0.55;

	//float N0 = 0.42;//0.44;
	//float Na = 0.485;//0.458;
	float N0 = 0.438;//0.44;
	float Na = 0.466;//0.458;
	float Nr = 1- Na - N0;

	float p_v = (N0 - 0.5*Bpm)/(B0 - 0.5*Bpm);
	float q_v = 1- p_v;
	cout << "p_v = " << p_v << endl;
	cout << "p'_v = " << p_v+0.5*q_v << endl;

	float p_kpm = (Na-0.5*q_v*B0) / ((p_v+0.5*q_v)*Bpm);

	cout << "p_kpm = " << p_kpm << endl;

	float checkNr = 0.5*q_v*B0 + p_v*(1-p_kpm)*Bpm + 0.5 * q_v*(1-p_kpm)*Bpm;
	cout << "Check = " << checkNr << endl;

	float Nr1 = 2 * Bpm*B0 * (0.5 * q_v * p_v + 0.25 * q_v*q_v) + B0*B0 * 0.5*q_v*q_v;
	float Na1 = Bpm*Bpm * (p_v*p_v + q_v*p_v + 0.25*q_v*q_v) + 2 * Bpm*B0 * (0.5*q_v*p_v + 0.25*q_v*q_v) + B0*B0*0.5*q_v*q_v;
	float N01 = Bpm*Bpm * (q_v*p_v + 0.75 * q_v*q_v) +  2 * Bpm*B0 * (p_v*p_v +q_v*p_v + 0.5*q_v*q_v) + B0*B0*(p_v*p_v + 2*p_v*q_v);
	cout << "Check1: Na = " << Na1 << "; Nr = " << Nr1 << "; N0 = " << N01 << endl;
	TF1 * fNr = new TF1("fNr;p_v","2*(1-[0])*[0] * (0.5* (1-x)*x + 0.25 * (1-x)*(1-x)) +0.5*[0]*[0]*(1-x)*(1-x)",0,1);
	TF1 * fNa = new TF1("fNa;p_v","(1-[0])*(1-[0])* (x*x+(1-x)*x + 0.25*(1-x)*(1-x)) + 2*(1-[0])*[0] * (0.5* (1-x)*x + 0.25*(1-x)*(1-x)) + [0]*[0]*0.5*(1-x)*(1-x)",0,1);
	TF1 * fN0 = new TF1("fN0;p_v","(1-[0])*(1-[0])* ((1-x)*x + 0.75*(1-x)*(1-x)) + 2*(1-[0])*[0] *(x*x+(1-x)*x + 0.5*(1-x)*(1-x))+ [0]*[0]*(x*x+2*x*(1-x))",0,1);
	//TF1 * fNa = new TF1("fNa","[1]*x*(1-[0])+0.5*(1-x)*[0]+0.5*(1-x)*[1]*(1-[0])",0,1);
	//TF1 * fN0 = new TF1("fN0","0.5*(1-x)*(1-[0])+x*[0]",0,1);
	//TF1 * fNr = new TF1("fNr","0.5*(1-x)*[0]+x*(1-[1])*(1-[0]) + 0.5*(1-[1])*(1-x)*(1-[0])",0,1);
	TF1 * fN0val = new TF1("fN0val","0.693",0,1);
	TF1 * fNaval = new TF1("fNaval","0.229",0,1);
	TF1 * fNrval = new TF1("fNrval","0.078",0,1);
	fN0val->SetLineStyle(2);
	fNaval->SetLineStyle(2);
	fNrval->SetLineStyle(2);
	string f = fNa->GetExpFormula() + "+" + fNr->GetExpFormula() + "+" + fN0->GetExpFormula();
	TF1 * fNt = new TF1("fNt",f.c_str(),0,1);
	//TF1 * fNt = new TF1("fNr","[1]*x*(1-[0])+0.5*(1-x)*[0]+0.5*(1-x)*[1]*(1-[0])+0.5*(1-x)*(1-[0])+x*[0]+0.5*(1-x)*([0])+x*(1-[1])*(1-[0]) + 0.5*(1-[1])*(1-x)*(1-[0])",0,1);
	p_kpm = 0.95;
	B0 = 0.54;
	fNa->SetLineColor(kBlue);
	fNr->SetLineColor(kRed);
	fN0->SetLineColor(kGreen);
	fNt->SetLineColor(kBlack);
	fNa->SetParameter(0, B0);
	fNa->SetParameter(1, p_kpm);
	fN0->SetParameter(0, B0);
	fN0->SetParameter(1, p_kpm);
	fNr->SetParameter(0, B0);
	fNr->SetParameter(1, p_kpm);
	fNt->SetParameter(0, B0);
	fNt->SetParameter(1, p_kpm);
	fNa->SetMaximum( 1);
	fNa->SetMinimum( 0);
	fNa->GetXaxis()->SetTitle("p_{v}");
	fNa->Draw();
	fN0->Draw("same");
	fNr->Draw("same");
	fNt->Draw("same");
	fN0val->Draw("same");
	fNaval->Draw("same");
	fNrval->Draw("same");

}


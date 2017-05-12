void _cvar()
{
	TCanvas * c1 = new TCanvas("c1", "Data-MC",0,0,500,500);
	TF1 * f1 = new TF1("f1","(2525.17*(1+x*x) + 4763.14*x+13.0821*(1-x*x))/(2525.17*(1+x*x) + 4763.14*x+13.0821*(1-x*x))",-1,1);
	f1->SetLineStyle(1);
	f1->Draw();
	TF1 * f2 = new TF1("f2","(2525.17*(1+x*x) + 4763.14*x+26*(1-x*x))/(2525.17*(1+x*x) + 4763.14*x+13.0821*(1-x*x))",-1,1);
	f2->SetLineStyle(2);
	f2->Draw("same");
	TF1 * f2 = new TF1("f2","(2525.17*(1+x*x) + 4763.14*x+39*(1-x*x))/(2525.17*(1+x*x) + 4763.14*x+13.0821*(1-x*x))",-1,1);
	f2->SetLineStyle(2);
	f2->SetLineColor(kGreen);
	f2->Draw("same");
	TF1 * f2 = new TF1("f2","(2515.17*(1+x*x) + 4763.14*x+39*(1-x*x))/(2525.17*(1+x*x) + 4763.14*x+13.0821*(1-x*x))",-1,1);
	f2->SetLineStyle(3);
	f2->SetLineColor(kGreen+1);
	f2->Draw("same");
	TF1 * f2 = new TF1("f2","(2525.17*(1+x*x) + 4763.14*x-13*(1-x*x))/(2525.17*(1+x*x) + 4763.14*x+13.0821*(1-x*x))",-1,1);
	f2->SetLineStyle(2);
	f2->SetLineColor(kGreen+2);
	f2->Draw("same");
	TF1 * f2 = new TF1("f2","(2515.17*(1+x*x) + 4753.14*x+36*(1-x*x))/(2525.17*(1+x*x) + 4763.14*x+13.0821*(1-x*x))",-1,1);
	f2->SetLineStyle(2);
	f2->SetLineColor(kGreen-1);
	f2->Draw("same");
	TF1 * f2 = new TF1("f2","(2515.17*(1+x*x) + 4753.14*x+2*(1-x*x))/(2525.17*(1+x*x) + 4763.14*x+13.0821*(1-x*x))",-1,1);
	f2->SetLineStyle(2);
	f2->SetLineColor(kGreen-2);
	f2->Draw("same");
}

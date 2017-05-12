void LRAsymmetry()
{
	
	TCanvas * c1 = new TCanvas("c1", "Data-MC",0,0,1000,500);
	c1->Divide(2,1);
	c1->cd(1);
	gROOT->ProcessLine(".x /users/flc/bilokin/Processors/Macros/top/asymmetry.C(\"TTBarProcessorLeft.root\",c1)");
	c1->cd(2);
	gROOT->ProcessLine(".x /users/flc/bilokin/Processors/Macros/top/asymmetry.C(\"TTBarProcessorRight.root\",c1)");
}

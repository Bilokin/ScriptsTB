#include "/users/flc/bilokin/Processors/Macros/top/basymmetry.C"
void LRBasymmetry()
{
	
	TCanvas * c1 = new TCanvas("c1", "Data-MC",0,0,1000,500);
	c1->Divide(2,1);
	c1->cd(1);
	basymmetry("TTBarProcessorLeft.root",c1);
//	gROOT->ProcessLine(".x /users/flc/bilokin/Processors/Macros/top/basymmetry.C(\"TTBarProcessorLeft.root\",c1)");
	c1->cd(2);
	basymmetry("TTBarProcessorRight.root",c1, false);
	//gROOT->ProcessLine(".x /users/flc/bilokin/Processors/Macros/top/basymmetry.C(\"TTBarProcessorRight.root\",c1)");
}

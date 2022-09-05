#include <cmath> 

void yield_cs(int syst,TString varExp){

	
	double brafrac_ntKp = 0.00102;
	double brafrac_ntKp_err = 0.000019;

	double lumi = 302.3;
	
	TString* var;
	if(varExp == "By"){var = new TString("Y");}
	if(varExp == "nMult"){var = new TString("Mult");}

	TFile *diff_f = new TFile(Form("./results/ntKp_%s_ratio.root",varExp.Data()),"read");
	TFile *eff_f = new TFile(Form("./FinalFiles/BPPPCorrYield%sNoTnP.root",var->Data()),"read");

	TMultiGraph *TMG_diff = (TMultiGraph*) diff_f->Get("TG");
	TH1D *TH_eff = (TH1D*) eff_f->Get("hInvEff");

	auto tl_diff = TMG_diff->GetListOfGraphs();


	auto TG_diff = (TGraphAsymmErrors *) tl_diff->At(0);
	
	const int _nBins = (TG_diff->GetN());
	TCanvas* c=new TCanvas();
	TLegend* leg_diff=new TLegend(0.7,0.7,0.9,0.9);
	TMultiGraph* m_diff=new TMultiGraph();
	double x[_nBins];
	double y[_nBins];
	double ex[_nBins];
	double ey[_nBins];
	
	for (int i=0;i<_nBins;i++){
		y[i]=TG_diff->GetY()[i] / (TH_eff->GetBinContent(i+1)*brafrac_ntKp*lumi);
		x[i]=TG_diff->GetX()[i];
		ey[i]=abs(y[i]) * sqrt(pow(TG_diff->GetErrorY(i)/TG_diff->GetY()[i],2) + pow(TH_eff->GetBinError(i+1)/TH_eff->GetBinContent(i+1),2) + pow(brafrac_ntKp_err/brafrac_ntKp,2) + pow(0.019,2));
		ex[i]=TG_diff->GetErrorX(i);
	}
/*
	if (syst==1){
		auto TG1_syst = (TGraphAsymmErrors *) tl1->At(1);
  		auto TG2_syst = (TGraphAsymmErrors *) tl2->At(1);
  		double x_syst[_nBins];
		double y_syst[_nBins];
		double ex_syst[_nBins];
		double ey_syst[_nBins];
	
		for (int i=0;i<_nBins;i++){
			y_syst[i]=TG1_syst->GetY()[i]/TG2_syst->GetY()[i];
			x_syst[i]=(TG1_syst->GetX()[i]+TG2_syst->GetX()[i])/2;
			ey_syst[i]=abs(y_syst[i])*sqrt(pow(TG1_syst->GetErrorY(i)/TG1_syst->GetY()[i],2)+pow(TG2_syst->GetErrorY(i)/TG2_syst->GetY()[i],2));
			ex_syst[i]=sqrt(pow(TG1_syst->GetErrorX(i),2)+pow(TG2_syst->GetErrorX(i),2))/2;
		}
		
		TGraph *g_diff_syst= new TGraphErrors (_nBins,x_syst,y_syst,ex_syst,ey_syst);
		g_diff_syst->SetLineColor(2);
		m_diff->Add(g_diff_syst);
		leg_diff->AddEntry(g_diff_syst,"Systematic Uncertainty", "e");
	}
*/
	
	TGraph *g_diff= new TGraphErrors (_nBins,x,y,ex,ey);
	g_diff->SetMarkerColor(1);
	g_diff->SetMarkerStyle(21);
	
	m_diff->Add(g_diff);
	m_diff->GetYaxis()->SetTitle("Differential Cross Section");
	
	double yield_max=0;
	for(int i = 0; i < _nBins; i++){
		if(y[i] > yield_max){
		yield_max = y[i];
		}
	}
	m_diff->GetYaxis()->SetRangeUser(0, yield_max * 1.5);
	m_diff->GetYaxis()->SetTitle(Form("d#sigma/d%s", varExp.Data()));
	if(varExp == "By"){
		 m_diff->GetXaxis()->SetTitle("Rapidity (y)");
		 m_diff->GetYaxis()->SetTitle("d#sigma/dy (m^{2})");
		 m_diff->GetXaxis()->SetLimits(-2.4 ,2.4);
	 }
	 if(varExp == "Bpt"){
		 m_diff->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
		 //m_diff->GetYaxis()->SetTitle("dN_{S}/dp_{T}");
		//if (tree == "ntKp"){ m_diff->GetXaxis()->SetLimits(0 ,80); }
		//if (tree == "ntphi"){ m_diff->GetXaxis()->SetLimits(0 ,60); }
	 }
	 if(varExp == "nMult"){
		 m_diff->GetXaxis()->SetTitle("Multiplicity (Mult)");
		 m_diff->GetYaxis()->SetTitle("d#sigma/dMult (m^{2})");
		 m_diff->GetXaxis()->SetLimits(0, 110);
	 }

	m_diff->Draw("ap");
	leg_diff->AddEntry(g_diff, "Statistical Uncertainty", "e");
	leg_diff->SetBorderSize(0);
	leg_diff->SetFillStyle(0);
	leg_diff->SetTextSize(0);
	leg_diff->Draw();
	c->SaveAs(Form("./results/ntKp_%s_diffcrosssectionplot.png",varExp.Data())); 
	return;
}

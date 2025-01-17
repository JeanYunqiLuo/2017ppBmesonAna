#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom.h"
#include <iostream>
#include <iomanip>
#include <fstream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"

#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
#include "TLine.h"
//#include "tnp_weight_lowptPbPb.h"



//#include "his.h"
using namespace std;

using std::cout;
using std::endl;

void BPRAA(){



	
	gStyle->SetOptStat(0);

	TCanvas * c = new TCanvas("c","c",600,600);

	c->cd();

	c->SetLeftMargin(0.16);


	TString InfileBP = "../../BP/EffAna/FinalFiles/BPPPCorrYieldPT.root";

	TFile * FileBP = new TFile(InfileBP.Data());
	TH1D * BPCross = (TH1D *) FileBP->Get("hPtSigma");
	BPCross->SetMarkerStyle(20);
	BPCross->SetMarkerSize(1);
	BPCross->SetMarkerColor(1);
	BPCross->SetLineColor(1);



	//B+ PbPb//
	
	const int NBins = 7;

	float BPXsecPPY[NBins];
	float BPXsecPPX[NBins] = {6,8.5,12.5,17.5,25,40,55};

	float BPXSecPPYErrUp[NBins];
	float BPXSecPPYErrDown[NBins];


	float BPXSecPPYErrUpRatio[NBins];
	float BPXSecPPYErrDownRatio[NBins];


	for(int i = 0; i < NBins; i++){

		BPXsecPPY[i] = BPCross->GetBinContent(i+1);
		BPXSecPPYErrUp[i] = BPCross->GetBinError(i+1);
		BPXSecPPYErrDown[i] = BPCross->GetBinError(i+1);
		BPXSecPPYErrUpRatio[i] = BPXSecPPYErrUp[i]/BPXsecPPY[i];
		BPXSecPPYErrDownRatio[i] = BPXSecPPYErrDown[i]/BPXsecPPY[i];
		
	}


	float BPXSecPPXErrUp[NBins] = {1,1.5,2.5,2.5,5,10,5};
	float BPXSecPPXErrDown[NBins] = {1,1.5,2.5,2.5,5,10,5};
	


	float BPXsecPbPbY[NBins] = {4.82132e+06/11.1,311668,270167,64384.4,208537/11.1,28700.6/11.1,7000.73/11.1};
	float BPXsecPbPbX[NBins] = {6,8.73,12.4,17.2,25,40,55};


	float BPXSecPbPbXErrUp[NBins] = {1,1.27,2.6,2.8,5,10,5};
	float BPXSecPbPbXErrDown[NBins] = {1,1.23,2.4,2.2,5,10,5};

	float BPXSecPbPbYErrUpRatio[NBins] = {0.278198,0.159,0.041,0.0654,0.0690334,0.104543,0.24575};
	float BPXSecPbPbYErrDownRatio[NBins] = {0.278198,0.145,0.0795,0.065,0.0690334,0.104543,0.24575};


	float BPXSecPbPbYErrUp[NBins];
	float BPXSecPbPbYErrDown[NBins];

	for(int i = 0; i < NBins; i++){

		BPXSecPbPbYErrUp[i] = BPXSecPbPbYErrUpRatio[i] * BPXsecPbPbY[i];
		BPXSecPbPbYErrDown[i] = BPXSecPbPbYErrDownRatio[i] * BPXsecPbPbY[i];

	}



	//Syst//
  TString errorFile = "../../2DMapSyst/OutFiles/BPError2D.root";
  TFile fError(errorFile);

	TH1D * TnPSyst = (TH1D *) fError.Get("TnPSyst");
	TH1D * BptSyst = (TH1D *) fError.Get("BptSyst");
	TH1D * BDTSyst = (TH1D *) fError.Get("BDTSyst");

  TString pdfErrorFile = "../../bp_pdf.root";
  TFile fPdfError(pdfErrorFile);
  TGraph* pdfSyst = (TGraph *) fPdfError.Get("bp_error");

  TString trackSelErrorFile = "../../syst_track_sel.root";
  TFile fTrackSelError(trackSelErrorFile);
  TGraph* trackSelSyst = (TGraph *) fTrackSelError.Get("bp_track_sel_error");

	float BPXSecPPYSystUp[NBins];
	float BPXSecPPYSystDown[NBins];

  // percent error
	float BPTrackingSyst[NBins] = {[0 ... NBins - 1] = 5};
	float BPMCDataSyst[NBins];
	float BPPDFSyst[NBins];
	float BPTrackSelSyst[NBins];
	float BPPtShapeSyst[NBins];
	float BPTnPSystDown[NBins];
	float BPTnPSystUp[NBins];

  // Get systematics from input files
  for (auto ibin = 0; ibin < NBins; ++ibin) {
    BPMCDataSyst[ibin] = BDTSyst->GetBinContent(ibin + 1);
    BPPtShapeSyst[ibin] = BptSyst->GetBinContent(ibin + 1);
    BPTnPSystDown[ibin] = TnPSyst->GetBinContent(ibin + 1);
    // TnP systematics are symmetric in the binned pT case
    BPTnPSystUp[ibin] = BPTnPSystDown[ibin];
    BPPDFSyst[ibin] = pdfSyst->GetY()[ibin];
    BPTrackSelSyst[ibin] = trackSelSyst->GetY()[ibin];
  }

  // RMS of all the errors
	float BPTotalSystDownRatio[NBins];
	float BPTotalSystUpRatio[NBins];

	for(int i = 0; i < NBins; i++){
		BPTotalSystDownRatio[i] =
      TMath::Sqrt(TMath::Power(BPTrackingSyst[i], 2) + TMath::Power(BPMCDataSyst[i], 2) +
                  TMath::Power(BPPDFSyst[i], 2) + TMath::Power(BPTrackSelSyst[i], 2) +
                  TMath::Power(BPPtShapeSyst[i], 2) + TMath::Power(BPTnPSystDown[i], 2)) / 100;
		BPTotalSystUpRatio[i] = BPTotalSystDownRatio[i];
	}

	for(int i = 0; i < NBins; i++){

		BPXSecPPYSystUp[i] = BPXsecPPY[i] * ( BPTotalSystUpRatio[i]);
		BPXSecPPYSystDown[i] = BPXsecPPY[i] * (BPTotalSystDownRatio[i] );
		//cout << "i = " << i << "     BPXSecPPYSystDown[i] = " << BPXSecPPYSystDown[i] << "  BPXSecPPYErrDown[i] =  " << BPXSecPPYErrDown[i] << endl;
		
	}


	//PbPb

	float BPXSecPbPbYSystUpRatio[NBins] = {0.3577,0.1404,0.1714,0.0775,0.0858,0.0715,0.1253};
	float BPXSecPbPbYSystDownRatio[NBins] = {0.3210,0.1359,0.1705,0.0761,0.0843,0.0699,0.1220};


	float BPXSecPbPbYSystUp[NBins];
	float BPXSecPbPbYSystDown[NBins];


	for(int i = 0; i < NBins; i++){

		BPXSecPbPbYSystDown[i] = (BPXSecPbPbYSystDownRatio[i]) * BPXsecPbPbY[i];
		BPXSecPbPbYSystUp[i] = (BPXSecPbPbYSystUpRatio[i]) * BPXsecPbPbY[i];

	}



	








	TH2D * HisEmpty = new TH2D("HisEmpty","",100,5,60,100,100.0,2000000);
	HisEmpty->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	HisEmpty->GetYaxis()->SetTitle("Cross Section (Or Production Yield)");
	HisEmpty->GetXaxis()->CenterTitle();
	HisEmpty->GetYaxis()->CenterTitle();
	HisEmpty->GetYaxis()->SetTitleOffset(1.8);
	HisEmpty->Draw();

	TGraphAsymmErrors *BPPPCrossGraph = new TGraphAsymmErrors(NBins, BPXsecPPX, BPXsecPPY,BPXSecPPXErrDown, BPXSecPPXErrUp,BPXSecPPYErrDown,BPXSecPPYErrUp);
	
	TGraphAsymmErrors *BPPbPbCrossGraph = new TGraphAsymmErrors(NBins, BPXsecPbPbX, BPXsecPbPbY,BPXSecPbPbXErrDown, BPXSecPbPbXErrUp,BPXSecPbPbYErrDown,BPXSecPbPbYErrUp);

	BPPbPbCrossGraph->SetLineColor(kGreen+2);
//	BPPbPbCrossGraph->SetFillColorAlpha(kGreen-9,0.5);
	BPPbPbCrossGraph->SetMarkerStyle(20);
	BPPbPbCrossGraph->SetMarkerSize(1);
	BPPbPbCrossGraph->SetMarkerColor(kGreen+2);

	BPPPCrossGraph->SetLineColor(kBlue+2);
//	BPPPCrossGraph->SetFillColorAlpha(kBlue-9,0.5);
	BPPPCrossGraph->SetMarkerStyle(21);
	BPPPCrossGraph->SetMarkerSize(1);
	BPPPCrossGraph->SetMarkerColor(kBlue+2);


	TGraphAsymmErrors *BPPPCrossGraphSyst  = new TGraphAsymmErrors(NBins, BPXsecPPX, BPXsecPPY, BPXSecPPXErrDown, BPXSecPPXErrUp, BPXSecPPYSystDown,BPXSecPPYSystUp);
  	TGraphAsymmErrors *BPPbPbCrossGraphSyst    = new TGraphAsymmErrors(NBins, BPXsecPbPbX, BPXsecPbPbY, BPXSecPbPbXErrDown, BPXSecPbPbXErrUp, BPXSecPbPbYSystDown,BPXSecPbPbYSystUp);


	BPPPCrossGraphSyst->SetFillColorAlpha(kBlue-9,0.5);
	BPPPCrossGraphSyst->SetLineColor(kBlue-9);
	BPPbPbCrossGraphSyst->SetFillColorAlpha(kGreen-9,0.5);
	BPPbPbCrossGraphSyst->SetLineColor(kGreen-9);



	BPPPCrossGraph->Draw("epsame");	
	BPPbPbCrossGraph->Draw("epsame");	

	BPPPCrossGraphSyst->Draw("5same");	
	BPPbPbCrossGraphSyst->Draw("5same");	
	

	TLegend* leg = new TLegend(0.45,0.65,0.75,0.85,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.040);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);
	leg->AddEntry(BPPbPbCrossGraph,"2018 PbPb 5.02 TeV","PL");
	leg->AddEntry(BPPPCrossGraph,"2017 pp 5.02 TeV","PL");
	leg->Draw("same");


	c->SaveAs("RAAPlots/BP/BPPbPbPPCross.png");

	c->SetLogy();

	c->SaveAs("RAAPlots/BP/BPPbPbPPCrossLog.png");






	//Fid and Non-Fid

	
	//Fid Non Fid

	TH2D * HisEmpty6 = new TH2D("HisEmpty6","",100,5,60,100,100.0,8000000);
	HisEmpty6->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	HisEmpty6->GetYaxis()->SetTitle("Cross Section (Or Production Yield)");
	HisEmpty6->GetXaxis()->CenterTitle();
	HisEmpty6->GetYaxis()->CenterTitle();
	HisEmpty6->GetYaxis()->SetTitleOffset(1.8);
	HisEmpty6->Draw();

	

	// TString InfileBP2 = "../Comparisons/NoFiducial/FinalFiles/BPPPCorrYieldPT.root";
	TString InfileBP2 = "../../../braa_nofid/BP/EffAna/FinalFiles/BPPPCorrYieldPT.root";


	TFile * FileBP2 = new TFile(InfileBP2.Data());
	TH1D * BPCross2 = (TH1D *) FileBP2->Get("CorrDiffHisBin");
	BPCross2->SetMarkerStyle(20);
	BPCross2->SetMarkerSize(1);
	BPCross2->SetMarkerColor(1);
	BPCross2->SetLineColor(1);


	


	float BPXsecPPY2[NBins];


	float BPXSecPPYErrUpRatio2[NBins];
	float BPXSecPPYErrDownRatio2[NBins];
	float BPXSecPPYErrUp2[NBins];
	float BPXSecPPYErrDown2[NBins];

	for(int i = 0; i < NBins; i++){

		BPXsecPPY2[i] = BPCross2->GetBinContent(i+1);
		BPXSecPPYErrUp2[i] = BPCross2->GetBinError(i+1);
		BPXSecPPYErrDown2[i] = BPCross2->GetBinError(i+1);
		BPXSecPPYErrUpRatio2[i] = BPXSecPPYErrUp2[i]/BPXsecPPY2[i];
		BPXSecPPYErrDownRatio2[i] = BPXSecPPYErrDown2[i]/BPXsecPPY2[i];
		
	}

	TGraphAsymmErrors *BPPPCrossGraph2 = new TGraphAsymmErrors(NBins, BPXsecPPX, BPXsecPPY2,BPXSecPPXErrDown, BPXSecPPXErrUp,BPXSecPPYErrDown2,BPXSecPPYErrUp2);


	BPPPCrossGraph2->SetLineColor(kGreen+2);
//	BPPPCrossGraph->SetFillColorAlpha(kBlue-9,0.5);
	BPPPCrossGraph2->SetMarkerStyle(20);
	BPPPCrossGraph2->SetMarkerSize(1);
	BPPPCrossGraph2->SetMarkerColor(kGreen+2);

	BPPPCrossGraph2->Draw("ep");
	BPPPCrossGraph->Draw("epSAME");	

	TLegend* leg5 = new TLegend(0.45,0.55,0.75,0.80,NULL,"brNDC");
	leg5->SetBorderSize(0);
	leg5->SetTextSize(0.040);
	leg5->SetTextFont(42);
	leg5->SetFillStyle(0);
	leg5->SetLineWidth(3);
	leg5->AddEntry(BPPPCrossGraph,"Fiducial Region Used","PL");
	leg5->AddEntry(BPPPCrossGraph2,"No Fiducial Region","PL");
	leg5->Draw("same");





    c->SaveAs("RAAPlots/BP/BPFidOrNotComp.png");





	TCanvas * c2 = new TCanvas("c","c",600,600);

	c2->cd();

	c2->SetLeftMargin(0.16);



	TH2D * HisEmptyRAA = new TH2D("HisEmptyRAA","",100,5,60,100,0,1.3);
	HisEmptyRAA->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
	HisEmptyRAA->GetYaxis()->SetTitle("RAA = #frac{1}{TAA} #frac{dN_{PbPb}/dp_{T}}{d #sigma_{pp}/d p_{T}}");
	HisEmptyRAA->GetXaxis()->CenterTitle();
	HisEmptyRAA->GetYaxis()->CenterTitle();
	HisEmptyRAA->GetYaxis()->SetTitleOffset(1.8);
	
	HisEmptyRAA->Draw();

	


	float BPRAAY[NBins];
//	float BPRAAX[NBins] = {8.73,12.4,17.2,27.3};
	float BPRAAX[NBins] = {6,8.73,12.4,17.2,25,40,55};


//	float BPRAAXErrUp[NBins] = {1.27,2.6,2.8,22.7};
//	float BPRAAXErrDown[NBins] = {1.73,2.4,2.2,7.3};

	float BPRAAXErrUp[NBins] = {1,1.27,2.6,2.8,5,10,5};
	float BPRAAXErrDown[NBins] = {1,1.23,2.4,2.2,5,10,5};


	float BPRAAYSystUp[NBins] ;
	float BPRAAYSystDown[NBins];
	float BPRAAYSystUpRatio[NBins] ;
	float BPRAAYSystDownRatio[NBins];



	float BPRAAYErrUp[NBins];
	float BPRAAYErrDown[NBins];
	float BPRAAYErrUpRatio[NBins];
	float BPRAAYErrDownRatio[NBins];


	//2015 Data Sets
	
	const int NBins2015 = 5;

	float BPRAAY2015[NBins2015] = {0.35,0.448,0.440,0.615,0.35};
	float BPRAAX2015[NBins2015] = {8.5,12.5,17.5,25,40};




	float BPRAAXErrUp2015[NBins2015] = {1.5,2.5,2.5,5,10};
	float BPRAAXErrDown2015[NBins2015] = {1.5,2.5,2.5,5,10};

	float BPRAAYErrUp2015[NBins2015] = {0.11,0.074,0.075,0.092,0.11};
	float BPRAAYErrDown2015[NBins2015] = {0.11,0.074,0.075,0.092,0.11};

	float BPRAAYSystUp2015[NBins2015] = {0.064,0.077,0.074,0.102,0.0539};
	float BPRAAYSystDown2015[NBins2015] = {0.064,0.077,0.074,0.102,0.0539};



	for(int i = 0; i < NBins; i++){

		BPRAAY[i] = BPXsecPbPbY[i] / BPXsecPPY[i];


		BPRAAYErrUp[i] = BPRAAY[i] * TMath::Sqrt(BPXSecPbPbYErrUpRatio[i] * BPXSecPbPbYErrUpRatio[i] + BPXSecPPYErrUpRatio[i] * BPXSecPPYErrUpRatio[i]);
		BPRAAYErrDown[i] = BPRAAY[i] * TMath::Sqrt(BPXSecPbPbYErrDownRatio[i] * BPXSecPbPbYErrDownRatio[i] + BPXSecPPYErrDownRatio[i] * BPXSecPPYErrDownRatio[i]);
		

		BPRAAYSystDown[i] = BPRAAY[i] * TMath::Sqrt(BPXSecPbPbYSystDownRatio[i] * BPXSecPbPbYSystDownRatio[i] + BPTotalSystDownRatio[i] * BPTotalSystDownRatio[i]);
		BPRAAYSystUp[i] = BPRAAY[i] * TMath::Sqrt(BPXSecPbPbYSystUpRatio[i] * BPXSecPbPbYSystUpRatio[i] + BPTotalSystUpRatio[i] * BPTotalSystUpRatio[i]);

    BPRAAYErrUpRatio[i] = BPRAAYErrUp[i]/BPRAAY[i];
    BPRAAYErrDownRatio[i] = BPRAAYErrDown[i]/BPRAAY[i];
    BPRAAYSystUpRatio[i] = BPRAAYSystUp[i]/BPRAAY[i];
    BPRAAYSystDownRatio[i] = BPRAAYSystDown[i]/BPRAAY[i];

		cout << "i = " << i <<  "    BPRAAY[i] = " << BPRAAY[i] <<
      "   BPRAAYErrUpRatio[i]  = "  <<  BPRAAYErrUpRatio[i] <<
      "  BPRAAYErrDownRatio[i] =   " <<  BPRAAYErrDownRatio[i] <<
      "   BPRAAYSystUpRatio[i]  = "  << BPRAAYSystUpRatio[i]   <<
      "  BPRAAYSystDownRatio[i] =   " << BPRAAYSystDownRatio[i]  << endl;


	}






	TGraphAsymmErrors *BPRAAGraph = new TGraphAsymmErrors(NBins, BPRAAX, BPRAAY,BPRAAXErrDown, BPRAAXErrUp,BPRAAYErrDown,BPRAAYErrUp);
	BPRAAGraph->SetName("BPRAAGraph");
	TGraphAsymmErrors *BPRAAGraphSyst = new TGraphAsymmErrors(NBins, BPRAAX, BPRAAY,BPRAAXErrDown, BPRAAXErrUp,BPRAAYSystDown,BPRAAYSystUp);



	TGraphAsymmErrors *BPRAAGraph2015 = new TGraphAsymmErrors(NBins2015, BPRAAX2015, BPRAAY2015,BPRAAXErrDown2015, BPRAAXErrUp2015,BPRAAYErrDown2015,BPRAAYErrUp2015);
	TGraphAsymmErrors *BPRAAGraphSyst2015 = new TGraphAsymmErrors(NBins2015, BPRAAX2015, BPRAAY2015,BPRAAXErrDown2015, BPRAAXErrUp2015,BPRAAYSystDown2015,BPRAAYSystUp2015);


	BPRAAGraph->SetLineColor(kRed+2);
//	BPRAAGraph->SetFillColorAlpha(kRed+2,0.5);
	BPRAAGraph->SetMarkerStyle(20);
	BPRAAGraph->SetMarkerSize(1);
	BPRAAGraph->SetMarkerColor(kRed+2);




	BPRAAGraph2015->SetLineColor(kBlue+2);
//	BPRAAGraph->SetFillColorAlpha(kRed+2,0.5);
	BPRAAGraph2015->SetMarkerStyle(21);
	BPRAAGraph2015->SetMarkerSize(1);
	BPRAAGraph2015->SetMarkerColor(kBlue+2);
	

	BPRAAGraphSyst->SetFillColorAlpha(kRed-9,0.5);
	BPRAAGraphSyst->SetLineColor(kRed-9);
	BPRAAGraphSyst2015->SetFillColorAlpha(kBlue-9,0.5);
	BPRAAGraphSyst2015->SetLineColor(kBlue-9);





	TLine * Unity = new TLine(5,1,60,1);
	Unity->SetLineWidth(2);
	Unity->SetLineStyle(2);
	Unity->SetLineColor(1);

	HisEmptyRAA->Draw();
	BPRAAGraph->Draw("ep");
	BPRAAGraphSyst->Draw("5same");
	Unity->Draw("SAME");
	
	c2->SaveAs("RAAPlots/BP/BPRAA.png");
	c2->SaveAs("RAAPlots/BP/BPRAA.pdf");


	HisEmptyRAA->Draw();
	BPRAAGraph->Draw("ep");
	BPRAAGraph->Draw("epSAME");
	BPRAAGraph2015->Draw("epSAME");
	BPRAAGraphSyst->Draw("5same");
	BPRAAGraphSyst2015->Draw("5same");


	TLegend* leg2 = new TLegend(0.30,0.75,0.60,0.90,NULL,"brNDC");
	leg2->SetBorderSize(0);
	leg2->SetTextSize(0.040);
	leg2->SetTextFont(42);
	leg2->SetFillStyle(0);
	leg2->SetLineWidth(3);
	leg2->AddEntry(BPRAAGraph,"2018 PbPb + 2017 pp","PL");
	leg2->AddEntry(BPRAAGraph2015,"2015 PbPb + 2015 pp","PL");
	leg2->Draw("same");



	Unity->Draw("SAME");




	c2->SaveAs("RAAPlots/BP/BPRAACompairson.png");
	c2->SaveAs("RAAPlots/BP/BPRAACompairson.pdf");
	

	TFile * fout = new TFile("OutFile/BPRAA.root","RECREATE");

	fout->cd();

	BPRAAGraph->Write();

	fout->Close();

  std::vector<float> globUncert(NBins, 0.035);
  // summary of errors (in ratio, not percent)
  std::vector<int> ptbins = {5, 7, 10, 15, 20, 30, 50, 60};
  std::vector<float> abscissae = {6.0, 8.5, 12.5, 17.5, 25, 40, 55};

  string outFile = "../../MakeFinalPlots/NominalPlots/RAA/dataSource/RAA_pt_Bp_New.txt";
  ofstream out;
  out.open(outFile);
  unsigned columnWidth = 14;
  out << std::left << std::setw(columnWidth) <<
    "ptmin" << std::setw(columnWidth) << "ptmax" << std::setw(columnWidth) <<
    "central_val" << std::setw(columnWidth) <<
    "statUp" << std::setw(columnWidth) << "statDown" << std::setw(columnWidth) <<
    "systUp" << std::setw(columnWidth) << "systDown" << std::setw(columnWidth) <<
    "glbUp" << std::setw(columnWidth) << "glbDown" << std::setw(columnWidth) <<
    "abscissae" << endl;
  for (auto i = 0; i < NBins; ++i ) {
    out << std::setw(columnWidth) <<
      ptbins[i] << std::setw(columnWidth) << ptbins[i + 1] << std::setw(columnWidth) <<
      setprecision(3) << BPRAAY[i] << std::setw(columnWidth) <<
      setprecision(2) << std::defaultfloat <<
      BPRAAYErrUpRatio[i] << std::setw(columnWidth) <<
      BPRAAYErrDownRatio[i] << std::setw(columnWidth) <<
      BPRAAYSystUpRatio[i] << std::setw(columnWidth) <<
      BPRAAYSystDownRatio[i] << std::setw(columnWidth) <<
      globUncert[i] << std::setw(columnWidth) <<
      globUncert[i] << std::setw(columnWidth) <<
      setprecision(3) << abscissae[i] << "\n";
  }
  out.close();
}

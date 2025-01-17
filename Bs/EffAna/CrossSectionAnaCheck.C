#include "TROOT.h"
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
#include <fstream>
//#include "tnp_weight_lowptPbPb.h"



//#include "his.h"
using namespace std;

using std::cout;
using std::endl;

void CrossSectionAnaCheck(int DoTnP, float Factor){

	const int NBins = 4;
	//const int NBins = 6;

	int TnP = 1;
	

//	double BRchain = 6.02061e-5;
	double BRchain = 3.1189e-5;

//	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);


	TString FileName;

	FileName = "../../SkimmedSamples/BsData.root";
	TFile * fin = new TFile(FileName.Data());
	fin->cd();

	TTree * EffInfoTree = (TTree * ) fin->Get("ntphi");

	int NEvents = EffInfoTree->GetEntries();

	const int NCand = 15;

	Int_t BsizeNew;
	Int_t runNew;
	Int_t lumiNew;
	Int_t evtNew;
	Float_t BmassNew[NCand];
	Float_t BptNew[NCand];
	Float_t ByNew[NCand];
	Float_t BEff[NCand];
	Float_t BEffErr[NCand];


	Float_t BEffInv[NCand];
	Float_t BEffInvErr[NCand];
	Float_t BEffInv1D[NCand];
	Float_t BEffInvErr1D[NCand];

	Float_t BEffInvFit[NCand];
	Float_t BEffInvErrFit[NCand];

	Float_t BEffInvBDTWeighted[NCand];
	Float_t BEffInvErrBDTWeighted[NCand];



	Float_t BEffInvUp[NCand];
	Float_t BEffInvErrUp[NCand];
	Float_t BEffInvDown[NCand];
	Float_t BEffInvErrDown[NCand];

	EffInfoTree->SetBranchAddress("Bsize",&BsizeNew);
	EffInfoTree->SetBranchAddress("Bmass",BmassNew);
	EffInfoTree->SetBranchAddress("By",ByNew);
	EffInfoTree->SetBranchAddress("Bpt",BptNew);
	

/*
	EffInfoTree->SetBranchAddress("BEffInv",BEffInv);
	EffInfoTree->SetBranchAddress("BEffInvErr",BEffInvErr);
	EffInfoTree->SetBranchAddress("BEff",BEff);
	EffInfoTree->SetBranchAddress("BEffInv1D",BEffInv1D);
	EffInfoTree->SetBranchAddress("BEffInvErr1D",BEffInvErr1D);

	EffInfoTree->SetBranchAddress("BEffInvFit",BEffInvFit);
	EffInfoTree->SetBranchAddress("BEffInvErrFit",BEffInvErrFit);

	EffInfoTree->SetBranchAddress("BEffInv",BEffInv);
	EffInfoTree->SetBranchAddress("BEffInvErr",BEffInvErr);

	EffInfoTree->SetBranchAddress("BEffInvBDTWeighted",BEffInvBDTWeighted);
	EffInfoTree->SetBranchAddress("BEffInvErrBDTWeighted",BEffInvErrBDTWeighted);



	EffInfoTree->SetBranchAddress("BEffInvUp",BEffInvUp);
	EffInfoTree->SetBranchAddress("BEffInvErrUp",BEffInvErrUp);

	EffInfoTree->SetBranchAddress("BEffInvDown",BEffInvDown);
	EffInfoTree->SetBranchAddress("BEffInvErrDown",BEffInvErrDown);
*/

	/*
	   Int_t Bmu1Type[NCand];
	   Int_t Bmu2Type[NCand];





	   TTree * MuonInfoTree = (TTree * ) fin->Get("MuonInfoTree");

	   MuonInfoTree->SetBranchAddress("Bmu1Type",Bmu1Type);
	   MuonInfoTree->SetBranchAddress("Bmu2Type",Bmu2Type);


	   Float_t Bmu1etaNew[NCand];
	   Float_t Bmu2etaNew[NCand];

	   Float_t Bmu1ptNew[NCand];
	   Float_t Bmu2ptNew[NCand];


	   EffInfoTree->SetBranchAddress("Bmu1eta",Bmu1etaNew);
	   EffInfoTree->SetBranchAddress("Bmu2eta",Bmu2etaNew);

	   EffInfoTree->SetBranchAddress("Bmu1pt",Bmu1ptNew);
	   EffInfoTree->SetBranchAddress("Bmu2pt",Bmu2ptNew);

*/


	double ptBins[NBins + 1];

	int Counts[NBins];
	double SumCounts[NBins];
	double SumCountsErr[NBins];
	double NewEff[NBins];
	double NewEffErr[NBins];

	double NewEffReal[NBins];
	double NewEffRealErr[NBins];



	double SumCountsUp[NBins];
	double SumCountsErrUp[NBins];

	double SumCountsDown[NBins];
	double SumCountsErrDown[NBins];


	double SumCountsEff[NBins];
	double SumCountsEffErr[NBins];



	double SumCountsSyst[NBins];
	double SumCountsSystErr[NBins];
	double NewEffSyst[NBins];
	double NewEffSystErr[NBins];


	//	double CorrectionFactor[NBins];

	double NewEffUp[NBins];
	double NewEffErrUp[NBins];
	double NewEffDown[NBins];
	double NewEffErrDown[NBins];


	double lumi = 302.3;
	std::vector<double> ptbinsvec;
	std::vector<double> corrfactvec;


	if(NBins == 1){


		ptbinsvec.push_back(10.0);
		ptbinsvec.push_back(50);
	}


	if(NBins == 3){

		ptbinsvec.push_back(5);
		ptbinsvec.push_back(15);
		ptbinsvec.push_back(20);
		ptbinsvec.push_back(50);

	}


	if(NBins == 4){

		ptbinsvec.push_back(7);
		ptbinsvec.push_back(10);
		ptbinsvec.push_back(15);
		ptbinsvec.push_back(20);
		ptbinsvec.push_back(50);

		/*
		   corrfactvec.push_back(1.24759);
		   corrfactvec.push_back(1.05256);
		   corrfactvec.push_back(1.02614);
		   corrfactvec.push_back(1.01174);
		   */


	}

	if(NBins == 6){


		ptbinsvec.push_back(5);
		ptbinsvec.push_back(7);		
		ptbinsvec.push_back(10);
		ptbinsvec.push_back(15);
		ptbinsvec.push_back(20);
		ptbinsvec.push_back(50);
		ptbinsvec.push_back(100);



	}



	if(NBins == 7){


		//ptbinsvec.push_back(3);
		ptbinsvec.push_back(5);
		ptbinsvec.push_back(7);		
		ptbinsvec.push_back(10);
		ptbinsvec.push_back(15);
		ptbinsvec.push_back(20);
		ptbinsvec.push_back(30);		

		ptbinsvec.push_back(50);
	//	ptbinsvec.push_back(100);
		ptbinsvec.push_back(60);



	}



	if(NBins == 9){

		ptbinsvec.push_back(2);
		ptbinsvec.push_back(3);
		ptbinsvec.push_back(5);
		ptbinsvec.push_back(7);		
		ptbinsvec.push_back(10);
		ptbinsvec.push_back(15);
		ptbinsvec.push_back(20);
		ptbinsvec.push_back(30);
		ptbinsvec.push_back(50);
		ptbinsvec.push_back(100);



	}

	if(NBins == 10){

		ptbinsvec.push_back(1);
		ptbinsvec.push_back(2);
		ptbinsvec.push_back(3);
		ptbinsvec.push_back(5);
		ptbinsvec.push_back(7);		
		ptbinsvec.push_back(10);
		ptbinsvec.push_back(15);
		ptbinsvec.push_back(20);
		ptbinsvec.push_back(30);
		ptbinsvec.push_back(50);
		ptbinsvec.push_back(100);



	}





	for(int i = 0; i < NBins + 1; i++){
		ptBins[i] =  ptbinsvec[i];
	}

	for(int i = 0; i < NBins; i++){
		Counts[i] = 0;
		SumCounts[i] = 0;
		SumCountsErr[i] = 0;
		SumCountsEff[i] = 0;
		SumCountsEffErr[i] = 0;
		SumCountsSyst[i] = 0;
		SumCountsSystErr[i] = 0;
		//	CorrectionFactor[i] = corrfactvec[i];
		SumCountsUp[i] = 0;
		SumCountsErrUp[i] = 0;
		SumCountsDown[i] = 0;
		SumCountsErrDown[i] = 0;




	}



	/*
	   const int NBins = 3;
	   double ptBins[NBins + 1] ={5,15,20,50};

	   int Counts[NBins]={0,0,0};
	   double SumCounts[NBins]={0,0,0};
	   double SumCountsErr[NBins]={0,0,0};
	   */


	/*
	   const int NBins = 4;
	   double ptBins[NBins + 1] ={5,10,15,20,50};

	   int Counts[NBins]={0,0,0,0};
	   double SumCounts[NBins]={0,0,0,0};
	   double SumCountsErr[NBins]={0,0,0,0};
	   */



	int EtaBin;
	int PtBin;




	double trgtnp1;
	double trktnp1;
	double muidtnp1;

	double trgtnp1systup;
	double trgtnp1systdown;
	double trgtnp1statup;
	double trgtnp1statdown;


	double trktnp1systup;
	double trktnp1systdown;
	double trktnp1statup;
	double trktnp1statdown;

	double muidtnp1systup;
	double muidtnp1systdown;
	double muidtnp1statup;
	double muidtnp1statdown;


	double tnptotal1;
	double tnptotal1up;
	double tnptotal1down;


	double tnptotal1systup;
	double tnptotal1systdown;
	double tnptotal1statup;
	double tnptotal1statdown;



	double trgtnp2;
	double trktnp2;
	double muidtnp2;

	double trgtnp2systup;
	double trgtnp2systdown;
	double trgtnp2statup;
	double trgtnp2statdown;


	double trktnp2systup;
	double trktnp2systdown;
	double trktnp2statup;
	double trktnp2statdown;

	double muidtnp2systup;
	double muidtnp2systdown;
	double muidtnp2statup;
	double muidtnp2statdown;


	double tnptotal2;
	double tnptotal2up;
	double tnptotal2down;

	double tnptotal2systup;
	double tnptotal2systdown;
	double tnptotal2statup;
	double tnptotal2statdown;






	double tnpabssystup;
	double tnpabssystdown;

	TFile * fin1DEff;

	if(DoTnP == 0) fin1DEff = new TFile("NewEff2DMaps/EffFineNoTnP.root");
//	if(DoTnP == 1) fin1DEff = new TFile("NewEff2DMaps/EffFineBDT.root");
	if(DoTnP == 1) fin1DEff = new TFile(Form("NewEff2DMaps/EffFineBDT_%.0f.root",Factor));

	fin1DEff->cd();

	TH2D * invEff2D = (TH2D *) fin1DEff->Get("invEff2D");

	int XBin;
	int YBin;

	for( int i = 0; i < NEvents; i++){

		EffInfoTree->GetEntry(i);
		//MuonInfoTree->GetEntry(i);


		for(int j = 0; j < BsizeNew; j++){


			for(int k = 0; k < NBins; k++){

				//	if((BptNew[j] > ptBins[k] && BptNew[j] < ptBins[k+1] && TMath::Abs(BmassNew[j] - 5.27932) < 0.08  && ((BptNew[j] > 7 && BptNew[j] < 10 && ByNew[j] > 1.5 )||(BptNew[j] > 10)) && (Bmu1Type > -0.1 && Bmu2Type > -0.1)))
				if(BptNew[j] > ptBins[k] && BptNew[j] < ptBins[k+1] && TMath::Abs(BmassNew[j] - 5.3663) < 0.08 && TMath::Abs(ByNew[j]) < 2.4 && ((BptNew[j] > 5 && BptNew[j] < 10 && abs(ByNew[j]) > 1.5 )||(BptNew[j] > 10)) )
				{



					XBin = invEff2D->GetXaxis()->FindBin( BptNew[j]);
					YBin = invEff2D->GetYaxis()->FindBin( TMath::Abs(ByNew[j]));
					BEffInv[j] = invEff2D->GetBinContent(XBin,YBin);
					BEffInvErr[j] = invEff2D->GetBinError(XBin,YBin);
					BEff[j] = 1.0/invEff2D->GetBinContent(XBin,YBin);

					//if(BptNew[j] < 10) cout  << "BptNew[j] = " << BptNew[j] << "  TMath::Abs(ByNew[j])  = " << TMath::Abs(ByNew[j])   << 	"   XBin = " << XBin << "  YBin = " << YBin << endl;


					BEffErr[j] = BEffInvErr[j]/(BEffInv[j] * BEffInv[j]);
					if(BEffInv[j] > 0){
						SumCounts[k] = SumCounts[k] + BEffInv[j];
						SumCountsErr[k] = SumCountsErr[k] + BEffInvErr[j] * BEffInvErr[j];
						SumCountsEff[k] = SumCountsEff[k] + BEff[j];
						SumCountsEffErr[k] = SumCountsEffErr[k] + BEffErr[j] * BEffErr[j];
						SumCountsSyst[k] = 	SumCountsSyst[k]  + BEffInvBDTWeighted[j];
						SumCountsSystErr[k] = 	SumCountsSystErr[k]  + BEffInvErrBDTWeighted[j] * BEffInvErrBDTWeighted[j];

						SumCountsUp[k] = SumCountsUp[k] + BEffInvUp[j];
						SumCountsErrUp[k] = SumCountsErrUp[k] + BEffInvErrUp[j] * BEffInvErrUp[j];

						SumCountsDown[k] = SumCountsDown[k] + BEffInvDown[j];
						SumCountsErrDown[k] = SumCountsErrUp[k] + BEffInvErrDown[j] * BEffInvErrDown[j];

						Counts[k] = Counts[k] + 1;

						//cout << "SumCounts = " << SumCounts[k] << endl;

						//cout << "SumCountsUp = " << SumCountsUp[k] << endl;
						//cout << "Candidate: " << "   Bpt = " << BptNew[j] << "   Efficiency =  " << BEffInv[j] << "   Efficiency Error  = " << BEffInvErr[j] << endl;
					}
				}

			}


		}

	}


	TH1D * hInvEff = new TH1D("hInvEff","",NBins,ptBins);


	hInvEff->GetXaxis()->SetTitle("B^{0}_{s} p_{T} (GeV/c)");
	hInvEff->GetYaxis()->SetTitle("<1/(Eff * Acc)>");
	hInvEff->GetYaxis()->SetTitleOffset(1.4);
	hInvEff->GetXaxis()->CenterTitle();
	hInvEff->GetYaxis()->CenterTitle();
	hInvEff->SetMarkerColor(1);
	hInvEff->SetLineColor(1);
	hInvEff->SetMarkerStyle(20);

	hInvEff->SetMinimum(0);


	TH1D * hInvEffSyst = new TH1D("hInvEffSyst","",NBins,ptBins);

	hInvEffSyst->GetXaxis()->SetTitle("B^{0}_{s} p_{T} (GeV/c)");
	hInvEffSyst->GetYaxis()->SetTitle("<1/(Eff * Acc)> - BDT Data-MC Weighted");
	hInvEffSyst->GetYaxis()->SetTitleOffset(1.4);
	hInvEffSyst->GetXaxis()->CenterTitle();
	hInvEffSyst->GetYaxis()->CenterTitle();
	hInvEffSyst->SetMarkerColor(1);
	hInvEffSyst->SetLineColor(2);
	hInvEffSyst->SetMarkerStyle(20);

	hInvEffSyst->SetMinimum(0);

	TH1D * hEff = new TH1D("hEff","",NBins,ptBins);


	hEff->GetXaxis()->SetTitle("B^{0}_{s} p_{T} (GeV/c)");
	hEff->GetYaxis()->SetTitle("<(Eff * Acc)>");
	hEff->GetYaxis()->SetTitleOffset(1.4);
	hEff->GetXaxis()->CenterTitle();
	hEff->GetYaxis()->CenterTitle();
	hEff->SetMarkerColor(1);
	hEff->SetLineColor(1);
	hEff->SetMarkerStyle(20);

	hEff->SetMinimum(0);


	TH1D * hInvEffUp = new TH1D("hInvEffUp","",NBins,ptBins);


	hInvEffUp->GetXaxis()->SetTitle("B^{0}_{s} p_{T} (GeV/c)");
	hInvEffUp->GetYaxis()->SetTitle("<1./(Eff * Acc)>");
	hInvEffUp->GetYaxis()->SetTitleOffset(1.4);
	hInvEffUp->GetXaxis()->CenterTitle();
	hInvEffUp->GetYaxis()->CenterTitle();
	hInvEffUp->SetMarkerColor(1);
	hInvEffUp->SetLineColor(1);
	hInvEffUp->SetMarkerStyle(20);
	hInvEffUp->SetMinimum(0);



	TH1D * hInvEffDown = new TH1D("hInvEffDown","",NBins,ptBins);


	hInvEffDown->GetXaxis()->SetTitle("B^{0}_{s} p_{T} (GeV/c)");
	hInvEffDown->GetYaxis()->SetTitle("<1/(Eff * Acc)>");
	hInvEffDown->GetYaxis()->SetTitleOffset(1.4);
	hInvEffDown->GetXaxis()->CenterTitle();
	hInvEffDown->GetYaxis()->CenterTitle();
	hInvEffDown->SetMarkerColor(1);
	hInvEffDown->SetLineColor(1);
	hInvEffDown->SetMarkerStyle(20);
	hInvEffDown->SetMinimum(0);



	for(int i = 0; i < NBins; i++){


		NewEff[i] = SumCounts[i]/Counts[i];
		NewEffErr[i] = TMath::Sqrt(SumCountsErr[i])/Counts[i];


		NewEffUp[i] = SumCountsUp[i]/Counts[i];
		NewEffErrUp[i] = TMath::Sqrt(SumCountsErrUp[i])/Counts[i];



		NewEffDown[i] = SumCountsDown[i]/Counts[i];
		NewEffErrDown[i] = TMath::Sqrt(SumCountsErrDown[i])/Counts[i];


		NewEffReal[i] = SumCountsEff[i]/Counts[i];
		NewEffRealErr[i] = TMath::Sqrt(SumCountsEffErr[i])/Counts[i];


		hInvEff->SetBinContent(i+1,NewEff[i]);
		hInvEff->SetBinError(i+1,NewEffErr[i]);

		hEff->SetBinContent(i+1,1/NewEff[i]);
		hEff->SetBinError(i+1,NewEffErr[i]/(NewEff[i] * NewEff[i]));


		NewEffSyst[i] = SumCountsSyst[i]/Counts[i];
		NewEffSystErr[i] = TMath::Sqrt(SumCountsSystErr[i])/Counts[i];


		hInvEffSyst->SetBinContent(i+1,	NewEffSyst[i]);
		hInvEffSyst->SetBinError(i+1, NewEffSystErr[i]);



		hInvEffUp->SetBinContent(i+1,NewEffUp[i]);
		hInvEffUp->SetBinError(i+1,NewEffErrUp[i]);


		hInvEffDown->SetBinContent(i+1,NewEffDown[i]);
		hInvEffDown->SetBinError(i+1,NewEffErrDown[i]);

		//	cout << "Real eff = " << SumCountsReal[i]/Counts[i] << endl;
		//cout << "Counts = " << Counts[i] << endl;
		cout << "Count =  " <<  Counts[i] << "   NewEff = " << NewEff[i] << "     NewEffErr = " << NewEffErr[i] << endl;
		cout << "Count =  " <<  Counts[i] << "   NewEffSyst = " << NewEffSyst[i] << "     NewEffSystErr = " << NewEffSystErr[i] << endl;



		cout << "-----------------------------------------------------------------------------------------------" << endl;

		cout << "   NewEff = " << NewEff[i] << "     NewEffErr = " << NewEffErr[i] << "  Fractional = " << NewEffErr[i]/NewEff[i] << endl;
		//	cout << "   NewEff = " << NewEffUp[i] << "     NewEffErr = " << NewEffErrUp[i] << "  Fractional = " << NewEffErrUp[i]/NewEffUp[i] << endl;



		//NewEffErr[i] = 0; //Remove Error on Efficiency Correction//
	}


	for(int i = 0 ; i < NBins; i++){

		cout << "--------------------------------------------------  Eff Systematics Uncertainties  -----------------------------------------------" << endl;


		cout << "i = " << i << "    NewEffUp =   "  << NewEffUp[i]  <<    "    NewEffDown =   " <<  NewEffDown[i] << endl;

		cout << "i = " << i << "    SystUp =   "  << (NewEffUp[i] - NewEff[i])/ NewEff[i] <<    "    NewEffDown =   " << (NewEff[i] - NewEffDown[i])/NewEff[i] << endl;
		cout << "-----------------------------------------------------------------------------------------------------------------------------------" << endl;
	}

	//return;

	hInvEff->SetMaximum(NewEff[0]*1.5);
	TCanvas *c = new TCanvas("c","c",600,600);

	c->cd();

	hInvEff->Draw("ep");

	cout << "OK" << endl;

	c->SaveAs(Form("EffFinal/ReAnaEff_%dBins.png",NBins));
	c->SaveAs(Form("EffFinal/pdf/ReAnaEff_%dBins.pdf",NBins));

	hEff->Draw("ep");

	c->SaveAs(Form("EffFinal/ReAnaEffReal_%dBins.png",NBins));
	c->SaveAs(Form("EffFinal/pdf/ReAnaEffReal_%dBins.pdf",NBins));


	/*

	   hInvEff->SetMarkerColor(2);
	   hInvEff->SetLineColor(2);
	   hInvEffSyst->SetMarkerColor(3);
	   hInvEffSyst->SetLineColor(3);



	   TLegend *leg = new TLegend(0.19,0.60,0.39,0.87,NULL,"brNDC");
	   leg->SetBorderSize(0);
	   leg->SetTextSize(0.04);
	   leg->SetTextFont(42);
	   leg->SetFillStyle(0);

	   leg->AddEntry(hInvEff,"Nominal Correction","pl");
	   leg->AddEntry(hInvEffSyst,"BDT MC-Data Weighted Correction","pl");

	   hInvEffSyst->SetMinimum(0);

	   hInvEffSyst->Draw("ep");
	   hInvEff->Draw("epSAME");
	   leg->Draw("SAME");

	   c->SaveAs(Form("EBDTWeightedComp_%dBins.png",NBins));

	   TH1D * SelEffSystRatio =  (TH1D * ) hInvEffSyst->Clone("SelEffSystRatio");
	   SelEffSystRatio->GetYaxis()->SetTitle("Syst Variation/Nominal");
	   SelEffSystRatio->Sumw2();
	   hInvEff->Sumw2();
	   SelEffSystRatio->Divide(hInvEff);

	   SelEffSystRatio->SetMaximum(3);

	   SelEffSystRatio->SetMinimum(0);


	   TLatex	*texChi = new TLatex(0.20,0.95, "BDT Weighted/Nominal");
	   texChi->SetNDC();
	   texChi->SetTextAlign(12);
	   texChi->SetTextSize(0.04);
	   texChi->SetTextFont(42);
	   texChi->SetTextColor(1);

	   SelEffSystRatio->GetYaxis()->SetTitleOffset(1.2);
	   SelEffSystRatio->SetMarkerColor(1);
	   SelEffSystRatio->SetLineColor(1);
	   SelEffSystRatio->Draw("ep");



	   TLine *l5 = new TLine(7,1,50,1);
	   l5->SetLineStyle(2);
	   l5->SetLineWidth(2);
	   l5->SetLineColor(2);

	   l5->Draw("SAME");
	   texChi->Draw("SAME");
	   c->SaveAs(Form("SystEffRatio_%dBins.png",NBins));


	//TnP Comparison//

	hInvEff->SetMinimum(5);
	hInvEff->SetMaximum(25);



	hInvEff->SetMarkerColor(1);
	hInvEff->SetLineColor(1);
	hInvEff->Draw("ep");


	hInvEffUp->SetMarkerColor(2);
	hInvEffUp->SetLineColor(2);
	hInvEffDown->SetMarkerColor(3);
	hInvEffDown->SetLineColor(3../RawYieldFits/ROOTfiles/);


	hInvEffUp->Draw("epSAME");
	hInvEffDown->Draw("epSAME");

	TLegend *legTnP = new TLegend(0.17,0.20,0.44,0.47,NULL,"brNDC");
	legTnP->SetBorderSize(0);
	legTnP->SetTextSize(0.04);
	legTnP->SetTextFont(42);
	legTnP->SetFillStyle(0);

	legTnP->AddEntry(hInvEff,"TnP - Total: Nominal Correction","pl");
	legTnP->AddEntry(hInvEffUp,"TnP - Total: Upper Bound of Scale Factor","pl");
	legTnP->AddEntry(hInvEffDown,"TnP - Total: Lower Bound of Scale Factor","pl");
	legTnP->Draw("SAME");



	c->SaveAs(Form("CheckSystNuno/%s/TnPFinal/TnPEffComp_%dBins_%d_%d.png",WeightName.Data(),NBins,CentMin,CentMax));



	TH1D * hInvUpRatio =	(TH1D *) hInvEffUp->Clone("hInvUpRatio");
	TH1D * hInvDownRatio =	(TH1D *) hInvEffDown->Clone("hInvDownRatio");

	hInvUpRatio->Sumw2();
	hInvEff->Sumw2();
	hInvUpRatio->Divide(hInvEff);

	hInvDownRatio->Sumw2();
	hInvEff->Sumw2();
	hInvDownRatio->Divide(hInvEff);

	hInvUpRatio->Draw("ep");
	hInvDownRatio->Draw("epSAME");


	TLegend *legTnPRatio = new TLegend(0.49,0.30,0.69,0.45,NULL,"brNDC");
	legTnPRatio->SetBorderSize(0);
	legTnPRatio->SetTextSize(0.04);
	legTnPRatio->SetTextFont(42);
	legTnPRatio->SetFillStyle(0);

	legTnPRatio->AddEntry(hInvUpRatio,"TnP - Total Up/Nominal","pl");
	legTnPRatio->AddEntry(hInvDownRatio,"TnP - Total Down/Nominal","pl");
	legTnPRatio->Draw("SAME");

	TLine *l6 = new TLine(7,1,50,1);
	l6->SetLineStyle(2);
	l6->SetLineWidth(1);
	l6->SetLineColor(1);

	l6->Draw("SAME");
CorrDiffHisBin


	c->SaveAs(Form("CheckSystNuno/%s/TnPFinal/SystTnPRatio_%dBins_%d_%d.png",WeightName.Data(),NBins,CentMin,CentMax));



	for(int i = 0; i < hInvEffUp->GetNbinsX();i++){

		cout << "i = " << i << "   Upper TnP Syst = " << hInvUpRatio->GetBinContent(i+1) - 1 << endl;
		cout << "i = " << i << "   Lower TnP Syst = " << 1 - hInvDownRatio->GetBinContent(i+1)  << endl;

	}

	*/

		//TnP Syst DONE//

	TFile * RawYield = new TFile("../RawYieldFits/ROOTfiles/yields_Bs_binned_pt.root");
	RawYield->cd();
	TH1D * hPt = (TH1D *) RawYield->Get("hPt");


	double RawCount;
	double RawCountErr;

	double CorrYield= 0;
	double CorrYieldErr= 0;

	double CorrYieldDiff[NBins];
	double CorrYieldDiffErr[NBins];

	for(int i = 0; i < NBins;i++){

		RawCount = hPt->GetBinContent(i+1);
		RawCountErr = hPt->GetBinError(i+1);

	//	cout << "RawCount = " << RawCount << "  RawCountErr = " << RawCountErr << " NewEff[i] =   " << NewEff[i] << "  NewEffErr[i] =  " << NewEffErr[i] << endl; 

		cout << "CORR YIELD PT:  " <<  RawCount *   NewEff[i]/(BRchain*2 * lumi) << endl;
		CorrYield = RawCount * (ptBins[i+1] - ptBins[i]) *  NewEff[i]  + CorrYield;
		CorrYieldErr = ((RawCountErr * (ptBins[i+1] - ptBins[i]) *  NewEff[i]) *(RawCountErr * (ptBins[i+1] - ptBins[i]) *  NewEff[i]) + (RawCount * (ptBins[i+1] - ptBins[i]) *  NewEffErr[i]) * (RawCount * (ptBins[i+1] - ptBins[i]) *  NewEffErr[i]))  + CorrYieldErr;

		cout << "PrintYield = " << RawCount* (ptBins[i+1] - ptBins[i]) *  NewEff[i] << endl;

	}



	TFile * foutCorr;
	if(DoTnP == 0)	foutCorr = new TFile("FinalFiles/BsPPCorrYieldPTNoTnP.root","RECREATE");
//	if(DoTnP == 1)	foutCorr = new TFile("FinalFiles/BsPPCorrYieldPT.root","RECREATE");
	if(DoTnP == 1)	foutCorr = new TFile(Form("FinalFiles/BsPPCorrYieldPT_%.0f.root",Factor),"RECREATE");


	TH1D * CorrDiffHis = new TH1D("hPtSigma","",NBins,ptBins);
	CorrDiffHis->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	CorrDiffHis->GetYaxis()->SetTitle("d #sigma/d p_{T} (pb GeV^{-1} c)");

	CorrDiffHis->GetYaxis()->SetTitleOffset(1.3);
	CorrDiffHis->GetXaxis()->CenterTitle();
	CorrDiffHis->GetYaxis()->CenterTitle();



	for(int i = 0; i < NBins;i++){
		RawCount = hPt->GetBinContent(i+1);
		RawCountErr = hPt->GetBinError(i+1);
		CorrYieldDiff[i] = (RawCount *  NewEff[i])/(BRchain*2* lumi);
		CorrYieldDiffErr[i] = TMath::Sqrt((RawCountErr *  NewEff[i]) *(RawCountErr  *  NewEff[i]) + (RawCount *  NewEffErr[i]) * (RawCount  *  NewEffErr[i]))/(BRchain*2* lumi);
		CorrDiffHis->SetBinContent(i+1,CorrYieldDiff[i]);
		CorrDiffHis->SetBinError(i+1,CorrYieldDiffErr[i]);

	}


//	CorrDiffHis->SetTitle("2017 B^{0}_{s} pp Cross Section");
	CorrDiffHis->SetTitle("(Preliminary) B^{0}_{s} #rightarrow J/#psi K^{+} p_{T} Differential Cross Section in pp");

	CorrDiffHis->SetMarkerColor(kBlack);
	CorrDiffHis->SetMarkerSize(1);
	CorrDiffHis->SetMarkerStyle(20);



	c->cd();
	c->SetLogy();
	CorrDiffHis->Draw("ep");

	//pt Binned Correction//


	TH1D * Eff1DHis = (TH1D * ) fin1DEff->Get("Eff1DHis");
	
	TH1D * CorrDiffHisBin = (TH1D * ) CorrDiffHis ->Clone("CorrDiffHisBin");

	CorrDiffHisBin->SetMarkerColor(kRed);
	CorrDiffHisBin->SetLineColor(kRed);	
	CorrDiffHisBin->SetMarkerSize(1);
	CorrDiffHisBin->SetMarkerStyle(20);

	float Eff1D[NBins];
	float Eff1DErr[NBins];

	for(int i = 0; i < NBins;i++){
		RawCount = hPt->GetBinContent(i+1);
		RawCountErr = hPt->GetBinError(i+1);
		Eff1D[i] = Eff1DHis->GetBinContent(i+1);
		Eff1DErr[i] = Eff1DHis->GetBinError(i+1);

		cout << "RawCount = " << RawCount << "  RawCountErr = " << RawCountErr << " Eff1D[i] =   " << Eff1D[i] << "  Eff1DErr[i] =  " << Eff1DErr[i] << endl; 


//		CorrYieldDiff[i] = (RawCount *  Eff1D[i])/(BRchain*2* lumi);
//		CorrYieldDiffErr[i] = TMath::Sqrt((RawCountErr *  Eff1D[i]) *(RawCountErr  *  Eff1D[i]) + (RawCount *  Eff1DErr[i]) * (RawCount  *  Eff1DErr[i]))/(BRchain*2* lumi);
		CorrYieldDiff[i] = (RawCount /  Eff1D[i])/(BRchain*2* lumi);
		CorrYieldDiffErr[i] = TMath::Sqrt((RawCountErr /  Eff1D[i]) *(RawCountErr  /  Eff1D[i]) + (RawCount /Eff1D[i] *  Eff1DErr[i]) * (RawCount /Eff1D[i] *  Eff1DErr[i]))/(BRchain*2* lumi);

		CorrDiffHisBin->SetBinContent(i+1,CorrYieldDiff[i]);
		CorrDiffHisBin->SetBinError(i+1,CorrYieldDiffErr[i]);

	}
	CorrDiffHisBin->Draw("epSAME");

/*
	//2015 B+ pp//

	const int NPt = 5;
	double Pt[NPt+1] = {7,10,15,20,30,50};
	double BPXSec[NPt] = {2610000,744000,197000,46500,5300};
	double StatErr[NPt] = {170000,29000,9000,2400,500};
	double SystErr[NPt] = {230000,59000,15000.3500,400};
	double TotalErr[NPt];


	TH1D * BPXSecHis = new TH1D("BPXSecHis","",NPt,Pt);
	BPXSecHis->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	BPXSecHis->GetYaxis()->SetTitle("d #sigma/d p_{T} (pb GeV^{-1} c)");
	BPXSecHis->GetYaxis()->SetTitleOffset(1.3);
	BPXSecHis->GetXaxis()->CenterTitle();
	BPXSecHis->GetYaxis()->CenterTitle();

	//BPXSecHis->SetTitle("2015 B^{0}_{s} pp Cross Section");
	BPXSecHis->SetMarkerColor(kGreen);
	BPXSecHis->SetLineColor(kGreen);
	
	BPXSecHis->SetMarkerSize(1);
	BPXSecHis->SetMarkerStyle(20);

	for(int i = 0; i < NPt; i++){

		TotalErr[i] = sqrt(StatErr[i] * StatErr[i] + SystErr[i] * SystErr[i]);

		BPXSecHis->SetBinContent(i+1,BPXSec[i]);
		BPXSecHis->SetBinError(i+1,TotalErr[i]);

	}

//	BPXSecHis->Draw("epSAME");

	TLegend *leg = new TLegend(0.30,0.65,0.75,0.85,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.040);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->SetLineWidth(3);

	leg->AddEntry(CorrDiffHis,"2017 pp","pl");
	leg->AddEntry(CorrDiffHisBin,"2017 pp Eff Binned Corr","pl");
	
//	leg->AddEntry(BPXSecHis,"2015 pp","pl");

	//DONE 2015 B+ pp//

	leg->Draw("SAME");

	c->SaveAs("FinalPlots/CorrectedYield9BinsBDT.png");

*/
	foutCorr->cd();
	CorrDiffHis->Write();
	CorrDiffHisBin->Write();
	foutCorr->Close();

	




}

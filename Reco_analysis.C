#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
//#include "edm4eic/ReconstructedParticleData.h"
#include <iostream>

int ican2 = 0;
void makeCanvas()  {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican2++ ), "", 900, 600);
    //can->SetTopMargin(0.04);
    //can->SetRightMargin(0.01);
}

double calc_Phi( TLorentzVector lv1, TLorentzVector lv2) {
    TLorentzVector lvPlus = lv1 + lv2;
    lv1.Boost(-lvPlus.BoostVector());
    lv2.Boost(-lvPlus.BoostVector());
    TLorentzVector lvMinus = lv1 - lv2;
    double Px = lvPlus.Px();
    double Py = lvPlus.Py();
    double Qx = lvMinus.Px();
    double Qy = lvMinus.Py();
    double PcrossQ = (Px*Qy) - (Py*Qx);
    double cosphi = (Px*Qx + Py*Qy) / (lvPlus.Pt()*lvMinus.Pt());
    double PairPhi = acos(cosphi);
    if ( PcrossQ > 0 ){
        return PairPhi - 3.141592;
    } else {
        return 3.141592 - PairPhi;
    }
    return PairPhi;
}


void Reco_analysis() {
    // Creating 1D histograms
    TH1F("h1", "ntuple", 100, -4, 4);

    TFile * fo = new TFile( "analysis_plots.root", "RECREATE" );

    //TH1F("name", "title;xlabel;ylabel;zlabel", n nbins, xmin, xmax)

    auto * mPt = new TH1F("mPt", "J/#Psi Transverse Momentum; Pair P_{T} (GeV/c); counts", 50, 0, 1);
    auto * mMass = new TH1F("mMass", "J/#Psi Mass; M_{ee} (GeV); counts", 100, 0, 5);
    auto * mPairPhi = new TH1F("mPairPhi", "e^{+}e^{-} Phi distribution;Phi (rad);# events", 50, -3.1415, 3.1415);
    auto * mCos2phivsPT = new TH2F("mCos2phivsPT", "cos2#phi distribution vs P_{T}", 100, -2, 2, 10, 0, 0.2);

    //Open file with the tree
    TFile *myFile = TFile::Open("/Users/samcorey/code/data/JpsiEIC_output_podio.root");
    //TFile *myFile = TFile::Open("/Users/samcorey/eic/podio_output.root");
    TTreeReader myReader("events", myFile);

    TTreeReaderArray<float> PxVals(myReader, "ReconstructedParticles.momentum.x");
    TTreeReaderArray<float> PyVals(myReader, "ReconstructedParticles.momentum.y");
    TTreeReaderArray<float> PzVals(myReader, "ReconstructedParticles.momentum.z");
    TTreeReaderArray<float> EVals(myReader, "ReconstructedParticles.energy");

    //TLorentzVector lv1, lv2, lv, lvn;

    int evt_number = 0;

    while (myReader.Next()) {
 
        TLorentzVector lv1, lv2, lv, lvn;

        if (PxVals.GetSize() > 1 ){
            lv1.SetPxPyPzE(PxVals[0], PyVals[0], PzVals[0], EVals[0]);
            lv2.SetPxPyPzE(PxVals[1], PyVals[1], PzVals[1], EVals[1]);

            lv = lv1+lv2;
            double PhiVal = calc_Phi(lv1, lv2);

            mMass->Fill( lv.M() );
            mPt->Fill( lv.Pt() );
            mPairPhi->Fill(PhiVal);

            mCos2phivsPT->Fill( cos(2*PhiVal), lv.Pt() );
        }
        evt_number++;
        
    }

auto * mPTcos2phimoments = mCos2phivsPT->ProfileY("mPTcos2phimoments", 1, -1);

fo -> cd();

makeCanvas();
mMass->SetLineColor(kBlack);
mMass->Draw("E1");
gPad->Print( "RecoPlots/plot_mMass.pdf" );
gPad->Print( "RecoPlots/plot_mMass.png" );

makeCanvas();
mPt->SetLineColor(kBlack);
gPad->SetLogy();
mPt->Draw("E1");
gPad->Print( "RecoPlots/plot_mPt.pdf" );
gPad->Print( "RecoPlots/plot_mPt.png" );

makeCanvas();
mPairPhi->SetLineColor(kBlack);
mPairPhi->Draw("E1");
gPad->Print( "RecoPlots/plot_mPairPhi.pdf" );
gPad->Print( "RecoPlots/plot_mPairPhi.png" );

makeCanvas();
mPTcos2phimoments->SetLineColor(kBlack);
mPTcos2phimoments->SetTitle("Reconstructed 2#phi Fourier coefficient vs. Pair P_{T}; Pair P_{T} (GeV/c); <2cos2#phi>");
mPTcos2phimoments->Draw("E1");
gPad->Print( "RecoPlots/plot_mPTcos2phimoments.pdf" );
gPad->Print( "RecoPlots/plot_mPTcos2phimoments.png" );

fo->Write();
}

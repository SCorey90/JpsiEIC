#include "HepMC3/GenEvent.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/ReaderRoot.h"
#include "HepMC3/Print.h"

#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TRandom3.h"

#include <random>
#include <vector>
#include <iostream>

using namespace HepMC3;

std::random_device global_rng;
TRandom3 rng(global_rng());
TRandom ptr_rng(global_rng());

int ican2 = 0;
void makeCanvas()  {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican2++ ), "", 900, 600);
    //can->SetTopMargin(0.04);
    can->SetRightMargin(0.35);
}

TLorentzVector gen_Jpsi( TH1F* Pt_hist ) {
    double jpsi_mass = rng.Gaus(3.1, 0.01);
    double jpsi_pt = Pt_hist->GetRandom();
    double jpsi_phi = rng.Uniform(-3.141592, 3.141592);
    double jpsi_eta = rng.Uniform(-6, 6);
    TLorentzVector lvJpsi;
    lvJpsi.SetPtEtaPhiM( jpsi_pt, jpsi_eta, jpsi_phi, jpsi_mass );

    return lvJpsi;
}

vector<TLorentzVector> decay_Jpsi( TH1F *cos2phi_hist , TLorentzVector lvJpsi ) {
    double elec_mass = 0.000511 ;
    double Jpsi_pt = lvJpsi.Pt();

    double cos2phi_val = cos2phi_hist->GetBinContent( cos2phi_hist->FindBin(Jpsi_pt) ) ;
    double phi_val = acos( cos2phi_val ) / 2 ;

    TLorentzVector PosLV, ElecLV ;
    double az_ang =  lvJpsi.Phi() + phi_val ;
    double pol_ang = acos( rng.Uniform(-1, 1) );
    double absP = sqrt((lvJpsi.M()*lvJpsi.M()/4) - elec_mass*elec_mass);
    PosLV.SetPxPyPzE( (absP*sin(pol_ang)*cos(az_ang)), (absP*sin(pol_ang)*sin(az_ang)), (absP*cos(pol_ang)), (sqrt(absP*absP + elec_mass*elec_mass)) );
    ElecLV.SetPxPyPzE( (-absP*sin(pol_ang)*cos(az_ang)), (-absP*sin(pol_ang)*sin(az_ang)), (-absP*cos(pol_ang)), (sqrt(absP*absP + elec_mass*elec_mass)) );

    PosLV.Boost(lvJpsi.BoostVector());
    ElecLV.Boost(lvJpsi.BoostVector());
    vector<TLorentzVector> pair_LVs{ PosLV, ElecLV };

    return pair_LVs;
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
}



/** Main */
int main(int argc, char **argv) {
    //arg setup
    if( argc<2 ) {
        std::cout << "Usage: " << argv[0] << "<output_hepmc3_file>" << std::endl;
        exit(-1);
    }

    WriterAscii text_output(argv[1]);

    //sample setup
    double p_photon = 1 ;

    auto *pt_csv = new TTree("pt_csv", "tree from csv");
    pt_csv->ReadFile("/Users/samcorey/code/JpsiEIC/qt2crosssection.csv", "x/D:y");
    pt_csv->Draw("y:x");

    auto *cos2phi_csv = new TTree("cos2phi_csv", "tree from csv");
    cos2phi_csv->ReadFile("/Users/samcorey/code/JpsiEIC/cos2phivspt.csv", "x/D:y");
    cos2phi_csv->Draw("y:x");

    auto *no_photon_cos2phi_csv = new TTree("no_photon_cos2phi_csv", "tree from csv");
    no_photon_cos2phi_csv->ReadFile("/Users/samcorey/code/JpsiEIC/no_photon_cos2phivspt.csv", "x/D:y");
    no_photon_cos2phi_csv->Draw("y:x");

    auto * pt_hist = new TH1F("pt_hist", "P_T hist to sample; P_{T}; counts", 50, 0, sqrt(0.16));
    double * sample_pt = pt_csv->GetV2();
    double * sample_dndp = pt_csv->GetV1();
    for (int n=0; n<50; n++) {
        for (int i=0; i<83; i++) {
            if ( abs(sqrt(sample_pt[i]) - pt_hist->GetBinCenter(n+1)) < 0.8 * pt_hist->GetBinWidth(n+1) ) {
                pt_hist->SetBinContent(n+1, sample_dndp[i]);
            }
        }
    }

    auto * mCos2phi_hist = new TH1F("cos2phi_hist", "cos2phi hist to sample", 100, 0, 0.2);
    double * sample_xvals = cos2phi_csv->GetV2();
    double * sample_cos2phi = cos2phi_csv->GetV1();
    for (int n=1; n<101; n++) {
        vector<double> pt_diffs;
        for (int i=0; i<166 ; i++) {
            pt_diffs.push_back( abs(mCos2phi_hist->GetBinCenter(n+1) - sample_xvals[i]) );
        }
        int min_index = distance( begin(pt_diffs), min_element( begin(pt_diffs), end(pt_diffs) ) );
        mCos2phi_hist->SetBinContent(n, sample_cos2phi[min_index]);
    }

    auto * mNo_Photon_Cos2phi_hist = new TH1F("no_photon_cos2phi_hist", "cos2phi hist to sample", 100, 0, 0.2);
    double * no_photon_xvals = no_photon_cos2phi_csv->GetV2();
    double * no_photon_cos2phi = no_photon_cos2phi_csv->GetV1();
    for (int n=1; n<101; n++) {
        vector<double> pt_diffs;
        for (int i=0; i<253 ; i++) {
            pt_diffs.push_back( abs(mNo_Photon_Cos2phi_hist->GetBinCenter(n+1) - no_photon_xvals[i]) );
        }
        int min_index_no_photon = distance( begin(pt_diffs), min_element( begin(pt_diffs), end(pt_diffs) ) );
        mNo_Photon_Cos2phi_hist->SetBinContent(n, no_photon_cos2phi[min_index_no_photon]);
    }

    //tests
    auto * mTest_Cos2phivsPT = new TH2F("mTest_Cos2phivsPT", "Sampled electron pair", 100, -2, 2, 50, 0, 0.2);

    //event loop
    int events_parsed = 0;

    for( int i = 0; i<3000; i++ ) {

        GenEvent evt;

        TLorentzVector JpsiVec = gen_Jpsi(pt_hist);

        double pcheck = rng.Uniform(0, 1);
        if ( pcheck < p_photon ) {
            vector<TLorentzVector> pair_LVs = decay_Jpsi( mCos2phi_hist , JpsiVec );
            TLorentzVector PosLV = pair_LVs[0];
            TLorentzVector ElecLV = pair_LVs[1];

            mTest_Cos2phivsPT->Fill( 2*cos(2*calc_Phi(PosLV, ElecLV)), (PosLV+ElecLV).Pt() );

            double px0 = JpsiVec.Px();
            double py0 = JpsiVec.Py();
            double pz0 = JpsiVec.Pz();
            double E0 = JpsiVec.E();

            double px1 = PosLV.Px();
            double py1 = PosLV.Py();
            double pz1 = PosLV.Pz();
            double E1 = PosLV.E();

            double px2 = ElecLV.Pz();
            double py2 = ElecLV.Py();
            double pz2 = ElecLV.Pz();
            double E2 = ElecLV.E();

            GenParticlePtr p0 = std::make_shared<GenParticle>( FourVector( px0, py0, pz0, E0  ), 443, 4 );
            GenParticlePtr p1 = std::make_shared<GenParticle>( FourVector( px1, py1, pz1, E1  ), -11,  1 );
            GenParticlePtr p2 = std::make_shared<GenParticle>( FourVector( px2, py2, pz2, E2  ), 11,  1 );

            GenVertexPtr v1 = std::make_shared<GenVertex>();
            v1->add_particle_out(p1);
            v1->add_particle_out(p2);
            v1->add_particle_in(p0);
            v1->set_status(4);
            evt.add_vertex(v1);

        } else {
            vector<TLorentzVector> pair_LVs = decay_Jpsi( mNo_Photon_Cos2phi_hist , JpsiVec );
            TLorentzVector PosLV = pair_LVs[0];
            TLorentzVector ElecLV = pair_LVs[1];

            mTest_Cos2phivsPT->Fill( 2*cos(2*calc_Phi(PosLV, ElecLV)), (PosLV+ElecLV).Pt() );

            double px0 = JpsiVec.Px();
            double py0 = JpsiVec.Py();
            double pz0 = JpsiVec.Pz();
            double E0 = JpsiVec.E();

            double px1 = PosLV.Px();
            double py1 = PosLV.Py();
            double pz1 = PosLV.Pz();
            double E1 = PosLV.E();
        
            double px2 = ElecLV.Pz();
            double py2 = ElecLV.Py();
            double pz2 = ElecLV.Pz();
            double E2 = ElecLV.E();

            GenParticlePtr p0 = std::make_shared<GenParticle>( FourVector( px0, py0, pz0, E0  ), 443, 4 );
            GenParticlePtr p1 = std::make_shared<GenParticle>( FourVector( px1, py1, pz1, E1  ), -11,  1 );
            GenParticlePtr p2 = std::make_shared<GenParticle>( FourVector( px2, py2, pz2, E2  ), 11,  1 );

            GenVertexPtr v1 = std::make_shared<GenVertex>();
	    v1->add_particle_out(p1);
            v1->add_particle_out(p2);
            v1->add_particle_in(p0);
            v1->set_status(4);
            evt.add_vertex(v1);
        }

        if( events_parsed == 0 ) {
            std::cout << "First event: " << std::endl;
            Print::listing(evt);
        }

        text_output.write_event(evt);
        ++events_parsed;

        if( events_parsed%100 == 0 ) {
            std::cout << "Event: " << events_parsed << std::endl;
        }
    }

    text_output.close();

    std::cout << "Events parsed and written: " << events_parsed << std::endl;

    makeCanvas();
    mTest_Cos2phivsPT->Draw("colz");
    gPad->Print( "/Users/samcorey/code/JpsiEIC/plots/plot_mTest_Cos2phivsPT.pdf" );
    gPad->Print( "/Users/samcorey/code/JpsiEIC/plots/plot_mTest_Cos2phivsPT.png" );

    return 0;
}

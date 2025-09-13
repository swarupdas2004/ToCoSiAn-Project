#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"

using namespace std;

void analysis2(){
    TFile* outFile=new TFile("/Users/swarupdas/icloud/Documents/ToCoSiAn/data/analysis_output1.root","RECREATE");
    TTree* tree=new TTree("dataTree","Pion+ + Pion+ --> Pion- (-1.5,1.5)GeV");
    tree->ReadFile("/Users/swarupdas/icloud/Documents/ToCoSiAn/data/event2.dat", "pid/I:px/D:py/D:pz/D");

    Int_t pid;
    Double_t px,py,pz;

    tree->SetBranchAddress("pid", &pid);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);

    TH1D* h_pT=new TH1D("h_pT","Transverse Momenta;p_{T} Mev;Entries",100,0,2500);
    TH1D* h_eta=new TH1D("h_eta","Pseudorapidity;#eta;Entries",100,-5,5);
    TH1D* h_phi=new TH1D("h_phi","Azimuthal Angle;#phi [rad];Entries",100,-3.2,3.2);
    TH1D* h_pid=new TH1D("h_pid","Particle IDs; PDG ID Code;Entries",100,-212,212);

    Long64_t nEntries=tree->GetEntries();

    cout<<"Processing "<<nEntries<<" particle collisons...\n";

    for(Long64_t i=1; i<=nEntries; i++){

        tree->GetEntry(i);

        double p= sqrt(px*px + py*py + pz*pz);
        double pT= sqrt(px*px + py*py);
        double phi=atan2(py,px);
        double eta=0.0;

        if (p>1e-10 && (p-abs(pT))>1e-10){
            eta=-log(tan(acos(pz/p)/2.0));
        }

        h_pT->Fill(pT);
        h_eta->Fill(eta);
        h_phi->Fill(phi);
        h_pid->Fill(pid);
    }

    gStyle->SetOptStat("emr");
    gStyle->SetOptFit(1111);

    TF1* fitFunc_pT=new TF1("fitFunc_pT","gaus",500,2000);
    fitFunc_pT->SetParameters(h_pT->GetMaximum(),h_pT->GetMean(),h_pT->GetStdDev());
    fitFunc_pT->SetLineColor(kMagenta);
    

    TF1* fitFunc_phi=new TF1("fitFunc_eta","pol0",-3.2,3.2);
    fitFunc_phi->SetLineColor(kGreen+7);
    

    TCanvas* c1=new TCanvas("c1","Analysis plots",1200,800);
    c1->Divide(2,2);

    c1->cd(1);
    h_pT->SetLineColor(kBlue);h_pT->SetFillColorAlpha(kBlue, 0.35);h_pT->SetFillStyle(3001);
    h_pT->Draw("HIST");
    h_pT->Fit(fitFunc_pT,"R+");
    fitFunc_pT->Draw("SAME");

    c1->cd(2);
    h_eta->SetLineColor(kRed);h_eta->SetFillColorAlpha(kRed, 0.35);h_eta->SetFillStyle(3001);
    h_eta->Draw("HIST");

    c1->cd(3);
    h_phi->SetLineColor(kGreen+2);h_phi->SetFillColorAlpha(kGreen+2, 0.35);h_phi->SetFillStyle(3001);
    h_phi->Draw("HIST");
    h_phi->Fit(fitFunc_phi,"R+");
    fitFunc_phi->Draw("SAME");

    c1->cd(4);
    h_pid->SetLineColor(kViolet);h_pid->SetFillColorAlpha(kViolet, 0.35);h_pid->SetFillStyle(3001);
    h_pid->Draw("HIST");

    c1->SaveAs("/Users/swarupdas/icloud/Documents/ToCoSiAn/data/Analysis_plots1.pdf");

    outFile->cd();
    h_pT->Write();
    h_eta->Write();
    h_phi->Write();
    h_pid->Write();

    outFile->Close();
    delete outFile;


    cout<<"=========Fit Results for pT Distribution========\n";
    cout<<"Fit Function : Gaussian\n";
    cout<<"Chi^2/ndf : "<<fitFunc_pT->GetChisquare()<<"/"<<fitFunc_pT->GetNDF()<<endl;
    cout<<"Mean : "<<fitFunc_pT->GetParameter(1)<<"+/-"<<fitFunc_pT->GetParError(1)<<" MeV\n";
    cout<<"Sigma : "<<fitFunc_pT->GetParameter(2)<<"+/-"<<fitFunc_pT->GetParError(2)<<" MeV\n";

    cout<<"=========Fit Results for phi Distribution=========\n";
    cout<<"Fit Function : Constant\n";
    cout<<"Constant Value : "<<fitFunc_phi->GetParameter(0)<<"+/-"<<fitFunc_phi->GetParError(0)<<"\n";

}
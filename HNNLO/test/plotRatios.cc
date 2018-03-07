#include "HelperFunctions.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TText.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TArrayI.h"


using namespace std;
using namespace HelperFunctions;


const unsigned int nfolders=3;
TString folderName[nfolders]={
  "Nominal",
  "mteqmb",
  "mtlarge"
};
TString orderName[3]={ "LO","NLO","NNLO" };

void extractRatios(int iorder){
  assert(iorder>=0 && iorder<=2);

  TFile* foutput = TFile::Open(orderName[iorder] + "_Ratios.root", "recreate");
  vector<TH1F*> hList[3];
  for (unsigned int f=0; f<nfolders; f++){
    for (unsigned int c=0; c<(iorder==2 ? 2 : 1); c++){
      TString cinput="./";
      cinput = cinput + folderName[f] + "/" + orderName[iorder] + "_NNPDFNominal_13TeV_125";
      if (iorder==2) cinput = cinput + Form("_Copy%i", c);
      cinput += ".root";

      vector<TH1F*> hList_tmp;
      TFile* finput = TFile::Open(cinput, "read");
      finput->cd();
      extractHistogramsFromDirectory(finput, hList_tmp);
      //for (TH1F*& hh:hList_tmp) hh->Rebin(2);

      foutput->cd();
      if (hList[f].empty()){
        for (TH1F*& hh:hList_tmp) hList[f].push_back(new TH1F(*hh));
      }
      else{
        for (unsigned int ih=0; ih<hList_tmp.size(); ih++){
          for (int ix=0; ix<=hList[f].at(ih)->GetNbinsX()+1; ix++){
            double c1=hList[f].at(ih)->GetBinContent(ix);
            double c2=hList_tmp.at(ih)->GetBinContent(ix);
            double e1=hList[f].at(ih)->GetBinError(ix);
            double e2=hList_tmp.at(ih)->GetBinError(ix);
            double c=0, e=0;
            if (c1>0. && c2>0.){
              c=(c1+c2)/2.;
              e=sqrt(pow(e1, 2)+pow(e2, 2))/2.;
            }
            else if (c2>0.){
              c=c2;
              e=e2;
            }
            else if (c1>0.){
              c=c1;
              e=e1;
            }
            hList[f].at(ih)->SetBinContent(ix, c);
            hList[f].at(ih)->SetBinError(ix, e);
          }
        }
      }
      finput->Close();
    }
  }
  for (unsigned int f=0; f<nfolders; f++){
    for (TH1F*& hh:hList[f]) hh->Scale(1./getHistogramIntegralAndError(hh, 0, hh->GetNbinsX()+1, false, nullptr));
  }
  for (unsigned int f=1; f<nfolders; f++){
    TDirectory* hfolder = foutput->mkdir(folderName[f]);
    hfolder->cd();
    for (unsigned int ih=0; ih<hList[f].size(); ih++){
      TH1F* hRatio = new TH1F(*(hList[f].at(ih)));
      TString hname = hRatio->GetName(); replaceString(hname, "distribution", "ratio");
      hRatio->SetName(hname);
      hRatio->SetTitle("");
      hRatio->Reset("ICESM");
      divideHistograms(hList[f].at(ih), hList[0].at(ih), hRatio, false);

      TString xtitle = hRatio->GetXaxis()->GetTitle();
      TString ytitle = hRatio->GetYaxis()->GetTitle();
      replaceString(xtitle, "y", "y(");
      replaceString(xtitle, "pt", "p_{T}(");
      replaceString(xtitle, "eta", "#eta(");
      xtitle += ")";
      ytitle = "Ratio";
      hRatio->GetXaxis()->SetTitle(xtitle);
      hRatio->GetYaxis()->SetTitle(ytitle);
      hRatio->GetYaxis()->SetTitleOffset(1.2);

      hfolder->WriteTObject(hRatio);
      delete hRatio;
    }
  }
  for (unsigned int f=0; f<nfolders; f++){ for (TH1F*& hh:hList[f]) delete hh; }
  foutput->Close();
}

void plotRatios(TString distribution){
  TFile* finput[3];
  TH1F* hRatio[3][2];

  for (unsigned int iorder=0; iorder<3; iorder++){
    finput[iorder]=TFile::Open(orderName[iorder]+"_Ratios.root", "read");
    for (unsigned int f=0; f<2; f++){
      hRatio[iorder][f] = (TH1F*) finput[iorder]->Get(folderName[f+1]+"/"+distribution+"_ratio");
    }
  }

  gStyle->SetOptStat(0);
  {
    int colors[100];
    Double_t Red[]    ={ 0.1, 0.1, 1.0 };
    Double_t Green[]  ={ 0.1, 1.0, 0.1 };
    Double_t Blue[]   ={ 1.0, 0.1, 0.1 };
    Double_t Length[] ={ 0.00, 0.50, 1.00 };
    int FI = TColor::CreateGradientColorTable(3, Length, Red, Green, Blue, 100);
    const unsigned int ncolors = gStyle->GetNumberOfColors();
    if (FI<0) cout << "Failed to set color palette" << endl;
    else{
      for (unsigned int ic=0; ic<100; ic++) colors[ic] = FI+ic;
      gStyle->SetPalette(100, colors);
    }
    cout << "Ncolors: " << ncolors << endl;
  }

  vector<TH1F*> histolist;
  vector<TString> strLabels;
  for (unsigned int f=0; f<2; f++){
    for (unsigned int iorder=0; iorder<3; iorder++){
      if (iorder==0 && distribution=="pt3456") continue;
      histolist.push_back(hRatio[iorder][f]);
      TString strLabel;
      if (f==0) strLabel = "m_{t}=m_{b}=4.75 GeV";
      else strLabel = "Large m_{t}, m_{b}=4.75 GeV";
      strLabel = strLabel + " (" + orderName[iorder] + ")";
      strLabels.push_back(strLabel);
    }
  }
  const unsigned int nhistos=histolist.size();
  {
    const unsigned int ncolors = gStyle->GetNumberOfColors();
    const unsigned int stepsize = ncolors/nhistos;
    for (unsigned int ih=0; ih<histolist.size(); ih++){
      int colorToUse = gStyle->GetColorPalette(ih*stepsize);
      histolist.at(ih)->SetMarkerColor(colorToUse);
      histolist.at(ih)->SetLineColor(colorToUse);
      histolist.at(ih)->SetLineWidth(2);
    }
  }
  TString canvasname = TString("c_")+distribution;
  TCanvas* canvas = new TCanvas(canvasname, "", 8, 30, 800, 800);
  canvas->cd();
  gStyle->SetOptStat(0);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetLeftMargin(0.17);
  canvas->SetRightMargin(0.05);
  canvas->SetTopMargin(0.07);
  canvas->SetBottomMargin(0.13);
  canvas->SetFrameFillStyle(0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameFillStyle(0);
  canvas->SetFrameBorderMode(0);

  TLegend* legend = new TLegend(0.38, 0.90-0.12/3.*nhistos, 0.70, 0.90);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->SetLineColor(1);
  legend->SetLineStyle(1);
  legend->SetLineWidth(1);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);

  TPaveText* pt = new TPaveText(0.15, 0.93, 0.85, 1, "brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  pt->SetTextFont(42);
  pt->SetTextSize(0.045);
  TText* text = pt->AddText(0.025, 0.45, "#font[61]{CMS}");
  text->SetTextSize(0.044);
  text = pt->AddText(0.165, 0.42, "#font[52]{Simulation}");
  text->SetTextSize(0.0315);
  TString cErgTev = Form("#font[42]{%i TeV}", 13);
  text = pt->AddText(0.99, 0.45, cErgTev);
  text->SetTextSize(0.0315);

  bool first=true;
  float minY=9e9, maxY=-9e9;
  for (unsigned int ih=0; ih<nhistos; ih++){
    TH1F* hh=histolist.at(ih);

    int binXlow = 1;
    int binXhigh = hh->GetNbinsX();
    if (distribution=="y3456") binXhigh=hh->GetXaxis()->FindBin(3.99);
    float minX, maxX;
    minX=hh->GetXaxis()->GetBinLowEdge(binXlow);
    maxX=hh->GetXaxis()->GetBinUpEdge(binXhigh);

    hh->GetXaxis()->SetRangeUser(minX, maxX);
    hh->GetXaxis()->SetNdivisions(505);
    hh->GetXaxis()->SetLabelFont(42);
    hh->GetXaxis()->SetLabelOffset(0.007);
    hh->GetXaxis()->SetLabelSize(0.04);
    hh->GetXaxis()->SetTitleSize(0.06);
    hh->GetXaxis()->SetTitleOffset(0.9);
    hh->GetXaxis()->SetTitleFont(42);
    hh->GetYaxis()->SetNdivisions(505);
    hh->GetYaxis()->SetLabelFont(42);
    hh->GetYaxis()->SetLabelOffset(0.007);
    hh->GetYaxis()->SetLabelSize(0.04);
    hh->GetYaxis()->SetTitleSize(0.06);
    hh->GetYaxis()->SetTitleOffset(1.1);
    hh->GetYaxis()->SetTitleFont(42);
    for (int ix=binXlow; ix<binXhigh; ix++){
      float bc = hh->GetBinContent(ix);
      if (bc!=0.){
        minY = std::min(bc, minY);
        maxY = std::max(bc, maxY);
      }
    }
    cout << "Min = " << minY << ", max = " << maxY << " after " << hh->GetName() << " (" << hh->GetTitle() << ")" << endl;
    legend->AddEntry(hh, strLabels.at(ih), "l");
  }
  for (TH1F*& hh:histolist){
    hh->GetYaxis()->SetRangeUser((minY<0. ? minY*1.2 : minY*0.8), maxY*1.2);
    if (first){
      hh->Draw("hist");
      first=false;
    }
    else{
      hh->Draw("histsame");
    }
  }
  legend->Draw("same");
  pt->Draw();
  canvas->RedrawAxis();
  canvas->Modified();
  canvas->Update();
  canvas->SaveAs(Form("%s%s", canvasname.Data(), ".png"));
  canvas->SaveAs(Form("%s%s", canvasname.Data(), ".pdf"));
  delete pt;
  delete legend;
  canvas->Close();

  for (unsigned int iorder=0; iorder<3; iorder++) finput[iorder]->Close();
}

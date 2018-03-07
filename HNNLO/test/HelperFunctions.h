#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include <ctime>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include "TClass.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TObject.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"

namespace HelperFunctions{

  float calculateEfficiencyError(
    float const sumW, float const sumWAll,
    float const sumWsq, float const sumWsqAll
  );

  template<typename T> bool checkVarNanInf(T const& var);

  template<typename T> void extractHistogramsFromDirectory(TDirectory* source, std::vector<T*>& histolist);

  template <typename T> void divideHistograms(T* hnum, T* hden, T*& hAssign, bool useEffErr);
  template<> void divideHistograms<TH1F>(TH1F* hnum, TH1F* hden, TH1F*& hAssign, bool useEffErr);
  template<> void divideHistograms<TH2F>(TH2F* hnum, TH2F* hden, TH2F*& hAssign, bool useEffErr);
  template<> void divideHistograms<TH3F>(TH3F* hnum, TH3F* hden, TH3F*& hAssign, bool useEffErr);

  template<typename T, typename U> void replaceString(T& strinput, U strTakeOut, U strPutIn);
  template<> void replaceString<TString, const TString>(TString& strinput, const TString strTakeOut, const TString strPutIn);
  template<> void replaceString<TString, const char*>(TString& strinput, const char* strTakeOut, const char* strPutIn);
  template<> void replaceString<std::string, const std::string>(std::string& strinput, const std::string strTakeOut, const std::string strPutIn);
  template<> void replaceString<std::string, const char*>(std::string& strinput, const char* strTakeOut, const char* strPutIn);

  template <typename T> double getHistogramIntegralAndError(T const* histo, int ix, int jx, bool useWidth, double* error=nullptr);
  template <typename T> double getHistogramIntegralAndError(T const* histo, int ix, int jx, int iy, int jy, bool useWidth, double* error=nullptr);
  template <typename T> double getHistogramIntegralAndError(T const* histo, int ix, int jx, int iy, int jy, int iz, int jz, bool useWidth, double* error=nullptr);

}

template<typename T> bool HelperFunctions::checkVarNanInf(T const& var){
  return !(std::isnan(var) || std::isinf(var));
}
template bool HelperFunctions::checkVarNanInf<short>(short const& var);
template bool HelperFunctions::checkVarNanInf<unsigned int>(unsigned int const& var);
template bool HelperFunctions::checkVarNanInf<int>(int const& var);
template bool HelperFunctions::checkVarNanInf<float>(float const& var);
template bool HelperFunctions::checkVarNanInf<double>(double const& var);

float HelperFunctions::calculateEfficiencyError(
  float const sumW, float const sumWAll,
  float const sumWsq, float const sumWsqAll
){
  float const& sumWp=sumW;
  float const& sumWsqp=sumWsq;
  float const sumWm = sumWAll-sumWp;
  float const sumWsqm = sumWsqAll-sumWsqp;
  float numerator, denominator;
  float ratio=0;
  if (sumWAll!=0.){
    numerator = sqrt(sumWsqp*pow(sumWm, 2) + sumWsqm*pow(sumWp, 2));
    denominator = pow(sumWAll, 2);
    ratio = numerator/denominator;
  }
  return ratio;
}

template<> void HelperFunctions::divideHistograms<TH1F>(TH1F* hnum, TH1F* hden, TH1F*& hAssign, bool useEffErr){
  if (hnum->GetNbinsX()!=hden->GetNbinsX()) return; const int nbinsx = hnum->GetNbinsX();
  for (int binx=0; binx<=nbinsx+1; binx++){
    float sumW = hnum->GetBinContent(binx);
    float sumWAll = hden->GetBinContent(binx);
    float sumWsq = pow(hnum->GetBinError(binx), 2);
    float sumWsqAll = pow(hden->GetBinError(binx), 2);
    float bincontent=0;
    float binerror=0;
    if (sumWAll!=0.) bincontent = sumW/sumWAll;
    if (useEffErr) binerror = calculateEfficiencyError(sumW, sumWAll, sumWsq, sumWsqAll);
    else binerror = bincontent*sqrt(sumWsq/pow(sumW, 2) + sumWsqAll/pow(sumWAll, 2));
    if (!checkVarNanInf(bincontent) || !checkVarNanInf(binerror)){
      bincontent=0;
      binerror=0;
    }
    hAssign->SetBinContent(binx, bincontent);
    hAssign->SetBinError(binx, binerror);
  }
}
template<> void HelperFunctions::divideHistograms<TH2F>(TH2F* hnum, TH2F* hden, TH2F*& hAssign, bool useEffErr){
  if (hnum->GetNbinsX()!=hden->GetNbinsX()) return; const int nbinsx = hnum->GetNbinsX();
  if (hnum->GetNbinsY()!=hden->GetNbinsY()) return; const int nbinsy = hnum->GetNbinsY();
  for (int binx=0; binx<=nbinsx+1; binx++){
    for (int biny=0; biny<=nbinsy+1; biny++){
      float sumW = hnum->GetBinContent(binx, biny);
      float sumWAll = hden->GetBinContent(binx, biny);
      float sumWsq = pow(hnum->GetBinError(binx, biny), 2);
      float sumWsqAll = pow(hden->GetBinError(binx, biny), 2);
      float bincontent=0;
      float binerror=0;
      if (sumWAll!=0.) bincontent = sumW/sumWAll;
      if (useEffErr) binerror = calculateEfficiencyError(sumW, sumWAll, sumWsq, sumWsqAll);
      else binerror = bincontent*sqrt(sumWsq/pow(sumW, 2) + sumWsqAll/pow(sumWAll, 2));
      if (!checkVarNanInf(bincontent) || !checkVarNanInf(binerror)){
        bincontent=0;
        binerror=0;
      }
      hAssign->SetBinContent(binx, biny, bincontent);
      hAssign->SetBinError(binx, biny, binerror);
    }
  }
}
template<> void HelperFunctions::divideHistograms<TH3F>(TH3F* hnum, TH3F* hden, TH3F*& hAssign, bool useEffErr){
  if (hnum->GetNbinsX()!=hden->GetNbinsX()) return; const int nbinsx = hnum->GetNbinsX();
  if (hnum->GetNbinsY()!=hden->GetNbinsY()) return; const int nbinsy = hnum->GetNbinsY();
  if (hnum->GetNbinsZ()!=hden->GetNbinsZ()) return; const int nbinsz = hnum->GetNbinsZ();
  for (int binx=0; binx<=nbinsx+1; binx++){
    for (int biny=0; biny<=nbinsy+1; biny++){
      for (int binz=0; binz<=nbinsz+1; binz++){
        float sumW = hnum->GetBinContent(binx, biny, binz);
        float sumWAll = hden->GetBinContent(binx, biny, binz);
        float sumWsq = pow(hnum->GetBinError(binx, biny, binz), 2);
        float sumWsqAll = pow(hden->GetBinError(binx, biny, binz), 2);
        float bincontent=0;
        float binerror=0;
        if (sumWAll!=0.) bincontent = sumW/sumWAll;
        if (useEffErr) binerror = calculateEfficiencyError(sumW, sumWAll, sumWsq, sumWsqAll);
        else binerror = bincontent*sqrt(sumWsq/pow(sumW, 2) + sumWsqAll/pow(sumWAll, 2));
        if (!checkVarNanInf(bincontent) || !checkVarNanInf(binerror)){
          bincontent=0;
          binerror=0;
        }
        hAssign->SetBinContent(binx, biny, binz, bincontent);
        hAssign->SetBinError(binx, biny, binz, binerror);
      }
    }
  }
}

template<typename T> void HelperFunctions::extractHistogramsFromDirectory(TDirectory* source, std::vector<T*>& histolist){
  // Copy all objects and subdirs of directory source as a subdir of the current directory
  TDirectory* target = gDirectory;
  source->ls();
  source->cd();
  // Loop on all entries of this directory
  TKey* key;
  TIter nextkey(source->GetListOfKeys());
  vector<TString> copiedKeys;
  while ((key = (TKey*) nextkey())){
    const char* classname = key->GetClassName();
    TClass* cl = gROOT->GetClass(classname);
    if (!cl) continue;
    if (cl->InheritsFrom(TDirectory::Class())){
      source->cd(key->GetName());
      TDirectory* subdir = gDirectory;
      source->cd();
      extractHistogramsFromDirectory<T>(subdir, histolist);
    }
    else if (cl->InheritsFrom(T::Class())){
      T* hist = (T*) source->Get(key->GetName());
      TString histname=hist->GetName();
      if ((histname=="Graph" && TString(key->GetName())!="Graph") || histname=="") histname=key->GetName(); // Holy jumping monkeys for fake rates
      hist->SetName(histname);
      bool alreadyCopied=false;
      for (auto& k:copiedKeys){
        if (k==key->GetName()){
          alreadyCopied=true;
          break;
        }
      }
      if (!alreadyCopied){
        if (hist){
          copiedKeys.push_back(key->GetName());
          histolist.push_back(hist);
        }
      }
    }
  }
  target->cd();
}
template void HelperFunctions::extractHistogramsFromDirectory<TH1F>(TDirectory* source, std::vector<TH1F*>& histolist);
template void HelperFunctions::extractHistogramsFromDirectory<TH2F>(TDirectory* source, std::vector<TH2F*>& histolist);
template void HelperFunctions::extractHistogramsFromDirectory<TH3F>(TDirectory* source, std::vector<TH3F*>& histolist);
template void HelperFunctions::extractHistogramsFromDirectory<TH1D>(TDirectory* source, std::vector<TH1D*>& histolist);
template void HelperFunctions::extractHistogramsFromDirectory<TH2D>(TDirectory* source, std::vector<TH2D*>& histolist);
template void HelperFunctions::extractHistogramsFromDirectory<TH3D>(TDirectory* source, std::vector<TH3D*>& histolist);
// Overloads for TGraph that can be used just like histograms
template void HelperFunctions::extractHistogramsFromDirectory<TGraphAsymmErrors>(TDirectory* source, std::vector<TGraphAsymmErrors*>& histolist);
template void HelperFunctions::extractHistogramsFromDirectory<TGraphErrors>(TDirectory* source, std::vector<TGraphErrors*>& histolist);
template void HelperFunctions::extractHistogramsFromDirectory<TGraph>(TDirectory* source, std::vector<TGraph*>& histolist);

template<> void HelperFunctions::replaceString<TString, const TString>(TString& strinput, const TString strTakeOut, const TString strPutIn){
  Ssiz_t ipos=strinput.Index(strTakeOut);
  if (ipos!=-1) strinput.Replace(ipos, strTakeOut.Length(), strPutIn);
}
template<> void HelperFunctions::replaceString<TString, const char*>(TString& strinput, const char* strTakeOut, const char* strPutIn){
  Ssiz_t ipos=strinput.Index(strTakeOut);
  if (ipos!=-1) strinput.Replace(ipos, strlen(strTakeOut), strPutIn);
}
template<> void HelperFunctions::replaceString<std::string, const std::string>(std::string& strinput, const std::string strTakeOut, const std::string strPutIn){
  std::string::size_type ipos=strinput.find(strTakeOut);
  if (ipos!=std::string::npos) strinput.replace(ipos, strTakeOut.length(), strPutIn);
}
template<> void HelperFunctions::replaceString<std::string, const char*>(std::string& strinput, const char* strTakeOut, const char* strPutIn){
  std::string::size_type ipos=strinput.find(strTakeOut);
  if (ipos!=std::string::npos) strinput.replace(ipos, strlen(strTakeOut), strPutIn);
}

template <typename T> double HelperFunctions::getHistogramIntegralAndError(T const* histo, int ix, int jx, bool useWidth, double* error){
  double res=0;
  double reserror=0;
  if (histo){
    if (!useWidth) res=histo->IntegralAndError(ix, jx, reserror, "");
    else{
      int xb[2]={ std::max(1, std::min(histo->GetNbinsX(), ix)), std::max(1, std::min(histo->GetNbinsX(), jx)) };

      res=histo->IntegralAndError(xb[0], xb[1], reserror, "width");

      double integralinside, integralerrorinside;
      integralinside=histo->IntegralAndError(xb[0], xb[1], integralerrorinside, "");

      double integraloutside, integralerroroutside;
      integraloutside=histo->IntegralAndError(ix, jx, integralerroroutside, "");

      res = res + integraloutside - integralinside;
      reserror = sqrt(std::max(0., pow(reserror, 2) + pow(integralerroroutside, 2) - pow(integralinside, 2)));
    }
  }
  if (error) *error=reserror;
  return res;
}
template <typename T> double HelperFunctions::getHistogramIntegralAndError(T const* histo, int ix, int jx, int iy, int jy, bool useWidth, double* error){
  double res=0;
  double reserror=0;
  if (histo){
    if (!useWidth) res=histo->IntegralAndError(ix, jx, iy, jy, reserror, "");
    else{
      int xb[2]={ std::max(1, std::min(histo->GetNbinsX(), ix)), std::max(1, std::min(histo->GetNbinsX(), jx)) };
      int yb[2]={ std::max(1, std::min(histo->GetNbinsY(), iy)), std::max(1, std::min(histo->GetNbinsY(), jy)) };

      res=histo->IntegralAndError(xb[0], xb[1], yb[0], yb[1], reserror, "width");

      double integralinside, integralerrorinside;
      integralinside=histo->IntegralAndError(xb[0], xb[1], yb[0], yb[1], integralerrorinside, "");

      double integraloutside, integralerroroutside;
      integraloutside=histo->IntegralAndError(ix, jx, iy, jy, integralerroroutside, "");

      res = res + integraloutside - integralinside;
      reserror = sqrt(std::max(0., pow(reserror, 2) + pow(integralerroroutside, 2) - pow(integralinside, 2)));
    }
  }
  if (error) *error=reserror;
  return res;
}
template <typename T> double HelperFunctions::getHistogramIntegralAndError(T const* histo, int ix, int jx, int iy, int jy, int iz, int jz, bool useWidth, double* error){
  double res=0;
  double reserror=0;
  if (histo){
    if (!useWidth) res=histo->IntegralAndError(ix, jx, iy, jy, iz, jz, reserror, "");
    else{
      int xb[2]={ std::max(1, std::min(histo->GetNbinsX(), ix)), std::max(1, std::min(histo->GetNbinsX(), jx)) };
      int yb[2]={ std::max(1, std::min(histo->GetNbinsY(), iy)), std::max(1, std::min(histo->GetNbinsY(), jy)) };
      int zb[2]={ std::max(1, std::min(histo->GetNbinsZ(), iz)), std::max(1, std::min(histo->GetNbinsZ(), jz)) };

      res=histo->IntegralAndError(xb[0], xb[1], yb[0], yb[1], zb[0], zb[1], reserror, "width");

      double integralinside, integralerrorinside;
      integralinside=histo->IntegralAndError(xb[0], xb[1], yb[0], yb[1], zb[0], zb[1], integralerrorinside, "");

      double integraloutside, integralerroroutside;
      integraloutside=histo->IntegralAndError(ix, jx, iy, jy, iz, jz, integralerroroutside, "");

      res = res + integraloutside - integralinside;
      reserror = sqrt(std::max(0., pow(reserror, 2) + pow(integralerroroutside, 2) - pow(integralinside, 2)));
    }
  }
  if (error) *error=reserror;
  return res;
}

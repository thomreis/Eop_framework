#ifndef TIMEBIN__
#define TIMEBIN__

#include <iostream>
#include <string>
#include <set>
#include <map>

#include "TH1F.h"
#include "TTree.h"

using namespace std;

class TimeBin
{
 public:

  //---ctors---
  TimeBin();
  //TimeBin(const UInt_t &runNumber, UShort_t &lumiBlock, UInt_t &time0, UInt_t &timef, int &weight);
  TimeBin(const TimeBin &bincopy);

  //---dtor---
  ~TimeBin();

  //---utils--
  void     AddEvent(const UInt_t &run, const UShort_t &ls, const UInt_t &t);
  void     AddEvent(const TimeBin& other);
  void     SetBinRanges(const UInt_t &runmin, const UInt_t &runmax, const UShort_t &lsmin, const UShort_t &lsmax, const UInt_t &timemin, const UInt_t &timemax);
  UInt_t   DeltaT() const {return timemax_-timemin_;};
  void     Reset();
  void     SetNev(const int &Nev_bin);
  bool     operator<(const TimeBin& other) const;
  void     BranchOutput(TTree* outtree);
  void     BranchInput(TTree* intree);
  void     operator=(const TimeBin& other);
  bool     Match(const UInt_t &run, const UShort_t &ls, const UInt_t &time) const;
  bool     Match(const UInt_t &run, const UShort_t &ls) const;
  void     FillHisto(double x) const {h_scale_->Fill(x);} ;
  bool     InitHisto( char* name, char* title, const int &Nbin, const double &xmin, const double &xmax);
  int      GetNev() const {return Nev_;};
  double   GetXminScale() const {return h_scale_->GetXaxis()->GetXmin();};
  double   GetXmaxScale() const {return h_scale_->GetXaxis()->GetXmax();};
  bool     Fit(TF1* fitfunc, string fitopt="QRL+", int nTrial=10, string TemplatePlotsFolder="");
  double   GetMean();
  double   GetMeanError();
  //double GetMean(double xmin, double xmax);
  //double GetMean(double evfraction);
  double   GetMedian();
  double   GetIntegral(const float &xmin, const float &xmax);
  double   GetBinWidth(const int &ibin);
  UInt_t   GetTimemin() {return timemin_;}
  UInt_t   GetTimemax() {return timemax_;}
  UInt_t   GetTime()    {return 0.5*(timemax_+timemin_);}
  UInt_t   GetIntlumimin() {return intlumimin_;}
  UInt_t   GetIntlumimax() {return intlumimax_;}
  UInt_t   GetIntlumi()    {return 0.5*(intlumimax_+intlumimin_);}
  //void   SaveAs(std::string outputfilename);  
  void     SetVariable(const std::string &variablename, const float &variablevalue);
  void     SetIntlumimin(const double &intlumi) {intlumimin_ = intlumi;}
  void     SetIntlumimax(const double &intlumi) {intlumimax_ = intlumi;}
  float    GetVariable(const std::string &variablename){return variablelist_[variablename];};
  void     PrintVariables();  
  void     UpdateNev();

 protected:
  UInt_t runmin_;
  UInt_t runmax_;
  UShort_t lsmin_;
  UShort_t lsmax_;
  UInt_t timemin_;
  UInt_t timemax_;
  double intlumimin_;
  double intlumimax_;
  int Nev_;
  TH1F* h_scale_;
  std::map<std::string,float> variablelist_;
};
#endif

#include "ICmanager.h"
#include "TMath.h"
#include "assert.h"
//#include "utils.h"

using namespace std;

ICmanager::ICmanager(CfgManager conf)
{
  InitIC(1);
  //-------------------------------------
  //load input IC
  if(conf.OptExist("Input.inputIC"))
    LoadIC( conf.GetOpt<std::vector<std::string> > ("Input.inputIC") );
  else
    cout<<"[WARNING]: no inputIC found in Input in cfg"<<endl;
  
}

ICmanager::ICmanager(const std::vector<std::string> &ICcfg)
{
  InitIC(1);
  LoadIC( ICcfg );
}

ICmanager::ICmanager()
{  
  InitIC(1);
}  

ICmanager::~ICmanager()
{ }

void ICmanager::LoadIC(TH2D* ICmap, const int &iz)
{
  timedependent_ICvalues_.clear();
  timedependent_ICvalues_.push_back( GetICFromTH2D(ICmap,iz) );
}

void ICmanager::LoadIC(const std::vector<std::string> &ICcfg)
{
  timedependent_ICvalues_.clear();

  string inputtype   = ICcfg.front();
  string filename    = ICcfg.back();
  if(inputtype=="txtIC")
  {
    cout<<"> Loading IC from txt file "<<filename<<endl;
    timedependent_ICvalues_.push_back( GetICFromtxt(filename) );
  }
  else
    if(inputtype=="txtICdictionary")
    {
      ifstream ICdictionary(filename.c_str());
      UInt_t runmin, runmax;
      UShort_t lsmin, lsmax;
      string ICfilename;
      while(!ICdictionary.eof())
      {
	if(!(ICdictionary  >> runmin >> lsmin >> runmax >> lsmax >> ICfilename))
	  continue;
	IOV thisIOV{runmin,lsmin,runmax,lsmax};
	cout<<"> Loading IC for IOV "<<runmin<<":"<<lsmin<<" - "<<runmax<<":"<<lsmax<<" from txt file "<<ICfilename<<endl;
	IOVlist_.push_back( thisIOV );
	timedependent_ICvalues_.push_back( GetICFromtxt(ICfilename) );
      }
      ICdictionary.close();
    }
    else
    {
      string objkey(inputtype);
      cout<<"> Loading IC from "<<filename<<"/"<<objkey<<endl;
      TFile* inICfile = new TFile(filename.c_str(),"READ");
      TH2D* ICmap = (TH2D*) inICfile->Get(objkey.c_str());
      int iz=0;
      if(objkey.find("EEm")!=std::string::npos || objkey.find("EEM")!=std::string::npos)
	iz=-1;
      else
	if(objkey.find("EEp")!=std::string::npos || objkey.find("EEP")!=std::string::npos)
	  iz=+1;
      timedependent_ICvalues_.push_back( GetICFromTH2D(ICmap,iz) );
      inICfile->Close();
    }
}

IC ICmanager::GetICFromtxt(const std::string &txtfilename)
{
  int ix, iy, iz;
  double icvalue, eic;
  IC icvalues;
  ifstream infile(txtfilename.c_str());
  while (!infile.eof()) 
  {
    infile >> ix >> iy >> iz >> icvalue >> eic ;
    //cout << ix <<"\t"<< iy <<"\t"<< iz <<"\t"<< icvalue << "\t"<< eic <<endl ;
    icvalues[ix][iy][iz] = icvalue;
  }
  return icvalues;
}

IC ICmanager::GetICFromTH2D(TH2D* ICmap, const int &iz)
{
  IC icvalues;
  bool toshift = (ICmap->GetXaxis()->GetXmin() == 1);
  for(int xbin=1; xbin<ICmap->GetNbinsX()+1; ++xbin)
  {
    int ix = (int) (ICmap->GetXaxis()->GetBinCenter(xbin) - 0.5*toshift);
    for(int ybin=1; ybin<ICmap->GetNbinsY()+1; ++ybin)
    {
      int iy = (int) (ICmap->GetYaxis()->GetBinCenter(ybin) - 0.5*toshift);
      if(iz==0) { //for barrel, for historical reason, in the th2f ix(ieta) and iy(iphi) are inverted 
	icvalues[iy][ix][iz] = ICmap->GetBinContent(xbin,ybin);
      }
      else
	icvalues[ix][iy][iz] = ICmap->GetBinContent(xbin,ybin);

      
    }
  }
  return icvalues;
}

void ICmanager::InitIC(Int_t ICvalue)
{
  IC icvalues;
  for(int iz=izmin_; iz<=izmax_; ++iz)
    for(int ix=ixmin_.at(iz); ix<=ixmax_.at(iz); ++ix)
      for(int iy=iymin_.at(iz); iy<=iymax_.at(iz); ++iy)
	icvalues[ix][iy][iz]=1.;

  timedependent_ICvalues_.push_back( icvalues );
}

Float_t ICmanager::GetIC(const Int_t &ix, const Int_t &iy, const Int_t &iz)
{
  //#ifdef DEBUG
  assert(iz>=izmin_     && iz<=izmax_); 
  assert(ix>=ixmin_.at(iz) && ix<=ixmax_.at(iz));
  assert(iy>=iymin_.at(iz) && iy<=iymax_.at(iz));
  //#endif
  return timedependent_ICvalues_.at(0)[ix][iy][iz];
}

Float_t ICmanager::GetIC(const Int_t &ix, const Int_t &iy, const Int_t &iz, const Int_t &iIOV)
{
  //#ifdef DEBUG
  assert(iz>=izmin_     && iz<=izmax_); 
  assert(ix>=ixmin_.at(iz) && ix<=ixmax_.at(iz));
  assert(iy>=iymin_.at(iz) && iy<=iymax_.at(iz));
  assert(iIOV<timedependent_ICvalues_.size());
  //#endif
  return timedependent_ICvalues_.at(iIOV)[ix][iy][iz];
}

int ICmanager::FindIOVNumber(const UInt_t &run, const UShort_t &ls)
{
  if(IOVlist_.size()==0)
    return 0;
  int iIOV=0;
  for(int iIOV=0; iIOV<IOVlist_.size(); ++iIOV)
    if(IOVlist_.at(iIOV).Contains(run,ls))
      return iIOV;
  //cout<<"[ERROR]: can't find the given (run,ls)=("<<run<<","<<ls<<") in the IOV list --> will use the first "<<endl;
  return -1;
}


int ICmanager::FindCloserIOVNumber(const UInt_t &run, const UShort_t &ls)
{
  int iIOV=FindIOVNumber(run,ls);
  if(iIOV>=0)
    return iIOV;

  if(IOVlist_.size()==0)
    return 0;

  if(IOVlist_.at(0).GreaterThan(run,ls))
    return 0;

  for(int iIOV=0; iIOV<IOVlist_.size()-1; ++iIOV)
    if(IOVlist_.at(iIOV).SmallerThan(run,ls) && IOVlist_.at(iIOV+1).GreaterThan(run,ls))
      return iIOV;

  return IOVlist_.size()-1;
}

TH2D* ICmanager::GetHisto(const int &iz, const char* name, const char* title)
{
  TH2D* ICmap;
  //if it is barrel, for historical reasons, I have to invert ix(ieta) and iy(iphi)
  if(iz==0)
    ICmap = new TH2D(name,title, Ny_.at(iz), iymin_.at(iz)-0.5, iymax_.at(iz)+0.5, Nx_.at(iz), ixmin_.at(iz)-0.5, ixmax_.at(iz)+0.5);
  else
    ICmap = new TH2D(name,title, Nx_.at(iz), ixmin_.at(iz)-0.5, ixmax_.at(iz)+0.5, Ny_.at(iz), iymin_.at(iz)-0.5, iymax_.at(iz)+0.5);

  ICmap->SetDirectory(0);

  for(int xbin=1; xbin<ICmap->GetNbinsX()+1; ++xbin)
    for(int ybin=1; ybin<ICmap->GetNbinsY()+1; ++ybin)
    {
      int ix,iy;
      if(iz==0)
      {
	ix=ICmap->GetYaxis()->GetBinCenter(ybin);
	iy=ICmap->GetXaxis()->GetBinCenter(xbin);
      }
      else
      {
	ix=ICmap->GetXaxis()->GetBinCenter(xbin);
	iy=ICmap->GetYaxis()->GetBinCenter(ybin);
      }
      ICmap->SetBinContent(xbin,ybin,this->GetIC(ix,iy,iz));
    }
  return ICmap;
}

double&  ICmanager::operator()(const Int_t &ix, const Int_t &iy, const Int_t &iz)
{
  return (timedependent_ICvalues_.at(0))[ix][iy][iz];
}

TH2D* ICmanager::GetPulledIC(TH2D* h2_ICpull, const int &iz)
{
  bool toshift = (h2_ICpull->GetXaxis()->GetXmin()==1);
  TH2D* pulledIC = new TH2D
    ("pulledIC","pulledIC",
     h2_ICpull->GetNbinsX(),h2_ICpull->GetXaxis()->GetXmin(),h2_ICpull->GetXaxis()->GetXmax(),
     h2_ICpull->GetNbinsY(),h2_ICpull->GetYaxis()->GetXmin(),h2_ICpull->GetYaxis()->GetXmax());
  pulledIC->SetDirectory(0);

  for(int xbin=1; xbin<h2_ICpull->GetNbinsX()+1; ++xbin)
    for(int ybin=1; ybin<h2_ICpull->GetNbinsY()+1; ++ybin)
    {
      int ix,iy;
      if(iz==0)
      {
	iy=h2_ICpull->GetXaxis()->GetBinCenter(xbin) - 0.5*toshift;
	ix=h2_ICpull->GetYaxis()->GetBinCenter(ybin) - 0.5*toshift;
      }
      else
      {   
	ix=h2_ICpull->GetXaxis()->GetBinCenter(xbin) - 0.5*toshift;
	iy=h2_ICpull->GetYaxis()->GetBinCenter(ybin) - 0.5*toshift;
      }
      double pull = h2_ICpull->GetBinContent(xbin,ybin);
      double oldIC = this->GetIC(ix,iy,iz);
      pulledIC->SetBinContent(xbin,ybin,oldIC*pull);
    }
  return pulledIC;
}

TH2D* ICmanager::PullIC(TH2D* h2_ICpull, const Int_t &iz)
{
  TH2D* pulledIC = this->GetPulledIC(h2_ICpull, iz);
  this->LoadIC(pulledIC, iz);
  return pulledIC;
}

void ICmanager::EtaringNormalizationEB()
{
  cout<<"ICmanager::EtaringNormalizationEB still to be validated"<<endl;
  /*  TGraphErrors* avgIC_vs_iEta = this->GetAvgICvsEtaEB();

  for(int ieta = ietamin_; ieta <= ietamax_; ++ieta)
  {
    double avgIC = avgIC_vs_iEta->Eval(1.*ieta);
    // normalize IC skipping bad channels and bad TTs
    for(int iphi = iphimin_; iphi <= iphimax_ ; ++iphi)
    {
      int index1D = fromIetaIphito1Dindex(ieta, iphi, Neta_, Nphi_, ietamin_, iphimin_);
      if( avgIC<0 ) continue;
      xtal_[index1D].IC /= avgIC;
    }
  }
  delete avgIC_vs_iEta;*/
}

void ICmanager::PrintSettings()
{
  cout<<"----------------------------------------------------------------------------------"<<endl;
  cout<<"> ICmanager settings:"<<endl;
  //cout<<">>> Neta="<<Neta_<<" in ["<<ietamin_<<","<<ietamax_<<"] and Nphi="<<Nphi_<<" in ["<<iphimin_<<","<<iphimax_<<"]"<<endl;
  cout<<"----------------------------------------------------------------------------------"<<endl;
}


TH1D* ICmanager::GetICspread( int nBins_spread, float spreadMin, float spreadMax)
{
  cout<<"ICmanager::GetICspread still to be validated"<<endl;
  /*  TH1D* h_ICspread = new TH1D("ICspread","ICspread",nBins_spread,spreadMin,spreadMax);
  for(int ieta = ietamin_; ieta <= ietamax_; ++ieta)
  {
    for(int iphi = iphimin_; iphi <= iphimax_ ; ++iphi)
    {
      int index1D = fromIetaIphito1Dindex(ieta, iphi, Neta_, Nphi_, ietamin_, iphimin_);
      h_ICspread->Fill(xtal_[index1D].IC);
    }
  }
  return h_ICspread;*/
  return 0;
}

TGraphErrors* ICmanager::GetICspreadvsEtaEB(int nBins_spread, float spreadMin, float spreadMax)
{
  cout<<"ICmanager::GetICspreadvsEtaEB still to be validated"<<endl;
  /*  TGraphErrors* g_ICspread_vs_iEta = new TGraphErrors();
  g_ICspread_vs_iEta->SetName("ICspread_vs_iEta");
  g_ICspread_vs_iEta->SetTitle("ICspread_vs_iEta");
  map<int,TH1F*> h_ICspread_vs_iEta;

  //Fill IC histos 
  for(int ieta = ietamax_; ieta >= ietamin_; --ieta)
  {
    if(ieta>0)
      h_ICspread_vs_iEta[ieta] = new TH1F(Form("h_ICspread_ieta%i",ieta),Form("h_ICspread_ieta%i",ieta),nBins_spread,spreadMin,spreadMax);
    for(int iphi = iphimin_; iphi <= iphimax_ ; ++iphi)
    {
      int index1D = fromIetaIphito1Dindex(ieta, iphi, Neta_, Nphi_, ietamin_, iphimin_);
      h_ICspread_vs_iEta[abs(ieta)]->Fill(xtal_[index1D].IC);
    }
  }

  //Fit IC histos and fill the graph
  TF1 fitfunc("fitfunc","gaus",spreadMin,spreadMax);
  for(int ieta = 1; ieta <= ietamax_; ++ieta)
  {
    double mean = h_ICspread_vs_iEta[ieta]->GetMean();
    double rms  = h_ICspread_vs_iEta[ieta]->GetRMS();
    fitfunc.SetParameter(1, mean);
    fitfunc.SetParameter(2, rms);
    h_ICspread_vs_iEta[ieta] -> Fit("fitfunc", "NQL", "", mean-3*rms, mean+3*rms);
    g_ICspread_vs_iEta -> SetPoint(ieta, ieta, fitfunc.GetParameter(2));
    g_ICspread_vs_iEta -> SetPointError(ieta, 0.5, fitfunc.GetParError(2));
  }
  
  for(auto h : h_ICspread_vs_iEta)
    if(h.second)
      delete h.second;

      return g_ICspread_vs_iEta;*/
  return 0;

}

TGraphErrors* ICmanager::GetAvgICvsEtaEB()
{
  cout<<"ICmanager::GetAvgICvsEtaEB still to be validated"<<endl;
  /*  TGraphErrors* avgIC_vs_iEta = new TGraphErrors();
  avgIC_vs_iEta->SetName("avgIC_vs_iEta");
  avgIC_vs_iEta->SetTitle("avgIC_vs_iEta");
  
  for(int ieta = ietamin_; ieta <= ietamax_; ++ieta)
  {
    float sumIC = 0.;
    int numIC = 0;
    
    // mean over phi corrected skipping dead channel
    for(int iphi = iphimin_; iphi <= iphimax_ ; ++iphi)
    {
      int index1D = fromIetaIphito1Dindex(ieta, iphi, Neta_, Nphi_, ietamin_, iphimin_);
      if(xtal_[index1D].status == 1)
      {
	sumIC += xtal_[index1D].IC;
	++numIC;
      }
      if(numIC!=0)
	avgIC_vs_iEta->SetPoint(ieta-ietamin_, ieta, sumIC/numIC);
      else
	avgIC_vs_iEta->SetPoint(ieta-ietamin_, ieta, 0.);
    }
  }

  return avgIC_vs_iEta;*/
  return 0;
}


TGraphErrors* ICmanager::GetICspreadvsPhiEB(int nBins_spread, float spreadMin, float spreadMax)
{

  cout<<"ICmanager::GetICspreadvsPhiEB still to be validated"<<endl;
  /*TGraphErrors* g_ICspread_vs_iphi = new TGraphErrors();
  g_ICspread_vs_iphi->SetName("ICspread_vs_iphi");
  g_ICspread_vs_iphi->SetTitle("ICspread_vs_iphi");
  map<int,TH1F*> h_ICspread_vs_iphi;

  //Fill IC histos 
  for(int iphi = iphimin_; iphi <= iphimin_; ++iphi)
  {
    h_ICspread_vs_iphi[iphi] = new TH1F(Form("h_ICspread_iphi%i",iphi),Form("h_ICspread_iphi%i",iphi),nBins_spread,spreadMin,spreadMax);
    for(int ieta = ietamin_; ieta <= ietamax_ ; ++ieta)
    {
      int index1D = fromIetaIphito1Dindex(ieta, iphi, Neta_, Nphi_, ietamin_, iphimin_);
      h_ICspread_vs_iphi[iphi]->Fill(xtal_[index1D].IC);
    }
  }

  //Fit IC histos and fill the graph
  TF1 fitfunc("fitfunc","gaus",spreadMin,spreadMax);
  for(int iphi = iphimin_; iphi <= iphimin_; ++iphi)
  {
    double mean = h_ICspread_vs_iphi[iphi]->GetMean();
    double rms  = h_ICspread_vs_iphi[iphi]->GetRMS();
    fitfunc.SetParameter(1, mean);
    fitfunc.SetParameter(2, rms);
    h_ICspread_vs_iphi[iphi] -> Fit("fitfunc", "NQL", "", mean-3*rms, mean+3*rms);
    g_ICspread_vs_iphi -> SetPoint(iphi-iphimin_, iphi, fitfunc.GetParameter(2));
    g_ICspread_vs_iphi -> SetPointError(iphi-iphimin_, 0.5, fitfunc.GetParError(2));
  }
  
  for(auto h : h_ICspread_vs_iphi)
    if(h.second)
      delete h.second;

  return g_ICspread_vs_iphi;
  */
  return 0;
}

TGraphErrors* ICmanager::GetAvgICvsPhiEB()
{
  cout<<"ICmanager::GetAvgICvsPhiEB still to be validated"<<endl;
  /*  TGraphErrors* avgIC_vs_iPhi = new TGraphErrors();
  avgIC_vs_iPhi->SetName("avgIC_vs_iPhi");
  avgIC_vs_iPhi->SetTitle("avgIC_vs_iPhi");
  
  for(int iphi = iphimin_; iphi <= iphimax_; ++iphi)
  {
    float sumIC = 0.;
    int numIC = 0;
    
    // mean over eta corrected skipping dead channel
    for(int ieta = ietamin_; ieta <= ietamax_ ; ++ieta)
    {
      int index1D = fromIetaIphito1Dindex(ieta, iphi, Neta_, Nphi_, ietamin_, iphimin_);
      if(xtal_[index1D].status == 1)
      {
	sumIC += xtal_[index1D].IC;
	++numIC;
      }
      if(numIC!=0)
	avgIC_vs_iPhi->SetPoint(iphi-iphimin_, iphi, sumIC/numIC);
      else
	avgIC_vs_iPhi->SetPoint(iphi-iphimin_, iphi, 0.);
    }
  }

  return avgIC_vs_iPhi;*/
  return 0;
}

TGraphErrors* ICmanager::GetPhiFoldProfileEB(int ietamin, int ietamax, int PhiPeriod )
{
  cout<<"ICmanager::GetPhiFoldProfileEB still to be validated"<<endl;
  /*  //checks on the inputs
  if(ietamin<ietamin_)
  {
    cout<<"[ERROR]: provided ietamin is lower than "<<ietamin_<<endl;
    return 0;
  }
  if(ietamax>ietamax_)
  {
    cout<<"[ERROR]: provided ietamax is higher than "<<ietamax_<<endl;
    return 0;
  }

  TGraphErrors* g_avgIC_vsPhiFold = new TGraphErrors();
  g_avgIC_vsPhiFold->SetName("g_avgIC_vsPhiFold");
  g_avgIC_vsPhiFold->SetTitle("g_avgIC_vsPhiFold");
    
  // define one histo per phi-region
  map<int,TH1F*> h_IC_vsPhiFold;
  for(int iregion=0; iregion<PhiPeriod; ++iregion)
    h_IC_vsPhiFold[iregion] = new TH1F(Form("h_IC_vsPhiFold_iregion%i",iregion),Form("h_IC_vsPhiFold_iregion%i",iregion),1000,0.,2.);

  // fill the histos
  for(int iphi=iphimin_; iphi<=iphimax_; ++iphi)
  {
    int iregion = (iphi-iphimin_) % PhiPeriod;
    for(int ieta=ietamin; ieta<=ietamax; ++ieta)
    {
      int index1D = fromIetaIphito1Dindex(ieta, iphi, Neta_, Nphi_, ietamin_, iphimin_);
      if(xtal_[index1D].status == 1)
	h_IC_vsPhiFold[iregion] -> Fill(xtal_[index1D].IC);
    }
  }

  //Fill the graph
  for(int iregion=0; iregion<PhiPeriod; ++iregion)
  {
    g_avgIC_vsPhiFold -> SetPoint(iregion, iregion, h_IC_vsPhiFold[iregion]->GetMean());
    g_avgIC_vsPhiFold -> SetPointError(iregion, 0.5, h_IC_vsPhiFold[iregion]->GetMeanError());
  }

  //delete histos
  for(auto h : h_IC_vsPhiFold)
    if(h.second)
      delete h.second;
  
      return g_avgIC_vsPhiFold;*/
  return 0;
}

void ICmanager::SupermoduleGapCorrectionEB(int ietamin, int ietamax, int PhiPeriod)
{
  /*  TGraphErrors* g_avgIC_vsPhiFold = this->GetPhiFoldProfileEB(ietamin,ietamax,PhiPeriod);
  //i don't want to bias the avg IC value
  double foldedICmean = TMath::Mean(g_avgIC_vsPhiFold->GetN(), g_avgIC_vsPhiFold->GetY(), NULL);

  // rescale the ICs
  for(int iphi=iphimin_; iphi<=iphimax_; ++iphi)
  {
    int iregion = (iphi-iphimin_) % PhiPeriod;
    for(int ieta=ietamin; ieta<=ietamax; ++ieta)
    {
      int index1D = fromIetaIphito1Dindex(ieta, iphi, Neta_, Nphi_, ietamin_, iphimin_);
      xtal_[index1D].IC *= foldedICmean/g_avgIC_vsPhiFold->Eval(1.*iregion);
    }
  }
  
  delete g_avgIC_vsPhiFold;
  */
}

TH2D* ICmanager::GetStatPrec(ICmanager* IC2)
{
  /*  //safety controls
  if(Neta_!=IC2->Neta_ || Nphi_!=IC2->Nphi_ || ietamin_!=IC2->ietamin_ || ietamax_!=IC2->ietamax_ || iphimin_!=IC2->iphimin_ || iphimax_!=IC2->iphimax_)
  {
    cout<<"[ERROR]::GetStatPrec incompatible ICs"<<endl;
    return 0;
  }

  TH2D* StatPrec = new TH2D("StatisticPrecision","Statistic Precision",Nphi_,iphimin_,iphimax_+1,Neta_,ietamin_,ietamax_+1);
  for(int xbin=1; xbin<StatPrec->GetNbinsX()+1; ++xbin)
    for(int ybin=1; ybin<StatPrec->GetNbinsY()+1; ++ybin)
    {
      int index = fromTH2indexto1Dindex(xbin, ybin, Nphi_, Neta_);
      double ICvalue1 = xtal_[index].IC;
      double ICvalue2 = IC2->xtal_[index].IC;
      if(xtal_[index].status==1 && IC2->xtal_[index].status==1 && ICvalue1!=0 && ICvalue2!=0)
	StatPrec->SetBinContent(xbin, ybin, (ICvalue1-ICvalue2)/(ICvalue1+ICvalue2) );
      else
	StatPrec->SetBinContent(xbin, ybin, 0.);
    }

  return StatPrec;
  */
  return 0;
}
/*
TGraphErrors* ICmanager::GetStatPrecvsEtaringEB(ICmanager* IC2)
{
  TH2D* StatPrec = this->GetStatPrec(IC2);
  TGraphErrors* g_StatPrecvsEtaring = new TGraphErrors();
  TH1D* h_StatPrec_etaRing;
  for(int ieta=1; ieta<=StatPrec->GetNbinsX(); ++ieta)
  {
    h_StatPrec_etaRing = StatPrec->Projection
  
}
*/

TH2D* GetICpull(TH2D* h2_numerator,TH2D* h2_denominator)
{
  if(h2_numerator->GetNbinsX() != h2_denominator->GetNbinsX() ||
     h2_numerator->GetNbinsY() != h2_denominator->GetNbinsY() ||
     h2_numerator->GetXaxis()->GetXmin() != h2_denominator->GetXaxis()->GetXmin() ||
     h2_numerator->GetXaxis()->GetXmax() != h2_denominator->GetXaxis()->GetXmax() ||
     h2_numerator->GetYaxis()->GetXmin() != h2_denominator->GetYaxis()->GetXmin() ||
     h2_numerator->GetYaxis()->GetXmax() != h2_denominator->GetYaxis()->GetXmax() )
  {
    cout<<"[ERROR]: GetICpull: histograms not compatible"<<endl;
    return NULL;
  }

  TH2D* ICpull = new TH2D
    ("ICpull","ICpull",
     h2_numerator->GetNbinsX(),h2_numerator->GetXaxis()->GetXmin(),h2_numerator->GetXaxis()->GetXmax(),
     h2_numerator->GetNbinsY(),h2_numerator->GetYaxis()->GetXmin(),h2_numerator->GetYaxis()->GetXmax());
  ICpull->SetDirectory(0);

  for(int xbin=1; xbin<h2_numerator->GetNbinsX()+1; ++xbin)
    for(int ybin=1; ybin<h2_numerator->GetNbinsY()+1; ++ybin)
    {
      double num = h2_numerator->GetBinContent(xbin,ybin);
      double den = h2_denominator->GetBinContent(xbin,ybin);
      if (den!=0)
	ICpull->SetBinContent(xbin,ybin,num/den);
      else
	ICpull->SetBinContent(xbin,ybin,0);
    }

  return ICpull;
}

//INDEXING UTILS!!!!
//there are the 4 different reference frames in usage, taking the barrel as example

//0 CMS reference system eta, phi

//1  ECALELF reference
//   XSeedSCEle=ieta in [-85,85] with 0 excluded for a total of 171 bins
//   YSeedSCEle=iphi in [1,360] for a total of 360 values

//2  TH2 bin numbering: for hystorical reasons phi is on the x-axis 
//   binx in [1,360] for a total of 360 bins --> actually the range is [1,361] because one has to account for the width of the last bin  
//   biny in [-85,85] for a total of 171 bins --> actually the range is [-85,86] because one has to account for the width of the last bin

//3  ix, iy auxiliary reference system 
//   ix in [0,359] for a total of 360 bins
//   iy in [0,170] for a total of 171 bins

//4  1-D index
//   360*171 bins

int fromIetaIphito1Dindex(const int &ieta, const int &iphi, const int &Neta, const int &Nphi, const int &ietamin, const int &iphimin)
{
  return (iphi - iphimin)+Nphi*(ieta - ietamin);
}

int fromTH2indexto1Dindex(const int &binx, const int &biny, const int &Nbinx, const int &Nbiny)
{
  return (binx - 1)+Nbinx*(biny - 1);
}

void from1DindextoIetaIphi(const int &index, int &ieta, int &iphi, const int &Neta, const int &Nphi, const int &ietamin, const int &iphimin)
{
  iphi = index % Nphi+iphimin;
  ieta = index / Nphi+ietamin;
}

void from1DindextoTH2index(const int &index, int &binx, int &biny, const int &Nbinx, const int &Nbiny)
{
  binx = index % Nbinx+1;
  biny = index / Nbinx+1;
}

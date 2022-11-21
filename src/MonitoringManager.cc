#include "MonitoringManager.h"
#include "histoFunc.h"
#include "FitUtils.h"

using namespace std;

MonitoringManager::MonitoringManager(CfgManager conf):
  calibrator(conf),
  conf_(conf),
  h_template_(0)
{
  label_ = conf.GetOpt<string> ("Input.label");  
  variablename_ = conf.GetOpt<string> ("LaserMonitoring.variable");
  SetScaleVariable(variablename_);
  last_accessed_bin_=timebins.end();
}

MonitoringManager::~MonitoringManager()
{
  if(h_template_)
    delete h_template_;
}

void MonitoringManager::SetScaleVariable(const string &variablename)
{
  if(variablename=="ICenergy_over_p")
  {
    cout<<">> SetScaleVariable: special keyword detected"<<endl;
    variabletype_=kICenergy_over_p;
  }
  else
    if(variablename=="ICMee")
    {
      cout<<">> SetScaleVariable: special keyword detected"<<endl;
      variabletype_=kICMee;
    }
    else
    {
      variabletype_=kregular;
      this -> AddVariable("Eop_monitoring_scale",variablename);//method of ECALELFInterface
    }
}

float MonitoringManager::GetScaleVariableValue(const int &iEle)
{
  if(variabletype_==kICenergy_over_p)
  {
    if(GetP(iEle)!=0)
      return GetICEnergy(iEle)/GetP(iEle);
    else
      return -999;
  }
  else
    if(variabletype_==kICMee)
      return GetMee()/91. *sqrt(GetICEnergy(0)*GetICEnergy(1)) /sqrt(GetEnergy(0)*GetEnergy(1));
    else
      return this -> GetVariableValue("Eop_monitoring_scale",iEle);//method of ECALELFInterface
}


TH1F* MonitoringManager::BuildTemplate()
{
  int Nbin   = conf_.GetOpt<int>   ("LaserMonitoring.BuildTemplate.Nbin");
  float xmin = conf_.GetOpt<float> ("LaserMonitoring.BuildTemplate.xmin");
  float xmax = conf_.GetOpt<float> ("LaserMonitoring.BuildTemplate.xmax");

  string templatename = Form("h_template(%i,%f,%f)",Nbin,xmin,xmax);
  
  gDirectory->cd(0);

  cout<<">> Building template \""<<variablename_<<">>"<<templatename<<"\"\n"
      <<"   with selection \""<<selection_str_<<"\"\n"
      <<">> It can take hours.."<<endl; 

  chain_->Draw( Form("%s>>%s",variablename_.c_str(),templatename.c_str()) , selection_str_.c_str() , "goff");
  h_template_=(TH1F*)gDirectory->Get("h_template");
  h_template_->SetName(("h_template_"+label_).c_str());
  h_template_->SetTitle(("h_template_"+label_).c_str());
  //curr_dir->cd();
  return h_template_;
}

void  MonitoringManager::fitScale()
{
  cout<<">> fitScale in function"<<endl;
  //Parse the cfg
  string xname = conf_.GetOpt<string>    ("LaserMonitoring.scaleFit.xname");
  string yname = conf_.GetOpt<string>    ("LaserMonitoring.scaleFit.yname");
  string yuncname = conf_.GetOpt<string> ("LaserMonitoring.scaleFit.yuncname");
  string func = conf_.GetOpt<string>      ("LaserMonitoring.scaleFit.func");
  string label = conf_.GetOpt<string>      ("LaserMonitoring.scaleFit.label");
  string fitopt  = "QRL+";
  int Ntrialfit  = 10;
  string TemplatePlotsFolder="";
  if(conf_.OptExist("LaserMonitoring.scaleFit.fitoptions"))
     fitopt = conf_.GetOpt<string> ("LaserMonitoring.scaleFit.fitoptions");
  if(conf_.OptExist("LaserMonitoring.scaleFit.Ntrialfit"))
    Ntrialfit = conf_.GetOpt<int> ("LaserMonitoring.scaleFit.Ntrialfit");
  if(conf_.OptExist("LaserMonitoring.scaleFit.fitplots_folder"))
    TemplatePlotsFolder=conf_.GetOpt<string> ("LaserMonitoring.scaleFit.fitplots_folder");

  //Build the TGraph
  cout<<">>> building the graph"<<endl;
  TGraphErrors* g_scale_vs_x = new TGraphErrors();
  g_scale_vs_x->SetName(label.c_str());
  g_scale_vs_x->SetTitle(label.c_str());
  for(std::vector<TimeBin>::iterator it_bin = timebins.begin(); it_bin<timebins.end(); ++it_bin)
  {
    double x=-999;
    double ex=0.;
    if(xname=="time")
    {
      x = it_bin -> GetTime();
      ex = (it_bin->GetTimemax()-it_bin->GetTimemin())/sqrt(12.);
    }

    double y = it_bin->GetVariable(yname);
    if (y>2. || y<0.) continue;
    double ey=0.;
    if(yuncname!="")
      ey = it_bin->GetVariable(yuncname);

    g_scale_vs_x -> SetPoint( g_scale_vs_x->GetN(), x, it_bin->GetVariable(yname));
    g_scale_vs_x -> SetPointError( g_scale_vs_x->GetN()-1, ex, ey);
  }

  //Build the fit function
  cout<<">>> building the fit function"<<endl;
  float xmin_fit=0;
  if(conf_.OptExist("LaserMonitoring.scaleFit.xmin_fit"))
    xmin_fit = conf_.GetOpt<float>   ("LaserMonitoring.scaleFit.xmin_fit");
  else
    if(xname=="time")
      xmin_fit = timebins.begin()->GetTimemin()-1; 
  float xmax_fit=0;
  if(conf_.OptExist("LaserMonitoring.scaleFit.xmax_fit"))
    xmax_fit = conf_.GetOpt<float>   ("LaserMonitoring.scaleFit.xmax_fit");
  else
    if(xname=="time")
      xmax_fit = timebins.begin()->GetTimemax()-1;
  cout<<">>> building fit function "<<func<<" defined in range ("<<xmin_fit<<","<<xmax_fit<<")"<<endl; 
  TF1* fitfunc = new TF1(Form("fitfunc_%s",label.c_str()), func.c_str(), xmin_fit, xmax_fit);
  
  //Fit the graph
  cout<<">>> fitting"<<endl; 
  bool isgoodfit = FitUtils::PerseverantFit(g_scale_vs_x, fitfunc, fitopt, Ntrialfit, TemplatePlotsFolder);

  //Add the resulting parameters to the timebin
  for(std::vector<TimeBin>::iterator it_bin = timebins.begin(); it_bin<timebins.end(); ++it_bin)
  {
    //cout<<"reading bin "<<it_bin-timebins.begin()<<endl;
    for( int iPar=0; iPar<fitfunc->GetNpar(); ++iPar)
    {
      it_bin->SetVariable(Form("scale_fit_%s_par%i",label.c_str(),iPar), fitfunc->GetParameter(iPar));
      it_bin->SetVariable(Form("scale_unc_fit_%s_par%i",label.c_str(),iPar), fitfunc->GetParError(iPar));
    }
  }

  delete fitfunc;
  delete g_scale_vs_x;
  
}

void  MonitoringManager::RunDivide()
{
  cout<<">> Running RunDivide"<<endl;
  int    Nevmax_bin     = conf_.GetOpt<int>          ("LaserMonitoring.RunDivide.Nevmax_bin");
  float  maxduration    = 60*60*conf_.GetOpt<float>  ("LaserMonitoring.RunDivide.maxduration");//It is provided in hours

  //Loop on events to build bins corresponding to single lumisections, they will be merged afterward
  //Exploit methods inherited from ECALELFInterface to access the ntuple content
  long Nentries = this->GetEntries();
  cout<<Nentries<<" total entries\n"<<endl;
  for(long ientry=0; ientry<Nentries; ++ientry)
  {
    this->GetEntry(ientry);
    if(ientry%100000==0)
      cout<<"reading entry "<<ientry<<"\r"<<std::flush;
    int w=0;
    if(this->isSelected(0)) ++w;
    if(this->isSelected(1)) ++w;
    if(w)
    {
      unsigned t = this->GetTime();
      auto run =   this->GetRunNumber();
      auto ls =    this->GetLS();
      auto bin_iterator = FindBin(run,ls);
      if(bin_iterator!=timebins.end())
	bin_iterator->AddEvent(run,ls,t);
      else
      {
	TimeBin newbin;
	newbin.SetBinRanges(run,run,ls,ls,t,t);
	newbin.SetNev(w);
	timebins.push_back(newbin);
	last_accessed_bin_ = timebins.end();
      }
    }
  }
  cout<<endl;

  //Merge the lumisections to create TimeBins with about the required number of events
  cout<<">> Merging lumisections"<<endl;
  std::sort(timebins.begin(), timebins.end());
  vector<TimeBin> ls_bins;
  timebins.swap(ls_bins);//switch the contents of timebins vector with the one of ls_bins vector(empty)
  auto lsbin_iterator = ls_bins.begin();
  TimeBin bufferbin(*lsbin_iterator);
  lsbin_iterator++;
  while(lsbin_iterator != ls_bins.end())
  {
    if(bufferbin.DeltaT() > maxduration)
    {
      timebins.push_back(bufferbin);
      bufferbin.Reset();
    }
    else
      if(bufferbin.GetNev() >= Nevmax_bin)
      {
	timebins.push_back(bufferbin);
	bufferbin.Reset();
      }
      else
      {
	bufferbin.AddEvent(*lsbin_iterator);
	lsbin_iterator++;
      }
  }

  if(bufferbin.GetNev() >= Nevmax_bin*0.5)
    timebins.push_back(bufferbin);

  std::sort(timebins.begin(), timebins.end());

  if(conf_.OptExist("LaserMonitoring.RunDivide.lumifilename"))
  {
    string  intlumi_vs_time_filename   = conf_.GetOpt<string>  ("LaserMonitoring.RunDivide.lumifilename");
    LoadIntegratedLuminosity(intlumi_vs_time_filename);
  }

  last_accessed_bin_ = timebins.end();

}

void  MonitoringManager::LoadIntegratedLuminosity(string intlumi_vs_time_filename)
{
  cout<<"Loading integrated lumi info from "<<intlumi_vs_time_filename<<endl;
  TFile* infile_intlumi_vs_time = new TFile(intlumi_vs_time_filename.c_str());
  TGraph* g_intlumi_vs_time = (TGraph*)infile_intlumi_vs_time->Get("g_intlumi_time");
  infile_intlumi_vs_time->Close();
  if(!g_intlumi_vs_time)
    cout<<"[ERROR]: cannot load g_intlumi_time from "<<intlumi_vs_time_filename<<endl;
  else
    for(std::vector<TimeBin>::iterator it_bin = timebins.begin(); it_bin<timebins.end(); ++it_bin)
    {
      double timemin = 1.*it_bin->GetTimemin(); 
      it_bin->SetIntlumimin( g_intlumi_vs_time ->Eval(timemin) ); 
      double timemax = 1.*it_bin->GetTimemax(); 
      it_bin->SetIntlumimax( g_intlumi_vs_time ->Eval(timemax) ); 
    }
}


void  MonitoringManager::SaveTimeBins(std::string outfilename, std::string writemethod)
{
  TFile* outfile = new TFile(outfilename.c_str(), writemethod.c_str());
  outfile->cd();
  TTree* outtree = new TTree(label_.c_str(), label_.c_str());
  TimeBin bin( *(timebins.begin()) );
  bin.BranchOutput(outtree);
  
  for(auto bincontent : timebins)
  {
    if(bincontent.GetNev() > 0)
    {
      bin=bincontent;
      outtree->Fill();
    }
  }

  outfile->cd();
  outtree->AutoSave();
  outfile->Close();
}

// create list of time bins from list of runs
// this function is meant to use by the automatic prompt calibration system
void  MonitoringManager::LoadTimeBins(std::vector<UInt_t>& runs, std::vector<UInt_t>& times, std::string option)
{
  cout<<">> Loading timebins"<<endl; 
  if(timebins.size()>0)
    if(option=="RELOAD")
      timebins.clear();
    else
    {
      cout<<"[WARNING]: timebins not loaded because already in memory"<<endl
	  <<"           if you want to overwrite call LoadTimeBins(inputfilename, objname, \"RELOAD\")"<<endl;
      return;
    }
  if(runs.size() < 2)
    {
      cout<<"[ERROR]: Not enough runs specified in LoadTimeBins (at least two required)"<<endl;
      return;
    }
    
  for(unsigned int i=1; i<runs.size(); ++i)
  {      
    TimeBin ibin;
    ibin.SetBinRanges(runs[i-1], runs[i], 0, std::numeric_limits<bool>::max(), times[i-1], times[i]);
    timebins.push_back(ibin);
  }
  std::sort(timebins.begin(), timebins.end());//It should be already ordered, just for security
  cout<<">> Loaded "<<timebins.size()<<" bins"<<endl;
  last_accessed_bin_=timebins.end();    
}

void  MonitoringManager::LoadTimeBins(string inputfilename, string objname, std::string option)
{
  cout<<">> Loading timebins"<<endl; 
  if(timebins.size()>0)
    if(option=="RELOAD")
      timebins.clear();
    else
    {
      cout<<"[WARNING]: timebins not loaded because already in memory"<<endl
	  <<"           if you want to overwrite call LoadTimeBins(inputfilename, objname, \"RELOAD\")"<<endl;
      return;
    }

  TFile* inputfile;
  TTree* intree;
  if(objname=="")
  {
    objname=label_;
    cout<<">> Loading TTree "<<objname<<" (default name) from file "<<inputfilename<<endl;
  }
  else
    cout<<">> Loading TTree "<<objname<<" from file "<<inputfilename<<endl;

  inputfile = new TFile(inputfilename.c_str(),"READ");
  intree = (TTree*) inputfile->Get(objname.c_str());
  
  if(!intree)
  {
    cout<<"[ERROR]: can't get tree "<<label_<<" in the file"<<endl;
    return;
  }
  TimeBin bin;
  bin.BranchInput(intree);
  Long64_t Nbins = intree->GetEntries();
  for(Long64_t ibin=0; ibin<Nbins; ++ibin)
  {
    intree->GetEntry(ibin);//tree entry is copied in bin
    TimeBin bincopy(bin);//perhaps not needed, but for security I avoid to make a mess with pointers
    timebins.push_back(bincopy);
  }
  
  std::sort(timebins.begin(), timebins.end());//It should be already ordered, just for security
  cout<<">> Loaded "<<timebins.size()<<" bins"<<endl;
  inputfile->Close();
  last_accessed_bin_=timebins.end();
}

bool MonitoringManager::BookHistos()
{
  int Nbin_histos = conf_.GetOpt<int>      ("LaserMonitoring.scaleMonitor.Nbin_histos");
  float xmin_histos = conf_.GetOpt<float>  ("LaserMonitoring.scaleMonitor.xmin_histos");
  float xmax_histos = conf_.GetOpt<float>  ("LaserMonitoring.scaleMonitor.xmax_histos");
  for(unsigned ibin=0; ibin<timebins.size(); ++ibin)
    if(!timebins.at(ibin).InitHisto( Form("Histo%i",ibin), Form("Histo%i",ibin), Nbin_histos, xmin_histos, xmax_histos))
      return false;

  return true;
}

//Loop over ECALELF tree to fill the timebins with the corresponding values
void  MonitoringManager::FillTimeBins()
{
  cout<<">> Filling timebin histos with variable "<<variablename_<<endl;
  if(!BookHistos())
  {
    cout<<">> Cannot book the histos... Maybe you have already filled them " << endl;
    return;
  }

  
  long Nentries = this->GetEntries();
  cout<<Nentries<<" total entries\n"<<endl;
  for(long ientry=0; ientry<Nentries; ++ientry)
  {
    this->GetEntry(ientry);
    if(ientry%100000==0)
      cout<<"reading entry "<<ientry<<"\r"<<std::flush;

    for(int iEle=0; iEle<2; ++iEle)
    {
      if(this->isSelected(iEle))
      {
	int iIOV = FindIOVNumber(GetRunNumber(),GetLS());
	//printf("iIOV=%i\tE=%.1f\tICEnergy=%.1f\tp=%.1f\tE/p=%.1f\tscale=%.1f",iIOV,GetEnergy(iEle),GetICEnergy(iEle),GetP(iEle),GetICEnergy(iEle)/GetP(iEle),GetScaleVariableValue(iEle));
	//getchar();
	auto bin_iterator = FindBin(this->GetRunNumber(),this->GetLS(),this->GetTime());
	if(bin_iterator!=timebins.end())
	  bin_iterator->FillHisto( GetScaleVariableValue(iEle) );
      }
    }
  }

  cout<<">> Histos filled"<<endl;

  //Updating Nev of the bins (just a precaution)
  for(std::vector<TimeBin>::iterator it_bin = timebins.begin(); it_bin<timebins.end(); ++it_bin)
    it_bin->UpdateNev();

}

std::vector<TimeBin>::iterator MonitoringManager::FindBin(const UInt_t &run, const UShort_t &ls)
{
  //cout<<"finding bin"<<endl;
  //usually events are in the same time bin of the previous iteration or in adjacent timebins so i start to look for them from there
  std::vector<TimeBin>::iterator it_end   = timebins.end();
  std::vector<TimeBin>::iterator it_begin = timebins.begin();

  if(last_accessed_bin_!=it_end)
  {
    if(last_accessed_bin_ -> Match(run,ls))
      return last_accessed_bin_;

    if(last_accessed_bin_>it_begin)
      if((last_accessed_bin_-1) -> Match(run,ls))
      {
	last_accessed_bin_--;
	return (last_accessed_bin_);
      }

    if(last_accessed_bin_<it_end-1)
       if((last_accessed_bin_+1) -> Match(run,ls))
       {
	 last_accessed_bin_++;
	 return (last_accessed_bin_);
       }
  }    
  //if I am here, unfortunately I have to perform a search through the entire set
  //cout<<"last_accessed_bin_=it_end"<<endl;
  for(std::vector<TimeBin>::iterator it_bin = it_begin; it_bin<it_end; ++it_bin)//not the smarter way considering that the bins are ordered --> can be improved
    if(it_bin -> Match(run,ls))
    {
      //cout<<"match with bin "<<it_bin-it_begin<<endl;
      last_accessed_bin_=it_bin;
      return (last_accessed_bin_);
    }

  //if I am here, I didn't found any match
  return it_end;
}


std::vector<TimeBin>::iterator MonitoringManager::FindBin(const UInt_t &run, const UShort_t &ls, const UInt_t &time)
{
  //cout<<"finding bin"<<endl;
  //usually events are in the same time bin of the previous iteration or in adjacent timebins so i start to look for them from there
  std::vector<TimeBin>::iterator it_end   = timebins.end();
  std::vector<TimeBin>::iterator it_begin = timebins.begin();

  if(last_accessed_bin_!=it_end)
  {
    //cout<<"last_accessed_bin_!=it_end"<<endl;
    if(last_accessed_bin_ -> Match(run,ls,time))
      return last_accessed_bin_;

    if(last_accessed_bin_>it_begin)
      if((last_accessed_bin_-1) -> Match(run,ls,time))
      {
	last_accessed_bin_--;
	return (last_accessed_bin_);
      }

    if(last_accessed_bin_<it_end-1)
       if((last_accessed_bin_+1) -> Match(run,ls,time))
       {
	 last_accessed_bin_++;
	 return (last_accessed_bin_);
       }
  }    
  //if I am here, unfortunately I have to perform a search through the entire set
  //cout<<"last_accessed_bin_=it_end"<<endl;
  for(std::vector<TimeBin>::iterator it_bin = it_begin; it_bin<it_end; ++it_bin)//not the smarter way considering that the bins are ordered --> can be improved
    if(it_bin -> Match(run,ls,time))
    {
      //cout<<"match with bin "<<it_bin-it_begin<<endl;
      last_accessed_bin_=it_bin;
      return (last_accessed_bin_);
    }

  //if I am here, I didn't found any match
  return it_end;
}
	     
void  MonitoringManager::RunTemplateFit(string scale)
{
  cout<<">> RunTemplateFit in function"<<endl;
  if(h_template_)
  {
    cout<<"[WARNING]: a template histogram is already loaded in memory  deleting it"<<endl;
    delete h_template_;
  }
  //Load the template histogram
  vector<string> templatename = conf_.GetOpt<vector<string> > (Form("LaserMonitoring.scaleMonitor.%s.template",scale.c_str()));
  string templatekeyname = templatename.at(0);
  string templatefilename = templatename.at(1);
  TFile* templatefile = new TFile(templatefilename.c_str(),"READ");
  h_template_ = (TH1F*) templatefile->Get(templatekeyname.c_str());
  h_template_ ->SetDirectory(0);
  templatefile->Close();

  //Build the TF1 from the template histogram
  float xmin_fit = conf_.GetOpt<float>  (Form("LaserMonitoring.scaleMonitor.%s.xmin_fit",scale.c_str()));
  float xmax_fit = conf_.GetOpt<float>  (Form("LaserMonitoring.scaleMonitor.%s.xmax_fit",scale.c_str()));
  string fitopt  = "QRL+";
  int Ntrialfit  = 10;
  string TemplatePlotsFolder="";
  if(conf_.OptExist(Form("LaserMonitoring.scaleMonitor.%s.fitoptions",scale.c_str())))
    fitopt = conf_.GetOpt<string> (Form("LaserMonitoring.scaleMonitor.%s.fitoptions",scale.c_str()));
  if(conf_.OptExist(Form("LaserMonitoring.scaleMonitor.%s.Ntrialfit",scale.c_str())))
    Ntrialfit = conf_.GetOpt<int> (Form("LaserMonitoring.scaleMonitor.%s.Ntrialfit",scale.c_str()));
  if(conf_.OptExist(Form("LaserMonitoring.scaleMonitor.%s.fitplots_folder",scale.c_str())))
    TemplatePlotsFolder=conf_.GetOpt<string> (Form("LaserMonitoring.scaleMonitor.%s.fitplots_folder",scale.c_str()));

  histoFunc* templateHistoFunc = new histoFunc(h_template_);
  TF1* fitfunc = new TF1("fitfunc",templateHistoFunc, xmin_fit, xmax_fit, 3, "histofunc");
  fitfunc -> SetParName(0, "Norm");
  fitfunc -> SetParName(1, "Scale factor");
  fitfunc -> SetLineWidth(1);
  fitfunc -> SetNpx(10000);
  fitfunc -> SetLineColor(kGreen + 2);

  double templateIntegral = h_template_->Integral(h_template_->GetXaxis()->FindBin(xmin_fit), h_template_->GetXaxis()->FindBin(xmax_fit));
  //Run the fits
  for(std::vector<TimeBin>::iterator it_bin = timebins.begin(); it_bin<timebins.end(); ++it_bin)
  {
    double binwidthRatio = it_bin->GetBinWidth(1) / h_template_->GetBinWidth(1); 
    if(xmin_fit <= it_bin->GetXminScale() || xmax_fit >= it_bin->GetXmaxScale())
      cout<<"[WARNING]: the fit range is wider or equal to the histogram range! Possible issues in the normalization or fitting..."<<endl;
    double xNorm = it_bin->GetIntegral(xmin_fit,xmax_fit) / templateIntegral * binwidthRatio;
    fitfunc -> FixParameter(0, xNorm);
    fitfunc -> SetParameter(1, 0.99);
    fitfunc -> FixParameter(2, 0.);
    
    //cout<<"prefit fitfunc integral = "<<fitfunc->Integral(xmin_fit,xmax_fit)<<endl;
    //cout<<"reading bin "<<it_bin-timebins.begin()<<endl;
    bool isgoodfit = it_bin->Fit(fitfunc,fitopt,Ntrialfit,TemplatePlotsFolder);
    float fitscale = fitfunc->GetParameter(1);
    if(isgoodfit && fitscale!=0)
    {
      it_bin->SetVariable("scale_"+scale,    1./fitscale);
      if (fitfunc->GetParError(1) / fitscale / fitscale > 0.001)  it_bin->SetVariable("scale_unc_"+scale, fitfunc->GetParError(1) / fitscale / fitscale);
      else  it_bin->SetVariable("scale_unc_"+scale, 0.001 );

    }
    else
    {
      it_bin->SetVariable("scale_"+scale,    -999);
      it_bin->SetVariable("scale_unc_"+scale, 0);
    }
    //cout<<"postfit fitfunc integral = "<<fitfunc->Integral(xmin_fit,xmax_fit)<<endl;
  }

  delete fitfunc;
  delete templateHistoFunc;
}
		 
	  
void  MonitoringManager::RunComputeMean(string scale)
{
  cout<<">> RunComputeMean in function"<<endl;
  for(std::vector<TimeBin>::iterator it_bin = timebins.begin(); it_bin<timebins.end(); ++it_bin)
  {
    //cout<<"reading bin "<<it_bin-timebins.begin()<<endl;
    if(it_bin->GetNev()==0)
    {
      cout<<"[ERROR]: bin "<<it_bin-timebins.begin()<<" is empty"<<endl;
      it_bin->SetVariable("scale_"+scale, -999.);
      it_bin->SetVariable("scale_unc_"+scale, 0.);
    }
    else
    {
      it_bin->SetVariable("scale_"+scale, it_bin->GetMean());
      it_bin->SetVariable("scale_unc_"+scale, it_bin->GetMeanError());
    }
  }
}

	  
void  MonitoringManager::RunComputeMedian(string scale)
{
  cout<<">> RunComputeMedian in function"<<endl;
  for(std::vector<TimeBin>::iterator it_bin = timebins.begin(); it_bin<timebins.end(); ++it_bin)
  {
    //cout<<"reading bin "<<it_bin-timebins.begin()<<endl;
    if(it_bin->GetNev()==0)
    {
      cout<<"[ERROR]: bin "<<it_bin-timebins.begin()<<" is empty"<<endl;
      it_bin->SetVariable("scale_"+scale, -999.);
      it_bin->SetVariable("scale_unc_"+scale, 0.);
    }
    else
    {
      it_bin->SetVariable("scale_"+scale, it_bin->GetMedian());
      it_bin->SetVariable("scale_unc_"+scale, it_bin->GetMeanError());
    }
  }
}

void  MonitoringManager::PrintScales()
{
  auto firstbin=timebins.begin();
  
}



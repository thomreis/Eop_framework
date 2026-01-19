#include "ECALELFInterface.h"
#include "TObjArray.h"
#include "TChainElement.h"
#include "TList.h"
#include "TFriendElement.h"
#include "TROOT.h"

using namespace std;

void ECALELFInterface::BranchSelected(TChain* chain)
{
  chain->SetBranchAddress("runNumber",          &runNumber_);
  chain->SetBranchAddress("lumiBlock",          &lumiBlock_);
  chain->SetBranchAddress("eventTime",          &eventTime_);
  chain->SetBranchAddress("eventNumber",        &eventNumber_);
  chain->SetBranchAddress("chargeEle",          chargeEle_);
  chain->SetBranchAddress("etaEle",             etaEle_);
  chain->SetBranchAddress("phiEle",             phiEle_);
  chain->SetBranchAddress("rawEnergySCEle",     rawEnergySCEle_);
  chain->SetBranchAddress("energy_ECAL_ele",    energySCEle_);
  chain->SetBranchAddress("invMass_ECAL_ele",   &Mee_);
  chain->SetBranchAddress("etaSCEle",           etaSCEle_);
  chain->SetBranchAddress("phiSCEle",           phiSCEle_);
  chain->SetBranchAddress("esEnergySCEle",      esEnergySCEle_);
  chain->SetBranchAddress("pAtVtxGsfEle",       pAtVtxGsfEle_);
  chain->SetBranchAddress("fbremEle",           fbremEle_);
  chain->SetBranchAddress("xSeedSC",            xSeed_);
  chain->SetBranchAddress("ySeedSC",            ySeed_);
}


void ECALELFInterface::BranchExtraCalib(TChain* chain)
{
  for(int i=0;i<2;++i)
  {
    ERecHit_[i]=0;
    XRecHit_[i]=0;  //ETA
    YRecHit_[i]=0;  //PHI
    ZRecHit_[i]=0;
    recoFlagRecHit_[i]=0;
    fracRecHit_[i]=0;
  }

  // ele1
  chain->SetBranchAddress("energyRecHitSCEle1",   &ERecHit_[0]);
  chain->SetBranchAddress("XRecHitSCEle1",        &XRecHit_[0]);
  chain->SetBranchAddress("YRecHitSCEle1",        &YRecHit_[0]);
  chain->SetBranchAddress("ZRecHitSCEle1",        &ZRecHit_[0]);
  chain->SetBranchAddress("recoFlagRecHitSCEle1", &recoFlagRecHit_[0]);
  chain->SetBranchAddress("fracRecHitSCEle1",     &fracRecHit_[0]);
  // ele2
  chain->SetBranchAddress("energyRecHitSCEle2",   &ERecHit_[1]);
  chain->SetBranchAddress("XRecHitSCEle2",        &XRecHit_[1]);
  chain->SetBranchAddress("YRecHitSCEle2",        &YRecHit_[1]);
  chain->SetBranchAddress("ZRecHitSCEle2",        &ZRecHit_[1]);
  chain->SetBranchAddress("recoFlagRecHitSCEle2", &recoFlagRecHit_[1]);
  chain->SetBranchAddress("fracRecHitSCEle2",     &fracRecHit_[1]);
}


ECALELFInterface::ECALELFInterface(CfgManager conf):
  eeRing_(0),
  selection_(0)
{
  //-------------------------------------
  //initialize chain and branch tree
  std::vector<std::string> treelist = conf.GetOpt<std::vector<std::string> >("Input.treelist");
  for(auto treename : treelist)
  { 
    ch_[treename] = new TChain(treename.c_str(),treename.c_str());
    std::vector<std::string> filelist = conf.GetOpt<std::vector<std::string> >(Form("Input.%s.filelist",treename.c_str()));
    for(auto filename : filelist)
      ch_[treename]->Add(filename.c_str());
    if(treename=="selected")
      BranchSelected(ch_[treename]);
    else
      if(treename=="extraCalibTree")
	BranchExtraCalib(ch_[treename]);
      else
	cerr<<"[WARNING]: unknown tree "<<treename<<endl;
  }

  auto Nentries = ch_[treelist.at(0)]->GetEntries(); 
  for(unsigned int nchain = 1; nchain < treelist.size(); ++nchain)
  {
    cout << ">>> Adding chain " << treelist.at(nchain) << " as friend to chain " << treelist.at(0) << endl;
    assert(Nentries == ch_[treelist.at(nchain)]->GetEntries());
    ch_[treelist.at(0)]->AddFriend(ch_[treelist.at(nchain)],"");
    //ch_[treelist.at(0)]->BuildIndex("runNumber","eventNumber");
  }
  chain_=ch_[treelist.at(0)];
  Ncurrtree_=1;

  //-------------------------------------
  //load event selection
  this->SetSelection( conf.GetOpt<string> ("Input.selection") );

  //-------------------------------------
  //initialize EEring
  string EEringsFileName = conf.GetOpt<string> ("Input.eeringsFileName");
  eeRing_ = new TEndcapRings(EEringsFileName);
}

ECALELFInterface::~ECALELFInterface()
{
  if(selection_)
    delete selection_;

  for(auto ch_iterator : ch_)
    if(ch_iterator.second)
      (ch_iterator.second)->Delete();

  if(eeRing_)
    delete eeRing_;

  for(auto customvariablesiterator : customvariablesmap_)
    if(customvariablesiterator.second)
      delete customvariablesiterator.second;

}

Long64_t ECALELFInterface::GetEntry(const Long64_t &entry)
{
  Long64_t i=chain_->GetEntry(entry);
  if(chain_->GetTreeNumber() != Ncurrtree_)
  {
    Ncurrtree_ = chain_->GetTreeNumber();
    selection_->UpdateFormulaLeaves();
    for(auto customvariablesiterator : customvariablesmap_)
      customvariablesiterator.second -> UpdateFormulaLeaves();
  }
  return i;
}

Bool_t  ECALELFInterface::isEB(const Int_t &i)
{
  if(fabs(etaSCEle_[i])<1.479)
    return true;
  return false;
}

Bool_t  ECALELFInterface::isEE(const Int_t &i)
{
  if(fabs(etaSCEle_[i])>1.479 && fabs(etaSCEle_[i])<2.5)
    return true;
  return false;
}

void ECALELFInterface::GetSeed(Int_t &ix, Int_t &iy, const Int_t &i)
{
  ix=xSeed_[i];
  iy=ySeed_[i];
}

int ECALELFInterface::GetietaSeed(const Int_t &i)
{
  if(this->isEB(i))
    return xSeed_[i];
  else
    return eeRing_->GetEndcapIeta(xSeed_[i], ySeed_[i], (int)(etaSCEle_[i]>0) );
}

int ECALELFInterface::GetiphiSeed(const Int_t &i)
{
  if(this->isEB(i))
    return ySeed_[i];
  else
    return eeRing_->GetEndcapIphi(xSeed_[i], ySeed_[i], (int)(etaSCEle_[i]>0) );
}

int ECALELFInterface::GetEERingSeed(const Int_t &i)
{
  if( ! this->isEB(i) )
    return eeRing_->GetEndcapRing(xSeed_[i], ySeed_[i], (int)(etaSCEle_[i]>0) );
  else
    return -999;
}
  
void ECALELFInterface::PrintEleSummary    (const Int_t &i)
{
  cout<<std::left
      <<std::setw(15)<<"event "  <<eventNumber_
      <<std::setw(15)<<"run "    <<runNumber_
      <<std::setw(15)<<"ls "     <<lumiBlock_
      <<std::setw(6) <<"ele "    <<i
      <<std::setw(13)<<"charge " <<chargeEle_[i]
      <<std::setw(13)<<"E "      <<energySCEle_[i]
      <<std::setw(13)<<"p "      <<pAtVtxGsfEle_[i]
      <<std::setw(10)<<"etaSC "  <<etaSCEle_[i]
      <<std::setw(13)<<"xSeed "  <<xSeed_[i]
      <<std::setw(13)<<"ySeed "  <<ySeed_[i]
      <<endl;
}

void ECALELFInterface::PrintRHEleSummary  (const Int_t &i)
{    
  cout<<"ele "<<i<<endl;
  cout<<std::left
      <<std::setw(10)<<"x "
      <<std::setw(10)<<"y "
      <<std::setw(10)<<"z "
      <<std::setw(10)<<"E "
      <<std::setw(15)<<"frac "
      <<std::setw(15)<<"flag "
      <<endl;
 
  for(unsigned iRH=0; iRH<ERecHit_[i]->size(); ++iRH)
    cout<<std::left
	<<std::setw(10)<<XRecHit_[i]
	<<std::setw(10)<<YRecHit_[i]
	<<std::setw(10)<<ZRecHit_[i]
	<<std::setw(10)<<ERecHit_[i]
	<<std::setw(15)<<fracRecHit_[i]
	<<std::setw(15)<<recoFlagRecHit_[i]
	<<endl;

}

void ECALELFInterface::SetSelection(string selection_str)
{
  if(selection_)
    delete selection_;
  selection_str_ = selection_str;
  selection_ = new TTreeFormula("selection", selection_str_.c_str(), chain_);
} 


void ECALELFInterface::AddSelection(string additional_selection_str)
{
  if(selection_)
    delete selection_;

  if(selection_str_=="")
    selection_str_ = additional_selection_str;
  else
  {
    selection_str_.insert(0,"(");
    selection_str_.append(") && ("+additional_selection_str+")");
  }

  selection_ = new TTreeFormula("selection", selection_str_.c_str(), chain_);
} 

void ECALELFInterface::PrintSettings()
{
  cout<<"----------------------------------------------------------------------------------"<<endl;
  cout<<"> ECALELFInterface settings:"<<endl;

  cout<<">>> MAIN CHAIN NAME: "<<chain_->GetName()<<endl; 
  TObjArray* filelist = chain_->GetListOfFiles();
  TIter next(filelist);
  TChainElement *element;
  cout<<">>> FILENAMES: "<<endl;
  while ((element = (TChainElement*)next())) 
    cout<<">>> \t"<<element->GetTitle()<<endl;

  TList* friendlist = chain_->GetListOfFriends();
  TIter nextf(friendlist);
  TFriendElement *fr;
  while ((fr = (TFriendElement*)nextf())) 
  {
    cout<<endl<<">>> FRIEND CHAIN NAME: "<<fr->GetTreeName()<<endl;
    TChain* fc = (TChain*)gROOT->FindObject(fr->GetTreeName());
    if(fc)
    {
      cout<<">>> FILENAMES: "<<endl;
      TObjArray* ffilelist = fc->GetListOfFiles();
      TIter ffnext(ffilelist);
      TChainElement *felement;
      while ((felement = (TChainElement*)ffnext()))
      cout<<">>> \t"<<felement->GetTitle()<<endl;
    }
  }

  cout<<"> APPLIED SELECTION: "<<selection_->GetExpFormula().Data()<<endl;
  cout<<"----------------------------------------------------------------------------------"<<endl;
}

void ECALELFInterface::AddVariable(const string &name, const string &expr)
{
  customvariablesmap_[name] = new TTreeFormula(name.c_str(), expr.c_str(), chain_);
}

double ECALELFInterface::GetVariableValue(const string &name, const Int_t &i)
{
  return customvariablesmap_[name] -> EvalInstance(i);
}

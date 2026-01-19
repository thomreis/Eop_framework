#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH2D.h"
#include "TObject.h"

using namespace std;

//#define DEBUG

void PrintUsage()
{
  cerr << ">>>>> usage:  UpdateIC --oldIC <objname> <filename> --newIC <objname> <filename>" << endl;
}

int main(int argc, char* argv[])
{
  //Parse input parameters
  TString oldFileName="";
  TString oldObjName="";
  TString newFileName="";
  TString newObjName="";

  for(int iarg=1; iarg<argc; ++iarg)
  {
    if(string(argv[iarg])=="--oldIC")
    {
      oldObjName=argv[iarg+1];
      oldFileName=argv[iarg+2];

    }
    if(string(argv[iarg])=="--newIC")
    {
      newObjName=argv[iarg+1];
      newFileName=argv[iarg+2];
    }
  }

  if(oldFileName=="")
    cout<<"[WARNING] no old ICs provided --> Assume oldIC equal to 1 by default"<<endl;
  else
    if(oldObjName=="" || newFileName=="" || newObjName=="")
    {
      PrintUsage();
      return -1;
    }

  //Load old ICs
  TH2D* oldIC=0;
  if(oldFileName!="")
  {
    cout<<">> Reading old ICs from "<<oldFileName.Data()<<"/"<<oldObjName.Data()<<endl;
    TFile oldFile(oldFileName.Data(),"READ");
    oldIC = (TH2D*) oldFile.Get(oldObjName.Data());
    oldIC->SetDirectory(0);
    oldFile.Close();
  }

  //Load new (incomplete) ICs
  cout<<">> Reading new numerator and denominator from "<<newFileName.Data()<<endl;
  TFile newFile(newFileName.Data(),"UPDATE");
  TH2D* numerator = (TH2D*)newFile.Get("numeratorEB");
  TH2D* denominator = (TH2D*)newFile.Get("denominatorEB");
  TH2D* ICpull = (TH2D*)numerator->Clone("ICpull");//the initial content does not matter 
  TH2D* temporaryIC = (TH2D*)newFile.Get("temporaryICEB");

  //Reset objects that are in general wrongly filled by hadd
  ICpull->Reset();
  temporaryIC->Reset();
  //Create and fill ICpull and newIC
  TH2D* newIC = (TH2D*)temporaryIC->Clone();
  newIC->SetName(newObjName.Data());
  newIC->SetTitle(newObjName.Data());

  ICpull->Divide(numerator,denominator,1.,1.);//ICpull[i]=numerator[i]/denominator[i]

  if(oldFileName!="")
    newIC->Multiply(oldIC,ICpull,1.,1.);
  else
  {
    newIC=(TH2D*)ICpull->Clone();
    newIC->SetName(newObjName.Data());
    newIC->SetTitle(newObjName.Data());
  }

  //Overwrite the newIC file
  cout<<">> Writing new ICs to "<<newFileName.Data()<<"/"<<newIC->GetName()<<endl;
  newFile.cd();
  ICpull->Write(ICpull->GetName(),TObject::kOverwrite);
  temporaryIC->Write(temporaryIC->GetName(),TObject::kOverwrite);
  newIC->Write(newIC->GetName());
  newFile.Close();

  if(oldIC)
    delete oldIC;  
  return 0;
}


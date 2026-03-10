#ifndef FITUTILS__
#define FITUTILS__

#include <iostream>
#include <string>

#include "TFitResult.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"

using namespace std;

namespace FitUtils
{

  inline bool PerseverantFit( TGraph* gr, TF1* fitfunc, string fitopt="RSQME+", int nTrial=10, string TemplatePlotsFolder="")
  {
    //gStyle->SetOptFit(110);
    TFitResultPtr rp;
    int fStatus = 1;
    int iTrial = 0;
    TCanvas c_fit;//canvas to store the result of the template fit if requested by cfg

    while( fStatus!=3 && iTrial<nTrial )
    {
      rp = gr -> Fit(fitfunc, fitopt.c_str());
      fStatus=rp->CovMatrixStatus();
      //cout<<"STATUS"<<fStatus<<endl;                                                                                                                     
      ++iTrial;
      if(fStatus==3)
        cout<<">>>>>>> Converged after "<<iTrial<<" trials"<<endl;
      //cout<<"\n"<<endl;                                                                                                                                    
    }

    if (fitfunc->GetChisquare() / fitfunc->GetNDF() < 1 || fStatus == 2)
      cout<<"[WARNING]: points very close to the line --> matrix not pos def"<<endl;
    else
    {
      if (fitfunc->GetChisquare() / fitfunc->GetNDF() > 1000 || fStatus != 3)
      {
        cout<<"[WARNING]: bad chisquare"<<endl;
	//here I could try to fit changing the minimization method or the minimization steps...
      }
      if(fStatus != 3)
	cout<<">>>>>>> NOT Converged"<<endl;
    }

    if(TemplatePlotsFolder!="")
    {
      gr->Draw("AP");
      c_fit.Print( Form("%s/fit_%s.png",TemplatePlotsFolder.c_str(), gr->GetName()) );
      c_fit.SaveAs( Form("%s/fit_%s.root",TemplatePlotsFolder.c_str(), gr->GetName()) );
    }

    if(fStatus!=2 && fStatus!=3)
      return false;
    return true;
  }
  
  inline bool PerseverantFit( TH1* h, TF1* fitfunc, string fitopt="QRL+", int nTrial=10, string TemplatePlotsFolder="")
  {
    TFitResultPtr rp;
    int fStatus;
    TCanvas c_template_fit;//canvas to store the result of the template fit if requested by cfg
#ifdef DEBUG
    std::cout << "entries: " << h->GetEntries() << std::endl;
#endif
    fitfunc -> SetParameter(1, 0.99);
    for(int nTrial = 0; nTrial < 10; ++nTrial)
    {
      c_template_fit.cd();
      rp = h -> Fit(fitfunc, fitopt.c_str());
      fStatus = rp;
      //std::cout << "nTrial: " << nTrial << "--> fstatus: " << fStatus << std::endl;

      if(fStatus != 4 && fitfunc->GetParError(1) != 0. && fitfunc->GetParameter(1) != 0.99)
      {
        //std::cout << "So I'm HERE goodfit!!!!! " << std::endl;
	if(TemplatePlotsFolder!="")
	{
	  c_template_fit.Print( Form("%s/fit_%s.png",TemplatePlotsFolder.c_str(), h->GetName()) );
	  c_template_fit.SaveAs( Form("%s/fit_%s.root",TemplatePlotsFolder.c_str(), h->GetName()) );
	}
	return true;
      }
    }
    if(TemplatePlotsFolder!="")
    {
      c_template_fit.Print( Form("%s/fit_%s.png",TemplatePlotsFolder.c_str(), h->GetName()) );
      c_template_fit.SaveAs( Form("%s/fit_%s.root",TemplatePlotsFolder.c_str(), h->GetName()) );
    }
    //std::cout << "OH NO --- NOT a goodfit " << std::endl;
    return false;
  }

}

#endif

//////////////////////////////////////////////////////////////////////////////
//
// File to simulate DIS process e+p->e'+jet+X 				  ////
//
/////////////////////////////////////////////////////////////////////////////
#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace Pythia8;
using namespace fastjet;
using namespace std;

// Define User information including event charge
class MyInfo: public PseudoJet::UserInfoBase {
public:
MyInfo(int event_id_in, double event_charge_in, int event_chargeType_in, string event_name_in) :
  _event_id(event_id_in),_event_charge(event_charge_in),_event_chargeType(event_chargeType_in),_event_name(event_name_in){}

int event_id() const {return _event_id;}
double event_charge() const {return _event_charge;}
int event_chargeType() const {return _event_chargeType;}
string event_name() const {return _event_name;}

protected:
int _event_id;
double _event_charge;
int _event_chargeType;
string _event_name;
};
// End user information

double pTjetMin = 10.0, pTjetMax = 25.0, Q2min = 10.0;
double pTeMin = 15.0, pTeMax = 20.0;
double ymin = 0.1, ymax = 0.9;
int   nrange = 10;
double qTMin = 0.0, qTMax = 5;
int nqT = 25;
float radius = 1.; // jet cone radius
float min_eta_jet = -5.0; // minimum of eta for jet
float max_eta_jet =  5.0; // maximum of eta for jet
float min_eta_photon = -5.0; // minimum of eta for photon
float max_eta_photon =  5.0; // maximum of eta for photon
double kTMin = 0.0;
double kTMax = 3.0;
int kTrange = 30;
double zMin = 0.0;
double zMax = 1.0;
double z_min_cut = 0.1; // for kT-distribution, cut on z_min
double z_max_cut = 0.5; // for kT-distribution, cut on z_max
int nz = 20;
double pthMin = 0.25; // cutoff for hadron pT in the jet
int nQjet = 100;
double QjetMin = -3.0;
double QjetMax = 3.0;
ostream & operator<<(ostream &, const PseudoJet &); 

int main() {
 
  // initate root histograms 
  TFile *file = new TFile("e_jet.root", "RECREATE");
  file->cd();

  TH1F *hqT = new TH1F("imbalance qT", Form("qT"), nqT, qTMin, qTMax);
  TH1F *hqTnormal = new TH1F("imbalance qT normalized to one", Form("qT"), nqT, qTMin, qTMax);
  TH1F *hqTG = new TH1F("imbalance qT Gaussian", Form("qT"), nqT, qTMin, qTMax);  

  TH1F *hpe = new TH1F("electron p_{T}^{e}", Form("pTe"), nrange, pTeMin, pTeMax);
  TH1F *hpj = new TH1F("jet p_{T}^{jet}", Form("pTj"), nrange, pTjetMin, pTjetMax);
  TH2F *heJetPT = new TH2F("pt e jet", "p_{T}^{e}p_{T}^{j} (GeV)", nrange, pTeMin, pTeMax, nrange, pTjetMin, pTjetMax);

  TH1F *hadronzh = new TH1F("hadron zh", Form("hadron zh"), nz, zMin, zMax);
  TH1F *hadronkT = new TH1F("hadron kT w.r.t jet", Form("hadron kT"), kTrange, kTMin, kTMax);
  TH1F *hQjet = new TH1F("jet charge", Form("jet charge"), nQjet, QjetMin, QjetMax);
  TH1F *hqTQpos = new TH1F("imbalance qT for positive jet charge", Form("qT positive charge, normalized"), nqT, qTMin, qTMax);
  TH1F *hqTQneu = new TH1F("imbalance qT for neutral jet charge", Form("qT neutral charge, normalized"), nqT, qTMin, qTMax);
  TH1F *hqTQneg = new TH1F("imbalance qT for negtive jet charge", Form("qT negative charge, normalized"), nqT, qTMin, qTMax);
  // end initiation

  // Generator Process selection
  Pythia pythia;
  
  // Shorthand for some public members
  Event &event = pythia.event;
  Settings &settings = pythia.settings;
  //  Info &info = pythia.info;
  ParticleData &particleData = pythia.particleData;
  
  // Set up FastJet jet finder, with anti-kT clustering 
  JetDefinition jet_def(antikt_algorithm, radius);
  vector <PseudoJet> fjInputs; /* Fastjet input */

  // Common settings for all the subruns
  pythia.readFile("main99.cmnd");

  // make charged hadron not decay, so we can observe them in final-state
  // pi0: 111 pi+: 211 K+: 321 p: 2212 n: 2112
  // it also turn off the anti-particle of the above
  int notDecay[4]= {111,211,321,2212};
  for (int iC = 0; iC < 4; ++iC) {
    particleData.mayDecay(notDecay[iC],false);
  }
  
  // Number of events per bin, generated and listed ones
  int nEvent = 1.e6;//1.e7;
  int nList = 0, nListJets = 3;
  int njet = 0;  // count total number of jets

  settings.parm("PhaseSpace:Q2Min",Q2min);
  pythia.init();

  // Begin event loop. Generate event, Skip if error. list first few.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    if (event[5].id() != 1 ) continue;  

    // scattered electron momentum
    Vec4 pProton = event[1].p();
    Vec4 peIn = event[4].p();
    Vec4 peOut = event[6].p();
    Vec4 pPhoton = peIn - peOut;
	
    // Q2, y
    double Q2 = -pPhoton.m2Calc();
    double y = (pProton * pPhoton) / (pProton * peIn);

    double pTel = event[6].pT();
    double pxel = event[6].px();
    double pyel = event[6].py();	 
    double eta_el = event[6].eta();

    // fastjet analysis
    fjInputs.resize(0);
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal() && event[i].id()!=11 && event[i].mother1()!=6) {
	PseudoJet particle( event[i].px(), event[i].py(), 
			    event[i].pz(), event[i].e() );
	// add particle id into fastjet
	particle.set_user_index(abs(event[i].id()));
        // add particle charge         
        particle.set_user_info( new MyInfo(event[i].id(), event[i].charge(),
                                 event[i].chargeType(), event[i].name() ) );
	fjInputs.push_back(particle);
      }
    if(fjInputs.size() == 0) continue;
    vector <PseudoJet> inclusiveJets, sortedJets;
    ClusterSequence clust_seq(fjInputs, jet_def);
    inclusiveJets = clust_seq.inclusive_jets(pTjetMin);
    sortedJets = sorted_by_pt(inclusiveJets);

    // analyze the jet
    for (int ijet = 0; ijet < sortedJets.size(); ijet++) {     
      if (sortedJets[ijet].eta() > max_eta_jet
	  || sortedJets[ijet].eta() < min_eta_jet) continue;
      if (y > ymax || y < ymin) continue;
      if (pTel > pTeMax || pTel < pTeMin) continue;
      njet++;
      
      double pTJet = sortedJets[ijet].pt();
      double pxJet = sortedJets[ijet].px();
      double pyJet = sortedJets[ijet].py();
      double pzJet = sortedJets[ijet].pz();
      double pJet = sqrt(pow(pxJet,2)+pow(pyJet,2)+pow(pzJet,2));
      
      //      cout << sortedJets.size() << '\t' <<  "pTJet:" << pTJet << endl;

      // compute the average transverse momentum of electron+jet
      double px_avg = pxel - pxJet;
      double py_avg = pyel - pyJet;
      double pT_avg = sqrt(pow(px_avg,2.0)+pow(py_avg,2.0))/2.0;
      
      // compute the imbalance qT
      double qx = pxel + pxJet;
      double qy = pyel + pyJet;
      double qT = sqrt(pow(qx,2.0)+pow(qy,2.0));

      if(qT>qTMax) continue;

      hqT->Fill(qT); /* d[sigma]/dq_T */
      hqTnormal->Fill(qT); /* d[sigma]/dq_T */
      hqTG->Fill(qT,1.0/qT); /* d[sigma]/(q_T*dq_T) */
      
      hpe->Fill(pTel);
      hpj->Fill(pTJet);
      heJetPT->Fill(pTel,pTJet);

      // analyze the constituents of a jet
      vector <PseudoJet> constituents = sortedJets[ijet].constituents();  
      for (int iPart = 0; iPart < constituents.size(); iPart++) { 
	// check charged hadrons: 211 pi+/-, 321 K+/-, 2212 proton/anti-proton
	if ( constituents[iPart].user_index() == 211 || constituents[iPart].user_index()==321
	     || constituents[iPart].user_index()==2212 ) {
	  double pth = constituents[iPart].pt();
	  double px_i = constituents[iPart].px();
	  double py_i = constituents[iPart].py();
	  double pz_i = constituents[iPart].pz();

	  double prodmag = pxJet*px_i + pyJet*py_i + pzJet*pz_i;
	  double crossx = pyJet*pz_i - pzJet*py_i;
	  double crossy = pzJet*px_i - pxJet*pz_i;
	  double crossz = pxJet*py_i - pyJet*px_i;
	  double crossmag = sqrt(pow(crossx,2)+pow(crossy,2)+pow(crossz,2));

	  double z_i = prodmag / pow(pJet,2);
	  double kT = crossmag / pJet;
	
	  double psi,zsel,hsel;
	  //		psi = ( (ptJet > pTMin && ptJet < pTMax)? 1 : 0 );
	  zsel = ( (z_i > z_min_cut && z_i < z_max_cut)? 1 : 0 );
	  hsel = ( (pth > pthMin)? 1 : 0 );
	  hadronkT -> Fill(kT, zsel*hsel);
	  hadronzh -> Fill(z_i, hsel);
	}   
      } // end hadron-in-jet loop


        // start the standard jet charge
        // \sum_i pt_i^kappa * Q_i / ptJ^kappa
        double Qjet = 0., kappa = 0.3;       
//        int Cflavor = 0;
        vector<PseudoJet> jevent = sortedJets[ijet].constituents();
        if(jevent.size() == 0) continue;
               
        for (int j = 0; j< jevent.size();++j) {
          if(jevent[j].user_info<MyInfo>().event_chargeType() && abs(jevent[j].user_info<MyInfo>().event_id())==321) {
            Qjet += jevent[j].user_info<MyInfo>().event_charge()*pow(jevent[j].pt(),kappa);            
//            Cflavor = Cflavor + 1;
          }
        }// end jet charge loop
        Qjet = Qjet / pow( sortedJets[ijet].perp(),kappa);

        
//        if(Cflavor==0) continue; // continue if no such flavor is found in the jet;
        hQjet->Fill(Qjet); //Fill in histogram
        double psel,nsel,qzsel; 
        psel = ( (Qjet >=  .25)? 1 : 0 );
        nsel = ( (Qjet <= -.25)? 1 : 0 );
        qzsel = (((Qjet > -.25)&&(Qjet < .25 ))? 1 : 0 );
        hqTQpos->Fill(qT,psel);
        hqTQneu->Fill(qT,qzsel);
        hqTQneg->Fill(qT,nsel);
        // end standard jet charge

     
    } // end jet loop   
  } // end of event loop. statistics
  pythia.stat();
  

  // start root analysis 
  // normalized to cross section
  float sigma = pythia.info.sigmaGen()*1.0e9; // total cross section in pb
  float normone = sigma/nEvent * nqT/(qTMax - qTMin);
  float norm2 = sigma/nEvent * nrange/(pTeMax - pTeMin);
  float norm3 = sigma/nEvent * nrange/(pTjetMax - pTjetMin);
  float normtwo = sigma/nEvent * nrange/(pTeMax - pTeMin) * nrange/(pTjetMax - pTjetMin);
  
  cout << "cross section check: " << sigma << endl;
  cout << "norm: " << normone << endl;
  hqT->Scale(normone);
  hqTG->Scale(normone);

  hpe->Scale(norm2);
  hpj->Scale(norm3);
  heJetPT->Scale(normtwo);
  
  double totcross,test1;
  totcross = hqTnormal->Integral("width");
  test1 = hqTnormal->Integral(0,26,"width");
  cout << "total cross section: " << totcross << " test: " << test1 << endl;
  
  hqTnormal->Scale(1/totcross);

  // hadron-in-jet distribution
  double normkT = njet * (kTMax - kTMin)/kTrange; 
  double normzh = njet * (zMax - zMin)/nz;
  double normQjet = njet * (QjetMax - QjetMin)/nQjet;

  hadronkT->Scale(1.0/normkT);
  hadronzh->Scale(1.0/normzh);
  
//  hQjet->Scale(1.0/normQjet);
  hqTQpos->Scale(1/totcross);
  hqTQneu->Scale(1/totcross);
  hqTQneg->Scale(1/totcross);
 
  double fracp,fracn,fracz;
  double errp,errn,errz;
  fracp = hqTQpos->IntegralAndError(0,nqT,errp,"width");
  fracz = hqTQneu->IntegralAndError(0,nqT,errz,"width");
  fracn = hqTQneg->IntegralAndError(0,nqT,errn,"width");
  cout << "fraction for postive charge: " << fracp  << " +/- " << errp << endl;
  cout << "fraction for zero    charge: " << fracz  << " +/- " << errz << endl;
  cout << "fraction for negtive charge: " << fracn  << " +/- " << errn << endl;

  double Qfracp=0.0,Qfracn=0.0,Qfracz=0.0,Qnorm=0.0,Qbin;
  double Qerrp=0.0,Qerrn=0.0,Qerrz=0.0;
  for (int j=1; j<=nQjet+1 ;++j){
    Qbin = (-3.0+ (j-1)*.06 -3.0 + j*.06 )*.5;
    Qnorm = Qnorm + hQjet->GetBinContent(j);
    if(Qbin >= .25){Qfracp=Qfracp + Qbin*hQjet->GetBinContent(j);Qerrp=Qerrp+pow(Qbin*hQjet->GetBinError(j),2);}
    if(Qbin <=-.25){Qfracn=Qfracn + Qbin*hQjet->GetBinContent(j);Qerrn=Qerrn+pow(Qbin*hQjet->GetBinError(j),2);}
    if(Qbin > -.25 && Qbin < .25){Qfracz=Qfracz + Qbin*hQjet->GetBinContent(j);Qerrz=Qerrz+pow(Qbin*hQjet->GetBinError(j),2);}
  } 
  cout << "Qnorm: " << Qnorm << endl; 
  cout << "postive charge: " << Qfracp/Qnorm << " +/- " << sqrt(abs(Qerrp))/Qnorm << endl;
  cout << "zero    charge: " << Qfracz/Qnorm << " +/- " << sqrt(abs(Qerrz))/Qnorm << endl;
  cout << "negtive charge: " << Qfracn/Qnorm << " +/- " << sqrt(abs(Qerrn))/Qnorm << endl;



  ofstream myfile;  
  myfile.open("charge.dat");
  myfile << "fraction for postive charge: " << fracp  << " +/- " << errp << endl;
  myfile << "fraction for zero    charge: " << fracz  << " +/- " << errz << endl;
  myfile << "fraction for negtive charge: " << fracn  << " +/- " << errn << endl;
  myfile << "postive charge: " << Qfracp/Qnorm << " +/- " << sqrt(Qerrp)/Qnorm << endl;
  myfile << "zero    charge: " << Qfracz/Qnorm << " +/- " << sqrt(Qerrz)/Qnorm << endl;
  myfile << "negtive charge: " << Qfracn/Qnorm << " +/- " << sqrt(Qerrn)/Qnorm << endl;
  myfile << "  charge " << "		" << "normalized events"  << "		" << " errors "  << endl; 
  for (int j=1; j<=nQjet+1 ;++j){
  myfile << -3.0+(j-1)*.06 << "		  " <<  hQjet->GetBinContent(j)/Qnorm << "          " << hQjet->GetBinError(j)/Qnorm  << endl;
  }

//  for (int j=1; j<=10 ;++j){
//    myfile << (hqTQneg->GetBinContent(j))/(hqTnormal->GetBinContent(j)) << "  " << (hqTQpos->GetBinContent(j))/(hqTnormal->GetBinContent(j)) <<endl;
//  }
  myfile << endl;
 
  for (int j=1; j <= nqT + 1; ++j){
    myfile << "{ "  <<qTMin + (j-1)*(qTMax-qTMin)/float(nqT) << ", "  << hqT->GetBinContent(j) << ", " << hqT->GetBinError(j) ;
    myfile << ", "  << normone*totcross*hqTQneg->GetBinContent(j) << ", " << normone*totcross*hqTQneg->GetBinError(j) ;
    myfile << ", "  << normone*totcross*hqTQneu->GetBinContent(j) << ", " << normone*totcross*hqTQneu->GetBinError(j) ;
    myfile << ", "  << normone*totcross*hqTQpos->GetBinContent(j) << ", " << normone*totcross*hqTQpos->GetBinError(j) << " }," << endl ;
  }

  myfile.close();

  file->Write();
  // end analysis 

  return 0;
}

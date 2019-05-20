#ifndef HighPtMuonAnalysis_GEScale
#define HighPtMuonAnalysis_GEScale

#include "GEScaleSyst.h"

using namespace std;

float GEScaleSyst::GetGEScaleKappa(int icopy, float eta, float phi)
{
  int the_eta = -999;
  int the_phi = -999;

  for(int ieta=0; ieta<neta; ++ieta) {
    if( eta >= etaBinEdge[ieta] && eta < etaBinEdge[ieta+1] ) {
      the_eta = ieta;
      break;
    }
  }

  for(int iphi=0; iphi<nphi; ++iphi) {
    if( phi >= phiBinEdge[iphi] && phi < phiBinEdge[iphi+1] ) {
      the_phi = iphi;
      break;
    }
  }

  if(the_eta < 0 || the_phi < 0) {
    cout << "ERROR   GEScaleSyst::GetGEScaleKappa   correction value not found... icopy= " << icopy << ", eta= " << eta << ", phi= " << phi << endl;
    return -1e9;
  }

  return _matrix[icopy][the_eta][the_phi];
}

float GEScaleSyst::GEScaleCorrPt(int icopy, float pt, float eta, float phi, int charge)
{
  if( _matrix.find(icopy) == _matrix.end() ) {
    cout << "ERROR   GEScaleSyst::GEScaleCorrPt   invalid icopy= " << icopy << endl;
    return -1e9;
  }
  if( pt < 0 ) {
    cout << "ERROR   GEScaleSyst::GEScaleCorrPt   invalid muon pt= " << pt << endl;
    return -1e9;
  }
  if( fabs(eta) > 2.4 ) {
    cout << "ERROR   GEScaleSyst::GEScaleCorrPt   invalid muon eta= " << eta << endl;
    return -1e9;
  }
  if( fabs(phi) > M_PI ) {
    cout << "ERROR   GEScaleSyst::GEScaleCorrPt   invalid muon phi= " << phi << endl;
    return -1e9;
  }
  if( fabs(charge) != 1 ) {
    cout << "ERROR   GEScaleSyst::GEScaleCorrPt   invalid muon charge= " << charge << endl;
    return -1e9;
  }

  float kappa = GetGEScaleKappa(icopy, eta, phi);

  if(verbose==1) {
    cout << "\n Start GE pT correction:" << endl;
    cout << "    initial pT= " << pt << ", eta= " << eta << ", phi= " << phi << ", q= " << charge << endl;
    cout << "    matrix for icopy= " << icopy << endl;
    for(int ieta=0; ieta<neta; ++ieta) {
      for(int iphi=0; iphi<nphi; ++iphi) {
        printf("%3.5f    ", _matrix[icopy][ieta][iphi]);
      }
      cout << endl;
    }
    cout << "\t --> kappa= " << kappa << endl;
  }

  float corr_pt = pt;
  corr_pt = corr_pt / 1000.;              // convert GeV to TeV
  if(verbose==1)  cout << "\t" << corr_pt << endl;
  corr_pt = charge * fabs(corr_pt);       // convert to signed pT
  if(verbose==1)  cout << "\t" << corr_pt << endl;
  corr_pt = 1./corr_pt;                   // convert to q/pT curvature
  if(verbose==1)  cout << "\t" << corr_pt << endl;
  corr_pt = corr_pt + kappa;              // apply correction
  if(verbose==1)  cout << "\t" << corr_pt << endl;
  corr_pt = 1./corr_pt;                   // revert to pT
  if(verbose==1)  cout << "\t" << corr_pt << endl;
  corr_pt = fabs(corr_pt);                // revert to pT
  if(verbose==1)  cout << "\t" << corr_pt << endl;
  corr_pt = corr_pt * 1000.;              // revert to GeV
  if(verbose==1)  cout << "\t" << corr_pt << endl;
  if(verbose==1)  cout << "\t --> diff= " << (corr_pt-pt) << endl;

  return corr_pt;
}

TLorentzVector GEScaleSyst::GEScaleCorrLvec(int icopy, float pt, float eta, float phi, int charge)
{
  TLorentzVector Lvec;
  float corr_pt = GEScaleCorrPt(icopy, pt, eta, phi, charge);

  if(corr_pt < 0) {
    cout << "ERROR   GEScaleSyst::GEScaleCorrLvec   invalid corrected pt= " << corr_pt << endl;
    Lvec.SetPtEtaPhiM(0, 0, 0, 0);
    return Lvec;
  }

  Lvec.SetPtEtaPhiM(corr_pt, eta, phi, mu_mass);

  if(verbose==1)  cout << "Final muon pT, eta, phi" << endl;
  if(verbose==1)  cout << "pT=  " << Lvec.Pt() << endl;
  if(verbose==1)  cout << "eta= " << Lvec.Eta() << endl;
  if(verbose==1)  cout << "phi= " << Lvec.Phi() << endl;
  if(verbose==1)  cout << "E=   " << Lvec.E() << endl;
  if(verbose==1)  cout << "m=   " << Lvec.M() << endl;

  return Lvec;
}

void GEScaleSyst::AddNewMatrix( int new_copy, std::map<int, std::map<int, float> > new_matrix ) {
  if(new_copy < 100) {
    cout << "ERROR   GEScaleSyst::AddNewMatrix   new_copy < 100 --> return" << endl;
    return;
  }

  for(int ieta=0; ieta<neta; ++ieta) {
    for(int iphi=0; iphi<nphi; ++iphi) {
      _matrix[new_copy][ieta][iphi] = new_matrix[ieta][iphi];
    }
  }
}


GEScaleSyst::GEScaleSyst()
{
  // icopy == 1, just for test
  _matrix[-1][0][0] =  0.2;  _matrix[-1][0][1] = 0.01;  _matrix[-1][0][2] = -0.2;
  _matrix[-1][1][0] =  0.2;  _matrix[-1][1][1] = 0.01;  _matrix[-1][1][2] = -0.2;
  _matrix[-1][2][0] =  0.2;  _matrix[-1][2][1] = 0.01;  _matrix[-1][2][2] = -0.2;
  _matrix[-1][3][0] = -0.2;  _matrix[-1][3][1] = 0.02;  _matrix[-1][3][2] =  0.2;
  _matrix[-1][4][0] = -0.2;  _matrix[-1][4][1] = 0.02;  _matrix[-1][4][2] =  0.2;
  _matrix[-1][5][0] = -0.2;  _matrix[-1][5][1] = 0.02;  _matrix[-1][5][2] =  0.2;

  
}

GEScaleSyst::~GEScaleSyst() {}

#endif

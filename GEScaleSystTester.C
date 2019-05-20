#include "GEScaleSyst.h"
#include "GEScaleSyst.cc"

void GEScaleSystTester( int icopy = -1, float pt = 1000., float eta = 0., float phi = 0., int q = 1 )
{
  GEScaleSyst *GE = new GEScaleSyst();
  GE->SetVerbose(1);

  TLorentzVector Lvec = GE->GEScaleCorrLvec( icopy, pt, eta, phi, q );

  delete GE;
}

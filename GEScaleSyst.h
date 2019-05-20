#ifndef HighPtMuonAnalysis_GEScale_H
#define HighPtMuonAnalysis_GEScale_H
// #pragma once
#include <iostream>
#include <cmath>

#include "TLorentzVector.h"

class GEScaleSyst {
  public:
    GEScaleSyst();
    virtual ~GEScaleSyst();

    float GEScaleCorrPt(int icopy, float pt, float eta, float phi, int charge);
    TLorentzVector GEScaleCorrLvec(int icopy, float pt, float eta, float phi, int charge);

    void SetVerbose(int _v) {
      verbose = _v;
    };

    // new_copy >= 100 for user defined matrices
    void AddNewMatrix( int new_copy, std::map<int, std::map<int, float> > new_matrix );

  private:
    const double mu_mass = 0.105658;

    float GetGEScaleKappa(int icopy, float eta, float phi);

    const int neta = 6;
    const double etaBinEdge[7] = { -2.4, -2.1, -1.2, 0, 1.2, 2.1, 2.4 };
    const int nphi = 3;
    const double phiBinEdge[4] = { -M_PI, -M_PI/3.0, M_PI/3.0, M_PI };
    // _matrix[copy #][eta bin #][phi bin #]
    // copy #:
    //   0 (original),
    //   1, 2, 3, ... (variation)
    // eta bin #:
    //   0: [-2.4, -2.1]
    //   1: [-2.1, -1.2]
    //   2: [-1.2, 0.0]
    //   3: [0.0, 1.2]
    //   4: [1.2, 2.1]
    //   5: [2.1, 2.4]
    // phi bin #:
    //   0: [-180, -60]
    //   1: [-60, 60]
    //   2: [60, 180]
    std::map<int, std::map<int, std::map<int, float> > > _matrix;

    int verbose;
};

#endif

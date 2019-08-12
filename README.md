# GEScaleSyst

``` C
float GEScaleCorrPt(int icopy, float pt, float eta, float phi, int charge, bool doOpp = false)
  - returns corrected pT

TLorentzVector GEScaleCorrLvec(int icopy, float pt, float eta, float phi, int charge, bool doOpp = false)
  - returns TLorentzVector with corrected pT
```
``` txt
icopy:
  - yr: 16, 17 (77 for 94X samples with correct tracker alignment), 18
  - yr0000, yr0001, ...: random copies using gaus(kappa, sigma)
  - yr1000, yr1001, ...: random copies passing sign contraint, restrict sign of bias if (k + 1 sigma)(k - 1 sigma) > 0
```

## Simple test
root -l 'GEScaleSystTester.C( icopy, pT, eta, phi, charge)'

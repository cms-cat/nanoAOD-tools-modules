/*  
 * Run3 Muon correction module adapted after the example in https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/.
 * See src/MuonScaleRe.cc file for techincal notes.
 */

#include <correction.h>
#include <TRandom3.h>
#include <string>

class MuonScaRe {
public:
  MuonScaRe(std::string json, double low_pt_threshold = 26);

  double pt_resol(double pt, double eta, double phi, float nL, uint64_t evtNumber, int lumiNumber);
  
  double pt_scale(bool is_data, double pt, double eta, double phi, int charge);

  double pt_resol_var(double pt_woresol, double pt_wresol, double eta, std::string updn);

  double pt_scale_var(double pt, double eta, double phi, int charge, std::string updn);

private:
  double get_k(double eta, std::string var);
  double get_std(double pt, double eta, float nL);
  double get_rndm(double eta, double phi, float nL, uint64_t evtNumber, int lumiNumber);

  std::unique_ptr<correction::CorrectionSet> cset;
  double low_pt_threshold;
  TRandom3 rnd;
};

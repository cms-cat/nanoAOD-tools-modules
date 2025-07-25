/*  
 * Run3 Muon correction module, imported from:
 * https://gitlab.cern.ch/cms-muonPOG/muonscarekit/-/blob/7e7dd87b26964cc16bea0cf9bc2a849bd7af870e/scripts/MuonScaRe.cc
 * with the following modifications:
 * -encapsulated in a class so that the correctionSet and RNG are kept 
 *  internally
 * -add function to optionally set the RNG seed, so that the user can implement 
 *  deterministic pseudo-random smearing. If desired, a proper seed should 
 *  be passed once for each muon before calling pt_resol the first time for 
 *  that mu.
 *  Note: setting the seed is left to the user since in the current 
 *  implementation get_rndm, which normally getsis called 3 times for each mu 
 *  (for var = norm, syst, stat), does not have access to enough sources of 
 *  entropy to create a decent seed internally.
 * -avoid trivial string manipulations in pt_scale
 */

#include <boost/math/special_functions/erf.hpp>
#include <string>
#include <PhysicsTools/NATModules/interface/MuonScaRe.h>

using namespace std;

MuonScaRe::MuonScaRe(string json) : cset(correction::CorrectionSet::from_file(json)){}
 

struct CrystalBall{
    double pi=3.14159;
    double sqrtPiOver2=sqrt(pi/2.0);
    double sqrt2=sqrt(2.0);
    double m;
    double s;
    double a;
    double n;
    double B;
    double C;
    double D;
    double N;
    double NA;
    double Ns;
    double NC;
    double F;
    double G;
    double k;
    double cdfMa;
    double cdfPa;
CrystalBall():m(0),s(1),a(10),n(10){
    init();
}
CrystalBall(double mean, double sigma, double alpha, double n)
    :m(mean),s(sigma),a(alpha),n(n){
    init();
}
void init(){
    double fa = fabs(a);
    double ex = exp(-fa*fa/2);
    double A  = pow(n/fa, n) * ex;
    double C1 = n/fa/(n-1) * ex; 
    double D1 = 2 * sqrtPiOver2 * erf(fa/sqrt2);
    B = n/fa-fa;
    C = (D1+2*C1)/C1;   
    D = (D1+2*C1)/2;   
    N = 1.0/s/(D1+2*C1); 
    k = 1.0/(n-1);  
    NA = N*A;       
    Ns = N*s;       
    NC = Ns*C1;     
    F = 1-fa*fa/n; 
    G = s*n/fa;    
    cdfMa = cdf(m-a*s);
    cdfPa = cdf(m+a*s);
}
double pdf(double x) const{ 
    double d=(x-m)/s;
    if(d<-a) return NA*pow(B-d, -n);
    if(d>a) return NA*pow(B+d, -n);
    return N*exp(-d*d/2);
}
double pdf(double x, double ks, double dm) const{ 
    double d=(x-m-dm)/(s*ks);
    if(d<-a) return NA/ks*pow(B-d, -n);
    if(d>a) return NA/ks*pow(B+d, -n);
    return N/ks*exp(-d*d/2);

}
double cdf(double x) const{
    double d = (x-m)/s;
    if(d<-a) return NC / pow(F-s*d/G, n-1);
    if(d>a) return NC * (C - pow(F+s*d/G, 1-n) );
    return Ns * (D - sqrtPiOver2 * erf(-d/sqrt2));
}
double invcdf(double u) const{
    if(u<cdfMa) return m + G*(F - pow(NC/u, k));
    if(u>cdfPa) return m - G*(F - pow(C-u/NC, -k) );
    return m - sqrt2 * s * boost::math::erf_inv((D - u/Ns )/sqrtPiOver2);
}
};

double MuonScaRe::get_rndm(double eta, float nL) {

    // obtain parameters from correctionlib
    double mean = cset->at("cb_params")->evaluate({abs(eta), nL, 0});
    double sigma = cset->at("cb_params")->evaluate({abs(eta), nL, 1});
    double n = cset->at("cb_params")->evaluate({abs(eta), nL, 2});
    double alpha = cset->at("cb_params")->evaluate({abs(eta), nL, 3});
    
    // instantiate CB and get random number following the CB
    CrystalBall cb(mean, sigma, alpha, n);
    /** Original implementation using time() as a seed; replaced by our own RNG
    TRandom3 rnd(time(0));
    double rndm = gRandom->Rndm();
    return cb.invcdf(rndm);
    **/
    return cb.invcdf(rng.Rndm());
}


double MuonScaRe::get_std(double pt, double eta, float nL) {

    // obtain paramters from correctionlib
    double param_0 = cset->at("poly_params")->evaluate({abs(eta), nL, 0});
    double param_1 = cset->at("poly_params")->evaluate({abs(eta), nL, 1});
    double param_2 = cset->at("poly_params")->evaluate({abs(eta), nL, 2});

    // calculate value and return max(0, val)
    double sigma = param_0 + param_1 * pt + param_2 * pt*pt;
    if (sigma < 0) sigma = 0;
    return sigma; 
}


double MuonScaRe::get_k(double eta, string var) {

    // obtain parameters from correctionlib
    double k_data = cset->at("k_data")->evaluate({abs(eta), var});
    double k_mc = cset->at("k_mc")->evaluate({abs(eta), var});

    // calculate residual smearing factor
    // return 0 if smearing in MC already larger than in data
    double k = 0;
    if (k_mc < k_data) k = sqrt(k_data*k_data - k_mc*k_mc);
    return k;
}


double MuonScaRe::pt_resol(double pt, double eta, float nL) {

    // load correction values
    double rndm = (double) get_rndm(eta, nL);
    double std = (double) get_std(pt, eta, nL);
    double k = (double) get_k(eta, "nom");

    // calculate corrected value and return original value if a parameter is nan
    double ptc = pt * ( 1 + k * std * rndm);
    if (isnan(ptc)) ptc = pt;
    if(ptc / pt > 2 || ptc / pt < 0.1 || ptc < 0) {
	ptc = pt;
    }	
    return ptc;
}

double MuonScaRe::pt_resol_var(double pt_woresol, double pt_wresol, double eta, string updn){

    double k = (double) get_k(eta, "nom");

    if (k==0) return pt_wresol;

    double k_unc = cset->at("k_mc")->evaluate({abs(eta), "stat"});

    double std_x_rndm = (pt_wresol / pt_woresol - 1) / k;

    double pt_var = pt_wresol;

    if (updn=="up"){
        pt_var = pt_woresol * (1 + (k+k_unc) * std_x_rndm);
    }
    else if (updn=="dn"){
        pt_var = pt_woresol * (1 + (k-k_unc) * std_x_rndm);
    }
    else {
        cout << "ERROR: updn must be 'up' or 'dn'" << endl;
    }
    if(pt_var / pt_woresol > 2 || pt_var / pt_woresol < 0.1 || pt_var < 0) {
      pt_var = pt_woresol;
    }

    return pt_var;
}

double MuonScaRe::pt_scale(bool is_data, double pt, double eta, double phi, int charge) {
  
    double a = cset->at(is_data?"a_data":"a_mc")->evaluate({eta, phi, "nom"});
    double m = cset->at(is_data?"m_data":"m_mc")->evaluate({eta, phi, "nom"});
    return 1. / (m/pt + charge * a);
}


double MuonScaRe::pt_scale_var(double pt, double eta, double phi, int charge, string updn) {

    double stat_a = cset->at("a_mc")->evaluate({eta, phi, "stat"});
    double stat_m = cset->at("m_mc")->evaluate({eta, phi, "stat"});
    double stat_rho = cset->at("m_mc")->evaluate({eta, phi, "rho_stat"});

    double unc = pt*pt*sqrt(stat_m*stat_m / (pt*pt) + stat_a*stat_a + 2*charge*stat_rho*stat_m/pt*stat_a);

    double pt_var = pt;
    
    if (updn=="up"){
        pt_var = pt + unc;
    }
    else if (updn=="dn"){
        pt_var = pt - unc;
    }

    return pt_var;
}

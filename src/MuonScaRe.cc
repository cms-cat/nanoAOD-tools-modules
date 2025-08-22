/*  
 * Muon correction module, imported from:
 * https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/773fd624aa67c3f742c85530e1a3979bce78b179/examples/MuonScaRe.cc
 * with the following modifications:
 * -encapsulated in a class so that the correctionSet and RNG are kept 
 *  internally
 * -avoid trivial string manipulations in pt_scale
 */

#include <PhysicsTools/NATModules/interface/MuonScaRe.h>
#include <boost/math/special_functions/erf.hpp>
//#include <cstdint>
//#include <cmath>
//#include <algorithm>

using namespace std;

MuonScaRe::MuonScaRe(string json, double threshold) : cset(correction::CorrectionSet::from_file(json)),
						      low_pt_threshold(threshold){}

class SeedSequence {
public:
    explicit SeedSequence(std::initializer_list<uint32_t> seeds)
        : m_seeds(seeds) {}

    template <typename Iter>
    void generate(Iter begin, Iter end) const {
        const size_t n = std::distance(begin, end);
	if (n == 0) return;

	const uint32_t mult = 0x9e3779b9;
	const uint32_t mix_const = 0x85ebca6b;

	std::vector<uint32_t> buffer(n, 0x8b8b8b8b);

	size_t s = m_seeds.size();

        size_t i = 0;

	for(; i < std::min(n, s); ++i) {
            buffer[i] = buffer[i] ^ (m_seeds[i] + mult * i);
	}
	for(; i < n; ++i) {
            buffer[i] = buffer[i] ^ (mult * i);
        }

	for (size_t k = 0; k < n; ++k) {
            uint32_t z = buffer[(k + n - 1) % n] ^ (buffer[k] >> 27);
	    buffer[k] = (z * mix_const) ^ (buffer[k] << 13);
        }

	std::copy(buffer.begin(), buffer.end(), begin);
    }

private:
    std::vector<uint32_t> m_seeds;
};

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

double MuonScaRe::get_rndm(double eta, double phi, float nL, int evtNumber, int lumiNumber) {
    // obtain parameters from correctionlib
    double mean = cset->at("cb_params")->evaluate({abs(eta), nL, 0});
    double sigma = cset->at("cb_params")->evaluate({abs(eta), nL, 1});
    double n = cset->at("cb_params")->evaluate({abs(eta), nL, 2});
    double alpha = cset->at("cb_params")->evaluate({abs(eta), nL, 3});
   
    // instantiate CB and get random number following the CB
    CrystalBall cb(mean, sigma, alpha, n);
    int64_t phi_seed = static_cast<int64_t>((phi / M_PI) * ((1LL << 31) - 1)) & 0xFFF;
    SeedSequence seq{static_cast<uint32_t>(evtNumber), static_cast<uint32_t>(lumiNumber), static_cast<uint32_t>(phi_seed)};
    uint32_t seed;
    seq.generate(&seed, &seed + 1);

    rnd.SetSeed(seed);
    return cb.invcdf(rnd.Rndm());
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


double MuonScaRe::pt_resol(double pt, double eta, double phi, float nL, int evtNumber, int lumiNumber) {
    // load correction values
    double rndm = (double) get_rndm(eta, phi, nL, evtNumber, lumiNumber);
    double std = (double) get_std(pt, eta, nL);
    double k = (double) get_k(eta, "nom");

    // calculate corrected value and return original value if a parameter is nan
    double ptc = pt * ( 1 + k * std * rndm);
    if (isnan(ptc)) ptc = pt;
    if(ptc / pt > 2 || ptc / pt < 0.1 || ptc < 0 || pt < low_pt_threshold || pt > 200){
	ptc = pt;
    }
    // TODO: Understand why for evts with pT < threshold the pt_corr is set to one
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
    if(pt_var / pt_woresol > 2 || pt_var / pt_woresol < 0.1 || pt_var < 0){
        pt_var = pt_woresol; 
    }

    return pt_var;
}

double MuonScaRe::pt_scale(bool is_data, double pt, double eta, double phi, int charge) {
        
    // use right correction
    double a = cset->at(is_data?"a_data":"a_mc")->evaluate({eta, phi, "nom"});
    double m = cset->at(is_data?"m_data":"m_mc")->evaluate({eta, phi, "nom"});
    if(pt < low_pt_threshold)
	    return pt;

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

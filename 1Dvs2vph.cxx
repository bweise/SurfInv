/* 1Dvs2vph
 * Berechnet aus einem 1D Vs,Vp,Dichte-Modell die Phasengeschwindigkeit
 * zu gegebenen Perioden.
 * Quellen: Haskell (1953), Dunkin (1965), Wathelet (2005),
 * Cercato (2007)
 */

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <fstream>
#include <netcdf>
#include <boost/math/tools/roots.hpp>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

typedef complex<double> dcomp;
const std::complex<double> i(0, 1.0);

bool verbose = 0; // set to 1 for more output
double tolerance = 0.0001; // Tolerance for phase velocity [m/s]
double mode_skip_it = 10.0;	// Number of additional iterations to check for mode skipping & factor to increase precision

// function of root of rayleigh velocity
double compute_fvr(double vp, double vs, double vr){
	double fvr = 4.0-4.0*(pow(vr,2)/pow(vs,2))+pow(vr,4)/pow(vs,4)-4.0*sqrt(1-pow(vr,2)/pow(vp,2))*sqrt(1.0-pow(vr,2)/pow(vs,2));
	return fvr;
}

// derivative of rayleigh velocity function
double compute_dfvr(double vp, double vs, double vr){
	double dfvr = -8.0*vr/pow(vs,2)+4.0*pow(vr,3)/pow(vs,4)+(4.0*vr*(pow(vp,2)+pow(vs,2)-2.0*pow(vr,2)))/(vp*vs*sqrt((vp-vr)*(vp+vr))*sqrt((vs-vr)*(vs+vr)));
	return dfvr;
}

// computation of Rayleigh velocity of homogenous half space using Newton Raphson method
double newton_vr(double vp, double vs){
	// variables to store rayleigh velocity
	double vrlast = vs-1;
	double vr;
	double diff=99999.0;
	// calculate Rayleigh velocity for homogenous half space
	while(diff>tolerance){ 
		double fvr = compute_fvr(vp, vs, vrlast);
		double dfvr = compute_dfvr(vp, vs, vrlast);
		vr = vrlast - fvr/dfvr;
		diff = sqrt(pow(vr-vrlast,2));
		vrlast = vr;
	}
	if (verbose == 1)
		cout << "vr: " << vr << "\t final diff: " << diff << "\n\n";
	return vr;
}

std::tuple<dcomp,dcomp,double,double,double,double,double,dcomp,dcomp,double> compute_util(double w, double c, double vp, double vs, double dn, bool botlay){
	// computes some constants for each layer (vertical wave numbers and some derived properties)
	dcomp hv, kv;
	double SH, CH, SK, CK, mh, mk, k = w/c;
	if (c < vp){
		mh = sqrt(pow(w/c,2)-pow(w/vp,2));
		hv = mh;
		if (botlay == 0){
			SH = (k/mh)*sinh(mh*dn);
			CH = cosh(mh*dn);
		}
	}
	else{
		mh = sqrt(pow(w/vp,2)-pow(w/c,2));
		hv = mh;
		if (botlay == 0){
			hv = hv*i;
			SH = (k/mh)*sin(mh*dn);
			CH = cos(mh*dn);
		}
	}
	if (c < vs){
		mk = sqrt(pow(w/c,2)-pow(w/vs,2));
		kv = mk;
		if (botlay == 0){
			SK = (k/mk)*sinh(mk*dn);
			CK = cosh(mk*dn);
		}
	}
	else{
		mk = sqrt(pow(w/vs,2)-pow(w/c,2));
		kv = mk;
		if (botlay == 0){
			kv = kv*i;
			SK = (k/mk)*sin(mk*dn);
			CK = cos(mk*dn);
		}
	}
	double gam = 2.0*pow(vs,2)/pow(c,2);
	dcomp hvnorm = hv/k;
	dcomp kvnorm = kv/k;
	double l = 2.0*pow(k,2)-pow(w/vs,2);
	
	if (verbose==1){
		cout << "hv: " << hv << " kv: " << kv << "\n"
			<< "hvnorm: " << hvnorm << " kvnorm: " << kvnorm << "\n"
			<< "SH: " << SH << " SK: " << SK << " CH: " << CH << " CK: " << CK << "\n"
			<< "gamma: " << gam << " l: " << l << "\n\n";
	}
	
	return std::make_tuple(hv, kv, SH, CH, SK, CK, gam, hvnorm, kvnorm, l);
}

std::tuple<double,double,double,double,double> compute_T(double w, double c, double vp, double vs, double mu){
	// computes layer matrix for bottom layer (an i denotes an imaginary subdeterminant e.g. iT1214)
	double k = w/c;
	auto util = compute_util(w, c, vp, vs, 99999, 1);
	dcomp hv = std::get<0>(util);
	dcomp kv = std::get<1>(util);
	double l = std::get<9>(util);
	
	dcomp fact = pow((-1.0)*pow(vs,2)/(2*mu*hv*kv*pow(w,2)),2);

	double T1212 = std::real(pow(mu,2)*kv*hv*(pow(l,2)-4.0*pow(k,2)*kv*hv)*fact);
	double T1213 = std::real(mu*pow(hv,2)*kv*(l-2.0*pow(k,2))*fact);
	double iT1214 = std::real(k*mu*hv*kv*(l-2.0*hv*kv)*fact);
	double T1224 = std::real(mu*hv*pow(kv,2)*(2.0*pow(k,2)-l)*fact);
	double T1234 = std::real(hv*kv*(pow(k,2)-hv*kv)*fact);
	
	if (verbose==1) {
		cout << "factor: " << fact << "\n"
			<< "T1212: " << T1212 << "\t"
			<< "T1213: " << T1213 << "\t"
			<< "iT1214: " << iT1214 << "\t"
			<< "T1224: " << T1224 << "\t"
			<< "T1234: " << T1234 << "\n\n";
	}
	return std::make_tuple(T1212,T1213,iT1214,T1224,T1234);
}

std::tuple<double,double,double,double,double,double,double,double,double,double,double,double,double,double,double> compute_G(double c, double dn, double w, double vp, double vs, double dens){
	// computes subdeterminants of G matrix (an i denotes an imaginary subdeterminant e.g. iG1214)
	auto kvert = compute_util(w, c, vp, vs, dn, 0);
	double SH = std::get<2>(kvert);
	double CH = std::get<3>(kvert);
	double SK = std::get<4>(kvert);
	double CK = std::get<5>(kvert);
	double gam = std::get<6>(kvert);
	dcomp hvnorm = std::get<7>(kvert);
	dcomp kvnorm = std::get<8>(kvert);
	
	double G1212 = std::real(2.0*gam * (1.0-gam) + (2.0*pow(gam,2)-2.0*gam+1.0) * CH*CK - (pow(1.0-gam,2) + pow(gam,2)*pow(hvnorm,2)*pow(kvnorm,2))*SH*SK);
	double G1213 = std::real((1.0/(dens*w*c))*(CH*SK - SH*CK*pow(hvnorm,2)));
	double iG1214 = std::real((1.0/(dens*w*c))*((1.0 - 2.0*gam)*(1.0-CK*CH) + (1.0 - gam - gam*pow(hvnorm,2)*pow(kvnorm,2))*SH*SK));
	double G1224 = std::real((1.0/(dens*w*c))*(pow(kvnorm,2)*CH*SK - SH*CK));
	double G1234 = std::real((-1.0/(pow(dens,2)*pow(w,2)*pow(c,2)))*(2.0*(1.0 - CH*CK) + (1.0 + pow(kvnorm,2)*pow(hvnorm,2))*SK*SH));
	double G1312 = std::real(dens*w*c*(pow(gam,2)*pow(kvnorm,2)*CH*SK - pow(1.0-gam,2)*SH*CK));
	double G1313 = std::real(CH*CK);
	double iG1314 = std::real((1.0 - gam)*SH*CK + gam*pow(kvnorm,2)*CH*SK);
	double G1324 = std::real((-1.0)*pow(kvnorm,2)*SH*SK);
	double iG1412 = std::real(dens*w*c*((3.0*pow(gam,2) - 2.0*pow(gam,3) - gam)*(1.0 - CH*CK)+(pow(1.0 - gam,3) - pow(gam,3)*pow(hvnorm,2)*pow(kvnorm,2))*SH*SK));
	double iG1413 = std::real((-1.0)*((1.0 - gam)*CH*SK + gam*pow(hvnorm,2)*SH*CK));
	double G1414 = std::real(1.0 - 2.0*gam*(1.0 - gam)*(1.0 - CH*CK) + (pow(1.0 - gam,2) + pow(gam,2)*pow(kvnorm,2)*pow(hvnorm,2))*SH*SK);
	double G2412 = std::real(dens*w*c*(pow(1.0 - gam,2)*CH*SK - pow(gam,2)*SH*CK*pow(hvnorm,2)));
	double G2413 = std::real((-1.0)*pow(hvnorm,2)*SH*SK);
	double G3412 = std::real((-1.0)*pow(dens,2)*pow(w,2)*pow(c,2)*(2.0*pow(gam,2)*pow(1.0 - gam,2)*(1.0 - CH*CK) + (pow(1.0 - gam,4)+pow(gam,4)*pow(hvnorm,2)*pow(kvnorm,2))*SH*SK));
	
	if (verbose==1){
		cout << "G-Komponenten:\n"
			<< "G1212: " << G1212 << " G1213: " << G1213 << " iG1214: " << iG1214 << " G1224: " << G1224 << " G1234: " << G1234 << "\n"
			<< "G1312: " << G1312 << " G1313: " << G1313 << " iG1314: " << iG1314 << " G1324: " << G1324 << "\n"
			<< "iG1412: " << iG1412 << " iG1413: " << iG1413 << " G1414: " << G1414 << "\n"
			<< "G2412: " << G2412 << " G2413: " << G2413 << "\n"
			<< "G3412: " << G3412 << "\n\n";
	}
		
	return std::make_tuple(G1212,G1213,iG1214,G1224,G1234,G1312,G1313,iG1314,G1324,iG1412,iG1413,G1414,G2412,G2413,G3412);
}

std::tuple<double,double,double,double,double> compute_R(double w, double c, double vp, double vs, double dn, double dens, std::tuple<double,double,double,double,double> T){
	// Recursive layer stacking from bottom to top layer (an i denotes an imaginary subdeterminant e.g. iR1214)
	double T1212 = std::get<0>(T);
	double T1213 = std::get<1>(T);
	double iT1214 = std::get<2>(T);
	double T1224 = std::get<3>(T);
	double T1234 = std::get<4>(T);
	auto G = compute_G(c,dn,w,vp,vs,dens);
	double G1212 = std::get<0>(G);
	double G1213 = std::get<1>(G);
	double iG1214 = std::get<2>(G);
	double G1224 = std::get<3>(G);
	double G1234 = std::get<4>(G);
	double G1312 = std::get<5>(G);
	double G1313 = std::get<6>(G);
	double iG1314 = std::get<7>(G);
	double G1324 = std::get<8>(G);
	double iG1412 = std::get<9>(G);
	double iG1413 = std::get<10>(G);
	double G1414 = std::get<11>(G);
	double G2412 = std::get<12>(G);
	double G2413 = std::get<13>(G);
	double G3412 = std::get<14>(G);
	
	double R1212 = T1212*G1212 + T1213*G1312 - 2.0*iT1214*iG1412 + T1224*G2412 + T1234*G3412;
	double R1213 = T1212*G1213 + T1213*G1313 - 2.0*iT1214*iG1413 + T1224*G2413 + T1234*G2412;
	double iR1214 = T1212*iG1214 + T1213*iG1314 + iT1214*(2.0*G1414-1.0) + T1224*iG1413 + T1234*iG1412;
	double R1224 = T1212*G1224 + T1213*G1324 - 2.0*iT1214*iG1314 + T1224*G1313 + T1234*G1312;
	double R1234 = T1212*G1234 + T1213*G1224 - 2.0*iT1214*iG1214 + T1224*G1213 + T1234*G1212;
	
	// Normalize R matrix components to +-100000
	double tmpmax, maxR = R1212;
	if(maxR<0)
		maxR = (-1.0)*maxR;
	tmpmax = R1213;
	if (tmpmax<0)
		tmpmax = (-1.0)*tmpmax;
	if (tmpmax>maxR)
		maxR = tmpmax;
	tmpmax = iR1214;
	if (tmpmax<0)
		tmpmax = (-1.0)*tmpmax;
	if (tmpmax>maxR)
		maxR = tmpmax;
	tmpmax = R1224;
	if (tmpmax<0)
		tmpmax = (-1.0)*tmpmax;
	if (tmpmax>maxR)
		maxR = tmpmax;
	tmpmax = R1234;
	if (tmpmax<0)
		tmpmax = (-1.0)*tmpmax;
	if (tmpmax>maxR)
		maxR = tmpmax;
	if (maxR>1.0e5){
		maxR = 1.0e5/maxR;
		R1212 = maxR*R1212;
		R1213 = maxR*R1213;
		iR1214 = maxR*iR1214;
		R1224 = maxR*R1224;
		R1234 = maxR*R1234;
	}
	
	if (verbose==1) {
		cout << "R1212: " << R1212 << "\n"
		<< "R1213: " << R1213 << "\n"
		<< "iR1214: " << iR1214 << "\n"
		<< "R1224: " << R1224 << "\n"
		<< "R1234: " << R1234 << "\n\n";
	}
	
	return std::make_tuple(R1212,R1213,iR1214,R1224,R1234);
}

double compute_R1212(double w, double c, std::vector<double> vp, std::vector<double> vs, double mu, std::vector<double> depth, std::vector<double> dens, int nlay){
	// Recursive layer stacking from bottom to top to get R1212
	std::tuple<double,double,double,double,double> R;
	for(int n=nlay-1;n>=0;n--){
		if (verbose==1){
			cout << "Schicht: " << n+1 << "\n";
			cout << "Geschwindigkeiten: " << vp[n] << "\t" << vs[n] << "\n";
		}
		if (n==nlay-1)
			R = compute_T(w, c, vp[n], vs[n], mu);
		else {
			double dn = depth[n+1]-depth[n];
			if (verbose==1)
				cout << "Schichtdicke: " << dn << "m\n";
			R = compute_R(w, c, vp[n], vs[n], dn, dens[n], R);
		}
	}
	return(std::get<0>(R));
}

class R1212_root{
	// Root finding functor for R1212. Keeps all other varibles constants and only changes c (phase velocity) to finds roots. Uses the TOMS algorithm from the boost libraries
	public:
		R1212_root(double w_, std::vector<double> vp_, std::vector<double> vs_, double mu_, std::vector<double> depth_, std::vector<double> dens_, int nlay_):w(w_),vp(vp_),vs(vs_),mu(mu_),depth(depth_),dens(dens_),nlay(nlay_) {};
		double operator()(const double c){
			return compute_R1212(w, c, vp, vs, mu, depth, dens, nlay);
		}
	private:
		double w;
		std::vector<double> vp;
		std::vector<double> vs;
		double mu;
		std::vector<double> depth;
		std::vector<double> dens;
		int nlay;
};

struct TerminationCondition{
	// checks whether root bracketing has sufficiently converged
	bool operator() (double min, double max){
    return abs(min - max) < tolerance;
	}
};

std::vector<std::pair<double, double>> = get_waypoints(east1, north1, east2, north2, model_res){
	const double R = 6371008.771 // Mean earth radius
	a1 = atan2(cos(north2)
}

int main(){
	
	// Read phase delay time observations
	NcFile dtpFile("/home/bweise/bmw/MATLAB/matgsdf-master/dt_usarray.nc", NcFile::read);
	NcDim nperiodsIn = dtpFile.getDim("NumberOfPeriods");
	NcDim nstatsIn = dtpFile.getDim("NumberOfStations");
	NcDim nraysIn = dtpFile.getDim("NumberOfRays");
	int nperiods = nperiodsIn.getSize();
	int nstats = nstatsIn.getSize();
	int nrays = nraysIn.getSize();
	
	std::vector<double> periods(nperiods);
	std::vector<double> mpn(nstats);
	std::vector<double> mpe(nstats);
	std::vector<double> mpz(nstats);
	std::vector<double> src_rcvr_cmb(nrays*2.0);
	std::vector<double> dtp(nrays*nperiods);
	
	NcVar periodsIn=dtpFile.getVar("Periods");
	periodsIn.getVar(periods.data());
	NcVar mpnIn=dtpFile.getVar("MeasPosX");
	mpnIn.getVar(mpn.data());
	NcVar mpeIn=dtpFile.getVar("MeasPosY");
	mpeIn.getVar(mpe.data());
	NcVar mpzIn=dtpFile.getVar("MeasPosZ");
	mpzIn.getVar(mpz.data());
	NcVar src_rcvr_cmbIn=dtpFile.getVar("SRCombinations");
	src_rcvr_cmbIn.getVar(src_rcvr_cmb.data());
	NcVar dtpIn=dtpFile.getVar("dtp");
	dtpIn.getVar(dtp.data());
	
	int * dtp_dummy;
	NcVarAtt dummy = dtpIn.getAtt("_FillValue");
	dummy.getValues(dtp_dummy);
	
	double line = 12;
	cout << src_rcvr_cmb[line] << "\t" << src_rcvr_cmb[line + nrays] << "\n";
	cout << dtp[line] << "\t" << dtp[line + nrays] << "\t" << dtp[line + 2*nrays] << "\t" << dtp[line + 3*nrays] << "\t" << dtp[line + 4*nrays] << "\t" << dtp[line + 5*nrays] << "\t" << dtp[line + 6*nrays] << "\t" << dtp[line + 7*nrays] << "\n";
	
	// Read density, vs, vp from nc file
	NcFile densFile("/home/bweise/bmw/WINTERC/dens_na.nc", NcFile::read);
	NcFile vpFile("/home/bweise/bmw/WINTERC/vp_na.nc", NcFile::read);
	NcFile vsFile("/home/bweise/bmw/WINTERC/vs_na.nc", NcFile::read);
	
	// Get dimensions of model
	NcDim nxIn=densFile.getDim("Northing");
	NcDim nyIn=densFile.getDim("Easting");
	NcDim nzIn=densFile.getDim("Depth");
	int NX = nxIn.getSize();
	int NY = nyIn.getSize();
	int NZ = nzIn.getSize();
	
	// Define data variables
	std::vector<double> depth(NZ);
	std::vector<double> north(NX);
	std::vector<double> east(NY);
	std::vector<double> dens_all(NX*NY*NZ);
	std::vector<double> vp_all(NX*NY*NZ);
	std::vector<double> vs_all(NX*NY*NZ);
	
	// Read model from nc file
	NcVar depthIn=densFile.getVar("Depth");
	depthIn.getVar(depth.data());
	NcVar northingIn=densFile.getVar("Northing");
	northingIn.getVar(north.data());
	NcVar eastingIn=densFile.getVar("Easting");
	eastingIn.getVar(east.data());
	NcVar densIn=densFile.getVar("Density");
	densIn.getVar(dens_all.data());
	NcVar vsIn=vsFile.getVar("Vs");
	vsIn.getVar(vs_all.data());
	NcVar vpIn=vpFile.getVar("Vp");
	vpIn.getVar(vp_all.data());
	
	// Generate additional "0" to start vector of depths
	depth.pop_back();
	depth.insert(depth.begin(), 0);
	int nlay = depth.size();	// number of layers
	
	// get minimum extent of model cell
	double diff_e = 99999999999, diff_n = 99999999999;
	for(int estep=1; estep<NY; estep++){
		double diff = abs(east[estep] - east[estep-1]);
		if (diff < diff_e)
			diff_e = diff;
	}
	for(int nstep=1; nstep<NY; nstep++){
		double diff = abs(north[estep] - north[nstep-1]);
		if (diff < diff_n)
			diff_n = diff;
	}
	double model_res = diff_n;
	if (diff_e < diff_n)
		model_res = diff_e;
	
	// Open output file, write header line
	ofstream resultfile;
	resultfile.open ("dispersion.out");
	resultfile << "# Easting [m] \t Northing [m] \t Period [s] \t Phase velocity [m/s] \t Differenz [m/s] \t Anzahl Iterationen";
	
	// Sort vector of periods
	std::sort(periods.begin(), periods.end());
	
	std::vector<double> dispersion;
	
	for (int estep = 0; estep<NY; estep++){
		for (int nstep = 0; nstep<NX; nstep++){
			std::vector<double> dens;
			std::vector<double> vs;
			std::vector<double> vp;	
			bool lvz=0;
			for (int n=0; n<nlay; n++){
				// sort velocities, densities into 1D models 
				dens.push_back(dens_all[estep*NX+nstep+n*NX*NY]);
				vs.push_back(vs_all[estep*NX+nstep+n*NX*NY]*1000);
				//cout << vs[n] << "\n";
				// check if there's a low velocity zone
				if (n>0 & vs[n]<vs[n-1])
					lvz=1;
				//cout << lvz << "\n";
				vp.push_back(vp_all[estep*NX+nstep+n*NX*NY]*1000);
			}
			
			if (vs[0]==0)
				// This still needs some work. Computation of dispersion curves does not work if top layer has vs=0
				continue;
			else{				
				// Calculation of velocity limits 
				std::vector<double> c_lim = {0,0};
				c_lim[0] = newton_vr(vp[0], vs[0])/1.05;
				c_lim[1] = newton_vr(vp[nlay-1], vs[nlay-1])*1.05;
				double vsmin = *std::min_element(vs.begin(),vs.end());
				
				// step ratio for root bracketing
				double stepratio = (vsmin - c_lim[0])/(2*vsmin);	
	
				if(verbose==1){
					cout << "cmin: " << c_lim[0] << "\t cmax: " << c_lim[1] << "\n";
					cout << "Anzahl Schichten: " << nlay << "\n";
					cout << "Easting: " << east[estep] << "\t Northing: " << north[nstep] << "\n";
	
					for(int n=0; n<nlay-1; n++)
						cout << "Schicht " << n+1 << " von " << depth[n] << " bis "
							<< depth[n+1] << " m mit vs = " << vs[n] << " m/s, vp = "
							<< vp[n] << " m/s und " << dens[n] << " kg/m3 Dichte.\n";
					cout << "Letzte Schicht ab " << depth[nlay-1] << " m: vs = "
						<< vs[nlay-1] << " m/s, vp = " << vp[nlay-1]
						<< " m/s, Dichte = " << dens[nlay-1] << " kg/m3.\n";
		
					cout << "Kreisfrequenzen:\t";
				}
				
				// Conversion of periods to angular frequencies
				std::vector<double> w;
				for(int n=0; n<periods.size(); n++){
					w.push_back(2*M_PI/periods[n]);
					if (verbose==1)
						cout << w[n] << " rad/s\n";
				}
				
				// Shear modulus bottom layer
				double mu = pow(vs[nlay-1],2)*dens[nlay-1];
				if (verbose==1)
					cout << "Schermodul unterste Schicht: " << mu << "\n"
						<< "Polarisation von R1212 für geringe Geschwindigkeit: \n\n";
				
				// Compute initial R1212 polarization for large period (2000 s)	below fundamental mode
				double R1212 = compute_R1212(2*M_PI/2000, c_lim[0], vp, vs, mu, depth, dens, nlay);
				bool pol0 = signbit(R1212);
				
				double c_last=c_lim[0];	//
				
				// Loop over all periods/frequencies
				for(int freq=0; freq<w.size(); freq++){
					double c0, c1 = c_last;	// stores brackets
					bool pol1 = pol0;	// initial polarization of R1212
					double precision = 1;
					
					if (lvz==1) // If there is a LVZ, we have to start at the lowest possible velocity
						c1 = c_lim[0];
					
					// Loop to find root brackets, breaks when sign of R1212 changes
					while (pol0==pol1){
						cnt:;
						c0 = c1;	// set lower bracket to value of last iteration's upper bracket
						c1 = c0 + c0*(stepratio/precision);	// increase upper bracket by step ratio
			
						if (verbose==1)
							cout << "Aktuelle Kreisfreq. & Geschwindigkeit: " << w[freq] << "\t" << c1 << "\n";
						
						// Check polarization of R1212 for the upper bracket	
						R1212 = compute_R1212(w[freq], c1, vp, vs, mu, depth, dens, nlay);
						pol1 = signbit(R1212);
						
						// If a sign change is found check for mode skipping
						if (pol0!=pol1 & (c1-c0)>(2*tolerance)){
							double c2 = (c1+c0)/2;	// set new speed between brackets
							double delta = (c2-c0)/mode_skip_it;
							if (delta < (mode_skip_it*tolerance))
								delta = tolerance;
							if (verbose==1){
								cout << "c0: " << c0 << "\n";
								cout << "c2: " << c2 << "\n";
								cout << "c1: " << c1 << "\n";
								cout << "precision: " << precision << "\n";
								cout << "delta: " << delta << "\n";
							}
							// check sign of R1212 between brackets
							while (tolerance<(c2-c0)){
								R1212 = compute_R1212(w[freq], c2, vp, vs, mu, depth, dens, nlay);
								bool pol2 = signbit(R1212);
								// if mode skipping detected increase precision (-> decrease step ratio) and return to bracket search
								if (pol2==pol1){
									precision = precision * mode_skip_it;
									c1 = c0;
									if (verbose==1)
										cout << "Error: Mode skipping detected!\n";
									goto cnt;
								}
								// "Downward" search along c-axis for mode skipping (10 runs per default)
								c2 = c2-delta;
								if (verbose==1)
									cout << "New c2: " << c2 << "\n";
							}
						}		
					}
					// If a sign change is found, brackets c0 & c2 are passed to an instance of R1212_root (-> Boost TOMS root finding algorithm)
					std::pair<double, double> brackets; // stores refinded root brackets
					boost::uintmax_t max_iter=500;	// Maximum number of TOMS iterations (500 is probably way to much...)
					R1212_root root(w[freq], vp, vs, mu, depth, dens, nlay);
					brackets = boost::math::tools::toms748_solve(root, c0, c1, TerminationCondition(), max_iter);
					if (lvz == 0 & (brackets.first + brackets.second)/2 < c_last) {
						c1 = c_lim[0];
						precision = precision * mode_skip_it;
						if (verbose==1)
							cout << "Error: non-monotonous shape of dispersion curve!\n";
						goto cnt;
					}
					c_last = (brackets.first + brackets.second)/2;
					
					dispersion.push_back(c_last);
							
					// Write output to file
					resultfile << "\n" << east[estep] << "\t" << north[nstep] << "\t" << (2*M_PI)/w[freq] << "\t" << c_last << "\t" << brackets.second-brackets.first << "\t" << max_iter;
				} // end loop over frequencies
			} 
		} // end of loop over northing
	} // end of loop over easting
	
	// close file
	resultfile.close();
	
	// loop over all rays, computes phase delays
	for (int ray=0; ray<nrays; ray++){
		for (int freq=0; freq<w.size(); freq++){
			
			if (dtp[ray,freq]==*dtp_dummy)
				// if there is only a dummy value we can skip this period
				continue;
			else {				
				auto wpts = get_waypoints(mpe[src_rcvr_cmp[ray,0]], mpn[src_rcvr_cmp[ray,0]], mpe[src_rcvr_cmp[ray,1]], mpn[src_rcvr_cmp[ray,1]], model_res);
			}
				
		} // end loop frequencies
	} // end loop rays
	
	// end program
	return 0;
}

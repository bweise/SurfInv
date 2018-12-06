/* 1Dvs2vph
 * Berechnet aus einem 1D Vs,Vp,Dichte-Modell die Phasengeschwindigkeit
 * zu gegebenen Perioden.
 * Quellen: Haskell (1953), Dunkin (1965), Wathelet (2005),
 * Cercato (2007)
 */

using namespace std;
#include <iostream>
#include <vector>
//#include <math.h>
#include <complex>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <fstream>
typedef complex<double> dcomp;
const std::complex<double> i(0, 1.0);
bool verbose = 0;

// Perioden in sec
int np = 100;
std::pair<double,double> plim = {1,100};
std::vector<double> periods;
int ip=0;
//std::vector<double> periods = {10};

int nk = 1000; //Anzahl Wellenzahlen zum durchprobieren (min. 2)

// Definition 1D Modell
/*std::vector<double> depth = {0,25};	// Tiefe Schichtgrenzen [m]
std::vector<double> vp = {1350,2000};	// Vp für Schichten [m/s]
std::vector<double> vs = {250,1000};	// Vs für Schichten [m/s]
std::vector<double> dens = {2400,2400}; // Dichten [kg/m3]*/
std::vector<double> depth = {0,5000,15000,30000,55000};	// Tiefe Schichtgrenzen [m]
std::vector<double> vp = {5190,6060,6930,7790,8660};	// Vp für Schichten [m/s]
std::vector<double> vs = {3000,3500,4000,4500,5000};	// Vs für Schichten [m/s]
std::vector<double> dens = {2400,2625,2850,3075,3300}; // Dichten [kg/m3]*/

double compute_fvr(double vp, double vs, double vr){
	double fvr = 4.0-4.0*(pow(vr,2)/pow(vs,2))+pow(vr,4)/pow(vs,4)-4.0*sqrt(1-pow(vr,2)/pow(vp,2))*sqrt(1.0-pow(vr,2)/pow(vs,2));
	return fvr;
}

double compute_dfvr(double vp, double vs, double vr){
	double dfvr = -8.0*vr/pow(vs,2)+4.0*pow(vr,3)/pow(vs,4)+(4.0*vr*(pow(vp,2)+pow(vs,2)-2.0*pow(vr,2)))/(vp*vs*sqrt((vp-vr)*(vp+vr))*sqrt((vs-vr)*(vs+vr)));
	return dfvr;
}

double newton_vr(double vp, double vs){
	double vrlast = vs-1;
	double vr;
	double diff=99999.0;
	while(diff>0.0001){
		//cout << "vr: " << vrlast << "\t diff: " << diff << "\n";
		double fvr = compute_fvr(vp, vs, vrlast);
		double dfvr = compute_dfvr(vp, vs, vrlast);
		//cout << "fvr: " << fvr << "\t dfvr: " << dfvr << "\n";
		vr = vrlast - fvr/dfvr;
		diff = sqrt(pow(vr-vrlast,2));
		vrlast = vr;
	}
	if (verbose == 1)
		cout << "vr: " << vr << "\t final diff: " << diff << "\n\n";
	return vr;
}

std::tuple<dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,double,dcomp,dcomp,double> compute_util(double w, double k, double vp, double vs, double dn, bool botlay){
	dcomp mh, mk, hv, kv, SH, CH, SK, CK;
	double c = w/k;
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

std::tuple<dcomp,dcomp,dcomp,dcomp,dcomp> compute_T(double w, double k, double vp, double vs, double mu){
	double c = w/k;
	auto util = compute_util(w, k, vp, vs, 99999, 1);
	dcomp hv = std::get<0>(util);
	dcomp kv = std::get<1>(util);
	dcomp l = std::get<9>(util);
	
	dcomp fact = pow((-1.0)*pow(vs,2)/(2*mu*hv*kv*pow(w,2)),2);

	dcomp T1212 = pow(mu,2)*kv*hv*(pow(l,2)-4.0*pow(k,2)*kv*hv)*fact;
	dcomp T1213 = mu*pow(hv,2)*kv*(l-2.0*pow(k,2))*fact;
	dcomp iT1214 = 1.0*k*mu*hv*kv*(l-2.0*hv*kv)*fact;
	dcomp T1224 = mu*hv*pow(kv,2)*(2.0*pow(k,2)-l)*fact;
	dcomp T1234 = hv*kv*(pow(k,2)-hv*kv)*fact;
	
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

std::tuple<dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp> compute_G(double k, double dn, double w, double vp, double vs, double dens){
	double c = w/k;
	auto kvert = compute_util(w, k, vp, vs, dn, 0);
	dcomp SH = std::get<2>(kvert);
	dcomp CH = std::get<3>(kvert);
	dcomp SK = std::get<4>(kvert);
	dcomp CK = std::get<5>(kvert);
	dcomp gam = std::get<6>(kvert);
	dcomp hvnorm = std::get<7>(kvert);
	dcomp kvnorm = std::get<8>(kvert);
	
	dcomp G1212 = 2.0*gam * (1.0-gam) + (2.0*pow(gam,2)-2.0*gam+1.0) * CH*CK - (pow(1.0-gam,2) + pow(gam,2)*pow(hvnorm,2)*pow(kvnorm,2))*SH*SK;
	dcomp G1213 = (1.0/(dens*w*c))*(CH*SK - SH*CK*pow(hvnorm,2));
	dcomp iG1214 = (1.0/(dens*w*c))*((1.0 - 2.0*gam)*(1.0-CK*CH) + (1.0 - gam - gam*pow(hvnorm,2)*pow(kvnorm,2))*SH*SK);
	dcomp G1224 = (1.0/(dens*w*c))*(pow(kvnorm,2)*CH*SK - SH*CK);
	dcomp G1234 = (-1.0/(pow(dens,2)*pow(w,2)*pow(c,2)))*(2.0*(1.0 - CH*CK) + (1.0 + pow(kvnorm,2)*pow(hvnorm,2))*SK*SH);
	dcomp G1312 = dens*w*c*(pow(gam,2)*pow(kvnorm,2)*CH*SK - pow(1.0-gam,2)*SH*CK);
	dcomp G1313 = CH*CK;
	dcomp iG1314 = 1.0*((1.0 - gam)*SH*CK + gam*pow(kvnorm,2)*CH*SK);
	dcomp G1324 = (-1.0)*pow(kvnorm,2)*SH*SK;
	dcomp iG1412 = 1.0*dens*w*c*((3.0*pow(gam,2) - 2.0*pow(gam,3) - gam)*(1.0 - CH*CK)+(pow(1.0 - gam,3) - pow(gam,3)*pow(hvnorm,2)*pow(kvnorm,2))*SH*SK);
	dcomp iG1413 = (-1.0)*1.0*((1.0 - gam)*CH*SK + gam*pow(hvnorm,2)*SH*CK);
	dcomp G1414 = 1.0 - 2.0*gam*(1.0 - gam)*(1.0 - CH*CK) + (pow(1.0 - gam,2) + pow(gam,2)*pow(kvnorm,2)*pow(hvnorm,2))*SH*SK;
	dcomp G2412 = dens*w*c*(pow(1.0 - gam,2)*CH*SK - pow(gam,2)*SH*CK*pow(hvnorm,2));
	dcomp G2413 = (-1.0)*pow(hvnorm,2)*SH*SK;
	dcomp G3412 = (-1.0)*pow(dens,2)*pow(w,2)*pow(c,2)*(2.0*pow(gam,2)*pow(1.0 - gam,2)*(1.0 - CH*CK) + (pow(1.0 - gam,4)+pow(gam,4)*pow(hvnorm,2)*pow(kvnorm,2))*SH*SK);
	
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

std::tuple<dcomp,dcomp,dcomp,dcomp,dcomp> compute_R(double w, double k, double vp, double vs, double dn, double dens, std::tuple<dcomp,dcomp,dcomp,dcomp,dcomp> T){
	dcomp T1212 = std::get<0>(T);
	dcomp T1213 = std::get<1>(T);
	dcomp iT1214 = std::get<2>(T);
	dcomp T1224 = std::get<3>(T);
	dcomp T1234 = std::get<4>(T);
	auto G = compute_G(k,dn,w,vp,vs,dens);
	dcomp G1212 = std::get<0>(G);
	dcomp G1213 = std::get<1>(G);
	dcomp iG1214 = std::get<2>(G);
	dcomp G1224 = std::get<3>(G);
	dcomp G1234 = std::get<4>(G);
	dcomp G1312 = std::get<5>(G);
	dcomp G1313 = std::get<6>(G);
	dcomp iG1314 = std::get<7>(G);
	dcomp G1324 = std::get<8>(G);
	dcomp iG1412 = std::get<9>(G);
	dcomp iG1413 = std::get<10>(G);
	dcomp G1414 = std::get<11>(G);
	dcomp G2412 = std::get<12>(G);
	dcomp G2413 = std::get<13>(G);
	dcomp G3412 = std::get<14>(G);
	
	dcomp R1212 = T1212*G1212 + T1213*G1312 - 2.0*iT1214*iG1412 + T1224*G2412 + T1234*G3412;
	dcomp R1213 = T1212*G1213 + T1213*G1313 - 2.0*iT1214*iG1413 + T1224*G2413 + T1234*G2412;
	dcomp iR1214 = T1212*iG1214 + T1213*iG1314 + iT1214*(2.0*G1414-1.0) + T1224*iG1413 + T1234*iG1412;
	dcomp R1224 = T1212*G1224 + T1213*G1324 - 2.0*iT1214*iG1314 + T1224*G1313 + T1234*G1312;
	dcomp R1234 = T1212*G1234 + T1213*G1224 - 2.0*iT1214*iG1214 + T1224*G1213 + T1234*G1212;
	
	double tmpmax, maxR = std::real(R1212);
	if(maxR<0)
		maxR = (-1.0)*maxR;
	tmpmax = std::real(R1213);
	if (tmpmax<0)
		tmpmax = (-1.0)*tmpmax;
	if (tmpmax>maxR)
		maxR = tmpmax;
	tmpmax = std::imag(iR1214);
	if (tmpmax<0)
		tmpmax = (-1.0)*tmpmax;
	if (tmpmax>maxR)
		maxR = tmpmax;
	tmpmax = std::real(R1224);
	if (tmpmax<0)
		tmpmax = (-1.0)*tmpmax;
	if (tmpmax>maxR)
		maxR = tmpmax;
	tmpmax = std::real(R1234);
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

int main()
{	
	while (ip<=np){
		double p = std::get<0>(plim)+ip*(std::get<1>(plim)-std::get<0>(plim))/np;
		periods.push_back(p);
		++ip;
	}
	
	
	ofstream resultfile;
	resultfile.open ("dispersion.out");
	resultfile << "# Output of 1Dvs2vph:\n# "
		<< nk << " wavenumbers at " << periods.size() << " Frequencies tested.\n"
		<<"# Frequency [1/s] \t Phase velocity [m/s] \t R1212 \n";
	int nlay = depth.size();
	
	std::vector<double> vs_sort = vs;
	std::sort(vs_sort.begin(), vs_sort.end());
	std::vector<double> c_lim = {0,0};
	c_lim[0] = newton_vr(vp[0], vs[0])/1.05;
	c_lim[1] = newton_vr(vp[nlay-1], vs[nlay-1])*1.05;
	
	if(verbose==1){
	cout << "cmin: " << c_lim[0] << "\t cmax: " << c_lim[1] << "\n";
	
	cout << "Anzahl Schichten: " << nlay << "\n";
		for(int n=0; n<nlay-1; n++)
			cout << "Schicht " << n+1 << " von " << depth[n] << " bis "
				<< depth[n+1] << " m mit vs = " << vs[n] << " m/s, vp = "
				<< vp[n] << " m/s und " << dens[n] << " kg/m3 Dichte.\n";
		cout << "Letzte Schicht ab " << depth[nlay-1] << " m: vs = "
			<< vs[nlay-1] << " m/s, vp = " << vp[nlay-1]
			<< " m/s, Dichte = " << dens[nlay-1] << " kg/m3.\n";
		
		cout << "Kreisfrequenzen:\t";
	}
	std::vector<double> w;
	for(int n=0; n<periods.size(); n++){
		w.push_back(2*M_PI/periods[n]);
		if (verbose==1)
			cout << w[n] << " rad/s\n";
	}
	
	double mu = pow(vs[nlay-1],2)*dens[nlay-1]; // Shear modulus unterste Schicht
	if (verbose==1)
		cout << "Schermodul unterste Schicht: " << mu << "\n\n";

	for(int freq=0; freq<w.size(); freq++){
		std::vector<double> k_lim = {w[freq]/c_lim[0],w[freq]/c_lim[1]};
		int kk;
		for(int kint=0; kint<nk; kint++){
			kk=kint;
			double k=k_lim[0]-kint*(k_lim[0]-k_lim[1])/(nk-1);
			if (verbose==1)
				cout << "Aktuelle Kreisfreq. & Wellenzahl: " << w[freq] << "\t" << k << "\n";
			
			std::tuple<dcomp,dcomp,dcomp,dcomp,dcomp> R;
			
			for(int n=nlay-1;n>=0;n--){
				if (verbose==1){
					cout << "Schicht: " << n+1 << "\n";
					cout << "Geschwindigkeiten: " << vp[n] << "\t" << vs[n] << "\n";
				}
				if (n==nlay-1){
					R = compute_T(w[freq], k, vp[n], vs[n], mu);
				}
				else {
					double dn = depth[n+1]-depth[n];
					if (verbose==1)
						cout << "Schichtdicke: " << dn << "m\n";
					R = compute_R(w[freq], k, vp[n], vs[n], dn, dens[n], R);
				}
			}
			if (freq==w.size()-1 & kk==nk-1){
				resultfile << w[freq]/(2*M_PI) << "\t" << w[freq]/k << "\t" << std::real(std::get<0>(R));
				resultfile.close();
			}
			else
				resultfile << w[freq]/(2*M_PI) << "\t" << w[freq]/k << "\t" << std::real(std::get<0>(R)) << "\n";
		}
	}
	
	return 0;
}


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
const std::complex<double> i(0, 1);

// Perioden in sec
int np = 1000;
std::pair<double,double> plim = {0.1,0.03333333};
std::vector<double> periods;
int ip=0;
//std::vector<double> periods = {0.1,0.5};

int nk = 10000; //Anzahl Wellenzahlen zum durchprobieren (min. 2)

// Definition 1D Modell
std::vector<double> depth = {0,25};	// Tiefe Schichtgrenzen [m]
std::vector<double> vp = {1350,2000};	// Vp für Schichten [m/s]
std::vector<double> vs = {250,1000};	// Vs für Schichten [m/s]
std::vector<double> dens = {2400,2400}; // Dichten [kg/m3]

double compute_fvr(double vp, double vs, double vr){
	double fvr = 4.0-4.0*(pow(vr,2)/pow(vs,2))+pow(vr,4)/pow(vs,4)-4.0*sqrt(1-pow(vr,2)/pow(vp,2))*sqrt(1.0-pow(vr,2)/pow(vs,2));
	return fvr;
}

double compute_dfvr(double vp, double vs, double vr){
	double dfvr = -8.0*vr/pow(vs,2)+4.0*pow(vr,3)/pow(vs,4)+(4.0*vr*(pow(vp,2)+pow(vs,2)-2.0*pow(vr,2)))/(vp*vs*sqrt((vp-vr)*(vp+vr))*sqrt((vs-vr)*(vs+vr)));
	return dfvr;
}

double newton_vr(double vp, double vs){
	double vrlast = vs-(vs*0.1);
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
	//cout << "final vr: " << vr << "\t final diff: " << diff << "\n\n\n";
	return vr;
}

std::tuple<dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,double,dcomp,dcomp,double> compute_util(double w, double k, double vp, double vs, double dn, bool botlay){
	dcomp mh, mk, hv, kv, SH, CH, SK, CK;
	if (w/k < vp | botlay == 1){
		mh = sqrt(pow(k,2)-pow(w/vp,2));
		hv = mh;
		if (botlay == 0){
		SH = (k/mh)*sinh(mh*dn);
		CH = cosh(mh*dn);
		}
	}
	else{
		mh = sqrt(pow(w/vp,2)-pow(k,2));
		hv = mh*i;
		SH = (k/mh)*sin(mh*dn);
		CH = cos(mh*dn);
	}
	if (w/k < vs | botlay == 1){
		mk = sqrt(pow(k,2)-pow(w/vs,2));
		kv = mk;
		if (botlay == 0){
		SK = (k/mk)*sinh(mk*dn);
		CK = cosh(mk*dn);
		}
	}
	else{
		mk = sqrt(pow(w/vs,2)-pow(k,2));
		kv = mk*i;
		SK = (k/mk)*sin(mk*dn);
		CK = cos(mk*dn);
	}
	double gam = 2.0*pow(vs,2)*pow(k/w,2);
	dcomp hvnorm = hv/k;
	dcomp kvnorm = kv/k;
	double l = 2.0*pow(k,2)-pow(w/vs,2);
	return std::make_tuple(mh, mk, hv, kv, SH, CH, SK, CK, gam, hvnorm, kvnorm, l);
}

std::tuple<dcomp,dcomp,dcomp,dcomp,dcomp> compute_T(double w, double k, double vp, double vs, double mu, double dens, double dn){
	auto util = compute_util(w, k, vp, vs, dn, 1);
	dcomp hv = std::get<2>(util);
	dcomp kv = std::get<3>(util);
	dcomp l = std::get<11>(util);
	//cout << "kvert: " << hv << "\t" << kv << "\n";
	dcomp T1212 = (pow(vs,4)/(4.0*pow(w,4)))*((pow(l,2)/(kv*hv))-4.0*pow(k,2));
	dcomp T1213 = (-1.0)/(4.0*dens*pow(w,2)*kv);
	dcomp T1214 = (i*pow(vs,2)*k/(4.0*dens*pow(w,4)))*((l/(kv*hv))-2.0);
	dcomp T1224 = (1/(4.0*pow(dens,2)*pow(w,4)))*((pow(k,2)/(kv*hv))-1.0);
	dcomp T1234 = (pow(vs,4)/(4.0*pow(mu,2)*pow(w,4)))*((pow(k,2)/(kv*hv))-1.0);
	/*cout << "T-Komponenten: " << T1212 << "\t"
		<< T1213 << "\t"
		<< T1214 << "\t"
		<< T1224 << "\t"
		<< T1234 << "\n\n\n";*/
	return std::make_tuple(T1212,T1213,T1214,T1224,T1234);
}

std::tuple<dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp> compute_G(double k, double dn, double w, double vp, double vs, double dens){
	auto kvert = compute_util(w, k, vp, vs, dn, 0);
	dcomp mh = std::get<0>(kvert);
	dcomp mk = std::get<1>(kvert);
	dcomp hv = std::get<2>(kvert);
	dcomp kv = std::get<3>(kvert);
	dcomp SH = std::get<4>(kvert);
	dcomp CH = std::get<5>(kvert);
	dcomp SK = std::get<6>(kvert);
	dcomp CK = std::get<7>(kvert);
	dcomp gam = std::get<8>(kvert);
	dcomp hvnorm = std::get<9>(kvert);
	dcomp kvnorm = std::get<10>(kvert);
	dcomp l = std::get<11>(kvert);
	/*cout << "kvert: " << hv << "\t" << kv << "\n"
		<< "kvert (normiert): " << hvnorm << "\t" << kvnorm << "\n"
		<< "SH: " << SH << "\t CH: " << CH << "\n"
		<< "SK: " << SK << "\t CK: " << CK << "\n"
		<< "gamma: " << gam << "\t" << l << "\n";*/
	dcomp G1212 = 2.0*gam*(1.0-gam)+(2.0*pow(gam,2)-2.0*gam+1.0)*CH*CK-(pow(1.0-gam,2)+pow(gam,2)*pow(hvnorm,2)*pow(kvnorm,2))*SH*SK;
	dcomp G1213 = (k/(dens*pow(w,2)))*(CH*SK-SH*CK*pow(kvnorm,2));
	dcomp iG1214 = (i*k/(dens*pow(w,2)))*((1.0-2.0*gam)*(1.0-CK*CH)+(1.0-gam-gam*pow(hvnorm,2)*pow(kvnorm,2))*SH*SK);
	dcomp G1224 = (k/(dens*pow(w,2)))*(pow(kvnorm,2)*CH*SK-SH*CK);
	dcomp G1234 = (pow((-1.0)*k,2)/(pow(dens,2)*pow(w,4)))*(2.0*(1.0-CH*CK)+(1.0+pow(kvnorm,2)*pow(hvnorm,2))*SK*SH);
	dcomp G1312 = (dens*pow(w,2)/k)*(pow(gam,2)*pow(kvnorm,2)*CH*SK-pow(1.0-gam,2)*SH*CK);
	dcomp G1313 = CH*CK;
	dcomp iG1314 = i*((1.0-gam)*SH*CK+gam*pow(kvnorm,2)*CH*SK);
	dcomp G1324 = (-1.0)*pow(kvnorm,2)*SH*SK;
	dcomp iG1412 = (i*dens*pow(w,2)/k)*((3.0*pow(gam,2)-2.0*pow(gam,3)-gam)*(1.0-CH*CK)+(pow(1.0-gam,3)-pow(gam,3)*pow(hvnorm,2)*pow(kvnorm,2))*SH*SK);
	dcomp iG1413 = (-1.0)*i*((1.0-gam)*CH*SK+gam*pow(hvnorm,2)*SH*CK);
	dcomp G1414 = 1.0-2.0*gam*(1.0-gam)*(1.0-CH*CK)+(pow(1.0-gam,2)+pow(gam,2)*pow(kvnorm,2)*pow(hvnorm,2))*SH*SK;
	dcomp G2412 = (dens*pow(w,2)/k)*(pow(1.0-gam,2)*CH*SK-pow(gam,2)*SH*CK*pow(hvnorm,2));
	dcomp G2413 = (-1.0)*pow(hvnorm,2)*SH*SK;
	dcomp G3412 = (((-1.0)*pow(dens,2)*pow(w,4))/pow(k,2))*(2.0*pow(gam,2)*pow(1.0-gam,2)*(1.0-CH*CK)+(pow(1.0-gam,4)+pow(gam,4)*pow(hvnorm,2)*pow(kvnorm,2))*SH*SK);
	/*cout << "G-Komponenten:\n" << G1212 << "\t" << G1213 << "\t" << iG1214 << "\t" << G1224 << "\t" << G1234 << "\n"
		<< G1312 << "\t" << G1313 << "\t" << iG1314 << "\t" << G1324 << "\n"
		<< iG1412 << "\t" << iG1413 << "\t" << "\n"
		<< G2412 << "\t" << G2413 << "\n"
		<< G3412 << "\t" << "\n";*/
	return std::make_tuple(G1212,G1213,iG1214,G1224,G1234,G1312,G1313,iG1314,G1324,iG1412,iG1413,G1414,G2412,G2413,G3412);
}

std::tuple<dcomp,dcomp,dcomp,dcomp,dcomp> compute_R(double w, double k, double vp, double vs, double dn, double dens, std::tuple<dcomp,dcomp,dcomp,dcomp,dcomp> T){
	dcomp T1212 = std::get<0>(T);
	dcomp T1213 = std::get<1>(T);
	dcomp T1214 = std::get<2>(T);
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
	dcomp R1212 = T1212*G1212+T1213*G1312-2.0*T1214*iG1412+T1224*G2412-T1234*G3412;
	dcomp R1213 = T1212*G1213+T1213*G1313-2.0*T1214*iG1413-T1224*G2413+T1234*G2412;
	dcomp R1214 = T1212*iG1214+T1213*iG1314+T1214*(2.0*G1414-1.0)-T1224*iG1413+T1234*iG1412;
	dcomp R1224 = T1212*G1224+T1213*G1324-2.0*T1214*iG1314+T1224*G1313+T1234*G1312;
	dcomp R1234 = T1212*G1234+T1213*G1224-2.0*T1214*iG1214+T1224*G1213+T1234*G1212;
	/*cout << "R-Komponenten: " << R1212 << "\t"
	<< R1213 << "\t"
	<< R1214 << "\t"
	<< R1224 << "\t"
	<< R1234 << "\n\n\n";*/
	return std::make_tuple(R1212,R1213,R1214,R1224,R1234);
}

int main()
{	
	while (ip<=np){
	double p = std::get<0>(plim)+ip*(std::get<1>(plim)-std::get<0>(plim))/np;
	periods.push_back(p);
	++ip;
	}
	
	
	ofstream resultfile;
	resultfile.open ("example.txt");
	resultfile << "# Output of 1Dvs2vph:\n# "
		<< nk << " wavenumbers at " << periods.size() << " Frequencies tested.\n"
		<<"# Frequency [1/s] \t Phase velocity [m/s] \t R1212 \n";
	int nlay = depth.size();
	
	std::vector<double> vs_sort = vs;
	std::sort(vs_sort.begin(), vs_sort.end());
	std::vector<double> c_lim = {0,vs_sort[nlay-1]-0.01};
	c_lim[0] = newton_vr(vp[0], vs[0])/1.05;
	//cout << c_lim[0] << "\t" << c_lim[1] << "\n";
	
	/*cout << "Anzahl Schichten: " << nlay << "\n";
	for(int n=0; n<nlay-1; n++)
		cout << "Schicht " << n+1 << " von " << depth[n] << " bis "
			<< depth[n+1] << " m mit vs = " << vs[n] << " m/s, vp = "
			<< vp[n] << " m/s und " << dens[n] << " kg/m3 Dichte.\n";
	cout << "Letzte Schicht ab " << depth[nlay-1] << " m: vs = "
		<< vs[nlay-1] << " m/s, vp = " << vp[nlay-1]
		<< " m/s, Dichte = " << dens[nlay-1] << " kg/m3.\n";*/
		
	//cout << "Kreisfrequenzen:\t";
	std::vector<double> w;
	for(int n=0; n<periods.size(); n++){
		w.push_back(2*M_PI/periods[n]);
		//cout << w[n] << " rad/s\n";
	}
	
	double mu = vs[nlay-1]*dens[nlay-1]; // Shear modulus unterste Schicht
	//cout << "Schermodul unterste Schicht: " << mu << "\n\n\n";

	for(int freq=0; freq<w.size(); freq++){
		std::vector<double> k_lim = {w[freq]/c_lim[0],w[freq]/c_lim[1]};
		int kk;
		for(int kint=0; kint<nk; kint++){
			kk=kint;
			double k=k_lim[0]-kint*(k_lim[0]-k_lim[1])/(nk-1);
			//cout << "Aktuelle Kreisfreq. & Wellenzahl: " << w[freq] << "\t" << k << "\n";
			
			std::tuple<dcomp,dcomp,dcomp,dcomp,dcomp> R;
			
			for(int n=nlay-1;n>=0;n--){
				//cout << "Schicht: " << n+1 << "\n";
				//cout << "Geschwindigkeiten: " << vp[n] << "\t" << vs[n] << "\n";
				if (n==nlay-1){
					R = compute_T(w[freq], k, vp[n], vs[n], mu, dens[n], 999);
				}
				else {
					double dn = depth[n+1]-depth[n];
					//cout << "Schichtdicke: " << dn << "m\n";
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


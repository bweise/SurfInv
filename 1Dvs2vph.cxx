/* 1Dvs2vph
 * Berechnet aus einem 1D Vs,Vp,Dichte-Modell die Phasengeschwindigkeit
 * zu gegebenen Perioden.
 * Quellen: Haskell (1953), Dunkin (1965), Wathelet (2005),
 * Cercato (2007)
 */

using namespace std;
#include <iostream>
#include <vector>
#include <math.h>
#include <complex>
#include <cmath>
#include <tuple>
#include <algorithm>
typedef complex<double> dcomp;
const std::complex<double> i(0, 1);

// Perioden in sec
//std::vector<double> periods = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
std::vector<double> periods = {0.1};

int nk = 2; //Anzahl Wellenzahlen zum durchprobieren (min. 2)

// Definition 1D Modell
std::vector<double> depth = {0,12500,25000,37500,50000};	// Tiefe Schichtgrenzen [m]
std::vector<double> vp = {5190,6060,6930,7790,8660};	// Vp für Schichten [m/s]
std::vector<double> vs = {3000,3500,4000,4500,5000};	// Vs für Schichten [m/s]
std::vector<double> dens = {2400,2625,2850,3075,3300}; // Dichten [kg/m3]

double compute_fvr(double vp, double vs, double vr){
	double fvr = 4-4*(pow(vr,2)/pow(vs,2))+pow(vr,4)/pow(vs,4)-4*sqrt(1-pow(vr,2)/pow(vp,2))*sqrt(1-pow(vr,2)/pow(vs,2));
	return fvr;
}

double compute_dfvr(double vp, double vs, double vr){
	double dfvr = -8*vr/pow(vs,2)+4*pow(vr,3)/pow(vs,4)+(4*vr*(pow(vp,2)+pow(vs,2)-2*pow(vr,2)))/(vp*vs*sqrt((vp-vr)*(vp+vr))*sqrt((vs-vr)*(vs+vr)));
	return dfvr;
}

double newton_vr(double vp, double vs){
	double vrlast = vs-(vs*0.1);
	double vr;
	double diff=99999.0;
	while(diff>0.0001){
		cout << "vr: " << vrlast << "\t diff: " << diff << "\n";
		double fvr = compute_fvr(vp, vs, vrlast);
		double dfvr = compute_dfvr(vp, vs, vrlast);
		cout << "fvr: " << fvr << "\t dfvr: " << dfvr << "\n";
		vr = vrlast - fvr/dfvr;
		diff = sqrt(pow(vr-vrlast,2));
		vrlast = vr;
	}
	cout << "final vr: " << vr << "\t final diff: " << diff << "\n\n\n";
	return vr;
}

std::pair<dcomp,dcomp> compute_kvert(double w, double k, double vp, double vs){
	dcomp hv = 2*pow(k,2)-(pow(w,2)/pow(vp,2));
	dcomp kv = 2*pow(k,2)-(pow(w,2)/pow(vs,2));
	return std::make_pair(sqrt(hv), sqrt(kv));
}

std::tuple<dcomp,dcomp,dcomp,dcomp,dcomp> compute_T(double w, double k, double vp, double vs, double mu){
	auto kvert = compute_kvert(w, k, vp, vs);
	dcomp hv = std::get<0>(kvert);
	dcomp kv = std::get<1>(kvert);
	cout << "kvert: " << std::get<0>(kvert) << "\t" << std::get<1>(kvert) << "\n";
	dcomp ln = pow(k,2)+pow(kv,2);
	dcomp t_factor = (-1.0)*pow(vs,2)/(2*mu*hv*kv*pow(w,2));
	//cout << ln << "\t" << t_factor << "\n";
	dcomp t11 = (2.0*i*mu*k*hv*kv)*t_factor;
	dcomp t12 = (mu*ln*kv)*t_factor;
	dcomp t13 = (hv*kv)*t_factor;
	dcomp t14 = (i*k*kv)*t_factor;
	dcomp t21 = ((-1.0)*mu*ln*kv)*t_factor;
	dcomp t22 = (2.0*i*mu*k*hv*kv)*t_factor;
	dcomp t23 = (i*k*hv)*t_factor;
	dcomp t24 = ((-1.0)*hv*kv)*t_factor;
	cout << "T-Komponenten:\n" << t11 << "\t" << t12 << "\t" << t13 << "\t" << t14 << "\n"
		<< t21 << "\t" << t22 << "\t" << t23 << "\t" << t24 << "\n";
	dcomp T1212 = t11*t22 - t12*t21;
	dcomp T1213 = t11*t23 - t13*t21;
	dcomp T1214 = t11*t24 - t14*t21;
	dcomp T1224 = t12*t24 - t14*t22;
	dcomp T1234 = t13*t24 - t14*t23;
	return std::make_tuple(T1212,T1213,T1214,T1224,T1234);
}

std::tuple<dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp> compute_G(double k, double dn, double w, double vp, double vs, double dens){
	
	auto kvert = compute_kvert(w, k, vp, vs);
	dcomp hv = std::get<0>(kvert);
	dcomp kv = std::get<1>(kvert);
	cout << "kvert: " << hv << "\t" << kv << "\n";
	
	dcomp hvnorm = hv/k;
	dcomp kvnorm = kv/k;
	cout << "kvert (normiert): " << hvnorm << "\t" << kvnorm << "\n";
	
	dcomp SH = 0;
	dcomp CH = 0;
	dcomp SK = 0;
	dcomp CK = 0;
	
	
	if (imag(hvnorm)==0){
		SH = 0.5*(1.0-exp((-2.0)*dn*hv))/hvnorm;
		CH = 0.5*(1.0+exp((-2.0)*dn*hv));
	}
	else{
		SH = sin((-1.0)*i*dn*hv)/hvnorm;
		CH = cos((-1.0)*i*dn*hv);
	}
	if (imag(kvnorm)==0){
		SK = 0.5*(1.0-exp((-2.0)*dn*kv))/kvnorm;
		CK = 0.5*(1.0+exp((-2.0)*dn*kv));
	}
	else{
		SK = sin((-1.0)*i*dn*kv)/kvnorm;
		cout << kv << " " << kvnorm << " " << dn << " " << i << " " << sin(i*dn*kv)/kvnorm << "\n";
		CK = cos((-1.0)*i*dn*kv);
	}
	cout << "SH: " << SH << "\t CH: " << CH << "\n"
		<< "SK: " << SK << "\t CK: " << CK << "\n";
	
	double gam = 2.0*pow(k,2)/pow((w/vp),2);
	double a1 = pow(gam,2)-2.0*gam+1.0;
	dcomp a2 = pow(hvnorm,2)*pow(kvnorm,2);
	double a3 = pow(gam,2)+a1;
	double a4 = 1.0-gam;
	dcomp a5 = pow(gam,2)*a2;
	dcomp expcorr = exp((-1.0)*hv*dn-kv*dn);
	double c1 = dens*pow(w,2)/k;
	double c2 = 1.0/c1;
	
	cout << "gamma: " << gam << "\n"
		<< "a1: " << a1 << "\t a2: " << a2 << "\t a3: " << a3 << "\t a4: " << a4 << "\t a5: " << a5 << "\n"
		<< "expcorr: " << expcorr << "\n"
		<< "c1: " << c1 << "\t c2: " << c2 << "\n";
	
	dcomp G1212 = a3*CH*CK-(a1+a5)*SH*SK-(a3-1)*expcorr;
	dcomp G1213 = c2*(CH*SK-pow(hvnorm,2)*SH*CK);
	cout << CH*SK << "\t" << SH*CK << "\n";
	dcomp iG1214 = i*c2*((a1-pow(gam,2))*(expcorr-CH*CK)+(a4-gam*a2)*SH*SK);
	dcomp G1224 = c2*(pow(kvnorm,2)*CH*SK-SH*CK);
	dcomp G1234 = pow(c2,2)*(2.0*CH*CK+(1.0+a2)*SH*SK);
	dcomp G1312 = c1*(pow(gam,2)*pow(kvnorm,2)*CH*SK-a1*SH*CK);
	dcomp iG1314 = i*(a4*SH*CK+gam*pow(kvnorm,2)*CH*SK);
	dcomp G1324 = pow(kvnorm,2)*SH*SK;
	dcomp iG1412 = i*c1*((a1-a4)*(a4-gam)*(expcorr-CH*CK)+(a4*a1-gam*a5)*SH*SK);
	dcomp iG1413 = i*(gam*pow(hvnorm,2)*SH*CH+a4*CH*SK);
	dcomp G1423 = CH*CK-G1212;
	dcomp G2412 = c1*(a1*CH*SK-pow(gam,2)*pow(hvnorm,2)*SH*CK);
	dcomp G2413 = pow(hvnorm,2)*SH*SK;
	dcomp G3412 = pow(c1,2)*(2.0*pow(gam,2)*a1*CH*CK+(pow(a1,2)+pow(gam,2)*a5)*SH*SK);
	cout << "G-Komponenten:\n" << G1212 << "\t" << G1213 << "\t" << iG1214 << "\t" << G1224 << "\t" << G1234 << "\n"
		<< G1312 << "\t" << iG1314 << "\t" << G1324 << "\n"
		<< iG1412 << "\t" << iG1413 << "\t" << G1423 << "\n"
		<< G2412 << "\t" << G2413 << "\n"
		<< G3412 << "\t" << "\n";
	return std::make_tuple(G1212,G1213,iG1214,G1224,G1234,G1312,iG1314,G1324,iG1412,iG1413,G1423,G2412,G2413,G3412,CH,CK,expcorr);
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
	dcomp iG1314 = std::get<6>(G);
	dcomp G1324 = std::get<7>(G);
	dcomp iG1412 = std::get<8>(G);
	dcomp iG1413 = std::get<9>(G);
	dcomp G1423 = std::get<10>(G);
	dcomp G2412 = std::get<11>(G);
	dcomp G2413 = std::get<12>(G);
	dcomp G3412 = std::get<13>(G);
	dcomp CH = std::get<14>(G);
	dcomp CK = std::get<15>(G);
	dcomp expcorr = std::get<16>(G);
	dcomp R1212 = T1212*G1212+(T1213*G1312-2.0*T1214*iG1412+T1224*G2412-T1234*G3412)/pow(w,1);
	dcomp R1213 = pow(w,2)*T1212*G1213+T1213*CH*CK-2.0*T1214*iG1413-T1224*G2413+T1234*G2412;
	dcomp R1214 = pow(w,2)*T1212*iG1214+T1213*iG1314+T1214*(2.0*G1423+expcorr)-T1224*iG1413+T1234*iG1412;
	dcomp R1224 = pow(w,2)*T1212*G1224+T1213*G1324-2.0*T1214*iG1314+T1224*CH*CK+T1234*G1312;
	dcomp R1234 = (-1.0)*pow(w,2)*T1212*G1234+T1213*G1224-2.0*T1214*iG1214+T1224*G1213+T1234*G1212;
	return std::make_tuple(R1212,R1213,R1214,R1224,R1234);
}

int main()
{	
	int nlay = depth.size();
	
	std::vector<double> vs_sort = vs;
	std::sort(vs_sort.begin(), vs_sort.end());
	std::vector<double> c_lim = {0,vs_sort[nlay-1]};
	c_lim[0] = newton_vr(vp[0], vs[0]);
	
	cout << "Anzahl Schichten: " << nlay << "\n";
	for(int n=0; n<nlay-1; n++)
		cout << "Schicht " << n+1 << " von " << depth[n] << " bis "
			<< depth[n+1] << " m mit vs = " << vs[n] << " m/s, vp = "
			<< vp[n] << " m/s und " << dens[n] << " kg/m3 Dichte.\n";
	cout << "Letzte Schicht ab " << depth[nlay-1] << " m: vs = "
		<< vs[nlay-1] << " m/s, vp = " << vp[nlay-1]
		<< " m/s, Dichte = " << dens[nlay-1] << " kg/m3.\n";
		
	cout << "Kreisfrequenzen:\t";
	std::vector<double> w;
	for(int n=0; n<periods.size(); n++){
		w.push_back(2*M_PI/periods[n]);
		cout << w[n] << " rad/s\n";
	}
	
	double mu = vs[nlay-1]*dens[nlay-1]; // Shear modulus unterste Schicht
	cout << "Schermodul unterste Schicht: " << mu << "\n\n\n";

	for(int freq=0; freq<w.size(); freq++){
		std::vector<double> k_lim = {w[freq]/c_lim[0],w[freq]/c_lim[1]};
		for(int kint=0; kint<nk; kint++){
			
			double k=k_lim[0]-kint*(k_lim[0]-k_lim[1])/(nk-1);
			cout << "Aktuelle Kreisfreq. & Wellenzahl: " << w[freq] << "\t" << k << "\n";
			
			std::tuple<dcomp,dcomp,dcomp,dcomp,dcomp> R;
			
			for(int n=nlay-1;n>=0;n--){
				cout << "Schicht: " << n+1 << "\n";
				cout << "Geschwindigkeiten: " << vp[n] << "\t" << vs[n] << "\n";
				
				if (n==nlay-1){
					R = compute_T(w[freq], k, vp[n], vs[n], mu);
					cout << "R-Komponenten: " << std::get<0>(R) << "\t"
						<< std::get<1>(R) << "\t"
						<< std::get<2>(R) << "\t"
						<< std::get<3>(R) << "\t"
						<< std::get<4>(R) << "\n\n\n";
				}
				else {
					double dn = depth[n+1]-depth[n];
					cout << "Schichtdicke: " << dn << "m\n";
					R = compute_R(w[freq], k, vp[n], vs[n], dn, dens[n], R);
					cout << "R-Komponenten: " << std::get<0>(R) << "\t"
						<< std::get<1>(R) << "\t"
						<< std::get<2>(R) << "\t"
						<< std::get<3>(R) << "\t"
						<< std::get<4>(R) << "\n\n\n";
				}
			}
			
		}
	}
	
	return 0;
}


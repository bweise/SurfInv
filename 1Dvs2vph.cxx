/*
 * 1Dvs2vph
 * Berechnet aus einem 1D Vs,Vp,Dichte-Modell die Phasengeschwindigkeit
 * zu gegebenen Perioden.
 * Quellen: Haskell (1953), Dunkin (1965), Wathelet (2005),
 * Cercato (2007), Fang (2015)
 */

using namespace std;
#include <iostream>
#include <vector>
#include <math.h>
#include <complex>
#include <cmath>
#include <tuple>
typedef complex<double> dcomp;
const std::complex<double> i(0, 1);

// Perioden in sec
//std::vector<double> periods = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
std::vector<double> periods = {0.1};

// Kreiswellenzahlen-Limits in rad/m
std::vector<double> k_lim = {20*M_PI/1000,20*M_PI/200};
int nk = 3; //Anzahl Wellenzahlen zum durchprobieren

// Definition 1D Modell
std::vector<double> depth = {0,25};	// Tiefe Schichtgrenzen [m]
std::vector<double> vs = {250,1000};	// Vs für Schichten [m/s]
std::vector<double> vp = {1350,2000};	// Vp für Schichten [m/s]
std::vector<double> dens = {2400,2400}; // Dichten [kg/m3]

std::pair<dcomp,dcomp> compute_kvert(double w, double k, double vp, double vs){
	dcomp hn = 2*pow(k,2)-(pow(w,2)/pow(vp,2));
	dcomp kn = 2*pow(k,2)-(pow(w,2)/pow(vs,2));
	return std::make_pair(sqrt(hn), sqrt(kn));
}

std::tuple<dcomp,dcomp,dcomp,dcomp,dcomp> compute_T(double w, double k, double vp, double vs, double mu, dcomp hn, dcomp kn){
	dcomp ln = pow(k,2)+pow(kn,2);
	dcomp t_factor = (-1.0)*pow(vs,2)/(2*mu*hn*kn*pow(w,2));
	//cout << ln << "\t" << t_factor << "\n";
	dcomp t11 = (2.0*i*mu*k*hn*kn)*t_factor;
	dcomp t12 = (mu*ln*kn)*t_factor;
	dcomp t13 = (hn*kn)*t_factor;
	dcomp t14 = (i*k*kn)*t_factor;
	dcomp t21 = ((-1.0)*mu*ln*kn)*t_factor;
	dcomp t22 = (2.0*i*mu*k*hn*kn)*t_factor;
	dcomp t23 = (i*k*hn)*t_factor;
	dcomp t24 = ((-1.0)*hn*kn)*t_factor;
	/*cout << t11 << "\t" << t12 << "\t" << t13 << "\t" << t14 << "\n"
		<< t21 << "\t" << t22 << "\t" << t23 << "\t" << t24 << "\n";*/
	dcomp T1212 = t11*t22 - t12*t21;
	dcomp T1213 = t11*t23 - t13*t21;
	dcomp T1214 = t11*t24 - t14*t21;
	dcomp T1224 = t12*t24 - t14*t22;
	dcomp T1234 = t13*t24 - t14*t23;
	return std::make_tuple(T1212,T1213,T1214,T1224,T1234);
}

std::tuple<dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp,dcomp> compute_G(dcomp hn, dcomp kn, double k, double dn, double w, double vp, double dens){
	
	dcomp hnn = hn/k;
	dcomp knn = kn/k;
	
	if (imag(hnn)==0){
		dcomp SH = 0.5*(1.0-exp((-2.0)*dn*hn))/hnn;
		dcomp CH = 0.5*(1.0+exp((-2.0)*dn*hn));
	}
	else{
		dcomp SH = sin((-1.0)*i*dn*hn)/hnn;
		dcomp CH = cos((-1.0)*i*dn*hn);
	}
	if (imag(knn)==0){
		dcomp SK = 0.5*(1.0-exp((-2.0)*dn*kn))/knn;
		dcomp CK = 0.5*(1.0+exp((-2.0)*dn*kn));
	}
	else{
		dcomp SK = sin((-1.0)*i*dn*kn)/knn;
		dcomp CK = cos((-1.0)*i*dn*kn);
	}
	
	double gam = 2.0*pow(k,2)/pow((w/vp),2);
	double a1 = pow(gam,2)-2.0*gam+1.0;
	dcomp a2 = pow(hnn,2)*pow(knn,2);
	double a3 = pow(gam,2)+a1;
	double a4 = 1.0-gam;
	dcomp a5 = pow(gam,2)*a2;
	dcomp expcorr = exp((-1.0)*hn*dn-kn*dn);
	double c1 = dens*pow(w,2)/k;
	double c2 = 1.0/c1;
	
	dcomp G1212 = a3*CH*CK-(a1+a5)*SH*SK-(a3-1)*expcorr;
	dcomp G1213 = c2*(CH*SK-pow(hnn,2)*SH*CK);
	dcomp iG1214 = i*c2*((a1-pow(gam,2))*(expcorr-CH*CK)+(a4-gam*a2)*SH*SK);
	dcomp G1224 = c2*(pow(knn,2)*CH*SK-SH*CK);
	dcomp G1234 = pow(c2,2)*(2.0*CH*CK+(1.0+a2)*SH*SK);
	dcomp G1312 = c1*(pow(gam,2)*pow(knn,2)*CH*SK-a1*SH*CK);
	dcomp iG1314 = i*(a4*SH*CK+gam*pow(knn,2)*CH*SK);
	dcomp G1324 = pow(knn,2)*SH*SK;
	dcomp iG1412 = i*c1*((a1-a4)*(a4-gam)*(expcorr-CH*CK)+(a4*a1-gam*a5)*SH*SK);
	dcomp iG1413 = i*(gam*pow(hnn,2)*SH*CH+a4*CH*SK);
	dcomp G1423 = CH*CK-G1212;
	dcomp G2412 = c1*(a1*CH*SK-pow(gam,2)*pow(hnn,2)*SH*CK);
	dcomp G2413 = pow(hnn,2)*SH*SK;
	dcomp G3412 = pow(c1,2)*(2*pow(gam,2)*a1*CH*CK+(pow(a1,2)+pow(gam,2)*a5)*SH*SK);
	
	return std::make_tuple(G1212,G1213,iG1214,G1224,G1234,G1312,iG1314,G1324,iG1412,iG1413,G1423,G2412,G2413,G3412);
}

int main()
{
	int nlay = depth.size();
	
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
		for(int kint=0; kint<=nk; kint++){
			
			double k=k_lim[0]+kint*(k_lim[1]-k_lim[0])/nk;
			cout << "Aktuelle Kreisfreq. & Wellenzahl: " << w[freq] << "\t" << k << "\n";
			
			for(int n=nlay-1;n>=0;n--){
				cout << "Schicht: " << n+1 << "\n";
				cout << "Geschwindigkeiten: " << vp[n] << "\t" << vs[n] << "\n";
				auto kvert = compute_kvert(w[freq], k, vp[n], vs[n]);
				cout << "kvert: " << std::get<0>(kvert) << "\t" << std::get<1>(kvert) << "\n\n\n";
				if (n==nlay-1){
					auto T = compute_T(w[freq], k, vp[n], vs[n], mu, std::get<0>(kvert), std::get<1>(kvert));
					cout << "T-Komponenten: " << std::get<0>(T) << "\t"
						<< std::get<1>(T) << "\t"
						<< std::get<2>(T) << "\t"
						<< std::get<3>(T) << "\t"
						<< std::get<4>(T) << "\n\n\n";
				}
				else{
					double dn = depth[n+1]-depth[n];
					cout << "Schichtdicke: " << dn << "m\n";
					auto G = compute_G(std::get<0>(kvert), std::get<1>(kvert),k,dn,w[freq],vp[n],dens[n]);
					cout << "Noch nicht programmiert!\n\n\n";
				}
			}
			
		}
	}
	
	return 0;
}


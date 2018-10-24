/*
 * 1Dvs2vph
 * Berechnet aus einem 1D Vs,Vp,Dichte-Modell die Phasengeschwindigkeit
 * zu gegebenen Perioden.
 * Quellen: Haskell (1953), Dunkin (1965)
 */

using namespace std;
#include <iostream>
#include <vector>
#include <math.h>
#include <complex>
#include <cmath>
#include <tuple>
typedef complex<double> dcomp;
dcomp i = sqrt(-1);

// Perioden in sec
std::vector<double> periods = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};

// Kreiswellenzahlen-Limits in rad/m
std::vector<double> k_lim = {2*M_PI/1000,2*M_PI/10};
int nk = 10; //Anzahl Wellenzahlen zum durchprobieren

// Definition 1D Modell
std::vector<double> depth = {0,25};	// Tiefe Schichtgrenzen [m]
std::vector<double> vs = {250,1000};	// Vs für Schichten [m/s]
std::vector<double> vp = {1350,2000};	// Vp für Schichten [m/s]
std::vector<double> dens = {2400,2400}; // Dichten [kg/m3]

std::tuple<dcomp,double,double,dcomp,double,dcomp,dcomp,double> compute_tn(double w, double k, double vp, double vs, double mu){
	double hn = sqrt(2*pow(k,2)-(pow(w,2)/pow(vp,2)));
	double kn = sqrt(2*pow(k,2)-(pow(w,2)/pow(vs,2)));
	double ln = pow(k,2)*+pow(kn,2);
	double t_factor = (-1)*pow(vs,2)/(2*mu*hn*kn*pow(w,2));
	dcomp tn11 = (2.0*i*mu*k*hn*kn)*t_factor;
	double tn12 = (mu*ln*kn)*t_factor;
	double tn13 = (hn*kn)*t_factor;
	dcomp tn14 = (i*k*kn)*t_factor;
	double tn21 = ((-1)*mu*ln*kn)*t_factor;
	dcomp tn22 = (2.0*i*mu*k*hn*kn)*t_factor;
	dcomp tn23 = (i*k*hn)*t_factor;
	double tn24 = ((-1)*hn*kn)*t_factor;
	return std::make_tuple(tn11,tn12,tn13,tn14,tn21,tn22,tn23,tn24);
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
		
	cout << "Kreisfrequenzen:\n";
	std::vector<double> w;
	for(int n=0; n<periods.size(); n++){
		w.push_back(2*M_PI/periods[n]);
		cout << w[n] << " rad/s\n";
	}
	
	double mu = vs[nlay-1]*dens[nlay-1]; // Shear modulus unterste Schicht
	
	for(int freq=0; freq<w.size(); freq++){
		for(int kint=0; kint<=nk; kint++){
			double k=k_lim[0]+kint*(k_lim[1]-k_lim[0])/nk;
			auto tn = compute_tn(w[freq], k, vp[nlay-1], vs[nlay-1], mu);
			
		}
	}
	
	return 0;
}


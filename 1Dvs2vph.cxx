/* 1Dvs2vph
 * Calculates the phase delay times between stations for given
 * earthquakes via computation of phase velocity maps (dispersion
 * curves) from a density, vp & vs model.
 * Uses boost, netcdf ans geographiclib libraries.
 * Sources: Haskell (1953), Dunkin (1965), Wathelet (2005),
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
#include <GeographicLib/TransverseMercatorExact.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/GeodesicLine.hpp>
#include <GeographicLib/Rhumb.hpp>
#include <GeographicLib/Constants.hpp>

using namespace std;
using namespace netCDF;
using namespace GeographicLib;

typedef complex<double> dcomp;
const std::complex<double> i(0, 1.0);

bool verbose = 0; // set to 1 for more output
double tolerance = 0.01; // Tolerance for phase velocity [m/s]
double length_tolerance = 1.0; // Tolerance for grat circle vs loxodrome length [m]
double mode_skip_it = 2.0;	// Number of additional iterations to check for mode skipping & factor to increase precision

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
	double vrlast = vs*0.99;
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

std::tuple<dcomp,dcomp,double,double,double,double,double,dcomp,dcomp,double,double,double> compute_util(double w, double c, double vp, double vs, double dn, bool botlay){
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
	
	/*if (verbose==1){
		cout << "hv: " << hv << " kv: " << kv << "\n"
			<< "hvnorm: " << hvnorm << " kvnorm: " << kvnorm << "\n"
			<< "SH: " << SH << " SK: " << SK << " CH: " << CH << " CK: " << CK << "\n"
			<< "gamma: " << gam << " l: " << l << "\n\n";
	}*/
	
	return std::make_tuple(hv, kv, SH, CH, SK, CK, gam, hvnorm, kvnorm, l, mh, mk);
}

std::tuple<dcomp, dcomp, double, double, double, double, double, double, double, double> compute_util_grads(double w, double vs, double vp, double c, double thck, bool botlay){
	auto util = compute_util(w, c, vp, vs, thck, botlay);
	
	double SH = std::get<2>(util);
	double CH = std::get<3>(util);
	double SK = std::get<4>(util);
	double CK = std::get<5>(util);
	double mh = std::get<10>(util);
	double mk = std::get<11>(util);
	
	dcomp hv_h, kv_k;
	double mkk, mhh, k = w/c;
	if(c<vp){
		mhh = pow(w,2)/(mh*pow(vp,3));
		hv_h = mhh;
	}
	else{
		mhh = ((-1.0)*pow(w,2))/(mh*pow(vp,3));
		hv_h = mhh;
		if(botlay==0)
			hv_h = i*mhh;
	}
		
	if(c<vs){
		mkk = pow(w,2)/(mk*pow(vs,3));
		kv_k = mkk;
	}
	else{
		mkk = ((-1.0)*pow(w,2))/(mk*pow(vs,3));
		kv_k = mkk;
		if(botlay==0)
			kv_k = i*mkk;
	}
		
	double hvnorm_h = (2.0*pow(c,2))/pow(vp,3);
	double kvnorm_k = (2.0*pow(c,2))/pow(vs,3);
	double gam_k = (4.0*vs)/pow(c,2);
	double l_k = 2.0*pow(w,2)/pow(vs,3);
	double CHH = (w*c*thck*SH)/pow(vp,3);
	double CKK = (w*c*thck*SK)/pow(vs,3);
	double SHH = (mhh/mh) * (k*thck*CH - SH);
	double SKK = (mkk/mk) * (k*thck*CK - SK);
	
	return std::make_tuple(hv_h, kv_k, hvnorm_h, kvnorm_k, gam_k, l_k, CHH, CKK, SHH, SKK);	
}

std::tuple<double,double,double,double,double> compute_T(double w, double c, double vp, double vs, double mu){
	// computes layer matrix for bottom layer (an i denotes an imaginary subdeterminant e.g. iT1214)
	double k = w/c;
	auto util = compute_util(w, c, vp, vs, 99999.0, 1);
	dcomp hv = std::get<0>(util);
	dcomp kv = std::get<1>(util);
	double l = std::get<9>(util);
	
	dcomp fact = pow((-1.0)*pow(vs,2)/(2*mu*hv*kv*pow(w,2)),2);

	double T1212 = std::real(pow(mu,2)*kv*hv*(pow(l,2)-4.0*pow(k,2)*kv*hv)*fact);
	double T1213 = std::real(mu*pow(hv,2)*kv*(l-2.0*pow(k,2))*fact);
	double iT1214 = std::real(k*mu*hv*kv*(l-2.0*hv*kv)*fact);
	double T1224 = std::real(mu*hv*pow(kv,2)*(2.0*pow(k,2)-l)*fact);
	double T1234 = std::real(hv*kv*(pow(k,2)-hv*kv)*fact);
	
	/*if (verbose==1) {
		cout << "factor: " << fact << "\n"
			<< "T1212: " << T1212 << "\t"
			<< "T1213: " << T1213 << "\t"
			<< "iT1214: " << iT1214 << "\t"
			<< "T1224: " << T1224 << "\t"
			<< "T1234: " << T1234 << "\n\n";
	}*/
	return std::make_tuple(T1212,T1213,iT1214,T1224,T1234);
}

std::tuple<double,double,double,double,double> compute_T_vs(double w, double c, double vp, double vs, double mu){
	auto T = compute_T(w, c, vp, vs, mu);
	double T1212 = std::get<0>(T);
	double iT1214 = std::get<2>(T);
	
	auto util = compute_util(w, c, vp, vs, 99999.0, 1);
	dcomp hv = std::get<0>(util);
	dcomp kv = std::get<1>(util);
	double l = std::get<9>(util);
	
	auto util_grads = compute_util_grads(w, vs, vp, c, 99999.0, 1);
	dcomp kv_k = std::get<1>(util_grads);
	double l_k = std::get<5>(util_grads);
	
	double gT1212 = std::real(4.0*T1212/vs + pow(vs,4)*(2.0*l*l_k*kv - pow(l,2)*kv_k)/(4.0*pow(w,4)*hv*pow(kv,2)));
	double gT1213 = std::real(kv_k*pow(vs,2)/(4.0*mu*pow(w,2)*pow(kv,2)));
	double igT1214 = std::real(2.0*iT1214/vs + pow(vs,4)*(l_k*kv - l*kv_k)/(4.0*mu*pow(w,3)*c*hv*pow(kv,2)));
	double gT1224 = 0.0;
	double gT1234 = std::real((-1.0)*kv_k*pow(vs,4)/(4.0*pow(mu,2)*pow(w,2)*pow(c,2)*hv*pow(kv,2)));
	
	return std::make_tuple(gT1212, gT1213, igT1214, gT1224, gT1234);
}

std::tuple<double,double,double,double,double> compute_T_vp(double w, double c, double vp, double vs, double mu){
	auto util = compute_util(w, c, vp, vs, 99999.0, 1);
	dcomp hv = std::get<0>(util);
	dcomp kv = std::get<1>(util);
	double l = std::get<9>(util);
	
	auto util_grads = compute_util_grads(w, vs, vp, c, 99999.0, 1);
	dcomp hv_h = std::get<0>(util_grads);
	
	double gT1212 = std::real((-1.0) * pow(vs,4)*pow(l,2)*hv_h/(4.0*pow(w,4)*pow(hv,2)*kv));
	double gT1213 = 0.0;
	double igT1214 = std::real((-1.0)*pow(vs,4)*l*hv_h/(4.0*mu*pow(w,3)*c*pow(hv,2)*kv));
	double gT1224 = std::real((-1.0)*hv_h*pow(vs,2)/(4.0*mu*pow(w,2)*pow(hv,2)));
	double gT1234 = std::real((-1.0)*hv_h*pow(vs,4)/(4.0*pow(mu,2)*pow(w,2)*pow(c,2)*pow(hv,2)*kv));
	
	return std::make_tuple(gT1212, gT1213, igT1214, gT1224, gT1234);
}

std::tuple<double,double,double,double,double> compute_T_rho(double w, double c, double vp, double vs, double mu){
	auto T = compute_T(w, c, vp, vs, mu);
	double T1213 = std::get<1>(T);
	double iT1214 = std::get<2>(T);
	double T1224 = std::get<3>(T);
	double T1234 = std::get<4>(T);
	
	double gT1212 = 0.0;
	double gT1213 = (-1.0)*T1213*pow(vs,2)/mu;
	double igT1214 = (-1.0)*iT1214*pow(vs,2)/mu;
	double gT1224 = (-1.0)*T1224*pow(vs,2)/mu;
	double gT1234 = (-2.0)*T1234*pow(vs,2)/mu;
	
	return std::make_tuple(gT1212, gT1213, igT1214, gT1224, gT1234);
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
	
	/*if (verbose==1){
		cout << "G-Komponenten:\n"
			<< "G1212: " << G1212 << " G1213: " << G1213 << " iG1214: " << iG1214 << " G1224: " << G1224 << " G1234: " << G1234 << "\n"
			<< "G1312: " << G1312 << " G1313: " << G1313 << " iG1314: " << iG1314 << " G1324: " << G1324 << "\n"
			<< "iG1412: " << iG1412 << " iG1413: " << iG1413 << " G1414: " << G1414 << "\n"
			<< "G2412: " << G2412 << " G2413: " << G2413 << "\n"
			<< "G3412: " << G3412 << "\n\n";
	}*/
		
	return std::make_tuple(G1212,G1213,iG1214,G1224,G1234,G1312,G1313,iG1314,G1324,iG1412,iG1413,G1414,G2412,G2413,G3412);
}

std::tuple<double,double,double,double,double,double,double,double,double,double,double,double,double,double,double> compute_G_vs(double c, double dn, double w, double vp, double vs, double dens){
	auto kvert = compute_util(w, c, vp, vs, dn, 0);
	double SH = std::get<2>(kvert);
	double CH = std::get<3>(kvert);
	double SK = std::get<4>(kvert);
	double CK = std::get<5>(kvert);
	double gam = std::get<6>(kvert);
	dcomp hvnorm = std::get<7>(kvert);
	dcomp kvnorm = std::get<8>(kvert);
	
	auto util_grads = compute_util_grads(w, vs, vp, c, dn, 0);
	double kvnorm_k = std::get<3>(util_grads);
	double gam_k = std::get<4>(util_grads);
	double CKK = std::get<7>(util_grads);
	double SKK = std::get<9>(util_grads);
	
	double gG1212 = std::real(2.0*gam_k*(1.0-2.0*gam)*(1.0-CH*CK)+(2.0*pow(gam,2)-2.0*gam+1)*CH*CKK-(2.0*gam_k*(gam-1.0)+2.0*gam*gam_k*pow(hvnorm,2)*pow(kvnorm,2)+pow(gam,2)*pow(hvnorm,2)*kvnorm_k)*SH*SK-(pow((1.0-gam),2)+pow(gam,2)*pow(hvnorm,2)*pow(kvnorm,2))*SH*SKK);
	double gG1213 = std::real((1.0/(dens*w*c))*(CH*SKK-pow(hvnorm,2)*SH*CKK));
	double igG1214 = std::real((1.0/(dens*w*c))*((-2.0)*gam_k*(1.0-CH*CK)-(1.0-2.0*gam)*CH*CKK)+(1.0/(dens*w*c))*(((-1.0)*gam_k-gam_k*pow(hvnorm,2)*pow(kvnorm,2)-gam*pow(hvnorm,2)*kvnorm_k)*SH*SK+(1.0-gam-gam*pow(hvnorm,2)*pow(kvnorm,2))*SH*SKK));
	double gG1224 = std::real((1.0/(dens*w*c))*(kvnorm_k*CH*SK+pow(kvnorm,2)*CH*SKK-SH*CKK));
	double gG1234 = std::real((1.0/pow((dens*w*c),2))*(-2.0*CH*CKK+(1.0+pow(hvnorm,2)*pow(kvnorm,2))*SH*SKK+pow(hvnorm,2)*kvnorm_k*SH*SK));
	double gG1312 = std::real(dens*w*c*(2.0*gam*gam_k*pow(kvnorm,2)*CH*SK+pow(gam,2)*kvnorm_k*CH*SK+pow(gam,2)*pow(kvnorm,2)*CH*SKK)+dens*w*c*(2.0*gam_k*(1.0-gam)*SH*CK-pow((1.0-gam),2)*SH*CKK));
	double gG1313 = std::real(CH*CKK);
	double igG1314 = std::real((-1.0)*gam_k*SH*CK+(1.0-gam)*SH*CKK+gam_k*pow(kvnorm,2)*CH*SK+gam*kvnorm_k*CH*SK+gam*pow(kvnorm,2)*CH*SKK);
	double gG1324 = std::real((-1.0)*kvnorm_k*SH*SK-pow(kvnorm,2)*SH*SKK);
	double igG1412 = std::real(dens*w*c*(gam_k*(-6.0*pow(gam,2)+6.0*gam-1.0)*(1.0-CH*CK)-(pow(gam,2)-gam)*(1.0-2.0*gam)*CH*CKK)+dens*w*c*(-3.0*pow((1.0-gam),2)*gam_k-3.0*pow(gam,2)*gam_k*pow(hvnorm,2)*pow(kvnorm,2)-pow(gam,3)*pow(hvnorm,2)*kvnorm_k)*SH*SK+dens*w*c*(pow((1.0-gam),3)-pow(gam,3)*pow(hvnorm,2)*pow(kvnorm,2))*SH*SKK);
	double igG1413 = std::real((-1.0)*((1.0-gam)*CH*SKK-gam_k*CH*SK+gam_k*pow(hvnorm,2)*SH*CK+gam*pow(hvnorm,2)*SH*CKK));
	double gG1414 = std::real(2.0*gam_k*(2.0*gam-1.0)*(1.0-CH*CK)-2.0*(pow(gam,2)-gam)*CH*CKK+(2.0*(gam-1)*gam_k+2.0*gam*gam_k*pow(hvnorm,2)*pow(kvnorm,2)+pow(gam,2)*pow(hvnorm,2)*kvnorm_k)*SH*SK+(pow((1.0-gam),2)+pow(gam,2)*pow(hvnorm,2)*pow(kvnorm,2))*SH*SKK);
	double gG2412 = std::real(dens*w*c*(2.0*(gam-1.0)*gam_k*CH*SK+pow((1.0-gam),2)*CH*SKK-2.0*gam*gam_k*pow(hvnorm,2)*SH*CK-pow(gam,2)*pow(hvnorm,2)*SH*CKK));
	double gG2413 = std::real((-1.0)*pow(hvnorm,2)*SH*SKK);
	double gG3412 = std::real((-1.0)*pow((dens*c*w),2)*(4.0*gam*gam_k*(1.0+2.0*pow(gam,2)-3.0*gam)*(1.0-CH*CK)-2.0*pow(gam,2)*pow((1.0-gam),2)*CH*CKK)-pow((dens*c*w),2)*((pow((1.0-gam),4)+pow(gam,4)*pow(hvnorm,2)*pow(kvnorm,2))*SH*SKK)-pow((dens*c*w),2)*((-4.0*gam_k*pow((1.0-gam),3)+4.0*pow(gam,3)*gam_k*pow(hvnorm,2)*pow(kvnorm,2)+pow(gam,4)*pow(hvnorm,2)*kvnorm_k)*SH*SK));
	return std::make_tuple(gG1212,gG1213,igG1214,gG1224,gG1234,gG1312,gG1313,igG1314,gG1324,igG1412,igG1413,gG1414,gG2412,gG2413,gG3412);
}

std::tuple<double,double,double,double,double,double,double,double,double,double,double,double,double,double,double> compute_G_vp(double c, double dn, double w, double vp, double vs, double dens){
	auto kvert = compute_util(w, c, vp, vs, dn, 0);
	double SH = std::get<2>(kvert);
	double SK = std::get<4>(kvert);
	double CK = std::get<5>(kvert);
	double gam = std::get<6>(kvert);
	dcomp hvnorm = std::get<7>(kvert);
	dcomp kvnorm = std::get<8>(kvert);
	
	auto util_grads = compute_util_grads(w, vs, vp, c, dn, 0);
	double hvnorm_h = std::get<2>(util_grads);
	double CHH = std::get<6>(util_grads);
	double SHH = std::get<8>(util_grads);
	
	double gG1212 = std::real((2.0*pow(gam,2)-2.0*gam+1.0)*CHH*CK-pow(gam,2)*hvnorm_h*pow(kvnorm,2)*SH*SK-(pow((1.0-gam),2)+pow(gam,2)*pow(hvnorm,2)*pow(kvnorm,2))*SHH*SK);
	double gG1213 = std::real((1.0/(dens*w*c))*(CHH*SK-hvnorm_h*SH*CK-pow(hvnorm,2)*SHH*CK));
	double igG1214 = std::real((1.0/(dens*w*c))*((2.0*gam-1.0)*CHH*CK+(1.0-gam-gam*pow(hvnorm,2)*pow(kvnorm,2))*SHH*SK-gam*hvnorm_h*pow(kvnorm,2)*SH*SK));
	double gG1224 = std::real((1.0/(dens*w*c))*(pow(kvnorm,2)*CHH*SK-SHH*CK));
	double gG1234 = std::real(((-1.0)/pow((w*dens*c),2))*(-2.0*CHH*CK+(1.0+pow(kvnorm,2)*pow(hvnorm,2))*SHH*SK+hvnorm_h*pow(kvnorm,2)*SH*SK));
	double gG1312 = std::real(dens*w*c*(pow(gam,2)*pow(kvnorm,2)*CHH*SK-pow((1.0-gam),2)*SHH*CK));
	double gG1313 = std::real(CHH*CK);
	double igG1314 = std::real((1.0-gam)*SHH*CK+gam*pow(kvnorm,2)*CHH*SK);
	double gG1324 = std::real((-1.0)*pow(kvnorm,2)*SHH*SK);
	double igG1412 = std::real(dens*w*c*(-1.0*gam*(gam-1.0)*(1.0-2.0*gam)*CHH*CK-pow(gam,3)*hvnorm_h*pow(kvnorm,2)*SH*SK)+dens*w*c*((pow((1.0-gam),3)-pow(gam,3)*pow(hvnorm,2)*pow(kvnorm,2))*SHH*SK));
	double igG1413 = std::real((gam-1.0)*CHH*SK-gam*hvnorm_h*SH*CK-gam*pow(hvnorm,2)*SHH*CK);
	double gG1414 = std::real(2.0*gam*(1.0-gam)*CHH*CK+(pow((1.0-gam),2)+pow(gam,2)*pow(hvnorm,2)*pow(kvnorm,2))*SHH*SK+pow(gam,2)*hvnorm_h*pow(kvnorm,2)*SH*SK);
	double gG2412 = std::real(dens*w*c*(pow((1.0-gam),2)*CHH*SK-pow(gam,2)*hvnorm_h*SH*CK-pow(gam,2)*pow(hvnorm,2)*SHH*CK));
	double gG2413 = std::real((-1.0)*(hvnorm_h*SH+pow(hvnorm,2)*SHH)*SK);
	double gG3412 = std::real((-1.0)*pow((dens*c*w),2)*(-2.0*pow(gam,2)*pow((1-gam),2)*CHH*CK)-pow((dens*c*w),2)*((pow((1-gam),4)+pow(gam,4)*pow(hvnorm,2)*pow(kvnorm,2))*SHH*SK+pow(gam,4)*hvnorm_h*pow(kvnorm,2)*SH*SK));
	
	return std::make_tuple(gG1212,gG1213,igG1214,gG1224,gG1234,gG1312,gG1313,igG1314,gG1324,igG1412,igG1413,gG1414,gG2412,gG2413,gG3412);
}

std::tuple<double,double,double,double,double,double,double,double,double,double,double,double,double,double,double> compute_G_rho(double c, double dn, double w, double vp, double vs, double dens){
	auto G = compute_G(c,dn,w,vp,vs,dens);
	double G1213 = std::get<1>(G);
	double iG1214 = std::get<2>(G);
	double G1224 = std::get<3>(G);
	double G1234 = std::get<4>(G);
	double G1312 = std::get<5>(G);
	double iG1412 = std::get<9>(G);
	double G2412 = std::get<12>(G);
	double G3412 = std::get<14>(G);
	
	double gG1212 = 0.0;
	double gG1213 = std::real((-1.0)*G1213/dens);
	double igG1214 = std::real((-1.0)*iG1214/dens);
	double gG1224 = std::real((-1.0)*G1224/dens);
	double gG1234 = std::real((-2.0)*G1234/dens);
	double gG1312 = std::real(G1312/dens);
	double gG1313 = 0.0;
	double igG1314 = 0.0;
	double gG1324 = 0.0;
	double igG1412 = std::real(iG1412/dens);
	double igG1413 = 0.0;
	double gG1414 = 0.0;
	double gG2412 = std::real(G2412/dens);
	double gG2413 = 0.0;
	double gG3412 = std::real(2.0*G3412/dens);
	return std::make_tuple(gG1212,gG1213,igG1214,gG1224,gG1234,gG1312,gG1313,igG1314,gG1324,igG1412,igG1413,gG1414,gG2412,gG2413,gG3412);
}
std::tuple<double,double,double,double,double> compute_R(double w, double c, double vp, double vs, double dn, double dens, std::tuple<double,double,double,double,double> T, int param){
	// Recursive layer stacking from bottom to top layer (an i denotes an imaginary subdeterminant e.g. iR1214)
	double T1212 = std::get<0>(T);
	double T1213 = std::get<1>(T);
	double iT1214 = std::get<2>(T);
	double T1224 = std::get<3>(T);
	double T1234 = std::get<4>(T);
	
	std::tuple<double,double,double,double,double,double,double,double,double,double,double,double,double,double,double> G;
	if(param==0){
		G = compute_G(c,dn,w,vp,vs,dens);
	}
	else if(param==1){
		G = compute_G_vs(c,dn,w,vp,vs,dens);
	}
	else if(param==2){
		G = compute_G_vp(c,dn,w,vp,vs,dens);
	}
	else if(param==3){
		G = compute_G_rho(c,dn,w,vp,vs,dens);
	}
	
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
	
	if(param==0){
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
	}
	
	/*if (verbose==1) {
		cout << "R1212: " << R1212 << "\n"
		<< "R1213: " << R1213 << "\n"
		<< "iR1214: " << iR1214 << "\n"
		<< "R1224: " << R1224 << "\n"
		<< "R1234: " << R1234 << "\n\n";
	}*/
	
	return std::make_tuple(R1212,R1213,iR1214,R1224,R1234);
}

double compute_R1212(double w, double c, std::vector<double> vp, std::vector<double> vs, double mu, std::vector<double> depth, std::vector<double> dens, int nlay, int param, int gradlay){
	// Recursive layer stacking from bottom to top to get R1212
	std::tuple<double,double,double,double,double> R;
	for(int n=nlay-1;n>=0;n--){
		/*if (verbose==1){
			cout << "Schicht: " << n+1 << "\n";
			cout << "Geschwindigkeiten: " << vp[n] << "\t" << vs[n] << "\n";
		}*/
		if (n==nlay-1){
			if(n==gradlay){
				if(param==1)
					R = compute_T_vs(w, c, vp[n], vs[n], mu);
				if(param==2)
					R = compute_T_vp(w, c, vp[n], vs[n], mu);
				else if(param==3)
					R = compute_T_rho(w, c, vp[n], vs[n], mu);
				}
			else
				R = compute_T(w, c, vp[n], vs[n], mu);
		}
		else {
			double dn = depth[n+1]-depth[n];
			/*if (verbose==1)
				cout << "Schichtdicke: " << dn << "m\n";*/
			if(n==gradlay)
				R = compute_R(w, c, vp[n], vs[n], dn, dens[n], R, param);
			else
				R = compute_R(w, c, vp[n], vs[n], dn, dens[n], R, 0);
		}
	}
	return(std::get<0>(R));
}

class R1212_root{
	// Root finding functor for R1212. Keeps all other varibles constants and only changes c (phase velocity) to finds roots. Uses the TOMS algorithm from the boost libraries
	public:
		R1212_root(double w_, std::vector<double> vp_, std::vector<double> vs_, double mu_, std::vector<double> depth_, std::vector<double> dens_, int nlay_):w(w_),vp(vp_),vs(vs_),mu(mu_),depth(depth_),dens(dens_),nlay(nlay_) {};
		double operator()(const double c){
			return compute_R1212(w, c, vp, vs, mu, depth, dens, nlay, 0, -999);
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

std::vector<vector<double>> get_gc_segments(double east0, double north0, double east1, double north1, double lon_centr) {
	// Approximates great circle path with loxodrome segments.
	double false_east = 500000.0; // false eating for utm coordinates
	
	// Set up of projections
	TransverseMercatorExact proj(Constants::WGS84_a(), Constants::WGS84_f(), Constants::UTM_k0());
	Geodesic geod(Constants::WGS84_a(), Constants::WGS84_f());
	Rhumb rhumb(Constants::WGS84_a(), Constants::WGS84_f());
	
	// Geographic coordinates for start & end point
	double lon0, lat0, lon1, lat1;
	proj.Reverse(lon_centr, east0-false_east, north0, lat0, lon0);
	proj.Reverse(lon_centr, east1-false_east, north1, lat1, lon1);
	std::vector<double> lon;
	std::vector<double> lat;
	std::vector<double> easting;
	std::vector<double> northing;
	lon.push_back(lon0);
	lon.push_back(lon1);
	lat.push_back(lat0);
	lat.push_back(lat1);
	easting.push_back(east0);
	easting.push_back(east1);
	northing.push_back(north0);
	northing.push_back(north1);
	
	// Great circle between start & end point
	const GeographicLib::GeodesicLine line = geod.InverseLine(lat[0], lon[0], lat[1], lon[1]);
	
	// Loxodrome distance between start & end point
	double azi_tmp, dist_lox;
	rhumb.Inverse(lat[0], lon[0], lat[1], lon[1], dist_lox, azi_tmp);
	
	// start with two points, check if distance difference is already smaller than tolerance
	int npts = 2;
	while(dist_lox-line.Distance() > length_tolerance) {
		// In each iteration keep only the starting point
		lon.erase(lon.begin()+1, lon.end());
		lat.erase(lat.begin()+1, lat.end());
		easting.erase(easting.begin()+1, easting.end());
		northing.erase(northing.begin()+1, northing.end());
		
		// calculate segment length
		double segment_length = line.Distance()/(npts+1);
		for(int pt = 1; pt <= npts; pt++) {
			// calculate point along great circle
			double lon_tmp, lat_tmp;
			line.Position(pt * segment_length, lat_tmp, lon_tmp);
			lat.push_back(lat_tmp);
			lon.push_back(lon_tmp);
			
			// transform to utm
			double x, y;
			proj.Forward(lon_centr, lat[pt], lon[pt], x, y);
			easting.push_back(x+false_east);
			northing.push_back(y);
		}
		// add end point to vectors
		lon.push_back(lon1);
		lat.push_back(lat1);
		easting.push_back(east1);
		northing.push_back(north1);
		
		//caculate loxodrome distance
		dist_lox = 0.0;
		double dist_tmp;
		for(int segment = 0; segment<easting.size()-1; segment++){
			rhumb.Inverse(lat[segment], lon[segment], lat[segment + 1], lon[segment + 1], dist_tmp, azi_tmp);
			dist_lox = dist_lox + dist_tmp;
		}
		
		// increase umber of points/segments
		npts = npts * 2;
	}
	
	// write eastings/northings to pts vector
	std::vector<vector<double>> pts;
	pts.push_back(easting);
	pts.push_back(northing);
	return pts;
}

/*std::vector<vector<double>>*/double get_t_segments(double east0, double north0, double east1, double north1, double event_e, double event_n, double lon_centr, std::vector<double> origin, double deast, double dnorth, std::vector<double> c, int ncells_north, int freq, int nperiods){//, std::vector<double> dcdvs, std::vector<double> dcdvp, std::vector<double> dcdrho, int nlay){
	// Computes relative phase delay for a station pair and a given earthquake
	double false_east = 500000.0; // false eating for utm coordinates
	origin[0] = origin[0] - deast;
	origin[1] = origin[1] - dnorth;
	
	// Set up coordinate transformations
	Geodesic geod(Constants::WGS84_a(), Constants::WGS84_f());
	TransverseMercatorExact proj(Constants::WGS84_a(), Constants::WGS84_f(), Constants::UTM_k0());
	double event_lat, event_lon;
	proj.Reverse(lon_centr, event_e-false_east, event_n, event_lat, event_lon);
	
	// Path segment is interpreted as line segment (not loxodrome)
	// because we don't care about its length (only orientation).
	// Check which model cells are crossed by inter-station path segment.
	double slope = (north1-north0)/(east1-east0);
	double intercept = north0-slope*east0;
	double ecell0 = floor((east0-origin[0])/deast);
	double ncell0 = floor((north0-origin[1])/dnorth);
	double ecell1 = floor((east1-origin[0])/deast);
	double ncell1 = floor((north1-origin[1])/dnorth);
	//std::vector<vector<double>> times_grads;
	//std::vector<double> times;
	double time = 0.0;
	double mid_e, mid_n, mid_lon, mid_lat, s12, az1, az2, se, sn, dist_segment_e, dist_segment_n;
	//std::vector<double> ecells,ncells,densgrad, vsgrad, vpgrad;
	//ecells.push_back(ecell0);
	//ncells.push_back(ncell0);
	
	while((abs(ecell0 - ecell1) > 0.0) | (abs(ncell0 - ncell1) > 0.0)) {
		double north_intercept, east_intercept, estep, nstep;
		if(east0 <= east1){
			east_intercept = origin[0]+(ecell0+1)*deast;
			estep = 1.0;
			if(slope < 0){
				north_intercept = ((origin[1]+ncell0*dnorth)-intercept)/slope;
				nstep = -1.0;
			}
			else {
				north_intercept = ((origin[1]+(ncell0+1)*dnorth)-intercept)/slope;
				nstep = 1.0;
			}
		}
		else{
			east_intercept = origin[0]+ecell0*deast;
			estep = -1.0;
			if(slope < 0){
				north_intercept = ((origin[1]+(ncell0+1)*dnorth)-intercept)/slope;
				nstep = 1.0;
			}
			else {
				north_intercept = ((origin[1]+ncell0*dnorth)-intercept)/slope;
				nstep = -1.0;
			}
		}
		double east_dist = abs(east0-east_intercept);
		double north_dist = abs(east0 - north_intercept);
		if(north_dist < east_dist){
			dist_segment_e = north_intercept-east0;
			dist_segment_n = (slope*north_intercept + intercept) - north0;
			mid_e = east0 + dist_segment_e/2.0;
			mid_n = north0 + dist_segment_n/2.0;
			north0 = slope*north_intercept + intercept;
			east0 = north_intercept;
			ncell0 = ncell0 + nstep;
		}
		else {
			dist_segment_e = east_intercept-east0;
			dist_segment_n = (slope*east_intercept+intercept)-north0;
			mid_e = east0 + dist_segment_e/2.0;
			mid_n = north0 + dist_segment_n/2.0;
			east0 = east_intercept;
			north0 = slope*east_intercept+intercept;
			ecell0 = ecell0 + estep;
		}
		//ecells.push_back(ecell0);
		//ncells.push_back(ncell0);
		proj.Reverse(lon_centr, mid_e-false_east, mid_n, mid_lat, mid_lon);
		geod.Inverse(event_lat, event_lon, mid_lat, mid_lon, s12, az1, az2);
		se = cos((90.0-az2)*M_PI/180.0)*(1.0/c[ecell0 * ncells_north + ncell0]);
		sn = sin((90.0-az2)*M_PI/180.0)*(1.0/c[ecell0 * ncells_north + ncell0]);
		time = time + (se * dist_segment_e + sn * dist_segment_n)*(-1.0);
		/*for(int n=0;n<nlay;n++){
			double tmp = dcdvs[ecell0 * ncells_north * nperiods * nlay + ncell0 * nperiods * nlay + freq*nlay + n];
			tmp = sqrt(2)/tmp;
			tmp = tmp*dist_segment_e+tmp*dist_segment_n;
			vsgrad.push_back(tmp);
			tmp = dcdvp[ecell0 * ncells_north * nperiods * nlay + ncell0 * nperiods * nlay + freq*nlay + n];
			tmp = sqrt(2)/tmp;
			tmp = tmp*dist_segment_e+tmp*dist_segment_n;
			vpgrad.push_back(tmp);
			tmp = dcdrho[ecell0 * ncells_north * nperiods * nlay + ncell0 * nperiods * nlay + freq*nlay + n];
			tmp = sqrt(2)/tmp;
			tmp = tmp*dist_segment_e+tmp*dist_segment_n;
			densgrad.push_back(tmp);
		}*/
	}
	//ecells.push_back(ecell1);
	//ncells.push_back(ncell1);
	dist_segment_e = east1-east0;
	dist_segment_n = north1-north0;
	mid_e = east0 + dist_segment_e/2.0;
	mid_n = north0 + dist_segment_n/2.0;
	proj.Reverse(lon_centr, mid_e-false_east, mid_n, mid_lat, mid_lon);
	geod.Inverse(event_lat, event_lon, mid_lat, mid_lon, s12, az1, az2);
	se = cos((90-az2)*M_PI/180)*(1/c[ecell1 * ncells_north + ncell1]);
	sn = sin((90-az2)*M_PI/180)*(1/c[ecell1 * ncells_north + ncell1]);
	time = time + (se * dist_segment_e + sn * dist_segment_n)*(-1.0);
	/*for(int n=0;n<nlay;n++){
		double tmp = dcdvs[ecell0 * ncells_north * nperiods * nlay + ncell0 * nperiods * nlay + freq*nlay + n];
		tmp = sqrt(2)/tmp;
		tmp = tmp*dist_segment_e+tmp*dist_segment_n;
		vsgrad.push_back(tmp);
		tmp = dcdvp[ecell0 * ncells_north * nperiods * nlay + ncell0 * nperiods * nlay + freq*nlay + n];
		tmp = sqrt(2)/tmp;
		tmp = tmp*dist_segment_e+tmp*dist_segment_n;
		vpgrad.push_back(tmp);
		tmp = dcdrho[ecell0 * ncells_north * nperiods * nlay + ncell0 * nperiods * nlay + freq*nlay + n];
		tmp = sqrt(2)/tmp;
		tmp = tmp*dist_segment_e+tmp*dist_segment_n;
		densgrad.push_back(tmp);
	}
	times.push_back(time);
	times_grads.push_back(times);
	times_grads.push_back(ecells);
	times_grads.push_back(ncells);
	times_grads.push_back(vsgrad);
	times_grads.push_back(vpgrad);
	times_grads.push_back(densgrad);*/
	return time; //times_grads;
}

int main(){
	
	// Read phase delay time observations
	NcFile dtpFile("./dt_usarray_win_neu_utm.nc", NcFile::read);
	NcDim nperiodsIn = dtpFile.getDim("NumberOfPeriods");
	NcDim nstatsIn = dtpFile.getDim("NumberOfStations");
	NcDim nsrcsIn = dtpFile.getDim("NumberOfRays");
	NcDim neventsIn = dtpFile.getDim("NumberOfEvents");
	NcDim nevents_per_srcIn = dtpFile.getDim("EventsPerSRC");
	int nperiods = nperiodsIn.getSize();
	int nstats = nstatsIn.getSize();
	int nsrcs = nsrcsIn.getSize();
	int nevents = neventsIn.getSize();
	int nevents_per_src = nevents_per_srcIn.getSize();
	
	std::vector<double> periods(nperiods);
	std::vector<double> mpn(nstats);
	std::vector<double> mpe(nstats);
	std::vector<double> mpz(nstats);
	std::vector<double> src_rcvr_cmb(nsrcs*2.0);
	std::vector<double> dtp(nsrcs*nevents_per_src*nperiods);
	std::vector<double> event_stat_cmb(nsrcs*nevents_per_src);
	std::vector<double> eventx(nevents);
	std::vector<double> eventy(nevents);
	
	NcVar periodsIn=dtpFile.getVar("Periods");
	periodsIn.getVar(periods.data());
	NcVar mpnIn=dtpFile.getVar("MeasPosX");
	mpnIn.getVar(mpn.data());
	NcVar mpeIn=dtpFile.getVar("MeasPosY");
	mpeIn.getVar(mpe.data());
	NcVar mpzIn=dtpFile.getVar("MeasPosZ");
	mpzIn.getVar(mpz.data());
	NcVar src_rcvr_cmbIn=dtpFile.getVar("StatComb");
	src_rcvr_cmbIn.getVar(src_rcvr_cmb.data());
	NcVar dtpIn=dtpFile.getVar("dtp");
	dtpIn.getVar(dtp.data());
	NcVar event_stat_cmbIn=dtpFile.getVar("EventStatComb");
	event_stat_cmbIn.getVar(event_stat_cmb.data());
	NcVar eventxIn=dtpFile.getVar("EventPosX");
	eventxIn.getVar(eventx.data());
	NcVar eventyIn=dtpFile.getVar("EventPosY");
	eventyIn.getVar(eventy.data());
	
	double dtp_dummy;
	double *dummy_pointer = &dtp_dummy;
	NcVarAtt dummy = dtpIn.getAtt("_FillValue");
	dummy.getValues(dummy_pointer);
	
	double lon_centr;
	double *lon_centr_pnt = &lon_centr;
	NcVarAtt lonc = mpnIn.getAtt("Central_meridian");
	lonc.getValues(lon_centr_pnt);
	
	// Conversion of periods to angular frequencies
	std::vector<double> w;
	for(int n=0; n<nperiods; n++){
		w.push_back(2.0*M_PI/periods[n]);
		if (verbose==1)
			cout << "Kreisfreq. " << n << ": " << w[n] << " rad/s\n";
	}
	
	// Read density, vs, vp from nc file
	NcFile densFile("./dens_na_neu_utm.nc", NcFile::read);
	NcFile vpFile("./vp_na_neu_utm.nc", NcFile::read);
	NcFile vsFile("./vs_na_neu_utm.nc", NcFile::read);
	
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
	
	// Shift vector of depths so that it starts from 0
	int nlay = depth.size();	// number of layers
	for(int n=nlay-1; n>=0; n--)
		depth[n] = depth[n] - depth[0];
	
	// get model cell sizes (east and northwards)
	double model_cell_east = east[1] - east[0];
	double model_cell_north = north[1] - north[0];
	std::vector<double> model_origin = {east[0], north[0]};
	
	// Open output files, write header lines
	ofstream resultfile;
	resultfile.open ("dispersion.out");
	resultfile << "# Easting [m] \t Northing [m] \t Period [s] \t Phase velocity [m/s] \t Difference [m/s] \t No. of iterations";
	
	ofstream delayfile;
	delayfile.open ("delays.out");
	delayfile << "#Easting_epi [m] \t Northing_epi [m] \t Event_num \t Easting_stat1 [m] \t Northing_stat1 [m] \t Easting_stat2 [m] \t Northing_stat2 [m] \t stat1_num \t stat2_num \t Period [s] \t Phase delay [s]";
	
	for(int freq=0; freq<nperiods; freq++){
		//Vectors to store gradients, dispersion curves
		std::vector<double> /*dcdrho, dcdvs, dcdvp,*/ dispersion;	
		for (int estep = 0; estep<NY; estep++){
			for (int nstep = 0; nstep<NX; nstep++){
				std::vector<double> dens;
				std::vector<double> vs;
				std::vector<double> vp;	
				bool lvz=0;
				for (int n=0; n<nlay; n++){
					// sort velocities, densities into 1D models 
					dens.push_back(dens_all[n+nlay*estep+NY*nlay*nstep]);
					vs.push_back(vs_all[n+nlay*estep+NY*nlay*nstep]);
					// check if there's a low velocity zone
					if (n>0 & vs[n]<vs[n-1])
						lvz=1;
					//cout << lvz << "\n";
					vp.push_back(vp_all[n+nlay*estep+NY*nlay*nstep]);
				}
			
				if (vs[0]<=0){
					// This still needs some work. Computation of dispersion curves does not work if top layer has vs=0
					dispersion.push_back(0.0);
					resultfile << "\n" << east[estep] << "\t" << north[nstep] << "\t" << (2.0*M_PI)/w[freq] << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0;
					continue;
				}
				else{				
					// Calculation of velocity limits 
					std::vector<double> c_lim={0,0};
					double vsmin = *std::min_element(vs.begin(),vs.end());
					double vsmax = *std::max_element(vs.begin(),vs.end());
					double vpmin = *std::min_element(vp.begin(),vp.end());
					double vpmax = *std::max_element(vp.begin(),vp.end());
					c_lim[0] = newton_vr(vpmin, vsmin)/1.05;
					c_lim[1] = newton_vr(vpmax, vsmax)*1.05;
				
					// step ratio for root bracketing
					double stepratio = (vsmin - c_lim[0])/(2.0*vsmin);	
	
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
					}
				
					// Shear modulus bottom layer
					double mu = pow(vs[nlay-1],2)*dens[nlay-1];
					if (verbose==1)
						cout << "Schermodul unterste Schicht: " << mu << "\n";
				
					// Compute initial R1212 polarization for large period below fundamental mode
					double R1212 = compute_R1212(w[nperiods]/10.0, c_lim[0], vp, vs, mu, depth, dens, nlay, 0, -999);
					bool pol0 = signbit(R1212);
				
					if(verbose==1)
						cout << "Polarisation von R1212 für geringe Geschwindigkeit: " << pol0 << "\n";
				
					double c_last=c_lim[0];	//
				
					// Loop over all periods/frequencies
					//for(int freq=0; freq<nperiods; freq++){
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
			
						if (verbose==1){
							cout << "c0: " << c0 << "\t" << "stepratio: " << stepratio << "\t" << "precision: " << precision << "\n";
							cout << "Aktuelle Kreisfreq. & Geschwindigkeit: " << w[freq] << "\t" << c1 << "\n";
						}
						
						// Check polarization of R1212 for the upper bracket	
						R1212 = compute_R1212(w[freq], c1, vp, vs, mu, depth, dens, nlay, 0, -999);
						pol1 = signbit(R1212);
						
						// If a sign change is found check for mode skipping
						if (pol0!=pol1 & (c1-c0)>(2.0*tolerance)){
							double c2 = (c1+c0)/2.0;	// set new speed between brackets
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
								R1212 = compute_R1212(w[freq], c2, vp, vs, mu, depth, dens, nlay, 0, -999);
								bool pol2 = signbit(R1212);
								// if mode skipping detected increase precision (-> decrease step ratio) and return to bracket search
								if (pol2==pol1){
									precision = precision * mode_skip_it;
									c1 = c0;
									if (verbose==1){
										cout << "Error: Mode skipping detected!\n";
									}
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
					if (lvz == 0 & (brackets.first + brackets.second)/2.0 < c_last) {
						c1 = c_lim[0];
						precision = precision * mode_skip_it;
						if (verbose==1)
							cout << "Error: non-monotonous shape of dispersion curve!\n";
						goto cnt;
					}
					c_last = (brackets.first + brackets.second)/2.0;
					
					dispersion.push_back(c_last);
							
					// Write output to file
					resultfile << "\n" << east[estep] << "\t" << north[nstep] << "\t" << (2.0*M_PI)/w[freq] << "\t" << c_last << "\t" << brackets.second-brackets.first << "\t" << max_iter;
					
					
					/*for(int n=0;n<nlay;n++){
						//Computation of Gradients
						double R_tmp = compute_R1212(w[freq], c_last, vp, vs, mu, depth, dens, nlay, 1, n);
						dcdvs.push_back(R_tmp);
						R_tmp = compute_R1212(w[freq], c_last, vp, vs, mu, depth, dens, nlay, 2, n);
						dcdvp.push_back(R_tmp);
						R_tmp = compute_R1212(w[freq], c_last, vp, vs, mu, depth, dens, nlay, 3, n);
						dcdrho.push_back(R_tmp);
					}*/
					
					
				//} // end loop over frequencies
				} 
			} // end of loop over northing
		} // end of loop over easting
	
		// close file
		resultfile << "\n";
		resultfile.close();
	
		/*ofstream gradfile;
		gradfile.open ("gradients.out");
		gradfile << "#Source receiver comb. \t period[s] \t ecell \t ncell \n #gradvs \t gradvp \t graddens";
		*/
	
		// loop over all rays, computes phase delays
		for (int src=0; src<nsrcs; src++){
			//cout << "src: \t" << src << "\n";
			std::vector<vector<double>> segments;
			segments = get_gc_segments(mpe[src_rcvr_cmb[src]], mpn[src_rcvr_cmb[src]], mpe[src_rcvr_cmb[src+nsrcs]], mpn[src_rcvr_cmb[src+nsrcs]], lon_centr);
			std::vector<double> seg_east = segments[0];
			std::vector<double> seg_north = segments[1];
			//for (int freq=0; freq<nperiods; freq++){
			//cout << "freq: \t" << freq << "\n";
			for(int event=0; event<nevents_per_src; event++){
				//cout << "event: \t" << event << "\n";
				if (dtp[freq*nsrcs*nevents_per_src+event*nsrcs+src]==dtp_dummy | event_stat_cmb[event*nsrcs+src]==dtp_dummy){
					// if there is only a dummy value we can skip this period
					delayfile << "\n" << dtp_dummy << "\t" << dtp_dummy << "\t" << event_stat_cmb[event*nsrcs+src] << "\t" << mpe[src_rcvr_cmb[src]] << "\t" << mpn[src_rcvr_cmb[src]] << "\t" << mpe[src_rcvr_cmb[src+nsrcs]] << "\t" << mpn[src_rcvr_cmb[src+nsrcs]] << "\t" << src_rcvr_cmb[src] << "\t" << src_rcvr_cmb[src+nsrcs] << "\t" << (2.0*M_PI)/w[freq] << "\t" << dtp_dummy;
					continue;
				}
				else {
					// loop over segments
					double time_total = 0.0;
					//std::vector<double> ecells, ncells, tgradvs, tgradvp, tgradrho;
					for(int seg=0; seg<seg_east.size()-1; seg++){
						/*std::vector<vector<double>>*/double time_segment = get_t_segments(seg_east[seg], seg_north[seg], seg_east[seg+1], seg_north[seg+1], eventy[event_stat_cmb[event*nsrcs+src]], eventx[event_stat_cmb[event*nsrcs+src]], lon_centr, model_origin, model_cell_east, model_cell_north, dispersion, NX, freq, nperiods);//, dcdvs, dcdvp, dcdrho, nlay);
						//std::vector<double> ts = time_segment[0];
						time_total = time_total + time_segment;//ts[0];
						/*std::vector<double> ectmp = time_segment[1];
						std::vector<double> nctmp = time_segment[2];
						for(int cell = 0; cell<ectmp.size(); cell++){
							ecells.push_back(ectmp[cell]);
							ncells.push_back(nctmp[cell]);
						}
						std::vector<double> vstmp = time_segment[3];
						std::vector<double> vptmp = time_segment[4];
						std::vector<double> denstmp = time_segment[5];
						for(int ngrads=0; ngrads < vstmp.size(); ngrads++){
							tgradvs.push_back(vstmp[ngrads]);
							tgradvp.push_back(vptmp[ngrads]);
							tgradrho.push_back(denstmp[ngrads]);
						}*/
					}
					delayfile << "\n" << eventy[event_stat_cmb[event*nsrcs+src]] << "\t" << eventx[event_stat_cmb[event*nsrcs+src]] << "\t" << event_stat_cmb[event*nsrcs+src] << "\t" << mpe[src_rcvr_cmb[src]] << "\t" << mpn[src_rcvr_cmb[src]] << "\t" << mpe[src_rcvr_cmb[src+nsrcs]] << "\t" << mpn[src_rcvr_cmb[src+nsrcs]] << "\t" << src_rcvr_cmb[src] << "\t" << src_rcvr_cmb[src+nsrcs] << "\t" << (2.0*M_PI)/w[freq] << "\t" << time_total;
					/*for(int cell=0; cell<ecells.size(); cell++){
						gradfile << "\n" << src << "\t" << (2.0*M_PI)/w[freq] << "\t" << ecells[cell] << "\t" << ncells[cell];
						for(int grad=0; grad<nlay; grad++){
							gradfile << "\n" << tgradvs[cell*nlay+grad] << "\t" << tgradvp[cell*nlay+grad] << "\t" << tgradrho[cell*nlay+grad];
						}
					}*/
				}
			}	
		} // end loop rays
	} // end loop frequencies
	
	delayfile << "\n";
	delayfile.close();
	
	/*gradfile << "\n";
	gradfile.close();*/

	// end program
	return 0;
}

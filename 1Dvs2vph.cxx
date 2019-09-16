/* 1Dvs2vph
 * Calculates the phase delay times between stations for given
 * earthquakes via computation of phase velocity maps (dispersion
 * curves) from a density, vp & vs model.
 * Uses boost, netcdf ans geographiclib libraries.
 * Sources: Haskell (1953), Dunkin (1965), Wathelet (2005),
 * Cercato (2007)
 */

#include <string>
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <functional>
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

string datafile = "dt_tiny.nc";
string dens_file = "dens_tiny.nc";
string vs_file = "vs_tiny.nc";
string vp_file = "vp_tiny.nc";
const double false_east = 500000.0; // false eating for utm coordinates
const bool calcgrads = 1; // set to 1 to calculate gradients.
const bool verbose = 0; // set to 1 for more output
const double tolerance = 0.01; // Tolerance for phase velocity [m/s]
const double length_tolerance = 1.0; // Tolerance for grat circle vs loxodrome length [m]
const double mode_skip_it = 2.0;	// Number of additional iterations to check for mode skipping & factor to increase precision

std::vector<std::vector<double>> read_data(const string &datafile, double &dtp_dummy, double &lon_centr){
	
	int nperiods, nstats, nsrcs, nevents, nevents_per_src;
	
	// Read phase delay time observations
	NcFile dtpFile(datafile, NcFile::read);
	NcDim nperiodsIn = dtpFile.getDim("NumberOfPeriods");
	NcDim nstatsIn = dtpFile.getDim("NumberOfStations");
	NcDim nsrcsIn = dtpFile.getDim("NumberOfRays");
	NcDim neventsIn = dtpFile.getDim("NumberOfEvents");
	NcDim nevents_per_srcIn = dtpFile.getDim("EventsPerSRC");
	nperiods = nperiodsIn.getSize();
	nstats = nstatsIn.getSize();
	nsrcs = nsrcsIn.getSize();
	nevents = neventsIn.getSize();
	nevents_per_src = nevents_per_srcIn.getSize();
	
	std::vector<double> periods(nperiods);
	std::vector<double> mpn(nstats);
	std::vector<double> mpe(nstats);
	std::vector<double> mpz(nstats);
	std::vector<double> src_rcvr_cmb(nsrcs*2);
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
	
	NcVarAtt dummy = dtpIn.getAtt("_FillValue");
	dummy.getValues(&dtp_dummy);
	
	NcVarAtt lonc = mpnIn.getAtt("Central_meridian");
	lonc.getValues(&lon_centr);
	
	std::vector<std::vector<double>> dtp_data(9);
	dtp_data[0] = periods;
	dtp_data[1] = mpn;
	dtp_data[2] = mpe;
	dtp_data[3] = mpz;
	dtp_data[4] = src_rcvr_cmb;
	dtp_data[5] = dtp;
	dtp_data[6] = event_stat_cmb;
	dtp_data[7] = eventx;
	dtp_data[8] = eventy;
	
	return dtp_data;
}

std::vector<std::vector<double>> read_model(const string &dens_file, const string &vs_file, const string &vp_file){
	// Read density, vs, vp from nc file
	NcFile densFile(dens_file, NcFile::read);
	NcFile vpFile(vp_file, NcFile::read);
	NcFile vsFile(vs_file, NcFile::read);
	
	// Get dimensions of model
	int NX, NY, NZ;
	NcDim nxIn=densFile.getDim("Northing");
	NcDim nyIn=densFile.getDim("Easting");
	NcDim nzIn=densFile.getDim("Depth");
	NX = nxIn.getSize();
	NY = nyIn.getSize();
	NZ = nzIn.getSize();
	
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
	
	std::vector<std::vector<double>> model(6);
	model[0] = depth;
	model[1] = north;
	model[2] = east;
	model[3] = dens_all;
	model[4] = vp_all;
	model[5] = vs_all;
	
	return model;
}

// computation of Rayleigh velocity of homogenous half space using Newton Raphson method
double newton_vr(const double &vp, const double &vs){
	// variables to store rayleigh velocity
	double vrlast = vs*0.99;
	double vr;
	double diff=99999.0;
	// calculate Rayleigh velocity for homogenous half space
	while(diff>tolerance){ 
		double fvr = 4.0-4.0*(pow(vrlast,2)/pow(vs,2))+pow(vrlast,4)/pow(vs,4)-4.0*sqrt(1-pow(vrlast,2)/pow(vp,2))*sqrt(1.0-pow(vrlast,2)/pow(vs,2));
		double dfvr = -8.0*vrlast/pow(vs,2)+4.0*pow(vrlast,3)/pow(vs,4)+(4.0*vrlast*(pow(vp,2)+pow(vs,2)-2.0*pow(vrlast,2)))/(vp*vs*sqrt((vp-vrlast)*(vp+vrlast))*sqrt((vs-vrlast)*(vs+vrlast)));
		vr = vrlast - fvr/dfvr;
		diff = sqrt(pow(vr-vrlast,2));
		vrlast = vr;
	}
	if (verbose == 1)
		cout << "vr: " << vr << "\t final diff: " << diff << "\n\n";
	return vr;
}

std::tuple<dcomp,dcomp,double,double,double,double,double,dcomp,dcomp,double,double,double> compute_util(const double &w, const double &c, const double &vp, const double &vs, const double &dn, const bool &botlay){
	// computes some constants for each layer (vertical wave numbers and some derived properties)
	dcomp hv, kv;
	double SH, CH, SK, CK, mh, mk;
	const double k = w/c;
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
	const double gam = 2.0*pow(vs,2)/pow(c,2);
	const dcomp hvnorm = hv/k;
	const dcomp kvnorm = kv/k;
	const double l = 2.0*pow(k,2)-pow(w/vs,2);
	
	return std::make_tuple(hv, kv, SH, CH, SK, CK, gam, hvnorm, kvnorm, l, mh, mk);
}

std::tuple<dcomp, dcomp, double, double, double, double, double, double, double, double> compute_util_grads(const double &w, const double &vs, const double &vp, const double &c, const double &thck, const bool &botlay){
	const auto util = compute_util(w, c, vp, vs, thck, botlay);
	
	const double SH = std::get<2>(util);
	const double CH = std::get<3>(util);
	const double SK = std::get<4>(util);
	const double CK = std::get<5>(util);
	const double mh = std::get<10>(util);
	const double mk = std::get<11>(util);
	
	dcomp hv_h, kv_k;
	double mkk, mhh;
	const double k = w/c;
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
		
	const double hvnorm_h = (2.0*pow(c,2))/pow(vp,3);
	const double kvnorm_k = (2.0*pow(c,2))/pow(vs,3);
	const double gam_k = (4.0*vs)/pow(c,2);
	const double l_k = 2.0*pow(w,2)/pow(vs,3);
	const double CHH = (w*c*thck*SH)/pow(vp,3);
	const double CKK = (w*c*thck*SK)/pow(vs,3);
	const double SHH = (mhh/mh) * (k*thck*CH - SH);
	const double SKK = (mkk/mk) * (k*thck*CK - SK);
	
	return std::make_tuple(hv_h, kv_k, hvnorm_h, kvnorm_k, gam_k, l_k, CHH, CKK, SHH, SKK);	
}

std::tuple<double,double,double,double,double> compute_T(const double &w, const double &c, const double &vp, const double &vs, const double &mu){
	// computes layer matrix for bottom layer (an i denotes an imaginary subdeterminant e.g. iT1214)
	const double k = w/c;
	const auto util = compute_util(w, c, vp, vs, 99999.0, 1);
	const dcomp hv = std::get<0>(util);
	const dcomp kv = std::get<1>(util);
	const double l = std::get<9>(util);
	
	const dcomp fact = pow((-1.0)*pow(vs,2)/(2*mu*hv*kv*pow(w,2)),2);

	const double T1212 = std::real(pow(mu,2)*kv*hv*(pow(l,2)-4.0*pow(k,2)*kv*hv)*fact);
	const double T1213 = std::real(mu*pow(hv,2)*kv*(l-2.0*pow(k,2))*fact);
	const double iT1214 = std::real(k*mu*hv*kv*(l-2.0*hv*kv)*fact);
	const double T1224 = std::real(mu*hv*pow(kv,2)*(2.0*pow(k,2)-l)*fact);
	const double T1234 = std::real(hv*kv*(pow(k,2)-hv*kv)*fact);
	
	return std::make_tuple(T1212,T1213,iT1214,T1224,T1234);
}

std::tuple<double,double,double,double,double> compute_T_vs(const double &w, const double &c, const double &vp, const double &vs, const double &mu){
	const auto T = compute_T(w, c, vp, vs, mu);
	const double T1212 = std::get<0>(T);
	const double iT1214 = std::get<2>(T);
	
	const auto util = compute_util(w, c, vp, vs, 99999.0, 1);
	const dcomp hv = std::get<0>(util);
	const dcomp kv = std::get<1>(util);
	const double l = std::get<9>(util);
	
	const auto util_grads = compute_util_grads(w, vs, vp, c, 99999.0, 1);
	const dcomp kv_k = std::get<1>(util_grads);
	const double l_k = std::get<5>(util_grads);
	
	const double gT1212 = std::real(4.0*T1212/vs + pow(vs,4)*(2.0*l*l_k*kv - pow(l,2)*kv_k)/(4.0*pow(w,4)*hv*pow(kv,2)));
	const double gT1213 = std::real(kv_k*pow(vs,2)/(4.0*mu*pow(w,2)*pow(kv,2)));
	const double igT1214 = std::real(2.0*iT1214/vs + pow(vs,4)*(l_k*kv - l*kv_k)/(4.0*mu*pow(w,3)*c*hv*pow(kv,2)));
	const double gT1224 = 0.0;
	const double gT1234 = std::real((-1.0)*kv_k*pow(vs,4)/(4.0*pow(mu,2)*pow(w,2)*pow(c,2)*hv*pow(kv,2)));
	
	return std::make_tuple(gT1212, gT1213, igT1214, gT1224, gT1234);
}

std::tuple<double,double,double,double,double> compute_T_vp(const double &w, const double &c, const double &vp, const double &vs, const double &mu){
	const auto util = compute_util(w, c, vp, vs, 99999.0, 1);
	const dcomp hv = std::get<0>(util);
	const dcomp kv = std::get<1>(util);
	const double l = std::get<9>(util);
	
	const auto util_grads = compute_util_grads(w, vs, vp, c, 99999.0, 1);
	const dcomp hv_h = std::get<0>(util_grads);
	
	const double gT1212 = std::real((-1.0) * pow(vs,4)*pow(l,2)*hv_h/(4.0*pow(w,4)*pow(hv,2)*kv));
	const double gT1213 = 0.0;
	const double igT1214 = std::real((-1.0)*pow(vs,4)*l*hv_h/(4.0*mu*pow(w,3)*c*pow(hv,2)*kv));
	const double gT1224 = std::real((-1.0)*hv_h*pow(vs,2)/(4.0*mu*pow(w,2)*pow(hv,2)));
	const double gT1234 = std::real((-1.0)*hv_h*pow(vs,4)/(4.0*pow(mu,2)*pow(w,2)*pow(c,2)*pow(hv,2)*kv));
	
	return std::make_tuple(gT1212, gT1213, igT1214, gT1224, gT1234);
}

std::tuple<double,double,double,double,double> compute_T_rho(const double &w, const double &c, const double &vp, const double &vs, const double &mu){
	const auto T = compute_T(w, c, vp, vs, mu);
	const double T1213 = std::get<1>(T);
	const double iT1214 = std::get<2>(T);
	const double T1224 = std::get<3>(T);
	const double T1234 = std::get<4>(T);
	
	const double gT1212 = 0.0;
	const double gT1213 = (-1.0)*T1213*pow(vs,2)/mu;
	const double igT1214 = (-1.0)*iT1214*pow(vs,2)/mu;
	const double gT1224 = (-1.0)*T1224*pow(vs,2)/mu;
	const double gT1234 = (-2.0)*T1234*pow(vs,2)/mu;
	
	return std::make_tuple(gT1212, gT1213, igT1214, gT1224, gT1234);
}

std::tuple<double,double,double,double,double,double,double,double,double,double,double,double,double,double,double> compute_G(const double &c, const double &dn, const double &w, const double &vp, const double &vs, const double &dens){
	// computes subdeterminants of G matrix (an i denotes an imaginary subdeterminant e.g. iG1214)
	const auto kvert = compute_util(w, c, vp, vs, dn, 0);
	const double SH = std::get<2>(kvert);
	const double CH = std::get<3>(kvert);
	const double SK = std::get<4>(kvert);
	const double CK = std::get<5>(kvert);
	const double gam = std::get<6>(kvert);
	const dcomp hvnorm = std::get<7>(kvert);
	const dcomp kvnorm = std::get<8>(kvert);
	
	const double G1212 = std::real(2.0*gam * (1.0-gam) + (2.0*pow(gam,2)-2.0*gam+1.0) * CH*CK - (pow(1.0-gam,2) + pow(gam,2)*pow(hvnorm,2)*pow(kvnorm,2))*SH*SK);
	const double G1213 = std::real((1.0/(dens*w*c))*(CH*SK - SH*CK*pow(hvnorm,2)));
	const double iG1214 = std::real((1.0/(dens*w*c))*((1.0 - 2.0*gam)*(1.0-CK*CH) + (1.0 - gam - gam*pow(hvnorm,2)*pow(kvnorm,2))*SH*SK));
	const double G1224 = std::real((1.0/(dens*w*c))*(pow(kvnorm,2)*CH*SK - SH*CK));
	const double G1234 = std::real((-1.0/(pow(dens,2)*pow(w,2)*pow(c,2)))*(2.0*(1.0 - CH*CK) + (1.0 + pow(kvnorm,2)*pow(hvnorm,2))*SK*SH));
	const double G1312 = std::real(dens*w*c*(pow(gam,2)*pow(kvnorm,2)*CH*SK - pow(1.0-gam,2)*SH*CK));
	const double G1313 = std::real(CH*CK);
	const double iG1314 = std::real((1.0 - gam)*SH*CK + gam*pow(kvnorm,2)*CH*SK);
	const double G1324 = std::real((-1.0)*pow(kvnorm,2)*SH*SK);
	const double iG1412 = std::real(dens*w*c*((3.0*pow(gam,2) - 2.0*pow(gam,3) - gam)*(1.0 - CH*CK)+(pow(1.0 - gam,3) - pow(gam,3)*pow(hvnorm,2)*pow(kvnorm,2))*SH*SK));
	const double iG1413 = std::real((-1.0)*((1.0 - gam)*CH*SK + gam*pow(hvnorm,2)*SH*CK));
	const double G1414 = std::real(1.0 - 2.0*gam*(1.0 - gam)*(1.0 - CH*CK) + (pow(1.0 - gam,2) + pow(gam,2)*pow(kvnorm,2)*pow(hvnorm,2))*SH*SK);
	const double G2412 = std::real(dens*w*c*(pow(1.0 - gam,2)*CH*SK - pow(gam,2)*SH*CK*pow(hvnorm,2)));
	const double G2413 = std::real((-1.0)*pow(hvnorm,2)*SH*SK);
	const double G3412 = std::real((-1.0)*pow(dens,2)*pow(w,2)*pow(c,2)*(2.0*pow(gam,2)*pow(1.0 - gam,2)*(1.0 - CH*CK) + (pow(1.0 - gam,4)+pow(gam,4)*pow(hvnorm,2)*pow(kvnorm,2))*SH*SK));
		
	return std::make_tuple(G1212,G1213,iG1214,G1224,G1234,G1312,G1313,iG1314,G1324,iG1412,iG1413,G1414,G2412,G2413,G3412);
}

std::tuple<double,double,double,double,double,double,double,double,double,double,double,double,double,double,double> compute_G_vs(const double &c, const double &dn, const double &w, const double &vp, const double &vs, const double &dens){
	const auto kvert = compute_util(w, c, vp, vs, dn, 0);
	const double SH = std::get<2>(kvert);
	const double CH = std::get<3>(kvert);
	const double SK = std::get<4>(kvert);
	const double CK = std::get<5>(kvert);
	const double gam = std::get<6>(kvert);
	const dcomp hvnorm = std::get<7>(kvert);
	const dcomp kvnorm = std::get<8>(kvert);
	
	const auto util_grads = compute_util_grads(w, vs, vp, c, dn, 0);
	const double kvnorm_k = std::get<3>(util_grads);
	const double gam_k = std::get<4>(util_grads);
	const double CKK = std::get<7>(util_grads);
	const double SKK = std::get<9>(util_grads);
	
	const double gG1212 = std::real(2.0*gam_k*(1.0-2.0*gam)*(1.0-CH*CK)+(2.0*pow(gam,2)-2.0*gam+1)*CH*CKK-(2.0*gam_k*(gam-1.0)+2.0*gam*gam_k*pow(hvnorm,2)*pow(kvnorm,2)+pow(gam,2)*pow(hvnorm,2)*kvnorm_k)*SH*SK-(pow((1.0-gam),2)+pow(gam,2)*pow(hvnorm,2)*pow(kvnorm,2))*SH*SKK);
	const double gG1213 = std::real((1.0/(dens*w*c))*(CH*SKK-pow(hvnorm,2)*SH*CKK));
	const double igG1214 = std::real((1.0/(dens*w*c))*((-2.0)*gam_k*(1.0-CH*CK)-(1.0-2.0*gam)*CH*CKK)+(1.0/(dens*w*c))*(((-1.0)*gam_k-gam_k*pow(hvnorm,2)*pow(kvnorm,2)-gam*pow(hvnorm,2)*kvnorm_k)*SH*SK+(1.0-gam-gam*pow(hvnorm,2)*pow(kvnorm,2))*SH*SKK));
	const double gG1224 = std::real((1.0/(dens*w*c))*(kvnorm_k*CH*SK+pow(kvnorm,2)*CH*SKK-SH*CKK));
	const double gG1234 = std::real((1.0/pow((dens*w*c),2))*(-2.0*CH*CKK+(1.0+pow(hvnorm,2)*pow(kvnorm,2))*SH*SKK+pow(hvnorm,2)*kvnorm_k*SH*SK));
	const double gG1312 = std::real(dens*w*c*(2.0*gam*gam_k*pow(kvnorm,2)*CH*SK+pow(gam,2)*kvnorm_k*CH*SK+pow(gam,2)*pow(kvnorm,2)*CH*SKK)+dens*w*c*(2.0*gam_k*(1.0-gam)*SH*CK-pow((1.0-gam),2)*SH*CKK));
	const double gG1313 = std::real(CH*CKK);
	const double igG1314 = std::real((-1.0)*gam_k*SH*CK+(1.0-gam)*SH*CKK+gam_k*pow(kvnorm,2)*CH*SK+gam*kvnorm_k*CH*SK+gam*pow(kvnorm,2)*CH*SKK);
	const double gG1324 = std::real((-1.0)*kvnorm_k*SH*SK-pow(kvnorm,2)*SH*SKK);
	const double igG1412 = std::real(dens*w*c*(gam_k*(-6.0*pow(gam,2)+6.0*gam-1.0)*(1.0-CH*CK)-(pow(gam,2)-gam)*(1.0-2.0*gam)*CH*CKK)+dens*w*c*(-3.0*pow((1.0-gam),2)*gam_k-3.0*pow(gam,2)*gam_k*pow(hvnorm,2)*pow(kvnorm,2)-pow(gam,3)*pow(hvnorm,2)*kvnorm_k)*SH*SK+dens*w*c*(pow((1.0-gam),3)-pow(gam,3)*pow(hvnorm,2)*pow(kvnorm,2))*SH*SKK);
	const double igG1413 = std::real((-1.0)*((1.0-gam)*CH*SKK-gam_k*CH*SK+gam_k*pow(hvnorm,2)*SH*CK+gam*pow(hvnorm,2)*SH*CKK));
	const double gG1414 = std::real(2.0*gam_k*(2.0*gam-1.0)*(1.0-CH*CK)-2.0*(pow(gam,2)-gam)*CH*CKK+(2.0*(gam-1)*gam_k+2.0*gam*gam_k*pow(hvnorm,2)*pow(kvnorm,2)+pow(gam,2)*pow(hvnorm,2)*kvnorm_k)*SH*SK+(pow((1.0-gam),2)+pow(gam,2)*pow(hvnorm,2)*pow(kvnorm,2))*SH*SKK);
	const double gG2412 = std::real(dens*w*c*(2.0*(gam-1.0)*gam_k*CH*SK+pow((1.0-gam),2)*CH*SKK-2.0*gam*gam_k*pow(hvnorm,2)*SH*CK-pow(gam,2)*pow(hvnorm,2)*SH*CKK));
	const double gG2413 = std::real((-1.0)*pow(hvnorm,2)*SH*SKK);
	const double gG3412 = std::real((-1.0)*pow((dens*c*w),2)*(4.0*gam*gam_k*(1.0+2.0*pow(gam,2)-3.0*gam)*(1.0-CH*CK)-2.0*pow(gam,2)*pow((1.0-gam),2)*CH*CKK)-pow((dens*c*w),2)*((pow((1.0-gam),4)+pow(gam,4)*pow(hvnorm,2)*pow(kvnorm,2))*SH*SKK)-pow((dens*c*w),2)*((-4.0*gam_k*pow((1.0-gam),3)+4.0*pow(gam,3)*gam_k*pow(hvnorm,2)*pow(kvnorm,2)+pow(gam,4)*pow(hvnorm,2)*kvnorm_k)*SH*SK));
	return std::make_tuple(gG1212,gG1213,igG1214,gG1224,gG1234,gG1312,gG1313,igG1314,gG1324,igG1412,igG1413,gG1414,gG2412,gG2413,gG3412);
}

std::tuple<double,double,double,double,double,double,double,double,double,double,double,double,double,double,double> compute_G_vp(const double &c, const double &dn, const double &w, const double &vp, const double &vs, const double &dens){
	const auto kvert = compute_util(w, c, vp, vs, dn, 0);
	const double SH = std::get<2>(kvert);
	const double SK = std::get<4>(kvert);
	const double CK = std::get<5>(kvert);
	const double gam = std::get<6>(kvert);
	const dcomp hvnorm = std::get<7>(kvert);
	const dcomp kvnorm = std::get<8>(kvert);
	
	const auto util_grads = compute_util_grads(w, vs, vp, c, dn, 0);
	const double hvnorm_h = std::get<2>(util_grads);
	const double CHH = std::get<6>(util_grads);
	const double SHH = std::get<8>(util_grads);
	
	const double gG1212 = std::real((2.0*pow(gam,2)-2.0*gam+1.0)*CHH*CK-pow(gam,2)*hvnorm_h*pow(kvnorm,2)*SH*SK-(pow((1.0-gam),2)+pow(gam,2)*pow(hvnorm,2)*pow(kvnorm,2))*SHH*SK);
	const double gG1213 = std::real((1.0/(dens*w*c))*(CHH*SK-hvnorm_h*SH*CK-pow(hvnorm,2)*SHH*CK));
	const double igG1214 = std::real((1.0/(dens*w*c))*((2.0*gam-1.0)*CHH*CK+(1.0-gam-gam*pow(hvnorm,2)*pow(kvnorm,2))*SHH*SK-gam*hvnorm_h*pow(kvnorm,2)*SH*SK));
	const double gG1224 = std::real((1.0/(dens*w*c))*(pow(kvnorm,2)*CHH*SK-SHH*CK));
	const double gG1234 = std::real(((-1.0)/pow((w*dens*c),2))*(-2.0*CHH*CK+(1.0+pow(kvnorm,2)*pow(hvnorm,2))*SHH*SK+hvnorm_h*pow(kvnorm,2)*SH*SK));
	const double gG1312 = std::real(dens*w*c*(pow(gam,2)*pow(kvnorm,2)*CHH*SK-pow((1.0-gam),2)*SHH*CK));
	const double gG1313 = std::real(CHH*CK);
	const double igG1314 = std::real((1.0-gam)*SHH*CK+gam*pow(kvnorm,2)*CHH*SK);
	const double gG1324 = std::real((-1.0)*pow(kvnorm,2)*SHH*SK);
	const double igG1412 = std::real(dens*w*c*(-1.0*gam*(gam-1.0)*(1.0-2.0*gam)*CHH*CK-pow(gam,3)*hvnorm_h*pow(kvnorm,2)*SH*SK)+dens*w*c*((pow((1.0-gam),3)-pow(gam,3)*pow(hvnorm,2)*pow(kvnorm,2))*SHH*SK));
	const double igG1413 = std::real((gam-1.0)*CHH*SK-gam*hvnorm_h*SH*CK-gam*pow(hvnorm,2)*SHH*CK);
	const double gG1414 = std::real(2.0*gam*(1.0-gam)*CHH*CK+(pow((1.0-gam),2)+pow(gam,2)*pow(hvnorm,2)*pow(kvnorm,2))*SHH*SK+pow(gam,2)*hvnorm_h*pow(kvnorm,2)*SH*SK);
	const double gG2412 = std::real(dens*w*c*(pow((1.0-gam),2)*CHH*SK-pow(gam,2)*hvnorm_h*SH*CK-pow(gam,2)*pow(hvnorm,2)*SHH*CK));
	const double gG2413 = std::real((-1.0)*(hvnorm_h*SH+pow(hvnorm,2)*SHH)*SK);
	const double gG3412 = std::real((-1.0)*pow((dens*c*w),2)*(-2.0*pow(gam,2)*pow((1-gam),2)*CHH*CK)-pow((dens*c*w),2)*((pow((1-gam),4)+pow(gam,4)*pow(hvnorm,2)*pow(kvnorm,2))*SHH*SK+pow(gam,4)*hvnorm_h*pow(kvnorm,2)*SH*SK));
	
	return std::make_tuple(gG1212,gG1213,igG1214,gG1224,gG1234,gG1312,gG1313,igG1314,gG1324,igG1412,igG1413,gG1414,gG2412,gG2413,gG3412);
}

std::tuple<double,double,double,double,double,double,double,double,double,double,double,double,double,double,double> compute_G_rho(const double &c, const double &dn, const double &w, const double &vp, const double &vs, const double &dens){
	const auto G = compute_G(c,dn,w,vp,vs,dens);
	const double G1213 = std::get<1>(G);
	const double iG1214 = std::get<2>(G);
	const double G1224 = std::get<3>(G);
	const double G1234 = std::get<4>(G);
	const double G1312 = std::get<5>(G);
	const double iG1412 = std::get<9>(G);
	const double G2412 = std::get<12>(G);
	const double G3412 = std::get<14>(G);
	
	const double gG1212 = 0.0;
	const double gG1213 = std::real((-1.0)*G1213/dens);
	const double igG1214 = std::real((-1.0)*iG1214/dens);
	const double gG1224 = std::real((-1.0)*G1224/dens);
	const double gG1234 = std::real((-2.0)*G1234/dens);
	const double gG1312 = std::real(G1312/dens);
	const double gG1313 = 0.0;
	const double igG1314 = 0.0;
	const double gG1324 = 0.0;
	const double igG1412 = std::real(iG1412/dens);
	const double igG1413 = 0.0;
	const double gG1414 = 0.0;
	const double gG2412 = std::real(G2412/dens);
	const double gG2413 = 0.0;
	const double gG3412 = std::real(2.0*G3412/dens);
	return std::make_tuple(gG1212,gG1213,igG1214,gG1224,gG1234,gG1312,gG1313,igG1314,gG1324,igG1412,igG1413,gG1414,gG2412,gG2413,gG3412);
}

std::tuple<double,double,double,double,double> compute_R(const double &w, const double &c, const double &vp, const double &vs, const double &dn, const double &dens, const std::tuple<double,double,double,double,double> &T, const int &param){
	// Recursive layer stacking from bottom to top layer (an i denotes an imaginary subdeterminant e.g. iR1214)
	const double T1212 = std::get<0>(T);
	const double T1213 = std::get<1>(T);
	const double iT1214 = std::get<2>(T);
	const double T1224 = std::get<3>(T);
	const double T1234 = std::get<4>(T);
	
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
	
	const double G1212 = std::get<0>(G);
	const double G1213 = std::get<1>(G);
	const double iG1214 = std::get<2>(G);
	const double G1224 = std::get<3>(G);
	const double G1234 = std::get<4>(G);
	const double G1312 = std::get<5>(G);
	const double G1313 = std::get<6>(G);
	const double iG1314 = std::get<7>(G);
	const double G1324 = std::get<8>(G);
	const double iG1412 = std::get<9>(G);
	const double iG1413 = std::get<10>(G);
	const double G1414 = std::get<11>(G);
	const double G2412 = std::get<12>(G);
	const double G2413 = std::get<13>(G);
	const double G3412 = std::get<14>(G);
	
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
	return std::make_tuple(R1212,R1213,iR1214,R1224,R1234);
}

double compute_R1212(const double &w, const double &c, const std::vector<double> &vp, const std::vector<double> &vs, const double &mu, const std::vector<double> &depth, const std::vector<double> &dens, const int &nlay, const int &param, const int &gradlay){
	// Recursive layer stacking from bottom to top to get R1212
	std::tuple<double,double,double,double,double> R;
	for(int n=nlay-1;n>=0;n--){
		if (n==nlay-1){
			if(n==gradlay){
				if(param==1)
					R = compute_T_vs(w, c, vp[n], vs[n], mu);
				else if(param==2)
					R = compute_T_vp(w, c, vp[n], vs[n], mu);
				else if(param==3)
					R = compute_T_rho(w, c, vp[n], vs[n], mu);
			}
			else{
				R = compute_T(w, c, vp[n], vs[n], mu);
			}
		}
		else {
			double dn = depth[n+1]-depth[n];
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
		R1212_root(const double &w_, const std::vector<double> &vp_, const std::vector<double> &vs_, const double &mu_, const std::vector<double> &depth_, const std::vector<double> &dens_, const int &nlay_):w(w_),vp(vp_),vs(vs_),mu(mu_),depth(depth_),dens(dens_),nlay(nlay_) {};
		double operator()(const double c){
			return compute_R1212(w, c, vp, vs, mu, depth, dens, nlay, 0, -999);
		}
	private:
		const double &w;
		const std::vector<double> &vp;
		const std::vector<double> &vs;
		const double &mu;
		const std::vector<double> &depth;
		const std::vector<double> &dens;
		const int &nlay;
};

struct TerminationCondition{
	// checks whether root bracketing has sufficiently converged
	bool operator() (const double &min, const double &max){
    return abs(min - max) < tolerance;
	}
};

std::vector<vector<double>> get_gc_segments(const double &east0, const double &north0, const double &east1, const double &north1, const double &lon_centr) {
	// Approximates great circle path with loxodrome segments.
	
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
		const double segment_length = line.Distance()/(npts+1);
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
	std::vector<vector<double>> pts(2);
	pts[0] = easting;
	pts[1] = northing;
	return pts;
}

std::vector<vector<double>> get_t_segments(double east0, double north0, const double &east1, const double &north1, const double &event_e, const double &event_n, const double &lon_centr, const std::vector<double> &origin, const double &deast, const double &dnorth, const std::vector<double> &c, const int &ncells_east, const std::vector<double> &dsdvs, const std::vector<double> &dsdvp, const std::vector<double> &dsdrho, const int &nlay){
	// Computes relative phase delay for a station pair and a given earthquake
	
	// Set up coordinate transformations
	Geodesic geod(Constants::WGS84_a(), Constants::WGS84_f());
	TransverseMercatorExact proj(Constants::WGS84_a(), Constants::WGS84_f(), Constants::UTM_k0());
	double event_lat, event_lon;
	proj.Reverse(lon_centr, event_e-false_east, event_n, event_lat, event_lon);
	
	// Path segment is interpreted as line segment (not loxodrome)
	// because we don't care about its length (only orientation).
	// Check which model cells are crossed by inter-station path segment.
	const double slope = (north1-north0)/(east1-east0);
	const double intercept = north0-slope*east0;
	double ecell0 = floor((east0-origin[0])/deast);
	double ncell0 = floor((north0-origin[1])/dnorth);
	const double ecell1 = floor((east1-origin[0])/deast);
	const double ncell1 = floor((north1-origin[1])/dnorth);
	std::vector<vector<double>> times_grads(4);
	std::vector<double> times(1);
	double time = 0.0;
	double mid_e, mid_n, mid_lon, mid_lat, s12, az1, az2, se, sn, dist_segment_e, dist_segment_n, dsedvs, dsndvs, dsedvp, dsndvp, dsedrho, dsndrho;
	std::vector<double> rhograd(dsdrho.size(),0.0), vsgrad(dsdrho.size(),0.0), vpgrad(dsdrho.size(),0.0);
	
	while((abs(ecell0 - ecell1) > 0.0) || (abs(ncell0 - ncell1) > 0.0)) {
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
		proj.Reverse(lon_centr, mid_e-false_east, mid_n, mid_lat, mid_lon);
		geod.Inverse(event_lat, event_lon, mid_lat, mid_lon, s12, az1, az2);
		se = cos((90.0-az2)*M_PI/180.0)*(1.0/c[ncell0 * ncells_east + ecell0]);
		sn = sin((90.0-az2)*M_PI/180.0)*(1.0/c[ncell0 * ncells_east + ecell0]);
		time = time + (se * dist_segment_e + sn * dist_segment_n)*(-1.0);
		
		if(calcgrads==1){
			for(int n=0;n<nlay;n++){
				dsedvs = cos((90.0-az2)*M_PI/180.0)*dsdvs[ncell0 * ncells_east * nlay + ecell0 * nlay + n];
				dsedvp = cos((90.0-az2)*M_PI/180.0)*dsdvp[ncell0 * ncells_east * nlay + ecell0 * nlay + n];
				dsedrho = cos((90.0-az2)*M_PI/180.0)*dsdrho[ncell0 * ncells_east * nlay + ecell0 * nlay + n];
				dsndvs = sin((90.0-az2)*M_PI/180.0)*dsdvs[ncell0 * ncells_east * nlay + ecell0 * nlay + n];
				dsndvp = sin((90.0-az2)*M_PI/180.0)*dsdvp[ncell0 * ncells_east * nlay + ecell0 * nlay + n];
				dsndrho = sin((90.0-az2)*M_PI/180.0)*dsdrho[ncell0 * ncells_east * nlay + ecell0 * nlay + n];
				vsgrad[ncell0 * ncells_east * nlay + ecell0 * nlay + n] = vsgrad[ncell0 * ncells_east * nlay + ecell0 * nlay + n] + (dsedvs * dist_segment_e + dsndvs * dist_segment_n)*(-1.0);
				vpgrad[ncell0 * ncells_east * nlay + ecell0 * nlay + n] = vpgrad[ncell0 * ncells_east * nlay + ecell0 * nlay + n] + (dsedvp * dist_segment_e + dsndvp * dist_segment_n)*(-1.0);
				rhograd[ncell0 * ncells_east * nlay + ecell0 * nlay + n] = rhograd[ncell0 * ncells_east * nlay + ecell0 * nlay + n] + (dsedrho * dist_segment_e + dsndrho * dist_segment_n)*(-1.0);
			}
		}
	}
	dist_segment_e = east1-east0;
	dist_segment_n = north1-north0;
	mid_e = east0 + dist_segment_e/2.0;
	mid_n = north0 + dist_segment_n/2.0;
	proj.Reverse(lon_centr, mid_e-false_east, mid_n, mid_lat, mid_lon);
	geod.Inverse(event_lat, event_lon, mid_lat, mid_lon, s12, az1, az2);
	se = cos((90-az2)*M_PI/180)*(1/c[ncell1 * ncells_east + ecell1]);
	sn = sin((90-az2)*M_PI/180)*(1/c[ncell1 * ncells_east + ecell1]);
	time = time + (se * dist_segment_e + sn * dist_segment_n)*(-1.0);
	
	if(calcgrads==1){
		for(int n=0;n<nlay;n++){
			dsedvs = cos((90.0-az2)*M_PI/180.0)*dsdvs[ncell1 * ncells_east * nlay + ecell1 * nlay + n];
			dsedvp = cos((90.0-az2)*M_PI/180.0)*dsdvp[ncell1 * ncells_east * nlay + ecell1 * nlay + n];
			dsedrho = cos((90.0-az2)*M_PI/180.0)*dsdrho[ncell1 * ncells_east * nlay + ecell1 * nlay + n];
			dsndvs = sin((90.0-az2)*M_PI/180.0)*dsdvs[ncell1 * ncells_east * nlay + ecell1 * nlay + n];
			dsndvp = sin((90.0-az2)*M_PI/180.0)*dsdvp[ncell1 * ncells_east * nlay + ecell1 * nlay + n];
			dsndrho = sin((90.0-az2)*M_PI/180.0)*dsdrho[ncell1 * ncells_east * nlay + ecell1 * nlay + n];
			vsgrad[ncell1 * ncells_east * nlay + ecell1 * nlay + n] = vsgrad[ncell1 * ncells_east * nlay + ecell1 * nlay + n] + (dsedvs * dist_segment_e + dsndvs * dist_segment_n)*(-1.0);
			vpgrad[ncell1 * ncells_east * nlay + ecell1 * nlay + n] = vpgrad[ncell1 * ncells_east * nlay + ecell1 * nlay + n] + (dsedvp * dist_segment_e + dsndvp * dist_segment_n)*(-1.0);
			rhograd[ncell1 * ncells_east * nlay + ecell1 * nlay + n] = rhograd[ncell1 * ncells_east * nlay + ecell1 * nlay + n] + (dsedrho * dist_segment_e + dsndrho * dist_segment_n)*(-1.0);
		}
	}
	
	times[0] = time;
	times_grads[0] = times;
	times_grads[1] = vsgrad;
	times_grads[2] = vpgrad;
	times_grads[3] = rhograd;
	return times_grads;
}

struct weighted_add {
    const double weight;
    double operator() (const double &aa, const double &bb) {
        return aa + bb*weight;
    }
    weighted_add(const double weight_) : weight(weight_) {}
};

std::vector<double> T2w(const std::vector<double> &periods){
	// Conversion of periods to angular frequencies
	std::vector<double> w(periods.size());
	for(int n=0; n<periods.size(); n++){
		w[n] = 2.0*M_PI/periods[n];
		if (verbose==1)
			cout << "Kreisfreq. " << n << ": " << w[n] << " rad/s\n";
	}
	return w;
}

class SurfaceWaves{
	public:
		SurfaceWaves(const double &lon_centr_, const double &dtp_dummy_, const std::vector<double> &depth_, const std::vector<double> &northing_, const std::vector<double> &easting_, const std::vector<double> &periods_, const std::vector<double> &mpn_, const std::vector<double> &mpe_, const std::vector<double> &mpz_, const std::vector<double> &src_rcvr_cmb_, const std::vector<double> &dtp_, const std::vector<double> &event_stat_cmb_, const std::vector<double> &eventx_, const std::vector<double> &eventy_):lon_centr(lon_centr_), dtp_dummy(dtp_dummy_), depth(depth_), northing(northing_), easting(easting_), periods(periods_), mpn(mpn_), mpe(mpe_), mpz(mpz_), src_rcvr_cmb(src_rcvr_cmb_), dtp(dtp_), event_stat_cmb(event_stat_cmb_), eventx(eventx_), eventy(eventy_), dens_grad(northing_.size()*easting_.size()*depth_.size()*periods_.size()), vs_grad(northing_.size()*easting_.size()*depth_.size()*periods_.size()), vp_grad(northing_.size()*easting_.size()*depth_.size()*periods_.size()), dtp_mod(dtp_.size()) {};
		void forward(const std::vector<double> &vs, const std::vector<double> &vp, const std::vector<double> &dens);
		void gradient();
		void residual();
	private:
		const double lon_centr, dtp_dummy;
		const std::vector<double> depth, northing, easting, periods, mpn, mpe, mpz, src_rcvr_cmb, dtp, event_stat_cmb, eventx, eventy;
		std::vector<double> dens_grad, vs_grad, vp_grad, dtp_mod, residual_dt;
};

void SurfaceWaves::residual(){
	residual_dt = std::transform(dtp_mod.begin(), dtp_mod.end(), dtp.begin(), dtp_mod.begin(), std::minus<double>());
	for(int ndt=0; ndt<dtp.size(); ndt++){
		residual_dt[ndt] = abs(residual_dt[ndt]);
	}
}

void SurfaceWaves::gradient(){}

void SurfaceWaves::forward(const std::vector<double> &vs, const std::vector<double> &vp, const std::vector<double> &dens){
	const int nperiods = periods.size();
	const int nstats = mpn.size();
	const int nsrcs = src_rcvr_cmb.size()/2;
	const int nevents = eventx.size();
	const int nevents_per_src = dtp.size()/(nsrcs*nperiods);
	const int NX = northing.size();
	const int NY = easting.size();
	const int NZ = depth.size();

	const double deast = easting[1] - easting[0];
	const double dnorth = northing[1] - northing[0];
	const std::vector<double> model_origin = {easting[0] - deast, northing[0] - dnorth};
	
	const std::vector<double> w = T2w(periods);
	
	// Open output files, write header lines
	ofstream resultfile;
	resultfile.open ("vph_map.out");
	resultfile << "# Easting [m] \t Northing [m] \t Period [s] \t Phase velocity [m/s] \t Difference [m/s] \t No. of iterations";
	
	ofstream delayfile;
	delayfile.open ("delays.out");
	delayfile << "#Event_num \t stat1_num \t stat2_num \t Period [s] \t Phase delay [s]";
	
	ofstream gradfile;
	gradfile.open ("vs_grad.out");
	gradfile << "#Period [s] \t Ncell \t Ecell \t Zcell \t vs_grad[s/(m/s)]";
	
	for(int freq=0; freq<nperiods; freq++){
		cout << "Period: " << periods[freq] << " s.";
		cout << "\n";
		//Vectors to store gradients, dispersion curves
		std::vector<double> dsdrho(dens.size()), dsdvs(vs.size()), dsdvp(vp.size()), vph_map(NX*NY);
		for (int nstep = 0; nstep<NX; nstep++){
			for (int estep = 0; estep<NY; estep++){
				std::vector<double> dens_1D(NZ);
				std::vector<double> vs_1D(NZ);
				std::vector<double> vp_1D(NZ);	
				bool lvz=0;
				for (int n=0; n<NZ; n++){
					// sort velocities, densities into 1D models 
					dens_1D[n] = dens[n+NZ*estep+NY*NZ*nstep];
					vp_1D[n] = vp[n+NZ*estep+NY*NZ*nstep];
					vs_1D[n] = vs[n+NZ*estep+NY*NZ*nstep];
					// check if there's a low velocity zone
					if (n>0 && vs_1D[n]<vs_1D[n-1]){
						lvz=1;
					}
				}
			
				if (vs_1D[0]<=0){
					// This still needs some work. Computation of dispersion curves does not work if top layer has vs=0
					vph_map[estep + NY*nstep] = 0.0;
					resultfile << "\n" << easting[estep] << "\t" << northing[nstep] << "\t" << (2.0*M_PI)/w[freq] << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0;
					continue;
				}
				else{				
					// Calculation of velocity limits
					const double vsmin = *std::min_element(vs_1D.begin(),vs_1D.end());
					const double vsmax = *std::max_element(vs_1D.begin(),vs_1D.end());
					const double vpmin = *std::min_element(vp_1D.begin(),vp_1D.end());
					const double vpmax = *std::max_element(vp_1D.begin(),vp_1D.end());
					const std::vector<double> c_lim{newton_vr(vpmin, vsmin)/1.05, newton_vr(vpmax, vsmax)*1.05};
					
					// step ratio for root bracketing
					const double stepratio = (vsmin - c_lim[0])/(2.0*vsmin);	
	
					if(verbose==1){
						cout << "cmin: " << c_lim[0] << "\t cmax: " << c_lim[1] << "\n";
						cout << "Anzahl Schichten: " << NZ << "\n";
						cout << "Easting: " << easting[estep] << "\t Northing: " << northing[nstep] << "\n";
					}
				
					// Shear modulus bottom layer
					const double mu = pow(vs_1D[NZ-1],2)*dens_1D[NZ-1];
					if (verbose==1)
						cout << "Schermodul unterste Schicht: " << mu << "\n";
						
					// Compute initial R1212 polarization for large period below fundamental mode
					double R1212 = compute_R1212(w[nperiods-1]/10.0, c_lim[0], vp_1D, vs_1D, mu, depth, dens_1D, NZ, 0, -999);
					const bool pol0 = signbit(R1212);
				
					if(verbose==1)
						cout << "Polarisation von R1212 für geringe Geschwindigkeit: " << pol0 << "\n";
				
					double c_last=c_lim[0];	//initial value for c to start search
				
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
						R1212 = compute_R1212(w[freq], c1, vp_1D, vs_1D, mu, depth, dens_1D, NZ, 0, -999);
						pol1 = signbit(R1212);
						
						// If a sign change is found check for mode skipping
						if (pol0!=pol1 && (c1-c0)>(2.0*tolerance)){
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
								R1212 = compute_R1212(w[freq], c2, vp_1D, vs_1D, mu, depth, dens_1D, NZ, 0, -999);
								const bool pol2 = signbit(R1212);
								// if mode skipping detected increase precision (-> decrease step ratio) and return to bracket search
								if (pol2==pol1){
									precision = precision * mode_skip_it;
									c1 = c0;
									if (verbose==1){
										cout << "Error: Mode skipping detected!\n";
									}
									goto cnt;
								}
								// "Downward" search along c-axis for mode skipping (2 runs per default)
								c2 = c2-delta;
								if (verbose==1)
									cout << "New c2: " << c2 << "\n";
							}
						}		
					}
					// If a sign change is found, brackets c0 & c2 are passed to an instance of R1212_root (-> Boost TOMS root finding algorithm)
					std::pair<double, double> brackets; // stores refinded root brackets
					boost::uintmax_t max_iter=500;	// Maximum number of TOMS iterations (500 is probably way to much...)
					R1212_root root(w[freq], vp_1D, vs_1D, mu, depth, dens_1D, NZ);
					brackets = boost::math::tools::toms748_solve(root, c0, c1, TerminationCondition(), max_iter);
					if (lvz == 0 && (brackets.first + brackets.second)/2.0 < c_last) {
						c1 = c_lim[0];
						precision = precision * mode_skip_it;
						if (verbose==1)
							cout << "Error: non-monotonous shape of dispersion curve!\n";
						goto cnt;
					}
					c_last = (brackets.first + brackets.second)/2.0;
					
					vph_map[estep + NY*nstep] = c_last;
					
					// Write output to file
					resultfile << "\n" << easting[estep] << "\t" << northing[nstep] << "\t" << (2.0*M_PI)/w[freq] << "\t" << c_last << "\t" << brackets.second-brackets.first << "\t" << max_iter;

					if(calcgrads==1){
						for(int n=0;n<NZ;n++){
							//Computation of Gradients
							double R_tmp = compute_R1212(w[freq], c_last, vp_1D, vs_1D, mu, depth, dens_1D, NZ, 1, n);
							dsdvs[n+NZ*estep+NY*NZ*nstep] = R_tmp*(-1.0/pow(c_last,2));
							R_tmp = compute_R1212(w[freq], c_last, vp_1D, vs_1D, mu, depth, dens_1D, NZ, 2, n);
							dsdvp[n+NZ*estep+NY*NZ*nstep] = R_tmp*(-1.0/pow(c_last,2));
							R_tmp = compute_R1212(w[freq], c_last, vp_1D, vs_1D, mu, depth, dens_1D, NZ, 3, n);
							dsdrho[n+NZ*estep+NY*NZ*nstep] = R_tmp*(-1.0/pow(c_last,2));
						}
					}
				} 
			} // end of loop over northing
		} // end of loop over easting
		
		// loop over all rays, computes phase delays
		for (int src=0; src<nsrcs; src++){
			std::vector<vector<double>> segments = get_gc_segments(mpe[src_rcvr_cmb[src]], mpn[src_rcvr_cmb[src]], mpe[src_rcvr_cmb[src+nsrcs]], mpn[src_rcvr_cmb[src+nsrcs]], lon_centr);
			std::vector<double> seg_east = segments[0];
			std::vector<double> seg_north = segments[1];
			std::vector<std::vector<double>>().swap(segments);
			for(int event=0; event<nevents_per_src; event++){ //loop over events for each station pair
				if (dtp[freq*nsrcs*nevents_per_src+event*nsrcs+src]==dtp_dummy | event_stat_cmb[event*nsrcs+src]==dtp_dummy){
					// if there is only a dummy value we can skip this period
					dtp_mod[src+nsrcs*event+freq*nsrcs*nevents_per_src] = dtp_dummy;
					delayfile << "\n" << event_stat_cmb[event*nsrcs+src] << "\t" << src_rcvr_cmb[src] << "\t" << src_rcvr_cmb[src+nsrcs] << "\t" << (2.0*M_PI)/w[freq] << "\t" << dtp_dummy;
					continue;
				}
				else {
					std::vector<double> vs_tmpgrd(NX*NY*NZ, 0.0), vp_tmpgrd(NX*NY*NZ, 0.0), dens_tmpgrd(NX*NY*NZ, 0.0);
					double time_total = 0.0;
					for(int seg=0; seg<seg_east.size()-1; seg++){ //loop over great circle segments
						//cout << "Period: " << periods[freq] << " s\t Station combination: " << src << " of " << nsrcs << ",\t Event: " << event << " of " << nevents_per_src << ",\t Segment: " << seg << "\n";
						std::vector<vector<double>> time_segment = get_t_segments(seg_east[seg], seg_north[seg], seg_east[seg+1], seg_north[seg+1], eventy[event_stat_cmb[event*nsrcs+src]], eventx[event_stat_cmb[event*nsrcs+src]], lon_centr, model_origin, deast, dnorth, vph_map, NY, dsdvs, dsdvp, dsdrho, NZ);
						std::vector<double> ts = time_segment[0];					
						time_total = time_total + ts[0];
						
						if(calcgrads==1){
							std::vector<double> tmp = time_segment[1];
							std::transform (vs_tmpgrd.begin(), vs_tmpgrd.end(), tmp.begin(), vs_tmpgrd.begin(), std::plus<double>());
							tmp = time_segment[2];
							std::transform (vp_tmpgrd.begin(), vp_tmpgrd.end(), tmp.begin(), vp_tmpgrd.begin(), std::plus<double>());
							tmp = time_segment[3];
							std::transform (dens_tmpgrd.begin(), dens_tmpgrd.end(), tmp.begin(), dens_tmpgrd.begin(), std::plus<double>());
						}
					}//end loop ever path segments
					dtp_mod[src+nsrcs*event+freq*nsrcs*nevents_per_src] = time_total;
					delayfile << "\n" << event_stat_cmb[event*nsrcs+src] << "\t" << src_rcvr_cmb[src] << "\t" << src_rcvr_cmb[src+nsrcs] << "\t" << (2.0*M_PI)/w[freq] << "\t" << time_total;
					if(calcgrads==1){
						const double residual = abs(dtp[freq*nsrcs*nevents_per_src+event*nsrcs+src]-time_total);
						const int model_length = NX*NY*NZ;
						const int element0 = freq*model_length;
						const int element1 = ((freq+1)*model_length)-1;
						std::transform(vs_grad.begin()+element0, vs_grad.begin()+element1, vs_tmpgrd.begin(), vs_grad.begin()+element0, weighted_add(residual));
						std::transform(vp_grad.begin()+element0, vp_grad.begin()+element1, vp_tmpgrd.begin(), vp_grad.begin()+element0, weighted_add(residual));
						std::transform(dens_grad.begin()+element0, dens_grad.begin()+element1, dens_tmpgrd.begin(), dens_grad.begin()+element0, weighted_add(residual));
					}
				}
			} //end loop over events
		} // end loop rays
	}// end loop frequencies
	
	resultfile << "\n";
	resultfile.close();
	
	delayfile << "\n";
	delayfile.close();
	
	int index = 0;
	if (calcgrads==1){
		for(int freq=0; freq<nperiods; freq++){
			for(int nstep=0; nstep<NX; nstep++){
				for(int estep=0; estep<NY; estep++){
					for(int n=0; n<NZ; n++){
						gradfile << "\n" << periods[freq] << "\t" << nstep << "\t" << estep << "\t" << n << "\t" << vs_grad[index];
						index++;
					}
				}
			}
		}		
	}
	
	gradfile << "\n";
	gradfile.close();
}

int main(){
	
	double lon_centr, dtp_dummy;
	
	std::vector<std::vector<double>> dtp_data = read_data(datafile, dtp_dummy, lon_centr);
	std::vector<double> periods = dtp_data[0];
	std::vector<double> mpn = dtp_data[1];
	std::vector<double> mpe = dtp_data[2];
	std::vector<double> mpz = dtp_data[3];
	std::vector<double> src_rcvr_cmb = dtp_data[4];
	std::vector<double> dtp = dtp_data[5];
	std::vector<double> event_stat_cmb = dtp_data[6];
	std::vector<double> eventx = dtp_data[7];
	std::vector<double> eventy = dtp_data[8];
	std::vector<std::vector<double>>().swap(dtp_data);
	
	std::vector<std::vector<double>> model = read_model(dens_file, vs_file, vp_file);
	std::vector<double> depth = model[0];
	std::vector<double> north = model[1];
	std::vector<double> east = model[2];
	std::vector<double> dens_all = model[3];
	std::vector<double> vp_all = model[4];
	std::vector<double> vs_all = model[5];
	std::vector<std::vector<double>>().swap(model);
	
	SurfaceWaves SW(lon_centr, dtp_dummy, depth, north, east, periods, mpn, mpe, mpz, src_rcvr_cmb, dtp, event_stat_cmb, eventx, eventy);
	SW.forward(vs_all, vp_all, dens_all);

	return 0;
}

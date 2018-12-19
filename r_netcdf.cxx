/*
 * r_netcdf.cxx
 * 
 * Copyright 2018 Bernhard Weise <bweise@eight>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */


#include <iostream>
#include <netcdf>
using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;


int main()
{
	
	static const int NX = 127; // Number of cells in nc file
	static const int NY = 127;
	static const int NZ = 31;
	
	double depth_arr[NZ];		// Define data variables
	double MPX[NX];
	double MPY[NY];
	double dens_arr[NX][NY][NZ];
	
	NcFile densFile("/home/bweise/bmw/WINTERC/dens_na.nc", NcFile::read);
	
	NcVar depthIn=densFile.getVar("Depth");
	depthIn.getVar(depth_arr);
	
	NcVar MPXIn=densFile.getVar("Northing");
	MPXIn.getVar(MPX);
	
	NcVar MPYIn=densFile.getVar("Easting");
	MPYIn.getVar(MPY);
	
	NcVar densIn=densFile.getVar("Density");
	densIn.getVar(dens_arr);
	
	std::vector<double> depth(depth_arr, depth_arr + sizeof depth_arr / sizeof depth_arr[0]);
	depth.insert ( depth.begin() , 0 );
	depth.pop_back();
	
	double arr2[NZ];
	
	for (int n=0;n<NZ;n++){
	arr2[n] = dens_arr[1][1][n];}
	
	std::vector<double> v2(arr2, arr2 + sizeof arr2 / sizeof arr2[0]);
	
	for (int n=0; n<v2.size(); n++){
		cout << v2[n] << "\n";}
	
	return 0;
}


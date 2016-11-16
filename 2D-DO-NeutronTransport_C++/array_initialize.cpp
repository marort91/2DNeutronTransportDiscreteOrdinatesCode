#include <vector>

//*******************************************************************************************************************************

void array_initialize( int Nx, int Ny, int ord, std::vector<std::vector<std::vector<double> > > &half_angular_flux_x, 
	                   std::vector<std::vector<std::vector<double> > > &half_angular_flux_y,
	                   std::vector<std::vector<std::vector<double> > > &angular_flux,
	                   std::vector<std::vector<double> > &scalar_flux,
	                   std::vector<std::vector<double> > &S, std::vector<std::vector<double> > &Q,
	                   std::vector<std::vector<double> > &sigt, std::vector<std::vector<double> > &sigs,
	                   std::vector<std::vector<double> > &nusigf )

//*******************************************************************************************************************************
//
//  ARRAY_INITIALIZE - Array Initialization for 2D Neutron Transport Solver
//
//  Licensing:
//             The MIT License (MIT)
//
//			   Copyright (c) [2016] [Mario I. Ortega]
//
//             Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
//             associated documentation files (the "Software"), to deal in the Software without restriction, including
//             without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
//             copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the 
//             following conditions:
//
//             The above copyright notice and this permission notice shall be included in all copies or substantial 
//             portions of the Software.
//
//             THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT 
//             LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
//             IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
//             WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
//             OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//  Author: 
//             Mario I. Ortega, University of California, Berkeley
//
//  Modified: 
//             2 May 2016
//
//  Parameters: 
//             Input, int Nx, mesh size in x-direction.
//			   Input, int Ny, mesh size in y-direction.
//			   Input, int ord, number of discrete ordinates.
//			   Input, vector<vector<vector<double> > > &half_angular_flux_x, pointer to angular flux half values in x-direction.
//			   Input, vector<vector<vector<double> > > &half_angular_flux_y, pointer to angular flux half values in y-direction.
//			   Input, vector<vector<double> > &angular_flux, pointer to angular flux.
//			   Input, vector<vector<double> > &S, pointer to scalar flux scattering source.
//			   Input, vector<vector<double> > &Q, pointer to external source.
{

	//Set up three-dimensional x-interval half-angular flux.
	half_angular_flux_x.resize(Nx);
	for ( int i = 0; i < Nx; i++ )
	{
		half_angular_flux_x[i].resize(Nx+1);

		for ( int j = 0; j < Nx+1; j++ )
		{
			half_angular_flux_x[i][j].resize(ord);
		}
	}

	// Set up three-dimensional y-interval half-angular flux.
	half_angular_flux_y.resize(Ny+1);
	for ( int i = 0; i < Ny+1; i++ )
	{
		half_angular_flux_y[i].resize(Ny);

		for ( int j = 0; j < Ny; j++ )
			half_angular_flux_y[i][j].resize(ord);
	}

	// Set up angular flux array.
	angular_flux.resize(Nx);
	for ( int i = 0; i < Nx; i++ )
	{
		angular_flux[i].resize(Ny);

		for ( int j = 0; j < Ny; j++ )
			angular_flux[i][j].resize(ord);
	}

	// Set up scalar flux array.
	scalar_flux.resize(Nx);
	for ( int i = 0; i < Nx; i++ )
	{
		scalar_flux[i].resize(Ny);
	}

	// Set up source and scattering scalar flux arrays.
	S.resize(Nx);
	Q.resize(Nx);
	sigt.resize(Nx);
	sigs.resize(Nx);
	nusigf.resize(Nx);
	for ( int i = 0; i < Ny; i++ )
	{
		S[i].resize(Ny);
		Q[i].resize(Ny);
		sigt[i].resize(Ny);
		sigs[i].resize(Ny);
		nusigf[i].resize(Ny);
	}

}

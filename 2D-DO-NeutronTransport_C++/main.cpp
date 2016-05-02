#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include <string>
#include "level_sym_quad.h"

using namespace std;

int main()
{

	// Spatial and Angular Variable Initialization
	int N = 2; int Nx = 2; int Ny = 2; int ord = N*(N+2)/2;
	double xL = 0; double xR = 1; double yB = 0; double yT = 1;
	std::vector<double> x; std::vector<double> y; double dx; double dy;

	dx = ( xR - xL )/(float)Nx;
	dy = ( yT - yB )/(float)Ny;

	// Angular and Scalar Flux Initialization
	std::vector<double> mu; std::vector<double> eta; std::vector<double> wi;
	std::vector<std::vector<std::vector<double> > > half_angular_flux_x;
	std::vector<std::vector<std::vector<double> > > half_angular_flux_y;
	
	//vector<vector<vector<double> > > half_angular_flux_x (Nx,vector<vector<double> > (Nx+1,vector <double>(ord)));
	//vector<vector<vector<double> > > half_angular_flux_y (Ny+1,vector<vector<double> > (Ny,vector <double>(ord)));
	//vector<vector<vector<double> > > angular_flux (Nx,vector<vector<double> > (Ny,vector <double>(ord)));
	
	//vector<vector<double> > scalar_flux (Nx, vector<double> (Ny));

	std::vector<std::vector<std::vector<double> > > angular_flux;
	std::vector<std::vector<double> > scalar_flux;

	// Source Initialization
	//vector<vector<double> > S (Nx, vector<double> (Ny));
	//vector<vector<double> > Q (Nx, vector<double> (Ny));
	std::vector<std::vector<double> > S; std::vector<std::vector<double> > Q;

	// Boundary Condition Flag
	//std::string bc = "Larsen-2D";
	int bc = 0;
	
	level_sym_quad( N, mu, eta, wi );
	array_initialize( Nx, Ny, ord, half_angular_flux_x, half_angular_flux_y, angular_flux, scalar_flux, S, Q );
	spatial_discretize( xL, xR, Nx, dx, yB, yT, Ny, dy, x, y );
	//set_boundary_condition( bc, Nx, Ny, ord, half_angular_flux_x, half_angular_flux_y );

	scatter_scalarflux_source(Ny,Ny,ord,angular_flux,scalar_flux,wi);

	int Q0 = 1;

	const_external_source_def(Q,Q0,Nx,Ny);

	//cout << half_angular_flux_x[3].size() << '\n';

	// Transport Sweep

	// Physical Constants
	double sigt = 1; double sigs = 0; 

	int itermax = 1;

	for ( int iter = 0; iter < itermax; iter++ )
	{

		for ( int k = 0; k < ord; k++ )
		{
			if ( mu[k] > 0 && eta[k] > 0 )
			{
				for ( int j = 0; j < Ny; j++ )
				{
					for ( int i = 0; i < Nx; i++ )
					{
						angular_flux[j][i][k] = ( 2*mu[k]*half_angular_flux_x[j][i][k]/dx + 2*eta[k]*half_angular_flux_y[j][i][k]/dy + Q[i][j] )/
												( 2*mu[k]/dx + 2*eta[k]/dy + sigt );
						cout << j << '	' << i << '	' << angular_flux[j][i][k] << '\n';
						half_angular_flux_x[j][i+1][k] = 2*angular_flux[j][i][k] - half_angular_flux_x[j][i][k];
					}

					for ( int m = 0; m < Nx; m++ )
					{
						half_angular_flux_y[j+1][m][k] = 2*angular_flux[j][m][k] - half_angular_flux_y[j][m][k];
					}
				}
			}
			else if ( mu[k] < 0 && eta[k] > 0 )
			{
				for ( int j = 0; j < Ny; j++ )
				{
					for ( int i = Nx-1; i >= 0; i-- )
					{
						angular_flux[j][i][k] = ( -2*mu[k]*half_angular_flux_x[j][i+1][k]/dx + 2*eta[k]*half_angular_flux_y[j][i][k]/dy + Q[i][j] )/
												( -2*mu[k]/dx + 2*eta[k]/dy + sigt );
						cout << j << '	' << i << '	' << angular_flux[j][i][k] << '\n';
						half_angular_flux_x[j][i][k] = 2*angular_flux[j][i][k] - half_angular_flux_x[j][i+1][k];
					}

					for ( int m = Nx-1; m >= 0; m-- )
					{
						half_angular_flux_y[j+1][m][k] = 2*angular_flux[j][m][k] - half_angular_flux_y[j][m][k];
					}
				}
			}
			else if ( mu[k] > 0 && eta[k] < 0 )
			{
				for ( int j = Ny-1; j >= 0; j-- )
				{
					for ( int i = 0; i < Nx; i++ )
					{
						angular_flux[j][i][k] = ( 2*mu[k]*half_angular_flux_x[j][i][k]/dx - 2*eta[k]*half_angular_flux_y[j+1][i][k]/dy + Q[i][j] )/
												( 2*mu[k]/dx - 2*eta[k]/dy + sigt );
						cout << j << '	' << i << '	' << angular_flux[j][i][k] << '\n';
						half_angular_flux_x[j][i+1][k] = 2*angular_flux[j][i][k] - half_angular_flux_x[j][i][k];						
					}

					for ( int m = 0; m < Nx; m++ )
					{
						half_angular_flux_y[j][m][k] = 2*angular_flux[j][m][k] - half_angular_flux_y[j+1][m][k];
					}
				}
			}
			else if ( mu[k] < 0 && eta[k] < 0 )
			{
				for ( int j = Ny-1; j >= 0; j-- )
				{
					for ( int i = Nx-1; i >= 0; i-- )
					{
						angular_flux[j][i][k] = ( -2*mu[k]*half_angular_flux_x[j][i+1][k]/dx - 2*eta[k]*half_angular_flux_y[j+1][i][k]/dy + Q[i][j] )/
												( -2*mu[k]/dx - 2*eta[k]/dy + sigt );
						cout << j << '	' << i << '	' << angular_flux[j][i][k] << '\n';
						half_angular_flux_x[j][i][k] = 2*angular_flux[j][i][k] - half_angular_flux_x[j][i+1][k];
					}

					for ( int m = 0; m < Nx; m++ )
					{
						half_angular_flux_y[j][m][k] = 2*angular_flux[j][m][k] - half_angular_flux_y[j+1][m][k];
					}
				}
			}
		}

	}

	for ( int i = 0; i < Nx; i++ )
	{
		for ( int j = 0; j < Ny; j++ )
		{
			for ( int k = 0; k < ord; k++ )
			{
				//cout << angular_flux[i][j][k] << '\n';
			}

			cout << '\n';
		}

	}

	//cout << mu[0] << '\n';

}
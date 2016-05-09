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
	//int N = 16; int Nx = 100; int Ny = 100; int ord = N*(N+2)/2;
	//double xL = 0; double xR = 1; double yB = 0; double yT = 1;

	int N; int Nx; int Ny; double xL; double xR; double yB; double yT;
				int bc; double sigt; double sigs; double tol; int ord;

	input_read( N, Nx, xL, xR, Ny, yB, yT, bc, sigt, sigs, tol );

	ord = N*(N+2)/2;

	std::vector<double> x; std::vector<double> y; double dx; double dy;
	double residual;

	dx = ( xR - xL )/(float)Nx;
	dy = ( yT - yB )/(float)Ny;

	// Angular and Scalar Flux Initialization
	std::vector<double> mu; std::vector<double> eta; std::vector<double> wi;
	std::vector<std::vector<std::vector<double> > > half_angular_flux_x;
	std::vector<std::vector<std::vector<double> > > half_angular_flux_y;

	std::vector<std::vector<std::vector<double> > > angular_flux;
	std::vector<std::vector<double> > scalar_flux;
	std::vector<std::vector<double> > scalar_flux_previous;

	// Source Initialization
	std::vector<std::vector<double> > S; std::vector<std::vector<double> > Q;

	// Boundary Condition Flag
	//std::string bc = "Larsen-2D";
	//int bc = 0;
	//int reflect_idx;
	
	level_sym_quad( N, mu, eta, wi );
	array_initialize( Nx, Ny, ord, half_angular_flux_x, half_angular_flux_y, angular_flux, scalar_flux, S, Q );
	spatial_discretize( xL, xR, Nx, dx, yB, yT, Ny, dy, x, y );

	int Q0 = 1;

	const_external_source_def(Q,Q0,Nx,Ny);

	// Transport Sweep

	// Physical Constants
	//double sigt = 1; double sigs = 0.99999; 

	//double tol = 1e-8;

	int itermax = 100000;

	for ( int iter = 0; iter < itermax; iter++ )
	{

		set_boundary_condition( bc, Nx, Ny, ord, 0, mu, eta, half_angular_flux_x, half_angular_flux_y );

		calculate_scalarflux( Ny,Ny,ord,angular_flux,scalar_flux,wi );

		scalar_flux_previous = scalar_flux;

		source_external_scattering( Q, scalar_flux, sigs, Q0, Nx, Ny );

		for ( int k = 0; k < ord; k++ )
		{
			if ( mu[k] > 0 && eta[k] > 0 )
			{
				if ( bc == 4)
				{
					set_boundary_condition( bc, Nx, Ny, ord, k, mu, eta, half_angular_flux_x, half_angular_flux_y );
				}

				for ( int j = 0; j < Ny; j++ )
				{
					for ( int i = 0; i < Nx; i++ )
					{
						angular_flux[j][i][k] = ( 2*mu[k]*half_angular_flux_x[j][i][k]/dx + 2*eta[k]*half_angular_flux_y[j][i][k]/dy + Q[i][j] )/
												( 2*mu[k]/dx + 2*eta[k]/dy + sigt );
						//cout << j << '	' << i << '	' << angular_flux[j][i][k] << '\n';
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
				if ( bc == 2 || bc == 4 )
				{
					set_boundary_condition( bc, Nx, Ny, ord, k, mu, eta, half_angular_flux_x, half_angular_flux_y );
				}

				for ( int j = 0; j < Ny; j++ )
				{
					for ( int i = Nx-1; i >= 0; i-- )
					{
						angular_flux[j][i][k] = ( -2*mu[k]*half_angular_flux_x[j][i+1][k]/dx + 2*eta[k]*half_angular_flux_y[j][i][k]/dy + Q[i][j] )/
												( -2*mu[k]/dx + 2*eta[k]/dy + sigt );
						//cout << j << '	' << i << '	' << angular_flux[j][i][k] << '\n';
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
				if ( bc == 2 || bc == 4 )
				{
					set_boundary_condition( bc, Nx, Ny, ord, k, mu, eta, half_angular_flux_x, half_angular_flux_y );
				}

				for ( int j = Ny-1; j >= 0; j-- )
				{
					for ( int i = 0; i < Nx; i++ )
					{
						angular_flux[j][i][k] = ( 2*mu[k]*half_angular_flux_x[j][i][k]/dx - 2*eta[k]*half_angular_flux_y[j+1][i][k]/dy + Q[i][j] )/
												( 2*mu[k]/dx - 2*eta[k]/dy + sigt );
						//cout << j << '	' << i << '	' << angular_flux[j][i][k] << '\n';
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
				if ( bc == 2 || bc == 4 )
				{
					set_boundary_condition( bc, Nx, Ny, ord, k, mu, eta, half_angular_flux_x, half_angular_flux_y );
				}

				for ( int j = Ny-1; j >= 0; j-- )
				{
					for ( int i = Nx-1; i >= 0; i-- )
					{
						angular_flux[j][i][k] = ( -2*mu[k]*half_angular_flux_x[j][i+1][k]/dx - 2*eta[k]*half_angular_flux_y[j+1][i][k]/dy + Q[i][j] )/
												( -2*mu[k]/dx - 2*eta[k]/dy + sigt );
						//cout << j << '	' << i << '	' << angular_flux[j][i][k] << '\n';
						half_angular_flux_x[j][i][k] = 2*angular_flux[j][i][k] - half_angular_flux_x[j][i+1][k];
					}

					for ( int m = 0; m < Nx; m++ )
					{
						half_angular_flux_y[j][m][k] = 2*angular_flux[j][m][k] - half_angular_flux_y[j+1][m][k];
					}
				}
			}
		}

		calculate_scalarflux( Ny,Ny,ord,angular_flux,scalar_flux,wi );

		residual = norm(Nx,Ny,scalar_flux,scalar_flux_previous);

		cout << iter << '\n';
		cout << "Residual:" << " " << residual << '\n';

		if ( residual < tol )
		{
			cout << residual << '\n';
			//cout << iter << '\n';
			break;
		}

	}

	output_write( Nx,Ny,scalar_flux );

}
#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include <string>
#include "level_sym_quad.h"

using namespace std;

int main()
{

	cout << '\n';
	cout << "2D Discrete Ordinates Radiation Transport Code" << '\n';
	cout << "Author: Mario I. Ortega" << '\n';
	cout << "Institution: University of California, Berkeley" <<'\n';
	cout << '\n';

	// Spatial and Angular Variable Initialization
	int N; int Nx; int Ny; double xL; double xR; double yB; double yT;
	int bc; double tol; int ord; int Egrp; double keff; double keff_previous;

	std::string srcfid;
	std::string sigtfid;
	std::string sigsfid;
	std::string nusigffid;
	std::string calc_mode;
	std::string chifid;

	input_read( N, Nx, xL, xR, Ny, yB, yT, bc, sigtfid, sigsfid, nusigffid, tol, srcfid, Egrp, calc_mode, chifid );

	cout << '\n';

	ord = N*(N+2)/2;

	std::vector<double> x; std::vector<double> y; double dx; double dy;
	double residual;
	double res_inner;

	dx = ( xR - xL )/(float)Nx;
	dy = ( yT - yB )/(float)Ny;

	// Angular and Scalar Flux Initialization
	std::vector<double> mu; std::vector<double> eta; std::vector<double> wi;
	std::vector<std::vector<std::vector<std::vector<double> > > > half_angular_flux_x;
	std::vector<std::vector<std::vector<std::vector<double> > > > half_angular_flux_y;

	std::vector<std::vector<std::vector<std::vector<double> > > > angular_flux;
	std::vector<std::vector<std::vector<double> > > scalar_flux;
	std::vector<std::vector<double> > group_scalar_flux;
	std::vector<std::vector<double> > group_scalar_flux_previous;
	std::vector<std::vector<std::vector<double> > > scalar_flux_outer;
	std::vector<std::vector<std::vector<double> > > scalar_flux_previous;
	std::vector<std::vector<std::vector<double> > > scalar_flux_previous_inner;
	std::vector<std::vector<std::vector<double> > > scalar_flux_outer_previous;

	// Source Initialization
	std::vector<std::vector<std::vector<double> > > S; std::vector<std::vector<std::vector<double> > > Q;
	std::vector<std::vector<std::vector<double> > > sigt; std::vector<std::vector<double> > sigs; 
	std::vector<std::vector<std::vector<double> > > nusigf;
	std::vector<double> chi;
	
	level_sym_quad( N, mu, eta, wi );
	array_initialize( Nx, Ny, ord, Egrp, half_angular_flux_x, half_angular_flux_y, angular_flux, scalar_flux, S, Q, sigt, sigs, nusigf, chi,
	group_scalar_flux );
	spatial_discretize( xL, xR, Nx, dx, yB, yT, Ny, dy, x, y );

	if ( calc_mode == "Source" )
	{
		cout << "Reading source file: " << srcfid << '\n';
		source_file_read( S, srcfid, Nx, Ny, calc_mode );
		keff = 0;
		cout << '\n';	
	}
	else
	{
		cout << "Criticality Calculation" << '\n';
		cout << "External Source Set to Zero" << '\n';
		set_extern_src_zero( Nx, Ny, Egrp, S);
		scalar_flux_critical_guess(Nx,Ny,Egrp,scalar_flux);
		//calculate_scalarflux( Ny, Ny, ord, Egrp, angular_flux, scalar_flux, wi );
		keff = 1;
		cout << '\n';
	}

	for ( int i = 0; i < Egrp; i++ )
	{
		for ( int j = 0; j < Nx; j++ )
		{
			for ( int k = 0; k < Ny; k++ )
			{

				cout << S[j][k][i] << " ";

			}

			cout << endl;
		}

		cout << endl;
	}

	cout << "Press enter to continue..." << '\n';
	system("read");
	
	xs_file_read( sigt, sigtfid, sigs, sigsfid, nusigf, nusigffid, chi, chifid, Nx, Ny, Egrp);

		for ( int i = 0; i < Egrp; i++ )
	{
		for ( int j = 0; j < Nx; j++ )
		{
			for ( int k = 0; k < Ny; k++ )
			{

				//cout << nusigf[j][k][i] << " ";

			}

		//	cout << endl;
		}

		//cout << endl;
	}

	cout << sigs[0][0] << endl;

	//for ( int e = 0; e < Egrp; e++ )
	//{
	//	cout << chi[e] << " ";
	//}

	cout << '\n';
	cout << "Beginning transport sweep! " << '\n';
	cout << '\n';

	// Transport Sweep

	int iterout = 1000;
	int iterin = 1000;

	std::vector<std::vector<double> > group_scalar_flux_prev;
	group_scalar_flux_prev.resize(Nx);
	for ( int i = 0; i < Ny; i++ )
	{
		group_scalar_flux_prev[i].resize(Ny);
	}

	//cout << "SEG FAULT! " << '\n';

	for ( int iter = 0; iter < iterout; iter++ )
	{

		scalar_flux_outer_previous = scalar_flux;

		for ( int Eiter = 0; Eiter < Egrp; Eiter++ )
		{
			//src_extrn_scalarflux( S, Q, scalar_flux, sigs, Nx, Ny, Egrp, keff, chi, nusigf );
			

			for ( int Eiiner = 0; Eiiner < iterin; Eiiner++)
			{
				//set_boundary_condition( bc, Nx, Ny, ord, Egrp, 0, mu, eta, half_angular_flux_x, half_angular_flux_y );
				//calculate_scalarflux( Ny, Ny, ord, Egrp, angular_flux, scalar_flux, wi );
				//if ( Eiiner == 0 )
				//{
				//for ( int i = 0; i < Nx; i ++)
				//{
				//	for ( int j = 0; j < Ny; j++ )
				//	{
				//		group_scalar_flux_prev[i][j] = scalar_flux[i][j][Eiter];
				//	}
				//}
				//}
				//else
				//{
				//	group_scalar_flux_prev = group_scalar_flux;
				//}

				//cout << "SEG FAULT 2! " << '\n';
				//calculate_group_scalarflux(Nx,Ny,ord,Eiter,angular_flux,group_scalar_flux,wi);
				//cout << "SEG FAULT 3! " << '\n';
				//group_src(S,Q,scalar_flux,group_scalar_flux,sigs,Nx,Ny,Eiter,keff,chi,nusigf,Egrp,calc_mode);
				//group_scalar_flux_previous = group_scalar_flux;
				//angular_flux_critical_guess(Nx,Ny,ord,Egrp,angular_flux);
				scalar_flux_previous = scalar_flux;
				src_extrn_scalarflux( S, Q, scalar_flux, sigs, Nx, Ny, Egrp, keff, chi, nusigf, calc_mode );

				for ( int k = 0; k < ord; k++ )
				{
					if ( mu[k] > 0 && eta[k] > 0 )
					{
						if ( bc == 4)
						{
							set_boundary_condition( bc, Nx, Ny, ord, Egrp, k, mu, eta, half_angular_flux_x, half_angular_flux_y );
						}

						for ( int j = 0; j < Ny; j++ )
						{
							for ( int i = 0; i < Nx; i++ )
							{
								angular_flux[j][i][k][Eiter] = ( 2*mu[k]*half_angular_flux_x[j][i][k][Eiter]/dx + 2*eta[k]*half_angular_flux_y[j][i][k][Eiter]/dy + Q[i][j][Eiter] )/
												( 2*mu[k]/dx + 2*eta[k]/dy + sigt[i][j][Eiter] );
								half_angular_flux_x[j][i+1][k][Eiter] = 2*angular_flux[j][i][k][Eiter] - half_angular_flux_x[j][i][k][Eiter];
							}

							for ( int m = 0; m < Nx; m++ )
							{
								half_angular_flux_y[j+1][m][k][Eiter] = 2*angular_flux[j][m][k][Eiter] - half_angular_flux_y[j][m][k][Eiter];
							}
						}
					}
					else if ( mu[k] < 0 && eta[k] > 0 )
					{
						if ( bc == 2 || bc == 4 )
						{
							set_boundary_condition( bc, Nx, Ny, ord, Egrp, k, mu, eta, half_angular_flux_x, half_angular_flux_y );
						}

				for ( int j = 0; j < Ny; j++ )
				{
					for ( int i = Nx-1; i >= 0; i-- )
					{
						angular_flux[j][i][k][Eiter] = ( -2*mu[k]*half_angular_flux_x[j][i+1][k][Eiter]/dx + 2*eta[k]*half_angular_flux_y[j][i][k][Eiter]/dy + Q[i][j][Eiter] )/
												( -2*mu[k]/dx + 2*eta[k]/dy + sigt[i][j][Eiter] );
						half_angular_flux_x[j][i][k][Eiter] = 2*angular_flux[j][i][k][Eiter] - half_angular_flux_x[j][i+1][k][Eiter];
					}

					for ( int m = Nx-1; m >= 0; m-- )
					{
						half_angular_flux_y[j+1][m][k][Eiter] = 2*angular_flux[j][m][k][Eiter] - half_angular_flux_y[j][m][k][Eiter];
					}
				}
			}
			else if ( mu[k] > 0 && eta[k] < 0 )
			{
				if ( bc == 2 || bc == 4 )
				{
					set_boundary_condition( bc, Nx, Ny, ord, Egrp, k, mu, eta, half_angular_flux_x, half_angular_flux_y );
				}

				for ( int j = Ny-1; j >= 0; j-- )
				{
					for ( int i = 0; i < Nx; i++ )
					{
						angular_flux[j][i][k][Eiter] = ( 2*mu[k]*half_angular_flux_x[j][i][k][Eiter]/dx - 2*eta[k]*half_angular_flux_y[j+1][i][k][Eiter]/dy + Q[i][j][Eiter] )/
												( 2*mu[k]/dx - 2*eta[k]/dy + sigt[i][j][Eiter] );
						half_angular_flux_x[j][i+1][k][Eiter] = 2*angular_flux[j][i][k][Eiter] - half_angular_flux_x[j][i][k][Eiter];						
					}

					for ( int m = 0; m < Nx; m++ )
					{
						half_angular_flux_y[j][m][k][Eiter] = 2*angular_flux[j][m][k][Eiter] - half_angular_flux_y[j+1][m][k][Eiter];
					}
				}
			}
			else if ( mu[k] < 0 && eta[k] < 0 )
			{
				if ( bc == 2 || bc == 4 )
				{
					set_boundary_condition( bc, Nx, Ny, ord, Egrp, k, mu, eta, half_angular_flux_x, half_angular_flux_y );
				}

				for ( int j = Ny-1; j >= 0; j-- )
				{
					for ( int i = Nx-1; i >= 0; i-- )
					{
						angular_flux[j][i][k][Eiter] = ( -2*mu[k]*half_angular_flux_x[j][i+1][k][Eiter]/dx - 2*eta[k]*half_angular_flux_y[j+1][i][k][Eiter]/dy + Q[i][j][Eiter] )/
												( -2*mu[k]/dx - 2*eta[k]/dy + sigt[i][j][Eiter] );
						half_angular_flux_x[j][i][k][Eiter] = 2*angular_flux[j][i][k][Eiter] - half_angular_flux_x[j][i+1][k][Eiter];
					}

					for ( int m = 0; m < Nx; m++ )
					{
						half_angular_flux_y[j][m][k][Eiter] = 2*angular_flux[j][m][k][Eiter] - half_angular_flux_y[j+1][m][k][Eiter];
					}
				}
			}

		}

		//calculate_group_scalarflux(Nx,Ny,ord,Eiter,angular_flux,group_scalar_flux,wi);
		calculate_scalarflux(Nx,Ny,ord,Egrp,angular_flux,scalar_flux,wi);
		res_inner = inner_norm(Nx,Ny,Eiter,scalar_flux,scalar_flux_previous);

		cout << "Energy Group: " << Eiter+1 << " " << "Inner Iteration: " << Eiiner+1 << " " << "Inner Error: " << " " << res_inner << '\n';

		if ( res_inner < tol )
		{
			cout << '\n' << endl;
			break;
		}

		}

	    }

	    //cout << "Seg Fault Here" << '\n';

		calculate_scalarflux( Nx,Ny,ord,Egrp,angular_flux,scalar_flux,wi );

		scalar_flux_outer = scalar_flux;

		//cout << "Seg Fault Here" << '\n';

		residual = norm(Nx,Ny,Egrp,scalar_flux,scalar_flux_outer_previous);

		//cout << "Seg Fault Here" << '\n';

		if ( calc_mode != "Source" )
		{
			keff_previous = keff;
			keff = kEigCalc(Nx,Ny,Egrp,scalar_flux,scalar_flux_outer_previous,nusigf,keff);
			//cout << "Seg Fault Here" << '\n';
			cout << "Outer Iteration: " << iter+1 << " " << "Residual: " << " " << residual << '\n';
			cout << "k-effective:" << keff << " " << " " << "Iteration Difference: " << " " << abs(keff-keff_previous) << '\n';
			cout << '\n';

			if ( residual < tol && abs(keff-keff_previous) < tol )
			//if ( residual < tol )
			{
				cout << '\n';
				cout << "Calculation converged with residual " << residual << '\n';
				cout << '\n';
				break;
			}
		}
		else
		{
			if ( residual < tol )
			{
				cout << '\n';
				cout << "Calculation converged with residual " << residual << '\n';
				cout << '\n';
				break;
			}
		}


	}

	output_write( Nx,Ny,Egrp,scalar_flux );

	cout << "Output file scalar_flux.out created in directory. " << '\n';
	cout << '\n';

}

#include <cmath>
#include <vector>
#include <iostream>

//*******************************************************************************************************************************

double kEigCalc( int Nx, int Ny, int Egrp, std::vector<std::vector<std::vector<double> > > scalar_flux, 
	             std::vector<std::vector<std::vector<double> > > scalar_flux_prev, 
	             std::vector<std::vector<std::vector<double> > > &nusigf, double keff)
{
	double sum_new = 0;
	double sum_prev = 0;

	for ( int e = 0; e < Egrp; e++ )
	{
		for ( int i = 0; i < Nx; i++ )
		{
			for ( int j = 0; j < Ny; j++ )
			{
				sum_prev = sum_prev + nusigf[i][j][e]*scalar_flux_prev[i][j][e];
				sum_new = sum_new + nusigf[i][j][e]*scalar_flux[i][j][e];
			}
		}
	}

	std::cout << "Sum New " << " " << sum_new << '\n';
	std::cout << "Sum Prev " << " " << sum_prev << '\n';

	return keff*(sum_new/sum_prev);
}
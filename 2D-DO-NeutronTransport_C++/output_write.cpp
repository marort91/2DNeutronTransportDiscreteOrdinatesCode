#include <iostream>
#include <fstream>
#include <vector>

void output_write( int Nx, int Ny, std::vector<std::vector<double> > &scalar_flux)
{
	std::ofstream outfile;
	outfile.open ( "scalar_flux.out");

	for ( int i = 0; i < Nx; i++ )
	{
		for ( int j = 0; j < Ny; j++ )
		{
			outfile << scalar_flux[i][j] << " ";

			if ( ( j + 1 ) % Nx == 0)
			{
				outfile << std::endl;
			}
		}
	}
}
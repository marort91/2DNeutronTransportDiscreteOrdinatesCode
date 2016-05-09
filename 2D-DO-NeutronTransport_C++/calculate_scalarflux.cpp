#include <vector>

//*******************************************************************************************************************************

void calculate_scalarflux( int Nx, int Ny, int ord, std::vector<std::vector<std::vector<double> > > &angular_flux,
	                             std::vector<std::vector<double> > &scalar_flux, std::vector<double> wi )
{
	double cumflux;
	
	for ( int i = 0; i < Nx; i++ )
	{
		for ( int j = 0; j < Ny; j++ )
		{
			cumflux = 0;

			for ( int k = 0; k < ord; k++ )
			{
				cumflux += 0.25*wi[k]*angular_flux[i][j][k];
			}

			scalar_flux[i][j] = cumflux;
		}
	}

}
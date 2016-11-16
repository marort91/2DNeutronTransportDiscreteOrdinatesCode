#include <cmath>
#include <vector>

//*******************************************************************************************************************************

double norm( int Nx, int Ny, int Egrp, std::vector<std::vector<std::vector<double> > > M1, std::vector<std::vector<std::vector<double> > > M2)
{
	double cum = 0;

	for ( int k = 0; k < Egrp; k++ )
	{
		for ( int i = 0; i < Nx; i++ )
		{
			for ( int j = 0; j < Ny; j++ )
			{
				cum += fabs( M1[i][j][k] - M2[i][j][k] ) * fabs( M1[i][j][k] - M2[i][j][k] );
			}
		}
	}

	return sqrt(cum);
}
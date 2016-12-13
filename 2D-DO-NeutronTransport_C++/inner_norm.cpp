#include <cmath>
#include <vector>

//*******************************************************************************************************************************

//double inner_norm( int Nx, int Ny, int Egrp, std::vector<std::vector<double> > M1, std::vector<std::vector<double> > M2)
//{
//	double cum = 0;
//
//		for ( int i = 0; i < Nx; i++ )
//		{
//			for ( int j = 0; j < Ny; j++ )
//			{
//				cum += fabs( M1[i][j] - M2[i][j] ) * fabs( M1[i][j] - M2[i][j] );
//			}
//		}
//
//	return sqrt(cum);
//}
double inner_norm( int Nx, int Ny, int Egrp, std::vector<std::vector<std::vector<double> > > M1, std::vector<std::vector<std::vector<double> > > M2)
{
	double cum = 0;

		for ( int i = 0; i < Nx; i++ )
		{
			for ( int j = 0; j < Ny; j++ )
			{
				cum += fabs( M1[i][j][Egrp] - M2[i][j][Egrp] ) * fabs( M1[i][j][Egrp] - M2[i][j][Egrp] );
			}
		}

	return sqrt(cum);
}
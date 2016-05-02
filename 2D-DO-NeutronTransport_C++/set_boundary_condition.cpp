#include <vector>
#include <string>

//*******************************************************************************************************************************

void set_boundary_condition( int bc, int Nx, int Ny, int ord, 
							 std::vector<std::vector<std::vector<double> > > &half_angular_flux_x,
	                         std::vector<std::vector<std::vector<double> > > &half_angular_flux_y )
{
	if ( bc == 1 )
	{
		for ( int i = 0; i < Nx+1; i++ )
		{
			for ( int j = 0; j < ord; j++ )
			{
				half_angular_flux_x[i][0][j] = 1.0;
			}
		}

		for ( int i = 0; i < Ny+1; i++ )
		{
			for ( int j = 0; j < ord; j++ )
			{
				half_angular_flux_y[0][i][j] = 1.0;
			}
		}
	}
}
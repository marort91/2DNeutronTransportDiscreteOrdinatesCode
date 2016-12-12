#include <vector>
#include <string>

//*******************************************************************************************************************************

void set_boundary_condition( int bc, int Nx, int Ny, int ord, int E, int k,
							 std::vector<double> &mu, std::vector<double> &eta,
							 std::vector<std::vector<std::vector<std::vector<double> > > > &half_angular_flux_x,
	                         std::vector<std::vector<std::vector<std::vector<double> > > > &half_angular_flux_y )
{
	int set_reflecting_BC( int ord, double muk, double etak, std::vector<double> &mu, std::vector<double> &eta );
	int reflect_idx;

	if ( bc == 1 )
	{
		for ( int e = 0; e < E; e++ )
		{
			for ( int i = 0; i < Nx; i++ )
			{
				for ( int j = 0; j < ord; j++ )
				{
					half_angular_flux_x[i][0][j][e] = 1.0;
				}
			}
		}

		for (int e = 0; e < E; e++ )
		{
			for ( int i = 0; i < Ny; i++ )
			{
				for ( int j = 0; j < ord; j++ )
				{
					half_angular_flux_y[0][i][j][e] = 1.0;
				}
			}
		}
	}
	else if ( bc == 2 || bc == 4 )
	{
		reflect_idx = set_reflecting_BC( ord, mu[k], eta[k], mu, eta );

		for ( int e = 0; e < E; e++ )
		{
			for ( int i = 0; i < Nx; i++ )
			{
				for ( int j = 0; j < Nx+1; j++ )
				{
					half_angular_flux_x[i][j][k][e] = half_angular_flux_x[i][j][reflect_idx][e];
				}
			}
		}
		for ( int e = 0; e < E; e++ )
		{
			for ( int i = 0; i < Ny+1; i++ )
			{
				for ( int j = 0; j < Ny; j++ )
				{
					half_angular_flux_y[i][j][k] = half_angular_flux_y[i][j][reflect_idx];
				}
			}
		}
	}
	else
	{
		for ( int e = 0; e < E; e++ )
		{
			for ( int i = 0; i < Nx; i++ )
			{
				for ( int j = 0; j < ord; j++ )
				{
					half_angular_flux_x[i][0][j][e] = 0.0;
				}
			}
		}

		for (int e = 0; e < E; e++ )
		{
			for ( int i = 0; i < Ny; i++ )
			{
				for ( int j = 0; j < ord; j++ )
				{
					half_angular_flux_y[0][i][j][e] = 0.0;
				}
			}
		}
	}
}

//*******************************************************************************************************************************

int set_reflecting_BC( int ord, double muk, double etak, std::vector<double> &mu, std::vector<double> &eta )
{
	int reflect_idx;

	if ( muk > 0 & etak > 0)
	{
		for ( int i = 0; i < ord; i++ )
		{
			if ( muk == -mu[i] & etak == -eta[i]) 
			{
				reflect_idx = i;
				break;
			}
		}
	}

	else if ( muk < 0 & etak > 0 )
	{
		for ( int i = 0; i < ord; i++ )
		{
			if ( muk == -mu[i] & etak == eta[i]) 
			{
				reflect_idx = i;
				break;
			}
		}
	}
	else if ( muk > 0 & etak < 0 )
	{
		for ( int i = 0; i < ord; i++ )
		{
			if ( muk == mu[i] & etak == -eta[i]) 
			{
				reflect_idx = i;
				break;
			}
		}
	}
	else if ( muk < 0 & etak < 0 )
	{
		for ( int i = 0; i < ord; i++ )
		{
			if ( muk == -mu[i] & etak == -eta[i]) 
			{
				reflect_idx = i;
				break;
			}
		}
	}
	
	return reflect_idx;
}
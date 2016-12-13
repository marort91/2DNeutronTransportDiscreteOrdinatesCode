#include <vector>
#include <iostream>

//*******************************************************************************************************************************

void src_extrn_scalarflux( std::vector<std::vector<std::vector<double> > > &S, 
						   std::vector<std::vector<std::vector<double> > > &Q,
                           std::vector<std::vector<std::vector<double> > > &scalar_flux,
						   std::vector<std::vector<double> > &sigs, int Nx, int Ny, int Egrp, double keff,
						   std::vector<double> &chi, std::vector<std::vector<std::vector<double> > > &nusigf, std::string calc_mode )
{

	double scatter_flux_temp;

	for ( int k = 0; k < Egrp; k++)
	{
		for ( int i = 0; i < Nx; i++ )
		{
			for ( int j = 0; j < Ny; j ++ )
			{
				scatter_flux_temp = 0;
				
				for ( int e = 0; e < Egrp; e++ )
				{
					if ( calc_mode != "Source" )
					{
						scatter_flux_temp = scatter_flux_temp + sigs[k][e]*scalar_flux[i][j][e] + (1/keff)*chi[e]*nusigf[i][j][e]*scalar_flux[i][j][e];	
					}
					else
					{
						scatter_flux_temp = scatter_flux_temp + sigs[k][e]*scalar_flux[i][j][e];	
					}
					
				}

				Q[i][j][k] = S[i][j][k] + scatter_flux_temp;
			}
		}
	}	
}
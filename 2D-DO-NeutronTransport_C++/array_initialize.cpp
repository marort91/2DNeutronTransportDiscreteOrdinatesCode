#include <vector>

//*******************************************************************************************************************************

void array_initialize( int Nx, int Ny, int ord, std::vector<std::vector<std::vector<double> > > &half_angular_flux_x, 
	                   std::vector<std::vector<std::vector<double> > > &half_angular_flux_y,
	                   std::vector<std::vector<std::vector<double> > > &angular_flux,
	                   std::vector<std::vector<double> > &scalar_flux,
	                   std::vector<std::vector<double> > &S, std::vector<std::vector<double> > &Q )
{

	// Set up three-dimensional x-interval half-angular flux.
	half_angular_flux_x.resize(Nx);
	for ( int i = 0; i < Nx; i++ )
	{
		half_angular_flux_x[i].resize(Nx+1);

		for ( int j = 0; j < ord; j++ )
			half_angular_flux_x[i][j].resize(ord);
	}

	// Set up three-dimensional y-interval half-angular flux.
	half_angular_flux_y.resize(Ny+1);
	for ( int i = 0; i < Ny + 1; i++ )
	{
		half_angular_flux_y[i].resize(Ny);

		for ( int j = 0; j < ord; j++ )
			half_angular_flux_y[i][j].resize(ord);
	}

	// Set up angular flux array.
	angular_flux.resize(Nx);
	for ( int i = 0; i < Ny; i++ )
	{
		angular_flux[i].resize(Ny);

		for ( int j = 0; j < ord; j++ )
			angular_flux[i][j].resize(ord);
	}

	// Set up scalar flux array.
	scalar_flux.resize(Nx);
	for ( int i = 0; i < Ny; i++ )
	{
		scalar_flux[i].resize(Ny);
	}

	// Set up source and scattering scalar flux arrays.
	S.resize(Nx);
	Q.resize(Ny);
	for ( int i = 0; i < Ny; i++ )
	{
		S[i].resize(Ny);
		Q[i].resize(Ny);
	}

}
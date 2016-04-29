#include <vector>

//*******************************************************************************************************************************

void spatial_discretize( double xL, double xR, int Nx, double dx, 
	                     double yB, double yT, int Ny, double dy, 
	                     std::vector<double> &x, std::vector<double> &y )
{
	//dx = (xR - xL)/Nx;
	//dy = (yT - yB)/Ny;

	x.resize(Nx);
	y.resize(Ny);

	for ( int i = 0; i < Nx; i++ )
	{
		x.at(i) = xL + dx/2 + i*dx;
	}

	for ( int i = 0; i < Ny; i++ )
	{
		y.at(i) = yB + dy/2 + i*dy;
	}
}
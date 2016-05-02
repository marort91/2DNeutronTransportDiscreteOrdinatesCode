#ifndef LEVEL_SYM_QUAD_H
#define LEVEL_SYM_QUAD_H

void level_sym_quad( int N, std::vector<double> &mu, std::vector<double> &eta, std::vector<double> &wi );
void array_initialize( int Nx, int Ny, int Nang, std::vector<std::vector<std::vector<double> > > &half_angular_flux_x, 
	                   std::vector<std::vector<std::vector<double> > > &half_angular_flux_y,
	                   std::vector<std::vector<std::vector<double> > > &angular_flux,
	                   std::vector<std::vector<double> > &scalar_flux, std::vector<std::vector<double> > &S, 
	                   std::vector<std::vector<double> > &Q );
void spatial_discretize( double xL, double xR, int Nx, double dx, 
	                     double yB, double yT, int Ny, double dy, 
	                     std::vector<double> &x, std::vector<double> &y );
double norm( int Nx, int Ny, std::vector<std::vector<double> > M1, std::vector<std::vector<double> > M2);
void set_boundary_condition( int bc, int Nx, int Ny, int ord, 
							 std::vector<std::vector<std::vector<double> > > &half_angular_flux_x,
	                         std::vector<std::vector<std::vector<double> > > &half_angular_flux_y );
void scatter_scalarflux_source( int Nx, int Ny, int ord, std::vector<std::vector<std::vector<double> > > &angular_flux,
	                            std::vector<std::vector<double> > &scalar_flux, std::vector<double> wi );
void const_external_source_def( std::vector<std::vector<double> > &Q, int Q0, int Nx, int Ny );

#endif 
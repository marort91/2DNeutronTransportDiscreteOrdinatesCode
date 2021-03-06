#ifndef LEVEL_SYM_QUAD_H
#define LEVEL_SYM_QUAD_H

void level_sym_quad( int N, std::vector<double> &mu, std::vector<double> &eta, std::vector<double> &wi );
void array_initialize( int Nx, int Ny, int ord, int Egrp, 
					   std::vector<std::vector<std::vector<std::vector<double> > > > &half_angular_flux_x, 
	                   std::vector<std::vector<std::vector<std::vector<double> > > > &half_angular_flux_y,
	                   std::vector<std::vector<std::vector<std::vector<double> > > > &angular_flux,
	                   std::vector<std::vector<std::vector<double> > > &scalar_flux,
	                   std::vector<std::vector<std::vector<double> > > &S, std::vector<std::vector<std::vector<double> > > &Q,
	                   std::vector<std::vector<std::vector<double> > > &sigt, std::vector<std::vector<double> > &sigs,
	                   std::vector<std::vector<std::vector<double> > > &nusigf,
	                   std::vector<double> &chi,
	                   std::vector<std::vector<double> > &group_scalar_flux );
void spatial_discretize( double xL, double xR, int Nx, double dx, 
	                     double yB, double yT, int Ny, double dy, 
	                     std::vector<double> &x, std::vector<double> &y );
double norm( int Nx, int Ny, int Egrp, std::vector<std::vector<std::vector<double> > > M1, std::vector<std::vector<std::vector<double> > > M2);
double inner_norm( int Nx, int Ny, int Egrp, std::vector<std::vector<std::vector<double> > > M1, std::vector<std::vector<std::vector<double> > > M2);
//double inner_norm( int Nx, int Ny, int Egrp, std::vector<std::vector<double> > M1, std::vector<std::vector<double> > M2);
void set_boundary_condition( int bc, int Nx, int Ny, int ord, int E, int k,
							 std::vector<double> &mu, std::vector<double> &eta,
							 std::vector<std::vector<std::vector<std::vector<double> > > > &half_angular_flux_x,
	                         std::vector<std::vector<std::vector<std::vector<double> > > > &half_angular_flux_y );
void calculate_scalarflux( int Nx, int Ny, int ord, int Egrp, 
						   std::vector<std::vector<std::vector<std::vector<double> > > > &angular_flux,
	                             std::vector<std::vector<std::vector<double> > > &scalar_flux, std::vector<double> wi );
//void const_external_source_def( std::vector<std::vector<double> > &Q, int Q0, int Nx, int Ny );
void source_external_scattering( std::vector<std::vector<std::vector<double> > > &Q,
                                 std::vector<std::vector<std::vector<double> > > &scalar_flux,
                                 std::vector<std::vector<double> > &sigs, int Q0, int Nx, int Ny, int Egrp );
void output_write( int Nx, int Ny, int Egrp, std::vector<std::vector<std::vector<double> > > &scalar_flux);
void input_read(int &N, int &Nx, double &xL, double &xR, int &Ny, double &yB, double &yT,
				int &bc, std::string &sigtfid, std::string &sigsfid, std::string &nusigffid, double &tol, 
				std::string &srcfid, int &Egrp, std::string &calc_mode, std::string &chifid);
void source_file_read( std::vector<std::vector<std::vector<double> > > &S, std::string fid, int Nx, int Ny, std::string calc_mode );
void src_extrn_scalarflux( std::vector<std::vector<std::vector<double> > > &S, 
						   std::vector<std::vector<std::vector<double> > > &Q,
                           std::vector<std::vector<std::vector<double> > > &scalar_flux,
						   std::vector<std::vector<double> > &sigs, int Nx, int Ny, int Egrp, double keff,
						   std::vector<double> &chi, std::vector<std::vector<std::vector<double> > > &nusigf,
						   std::string calc_mode );
void xs_file_read( std::vector<std::vector<std::vector<double> > > &sigt, std::string sigtfid, 
	               std::vector<std::vector<double> > &sigs, std::string sigsfid, 
	               std::vector<std::vector<std::vector<double> > > &nusigf, std::string nusigffid,
	               std::vector<double> &chi, std::string chifid,  int Nx, int Ny, int E );
int set_reflecting_BC( int ord, double muk, double etak, std::vector<double> &mu, std::vector<double> &eta );
void set_extern_src_zero( int Nx, int Ny, int Egrp, std::vector<std::vector<std::vector<double> > > &S );
void angular_flux_critical_guess( int Nx, int Ny, int ord, int Egrp, 
								  std::vector<std::vector<std::vector<std::vector<double> > > > &angular_flux );
double kEigCalc( int Nx, int Ny, int Egrp, std::vector<std::vector<std::vector<double> > > scalar_flux, 
	             std::vector<std::vector<std::vector<double> > > scalar_flux_prev, 
	             std::vector<std::vector<std::vector<double> > > &nusigf, double keff);
void calculate_group_scalarflux( int Nx, int Ny, int ord, int E, 
						   std::vector<std::vector<std::vector<std::vector<double> > > > &angular_flux,
	                             std::vector<std::vector<double> > &group_scalar_flux, std::vector<double> wi );
void group_src( std::vector<std::vector<std::vector<double> > > &S, 
						   std::vector<std::vector<std::vector<double> > > &Q,
						   std::vector<std::vector<std::vector<double> > > &scalar_flux,
                           std::vector<std::vector<double> >  &group_scalar_flux,
						   std::vector<std::vector<double> > &sigs, int Nx, int Ny, int E, double keff,
						   std::vector<double> &chi, std::vector<std::vector<std::vector<double> > > &nusigf, int Egrp,
						   std::string calc_mode );
void scalar_flux_critical_guess( int Nx, int Ny, int Egrp, 
								  std::vector<std::vector<std::vector<double> > > &scalar_flux );

#endif 
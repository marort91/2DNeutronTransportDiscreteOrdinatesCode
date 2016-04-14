%% Two-dimension, One-speed Discrete Ordinates Neutron Transport Equation Solver
% Solves the 2D neutron transport code using the discrete ordinates method.

clc, clear, clf

%% Material Properies
sigt = 1.0;
sigs0 = 0.0;
siga = sigt - sigs0;

%% Geometry and Angular Discretization
xL = 0; xR = 1;
yB = 0; yT = 1;

Nx = 10; Ny = 10;
Nang = 16;

[mu,eta,wi] = level_sym_table(Nang);

dx = (xR - xL)/Nx;
dy = (yT - yB)/Ny;

%% Generation of Angular Flux Array
angular_flux = zeros(Nx,Ny,Nang*(Nang+2)/2);

%% Source Generation
S = 1.*ones(Nx,Ny);
Q = zeros(Nx,Ny);

%% Calculation Parameters
maxiter = 1e3; %Maximum number of sweeps allowed
tol = 1e-8;

%% Boundary Conditions
bc = 1;

%% Code Parameters to be Printed to Terminal
fprintf('S%i calculation with scattering cross section = %f \n',Nang,sigs0);
fprintf('Maximum number of iterations: %i \n',maxiter);
fprintf('Scalar flux two-norm tolerance: %f \n',tol);
fprintf('\n');

%% Transport Sweep
for iter = 1:maxiter
    
    angular_flux_prev = angular_flux;
    
    scalar_flux = zeros(Nx,Ny);
    
    for k = 1:Nang*(Nang+2)/2
        
        scalar_flux = scalar_flux + 0.25.*wi(k).*angular_flux(:,:,k);
        
    end
    
    Q = ( S + sigs0.*scalar_flux );
    
    for l = 1:Nang*(Nang+2)/2
        
        if ( bc == 1 )
        
            angular_flux_half_x = zeros(Nx,Nx+1);
            angular_flux_half_y = zeros(Ny+1,Ny);
            
        end
        
        if ( mu(l) > 0 && eta(l) > 0 )
            
            for j = 1:Ny
                
                for i = 1:Nx
                    
                    angular_flux(i,j,l) = ( 2*mu(l)*angular_flux_half_x(j,i)/dx + ...
                        2*eta(l)*angular_flux_half_y(j,i)/dy + Q(i,j) )/...
                        ( 2*mu(l)/dx + 2*eta(l)/dy + sigt );
                    
                    angular_flux_half_x(j,i+1) = 2*angular_flux(i,j,l) - angular_flux_half_x(j,i);
                    
                end
                
                for m = 1:Nx
                    
                    angular_flux_half_y(j+1,m) = 2*angular_flux(m,j,l) - angular_flux_half_y(j,m);
                    
                end
                
            end
            
        elseif ( mu(l) < 0 && eta(l) > 0 )
            
            for j = 1:Ny
                
                for i = Nx:-1:1
            
                    angular_flux(i,j,l) = ( -2*mu(l)*angular_flux_half_x(j,i+1)/dx + ...
                        2*eta(l)*angular_flux_half_y(j,i)/dy + Q(i,j) )/...
                        ( -2*mu(l)/dx + 2*eta(l)/dy + sigt );
                    
                    angular_flux_half_x(j,i) = 2*angular_flux(i,j,l) - angular_flux_half_x(j,i+1);
                    
                end
                
                for m = Nx:-1:1
                    
                    angular_flux_half_y(j+1,m) = 2*angular_flux(m,j,l) - angular_flux_half_y(j,m);
                    
                end
                
            end
                    
        elseif ( mu(l) > 0 && eta(l) < 0 )
            
            for j = Ny:-1:1
                
                for i = 1:Nx
                    
                    angular_flux(i,j,l) = ( 2*mu(l)*angular_flux_half_x(j,i)/dx - ...
                        2*eta(l)*angular_flux_half_y(j+1,i)/dy + Q(i,j) )/...
                        ( 2*mu(l)/dx - 2*eta(l)/dy + sigt );
                    
                    angular_flux_half_x(j,i+1) = 2*angular_flux(i,j,l) - angular_flux_half_x(j,i);
                    
                end
                
                for m = 1:Nx
                    
                    angular_flux_half_y(j,m) = 2*angular_flux(m,j,l) - angular_flux_half_y(j+1,m);
                    
                end
                
            end
            
        elseif ( mu(l) < 0 && eta(l) < 0 )
            
            for j = Ny:-1:1
                
                for i = Nx:-1:1
                    
                    angular_flux(i,j,l) = ( -2*mu(l)*angular_flux_half_x(j,i+1)/dx - ...
                        2*eta(l)*angular_flux_half_y(j+1,i)/dy + Q(i,j) )/...
                        ( -2*mu(l)/dx - 2*eta(l)/dy + sigt );
                    
                    angular_flux_half_x(j,i) = 2*angular_flux(i,j,l) - angular_flux_half_x(j,i+1);
                    
                end
                
                for m = 1:Nx
                    
                    angular_flux_half_y(j,m) = 2*angular_flux(m,j,l) - angular_flux_half_y(j+1,m);
                    
                end
                
            end
                    
            
        else
            
            error('Angular discretization incorrect \n')
            
        end
        
    end
    
    scalar_flux_prev = zeros(Nx,Ny);
    scalar_flux_new = zeros(Nx,Ny);
    
    for k = 1:Nang*(Nang+2)/2
        
        scalar_flux_prev = scalar_flux_prev + 0.25.*wi(k).*angular_flux_prev(:,:,k);
        scalar_flux_new = scalar_flux_new + 0.25.*wi(k).*angular_flux(:,:,k);
        
    end
    
    residual = norm(scalar_flux_new-scalar_flux_prev);
    fprintf('Residual: %f     Iteration: %i \n',residual,iter-1);
    
    if ( residual < tol )
        
        break
        
    end
       
end
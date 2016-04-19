%% Two-dimension, One-speed Discrete Ordinates Neutron Transport Equation Solver
% Solves the 2D neutron transport code using the discrete ordinates method.

clc, clear, clf, close all

%% Material Properies
sigt = 1.0;
sigs0 = 0.99999;
siga = sigt - sigs0;

%% Geometry and Angular Discretization
xL = 0; xR = 1;
yB = 0; yT = 1;

Nx = 1000; Ny = 1000;
Nang = 16;

[mu,eta,wi] = level_sym_table(Nang);

dx = (xR - xL)/(Nx);
dy = (yT - yB)/(Ny);

x = xL+dx/2:dx:xR-dx/2;
y = yB+dy/2:dy:yT-dy/2;

[Xarr,Yarr] = meshgrid(x,y);

%% Generation of Angular Flux Array
angular_flux = zeros(Nx,Ny,Nang*(Nang+2)/2);

%% Source Generation
S = 1.*ones(Nx,Ny);
Q = zeros(Nx,Ny);

%% Calculation Parameters
maxiter = 1e3; %Maximum number of sweeps allowed
tol = 1e-8;

%% Boundary Conditions
bc = 3;
%bc = 'Larsen2D-Benchmark';

if ( bc == 1 )
    
    fprintf('Vacuum Boundary Conditions \n');
    
elseif ( strcmp(bc,'Larsen2D-Benchmark') == 1 )
    
    fprintf('Larsen 2D-benchmark testing, unit angular fluxes imposed on left and bottom boundaries \n');
    
elseif ( bc == 3 )
    
    fprintf('Reflective boundary conditions on right and top boundaries \n');
    
end

%% Code Parameters to be Printed to Terminal
fprintf('S%i calculation with scattering cross section = %f \n',Nang,sigs0);
fprintf('Maximum number of iterations: %i \n',maxiter);
fprintf('Scalar flux two-norm tolerance: %e \n',tol);
fprintf('\n');

%% Transport Sweep
for iter = 1:maxiter
    
    angular_flux_prev = angular_flux;
    
    scalar_flux = zeros(Nx,Ny);
    
    for k = 1:Nang*(Nang+2)/2
        
        scalar_flux = scalar_flux + 0.25.*wi(k).*angular_flux(:,:,k);
        
    end
    
    Q = ( S + sigs0.*scalar_flux );
    
    if ( strcmp(bc,'Larsen2D-Benchmark') == 1 )
         
        angular_flux_half_x = zeros(Nx,Nx+1,Nang*(Nang+2)/2);
        angular_flux_half_x(:,1,:) = 1;
        angular_flux_half_y = zeros(Ny+1,Ny,Nang*(Nang+2)/2);
        angular_flux_half_y(1,:,:) = 1;
            
    elseif ( bc == 1 || bc == 3 )
                
        angular_flux_half_x = zeros(Nx,Nx+1,Nang*(Nang+2)/2);
        angular_flux_half_y = zeros(Ny+1,Ny,Nang*(Nang+2)/2);
                
    else
                
        error('No boundary condition implemented');
                
    end
    
    for l = 1:Nang*(Nang+2)/2
        
        if ( mu(l) > 0 && eta(l) > 0 )
            
            for j = 1:Ny
                
                for i = 1:Nx
                    
                    angular_flux(j,i,l) = ( 2*mu(l)*angular_flux_half_x(j,i,l)/dx + ...
                        2*eta(l)*angular_flux_half_y(j,i,l)/dy + Q(i,j) )/...
                        ( 2*mu(l)/dx + 2*eta(l)/dy + sigt );
                    
                    angular_flux_half_x(j,i+1,l) = 2*angular_flux(j,i,l) - angular_flux_half_x(j,i,l);
                    
                end
                
                for m = 1:Nx
                    
                    angular_flux_half_y(j+1,m,l) = 2*angular_flux(j,m,l) - angular_flux_half_y(j,m,l);
                    
                end
                
            end
            
        elseif ( mu(l) < 0 && eta(l) > 0 )
            
            if ( bc == 3 )
                
                loc = find ( mu(l) == -mu & eta(l) == eta );
                
                angular_flux_half_x(:,:,l) = angular_flux_half_x(:,:,loc);
                
            end
            
            for j = 1:Ny
                
                for i = Nx:-1:1
            
                    angular_flux(j,i,l) = ( -2*mu(l)*angular_flux_half_x(j,i+1,l)/dx + ...
                        2*eta(l)*angular_flux_half_y(j,i,l)/dy + Q(i,j) )/...
                        ( -2*mu(l)/dx + 2*eta(l)/dy + sigt );
                    
                    angular_flux_half_x(j,i,l) = 2*angular_flux(j,i,l) - angular_flux_half_x(j,i+1,l);
                    
                end
                
                for m = Nx:-1:1
                    
                    angular_flux_half_y(j+1,m,l) = 2*angular_flux(j,m,l) - angular_flux_half_y(j,m,l);
                    
                end
                
            end
                    
        elseif ( mu(l) > 0 && eta(l) < 0 )
            
            if ( bc == 3 )
                
                loc = find( eta(l) == -eta & mu(l) == mu );
                
                angular_flux_half_y(:,:,l) = angular_flux_half_y(:,:,loc);
                
            end
            
            for j = Ny:-1:1
                
                for i = 1:Nx
                    
                    angular_flux(j,i,l) = ( 2*mu(l)*angular_flux_half_x(j,i,l)/dx - ...
                        2*eta(l)*angular_flux_half_y(j+1,i,l)/dy + Q(i,j) )/...
                        ( 2*mu(l)/dx - 2*eta(l)/dy + sigt );
                    
                    angular_flux_half_x(j,i+1,l) = 2*angular_flux(j,i,l) - angular_flux_half_x(j,i,l);
                    
                end
                
                for m = 1:Nx
                    
                    angular_flux_half_y(j,m,l) = 2*angular_flux(j,m,l) - angular_flux_half_y(j+1,m,l);
                    
                end
                
            end
            
        elseif ( mu(l) < 0 && eta(l) < 0 )
            
            if ( bc == 3 )
                
                loc = find ( mu(l) == -mu & eta(l) == -eta );
                
                angular_flux_half_x(:,:,l) = angular_flux_half_x(:,:,loc);
                angular_flux_half_y(:,:,l) = angular_flux_half_y(:,:,loc);
                
            end
            
            for j = Ny:-1:1
                
                for i = Nx:-1:1
                    
                    angular_flux(j,i,l) = ( -2*mu(l)*angular_flux_half_x(j,i+1,l)/dx - ...
                        2*eta(l)*angular_flux_half_y(j+1,i,l)/dy + Q(i,j) )/...
                        ( -2*mu(l)/dx - 2*eta(l)/dy + sigt );
                    
                    angular_flux_half_x(j,i,l) = 2*angular_flux(j,i,l) - angular_flux_half_x(j,i+1,l);
                    
                end
                
                for m = 1:Nx
                    
                    angular_flux_half_y(j,m,l) = 2*angular_flux(j,m,l) - angular_flux_half_y(j+1,m,l);
                    
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
    fprintf('Residual: %e     Iteration: %i \n',residual,iter-1);
    
    if ( residual < tol )
        
        break
        
    end
       
end

int = 1;

if ( strcmp(bc,'Larsen2D-Benchmark') == 1 )

    f = @(x,y) exp(-sigt*min(x/mu(int),y/eta(int)));
    Z = f(Xarr,Yarr);

    figure(1)
    mesh(x,y,Z);

    figure(2)
    mesh(x,y,angular_flux(:,:,int));

    max(max(abs(Z-angular_flux(:,:,int))))
    
else
    
    figure(1)
    mesh(x,y,scalar_flux);
    xlabel('x-direction');
    ylabel('y-direction');
    
end
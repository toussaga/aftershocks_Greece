function [P]=p2d(P_old, Beta, k_mu, Q, nx, ny, dx, dy, dt, Bound_type_hor, Bound_val_hor, Bound_type_ver, Bound_val_ver)

%   1999-2002    Boris Kaus
%   2015         Gunnar Jansen & Stephen A. Miller 
%   2022         Last modified: Thanushika Gunatilake

%Transpose input arrays
P_old = P_old';
k_mu  = k_mu';
Beta  = Beta'; 
Q     = Q';

%Initialise working arrays:
A                   =  sparse(nx*ny,nx*ny);
aa                  =  zeros(nx*ny,1);bb=aa;cc=aa;dd=aa;ee=aa;
P_vec               =  zeros(nx*ny,1);
k_mu_vec            =  zeros(nx*ny,1);
rhs                 =  zeros(nx*ny,1); %1 = spalte
P                   =  zeros(ny,nx);

% Check the timestep and give a warning if the timestep is too large
kappa               =  k_mu./Beta;
dt_opt              = (2*sqrt(dx^2+dy^2))^2/(4*pi^2*max(max(kappa)));

k_mu_vec               =  k_mu(:);
Q_vec                  =  Q(:);
beta_vec               =  Beta(:);

rhs                     =   -P_old(:) - Q_vec.*dt./beta_vec;

ii                      =   [ny+1:ny*nx-ny];
cc(1:end)=0;cc(ii)      =   dt./beta_vec(ii).*1/2.*(k_mu_vec(ii+1 )+k_mu_vec(ii))./dy^2;
aa(1:end)=0;aa(ii)      =   dt./beta_vec(ii).*1/2.*(k_mu_vec(ii-1 )+k_mu_vec(ii))./dy^2;
dd(1:end)=0;dd(ii)      =   dt./beta_vec(ii).*1/2.*(k_mu_vec(ii-ny)+k_mu_vec(ii))./dx^2;
ee(1:end)=0;ee(ii)      =   dt./beta_vec(ii).*1/2.*(k_mu_vec(ii+ny)+k_mu_vec(ii))./dx^2;
bb(1:end)=1;bb(ii)      =   -1 -aa(ii)-cc(ii)-dd(ii)-ee(ii);

if      strcmp(lower(Bound_type_hor{1}),    'dirichlet' )       % constant head left boundary
    ii                  =   1:ny;   bb(ii)=1; aa(ii)=0; cc(ii)=0; dd(ii)=0; ee(ii)=0; 
    rhs(ii)             =   Bound_val_hor(:,1);     
elseif  strcmp(lower(Bound_type_hor{1}),    'neumann'   )       % constant flux left boundary
    ii                  =   1:ny;   bb(ii)=-1/dx; aa(ii)=0; cc(ii)=0; dd(ii)=0; ee(ii)=1/dx;
    rhs(ii)             =   Bound_val_hor(:,1); 
elseif  strcmp(lower(Bound_type_hor{1}),    'periodic'   )       % periodic left boundary: watch out it has to be set after the matrix is formed!!!
    ii                  =   1:ny;   bb(ii)=-1/dx; aa(ii)=0; cc(ii)=0; dd(ii)=0; ee(ii)=1/dx;
else
    error('The left boundary condition you specified is not implemented yet')
end


if      strcmp(lower(Bound_type_hor{2}),    'dirichlet' )       % constant head right boundary
    ii                  =   nx*ny-ny:nx*ny;   bb(ii)=1; aa(ii)=0; cc(ii)=0; dd(ii)=0; ee(ii)=0; 
    rhs(ii)             =   Bound_val_hor(:,2); 
elseif  strcmp(lower(Bound_type_hor{2}),    'neumann'   )       % constant flux right boundary
    ii                  =   nx*ny-ny:nx*ny;   bb(ii)= 1/dx; aa(ii)=0; cc(ii)=0; dd(ii)=-1/dx; ee(ii)=0;
    rhs(ii)             =   Bound_val_hor(:,2); 
elseif  strcmp(lower(Bound_type_hor{2}),    'periodic'   )       % periodic right boundary: watch out it has to be set after the matrix is formed!!!
     ii                  =   nx*ny-ny:nx*ny;    bb(ii)=-1/dx; aa(ii)=0; cc(ii)=0; dd(ii)=0; ee(ii)=1/dx;
else
    error('The right boundary condition you specified is not implemented yet')
end


if      strcmp(lower(Bound_type_ver{1}),    'dirichlet' )       % constant head lower boundary
    ii                  =   1:ny:nx*ny;    bb(ii)=1; aa(ii)=0; cc(ii)=0; dd(ii)=0; ee(ii)=0; 
    rhs(ii)             =   Bound_val_ver(:,1); 
elseif  strcmp(lower(Bound_type_ver{1}),    'neumann'   )       % constant flux lower boundary
    ii                  =   1:ny:nx*ny;    bb(ii)=-1/dy; aa(ii)=0; cc(ii)= 1/dy; dd(ii)=0; ee(ii)=0;
    rhs(ii)             =   Bound_val_ver(:,1); 
else
    error('The lower boundary condition you specified is not implemented yet')
end


if      strcmp(lower(Bound_type_ver{2}),    'dirichlet' )       % constant head upper boundary
    ii                  =   ny:ny:nx*ny;   bb(ii)=1; aa(ii)=0; cc(ii)=0; dd(ii)=0; ee(ii)=0; 
    rhs(ii)             =   Bound_val_ver(:,2); 
elseif  strcmp(lower(Bound_type_ver{2}),    'neumann'   )       % constant flux upper boundary
    ii                  =   ny:ny:nx*ny;   bb(ii)= 1/dy; aa(ii)=-1/dy; cc(ii)=0; dd(ii)=0; ee(ii)=0;
    rhs(ii)             =   Bound_val_ver(:,2); 
else
    error('The upper boundary condition you specified is not implemented yet')
end

A                   =   spdiags([ee cc bb aa dd], [-ny -1:1 ny], nx*ny, nx*ny);
A                   =   A';

% If there ar periodic horizontal BC's, we have to set them here
if  strcmp(lower(Bound_type_hor{1}),    'periodic'   )    
    ii                      =   1:ny;
    for i=1:length(ii)
        A(ii(i),:) = 0;
    end
    
    for i=1:length(ii)
        ii_normal   = ii(i);
        ii_periodic = ny*(nx-1)+ii(i);
        A(ii_normal,ii_normal   + 1 )  =   dt./beta_vec(ii_normal).*1/2.*(k_mu_vec(ii_normal   + 1 )  + k_mu_vec(ii_normal))./dy^2;
        A(ii_normal,ii_periodic - 1 )  =   dt./beta_vec(ii_normal).*1/2.*(k_mu_vec(ii_periodic - 1 )  + k_mu_vec(ii_normal))./dy^2;
        A(ii_normal,ii_periodic - ny)  =   dt./beta_vec(ii_normal).*1/2.*(k_mu_vec(ii_periodic - ny)  + k_mu_vec(ii_normal))./dx^2;
        A(ii_normal,ii_normal   + ny)  =   dt./beta_vec(ii_normal).*1/2.*(k_mu_vec(ii_normal   + ny)  + k_mu_vec(ii_normal))./dx^2;
        
        A(ii_normal,ii_normal)         =   -1 - A(ii_normal,ii_normal+1) -  A(ii_normal,ii_periodic-1) - A(ii_normal,ii_periodic-ny) - A(ii_normal,ii_normal+ny);        
        A(ii_normal,:) = -A(ii_normal,:);
    end
end

if  strcmp(lower(Bound_type_hor{2}),    'periodic'   )    
    ii                      =   nx*ny-ny+1:nx*ny-1;
    for i=1:length(ii)
        A(ii(i),:) = 0;
    end
    
    for i=1:length(ii)
        ii_normal   = ii(i);
        ii_periodic = ii(i) - ny*(nx-1);
        
        A(ii_normal,ii_periodic + 1 ) =   dt./beta_vec(ii_normal).*1/2.*(k_mu_vec(ii_periodic +1 )   + k_mu_vec(ii_normal))./dy^2;
        A(ii_normal,ii_normal   - 1 ) =   dt./beta_vec(ii_normal).*1/2.*(k_mu_vec(ii_normal   -1 )   + k_mu_vec(ii_normal))./dy^2;
        A(ii_normal,ii_normal   - ny) =   dt./beta_vec(ii_normal).*1/2.*(k_mu_vec(ii_normal   -ny)   + k_mu_vec(ii_normal))./dx^2;
        A(ii_normal,ii_periodic + ny) =   dt./beta_vec(ii_normal).*1/2.*(k_mu_vec(ii_periodic +ny)   + k_mu_vec(ii_normal))./dx^2;
        A(ii_normal,ii_normal       ) =   -1 - A(ii_normal,ii_periodic+1) -  A(ii_normal,ii_normal-1) - A(ii_normal,ii_normal-ny) - A(ii_normal,ii_periodic+ny);


    end
end

P_vec = A\rhs;

P(find(P_vec)) = P_vec(find(P_vec)); % reshape()
P=P';

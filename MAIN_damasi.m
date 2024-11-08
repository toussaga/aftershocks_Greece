%% Clear version of Pressure Diffusion in 2D - Damasi, Greece
%% Gaëlle Toussaint (based on Thanushika Gunatilake & Stephen A. Miller) 
%% Require 3D-Simulation-Visualization MATLAB package (https://github.com/TerdikGyorgy/3D-Simulation-Visualization)

clear all
close all

%% Load Data
Damasi_data = readtable('Kassaras_sup/1-s2.0-S0264370722000023-mmc5.xlsx');
Mcompl = find(Damasi_data.Magnitude >= 2.3);
Damasi_data = Damasi_data(Mcompl,:);
save('Damasi_data.mat','Damasi_data');
load('Damasi_data.mat') 
Lat = Damasi_data.Latitude;
Long= Damasi_data.Longitude;
M = Damasi_data.Magnitude;
depth = Damasi_data.Depth;
yr = datetime(Damasi_data.Year,Damasi_data.Month,Damasi_data.Day,Damasi_data.Hour,Damasi_data.Min, Damasi_data.Sec);
absoluteDays = days(yr-yr(1));
relativeDays = diff(yr);
relativeDays=[0;days(relativeDays)];
day = relativeDays;
data_days_check_plot=absoluteDays;

data_subset_damasi;
subDam4 = Damasi_data(af4,:);
%% Domain parameters
L                       = 20*1000;      %   Length of domain  [m]
H                       = 15*1000;      %   Height of domain  [m]
nx                      = 400;          %   Number of gridpoints in x-direction
ny                      = 300;          %   Number of gridpoints in y-direction
dy                      = H/ny;         %   Spacing between gridpoints in y-direction
dx                      = L/nx;         %   Spacing between gridpoints in x-direction
[x2d,y2d]               = meshgrid(0:dx:(nx-1)*dx,(ny-1)*dy:-dy:0); %Computational domain
x2d                     = x2d'./1000;
y2d                     = y2d'./1000;
material                = ones(nx,ny);  %   Model geometry
beta_fluid              = 1e-10;        %   Fluid compressibility       [1/Pa]
beta_phi                = 1e-8;         %   Pore compressibility        [1/Pa]
porosity_rock           = 0.03;         %   Rock porosity
rho_fluid               = 1000;         %   Density of fluid            [kg/m3]
rho_rock                = 2700;         %   Density of rocks            [kg/m3]
g                       = 10;           %   Gravitational acceleration  [m2/s]
Geotherm_grad           = 30;           %   Geothermal gradient         [K/km]
Stress_state            = 2;            %   Stress state: 1-Extension, 2-Compression, 3-Transtension
num_faults              = 1;            %   Specify the number of discrete faults in the model (current maximum 5)
Dir_name                = '2D_output_1';%   Directory where the data should be saved
Day                     = 3600*24;      %   Number of seconds per day
dt                      = Day;          %   Timestep in seconds
tend                    = 58;          %   Overall sumulated time in days (%tend=1000,Save_timestep=2000;)
tot_time                = (tend*Day);   %   Total nunmber of timesteps
Num_timestep            = 0;            %   number of time steps simulated
Save_timestep           = 5;            % 0.01*2000=20Bilder, How many timesteps are there between 2 savings?  0.05*20 = 1Bild pro Tag  
Save_data               = true;         %   Save the data or not?
Save_tiff_pictures      = true;         %   Save pictures as tiff?
error_level             = 1e-2;
error_iterations        = realmax;
n_iter                  = 0;

%% Faults
Faults(1).x0            = 12400;        %   x-coordinate of faults upper limit
Faults(1).y0            = -4950;        %   y-coordinate of faults upper limit
Faults(1).yend          = -8500;        %   y-coordinate of faults lower limit (x-coordinate is calculated)
Faults(1).thickness     = 0.1*1000;     %   Thickness of the fault [m]
Faults(1).dip           = 45.3;         %   Dip of the fault. Measured anti-clock wise from the x-axis [deg]
Faults(1).azimuth       = 342.6;        %   Azimuth of the fault. Measured clock wise from the y-axis [deg]

Faults(2).x0            = 6000;         %   x-coordinate of faults upper limit
Faults(2).y0            = -7700;        %   y-coordinate of faults upper limit
Faults(2).yend          = -12000;       %   y-coordinate of faults lower limit (x-coordinate is calculated)
Faults(2).thickness     = 0.1*1000;     %   Thickness of the fault [m]
Faults(2).dip           = 30;           %   Dip of the fault. Measured anti-clock wise from the x-axis [deg]
Faults(2).azimuth       = 54.9;         %   Azimuth of the fault. Measured clock wise from the y-axis [deg]

Faults(3).x0            = 7500;         %   x-coordinate of faults upper limit
Faults(3).y0            = -5000;        %   y-coordinate of faults upper limit
Faults(3).yend          = -13000;       %   y-coordinate of faults lower limit (x-coordinate is calculated)
Faults(3).thickness     = 0.1*1000;     %   Thickness of the fault [m]
Faults(3).dip           = 38;           %   Dip of the fault. Measured anti-clock wise from the x-axis [deg]
Faults(3).azimuth       = 32;           %   Azimuth of the fault. Measured clock wise from the y-axis [deg]

Faults(4).x0            = 15700;        %   x-coordinate of faults upper limit
Faults(4).y0            = -6400;        %   y-coordinate of faults upper limit
Faults(4).yend          = -10500;       %   y-coordinate of faults lower limit (x-coordinate is calculated)
Faults(4).thickness     = 0.1*1000;     %   Thickness of the fault [m]
Faults(4).dip           = 40.5;         %   Dip of the fault. Measured anti-clock wise from the x-axis [deg]
Faults(4).azimuth       = 204.2;        %   Azimuth of the fault. Measured clock wise from the y-axis [deg]

Faults(5).x0            = 9350;         %   x-coordinate of faults upper limit
Faults(5).y0            = -6050;        %   y-coordinate of faults upper limit
Faults(5).yend          = -7900;        %   y-coordinate of faults lower limit (x-coordinate is calculated)
Faults(5).thickness     = 0.1*1000;     %   Thickness of the fault [m]
Faults(5).dip           = 23.1;         %   Dip of the fault. Measured anti-clock wise from the x-axis [deg]
Faults(5).azimuth       = 54.8;         %   Azimuth of the fault. Measured clock wise from the y-axis [deg]

Faults(6).x0            = 12650;        %   x-coordinate of faults upper limit
Faults(6).y0            = -5250;        %   y-coordinate of faults upper limit
Faults(6).yend          = -8500;        %   y-coordinate of faults lower limit (x-coordinate is calculated)
Faults(6).thickness     = 0.1*1000;     %   Thickness of the fault [m]
Faults(6).dip           = 40;           %   Dip of the fault. Measured anti-clock wise from the x-axis [deg]
Faults(6).azimuth       = 347.6;     %   Azimuth of the fault. Measured clock wise from the y-axis [deg]

%% Boundary Conditions
Bound_type_hor          =   {'Neumann'  ,'Neumann'  };
Bound_val_hor           =   [0   , 0];
Bound_type_ver          =   {'Neumann'  ,'Dirichlet'};
Bound_val_ver           =   [0              0];

%% Generates random fault planes on the domain
n=nx*ny; 
kappa = 0.2;
gamm = 0; 
Mu=[0 0 1];

Y = Random_FB4(kappa,gamm,Mu,n);
r = sqrt(Y(:,1).^2+Y(:,2).^2+Y(:,3).^2);
theta = atan2d(Y(:,2),Y(:,1));
phi = acosd(Y(:,3)./r)-90;
for i=1:length(theta)
    if theta(i)<0
        theta(i)=theta(i)+360;
    end
end
planes = [theta phi];

%% Material
Material_mat(1,:)       =  [   60       0       0.01    1e-14      35    0.1];      % Matrix
Material_mat(2,:)       =  [   60       0       0.05    1e-11      35    0.1];      % Fault 1 

P_input = [0 3.1e-5 2.2e-5 2.2e-5 10e-5]; 

mat_matrix=1;
mat_fault2=2;
mat_fault3=3;
mat_fault4=4;
mat_fault5=5;

%% Create Domain
for i=1:nx
    for j=1:ny
        material(i,j)=mat_matrix;
        P_matrix(i,j)=mat_matrix;
    end
end

%% Set fault indentifier, point source assignment and stress state
for k=1:num_faults
  if Faults(k).y0 <= Faults(k).yend
   error(['The end of fault %i should be deeper than the start!', k])
  end
  x0 = Faults(k).x0;
  y0 = Faults(k).y0;
  if Faults(k).azimuth<=180
      xend = x0-abs(Faults(k).yend-Faults(k).y0) / tan(-Faults(k).dip*pi/180.0);
  else
      xend = x0-abs(Faults(k).yend-Faults(k).y0) / tan(Faults(k).dip*pi/180.0);
  end
  Faults(k).xend = xend;
  yend = Faults(k).yend;
  for i=1:nx
   for j=1:ny  
       xpos = i*dx;
       ypos = -(ny-j)*dy;
       dist_to_fault = abs((yend-y0)*xpos-(xend-x0)*ypos+xend*y0-yend*x0)/sqrt((yend-y0)^2+(xend-x0)^2);
       if (dist_to_fault <= Faults(k).thickness && ypos >= Faults(k).yend && ypos <= Faults(k).y0)
            material(i,j)=1+k;
       end
   end
  end
end

% Assign Q0 to point where main aftershocks are observed
 for i=2:nx-1
   for j=2:ny-1  
       if ((p4_final(1,2)-cs1_final(4,2))*1000 >= (i-1)*dx && (p4_final(1,2)-cs1_final(4,2))*1000 <= (i)*dx && -depth(af4(1))*1000 >= -(ny-(j))*dy && -depth(af4(1))*1000 <= -(ny-(j+1))*dy)
           P_matrix(i,j)=2;
       end
   end
 end

 for i=2:nx-1
   for j=2:ny-1  
       if ((p4_final(66,2)-cs1_final(4,2))*1000 >= (i-1)*dx && (p4_final(66,2)-cs1_final(4,2))*1000 <= (i)*dx && -depth(af4(66))*1000 >= -(ny-(j))*dy && -depth(af4(66))*1000 <= -(ny-(j+1))*dy)
           P_matrix(i,j)=3;
       end
   end
 end

 for i=2:nx-1
   for j=2:ny-1  
       if ((p4_final(112,2)-cs1_final(4,2))*1000 >= (i-1)*dx && (p4_final(112,2)-cs1_final(4,2))*1000 <= (i)*dx && -depth(af4(112))*1000 >= -(ny-(j))*dy && -depth(af4(112))*1000 <= -(ny-(j+1))*dy)
           P_matrix(i,j)=4;
       end
   end
 end

 for i=2:nx-1
   for j=2:ny-1  
       if ((p4_final(135,2)-cs1_final(4,2))*1000 >= (i-1)*dx && (p4_final(135,2)-cs1_final(4,2))*1000 <= (i)*dx && -depth(af4(135))*1000 >= -(ny-(j))*dy && -depth(af4(135))*1000 <= -(ny-(j+1))*dy)
           P_matrix(i,j)=5;
       end
   end
 end
%% Stress state
P_fluid_hydro   =   rho_fluid*g*(y2d+.1)*1000;  % Hydrostatic fluid pressure
Sigma_1         =   rho_rock*g*(y2d+.1)*1000;   % Lithostatic rock pressure
phi = 0.64;
a1 = 0.58;

if Stress_state==1                              
    Sigma_3 = Sigma_1;
    Sigma_1 = 2.*Sigma_3;                      
elseif Stress_state==2                         
    Sigma_3 = a1*Sigma_1;
    Sigma_2 = (phi+(1-phi)*a1)*Sigma_1;
    Sigma_1 = Sigma_1;  
elseif Stress_state==3                        
    Sigma_3 = Sigma_1;
    Sigma_1 = 1.5*Sigma_3;                                        
else
    error('Stress state unknown')
end

%% Variables
P_source = zeros(nx,ny);   
P_source_1 = zeros(nx,ny);
P_source_2 = zeros(nx,ny);  % Fluid source pressure field (impermeable source is handled in initial conditions)
P_source_3 = zeros(nx,ny);
P_source2 = zeros(nx,ny); 
P_source3 = zeros(nx,ny); 
P_source4 = zeros(nx,ny); 
 
parameter_damasi;
P_old = zeros(nx,ny);               % Solution field
kc0 = ones(nx,ny);                  % Permeability matrix (m^2)
kc = ones(nx,ny);                  % Permeability matrix (m^2)
sig_star = zeros(nx,ny);            % Stress norm factor matrix
alpha = zeros(nx,ny);               % Permeability recovery factor matrix
Temp = y2d.*Geotherm_grad + 0;      % Temperature, assuming a constant gradient and 0 degrees on top
Fi = porosity_rock.*ones(nx,ny);    % Porosity matrix
Beta = Fi.*(beta_phi + beta_fluid); % beta coefficient used for the diffusion equation.
Visc = ones(nx,ny)*10^(-3)-0.9*(y2d./max(max(y2d))).*ones(nx,ny)*10^(-3); %visc in Pa-s (20 dec C=1e-3 Pa-s)
theal = zeros(nx,ny);               % Healing time. Basically stores the number of timesteps since healing became active.
trig_tog = zeros(nx,ny);            % Stores if the point has already been triggered.
plot_trig = zeros(nx,ny);           % Stores the points that were triggered in the current timestep
parameter_damasi;

%% Mechanic (Plot initial failure condition)
fric_fail = 0.6; 
eulera = 15;  %azimut of S3
eulerb = 90; %if S1 is vertical
eulerc = 0;
for i=1:nx
   for j=1:ny
      stress_tensor = [Sigma_1(i,j)-(P_old(i,j) + P_fluid_hydro(i,j)) 0 0; 0 Sigma_2(i,j)-(P_old(i,j) + P_fluid_hydro(i,j)) 0; 0 0 Sigma_3(i,j)-(P_old(i,j) + P_fluid_hydro(i,j))];
      [normS(i,j) shearS(i,j)]=shearstress_on_fracture(stress_tensor,eulera,eulerb,eulerc,azimuth(i,j),Dip_angle(i,j)); 
   end
end

%% Simulation Begin
time_count=1;           
time = 0;               

k_mean = NaN(floor(tot_time/dt)+1,1);
p_mean = NaN(floor(tot_time/dt)+1,1);
Q = NaN(floor(tot_time/dt)+1,1);
time_trig = zeros(floor(tot_time/dt)+1,1);
P_source0 = P_source;
P_old = P_source;

while time<tot_time
  count_trig=0;
  time =time+dt;
  n_iter= 0;

  if time_count==int16(3/dt*Day) || time_count==int16(12/dt*Day) || time_count==int16(24/dt*Day)
    theal = 0.*theal;
    theal2 = 0.*theal2;
  end
  
  for i=1:nx
    for j=1:ny
        stress_tensor = [Sigma_1(i,j)-(P_old(i,j) + P_fluid_hydro(i,j)) 0 0; 0 Sigma_2(i,j)-(P_old(i,j) + P_fluid_hydro(i,j)) 0; 0 0 Sigma_3(i,j)-(P_old(i,j) + P_fluid_hydro(i,j))];
        [sigma_n(i,j) tau_a(i,j)]=shearstress_on_fracture(stress_tensor,eulera,eulerb,eulerc,azimuth(i,j),Dip_angle(i,j)); 
        
        if (material(i,j)==1 || material(i,j)==2 )
            if ((tau_a(i,j)./sigma_n(i,j)>fric_fail||tau_a(i,j)./sigma_n(i,j)<-fric_fail)&&trig_tog(i,j)==0)
                trig_tog(i,j)=1;
                plot_trig(i,j)= plot_trig(i,j)+1;
                count_trig=count_trig+1;
                theal(i,j)=0;
            end
            theal(i,j)=theal(i,j)+dt;
            theal2(i,j)=theal(i,j)/24/3600;

        if (trig_tog(i,j)==1)

        if time_count<=int16(2/dt*Day)
            zeta = 0.6;
            P_source(i,j)=P_source0(i,j).*exp(-zeta.*theal2(i,j));
            kc(i,j)=kc0(i,j)+ (kc0(i,j)*100)*exp(-alpha(i,j).*theal2(i,j));
            P_source2(i,j) = P_source(i,j)+P_source_1(i,j);

        elseif time_count>int16(2/dt*Day) && time_count<=int16(11/dt*Day)
            P_source(i,j)=P_source2(i,j).*exp(-zeta.*theal2(i,j));
            kc(i,j)=kc0(i,j)+ (kc0(i,j)*100)*exp(-alpha(i,j).*theal2(i,j));
            P_source3(i,j) = P_source(i,j)+P_source_2(i,j);

        elseif time_count>int16(11/dt*Day) && time_count<=int16(23/dt*Day)
            P_source(i,j)=P_source3(i,j).*exp(-zeta.*theal2(i,j));
            kc(i,j)=kc0(i,j)+ (kc0(i,j)*100)*exp(-alpha(i,j).*theal2(i,j));
            P_source4(i,j) = P_source(i,j)+P_source_3(i,j);
        else
            P_source(i,j)=P_source4(i,j).*exp(-zeta.*theal2(i,j));
            kc(i,j)=kc0(i,j)+ (kc0(i,j)*100)*exp(-alpha(i,j).*theal2(i,j));
        end

        end
       end
        
    end
  end
    K                   =   kc.*exp(-sigma_n./1e6./sig_star);         
    k_mu                =   (1./Visc).*K;
    [P_new]             =   p2d(P_old, Beta, k_mu, P_source, nx, ny, dx, dy, dt, Bound_type_hor, Bound_val_hor, Bound_type_ver, Bound_val_ver);
    error_iterations    =   max(max(abs(P_old-P_new)))./max(max(abs(P_new)));
    dP_dt               =   (P_new-P_old)/(dt/Day);
    P_old               =   P_new;
    n_iter              =   n_iter + 1;
    Q                   =   (P_source./Beta);
    disp(['The error of iterations # ', num2str(n_iter),' = ',num2str(count_trig)]);

  time_count=time_count+1;
end

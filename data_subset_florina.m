%% Define cross-section
    % Coordinates of cross-section
x0 = 7.5;
%x0 = 10.5;
y0 = 3;
yend = 15;
strike = 355-180-90;
%strike = 233 - 180 
xend = x0 + (yend - y0) / tan(-strike*pi/180.0);

xd = [x0 xend];
yd = [y0 yend];

    % Ordinate a the origin
b = y0-(-tand(strike).*x0);

    % Degrees to km
lat0 = 40.68;
long0 = 21.27;
northing = (Lat-lat0).*111.11;
easting = (Long-long0).*111.11.*cosd(Lat);
northing_ms = (main_shock.Latitude-lat0)*111.11;
easting_ms = (main_shock.Longitude-long0).*111.11.*cosd(main_shock.Latitude);

%% Rotation of data
    % Defining rotation matrix
rot_mat = zeros(2,2);
rot_mat(1,1) = cosd(355-180);
rot_mat(1,2) = -sind(355-180);
rot_mat(2,1) = sind(355-180);
rot_mat(2,2) = cosd(355-180);

    % Rotate data
p1_tr    = [easting northing-b];              % translation to get origin at (0,0)
p1_rot   = rot_mat*p1_tr'; p1_rot = p1_rot';  % rotation
p1_final = [p1_rot(:,1) p1_rot(:,2)+b(1)];    % translation back to original 

ms_tr    = [easting_ms northing_ms-b];              % translation to get origin at (0,0)
ms_rot   = rot_mat*ms_tr'; ms_rot = ms_rot';  % rotation
ms_final = [ms_rot(:,1) ms_rot(:,2)+b(1)];    % translation back to original 


    % Origin of cross-sections
cs1_tr = [xd(1,:)' (yd(1,:)-b)'];
cs1_rot = rot_mat*cs1_tr'; cs1_rot = cs1_rot';
cs1_final = [cs1_rot(:,1) cs1_rot(:,2)+b'];

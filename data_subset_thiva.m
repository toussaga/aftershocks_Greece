%% Create data subsets
    % Coordinates of cross-sections
num_section = 4;

Lx = [1 15];
ymin = Lx.*tand(-30)+6.5;
ymax = Lx.*tand(-30)+13.5;
Ly = [ymin; ymax];

Lx2 = [12 20];
ymin2 = Lx2.*tand(-30)+9.5;
ymax2 = Lx2.*tand(-30)+16.5;
Ly2 = [ymin2; ymax2];

x1d = linspace(1.8,7,3);
ly = [1 10];
x22d = (ly(2)-ly(1)+x1d.*tand(60))./tand(60);
lx = [x1d(1) x22d(1); x1d(2) x22d(2); x1d(3) x22d(3); NaN NaN];

for i=1:num_section
    [xi(i),yi(i)] = linexline(Lx,Ly(1,:),lx(i,:),ly,0);
    [xj(i),yj(i)] = linexline(Lx,Ly(2,:),lx(i,:),ly,0);
end 
[xi(4),yi(4)] = linexline(Lx2,Ly2(1,:),[14.5 (ly(2)-ly(1)+14.5*tand(60))./tand(60)],ly,0);
[xj(4),yj(4)] = linexline(Lx2,Ly2(2,:),[14.5 (ly(2)-ly(1)+14.5*tand(60))./tand(60)],ly,0);
xd = [xi; xj];
yd = [yi; yj];

    % Spacing between sections
for i =2:num_section
    sp(i) = sqrt((xd(1,i)-xd(1,i-1)).^2+(yd(1,i)-yd(1,i-1)).^2)./2+0.165;
end
sp(1) = [];

    % Ordinate a the origin
for i=1:num_section
    b(i) = ly(1)-(tand(60).*xd(1,i));
end

    % Degrees to km
lat0 = 38.27;
long0 = 23.25;
northing = (Lat-lat0)*111.11;
easting = (Long-long0).*111.11.*cosd(Lat);

    % Data subset
af1 = find(inpolygon(easting,northing,[0 xd(1,1)+sp(1) xd(2,1)+sp(1) 0 0], [yd(1,1) yd(1,1) yd(2,1) yd(2,1) yd(1,1)]));
af2 = find(inpolygon(easting,northing,[xd(1,2)-sp(2) xd(1,2)+sp(2) xd(2,2)+sp(2) xd(2,2)-sp(2) xd(1,2)-sp(2)], [yd(1,2) yd(1,2) yd(2,2) yd(2,2) yd(1,2)]));
af3 = find(inpolygon(easting,northing,[xd(1,3)-sp(2) xd(1,3)+sp(3) xd(2,3)+sp(3) xd(2,3)-sp(3) xd(1,3)-sp(2)], [yd(1,3) yd(1,3) yd(2,3) yd(2,3) yd(1,3)]));
af4 = find(inpolygon(easting,northing,[xd(1,4)-sp(3) 28 28 xd(2,4)-sp(3) xd(1,4)-sp(3)], [yd(1,4) yd(1,4) yd(2,4) yd(2,4) yd(1,4)]));

af2 = sort([af2; 246; 942; 1011]);
    
% relative days
absoluteDays1 = days(yr(af1)-yr(af1(1)));
ev1 = 1:1:length(af1);
absoluteDays2 = days(yr(af2)-yr(af2(1)));
ev2 = 1:1:length(af2);
absoluteDays3 = days(yr(af3)-yr(af3(1)));
ev3 = 1:1:length(af3);
absoluteDays4 = days(yr(af4)-yr(af4(1)));
ev4 = 1:1:length(af4);


%% Rotation of data subsets
    % Defining rotation matrix
rot_mat = zeros(2,2);
rot_mat(1,1) = cosd(30);
rot_mat(1,2) = -sind(30);
rot_mat(2,1) = sind(30);
rot_mat(2,2) = cosd(30);

    % Rotate data
p1_tr = [easting(af1) northing(af1)-b(1)]; % translation to get origin at (0,0)
p1_rot = rot_mat*p1_tr'; p1_rot = p1_rot'; % rotation
p1_final = [p1_rot(:,1) p1_rot(:,2)+b(1)]; % translation back to original 

p2_tr = [easting(af2) northing(af2)-b(2)]; % translation to get origin at (0,0)
p2_rot = rot_mat*p2_tr'; p2_rot = p2_rot'; % rotation
p2_final = [p2_rot(:,1) p2_rot(:,2)+b(2)]; % translation back to original 

p3_tr = [easting(af3) northing(af3)-b(3)]; % translation to get origin at (0,0)
p3_rot = rot_mat*p3_tr'; p3_rot = p3_rot'; % rotation
p3_final = [p3_rot(:,1) p3_rot(:,2)+b(3)]; % translation back to original

p4_tr = [easting(af4) northing(af4)-b(4)]; % translation to get origin at (0,0)
p4_rot = rot_mat*p4_tr'; p4_rot = p4_rot'; % rotation
p4_final = [p4_rot(:,1) p4_rot(:,2)+b(4)]; % translation back to original

    % Origin of cross-sections
cs1_tr = [xd(1,:)' (yd(1,:)-b)'];
cs1_rot = rot_mat*cs1_tr'; cs1_rot = cs1_rot';
cs1_final = [cs1_rot(:,1) cs1_rot(:,2)+b'];

cs2_tr = [xd(2,:)' (yd(2,:)-b)'];
cs2_rot = rot_mat*cs2_tr'; cs2_rot = cs2_rot';
cs2_final = [cs2_rot(:,1) cs2_rot(:,2)+b'];

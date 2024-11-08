%% Create data subsets For Damasi sequence
% Coordinates of cross-sections
num_section = 6;

Lx = [1 45];
ymin = Lx.*tand(-40)+26;
ymax = Lx.*tand(-40)+54;
Ly = [ymin; ymax];

x1d = linspace(-15,26,num_section);
ly = [0 40];
x2d = (ly(2)-ly(1)+x1d.*tand(50))./tand(50);
%lx = [x1d(1) x2d(1); x1d(2) x2d(2); x1d(3) x2d(3); x1d(4) x2d(4); NaN NaN];

for i=1:num_section
    lx(i,:) = [x1d(i) x2d(i)];
    [xi(i),yi(i)] = linexline(Lx,Ly(1,:),lx(i,:),ly,0);
    hold on
    [xj(i),yj(i)] = linexline(Lx,Ly(2,:),lx(i,:),ly,0);
    hold on
end 
xd = [xi; xj];
yd = [yi; yj];

    % Spacing between sections
for i =2:num_section
    sp(i) = sqrt((xd(1,i)-xd(1,i-1)).^2+(yd(1,i)-yd(1,i-1)).^2)./2;
end
sp(1) = [];

    % Ordinate a the origin
for i=1:num_section
    b(i) = ly(1)-(tand(50).*xd(1,i));
end

    % Degrees to km
lat0 = 39.55;
long0 = 21.85;
northing = (Lat-lat0)*111.11;
easting =  (Long-long0).*111.11.*cosd(Lat);

    % Data subset
af1 = find(inpolygon(easting,northing,[0 xd(1,1)+sp(1) xd(2,1)+sp(1) 0 0], [yd(1,1) yd(1,1) yd(2,1) yd(2,1) yd(1,1)]));
af2 = find(inpolygon(easting,northing,[xd(1,2)-sp(2) xd(1,2)+sp(2) xd(2,2)+sp(2) xd(2,2)-sp(2) xd(1,2)-sp(2)], [yd(1,2) yd(1,2) yd(2,2) yd(2,2) yd(1,2)]));
af3 = find(inpolygon(easting,northing,[xd(1,3)-sp(3) xd(1,3)+sp(3) xd(2,3)+sp(3) xd(2,3)-sp(3) xd(1,3)-sp(3)], [yd(1,3) yd(1,3) yd(2,3) yd(2,3) yd(1,3)]));
af4 = find(inpolygon(easting,northing,[xd(1,4)-sp(4) xd(1,4)+sp(4) xd(2,4)+sp(4) xd(2,4)-sp(4) xd(1,4)-sp(4)], [yd(1,4) yd(1,4) yd(2,4) yd(2,4) yd(1,4)]));
af5 = find(inpolygon(easting,northing,[xd(1,5)-sp(5) xd(1,5)+sp(5) xd(2,5)+sp(5) xd(2,5)-sp(5) xd(1,5)-sp(5)], [yd(1,5) yd(1,5) yd(2,5) yd(2,5) yd(1,5)]));
af6 = find(inpolygon(easting,northing,[xd(1,6)-sp(5) max(easting) max(easting) xd(2,6)-sp(5) xd(1,6)-sp(5)], [yd(1,6) yd(1,6) yd(2,6) yd(2,6) yd(1,6)]));

af2 = af2(2:end);
    
% relative days
absoluteDays1 = days(yr(af1)-yr(af1(1)));
ev1 = 1:1:length(af1);
absoluteDays2 = days(yr(af2)-yr(af2(1)));
ev2 = 1:1:length(af2);
absoluteDays3 = days(yr(af3)-yr(af3(1)));
ev3 = 1:1:length(af3);
absoluteDays4 = days(yr(af4)-yr(af4(1)));
ev4 = 1:1:length(af4);
absoluteDays5 = days(yr(af5)-yr(af5(1)));
ev5 = 1:1:length(af5);
absoluteDays6 = days(yr(af6)-yr(af6(1)));
ev6 = 1:1:length(af6);


%% Rotation of data subsets
    % Defining rotation matrix
rot_mat = zeros(2,2);
rot_mat(1,1) = cosd(40);
rot_mat(1,2) = -sind(40);
rot_mat(2,1) = sind(40);
rot_mat(2,2) = cosd(40);

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

p5_tr = [easting(af5) northing(af5)-b(5)]; % translation to get origin at (0,0)
p5_rot = rot_mat*p5_tr'; p5_rot = p5_rot'; % rotation
p5_final = [p5_rot(:,1) p5_rot(:,2)+b(5)]; % translation back to original

p6_tr = [easting(af6) northing(af6)-b(6)]; % translation to get origin at (0,0)
p6_rot = rot_mat*p6_tr'; p6_rot = p6_rot'; % rotation
p6_final = [p6_rot(:,1) p6_rot(:,2)+b(6)]; % translation back to original

    % Origin of cross-sections
cs1_tr = [xd(1,:)' (yd(1,:)-b)'];
cs1_rot = rot_mat*cs1_tr'; cs1_rot = cs1_rot';
cs1_final = [cs1_rot(:,1) cs1_rot(:,2)+b'];

cs2_tr = [xd(2,:)' (yd(2,:)-b)'];
cs2_rot = rot_mat*cs2_tr'; cs2_rot = cs2_rot';
cs2_final = [cs2_rot(:,1) cs2_rot(:,2)+b'];
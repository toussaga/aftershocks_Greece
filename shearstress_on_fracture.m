function [n s]=shearstress_on_fracture(T,a,b,c,az,dip)
% SHEARSTRESS_ON_FRACTURE compute the shear and normal stress on a fracture
% [n s]=shearstress_on_fracture(T,a,b,c,az,dip)
% INPUT:
% T : stress tensor
% a,b,c : Angle d'Euler [°] for the stress Tensor orientation
%   pex: S1 vertical, S3 vers 80°: a=80; b=90; c=0;
% az : azimuth of the fracture
% dip : dip of the fracture
% OUTPUT:
% n : normal stress on the fracture
% s : shear stress on the fracture

% calculate conversion matrix
Rlmn=[  cosd(a)*cosd(b)                             sind(a)*cosd(b)                             -sind(b)
    cosd(a)*sind(b)*sind(c)-sind(a)*cosd(c)     sind(a)*sind(b)*sind(c)+cosd(a)*cosd(c)     cosd(b)*sind(c)
    cosd(a)*sind(b)*cosd(c)+sind(a)*sind(c)     sind(a)*sind(b)*cosd(c)-cosd(a)*sind(c)     cosd(b)*cosd(c)];

% compute stresses in the geographical referentiel  [x to north, y to east, z down]
Txyz=Rlmn'*T*Rlmn;

% calculate conversion matrice for fracture local referentiel 
a=az; b=-dip; c=0;
Ruvw=[  cosd(a)*cosd(b)                             sind(a)*cosd(b)                             -sind(b)
    cosd(a)*sind(b)*sind(c)-sind(a)*cosd(c)     sind(a)*sind(b)*sind(c)+cosd(a)*cosd(c)     cosd(b)*sind(c)
    cosd(a)*sind(b)*cosd(c)+sind(a)*sind(c)     sind(a)*sind(b)*cosd(c)-cosd(a)*sind(c)     cosd(b)*cosd(c)];

% stresses in local fracture referentiel  (total stress)
Tuvw=Ruvw*Txyz*Ruvw';
s=sqrt(Tuvw(1,3)^2+Tuvw(2,3)^2);
n=Tuvw(3,3);
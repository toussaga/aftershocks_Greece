%% Material parameter and input Matrix
P_source(P_matrix==1)= P_input(1);
P_source(P_matrix==2)= P_input(2);
P_source_1(P_matrix==3)= P_input(3);    
P_source_2(P_matrix==4)= P_input(4);
P_source_3(P_matrix==5)= P_input(5);
P_source_4(P_matrix==6)= P_input(6);
P_source_5(P_matrix==7)= P_input(7);
P_source_6(P_matrix==8)= P_input(8);

kc0(material==1)=Material_mat(1,4);sig_star(material==1)=Material_mat(1,5);
kc0(material==2)=Material_mat(2,4);sig_star(material==2)=Material_mat(2,5);
kc0(material==3)=Material_mat(3,4);sig_star(material==3)=Material_mat(3,5);

alpha(material==1)=Material_mat(1,6);
alpha(material==2)=Material_mat(2,6);
alpha(material==3)=Material_mat(3,6);

Dip_angle = NaN(nx,ny);
azimuth = NaN(nx,ny);

for i=1:1:nx 
    for j=1:1:ny
        k = (i-1)*ny+j;
        azimuth(i,j)   = planes(k,1);
        Dip_angle(i,j) = planes(k,2);
    end
end

azimuth(material==2)=Faults(1).azimuth;
azimuth(material==3)=Faults(2).azimuth;
Dip_angle(material==2)=Faults(1).dip;
Dip_angle(material==3)=Faults(2).dip;

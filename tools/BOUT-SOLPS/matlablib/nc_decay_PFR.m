%%%%%%%used for decay profiles from nc file in PFR

clear 
clc

addpath('~/bout/tools/matlablib');
%exist('b2fstati','file')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read in interpolation data
%fid = netcdf.open('d3d_158772.04205_174_v1_260_0.85psi_bout.grd.nc', 'nc_write');
%fid = netcdf.open('d3d_158772.04205_174_v1_516_bout.grd.nc', 'nc_write');
%fid = netcdf.open('d3d_158772.04205_174_v1_516x96y_bout.grd.nc', 'nc_write');
fid = netcdf.open('cmod-1100212023_260x64y_v1_bout.grd.nc', 'nc_write');
info_file(strcat('cmod-1100212023_260x64y_v1_bout.grd.nc'))
inp_var_list = {'nx', 'ny', 'ixseps1', 'Rxy', 'Zxy','Niexp', 'Neexp', 'Tiexp', 'Teexp', 'pressure_s'};
%inp_var_list = {'nx', 'ny', 'ixseps1', 'Rxy', 'Zxy','pressure'};
for i=1 : length(inp_var_list)

    var_name = char(inp_var_list(i));
    vid = netcdf.inqVarID( fid, var_name);
    command = [var_name '=netcdf.getVar(fid, vid);'];
    eval(command);
    command = [var_name '=double(' var_name ');'];
    eval(command);

end
clear command inp_var_list var_name vid i fid

Rxy=Rxy';
Zxy=Zxy';
Niexp=Niexp';
Neexp=Neexp';
Tiexp=Tiexp';
Teexp=Teexp';
pressure_s=pressure_s';
% pressure_s=pressure';

S1 = Rxy(nx,:);
[ypeak]=find(S1==max(max(S1)));
clear S1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fit divertor region, ndiver(divertor poloidal grid number), num count
%number 64:4, 96:18, 128:32
ndiver=4;
num=5;
[X,Y] = meshgrid(1:nx,[1:ndiver,ny-ndiver+1:ny]);
X1=1:ixseps1;
X1=X1';

%%%%%%%%%%%%%%%%% ni
Niexp([1:ixseps1-1],[1:ndiver,ny-ndiver+1:ny])=Niexp(nx,39)/2;
node1=ixseps1;
node2=fix(ixseps1*(num-1)/num);
while (node2>=1)
    Xi=X(:,[1:node2,node1:nx]);
    Yi=Y(:,[1:node2,node1:nx]);
    Nitmp=Niexp';
    Nitmp = interp2(Xi,Yi,Nitmp([1:ndiver,ny-ndiver+1:ny],[1:node2,node1:nx]),X,Y); %X,Y,Niexp should be same type.
    Niexp(:,[1:ndiver,ny-ndiver+1:ny]) = Nitmp';
    node1=fix(1/2*node1+1/2*node2);
    node2=fix(2*node2-node1);
end
    Xi=X(:,[1:1,node1:nx]);
    Yi=Y(:,[1:1,node1:nx]);
    Nitmp=Niexp';
    Nitmp = interp2(Xi,Yi,Nitmp([1:ndiver,ny-ndiver+1:ny],[1:1,node1:nx]),X,Y); %X,Y,Niexp should be same type.
    Niexp(:,[1:ndiver,ny-ndiver+1:ny]) = Nitmp';

for i=[1:ndiver,ny-ndiver+1:ny]
   [p,~,mu] = polyfit(X1,Niexp(1:ixseps1,i),6);
   nitmp = polyval(p,X1,[],mu);
   Niexp(1:ixseps1,i) = nitmp/1.2;
end

%%%%%%%%%%%%%%%%% ne
Neexp([1:ixseps1-1],[1:ndiver,ny-ndiver+1:ny])=Neexp(nx,39)/2;
node1=ixseps1;
node2=fix(ixseps1*(num-1)/num);
while (node2>=1)
    Xi=X(:,[1:node2,node1:nx]);
    Yi=Y(:,[1:node2,node1:nx]);
    Netmp=Neexp';
    Netmp = interp2(Xi,Yi,Netmp([1:ndiver,ny-ndiver+1:ny],[1:node2,node1:nx]),X,Y); %X,Y,Niexp should be same type.
    Neexp(:,[1:ndiver,ny-ndiver+1:ny]) = Netmp';
    node1=fix(1/2*node1+1/2*node2);
    node2=fix(2*node2-node1);
end
    Xi=X(:,[1:1,node1:nx]);
    Yi=Y(:,[1:1,node1:nx]);
    Netmp=Neexp';
    Netmp = interp2(Xi,Yi,Netmp([1:ndiver,ny-ndiver+1:ny],[1:1,node1:nx]),X,Y); %X,Y,Niexp should be same type.
    Neexp(:,[1:ndiver,ny-ndiver+1:ny]) = Netmp';
    
for i=[1:ndiver,ny-ndiver+1:ny]
   [p,~,mu] = polyfit(X1,Neexp(1:ixseps1,i),6);
   netmp = polyval(p,X1,[],mu);
   Neexp(1:ixseps1,i) = netmp/1.2;
end

%%%%%%%%%%%%%%%%% Ti
 Tiexp([1:ixseps1-1],[1:ndiver,ny-ndiver+1:ny])=Tiexp(nx,39)/2;
node1=ixseps1;
node2=fix(ixseps1*(num-1)/num);
while (node2>=1)
    Xi=X(:,[1:node2,node1:nx]);
    Yi=Y(:,[1:node2,node1:nx]);
    Titmp=Tiexp';
    Titmp = interp2(Xi,Yi,Titmp([1:ndiver,ny-ndiver+1:ny],[1:node2,node1:nx]),X,Y); %X,Y,Niexp should be same type.
    Tiexp(:,[1:ndiver,ny-ndiver+1:ny]) = Titmp';
    node1=fix(1/2*node1+1/2*node2);
    node2=fix(2*node2-node1);
end
    Xi=X(:,[1:1,node1:nx]);
    Yi=Y(:,[1:1,node1:nx]);
    Titmp=Tiexp';
    Titmp = interp2(Xi,Yi,Titmp([1:ndiver,ny-ndiver+1:ny],[1:1,node1:nx]),X,Y); %X,Y,Niexp should be same type.
    Tiexp(:,[1:ndiver,ny-ndiver+1:ny]) = Titmp';
    
for i=[1:ndiver,ny-ndiver+1:ny]
   [p,~,mu] = polyfit(X1,Tiexp(1:ixseps1,i),6);
   titmp = polyval(p,X1,[],mu);
   Tiexp(1:ixseps1,i) = titmp/1.2;
end

%%%%%%%%%%%%%%%%% Te
 Teexp([1:ixseps1-1],[1:ndiver,ny-ndiver+1:ny])=Teexp(nx,39)/2;
node1=ixseps1;
node2=fix(ixseps1*(num-1)/num);
while (node2>=1)
    Xi=X(:,[1:node2,node1:nx]);
    Yi=Y(:,[1:node2,node1:nx]);
    Tetmp=Teexp';
    Tetmp = interp2(Xi,Yi,Tetmp([1:ndiver,ny-ndiver+1:ny],[1:node2,node1:nx]),X,Y); %X,Y,Niexp should be same type.
    Teexp(:,[1:ndiver,ny-ndiver+1:ny]) = Tetmp';
    node1=fix(1/2*node1+1/2*node2);
    node2=fix(2*node2-node1);
end
    Xi=X(:,[1:1,node1:nx]);
    Yi=Y(:,[1:1,node1:nx]);
    Tetmp=Teexp';
    Tetmp = interp2(Xi,Yi,Tetmp([1:ndiver,ny-ndiver+1:ny],[1:1,node1:nx]),X,Y); %X,Y,Niexp should be same type.
    Teexp(:,[1:ndiver,ny-ndiver+1:ny]) = Tetmp';
    
for i=[1:ndiver,ny-ndiver+1:ny]
   [p,~,mu] = polyfit(X1,Teexp(1:ixseps1,i),6);
   tetmp = polyval(p,X1,[],mu);
   Teexp(1:ixseps1,i) = tetmp/1.2;
end

%%%%%%%%%%%%%%%%% pressure
pressure_s([1:ixseps1-1],[1:ndiver,ny-ndiver+1:ny])=pressure_s(nx,39)/2;
node1=ixseps1;
node2=fix(ixseps1*(num-1)/num);
while (node2>=1)
    Xi=X(:,[1:node2,node1:nx]);
    Yi=Y(:,[1:node2,node1:nx]);
    prtmp=pressure_s';
    prtmp = interp2(Xi,Yi,prtmp([1:ndiver,ny-ndiver+1:ny],[1:node2,node1:nx]),X,Y); %X,Y,Niexp should be same type.
    pressure_s(:,[1:ndiver,ny-ndiver+1:ny]) = prtmp';
    node1=fix(1/2*node1+1/2*node2);
    node2=fix(2*node2-node1);
end
    Xi=X(:,[1:1,node1:nx]);
    Yi=Y(:,[1:1,node1:nx]);
    prtmp=pressure_s';
    prtmp = interp2(Xi,Yi,prtmp([1:ndiver,ny-ndiver+1:ny],[1:1,node1:nx]),X,Y); %X,Y,Niexp should be same type.
    pressure_s(:,[1:ndiver,ny-ndiver+1:ny]) = prtmp';
    
for i=[1:ndiver,ny-ndiver+1:ny]
   [p,~,mu] = polyfit(X1,pressure_s(1:ixseps1,i),6);
   prtmp = polyval(p,X1,[],mu);
   pressure_s(1:ixseps1,i) = prtmp/1.15;
end
% for i=[1:ndiver-1,ny-ndiver+1:ny]
%    pressure_s(:,i) = pressure_s(:,ndiver);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Interpolation adjust data, Pmax=280,Sep=326 for grid 516, Pmax=149, Sep=169 for grid 260,

% Pmax=ixseps1-20, Sep=ixseps1, Nsol=0;
% [X,Y] = meshgrid(1:nx,1:ny);
% Xi=X(:,[1:Pmax,Sep+Nsol:nx]);
% Yi=Y(:,[1:Pmax,Sep+Nsol:nx]);
%  
% %%%%%%%%%%%%%%%%% ni
% Nitmp=Niexp';
% Nitmp = interp2(Xi,Yi,Nitmp(:,[1:Pmax,Sep+Nsol:nx]),X,Y,'spline'); %X,Y,Niexp should be same type.
% Niexp = Nitmp';
% %%%%%%%%%%%%%%%%% ne
% Netmp=Neexp';
% Netmp = interp2(Xi,Yi,Netmp(:,[1:Pmax,Sep+Nsol:nx]),X,Y,'spline'); %X,Y,Niexp should be same type.
% Neexp = Netmp';
% %%%%%%%%%%%%%%%%% Ti
% Titmp=Tiexp';
% Titmp = interp2(Xi,Yi,Titmp(:,[1:Pmax,Sep:nx]),X,Y,'spline'); %X,Y,Niexp should be same type.
% Tiexp = Titmp';
% %%%%%%%%%%%%%%%%% Te
% Tetmp=Teexp';
% Tetmp = interp2(Xi,Yi,Tetmp(:,[1:Pmax,Sep:nx]),X,Y,'spline'); %X,Y,Niexp should be same type.
% Teexp = Tetmp';
%%%%%%%%%%%%%%%%% pressure
% ptmp=pressure_s';
% ptmp = interp2(Xi,Yi,ptmp(:,[1:Pmax,Sep:nx]),X,Y,'spline'); %keep poloidal constant after interpolate
% pressure_s = ptmp';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%use griddata function do 2D smooth
Pmax=ixseps1-1, Sep=ixseps1+10, Nsol=0;
%Pmax=ixseps1-1, Sep=ixseps1+20, Nsol=0;
[X,Y] = meshgrid(1:nx,1:ny);
Xi=X(:,[1:Pmax,Sep+Nsol:nx]);
Yi=Y(:,[1:Pmax,Sep+Nsol:nx]);
ptmp=pressure_s';
ptmp =griddata(Xi,Yi,ptmp(:,[1:Pmax,Sep:nx]),X,Y,'v4'); %with poloidal decay after interpolate,must use interp2 together adjust the discontinuty point
%ptmp = interp2(Xi,Yi,ptmp(:,[1:Pmax,Sep:nx]),X,Y,'spline'); 
pressure_s = ptmp';

Pmax=ixseps1-1, Sep=ixseps1+20, Nsol=0;
%Pmax=ixseps1-1, Sep=ixseps1+20, Nsol=0;
[X,Y] = meshgrid(1:nx,1:ny);
Xi=X(:,[1:Pmax,Sep+Nsol:nx]);
Yi=Y(:,[1:Pmax,Sep+Nsol:nx]);
%%%%%%%%%%%%%%%%% Ti
Titmp=Tiexp';
Titmp = griddata(Xi,Yi,Titmp(:,[1:Pmax,Sep:nx]),X,Y,'v4'); %with poloidal decay after interpolate,must use interp2 together adjust the discontinuty point
Tiexp = Titmp';
%%%%%%%%%%%%%%%%% Te
Tetmp=Teexp';
Tetmp = griddata(Xi,Yi,Tetmp(:,[1:Pmax,Sep:nx]),X,Y,'v4'); %with poloidal decay after interpolate,must use interp2 together adjust the discontinuty point
Teexp = Tetmp';
%%%%%%%%%%%%%%%%% Ni
Nitmp=Niexp';
Nitmp = griddata(Xi,Yi,Nitmp(:,[1:Pmax,Sep:nx]),X,Y,'v4'); %with poloidal decay after interpolate,must use interp2 together adjust the discontinuty point
Niexp = Nitmp';
%%%%%%%%%%%%%%%%% Ne
Netmp=Neexp';
Netmp = griddata(Xi,Yi,Netmp(:,[1:Pmax,Sep:nx]),X,Y,'v4'); %with poloidal decay after interpolate,must use interp2 together adjust the discontinuty point
Neexp = Netmp';

%%%%%%%%%%%%%%%%%%%%%%%%use interp2 function do 2D smooth in PFR only radially no poloidal effect
%Pmax=ixseps1-20, Sep=ixseps1-2, Nsol=0;
Pmax=ixseps1-20, Sep=ixseps1, Nsol=0;
[X,Y] = meshgrid(1:nx,1:ny);
Xi=X(:,[1:Pmax,Sep+Nsol:nx]);
Yi=Y(:,[1:Pmax,Sep+Nsol:nx]);
ptmp=pressure_s';
ptmp = interp2(Xi,Yi,ptmp(:,[1:Pmax,Sep:nx]),X,Y,'spline'); %X,Y,Niexp should be same type.
pressure_s = ptmp'; 
%%%%%%%%%%%%%%%% Ti
Titmp=Tiexp';
Titmp = interp2(Xi,Yi,Titmp(:,[1:Pmax,Sep:nx]),X,Y,'spline'); %X,Y,Niexp should be same type.
Tiexp = Titmp';
%%%%%%%%%%%%%%%%% Te
Tetmp=Teexp';
Tetmp = interp2(Xi,Yi,Tetmp(:,[1:Pmax,Sep:nx]),X,Y,'spline'); %X,Y,Niexp should be same type.
Teexp = Tetmp';
%%%%%%%%%%%%%%%%% Ni
Nitmp=Niexp';
Nitmp = interp2(Xi,Yi,Nitmp(:,[1:Pmax,Sep:nx]),X,Y,'spline'); %with poloidal decay after interpolate,must use interp2 together adjust the discontinuty point
Niexp = Nitmp';
%%%%%%%%%%%%%%%%% Ne
Netmp=Neexp';
Netmp = interp2(Xi,Yi,Netmp(:,[1:Pmax,Sep:nx]),X,Y,'spline'); %with poloidal decay after interpolate,must use interp2 together adjust the discontinuty point
Neexp = Netmp';


% 
% Neexp = pressure_s./(Tiexp+Teexp)/16.02;
% Niexp = Neexp;

clear Netmp Nitmp Tetmp Titmp Nimptmp ertmp Nsol Xi Yi X Y


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% figure(1);
% clf;
% pcolor(Niexp);
% colormap('default');
% colorbar;
% shading interp;
% xlabel('grid Y');
% P=strcat(' grid X  ');
% ylabel(P);
% legend('location','NorthWest', 'ni')
% xlabel('grid Y','fontsize',18,'fontname','times');
% ylabel('grid X','fontsize',18,'fontname','times');
% set(gca,'fontweight','bold','fontsize',18,'fontname','times')
% set(gca,'linewidth',2)
%xlim([0 80])
%set(gca,'xtick',[0:5:80]);
%set(gca,'xminortick','on')
%set(gca,'xticklabel',{0:5:80});
% 
% figure(2);
% clf;
% surface(Niexp);
% colorbar('hsv')
% shading interp;
% grid on
% legend('location','NorthWest', 'ni')
% xlabel('Y','fontsize',18,'fontname','times');
% ylabel('X','fontsize',18,'fontname','times');
% set(gca,'fontweight','bold','fontsize',18,'fontname','times')
% set(gca,'linewidth',2)
% view(-260,36);
% 
% figure(3);
% clf;
% pcolor(Tiexp);
% colormap('default');
% colorbar;
% shading interp;
% xlabel('grid Y');
% P=strcat(' grid X  ');
% ylabel(P);
% legend('location','NorthWest', 'Ti')
% xlabel('grid Y','fontsize',18,'fontname','times');
% ylabel('grid X','fontsize',18,'fontname','times');
% set(gca,'fontweight','bold','fontsize',18,'fontname','times')
% set(gca,'linewidth',2)
% 
% figure(4);
% clf;
% surface(Teexp);
% colorbar('hsv')
% shading interp;
% grid on
% legend('location','NorthWest', 'Ti')
% xlabel('Y','fontsize',18,'fontname','times');
% ylabel('X','fontsize',18,'fontname','times');
% set(gca,'fontweight','bold','fontsize',18,'fontname','times')
% set(gca,'linewidth',2)
% view(-260,36);
% 
figure(5);
plot(Niexp(:,55)),hold
plot(Niexp(:,1))
plot(Niexp(:,10))
plot(Niexp(:,20))
plot(Niexp(:,30))
plot(Niexp(:,40))
plot(Niexp(:,50))
legend('location','Northwest','ni,midplane','ni1','ni10','ni20','ni30', 'ni40','ni50')
xlabel('grid X','fontsize',18,'fontname','times');
ylabel('n_{\iti} (10^{20} m^{-3})','fontsize',18,'fontname','times');
set(gca,'fontweight','bold','fontsize',18,'fontname','times')
set(gca,'linewidth',2)

figure(6)
plot(Tiexp(:,55)),hold
plot(Tiexp(:,1))
plot(Tiexp(:,10))
plot(Tiexp(:,20))
plot(Tiexp(:,30))
plot(Tiexp(:,40))
plot(Tiexp(:,50))
legend('location','Northeast','Ti,midplane','Ti1','Ti10','Ti20','Ti30', 'Ti40','Ti50')
xlabel('grid X','fontsize',18,'fontname','times');
ylabel('T_{\iti} (eV)','fontsize',18,'fontname','times');
set(gca,'fontweight','bold','fontsize',18,'fontname','times')
set(gca,'linewidth',2)

figure(7)
plot(Teexp(:,55)),hold
plot(Teexp(:,1))
plot(Teexp(:,10))
plot(Teexp(:,20))
plot(Teexp(:,30))
plot(Teexp(:,40))
plot(Teexp(:,50))
legend('location','Northeast','Te,midplane','Te1','Te10','Te20','Te30', 'Te40','Te50')
xlabel('grid X','fontsize',18,'fontname','times');
ylabel('T_{\ite} (eV)','fontsize',18,'fontname','times');
set(gca,'fontweight','bold','fontsize',18,'fontname','times')
set(gca,'linewidth',2)

figure(8);
plot(Neexp(:,55)),hold
plot(Neexp(:,1))
plot(Neexp(:,10))
plot(Neexp(:,20))
plot(Neexp(:,30))
plot(Neexp(:,40))
plot(Neexp(:,50))
legend('location','Northwest','ne,midplane','ne1','ne10','ne20','ne30', 'ne40','ne50')
xlabel('grid X','fontsize',18,'fontname','times');
ylabel('n_{\ite} (10^{20} m^{-3})','fontsize',18,'fontname','times');
set(gca,'fontweight','bold','fontsize',18,'fontname','times')
set(gca,'linewidth',2)

figure(9);
clf;
surface(pressure_s);
colorbar('hsv')
shading interp;
grid on
legend('location','NorthWest', 'pressure_s')
xlabel('Y','fontsize',18,'fontname','times');
ylabel('X','fontsize',18,'fontname','times');
set(gca,'fontweight','bold','fontsize',18,'fontname','times')
set(gca,'linewidth',2)
view(-260,36);

figure(10);
plot(pressure_s(:,55)),hold
plot(pressure_s(:,1))
plot(pressure_s(:,2))
plot(pressure_s(:,3))
plot(pressure_s(:,4))
plot(pressure_s(:,63))
plot(pressure_s(:,64))
legend('location','Northwest','pr,midplane','pr1','pr2','pr3','pr4', 'pr63','pr64')
xlabel('grid X','fontsize',18,'fontname','times');
ylabel('pressure_s (Pa)','fontsize',18,'fontname','times');
set(gca,'fontweight','bold','fontsize',18,'fontname','times')
set(gca,'linewidth',2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% save data
%fid = netcdf.open('d3d_158772.04205_174_v1_516_bout.grd.nc', 'nc_write');
fid = netcdf.open('cmod-1100212023_260x64y_v1_bout.grd.nc', 'nc_write');
%fid = netcdf.open('d3d_158772.04205_174_v1_516x96y_bout.grd.nc', 'nc_write');
vid = netcdf.inqVarID( fid, 'Niexp');
netcdf.putVar(fid, vid, Niexp');

vid = netcdf.inqVarID( fid, 'Neexp');
netcdf.putVar(fid, vid, Neexp');

vid = netcdf.inqVarID( fid, 'Tiexp');
netcdf.putVar(fid, vid, Tiexp');

vid = netcdf.inqVarID( fid, 'Teexp');
netcdf.putVar(fid, vid, Teexp');

% vid = netcdf.inqVarID( fid, 'N_imp');
% netcdf.putVar(fid, vid, N_imp');

vid = netcdf.inqVarID( fid, 'pressure_s');
netcdf.putVar(fid, vid, pressure_s');

% vid = netcdf.inqVarID( fid, 'E_r');
% netcdf.putVar(fid, vid, E_r');

netcdf.close(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adjust inner divertor profile near target
%for j=1:2
%    for i=1:nx
%Niexp(i,j)=Niexp(i,3);
%Neexp(i,j)=Neexp(i,3);
%N_imp(i,j)=N_imp(i,3);
%Tiexp(i,j)=Tiexp(i,3);
%Teexp(i,j)=Teexp(i,3);
%    end
%end

 

 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %fit divertor region to 2 order continuty on 2D
% ndiver=36;
% num=10;
% [X,Y] = meshgrid(1:nx,[1:ndiver,ny-ndiver+1:ny]);
% 
% %%%%%%%%%%%%%%%%% ni
% node1=ixseps1;
% node2=fix(ixseps1*(num-1)/num);
% while (node2>=1)
%     Xi=X(:,[1:node2,node1:nx]);
%     Yi=Y(:,[1:node2,node1:nx]);
%     Nitmp=Niexp';
%     Nitmp = interp2(Xi,Yi,Nitmp([1:ndiver,ny-ndiver+1:ny],[1:node2,node1:nx]),X,Y,'spline'); %X,Y,Niexp should be same type.
%     Niexp(:,[1:ndiver,ny-ndiver+1:ny]) = Nitmp';
%     node1=fix(1/2*node1+1/2*node2);
%     node2=fix(2*node2-node1);
% end
%     Xi=X(:,[1:1,node1:nx]);
%     Yi=Y(:,[1:1,node1:nx]);
%     Nitmp=Niexp';
%     Nitmp = interp2(Xi,Yi,Nitmp([1:ndiver,ny-ndiver+1:ny],[1:1,node1:nx]),X,Y,'spline'); %X,Y,Niexp should be same type.
%     Niexp(:,[1:ndiver,ny-ndiver+1:ny]) = Nitmp';
% 
% %%%%%%%%%%%%%%%%% ne
% node1=ixseps1;
% node2=fix(ixseps1*(num-1)/num);
% while (node2>=1)
%     Xi=X(:,[1:node2,node1:nx]);
%     Yi=Y(:,[1:node2,node1:nx]);
%     Netmp=Neexp';
%     Netmp = interp2(Xi,Yi,Netmp([1:ndiver,ny-ndiver+1:ny],[1:node2,node1:nx]),X,Y,'spline'); %X,Y,Niexp should be same type.
%     Neexp(:,[1:ndiver,ny-ndiver+1:ny]) = Netmp';
%     node1=fix(1/2*node1+1/2*node2);
%     node2=fix(2*node2-node1);
% end
%     Xi=X(:,[1:1,node1:nx]);
%     Yi=Y(:,[1:1,node1:nx]);
%     Netmp=Neexp';
%     Netmp = interp2(Xi,Yi,Netmp([1:ndiver,ny-ndiver+1:ny],[1:1,node1:nx]),X,Y,'spline'); %X,Y,Niexp should be same type.
%     Neexp(:,[1:ndiver,ny-ndiver+1:ny]) = Netmp';
%     
% %%%%%%%%%%%%%%%%% n_imp
% node1=ixseps1;
% node2=fix(ixseps1*(num-1)/num);
% while (node2>=1)
%     Xi=X(:,[1:node2,node1:nx]);
%     Yi=Y(:,[1:node2,node1:nx]);
%     Nimptmp=N_imp';
%     Nimptmp = interp2(Xi,Yi,Nimptmp([1:ndiver,ny-ndiver+1:ny],[1:node2,node1:nx]),X,Y,'spline'); %X,Y,Niexp should be same type.
%     N_imp(:,[1:ndiver,ny-ndiver+1:ny]) = Nimptmp';
%     node1=fix(1/2*node1+1/2*node2);
%     node2=fix(2*node2-node1);
% end
%     Xi=X(:,[1:1,node1:nx]);
%     Yi=Y(:,[1:1,node1:nx]);
%     Nimptmp=N_imp';
%     Nimptmp = interp2(Xi,Yi,Nimptmp([1:ndiver,ny-ndiver+1:ny],[1:1,node1:nx]),X,Y,'spline'); %X,Y,Niexp should be same type.
%     N_imp(:,[1:ndiver,ny-ndiver+1:ny]) = Nimptmp';
%     
% %%%%%%%%%%%%%%%%% Ti
% node1=ixseps1;
% node2=fix(ixseps1*(num-1)/num);
% while (node2>=1)
%     Xi=X(:,[1:node2,node1:nx]);
%     Yi=Y(:,[1:node2,node1:nx]);
%     Titmp=Tiexp';
%     Titmp = interp2(Xi,Yi,Titmp([1:ndiver,ny-ndiver+1:ny],[1:node2,node1:nx]),X,Y,'spline'); %X,Y,Niexp should be same type.
%     Tiexp(:,[1:ndiver,ny-ndiver+1:ny]) = Titmp';
%     node1=fix(1/2*node1+1/2*node2);
%     node2=fix(2*node2-node1);
% end
%     Xi=X(:,[1:1,node1:nx]);
%     Yi=Y(:,[1:1,node1:nx]);
%     Titmp=Tiexp';
%     Titmp = interp2(Xi,Yi,Titmp([1:ndiver,ny-ndiver+1:ny],[1:1,node1:nx]),X,Y,'spline'); %X,Y,Niexp should be same type.
%     Tiexp(:,[1:ndiver,ny-ndiver+1:ny]) = Titmp';
% 
% %%%%%%%%%%%%%%%%% Te
% node1=ixseps1;
% node2=fix(ixseps1*(num-1)/num);
% while (node2>=1)
%     Xi=X(:,[1:node2,node1:nx]);
%     Yi=Y(:,[1:node2,node1:nx]);
%     Tetmp=Teexp';
%     Tetmp = interp2(Xi,Yi,Tetmp([1:ndiver,ny-ndiver+1:ny],[1:node2,node1:nx]),X,Y,'spline'); %X,Y,Niexp should be same type.
%     Teexp(:,[1:ndiver,ny-ndiver+1:ny]) = Tetmp';
%     node1=fix(1/2*node1+1/2*node2);
%     node2=fix(2*node2-node1);
% end
%     Xi=X(:,[1:1,node1:nx]);
%     Yi=Y(:,[1:1,node1:nx]);
%     Tetmp=Teexp';
%     Tetmp = interp2(Xi,Yi,Tetmp([1:ndiver,ny-ndiver+1:ny],[1:1,node1:nx]),X,Y,'spline'); %X,Y,Niexp should be same type.
%     Teexp(:,[1:ndiver,ny-ndiver+1:ny]) = Tetmp';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%smooth profiles on x direction
X=1:nx;
X=X';
% for i=1 : ny
%     [p,~,mu] = polyfit(X,Niexp(1:nx,i),20);
%     nitmp = polyval(p,X,[],mu);
%     Niexp(1:nx,i) = nitmp;
% end
% 
%  for i=1 : ny
%     [p,~,mu] = polyfit(X,Neexp(1:nx,i),20);
%     netmp = polyval(p,X,[],mu);
%     Neexp(1:nx,i) = netmp;
%  end
% 
%  for i=1 : ny
%    [p,~,mu] = polyfit(X,Tiexp(1:nx,i),20);
%     titmp = polyval(p,X,[],mu);
%     Tiexp(1:nx,i) = titmp;
% end
%  
% for i=1 : ny
%    [p,~,mu] = polyfit(X,Teexp(1:nx,i),20);
%     tetmp = polyval(p,X,[],mu);
%     Teexp(1:nx,i) = tetmp;
% end
% 
% for i=1 : ny
%     [p,~,mu] = polyfit(X,N_imp(1:nx,i),6);
%     nimptmp = polyval(p,X,[],mu);
%     N_imp(1:nx,i) = nimptmp;
% end
% 

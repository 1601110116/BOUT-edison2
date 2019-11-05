clear all
clc
addpath('~/bout/tools/matlablib');
dirrun='/global/cscratch1/sd/binchen/nonlinear/cmod/CMod-high-current-1237';
dir=strcat(dirrun,'/data');
info_file(strcat(dir,'/BOUT.dmp.0.nc'));

Var='vexb_x';
vexb_x=import_data_netcdf(strcat(dir),Var,737,1);
disp(Var)

clear sa
sa=strcat(Var,'.mat');
save(sa,Var,'-v7.3');

Var='vbtild_x';
vbtild_x=import_data_netcdf(strcat(dir),Var,737,1);
disp(Var)

clear sa
sa=strcat(Var,'.mat');
save(sa,Var,'-v7.3');
% 
Var='Ni';
Ni=import_data_netcdf(strcat(dir),Var,737,1);
disp(Var)

clear sa
sa=strcat(Var,'.mat');
save(sa,Var,'-v7.3');

Var='N0';
N0=import_dmp(dir, Var);
disp(Var)

clear sa
sa=strcat(Var,'.mat');
save(sa,Var,'-v7.3');

Var='Ti';
Ti=import_data_netcdf(strcat(dir),Var,737,1);
disp(Var)

clear sa
sa=strcat(Var,'.mat');
save(sa,Var,'-v7.3'); 

Var='Ti0';
Ti0=import_dmp(dir, Var);
disp(Var)

clear sa
sa=strcat(Var,'.mat');
save(sa,Var,'-v7.3');

Var='Te';
Te=import_data_netcdf(strcat(dir),Var,737,1);
disp(Var)

clear sa
sa=strcat(Var,'.mat');
save(sa,Var,'-v7.3');
% 
Var='Te0';
Te0=import_dmp(dir, Var);
disp(Var)

clear sa
sa=strcat(Var,'.mat');
save(sa,Var,'-v7.3');

Var='N_imp0';
N_imp0=import_dmp(dir, Var);
disp(Var)

clear sa
sa=strcat(Var,'.mat');
save(sa,Var,'-v7.3');

% Var='heatflux_par_i';
% heatflux_par_i=import_dmp(dir, Var);
% disp(Var)
% 
% clear sa
% sa=strcat(Var,'.mat');
% save(sa,Var,'-v7.3');
% 
% Var='heatflux_par_e';
% heatflux_par_e=import_dmp(dir, Var);
% disp(Var)
% 
% clear sa
% sa=strcat(Var,'.mat');
% save(sa,Var,'-v7.3');

Var='heatflux_par_flutter_e';
heatflux_par_flutter_e=import_data_netcdf(strcat(dir),Var,737,1);
disp(Var)

clear sa
sa=strcat(Var,'.mat');
save(sa,Var,'-v7.3');

Var='heatflux_par_flutter_i';
heatflux_par_flutter_i=import_data_netcdf(strcat(dir),Var,737,1);
disp(Var)

clear sa
sa=strcat(Var,'.mat');
save(sa,Var,'-v7.3');

%%

% loaddata
% 
% Var='vexb_x'
% S=load(strcat(Var,'.mat'));
% eval([Var,'=S.',Var,';']);
% 
%  Var='vbtild_x'
%  S=load(strcat(Var,'.mat'));
%  eval([Var,'=S.',Var,';']);
% 
% Var='Ni'
% S=load(strcat(Var,'.mat'));
% eval([Var,'=S.',Var,';']);
% % 
% Var='N0'
% S=load(strcat(Var,'.mat'));
% eval([Var,'=S.',Var,';']);
% % 
%  Var='Ti'
%  S=load(strcat(Var,'.mat'));
%  eval([Var,'=S.',Var,';']);
% % 
%  Var='Ti0'
%  S=load(strcat(Var,'.mat'));
%  eval([Var,'=S.',Var,';']);
% % 
%  Var='Te0'
%  S=load(strcat(Var,'.mat'));
%  eval([Var,'=S.',Var,';']);
% % 
%  Var='Te'
%  S=load(strcat(Var,'.mat'));
%  eval([Var,'=S.',Var,';']);
% % 
%  Var='N_imp0';
%  S=load(strcat(Var,'.mat'));
%  eval([Var,'=S.',Var,';']);
%  
Var='heatflux_par_e';
S=load(strcat(Var,'.mat'));
eval([Var,'=S.',Var,';']);

Var='heatflux_par_i';
S=load(strcat(Var,'.mat'));
eval([Var,'=S.',Var,';']);
% % 
% Var='heatflux_par_flutter_e';
% S=load(strcat(Var,'.mat'));
% eval([Var,'=S.',Var,';']);
%  
% Var='heatflux_par_flutter_i';
% S=load(strcat(Var,'.mat'));
% eval([Var,'=S.',Var,';']);

clear S Var
% 
Nbar = import_dmp(dir,'Nbar')
Tbar = import_dmp(dir,'Tbar')
Tibar = import_dmp(dir,'Tibar')
Tebar = import_dmp(dir,'Tebar')
Lbar = import_dmp(dir,'Lbar')
Va = import_dmp(dir,'Va')
% 
%load data from grid file
dir=strcat(dirrun);
fid = netcdf.open(strcat(dir,'/cmod-1100212023_260x64y_v1_bout.grd.nc'), 'nc_write');
% 
inp_var_list = {'nx', 'ny', 'ixseps1','Rxy','Zxy','psixy', 'psi_axis', 'psi_bndry'}
% 
for i=1 : length(inp_var_list)
     var_name = char(inp_var_list(i));
     vid = netcdf.inqVarID( fid, var_name);
     command = [var_name '=netcdf.getVar(fid, vid);'];
     eval(command);
     command = [var_name '=double(' var_name ');'];
     eval(command);
end
Rxy=Rxy';
Zxy=Zxy';

clear command fid vid

 psi_norm = (psixy - psi_axis)/(psi_bndry - psi_axis);
 R=psi_norm(1,:);
 nz=65
 num=4
 I=size(Ni)
 t=I(4)
 Tepara1 = 1.027354e-07
 pinch=0.022 %pinch is averaged Bp/B near separatrix in divertor entrance
 
%%

disp('sec1')
%%
%kb=1.3806505*1e-23;
       
nit=(Ni)*Nbar*1e20;
Tit=(Ti)*Tibar*1.602176565e-19;

nit0=(N0)*Nbar*1e20;
Tit0=(Ti0)*Tibar*1.602176565e-19;

net=(Ni)*Nbar*1e20;
Tet=(Te)*Tebar*1.602176565e-19;

net0=(N0+N_imp0)*Nbar*1e20;
Tet0=(Te0)*Tebar*1.602176565e-19;

qir_r=zeros(nx,ny,nz,t);
qir_p=zeros(nx,ny,nz,t);
qer_r=zeros(nx,ny,nz,t);
qer_p=zeros(nx,ny,nz,t);

qir_r_xt=zeros(ny,t);
qir_p_xt=zeros(ny,t);
qer_r_xt=zeros(ny,t);
qer_p_xt=zeros(ny,t);
% 

qi_par_yt=zeros(ny,t);
qe_par_yt=zeros(ny,t);

clear dTit dTet r dr N0 Te0 Ti0 Ti Te Ni N_imp0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate paralle total heatflux

for k=1:t
     for ix=ixseps1:nx
        for i=1:ny
        for j=1:nz-1
        
        dr=sqrt((Rxy(ix,i)-Rxy(ix-1,i))^2+(Zxy(ix,i)-Zxy(ix-1,i))^2);
    	dz=2*pi*Rxy(ix,i)/5/(nz-1);
        heatflux_par_i(ix,i,j,k)=heatflux_par_i(ix,i,j,k)*Va*Nbar*1e20*Tibar*1.602176565e-19*pinch*dr*dz;   
        heatflux_par_e(ix,i,j,k)=heatflux_par_e(ix,i,j,k)*Va*Nbar*1e20*Tebar*1.602176565e-19*pinch*dr*dz;

        end
        end
    end
end
clear dr dz
heatflux_par_i(:,:,65,:)=0.0;
heatflux_par_e(:,:,65,:)=0.0;

for k=1:t
    for i=4:ny-3
    qi_par_yt(i,k)=sum(sum(heatflux_par_i(ixseps1+1:260,i,:,k)));
    qe_par_yt(i,k)=sum(sum(heatflux_par_e(ixseps1+1:260,i,:,k)));
    end
end
qi_par_yt=qi_par_yt*5;
qe_par_yt=qe_par_yt*5;

%clear heatflux_par_i heatflux_par_e

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate radial total heatflux

for k=1:t
     for ix=ixseps1:nx
        for i=5:ny-4
        for j=1:nz-1

        dy=sqrt((Rxy(ix,i+1)-Rxy(ix,i))^2+(Zxy(ix,i+1)-Zxy(ix,i))^2);
    	dz= 2*pi*Rxy(ix,i)/5/(nz-1);
        qir_r(ix,i,j,k)=(nit(ix,i,j,k)+nit0(ix,i))*(Tit(ix,i,j,k)+Tit0(ix,i))*vexb_x(ix,i,j,k)*Va*dy*dz;
        qir_p(ix,i,j,k)=-1.5*heatflux_par_flutter_i(ix,i,j,k)*vbtild_x(ix,i,j,k)*Va*(Nbar*1e20)*Tibar*1.602176565e-19*dy*dz;   
        qer_r(ix,i,j,k)=(net(ix,i,j,k)+net0(ix,i))*(Tet(ix,i,j,k)+Tet0(ix,i))*vexb_x(ix,i,j,k)*Va*dy*dz;
        qer_p(ix,i,j,k)=-1.5*heatflux_par_flutter_e(ix,i,j,k)*vbtild_x(ix,i,j,k)*Va*(Nbar*1e20)*Tebar*1.602176565e-19*dy*dz;

	    end
        end
    end
end
%%%%%calculate different component of the radial EXB flux
        qir_r_cross=zeros(nx,ny,nz,t);
        qer_r_cross=zeros(nx,ny,nz,t);
        qir_r_cov=zeros(nx,ny,nz,t);
        qer_r_cov=zeros(nx,ny,nz,t);
        qir_r_con=zeros(nx,ny,nz,t);
        qer_r_con=zeros(nx,ny,nz,t);
        qer_r_equ=zeros(nx,ny,nz,t);
for k=1:t
     for ix=ixseps1:nx
        for i=5:ny-4
        for j=1:nz-1

        dy=sqrt((Rxy(ix,i+1)-Rxy(ix,i))^2+(Zxy(ix,i+1)-Zxy(ix,i))^2);
    	dz= 2*pi*Rxy(ix,i)/5/(nz-1);
        qir_r_cross(ix,i,j,k)=(nit(ix,i,j,k))*(Tit(ix,i,j,k))*vexb_x(ix,i,j,k)*Va*dy*dz;
        qer_r_cross(ix,i,j,k)=(net(ix,i,j,k))*(Tet(ix,i,j,k))*vexb_x(ix,i,j,k)*Va*dy*dz;
        qir_r_cov(ix,i,j,k)=(nit(ix,i,j,k))*(Tit0(ix,i))*vexb_x(ix,i,j,k)*Va*dy*dz;
        qer_r_cov(ix,i,j,k)=(net(ix,i,j,k))*(Tet0(ix,i))*vexb_x(ix,i,j,k)*Va*dy*dz;
        qir_r_con(ix,i,j,k)=(nit0(ix,i))*(Tit(ix,i,j,k))*vexb_x(ix,i,j,k)*Va*dy*dz;
        qer_r_con(ix,i,j,k)=(net0(ix,i))*(Tet(ix,i,j,k))*vexb_x(ix,i,j,k)*Va*dy*dz;
        qer_r_equ(ix,i,j,k)=(net0(ix,i))*(Tet0(ix,i))*vexb_x(ix,i,j,k)*Va*dy*dz;
	    end
        end
    end
end

%clear nit nit0 net net0 Tit Tet0 vexb_x kappa_par_i kappa_par_e vbtild_x dTitdr dTetdr
qir_r(:,:,65,:)=0.0;
qir_p(:,:,65,:)=0.0;
qer_r(:,:,65,:)=0.0;
qer_p(:,:,65,:)=0.0;
qir_r=5*qir_r;
qer_r=5*qer_r;
qir_p=5*qir_p;
qer_p=5*qer_p;

        qir_r_cross_xt=zeros(nx,t);
        qer_r_cross_xt=zeros(nx,t);
        qir_r_cov_xt=zeros(nx,t);
        qer_r_cov_xt=zeros(nx,t);
        qir_r_con_xt=zeros(nx,t);
        qer_r_con_xt=zeros(nx,t);
        qer_r_equ_xt=zeros(nx,t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculate heat flux cross flux surface outside separatrix
for k=1:t
    for ix=ixseps1:nx
    qir_r_xt(ix,k)=sum(sum(qir_r(ix,num+1:ny-num,:,k)));
    qir_p_xt(ix,k)=sum(sum(qir_p(ix,num+1:ny-num,:,k)));
    qer_r_xt(ix,k)=sum(sum(qer_r(ix,num+1:ny-num,:,k)));
    qer_p_xt(ix,k)=sum(sum(qer_p(ix,num+1:ny-num,:,k)));
    
    qir_r_cross_xt(ix,k)=sum(sum(qir_r_cross(ix,num+1:ny-num,:,k)));
    qer_r_cross_xt(ix,k)=sum(sum(qer_r_cross(ix,num+1:ny-num,:,k)));
    qir_r_cov_xt(ix,k)=sum(sum(qir_r_cov(ix,num+1:ny-num,:,k)));
    qer_r_cov_xt(ix,k)=sum(sum(qer_r_cov(ix,num+1:ny-num,:,k)));
    qir_r_con_xt(ix,k)=sum(sum(qir_r_con(ix,num+1:ny-num,:,k)));
    qer_r_con_xt(ix,k)=sum(sum(qer_r_con(ix,num+1:ny-num,:,k)));
    qer_r_equ_xt(ix,k)=sum(sum(qer_r_equ(ix,num+1:ny-num,:,k)));
    end
end

        qir_r_cross_xt=5*qir_r_cross_xt;
        qer_r_cross_xt=5*qer_r_cross_xt;
        qir_r_cov_xt=5*qir_r_cov_xt;
        qer_r_cov_xt=5*qer_r_cov_xt;
        qir_r_con_xt=5*qir_r_con_xt;
        qer_r_con_xt=5*qer_r_con_xt;
        qer_r_equ_xt=5*qer_r_equ_xt;

q_r_sep_i=qir_r_xt(ixseps1,:)+qir_p_xt(ixseps1,:);
q_r_sep_e=qer_r_xt(ixseps1,:)+qer_p_xt(ixseps1,:);
q_t_sol_i=qir_r_xt(nx-2,:)+qir_p_xt(nx-2,:)+qi_par_yt(ny-num+1,:)-qi_par_yt(num,:);
q_t_sol_e=qer_r_xt(nx-2,:)+qer_p_xt(nx-2,:)+qe_par_yt(ny-num+1,:)-qe_par_yt(num,:);

figure(1)
clf
set(gca,'position',[0.15 0.2 0.7 0.75])
plot(q_r_sep_i/1e6,'r','linewidth',2),hold
plot(q_r_sep_e/1e6,'b','linewidth',2)
plot(q_t_sol_i/1e6,'--r','linewidth',2)
plot(q_t_sol_e/1e6,'--b','linewidth',2)
xlim([0 600])
xlabel ('Alfven time \tau_{\alpha}','fontweight','bold','fontsize',14,'fontname','times')
ylabel('Heat flux q (MW)','fontweight','bold','fontsize',14,'fontname','times');
legend('location','NorthWest', 'q_i across sep.','q_e across sep.', 'q_i into divertor&wall','q_e into divertor&wall')
set(gca,'fontweight','bold','fontsize',14,'fontname','arial')
set(gca,'linewidth',1.5)

figure(2)
clf
set(gca,'position',[0.15 0.2 0.7 0.75])
plot((q_r_sep_i+q_r_sep_e)/1e6,'r','linewidth',2.5),hold
plot((q_t_sol_i+q_t_sol_e)/1e6,'b','linewidth',2.5)
%ylim([1 1e5])
xlabel ('Alfven time \tau_{\alpha}','fontweight','bold','fontsize',14,'fontname','times')
ylabel('Heat flux q (MW)','fontweight','bold','fontsize',14,'fontname','times');
legend('location','NorthWest', 'q_t across sep.','q_t into divertor&wall')
set(gca,'fontweight','bold','fontsize',14,'fontname','arial')
set(gca,'linewidth',1.5)

figure(3)
clf
set(gca,'position',[0.15 0.2 0.7 0.75])
plot((qer_r_xt(nx-2,:)+qer_p_xt(nx-2,:)+qir_r_xt(nx-2,:)+qir_p_xt(nx-2,:))/1e6,'r','linewidth',2.5),hold
plot((qi_par_yt(ny-num+1,:)-qi_par_yt(num,:)+qe_par_yt(ny-num+1,:)-qe_par_yt(num,:))/1e6,'b','linewidth',2.5)
plot((q_r_sep_i+q_r_sep_e)/1e6-(q_t_sol_i+q_t_sol_e)/1e6,'g','linewidth',2.5)
%ylim([1 1e5])
xlabel ('Alfven time \tau_{\alpha}','fontweight','bold','fontsize',14,'fontname','times')
ylabel('Heat flux q (MW)','fontweight','bold','fontsize',14,'fontname','times');
legend('location','NorthWest', 'q_t onto wall','q_t into divertor','q_t stay in VV')
set(gca,'fontweight','bold','fontsize',14,'fontname','arial')
set(gca,'linewidth',1.5)

figure(8)
clf
set(gca,'position',[0.15 0.2 0.7 0.75])
plot(q_r_sep_e/1e6,'b','linewidth',2),hold
plot(q_t_sol_e/1e6,'--b','linewidth',2)
plot((qer_r_xt(nx-2,:)+qer_p_xt(nx-2,:))/1e6,'--r','linewidth',2)
plot((qe_par_yt(ny-num+1,:)-qe_par_yt(num,:))/1e6,'r','linewidth',2)
% plot((q_r_sep_e)/1e6-(q_t_sol_e)/1e6,'g','linewidth',2)
%ylim([1 1e5])
xlabel ('Alfven time \tau_{\alpha}','fontweight','bold','fontsize',14,'fontname','times')
ylabel('Heat flux q (MW)','fontweight','bold','fontsize',14,'fontname','times');
legend('location','NorthWest', 'q_e across sep.', 'q_e into divertor&wall','q_e onto wall','q_e into divertor','q_e stay in VV')
set(gca,'fontweight','bold','fontsize',14,'fontname','arial')
set(gca,'linewidth',1.5)

figure(4)
clf
set(gca,'position',[0.15 0.2 0.7 0.75])
plot((qir_r_xt(nx-2,:)+qir_p_xt(nx-2,:))/1e6,'r','linewidth',2.5),hold
plot((qi_par_yt(ny-num+1,:)-qi_par_yt(num,:))/1e6,'b','linewidth',2.5)
%ylim([1 1e5])
xlabel ('Alfven time \tau_{\alpha}','fontweight','bold','fontsize',14,'fontname','times')
ylabel('Heat flux q (MW)','fontweight','bold','fontsize',14,'fontname','times');
legend('location','NorthWest', 'q_i onto wall','q_i into divertor')
set(gca,'fontweight','bold','fontsize',14,'fontname','arial')
set(gca,'linewidth',1.5)

figure(5)
clf
set(gca,'position',[0.15 0.2 0.7 0.75])
plot(qir_r_xt(ixseps1,:)/1e6,'--black','linewidth',2),hold
plot(qir_r_cross_xt(ixseps1,:)/1e6,'r','linewidth',2)
plot(qir_r_cov_xt(ixseps1,:)/1e6,'b','linewidth',2)
plot(qir_r_con_xt(ixseps1,:)/1e6,'g','linewidth',2)
%ylim([1 1e5])
xlabel ('Alfven time \tau_{\alpha}','fontweight','bold','fontsize',14,'fontname','times')
ylabel('Heat flux q (MW)','fontweight','bold','fontsize',14,'fontname','times');
h=legend('location','NorthWest', 'q_i EXB total','$$q_i EXB \tilde{n_i}\tilde{T_i}$$', 'q_i EXB convection','q_i EXB conduction')
set(h,'interpreter','latex','fontweight','bold','FontSize',14,'fontname','arial');
set(gca,'fontweight','bold','fontsize',14,'fontname','arial')
set(gca,'linewidth',1.5)

figure(6)
clf
set(gca,'position',[0.15 0.2 0.7 0.75])
plot(qer_r_xt(ixseps1,:)/1e6,'--black','linewidth',2),hold
plot(qer_r_cross_xt(ixseps1,:)/1e6,'r','linewidth',2)
plot(qer_r_cov_xt(ixseps1,:)/1e6,'b','linewidth',2)
plot(qer_r_con_xt(ixseps1,:)/1e6,'g','linewidth',2)
%plot(qer_r_equ_xt(ixseps1,:)/1e6,'m','linewidth',2)
%ylim([1 1e5])
xlabel ('Alfven time \tau_{\alpha}','fontweight','bold','fontsize',14,'fontname','times')
ylabel('Heat flux q (MW)','fontweight','bold','fontsize',14,'fontname','times');
h=legend('location','NorthWest', 'q_e EXB total','$$q_e EXB \tilde{n_i}\tilde{T_i}$$', 'q_e EXB convection','q_e EXB conduction','q_e EXB equilibrium')
set(h,'interpreter','latex','FontSize',14);
set(gca,'fontweight','bold','fontsize',14,'fontname','arial')
set(gca,'linewidth',1.5)

figure(7)
clf
set(gca,'position',[0.15 0.2 0.7 0.75])
plot(qer_r_xt(ixseps1,:)/1e6,'black','linewidth',2),hold
plot(qer_p_xt(ixseps1,:)/1e6,'r','linewidth',2)
%ylim([1 1e5])
xlabel ('Alfven time \tau_{\alpha}','fontweight','bold','fontsize',14,'fontname','times')
ylabel('Heat flux q (MW)','fontweight','bold','fontsize',14,'fontname','times');
h=legend('location','NorthWest', 'q_e EXB','q_e magnetic flutter')
set(h,'interpreter','latex','FontSize',14);
set(gca,'fontweight','bold','fontsize',14,'fontname','arial')
set(gca,'linewidth',1.5)

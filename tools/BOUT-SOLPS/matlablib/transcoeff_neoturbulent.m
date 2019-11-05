clear all
clc
addpath('~/bout/tools/matlablib');

%loaddata
% 
Var='R'
S=load(strcat(Var,'.mat'));
eval([Var,'=S.',Var,';']);
%%%%%%%load calculated neoclassical transport coefficient 
Var='dneo'
S=load(strcat(Var,'.mat'));
eval([Var,'=S.',Var,';']);
% 
Var='xii'
S=load(strcat(Var,'.mat'));
eval([Var,'=S.',Var,';']);
% 
 Var='xie'
 S=load(strcat(Var,'.mat'));
 eval([Var,'=S.',Var,';']);
%%%%%%%%%%load transport coefficient calculated by BOUT++ turbulent simulation
Var='Di_x_WO_flutter'
S=load(strcat(Var,'.mat'));
Di_x=S.Di_x;
%
Var='kaii_x_WO_flutter'
S=load(strcat(Var,'.mat'));
kaii_x=S.kaii_x;
%
 Var='kaie_x_WO_flutter'
 S=load(strcat(Var,'.mat'));
kaie_x=S.kaie_x;


%%%%%%%%%%%%%%%%%%smooth the neoclassical transport coefficient near separatrix because of large peak at sep.
ixseps1=179

dneotmp=(1:230)';
dneotmp(1:ixseps1-10)=dneo(1:ixseps1-10);
dneotmp(ixseps1-9:230)=dneo(ixseps1+21:260);
Rtmp=(1:230)';
Rtmp(1:ixseps1-10)=R(1:ixseps1-10);
Rtmp(ixseps1-9:230)=R(ixseps1+21:260);
dneo=interp1(Rtmp,dneotmp,R,'PCHIP');

xietmp=(1:230)';
xietmp(1:ixseps1-10)=xie(1:ixseps1-10);
xietmp(ixseps1-9:230)=xie(ixseps1+21:260);
Rtmp=(1:230)';
Rtmp(1:ixseps1-10)=R(1:ixseps1-10);
Rtmp(ixseps1-9:230)=R(ixseps1+21:260);
xie=interp1(Rtmp,xietmp,R,'PCHIP');

xiitmp=(1:230)';
xiitmp(1:ixseps1-10)=xii(1:ixseps1-10);
xiitmp(ixseps1-9:230)=xii(ixseps1+21:260);
Rtmp=(1:230)';
Rtmp(1:ixseps1-10)=R(1:ixseps1-10);
Rtmp(ixseps1-9:230)=R(ixseps1+21:260);
xii=interp1(Rtmp,xiitmp,R,'PCHIP');

Di_x=Di_x+dneo;
kaie_x=kaie_x+xie;
kaii_x=kaii_x+xii;

clear S Var dneotmp xietmp xiitmp Rtmp

%%%%%%%%%%%save sum of BOUT++ turbulent transport coeff. and smoothed neoclassical transport coeff. to file
fid=fopen('transport.txt','wt');

fprintf(fid,'%s\n%s%s\n','&TRANSPORT','ndata(1, 1 , 1 )= 87' ,',');

 count=0
 for i=150 :2: 250
   count=count+1; 
    fprintf(fid,'%s%u%s%f%s%u%s%f%s\n','tdata(1,',count,', 1, 1 )=',R(i),', tdata(2,',count,', 1, 1)=',Di_x(i),',');
 end

fprintf(fid,'%s%s\n','ndata(1, 3 , 1 )= 87' ,',');
 
count=0
 for i=150 :2: 250
     count=count+1;
	fprintf(fid,'%s%u%s%f%s%u%s%f%s\n','tdata(1,',count,', 3, 1 )=',R(i),', tdata(2,',count,', 3, 1)=',kaii_x(i),',');
end

fprintf(fid,'%s%s\n','ndata(1, 4 , 1 )= 87' ,',');

count=0
 for i=150 :2: 250
count=count+1;
fprintf(fid,'%s%u%s%f%s%u%s%f%s\n','tdata(1,',count,', 4, 1 )=',R(i),', tdata(2,',count,', 4, 1)=',kaie_x(i),',');
end

fclose(fid)



% all the transport data relialbe
% count=0
%  for i=1 :3: 260
%    count=count+1; 
%     fprintf(fid,'%s%u%s%f%s%u%s%f%s\n','tdata(1,',count,', 1, 1 )=',R(i),', tdata(2,',count,', 1, 1)=',Di_x(i),',');
%  end
% 
% fprintf(fid,'%s%s\n','ndata(1, 3 , 1 )= 87' ,',');
%  
% count=0
%  for i=1 :3: 260
%      count=count+1;
% 	fprintf(fid,'%s%u%s%f%s%u%s%f%s\n','tdata(1,',count,', 3, 1 )=',R(i),', tdata(2,',count,', 3, 1)=',kaii_x(i),',');
% end
% 
% fprintf(fid,'%s%s\n','ndata(1, 4 , 1 )= 87' ,',');
% 
% count=0
%  for i=1 :3: 260
% count=count+1;
% fprintf(fid,'%s%u%s%f%s%u%s%f%s\n','tdata(1,',count,', 4, 1 )=',R(i),', tdata(2,',count,', 4, 1)=',kaie_x(i),',');
% end
% 
% fclose(fid)

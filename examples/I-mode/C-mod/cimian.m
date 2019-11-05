a='cmod1.grd.nc';
RXY=ncread(a,'Rxy');
ZXY=ncread(a,'Zxy');
psi=ncread(a,'psixy');
axis=ncread(a,'psi_axis');
bndry=ncread(a,'psi_bndry');
pol_angle=ncread(a,'pol_angle');
hthe=ncread(a,'hthe');
Bpxy=ncread(a,'Bpxy');
x=(psi-axis)./(bndry-axis);
%contour(RXY,ZXY,x,20,'ShowText','on');
%surf(RXY(:,1:68),ZXY(:,1:68),x(:,1:68),'edgecolor','none')
contour(RXY(:,:),ZXY(:,:),x(:,:),10,'ShowText','on');
hold on;
RXY1=[RXY(64,:)',RXY(1,:)']';
ZXY1=[ZXY(64,:)',ZXY(1,:)']';
x1=[x(64,:)',x(1,:)']';
contour(RXY1(:,:),ZXY1(:,:),x1(:,:),10,'ShowText','off');
xlabel('R/m');
ylabel('Z/m');



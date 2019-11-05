function plotbis(bi,bi_phase,flim,nfft)
%for bispecsum

if nargin<4
    nfft=1000;
end

if nargin<3
    flim=200e3;
end

[fx,fy]=bisaxis(nfft);
index1=findout(fx(1,:),0):findout(fx(1,:),flim);
index2=findout(fy(:,1),-flim):findout(fy(:,1),flim);

figure
pcolor(fx(index2,index1),fy(index2,index1),full(bi(index2,index1)));
shading interp
axis equal
axis tight
colorbar2
xlabel2('f_1(Hz)')
ylabel2('f_2(Hz)')
title2('bi^{2}(f_1,f_2)')
set(gca,'layer','top')

figure
pcolor(fx(index2,index1),fy(index2,index1),phaseturn(full(bi_phase(index2,index1)))/pi);
shading interp
axis equal
axis tight
colorbar2
xlabel2('f_1(Hz)')
ylabel2('f_2(Hz)')
title2('\Delta\theta(f_1,f_2)/\pi')
set(gca,'layer','top')
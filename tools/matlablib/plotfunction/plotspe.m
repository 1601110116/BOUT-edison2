function [SS,KK,FF]=plotspe(S,dr,frange,N,nfft,Fs)
%plot s(k/f),no s(k,f)
%dr: cm

if nargin<4;  N=15;end
if nargin<5;  nfft=1000;end
if nargin<6;  Fs=1e6;end

[K,F]=meshgrid((-N:N)*2*pi/(dr*(2*N+1)),(0:nfft/2)*Fs/nfft);
findex=(floor(frange(1)*nfft/Fs)+1):(floor(frange(2)*nfft/Fs)+1);
KK=K(findex,:);
FF=F(findex,:);

tem=sum(S,2);
tem=tem(:,ones(1,2*N+1));
S2=S./tem;
S1=S(findex,:);
SS{1}=S1;
SS{2}=S2(findex,:);

figure
[C,h]=contourf(KK,FF,S1/max(max(S1)),35);


set(h,'linestyle','none')

colorbar4
% set(gca,'layer','top')
xlabel2('k_{r}(cm^{-1})')
ylabel2('Frequency(kHz)')

% title2('log(s(k_{r},f))')
title2('s(k_{r},f)')

% figure
% contour(KK,FF,S2(findex,:),25);
% 
% colorbar
% set(gca,'layer','top')
% xlabel2('k(cm^{-1})')
% ylabel2('Frequency(Hz)')
% title2('s(k/f)')
function [SS,KK1,KK2]=plotspe2b(S,dp,dr,fre,N,nfft,Fs)
%dp:   the poloidal distance    cm
%dr:   the radial distance      cm

if nargin<5;  N=15;end
if nargin<6;  nfft=1000;end
if nargin<7;  Fs=1e6;end

if length(N)==1
    N1=N;
    N2=N;
elseif length(N)==2
    N1=N(1);
    N2=N(2);
end

temp=S;
S=temp([1:N1-1 N1+1:2*N1+1],[1:N2-1 N2+1:2*N2+1],:);

[KK1,KK2]=meshgrid(([-N2:-1 1:N2])*2*pi/(dr*(2*N2+1)),([-N1:-1 1:N1])*2*pi/(dp*(2*N1+1)));

findex=round(fre*nfft/Fs)+1;
figure
SS=squeeze(sum(S(:,:,findex),3)/length(findex));

contourf(KK1,KK2,SS/sum(sum(SS)));

linecross
colorbar2
xlabel2('k_{r}(cm^{-1})')
ylabel2('k_{\phi}(cm^{-1})')
if length(findex)==1
    title2(['f=' num2str((findex-1)*Fs/(1000*nfft)) 'kHz']);
else
    title2(['f=' num2str((findex(1)-1)*Fs/(1000*nfft)) '-' num2str((findex(end)-1)*Fs/(1000*nfft)),'kHz']);
end
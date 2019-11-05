function [bisum,n,bi,bi_phase] = bispecmode(x,y,z,nfft,maxgrid,nmrange,overlap,SideFlag)
% BISPECMODE;This function is used to calculate the bispectrum in the
% k-space or n mode number;
% SYNTAX
% [bisum,n,bi,bi_phase] = bispecmode(x,y,z,nfft,maxgrid,nmrange,overlap,SideFlag)
%
% INPUT PARAMETERS
% - x,y,z: 2D array of the same size. [m,n]=size(x). Calcualte the nolinear 
% interaction in toroidal mode number n space, and statistically averaged by
% the bicoherence on the same flux at different poloidal positions. For the
% auto bispectrum, x=y=z. For cross bispectrum study, please see the work
% by C. Holland in reference[2] below.
% - nfft: grid number in one period
% - maxgrid: total grid number in toroidal direction
% - nmrange: forget it, set it as default value;
% - overlap: set to be 0; 
%
% OUTPUT PARAMETERS
% bisum: total bicoherence
% n: mode nubmer 
% bi: 2D bicoherence or squared normalized bispectrum
% bi_phase: some times get interest meaning.
% 
% DESCRIPTION
% Exhaustive and long description of functionality of bispecmode.
%
% Examples:
% description of example for bispecmode
% >> x=squeeze(psi(300,20:50,1:65,300));
% >> [bisum,n,bi,bi_phase] = bispecmode(x,x,x,65,65*5,[1,33],0,1);
%
% See also:
% bispecfrq, bifrqline, bimodeline
%
% References:
% [1] Young C. Kim and Edward J. Powers. Digital Bispectral Analysis and Its
% Applications to Nonlinear Wave Interactions. IEEE. Transactions on Plasma
% Science PS-7(2), 120–131 (1979).
% [2] C. Holland, G. R. Tynan, P. H. Diamond et al. Evidence for Reynolds-stress
% driven shear flows using bispectral analysis: theory and experiment. Plasma
% Phys. Control. Fusion 44(5A), A453–A457 (2002).
% [3]The function is first coded by Adi Liu in USTC, China

% Copyright (c) ASIPP 1985-2015 Defeng Kong
% All Rights Reserved.
% email: dfkong@ipp.ac.cn
% $Revision: 1.0$ Created on: 28-May-2015 11:32:03


narginchk(4,8);
indsize=size(x);

if length(indsize)~=2,
error('x should be a 2D array with m=poloidal and n=torodial');
end
if nargin < 4, 
    nfft = indsize(2); 
end
if nargin < 5, 
    maxgrid = nfft; 
    warning('ZPERIOD=1, #Fraction of a torus to simulate in BOUT.inp');
end
if nargin < 6, 
    nmrange=[1,floor(nfft/2)];
end
if nargin < 7, 
    overlap = 0; 
    warning('In most of the time , overlap should be set 0');
end
if nargin < 8, 
    SideFlag = 1; 
end

if mod(indsize(2),2)==1;
    x=x(:,1:end-1);
    y=y(:,1:end-1);
    z=z(:,1:end-1);
    nfft=nfft-1;
end

xfft = separate(x, nfft, overlap);
yfft = separate(y, nfft, overlap);
zfft = separate(z, nfft, overlap);

n1=max(floor(nmrange(1)),1);
n2=floor(nmrange(2));

switch SideFlag
    case 1
        [bi,bi_phase,bisum] = bi2_f1Largerf2(xfft,yfft,zfft, [n1 n2], nfft);
    case -1
        [bi,bi_phase,bisum] = bi2_f2Largerf1(xfft,yfft,zfft, [n1 n2], nfft);
    case 0
        [bi,bi_phase,bisum] = bi2_whole(xfft,yfft,zfft, [n1 n2], nfft);
end
if mod(indsize(2),2)==1;
    n=(n1:n2)*maxgrid/(nfft+1);
else
    n=(n1:n2)*maxgrid/nfft;
end
if nargout==0
    plotbi(bi,n);
end

%--------------------------------------------------------------------------
function [bi,bi_phase,bisum] = bi2_f1Largerf2(xfft,yfft,zfft, range, nfft)
bi=zeros(nfft+1,nfft/2+1);
bi_phase=zeros(nfft+1,nfft/2+1);
bisum = zeros(range(2)-range(1)+1, 1);
for ii = range(1) : range(2)
    index = ii - range(1) +1;
    list=[max(fix(ii/2),1):ii-1 ii+1:nfft/2];
    num=length(list);
    for h=1:num
        jj=list(h);
        choo=(ii>jj)*yfft(abs(ii -jj),:)+(ii<jj)*conj(yfft(abs(jj - ii),:));
        [temp,temp_phase] = bi20(xfft(jj,:), choo, zfft(ii,:));
        bi(ii-jj+1+nfft/2,jj+1)=temp;
        bi_phase(ii-jj+1+nfft/2,jj+1)=temp_phase;
        bisum(index)=bisum(index)+temp;
    end
    bisum(index)=bisum(index)/num;
end
bi=sparse(bi);
bi_phase=sparse(bi_phase);

function [bi,bi_phase,bisum] = bi2_f2Largerf1(xfft,yfft,zfft, range, nfft)
bi=zeros(nfft/2+1,nfft+1);
bi_phase=zeros(nfft/2+1,nfft+1);
bisum = zeros(range(2)-range(1)+1, 1);
for ii = range(1) : range(2)
    index = ii - range(1) +1;
    list=[max(fix(ii/2),1):ii-1 ii+1:nfft/2];
    num=length(list);
    for h=1:num
        jj=list(h);
        choo=(ii>jj)*xfft(abs(ii-jj),:)+(ii<jj)*conj(xfft(abs(jj - ii),:));
        [temp,temp_phase] = bi2(choo, yfft(jj,:),  zfft(ii,:));
        bi(jj+1,ii-jj+1+nfft/2)=temp;        
        bi_phase(jj+1,ii-jj+1+nfft/2)=temp_phase;
        bisum(index)=bisum(index)+temp;
    end
    bisum(index)=bisum(index)/num;
end
bi=sparse(bi);
bi_phase=sparse(bi_phase);

function [bi,bi_phase,bisum] = bi2_whole(xfft, yfft, zfft, range, nfft)
bi= zeros(nfft+1, nfft+1);
bi_phase= zeros(nfft+1, nfft+1);
bisum = zeros(range(2)-range(1)+1, 1);
for ii = range(1) : range(2)
    index = ii - range(1) +1;
    list=nfft/2:-1:ii-nfft/2;
    list(list==ii)=[];
    list(list==0)=[];
    num=length(list);
    for h=1:num
        jj=list(h);%f2(yaxis)
        kk=ii-jj;%f1(xaxis)
        choo1=(kk>0)*xfft(abs(kk),:)+(kk<0)*conj(xfft(abs(kk),:));
        choo2=(jj>0)*yfft(abs(jj),:)+(jj<0)*conj(yfft(abs(jj),:));
        [temp,temp_phase] = bi20(choo1, choo2, zfft(ii,:));
        bi(jj+1+nfft/2,kk+1+nfft/2) = temp;
        bi_phase(jj+1+nfft/2,kk+1+nfft/2) = temp_phase;
        bisum(index)=bisum(index)+temp;
    end
    bisum(index)=bisum(index)/num;
end
bi=sparse(bi);
bi_phase=sparse(bi_phase);

function y = separate(x, nfft, overlap)
xind=size(x);
x0=reshape(x',xind(1)*xind(2),1);
xmat=divmat(x0,nfft,hanning(nfft),overlap);
temp=fftmat(xmat);
y=temp(2:nfft/2+1,:);

function y=fftmat(x,win)
nfft=size(x,1);
if nargin==1
    win=hanning(nfft);
end
temp=fft(x,nfft);
y=temp/norm(win);

function smat=divmat(s,n,win,overlap)
%divide a vector s to a nxm matrix
%always n refer to the length of Forier window(NFFT)
%       m refer to the number of realizations under the given NFFT
%       
%the default of overlap is n/2;
%--------------------------------------------------------------------------
%               best example of how to use function 'ones'!!!!!
%--------------------------------------------------------------------------
if nargin<4
    overlap=floor(n/2);
end
if nargin<3
    win=hanning(n);
elseif strcmp(win,'none')
    win=rectwin(n);
end

%s=detrend(s);
s=s(:);
len=length(s);
m=floor((len-overlap)/(n-overlap));
colindex=1 + (0:(m-1))*(n-overlap);
rowindex=(1:n)';

smat=s(rowindex(:,ones(1,m))+colindex(ones(n,1),:)-1);
smat=win(:,ones(1,m)).*detrend(smat);

function [y,an]= bi20(X, Y, Z)
temp1 = X.*Y;
temp2 = conj(Z)*temp1.';
an=angle(temp2);
y=abs(temp2)^2/((conj(X)*X.')*(conj(Y)*Y.')*(conj(Z)*Z.'));

function plotbi(bi,n)
% this function is used to give a contourf plot of bi
n1=[0,n];
n2=linspace(-n(end),n(end),2*length(n)+1);
[N1,N2]=meshgrid(n1,n2);
figure;
[~,h]=contourf(N1,N2,bi,30);
    set(h,'linestyle','none');
    h_axes=findobj(gcf,'type','axes');
    set(h_axes,'LineWidth',2);
    set(h_axes,'fontsize',14);
    set(h_axes,'fontweight','bold');
    set(h_axes,'XMinorTick','on');
    set(h_axes,'YMinorTick','on');
    set(h_axes,'layer','top');
xlabel('n_1','interpreter','tex','fontname','arial','fontsize',16);
ylabel('n_2','interpreter','tex','fontname','arial','fontsize',16);
title('Bicoherence','interpreter','tex','fontname','arial','fontsize',16);



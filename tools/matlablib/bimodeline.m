function [bis,bis_phase,nf] = bimodeline(x, y, z, nfreq, nfft, maxgrid, overlap)
% BIMODELINE: calculate the bicoherence between the mode with mode number
% 'nfreq' and other modes of different mode numbers
% SYNTAX
% [bis,bis_phase,nf] = bimodeline(x, y, z, nfreq, nfft, maxgrid, overlap)
%
% INPUT PARAMETERS
% - x,y,z: 2D array of the same size. [m,n]=size(x). Calcualte the nolinear 
% interaction in toroidal mode number n space, and statistically averaged by
% the bicoherence on the same flux at different poloidal positions. For the
% auto bispectrum, x=y=z. For cross bispectrum study, please see the work
% by C. Holland in reference[2] below.
% - nfreq: mode number you interest,ex. nfreq=5
% - nfft: grid number in one period
% - maxgrid: total grid number in toroidal direction
% - overlap: set to be 0; 
%
% OUTPUT PARAMETERS
% - bis: bicoherence 
% - bis_phase: cross phase of the bispectrum
% - nf: real mode number index
%
% DESCRIPTION
% Exhaustive and long description of functionality of bimodeline.
%
% Examples:
% description of example for bimodeline
% >> x=squeeze(psi(300,20:50,1:65,300));
% >> [bis,bis_phase,nf] = bimodeline(psi0, psi0, psi0, 5, 65, 65*5,0);
%
% See also:
% bispecfrq, bispecmode, bifrqline
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
% $Revision: 1.0$ Created on: 01-Jun-2015 15:49:30
narginchk(4,7);
indsize=size(x);

if length(indsize)~=2,
error('x should be a 2D array with m=poloidal and n=torodial');
end

if nargin < 5, 
    nfft = indsize(2); 
end
if nargin < 6, 
    maxgrid = nfft; 
    warning('ZPERIOD=1, #Fraction of a torus to simulate in BOUT.inp');
end
if nargin < 7, 
    overlap = 0; 
    warning('In most of the time , overlap should be set 0');
end

if mod(indsize(2),2)==1;
    x=x(:,1:end-1);
    y=y(:,1:end-1);
    z=z(:,1:end-1);
    Fs=maxgrid-maxgrid/nfft;
    nfft=nfft-1;  
end

xind=size(x);
x0=reshape(x',xind(1)*xind(2),1);
temp=fftmat(divmat(x0,nfft,hanning(nfft),overlap));
xfft=temp(2:nfft/2+1,:);

yind=size(y);
y0=reshape(y',yind(1)*yind(2),1);
temp=fftmat(divmat(y0,nfft,hanning(nfft),overlap));
yfft=temp(2:nfft/2+1,:);

zind=size(z);
z0=reshape(z',zind(1)*zind(2),1);
temp=fftmat(divmat(z0,nfft,hanning(nfft),overlap));
zfft=temp(2:nfft/2+1,:);

n=max(floor(nfreq*nfft/Fs),1);
if length(nfreq)==1
     [bis,bis_phase,nf] = bi2_f1Largerf2(xfft,yfft,zfft, n, nfft,Fs/nfft);
else
    bis=cell(length(n),1);
    bis_phase=cell(length(n),1);
    nf=cell(length(n),1);
    for i=1:length(n)
        [bi,bi_phase,fr] = bi2_f1Largerf2(xfft,yfft,zfft, n(i), nfft,Fs/nfft);
        bis{i}=bi;
        bis_phase{i}=bi_phase;
        nf{i}=fr;
    end
end

function [bi,bi_phase,f] = bi2_f1Largerf2(xfft,yfft,zfft,ii , nfft,K)
    list=[max(fix(ii/2),1):ii-1 ii+1:nfft/2];
    num=length(list);
    bi=zeros(num,1);
    bi_phase=zeros(num,1);
    f=zeros(num,1);
    for h=1:num
        jj=list(h);
        choo=(ii>jj)*yfft(abs(ii -jj),:)+(ii<jj)*conj(yfft(abs(jj - ii),:));
        [temp,temp_phase] = bi2(xfft(jj,:), choo, zfft(ii,:));
        bi(h)=temp;
        bi_phase(h)=temp_phase;
        f(h)=jj*K;
    end
    
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

function [y,an]= bi2(X, Y, Z)
temp1 = X.*Y;
temp2 = conj(Z)*temp1.';
an=angle(temp2);
y=abs(temp2)^2/((conj(temp1)*temp1.')*(conj(Z)*Z.'));



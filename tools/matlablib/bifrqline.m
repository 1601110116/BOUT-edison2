function [bis,bis_phase,f] = bifrqline(x, y, z, freq, nfft, Fs, overlap)
% BIFRQLINE: calculate the bicoherence between the mode with frequency
% 'freq' and other modes of different frequency
% SYNTAX
% [bis,bis_phase,f] = bifrqline(x, y, z, freq, nfft, Fs, overlap)
%
% INPUT PARAMETERS
% - x,y,z: 1D array of the same size. [m,n]=size(x). Calcualte the nolinear 
% interaction in frequency space, and statistically averaged by
% the bicoherence on the same flux at different poloidal positions. For the
% auto bispectrum, x=y=z. For cross bispectrum study, please see the work
% by C. Holland in reference[2] below.
% - freq: the frequency of the mode you are interest;
% - nfft: number point of one ensemble
% - Fs: Sampling rate
%
% OUTPUT PARAMETERS
% - bis: bicoherence 
% - bis_phase: cross phase of the bispectrum
%
% DESCRIPTION
% Exhaustive and long description of functionality of bifrqline.
%
% Examples:
% description of example for bifrqline
% >> t=0:1E-6:0.1; % Fs=1e6;
% >> mkarr=cos(2*pi*20e3*t+0.3)+cos(2*pi*50e3*t+0.5)+0.5*cos(2*pi*70e3*t+0.4)+cos(2*pi*20e3*t+0.3).*cos(2*pi*50e3*t+0.5)+rand(1,length(t));
% >> [bis,bis_phase,f] = bifrqline(mkarr, mkarr, mkarr, 20e3, 1000, 1e6, 500);
%
% See also:
% bispecfrq, bispecmode, bimodeline
%
% References:
% [1] Young C. Kim and Edward J. Powers. Digital Bispectral Analysis and Its
% Applications to Nonlinear Wave Interactions. IEEE. Transactions on Plasma
% Science PS-7(2), 120–131 (1979).
% [2] C. Holland, G. R. Tynan, P. H. Diamond et al. Evidence for Reynolds-stress
% driven shear flows using bispectral analysis: theory and experiment. Plasma
% Phys. Control. Fusion 44(5A), A453–A457 (2002).
% [3]The function is first coded by Adi Liu in USTC, China

% Copyright (c) ASIPP 1958-2015 Defeng Kong
% All Rights Reserved.
% email: dfkong@ipp.ac.cn
% $Revision: 1.0$ Created on: 01-Jun-2015 15:21:43

narginchk(4,8);
if nargin < 5, nfft = 100; end
if nargin < 6, Fs = 5.7e6; end
if nargin < 7, overlap = floor(nfft*3/4); end

if mod(nfft,2)==1;
    nfft=nfft-1;
    warning('nfft should be an even number');
end

temp=fftmat(divmat(x,nfft,hanning(nfft),overlap));
xfft=temp(2:nfft/2+1,:);
temp=fftmat(divmat(y,nfft,hanning(nfft),overlap));
yfft=temp(2:nfft/2+1,:);
temp=fftmat(divmat(z,nfft,hanning(nfft),overlap));
zfft=temp(2:nfft/2+1,:);

n=max(floor(freq*nfft/Fs),1);

if length(freq)==1
     [bis,bis_phase,f] = bi2_f1Largerf2(xfft,yfft,zfft, n, nfft,Fs/nfft);
else
    bis=cell(length(n),1);
    bis_phase=cell(length(n),1);
    f=cell(length(n),1);
    for i=1:length(n)
        [bi,bi_phase,fr] = bi2_f1Largerf2(xfft,yfft,zfft, n(i), nfft,Fs/nfft);
        bis{i}=bi;
        bis_phase{i}=bi_phase;
        f{i}=fr;
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

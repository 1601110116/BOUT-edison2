function [fx,fy]=bisaxis(nfft, Fs, SideFlag)
%get the grid matrixes of axis

if nargin < 1, nfft = 1000; end
if nargin < 2, Fs = 1e6; end
if nargin < 3, SideFlag = 1; end

switch SideFlag
    case 1
        xaxis=linspace(0,Fs/2,nfft/2+1);
        yaxis=linspace(-Fs/2,Fs/2,nfft+1);
%!!!! in bi(i,j) i infer to y-axis and j infer to x-axis,it is mistakable here. 
    case -1
        yaxis=linspace(0,Fs/2,nfft/2+1);
        xaxis=linspace(-Fs/2,Fs/2,nfft+1);
    case 0
        xaxis=linspace(-Fs/2,Fs/2,nfft+1);
        yaxis=linspace(-Fs/2,Fs/2,nfft+1);
end

[fx,fy]=meshgrid(xaxis,yaxis);
function y=fftmat(x,win)
nfft=size(x,1);
if nargin==1
    win=hanning(nfft);
end
temp=fft(x,nfft);
y=temp/norm(win);
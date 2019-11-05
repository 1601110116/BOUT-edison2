function smat=divmat(s,n,win,overlap)
%divide a vector s to a n¡Ám matrix
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
%smat=win(:,ones(1,m)).*smat;
function smat=divmat2(s,n,overlap)
%divide a vector s to a n¡Ám matrix

if nargin<3
    overlap=0;
end

s=s(:);
len=length(s);
m=floor((len-overlap)/(n-overlap));
colindex=1 + (0:(m-1))*(n-overlap);
rowindex=(1:n)';
smat=s(rowindex(:,ones(1,m))+colindex(ones(n,1),:)-1);

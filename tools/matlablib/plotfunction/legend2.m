function legend2(varargin)

len=length(varargin);
cellnew=cell(1,len);
for i=1:len
    temp=varargin{i};
    cellnew{i}=['$$' temp '$$'];
end
h=legend(cellnew);
set(h,'interpreter','latex','edgecolor','w','fontsize',14);
function ytick2(h_axe)

if nargin==0
    h_axe=gca;
end

xarray=get(h_axe,'ytick')/1e3;
xlab={num2str(xarray')};
set(h_axe,'yticklabel',xlab);
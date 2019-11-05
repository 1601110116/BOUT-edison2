function xtick2(h_axe)

if nargin==0
    h_axe=gca;
end

xarray=get(h_axe,'xtick')/1e3;
xlab={num2str(xarray')};
set(h_axe,'xticklabel',xlab);
set(h_axe,'fontsize',13);
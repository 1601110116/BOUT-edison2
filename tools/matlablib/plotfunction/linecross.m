function linecross(h)

if nargin==0
    h=gcf;
end

h_axe=findobj(gcf,'type','axes','-and','tag','');
yl=get(h_axe,'ylim');
xl=get(h_axe,'xlim');
line([0 0],yl,'linestyle','--','color','k','linewidth',1.5);
line(xl,[0 0],'linestyle','--','color','k','linewidth',1.5);
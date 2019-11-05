function xlim2(range)

h=findobj(gcf,'type','axes');
for i=1:length(h)
    set(h(i),'xlim',range);
end
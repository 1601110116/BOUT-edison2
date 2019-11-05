function ylim2(range)

h=findobj(gcf,'type','axes');
for i=1:length(h)
    set(h(i),'ylim',range);
end
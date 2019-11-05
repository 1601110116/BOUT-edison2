sta={'log','linear'};

h=findobj(gcf,'type','axes');
for i=1:length(h)
    temp=get(h(i),'yscale');
    k=strmatch(temp, strvcat('log','linear'));
    set(h(i),'yscale',sta{3-k});
end
clear h i k sta temp
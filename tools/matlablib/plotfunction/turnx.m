sta={'log','linear'};

h=findobj(gcf,'type','axes');
for i=1:length(h)
    temp=get(h(i),'xscale');
    k=strmatch(temp, strvcat('log','linear'));
    set(h(i),'xscale',sta{3-k});
end
clear h i k sta temp
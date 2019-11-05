function copyaxe(axe1,axe2)

hh=get(axe1,'children');

yli=get(axe1,'ylim');
xli=get(axe1,'xlim');

for i=1:length(hh)
    copyobj(hh(end-i+1),axe2);
end

set(axe2,'xlim',xli);
set(axe2,'ylim',yli);
gca=axe2;
box on
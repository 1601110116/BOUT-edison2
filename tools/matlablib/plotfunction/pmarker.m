function pmarker(hh)

% if nargin<=1
%     h=gco;
% end
% 
% if nargin==0
%     co=get(h,'color');
% end
% 
% if ~strcmp(get(h,'type'),'line')
%     error('wrong handle!')
% end
% 
% set(h,'markerfacecolor',co,'markeredgecolor',co);

if nargin==1
    hh=gco;
elseif nargin==0
    hh=findobj(gcf,'type','line','-not','marker','none');
end

for i=1:length(hh)
    %co=get(hh(i),'color');
    set(hh(i),'markerfacecolor','w');
end

    
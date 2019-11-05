function legend3(varargin)

% len=length(varargin);
% cellnew=cell(1,len);
% for i=1:len
%     temp=varargin{i};
%     cellnew{i}=['$$' temp '$$'];
% end


h=legend(varargin);
set(h,'interpreter','latex','edgecolor','w','fontname','arial','fontsize',13);
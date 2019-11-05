function lineall(n, fig)
% LINEALL plot lines on all stackplot obj.
% SYNTAX
% LINEALL(n, fighandle) : plot n lines on figure, by default n=1; 
% LINEALL(n) : handle can be the handle of a figure, by default
%                   handle=gcf.
%
% Examples:
% plot 2 lines on current stackplot obj.
% >> lineall(gcf, 2);
%
% See also:
% setptr, moveptr
%
% References:
% [1] A. Einstein, Die Grundlage der allgemeinen Relativitaetstheorie,
%     Annalen der Physik 49, 769-822 (1916).

% Copyright (c) CAS Key Laboratory of Basic Plasma Physics, USTC 1958-2012
% Author: lantao
% Email: lantao@ustc.edu.cn
% All Rights Reserved.
% $Revision: 1.0$ Created on: 11-Sep-2012 15:06:26

error(nargchk(0,2,nargin));   % check input arguments

if nargin<2, fig = gcf; end  % default = gcf
if nargin<1, n = 1; end   % default = 1

if ~ishandle(fig)
    error('first input argument should be figure handle!');
end

% linestyle
linestyle = '--';

% linecolor
linecolor = 'r';

% linewidth
linewidth = 2.0;

% tag
tag = 'lineall';

try
    figure(fig); % focus on fig
    
    
    % plot
    allaxes = findallaxes(fig);
    for ii=1:n
        [x,y] = ginput(1);
        hg = hggroup;
        for jj=1:length(allaxes)
            XX = [x x];
            YY = ylim(allaxes(jj));
            line(XX, YY, 'linestyle', linestyle, 'color',linecolor, 'linewidth', linewidth,'parent', allaxes(jj), 'tag', tag);
        end
    end

catch  E1

    rethrow(E1);
end


function Axes = findallaxes(fh)
temp = findobj(fh,'Type','axes');
Axes = temp;
% remove the colorbar axes
for ii=1:length(Axes)
    if strcmpi(get(Axes(ii),'Tag'), 'Colorbar')
        Axes(ii)= [];
    end
end
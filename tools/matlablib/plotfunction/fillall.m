function fillall(X,c)
% set a color field to a figure area
% INPUT:
%   X:  the color field of X axis
%   c:  color,such as 'r' 
% 2010/1/14  D.F.KONG

h_axes=findobj(gcf,'type','axes');
for i=1:length(h_axes)

    Y=get(h_axes(i),'ylim');
    axes(h_axes(i));
    fill2([X(1),X(2),X(2),X(1)],[Y(1),Y(1),Y(2),Y(2)],c);

end
    


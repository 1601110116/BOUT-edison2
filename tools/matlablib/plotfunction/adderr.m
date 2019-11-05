function [theta_err,x]=adderr(N,figure_h)
%N is the  realization number
%phase error for function 'PC1' and 'PC2'

if nargin==1
    figure_h=gcf;
end

h_axes=findobj(figure_h,'type','axes');
h_plot=findobj('parent',h_axes(2),'-not','color','k');
co=get(h_plot,'ydata');
x=get(h_plot,'xdata');

if iscell(co)
    error('wrong figure!!!')
else
    theta_err=sqrt(1-co.^2)./(abs(co)*sqrt(2*N));
end

h_th=findobj('parent',h_axes(1),'type','line','-not','color','k');
y=get(h_th,'ydata');
x=get(h_th,'xdata');

set(h_axes(1),'nextplot','add');
set(figure_h,'currentaxes',h_axes(1))
plot(x,y+theta_err/pi,'--r',x,y-theta_err/pi,'--r')


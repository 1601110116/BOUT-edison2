function adderr2(N,figure_h)
%N is the  realization number
%error for function 'plotbis2'('bispeclin')

if nargin==1
    figure_h=gcf;
end

h_axes=findobj(figure_h,'type','axes');

if length(h_axes)~=2
    error('wrong figure!')
end

h_plot=findobj(h_axes(2),'type','line');
co=get(h_plot,'ydata');
x=get(h_plot,'xdata');
if iscell(co)
    error('wrong figure!')
else
   co_err=(1-co)/N;
end
set(h_axes(2),'nextplot','add');
set(figure_h,'currentaxes',h_axes(2))
plot(x,co_err,'--','linewidth',2.5,'color',[0.502 0.502 0.502])

% theta_err=(1./co-1)/N;
% h_th=findobj('parent',h_axes(1),'type','line');
% y=get(h_th,'ydata');
% x=get(h_th,'xdata');
% 
% set(h_axes(1),'nextplot','add');
% set(figure_h,'currentaxes',h_axes(1))
% plot(x,y+theta_err/pi,'--r',x,y-theta_err/pi,'--r')
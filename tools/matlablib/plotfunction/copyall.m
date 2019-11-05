function copyall(fobj1,fobj2)
%
%copy all plots in figure fobj1 to figure fobj2
%fobj1 is just the number of figure
%i.e.    copyall(1,2)

axh1=findobj(fobj1,'type','axes');
axh2=findobj(fobj2,'type','axes');

if length(axh1)~=length(axh2)
    error('different axes number between two figures')
else
    for i=1:length(axh1)
        plot_h=findobj(axh1(i),'type','line');
        copyobj(plot_h,axh2(i));
    end
end
function title3(str)
% latex

if ischar(str)
    str_new=cstr(str);
    title(str_new,'interpreter','latex','fontsize',14)

elseif iscellstr(str)
    len=length(str);
    for i=1:len
        str_new{i}=cstr(str{i});
    end
   title(str_new,'interpreter','latex','fontsize',14)
end




    function y=cstr(x)
        x(x==' ')='~';
    y=['$$' x '$$'];
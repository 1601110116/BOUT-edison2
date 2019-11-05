function text2(varargin)

if mod(length(varargin),2)==1
    text(varargin{:},'Fontname','arial','fontsize',16);
else
    num=varargin{end};
     text(varargin{1:end-1},'Fontname','arial','fontsize',num);
end
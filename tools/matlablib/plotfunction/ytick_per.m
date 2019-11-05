function ytick_per(h)

if nargin==0
    h=gca;
end

y=get(h,'ytick');
y=y*100;
y=y';

temp=cell(length(y),1);
for i=1:length(y)
    temp{i}=[num2str(y(i)) '%'];
end
set(h,'yticklabel',temp)
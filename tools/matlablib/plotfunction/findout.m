function index=findout(y,y0,num)

if nargin==2
    num=1;
end

index=zeros(num,1);
temp=abs(y-y0);

if ~isempty(find(temp==0, 1))
    index=find(temp==0,num);
else    
    temp(temp==0)=max(temp);
    haha=find(temp==min(temp));
    index(1)=haha(1);
    if num>=2
        for i=2:num
            temp=temp-temp(index(i-1));
            temp(temp==0)=max(temp);
            haha=find(temp==min(temp));
            index(i)=haha(1);
        end
    end
end
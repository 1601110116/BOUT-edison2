function subplot2(m,n)
figure;
for ii=1:n*(m-1)
    subplot(m,n,ii);
    myplot(m,n,ii);
end
for ii=n*(m-1)+1:m*n
    subplot(m,n,ii);
end

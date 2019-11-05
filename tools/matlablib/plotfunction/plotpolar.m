function plotpolar(mat2,mat1,err)
%assistant routine for 'polar2'

index=find(mat1==0);
mat1(index)=[];
mat2(index)=[];
polar2(err,mat2,mat1,'*r');
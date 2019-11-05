//locate the position of x in vector xvec
int locate (BoutReal *xvec, int n, BoutReal x)
{	
	int i,jl,jm,ju;
	bool ascnd;
	ascnd = (xvec[n-1] >= xvec[0]);
	jl=0;
	ju=n;
	while (ju-jl > 1) {
		jm = (ju+jl) >> 1;
		if (x >= xvec[jm] == ascnd)
			jl=jm;
		else
			ju=jm;
	}
	if (x == xvec[0])
		i=0;
	else if (x == xvec[n-1])
		i=n-1;
	else
		i=jl;
	return i;
}

BoutReal interp1dlin ( BoutReal *xarray, int nx, BoutReal *fx, BoutReal x) {
	int i;
	BoutReal f;
	i = locate(xarray, nx, x);
	f = fx[i]+(fx[i+1]-fx[i])*(x-xarray[i])/(xarray[i+1]-xarray[i]);
	return f;
}

//Bilinear interpolation
BoutReal interp2dlin ( BoutReal *xarray, int nx, BoutReal *yarray, int ny, Field2D &fxy, BoutReal x, BoutReal y)
{
	int i,j;
	BoutReal t,u,f;
	i = locate (xarray, nx, x);
	j = locate (yarray, ny, y);
	t = (x-xarray[i])/(xarray[i+1]-xarray[i]);
	u = (y-yarray[j])/(yarray[j+1]-yarray[j]);
	f = (1.-t)*(1.-u)*fxy[i][j] + t*(1.-u)*fxy[i+1][j] + (1.-t)*u*fxy[i][j+1] + t*u*fxy[i+1][j+1];
	return f;
}

BoutReal interpqua ( BoutReal *xarray, int nx, BoutReal *yarray, BoutReal x)
{
	int i,ileft;
	BoutReal x1,x2,x3,y1,y2,y3,y;
	i = locate (xarray, nx, x);
	if(i==0)
		ileft=0;
	else if(i > nx-3)
		ileft=nx-3;
	else if(fabs(x-xarray[i])<fabs(xarray[i+1]-x))
		ileft=i-1;
	else
		ileft=i;
	x1=xarray[ileft];x2=xarray[ileft+1];x3=xarray[ileft+2];
	y1=yarray[ileft];y2=yarray[ileft+1];y3=yarray[ileft+2];
	y = (x-x2)*(x-x3)*y1/(x1-x2)/(x1-x3)+(x-x1)*(x-x3)*y2/(x2-x1)/(x2-x3)+(x-x1)*(x-x2)*y3/(x3-x1)/(x3-x2);
	return y;
}


BoutReal interp2dqua ( BoutReal *xarray, int nx, BoutReal *yarray, int ny, Field2D &fxy, BoutReal x, BoutReal y)
{
	int i,j,ileft,jleft;
	BoutReal x1,x2,x3,y1,y2,y3,tmpy[3],yy;
	i = locate (xarray, nx, x);
	if(i==0)
		ileft=0;
	else if(i > nx-3)
		ileft=nx-3;
	else if(fabs(x-xarray[i])<fabs(xarray[i+1]-x))
		ileft=i-1;
	else
		ileft=i;
	j = locate (yarray, ny, y);
	if(j==0)
		jleft=0;
	else if(j > ny-3)
		jleft=ny-3;
	else if(fabs(y-yarray[j])<fabs(yarray[j+1]-y))
		jleft=j-1;
	else
		jleft=j;
	x1=xarray[ileft];x2=xarray[ileft+1];x3=xarray[ileft+2];
	y1=yarray[jleft];y2=yarray[jleft+1];y3=yarray[jleft+2];
	for(int k=0;k<3;k++)
		tmpy[k]=(x-x2)*(x-x3)*fxy[ileft][jleft+k]/(x1-x2)/(x1-x3) +
			(x-x1)*(x-x3)*fxy[ileft+1][jleft+k]/(x2-x1)/(x2-x3) +
			(x-x1)*(x-x2)*fxy[ileft+2][jleft+k]/(x3-x1)/(x3-x2);
		/*tmpy[k]=(x-x2)*(x-x3)*(*((BoutReal*)fxy + ny*ileft + jleft+k))/(x1-x2)/(x1-x3) +
			(x-x1)*(x-x3)*(*((BoutReal*)fxy + ny*(ileft+1) + jleft+k))/(x2-x1)/(x2-x3) +
			(x-x1)*(x-x2)*(*((BoutReal*)fxy + ny*(ileft+2) + jleft+k))/(x3-x1)/(x3-x2);*/
	yy = (y-y2)*(y-y3)*tmpy[0]/(y1-y2)/(y1-y3)+(y-y1)*(y-y3)*tmpy[1]/(y2-y1)/(y2-y3)+(y-y1)*(y-y2)*tmpy[2]/(y3-y1)/(y3-y2);
	return yy;
}

BoutReal interpcub ( BoutReal *xarray, int nx, BoutReal *yarray, BoutReal x)
{
	int i,ileft;
	BoutReal x1,x2,x3,x4,y1,y2,y3,y4,y;
	i = locate (xarray, nx, x);
	if(i==0)
		ileft=0;
	else if(i > nx-3)
		ileft=nx-4;
	else
		ileft=i-1;
	x1=xarray[ileft];x2=xarray[ileft+1];x3=xarray[ileft+2];x4=xarray[ileft+3];
	y1=yarray[ileft];y2=yarray[ileft+1];y3=yarray[ileft+2];y4=yarray[ileft+3];
	y =	(x-x2)*(x-x3)*(x-x4)*y1/(x1-x2)/(x1-x3)/(x1-x4)+
		(x-x1)*(x-x3)*(x-x4)*y2/(x2-x1)/(x2-x3)/(x2-x4)+
		(x-x1)*(x-x2)*(x-x4)*y3/(x3-x1)/(x3-x2)/(x3-x4)+
		(x-x1)*(x-x2)*(x-x3)*y4/(x4-x1)/(x4-x2)/(x4-x3);
	return y;
}


BoutReal interp2dcub ( BoutReal *xarray, int nx, BoutReal *yarray, int ny, Field2D &fxy, BoutReal x, BoutReal y)
{
	int i,j,ileft,jleft;
	BoutReal x1,x2,x3,x4,y1,y2,y3,y4,tmpy[4],yy;
	i = locate (xarray, nx, x);
	if(i==0)
		ileft=0;
	else if(i > nx-3)
		ileft=nx-4;
	else
		ileft=i-1;
	j = locate (yarray, ny, y);
	if(j==0)
		jleft=0;
	else if(j > ny-3)
		jleft=ny-4;
	else
		jleft=j-1;
	x1=xarray[ileft];x2=xarray[ileft+1];x3=xarray[ileft+2];x4=xarray[ileft+3];
	y1=yarray[jleft];y2=yarray[jleft+1];y3=yarray[jleft+2];y4=yarray[jleft+3];
	for(int k=0;k<4;k++)
		tmpy[k]=(x-x2)*(x-x3)*(x-x4)*fxy[ileft][jleft+k]/(x1-x2)/(x1-x3)/(x1-x4) +
			(x-x1)*(x-x3)*(x-x4)*fxy[ileft+1][jleft+k]/(x2-x1)/(x2-x3)/(x2-x4) +
			(x-x1)*(x-x2)*(x-x4)*fxy[ileft+2][jleft+k]/(x3-x1)/(x3-x2)/(x3-x4) +
			(x-x1)*(x-x2)*(x-x3)*fxy[ileft+3][jleft+k]/(x4-x1)/(x4-x2)/(x4-x3);
		/*tmpy[k]=(x-x2)*(x-x3)*(x-x4)*(*((BoutReal*)fxy + ny*ileft + jleft+k))/(x1-x2)/(x1-x3)/(x1-x4) +
			(x-x1)*(x-x3)*(x-x4)*(*((BoutReal*)fxy + ny*(ileft+1) + jleft+k))/(x2-x1)/(x2-x3)/(x2-x4) +
			(x-x1)*(x-x2)*(x-x4)*(*((BoutReal*)fxy + ny*(ileft+2) + jleft+k))/(x3-x1)/(x3-x2)/(x3-x4) +
			(x-x1)*(x-x2)*(x-x3)*(*((BoutReal*)fxy + ny*(ileft+3) + jleft+k))/(x4-x1)/(x4-x2)/(x4-x3);*/
	yy = (y-y2)*(y-y3)*tmpy[0]/(y1-y2)/(y1-y3)+(y-y1)*(y-y3)*tmpy[1]/(y2-y1)/(y2-y3)+(y-y1)*(y-y2)*tmpy[2]/(y3-y1)/(y3-y2);
	yy = 	(y-y2)*(y-y3)*(y-y4)*tmpy[0]/(y1-y2)/(y1-y3)/(y1-y4)+
		(y-y1)*(y-y3)*(y-y4)*tmpy[1]/(y2-y1)/(y2-y3)/(y2-y4)+
		(y-y1)*(y-y2)*(y-y4)*tmpy[2]/(y3-y1)/(y3-y2)/(y3-y4)+
		(y-y1)*(y-y2)*(y-y3)*tmpy[3]/(y4-y1)/(y4-y2)/(y4-y3);
	return yy;
}

BoutReal interp3d ( BoutReal *xarray, int nx, BoutReal *yarray, int ny, BoutReal *zarray, int nz, Field3D &fxyz, BoutReal x, BoutReal y, BoutReal z) {
	int i,j,k;
	BoutReal xd,yd,zd,c00,c10,c01,c11,c0,c1,c;
	i = locate (xarray, nx, x);
	if(i==nx-1)
		{i=nx-2;}
	j = locate (yarray, ny, y);
	if(j==ny-1)
		{j=ny-2;}
	k = locate (zarray, nz, z);
	if(k==nz-1)
		{k=nz-2;}
	xd = (x-xarray[i])/(xarray[i+1]-xarray[i]);
	yd = (y-yarray[j])/(yarray[j+1]-yarray[j]);
	zd = (z-zarray[k])/(zarray[k+1]-zarray[k]);
	c00=fxyz[i][j][k]*(1.0-xd)     + fxyz[i+1][j][k]*xd;
	c10=fxyz[i][j+1][k]*(1.0-xd)   + fxyz[i+1][j+1][k]*xd;
	c01=fxyz[i][j][k+1]*(1.0-xd)   + fxyz[i+1][j][k+1]*xd;
	c11=fxyz[i][j+1][k+1]*(1.0-xd) + fxyz[i+1][j+1][k+1]*xd;
	c0 =c00*(1.0-yd)+c10*yd;
	c1 =c01*(1.0-yd)+c11*yd;
	c  =c0*(1.0-zd) +c1*zd;
	return c;
}

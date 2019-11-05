#include <iostream>
using namespace std;
#include <cmath>

//this subroutine generate a knot sequnce knot(n+k) for data points x(n), using not-a-kont conditions.
void cdbbsnak(BoutReal *x, int n, BoutReal *knot, int nk) {
	int k;
	k=nk-n;
	for (int i=0;i<k;i++)
		{knot[i]=x[0];knot[n+i]=x[n-1];}
	if (k % 2 == 0)
		for (int i=k;i<n;i++)
			knot[i]=x[i-k/2];
	else
		for (int i=k;i<n;i++)
			knot[i]=0.5*(x[i-(k-1)/2]+x[i-1-(k-1)/2]);
}

//given knot sequent knot, point x which located knot(left)<=x<knot(left+1), output all the B-value val from 1 order to k order.
//i.e. val(1) for 1 order, val(2:3) for 2 order, val((j-1)*j/2+1 : j*(j+1)/2) for j order, val((k-1)*k/2+1 : k*(k+1)/2) for k order.

void cdbbsvalvec(BoutReal *knot, int k, BoutReal x, int left, BoutReal *val)
{	int n1,n2,i,j;
	BoutReal deltar[k],deltal[k],term,temp;

	/*if (size(val) /= k*(k+1)/2) then
		write (*,*) 'error: dimension of val is wrong, cdbbsvalvec'
		STOP 'program terminated by cdbbsvalvec'
	end if*/
	if (knot[left] >= knot[left+1])
		{cout << "error: knot(left) < knot(left+1) is not satisfied, cdbbsvalvec\n";}
	
	val[0]=1.0;
	n1=0;
	for ( j=0;j<k-1;j++) {
		deltar[j] = knot[left+j+1] - x;
		deltal[j] = x - knot[left-j];
		temp=0.0;
		n2=n1+j+1;
		for (i=0;i<=j;i++) {
			term = val[n1+i] / ( deltar[i] + deltal[j-i] );
			val[n2+i] = temp + deltar[i] * term;
			temp = deltal[j-i] * term;
			}
		val[n2+j+1] = temp;
		n1=n2;
		}
}

//cdblocate.cxx is revised from cdblocate.f90
//cdblocate.f90 is revised from inrlocate.f90
//inrlocate can find the left closed and right open interval [ ), in which x located, for an strictly increase or decrease sequence xx.
//But in case of B-spline, the knot sequence may not strictly increase, so inrlocate should be revised.
int cdblocate(BoutReal *xx, int n, BoutReal x)
{	int i,jl,jm,ju;
	bool ascnd;
	ascnd = ( xx[n-1] >= xx[1] );
	jl=-1;
	ju=n;
	do {
		if (ju-jl <= 1) break;
		jm=(ju+jl)/2;
		if (ascnd && (x >= xx[jm]))
			jl=jm;
		else
			ju=jm;
	} while( true );
	if (x == xx[0])
		{for (i=1;i<n;i++)
			if(xx[i]==x)
				continue;
			else
				break;
		return i-1;}
	else if (x == xx[n-1])
		{for (i=n-2;i>=0;i--)
			if(xx[i]==x)
				continue;
			else
				break;
		return i;}
	else
		return jl;
}

//inrbandec.cxx is revised from inrbandec.f90
//inrbandec.f90 is revised from Numerical Recipes code bandec.f90
//But it is Independent on Numerical Recipes files nrtype.f90 and nrutil.f90
//a(n,m1+m2+1),al(n,m1),indx(n)

void inrbandec(BoutReal **a, int n, int m1, int m2, BoutReal **al, int *indx, BoutReal d)
{	
	const BoutReal TINY=1.0e-20;
	int i,j,k,l,mm;
	BoutReal dum, swaptemp(m1+m2+1);
	mm=m1+m2+1;
	l=m1;
	for (int i=0;i<m1;i++)
	{	//for (int j=m1-i;j<mm;j++) a[i][j-l]=a[i][j];
		for (int j=m1-i;j<mm;j++) (*((BoutReal*)a+i*mm+j-l))=(*((BoutReal*)a+i*mm+j));
		l--;
		//for (int j=mm-l-1;j<mm;j++) a[i][j]=0.0;
		for (int j=mm-l-1;j<mm;j++) (*((BoutReal*)a+i*mm+j))=0.0;
	}
//	cout << "eoshift:\n";
//	for (int i=0;i<n;i++)
//		for (int j=0;j<m1+m2+1;j++)
//			cout << i+1 << ", " << j+1 << ", " << *((BoutReal*)a+i*mm+j) << "\n";

	d=1.0;
	for (i=0;i<n;i++)
		for (j=0;j<m1;j++)
			//al[i][j]=0.0;
			(*((BoutReal*)al+i*m1+j))=0.0;
	for (k=0;k<n;k++)
	{	l=min(m1+k,n-1);
		//dum=a[k][0];
		dum=(*((BoutReal*)a+k*mm));
		i=k;
		for (j=k+1;j<l;j++) {
			//if (abs(a[j][0]) > abs(dum)) {
			if (abs((*((BoutReal*)a+j*mm))) > abs(dum)) {
				//dum=a[j][0];
				dum=(*((BoutReal*)a+j*mm));
				i=j;
			}
		}
		indx[k]=i;
		//cout << k << ", " << indx[k] << ", " << i <<"\n";
		//if (dum == 0.0) a[k][0]=TINY;
		if (dum == 0.0) (*((BoutReal*)a+k*mm))=TINY;
		if (i != k) {
			d = -d;
			//for (j=0;j<mm;j++) swap (a[k][j],a[i][j]);
			for (j=0;j<mm;j++) swap ((*((BoutReal*)a+k*mm+j)),(*((BoutReal*)a+i*mm+j)));
		}
		for (i=k+1;i<=l;i++) {
			//dum=a[i][0]/a[k][0];
			dum=(*((BoutReal*)a+i*mm)) / (*((BoutReal*)a+k*mm));
			//al[k][i-k-1]=dum;
			(*((BoutReal*)al+k*m1+i-k-1))=dum;
			//for (j=1;j<mm;j++) a[i][j-1]=a[i][j]-dum*a[k][j];
			for (j=1;j<mm;j++) (*((BoutReal*)a+i*mm+j-1))=(*((BoutReal*)a+i*mm+j))-dum*(*((BoutReal*)a+k*mm+j));
			//a[i][mm-1]=0.0;
			(*((BoutReal*)a+i*mm+mm-1))=0.0;
		}
	}
}

	//inrbanbks.cxx is revised from inrbanbks.f90
	//inrbanbks.f90 is revised from Numerical Recipes code banbks.f90
	//But it is Independent on Numerical Recipes files nrtype.f90 and nrutil.f90
//a(n,m1+m2+1), al(n,m1), indx(n), b(n)
void inrbanbks(BoutReal **a, int n, int m1, int m2, BoutReal **al, int  *indx, BoutReal *b)
{	int i,k,l,mdum,mm;
	BoutReal swapdum;
	mm = m1+m2+1;
	mdum = m1;
	/*do k=1,n
		l=min(n,m1+k)
		i=indx(k)
		!if (i /= k) call swap(b(i),b(k))
		if (i /=k) then
			swapdum=b(i)
			b(i)=b(k)
			b(k)=swapdum
		end if
		b(k+1:l)=b(k+1:l)-al(k,1:l-k)*b(k)
	end do*/
	for (k=0;k<n;k++) {
		l=min(n-1,m1+k);
		i=indx[k];
		if (i != k )
			{swapdum=b[i];b[i]=b[k];b[k]=swapdum;}
		//for(int j=1;j<=l-k;j++) b[k+j]=b[k+j]-al[k][j-1]*b[k];
		for(int j=1;j<=l-k;j++) b[k+j]=b[k+j]-(*((BoutReal*)al+k*m1+j-1))*b[k];
	}
	/*do i=n,1,-1
		l=min(mm,n-i+1)
		b(i)=(b(i)-dot_product(a(i,2:l),b(1+i:i+l-1)))/a(i,1)
	end do*/
	for (i=n-1;i>=0;i--) {
		l=min(mm,n-i);
		swapdum=0.0;
		for (int j=1;j<l;j++)
			//swapdum = swapdum + a[i][j]*b[i+j];
			swapdum = swapdum + (*((BoutReal*)a+i*mm+j))*b[i+j];
		//b[i]=(b[i]-swapdum)/a[i][0];
		b[i]=(b[i]-swapdum)/(*((BoutReal*)a+i*mm));
	}
}

void cdbbscoef(BoutReal *x, BoutReal *y, int n, BoutReal *knot, int k, BoutReal *bcoef)
{
	void cdbbsvalvec(BoutReal *knot, int k, BoutReal x, int left, BoutReal *val);
	void inrbandec(BoutReal **a, int n, int m1, int m2, BoutReal **al, int *indx, BoutReal d);
	void inrbanbks(BoutReal **a, int n, int m1, int m2, BoutReal **al, int  *indx, BoutReal *b);
	
	int left,right,n0;
	int indx[n];
	BoutReal val[k*(k+1)/2],b,a[n][2*k-1], al[n][k-1];

	left=k-1;
	n0=k*(k-1)/2;
	for (int i=0;i<n;i++)
		for (int j=0;j<2*k-1;j++)
			a[i][j]=0.0;
	for (int i=0;i<n;i++) {
		left=max(left,i);
		right = min ( i + k, n );
		if(x[i] < knot[left])
			cout << "error, x[i] >= knot[i] is not satisfied, i=" << i ;
		
		while ( knot[left+1] <= x[i] ) {
			left = left + 1;
			if ( left < right )
				continue;
			if ( knot[left] < x[i] )
      				cout << "error, x[i] < knot[i+k] is not satisfied, i=" << i;
      			left = left - 1;
			break;
		}
		cdbbsvalvec(knot,k,x[i],left,val);
		for (int j=0;j<k;j++)
			a[i][left-i+j]=val[n0+j];
	}
	inrbandec((BoutReal**)a,n,k-1,k-1,(BoutReal**)al,(int*)indx,b);
	for (int i=0;i<n;i++)
		bcoef[i]=y[i];
	inrbanbks((BoutReal**)a,n,k-1,k-1,(BoutReal**)al,(int*)indx,(BoutReal*)bcoef);
}

//function bvalue ( t, bcoef, n, k, x, jderiv )
//#include "cdblocate.cxx"
//knot[n+k],bcoef[n]
BoutReal cdbbsval(BoutReal *knot, int nk, BoutReal *bcoef, int n, BoutReal x, int jderiv)
{
	int cdblocate(BoutReal *xx, int n, BoutReal x);
	int k,i,ilo,j,jc,jcmax,jcmin,jj;
	k=nk-n;
	BoutReal aj[k], dl[k], dr[k];

	for (i=0;i<k;i++)
		{aj[i]=0.0; dl[i]=0.0; dr[i]=0.0;}

	if ( k <= jderiv )
		return 0.0;

	i=cdblocate(knot, nk, x);

	if ( k <= 1 )
		//cdbbsval = bcoef(i)
		return bcoef[i];

	//jcmin = 1
	jcmin=0;
	//if ( k <= i ) then
	if ( k <= i+1 )
		//do j = 1, k-1
		//	dl(j) = x - knot(i+1-j)
		//end do
		for (j = 0;j< k-1;j++)
			dl[j] = x - knot[i-j];
	else {
		//jcmin = 1 - ( i - k )
		jcmin = k-1-i;
		//do j = 1, i
		//	dl(j) = x - knot(i+1-j)
		//end do
		for (j=0;j<=i;j++)
			dl[j] = x - knot [i-j];
		//do j = i, k-1
		//	aj(k-j) = 0.0_sp
		//	dl(j) = dl(i)
		//end do
		for (j=i;j<k-1;j++)
			{aj[k-j-2]=0.0; dl[j]=dl[i];}
	}

	//jcmax = k
	jcmax = k-1;
	//if ( n < i ) then
	if (n < i+1) {
		//jcmax = k + n - i
		jcmax = k + n - i-2;
		//do j = 1, k + n - i
		//	dr(j) = knot(i+j) - x
		//end do
		for (j=0;j<k+n-i-1;j++)
			dr[j] = knot[i+j+1]-x;
		//do j = k+n-i, k-1
		//	aj(j+1) = 0.0_sp
		//	dr(j) = dr(k+n-i)
		//end do
		for (j=k+n-i-2;j<k-1;j++)
			{aj[j+1]=0.0; dr[j]=dr[k+n-i-2];}
	}
	else
		//do j = 1, k-1
		//	dr(j) = knot(i+j) - x
		//end do
		for (j=0;j<k-1;j++)
			dr[j]= knot[i+j+1]-x;


	//do jc = jcmin, jcmax
	//	aj(jc) = bcoef(i-k+jc)
	//end do
	for (jc=jcmin;jc<=jcmax;jc++)
		aj[jc] = bcoef[i-k+jc+1];

//Difference the coefficients JDERIV times.

	/*do j = 1, jderiv
		ilo = k - j
		do jj = 1, k - j
			aj(jj) = ( ( aj(jj+1) - aj(jj) ) / ( dl(ilo) + dr(jj) ) ) * real ( k - j, kind = sp )
			ilo = ilo - 1
		end do
	end do*/
	for (j=0;j<jderiv;j++) {
		ilo = k-j-2;
		for (jj=0;jj<k-j-1;jj++)
			{aj[jj] = (aj[jj+1]-aj[jj])/(dl[ilo]+dr[jj])*BoutReal(k-j-1); ilo=ilo-1;}
	}
//Compute value at X in (T(I),T(I+1)) of JDERIV-th derivative,
//given its relevant B-spline coefficients in AJ(1),...,AJ(K-JDERIV).

	/*do j = jderiv+1, k-1
		ilo = k-j
		do jj = 1, k-j
			aj(jj) = ( aj(jj+1) * dl(ilo) + aj(jj) * dr(jj) ) / ( dl(ilo) + dr(jj) )
			ilo = ilo - 1
		end do
	end do*/
	for (j = jderiv; j< k-1;j++) {
		ilo=k-j-2;
		for (jj=0;jj<k-j-1;jj++)
			{aj[jj]=(aj[jj+1]*dl[ilo]+aj[jj]*dr[jj])/(dl[ilo]+dr[jj]); ilo=ilo-1;}
	}

	//cdbbsval = aj(1)
	return aj[0];
}

//输入x，y，函数值fxy，以及节点序列knotx，knoty，输出系数coef2d。
//x(nx),y(ny),fxy(nx,ny),knotx(nx+kx),knoty(ny+ky),coef2d(nx,ny)
//#include "cdbbsvalvec.cxx"
//#include "inrbandec.cxx"
//#include "inrbanbks.cxx"
void cdbbscoef2d(BoutReal *x, BoutReal *y, Field2D &fxy, int nx, int ny, BoutReal *knotx, BoutReal *knoty, int kx, int ky, BoutReal **coef2d) {
	void cdbbsvalvec(BoutReal *knot, int k, BoutReal x, int left, BoutReal *val);
	void inrbandec(BoutReal **a, int n, int m1, int m2, BoutReal **al, int *indx, BoutReal d);
	void inrbanbks(BoutReal **a, int n, int m1, int m2, BoutReal **al, int  *indx, BoutReal *b);
	
	int i,j;
	BoutReal coef[nx],coefy[ny],temp[nx][ny];
	int left,n0,right,indxx[nx],indxy[ny];
	BoutReal b,alx[nx][kx-1],aly[ny][ky-1],valx[kx*(kx+1)/2],valy[ky*(ky+1)/2],ax[nx][2*kx-1],ay[ny][2*ky-1];	

	left=kx-1;
	n0=kx*(kx-1)/2;
	for (int i=0;i<nx;i++)
		for (int j=0;j<2*kx-1;j++)
			ax[i][j]=0.0;
	for (int i=0;i<nx;i++) {
		left=max(left,i);
		right = min ( i + kx, nx );
		if(x[i] < knotx[left])
			cout << "error, x[i] >= knotx[i] is not satisfied, i=" << i ;
		
		while ( knotx[left+1] <= x[i] ) {
			left = left + 1;
			if ( left < right )
				continue;
			if ( knotx[left] < x[i] )
      				cout << "error, x[i] < knotx[i+kx] is not satisfied, i=" << i;
      			left = left - 1;
			break;
		}
		cdbbsvalvec(knotx,kx,x[i],left,valx);
		for (int j=0;j<kx;j++)
			ax[i][left-i+j]=valx[n0+j];
	}
	inrbandec((BoutReal**)ax,nx,kx-1,kx-1,(BoutReal**)alx,(int*)indxx,b);
	for ( i=0;i<ny;i++) {
		for ( j=0;j<nx;j++)
			{coef[j]=fxy[j][i];}
			//{coef[j]=(*((BoutReal*)fxy+j*ny+i));}
		inrbanbks((BoutReal**)ax,nx,kx-1,kx-1,(BoutReal**)alx,(int*)indxx,(BoutReal*)coef);
		for (j=0;j<nx;j++)
			{temp[j][i]=coef[j];}
	}
	left=ky-1;
	n0=ky*(ky-1)/2;
	for (int i=0;i<ny;i++)
		for (int j=0;j<2*ky-1;j++)
			ay[i][j]=0.0;
	for (int i=0;i<ny;i++) {
		left=max(left,i);
		right = min ( i + ky, ny );
		if(y[i] < knoty[left])
			cout << "error, y[i] >= knoty[i] is not satisfied, i=" << i ;
		
		while ( knoty[left+1] <= y[i] ) {
			left = left + 1;
			if ( left < right )
				continue;
			if ( knoty[left] < y[i] )
      				cout << "error, y[i] < knoty[i+k] is not satisfied, i=" << i;
      			left = left - 1;
			break;
		}
		cdbbsvalvec(knoty,ky,y[i],left,valy);
		for (int j=0;j<ky;j++)
			ay[i][left-i+j]=valy[n0+j];
	}
	inrbandec((BoutReal**)ay,ny,ky-1,ky-1,(BoutReal**)aly,(int*)indxy,b);
	for (i=0;i<nx;i++) {
		for (j=0;j<ny;j++)
			{coefy[j]=temp[i][j];}
		inrbanbks((BoutReal**)ay,ny,ky-1,ky-1,(BoutReal**)aly,(int*)indxy,(BoutReal*)coefy);
		for (j=0;j<ny;j++)
			//{coef2d[i][j]=coefy[j];}
			{(*((BoutReal*)coef2d+i*ny+j))=coefy[j];}
	}
}	

void cdbbscoef2d(BoutReal *x, BoutReal *y, BoutReal **fxy, int nx, int ny, BoutReal *knotx, BoutReal *knoty, int kx, int ky, BoutReal **coef2d) {
	void cdbbsvalvec(BoutReal *knot, int k, BoutReal x, int left, BoutReal *val);
	void inrbandec(BoutReal **a, int n, int m1, int m2, BoutReal **al, int *indx, BoutReal d);
	void inrbanbks(BoutReal **a, int n, int m1, int m2, BoutReal **al, int  *indx, BoutReal *b);
	
	int i,j;
	BoutReal coef[nx],coefy[ny],temp[nx][ny];
	int left,n0,right,indxx[nx],indxy[ny];
	BoutReal b,alx[nx][kx-1],aly[ny][ky-1],valx[kx*(kx+1)/2],valy[ky*(ky+1)/2],ax[nx][2*kx-1],ay[ny][2*ky-1];	

	left=kx-1;
	n0=kx*(kx-1)/2;
	for (int i=0;i<nx;i++)
		for (int j=0;j<2*kx-1;j++)
			ax[i][j]=0.0;
	for (int i=0;i<nx;i++) {
		left=max(left,i);
		right = min ( i + kx, nx );
		if(x[i] < knotx[left])
			cout << "error, x[i] >= knotx[i] is not satisfied, i=" << i ;
		
		while ( knotx[left+1] <= x[i] ) {
			left = left + 1;
			if ( left < right )
				continue;
			if ( knotx[left] < x[i] )
      				cout << "error, x[i] < knotx[i+kx] is not satisfied, i=" << i;
      			left = left - 1;
			break;
		}
		cdbbsvalvec(knotx,kx,x[i],left,valx);
		for (int j=0;j<kx;j++)
			ax[i][left-i+j]=valx[n0+j];
	}
	inrbandec((BoutReal**)ax,nx,kx-1,kx-1,(BoutReal**)alx,(int*)indxx,b);
	for ( i=0;i<ny;i++) {
		for ( j=0;j<nx;j++)
			//{coef[j]=fxy[j][i];}
			{coef[j]=(*((BoutReal*)fxy+j*ny+i));}
		inrbanbks((BoutReal**)ax,nx,kx-1,kx-1,(BoutReal**)alx,(int*)indxx,(BoutReal*)coef);
		for (j=0;j<nx;j++)
			{temp[j][i]=coef[j];}
	}
	left=ky-1;
	n0=ky*(ky-1)/2;
	for (int i=0;i<ny;i++)
		for (int j=0;j<2*ky-1;j++)
			ay[i][j]=0.0;
	for (int i=0;i<ny;i++) {
		left=max(left,i);
		right = min ( i + ky, ny );
		if(y[i] < knoty[left])
			cout << "error, y[i] >= knoty[i] is not satisfied, i=" << i ;
		
		while ( knoty[left+1] <= y[i] ) {
			left = left + 1;
			if ( left < right )
				continue;
			if ( knoty[left] < y[i] )
      				cout << "error, y[i] < knoty[i+k] is not satisfied, i=" << i;
      			left = left - 1;
			break;
		}
		cdbbsvalvec(knoty,ky,y[i],left,valy);
		for (int j=0;j<ky;j++)
			ay[i][left-i+j]=valy[n0+j];
	}
	inrbandec((BoutReal**)ay,ny,ky-1,ky-1,(BoutReal**)aly,(int*)indxy,b);
	for (i=0;i<nx;i++) {
		for (j=0;j<ny;j++)
			{coefy[j]=temp[i][j];}
		inrbanbks((BoutReal**)ay,ny,ky-1,ky-1,(BoutReal**)aly,(int*)indxy,(BoutReal*)coefy);
		for (j=0;j<ny;j++)
			//{coef2d[i][j]=coefy[j];}
			{(*((BoutReal*)coef2d+i*ny+j))=coefy[j];}
	}
}	


//knotx(nx+kx),knoty(ny+ky),coef2d(nx,ny),x(nx),y(ny)
//#include "cdbbsval.cxx"
BoutReal cdbbsval2d(BoutReal *knotx, int nkx, BoutReal *knoty, int nky, BoutReal **coef2d, int nx, int ny, BoutReal x, BoutReal y, int derivx, int derivy) {
	int cdblocate(BoutReal *xx, int n, BoutReal x);
	BoutReal cdbbsval(BoutReal *knot, int nk, BoutReal *bcoef, int n, BoutReal x, int jderiv);
	
	int i,j,lefty,kx,ky;
	kx=nkx-nx;
	ky=nky-ny;
	BoutReal temp[ky],temp1[nx],temp2[2*ky-1];
	
	lefty=cdblocate(knoty, ny+ky, y);
	//do i=1,ky
	//	temp(i)=cdbbsval(knotx,kx,coef2d(1:nx,lefty-ky+i),x,derivx)
	//end do
	for (i=0;i<ky;i++) {
		for (j=0;j<nx;j++)
			//{temp1[j]=coef2d[j][lefty-ky+i+1];}
			{temp1[j]=(*((BoutReal*)coef2d+j*ny+lefty-ky+i+1));}
		temp[i] = cdbbsval(knotx, nx+kx, temp1, nx, x, derivx);
	}
	//cdbbsval2d = cdbbsval(knoty(lefty-ky+1:lefty+ky),ky,temp,y,derivy)
	for (i=0;i<2*ky;i++)
		{temp2[i] = knoty[lefty-ky+i+1];}
	return cdbbsval(temp2, 2*ky, temp, ky, y, derivy);
}

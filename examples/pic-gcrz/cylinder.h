class scal1D;
class scalRZ;
class vecRZ;

class cylinder {
	friend class vecRZ;
	public:
		BoutReal r,p,z;
		BoutReal dot(cylinder vec);
		cylinder cross(cylinder vec);
};

class mesh1D {
	friend class scal1D;
	protected:
		int n,k;
		BoutReal *grid,*knot;
	public:
		mesh1D();
		mesh1D(int,BoutReal*,int);
		void init(int,BoutReal*,int);
};

class scal1D {
	private:
		BoutReal *bscoef;
	public:
		scal1D();
		scal1D(mesh1D);
		void init(mesh1D);
		BoutReal *f;
		void coef(mesh1D);
		BoutReal val(BoutReal,int,mesh1D);
};

class meshRZ {
	friend class scalRZ;
	friend class vecRZ;
	protected:
		
	public:
		int nr,nz,kr,kz;
		meshRZ();
		BoutReal *gridr,*gridz,*knotr,*knotz;
		meshRZ(int, int, int, int);
		meshRZ(int, BoutReal*, int, BoutReal*, int, int);
		void init(int, int, int, int);
		void init(int, BoutReal*, int, BoutReal*, int, int);	
};

class scalRZ {
	friend class vecRZ;
	private:
		
	public:
		scalRZ();
		//~scalRZ();
		scalRZ(meshRZ);
		void init(meshRZ);
		BoutReal **f;
		void coef(meshRZ);
		BoutReal val(BoutReal,BoutReal,int,int,meshRZ);
		void clearf();
		void clearcoef();
		BoutReal **bscoef;	
};

class vecRZ {
	public:
		scalRZ vr,vp,vz;
		vecRZ();
		//~vecRZ();
		vecRZ(meshRZ);
		void init(meshRZ);
		void coef(meshRZ);
		cylinder val(BoutReal,BoutReal,int,int,meshRZ);
		void clearf();
		void clearcoef();
};

class particle {
	public:
		BoutReal y[4],mu,psi,theta;
		int region;
		//void locate();
		//void map();
};

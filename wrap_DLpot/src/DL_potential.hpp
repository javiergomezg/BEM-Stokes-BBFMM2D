//#ifndef __DL_potential_hpp__
//#define __DL_potential_hpp__

#include"header/BBFMM2D.hpp"
#define REAL double

using namespace std;
using namespace Eigen;

class myKernelTxxx: public kernel_Base {
public:
    virtual double kernel_Func(Point r0, Point r1){
        //implement your own kernel here
	double dx = (r1.x - r0.x);
	double dy = (r1.y - r0.y);
        double rSquare	= dx*dx + dy*dy;
	double r = pow(rSquare,0.5);
	if(r > 10e-16){
		return -4.*dx*dx*dx/pow(r,4.);
	}
    }
};

class myKernelTxyx: public kernel_Base { //sirve con Txxy y Tyxx
public:
    virtual double kernel_Func(Point r0, Point r1){
        //implement your own kernel here
	double dx = (r1.x - r0.x);
	double dy = (r1.y - r0.y);
        double rSquare	= dx*dx + dy*dy;
	double r = pow(rSquare,0.5);
	if(r > 10e-16){
	//if(r0.panel == r1.panel){
		return -4.*dx*dy*dx/pow(r,4.);
	}
    }
};

class myKernelTyxy: public kernel_Base { //sirve con Txyy y Tyyx
public:
    virtual double kernel_Func(Point r0, Point r1){
        //implement your own kernel here
	double dx = (r1.x - r0.x);
	double dy = (r1.y - r0.y);
        double rSquare	= dx*dx + dy*dy;
	double r = pow(rSquare,0.5);
	if(r > 10e-16){
	//if(r0.panel == r1.panel){
		return -4.*dy*dx*dy/pow(r,4.);
	}
    }
};

class myKernelTyyy: public kernel_Base {
public:
    virtual double kernel_Func(Point r0, Point r1){
        //implement your own kernel here
	double dx = (r1.x - r0.x);
	double dy = (r1.y - r0.y);
        double rSquare	= dx*dx + dy*dy;
	double r = pow(rSquare,0.5);
	if(r > 10e-16){
	//if(r0.panel == r1.panel){
		return -4.*dy*dy*dy/pow(r,4.);
	}
    }
};


void bbfmm_DL_potential_cpp(REAL *xsrc, int xsrcSize,
			REAL *ysrc, int ysrcSize,
			REAL *xtar, int xtarSize,
			REAL *ytar, int ytarSize,
			REAL *Rp, int n_pan,
			REAL *wc, int n_wc,
			REAL *ux, int n_ux,
			REAL *uy, int n_uy,
			REAL *nx, int n_nx,
			REAL *ny, int n_ny,
			REAL *bfmm, int n_bfmm,
			int nchebappol);

void bbfmm_DL_potential_cpp(REAL *xsrc, int xsrcSize,
			REAL *ysrc, int ysrcSize,
			REAL *xtar, int xtarSize,
			REAL *ytar, int ytarSize,
			REAL *Rp, int n_pan,
			REAL *wc, int n_wc,
			REAL *ux, int n_ux,
			REAL *uy, int n_uy,
			REAL *nx, int n_nx,
			REAL *ny, int n_ny,
			REAL *bfmm, int n_bfmm,
			int nchebappol){

	int nq = xsrcSize/n_pan;
	unsigned long N = xsrcSize + xtarSize;
	unsigned long n_src = xsrcSize;
	unsigned m = 1;
	vector<Point> location;
	unsigned short nChebNodes = nchebappol;

	//cout << "Digite los nodos de chebyshev para aproximar el kernel " <<endl;
	//cin >> nChebNodes;

	REAL pi = acos(-1.), ddpi = 4.*pi;

	for(int i=0; i < N; i++){
		if(i < xsrcSize){
			location.push_back(Point(i/nq, xsrc[i],ysrc[i],false));
			//cout<<i<<"; panel: "<< location[i].panel<<endl;
			
		}
		else{
			location.push_back(Point(i-xsrcSize, xtar[i-xsrcSize],ytar[i-xsrcSize],true));
			//cout<<i<<"  "<<i-xsrcSize <<"; panel: "<< location[i].panel<<endl;
		}
	}

	REAL *beta11;
	beta11 = new double[N*m];

	REAL *beta12;
	beta12 = new double[N*m];

	REAL *beta21;
	beta21 = new double[N*m];

	REAL *beta22;
	beta22 = new double[N*m];

	
	for(int i = 0; i < xsrcSize; i++){
		beta11[i] = (-1./ddpi)*Rp[i/nq]*wc[i]*ux[i/nq]*nx[i/nq];
		//cout << "b11_"<<i<<" = "<<beta11[i]<< endl;
	}

	for(int i= xsrcSize; i < N; i++){
		beta11[i] = 0.0;
		//cout << "b11_"<<i<<" = "<<beta11[i]<< endl;
	}

	
	for(int i = 0; i < xsrcSize; i++){
		beta12[i] = (-1./ddpi)*Rp[i/nq]*wc[i]*ux[i/nq]*ny[i/nq];
		//cout << "b12_"<<i<<" = "<<beta12[i]<< endl;
	}

	for(int i= xsrcSize; i < N; i++){
		beta12[i] = 0.0;
		//cout << "b12_"<<i<<" = "<<beta12[i]<< endl;
	}

	
	for(int i = 0; i < xsrcSize; i++){
		beta21[i] = (-1./ddpi)*Rp[i/nq]*wc[i]*uy[i/nq]*nx[i/nq];
		//cout << "b21_"<<i<<" = "<<beta21[i]<< endl;
	}

	for(int i= xsrcSize; i < N; i++){
		beta21[i] = 0.0;
		//cout << "b21_"<<i<<" = "<<beta21[i]<< endl;
	}

	for(int i = 0; i < xsrcSize; i++){
		beta22[i] = (-1./ddpi)*Rp[i/nq]*wc[i]*uy[i/nq]*ny[i/nq];
		//cout << "b22_"<<i<<" = "<<beta22[i]<< endl;
	}

	for(int i= xsrcSize; i < N; i++){
		beta22[i] = 0.0;
		//cout << "b22_"<<i<<" = "<<beta22[i]<< endl;
	}
	

	H2_2D_Tree Atree(nChebNodes, beta11, location, N, m);// Build the fmm tree;
	
	REAL *potentialbxxx; REAL *potentialbxyx;
	potentialbxxx = new double[N*m]; potentialbxyx = new double[N*m];

	myKernelTxxx bxxx;
    	bxxx.calculate_Potential(Atree, potentialbxxx);

	myKernelTxyx bxyx;
    	bxyx.calculate_Potential(Atree, potentialbxyx);

	//Atree.~H2_2D_Tree(); 	

	H2_2D_Tree Btree(nChebNodes, beta12, location, N, m);// Build the fmm tree;
	
	REAL *potentialbxxy; REAL *potentialbxyy;	
	potentialbxxy = new double[N*m]; potentialbxyy = new double[N*m];

	myKernelTxyx bxxy;
    	bxxy.calculate_Potential(Btree,potentialbxxy);

	myKernelTyxy bxyy;
    	bxyy.calculate_Potential(Btree,potentialbxyy);

	//Btree.~H2_2D_Tree();

	H2_2D_Tree Ctree(nChebNodes, beta21, location, N, m);// Build the fmm tree;
	
	REAL *potentialbyxx; REAL *potentialbyyx;
	potentialbyxx = new double[N*m]; potentialbyyx = new double[N*m];	

	myKernelTxyx byxx;
    	byxx.calculate_Potential(Ctree,potentialbyxx);

	myKernelTyxy byyx;
    	byyx.calculate_Potential(Ctree,potentialbyyx);

	//Ctree.~H2_2D_Tree();

	H2_2D_Tree Dtree(nChebNodes, beta22, location, N, m);// Build the fmm tree;
	
	REAL *potentialbyxy; REAL *potentialbyyy;	
	potentialbyxy = new double[N*m]; potentialbyyy = new double[N*m];	

	myKernelTyxy byxy;
    	byxy.calculate_Potential(Dtree,potentialbyxy);

	myKernelTyyy byyy;
    	byyy.calculate_Potential(Dtree,potentialbyyy);

	//Dtree.~H2_2D_Tree();

	REAL bxfmm[n_pan], byfmm[n_pan];
	
	for(int i = xsrcSize; i < N; i++){
		//bxfmm[i-xsrcSize] = potentialbxxx[i];
		bxfmm[i-xsrcSize] = ux[i-xsrcSize] + potentialbxxx[i] + potentialbxxy[i] + potentialbyxx[i] + potentialbyxy[i];
		bfmm[i-xsrcSize] = bxfmm[i-xsrcSize];
		//cout << "bxfmm_" << i - n_src << " = " << bxfmm[i-n_src] << endl;
	}

	cout <<endl;


	for(int i = xsrcSize; i < N; i++){
		byfmm[i-xsrcSize] = uy[i-xsrcSize] + potentialbxyx[i] + potentialbxyy[i] + potentialbyyx[i] + potentialbyyy[i];
		bfmm[i-xsrcSize + n_pan] = byfmm[i-xsrcSize];
		//cout << "byfmm_" << i << " = " << byfmm[i-n_src] << endl;
	}

	delete[] beta11, beta12, beta21, beta22;
	delete[] potentialbxxx, potentialbxxy, potentialbyxx, potentialbyxy;
	delete[] potentialbxyx, potentialbxyy, potentialbyyx, potentialbyyy; 

}


//#endif //__DL_potential_hpp__

#ifndef __ufield_fmm_hpp__
#define __ufielf_fmm_hpp__

#include"header/BBFMM2D.hpp"
#define REAL double

using namespace std;
using namespace Eigen;

// ##############Single layer potential kernel################

class myKernelPx: public kernel_Base {
public:
    virtual double kernel_Func(Point r0, Point r1){
        //implement your own kernel here
	double dx = -(r1.x - r0.x);
	double dy = -(r1.y - r0.y);
        double rSquare	= dx*dx + dy*dy;
	double r = pow(rSquare,0.5);
	if(r > 10e-16){
	//if(r0.panel != r1.panel){
	        return dx/rSquare;
	}
    }
};

class myKernelPy: public kernel_Base { //sirve con Gyx
public:
    virtual double kernel_Func(Point r0, Point r1){
        //implement your own kernel here
	double dx = -(r1.x - r0.x);
	double dy = -(r1.y - r0.y);
        double rSquare	= dx*dx + dy*dy;
	double r = pow(rSquare,0.5);
	if(r > 10e-16){
	//if(r0.panel != r1.panel){
		return dy/rSquare;
	}
    }
};

// ##############Single layer potential kernel################

// ##############Double layer potential kernel################

class myKernelPIxx: public kernel_Base {
public:
    virtual double kernel_Func(Point r0, Point r1){
        //implement your own kernel here
	double dx = (r1.x - r0.x);
	double dy = (r1.y - r0.y);
        double rSquare	= dx*dx + dy*dy;
	double r = pow(rSquare,0.5);
	if(r > 10e-16){
		return 4.*(-(1./rSquare) + 2.*dx*dx/pow(r,4.));
	}
    }
};

class myKernelPIxy: public kernel_Base { //sirve con PIyx
public:
    virtual double kernel_Func(Point r0, Point r1){
        //implement your own kernel here
	double dx = (r1.x - r0.x);
	double dy = (r1.y - r0.y);
        double rSquare	= dx*dx + dy*dy;
	double r = pow(rSquare,0.5);
	if(r > 10e-16){
	//if(r0.panel == r1.panel){
		return 4.*(2.*dx*dy/pow(r,4.));
	}
    }
};

class myKernelPIyy: public kernel_Base { 
public:
    virtual double kernel_Func(Point r0, Point r1){
        //implement your own kernel here
	double dx = (r1.x - r0.x);
	double dy = (r1.y - r0.y);
        double rSquare	= dx*dx + dy*dy;
	double r = pow(rSquare,0.5);
	if(r > 10e-16){
	//if(r0.panel == r1.panel){
		return 4.*(-(1./rSquare) + 2.*dy*dy/pow(r,4.));
	}
    }
};

// ##############Double layer potential kernel################

void fmm_pfield_cpp(REAL *xsrc, int xsrcSize,
			REAL *ysrc, int ysrcSize,
			REAL *xtar, int xtarSize,
			REAL *ytar, int ytarSize,
			REAL *Rp, int n_pan,
			REAL *wc, int n_wc,
			REAL *ux, int n_ux,
			REAL *uy, int n_uy,
			REAL *nx, int n_nx,
			REAL *ny, int n_ny,
			REAL *fx, int n_fx,
			REAL *fy, int n_fy,
			REAL Re,
			REAL *P_f, int n_P_f,
			int nchebappol);

void fmm_pfield_cpp(REAL *xsrc, int xsrcSize,
			REAL *ysrc, int ysrcSize,
			REAL *xtar, int xtarSize,
			REAL *ytar, int ytarSize,
			REAL *Rp, int n_pan,
			REAL *wc, int n_wc,
			REAL *ux, int n_ux,
			REAL *uy, int n_uy,
			REAL *nx, int n_nx,
			REAL *ny, int n_ny,
			REAL *fx, int n_fx,
			REAL *fy, int n_fy,
			REAL Re,
			REAL *P_f, int n_P_f,
			int nchebappol){

	int nq = xsrcSize/n_pan;
	unsigned long N = xsrcSize + xtarSize;
	unsigned long n_src = xsrcSize;
	unsigned m = 1;
	vector<Point> location;
	unsigned short nChebNodes = nchebappol;

	//cout << "Digite los nodos de chebyshev para aproximar el kernel " <<endl;
	//cin >> nChebNodes;

	double pi = acos(-1.), dpi = 2.*pi, ddpi = 4.*pi, octpi = 8.*pi;

	for(int i=0; i < N; i++){
		if(i < xsrcSize){
			location.push_back(Point(i/nq, xsrc[i],ysrc[i], false));
			//cout<<i<<"; panel: "<< location[i].panel<<endl;
			
		}
		else{
			location.push_back(Point(i-xsrcSize, xtar[i-xsrcSize],ytar[i-xsrcSize], true));
			//cout<<i<<"  "<<i-xsrcSize <<"; panel: "<< location[i].panel<<endl;
		}
	}

	//### SLP_precalc
	double* gamma1;
	gamma1 = new double[N*m];

	double* gamma2;
	gamma2 = new double[N*m];


	for(unsigned long i = 0; i < n_src; i++){
		gamma1[i] = ((-1./octpi)*wc[i]*fx[i/nq]*Rp[i/nq]);
	}

	for(unsigned long i = 0; i < n_src; i++){
		gamma2[i] = ((-1./octpi)*wc[i]*fy[i/nq]*Rp[i/nq]);
	}

	for(unsigned long i = n_src; i < N; i++){
		gamma1[i] = 0.0;
	}

	for(unsigned long i = n_src; i < N; i++){
		gamma2[i] = 0.0;
	}

	H2_2D_Tree Atree(nChebNodes, gamma1, location, N, m);
	
	double* potentialPx;
	potentialPx = new double[N*m]; 

	myKernelPx Px;
    	Px.calculate_Potential(Atree, potentialPx);


	H2_2D_Tree Btree(nChebNodes, gamma2, location, N, m);
	
	double* potentialPy; 
	potentialPy = new double[N*m];

	myKernelPy Py;
    	Py.calculate_Potential(Btree, potentialPy);

	//Atree.~H2_2D_Tree(); Btree.~H2_2D_Tree();
	delete[] gamma1, gamma2;

	//### SLP_precalc

	//### DLP_precalc

	REAL *beta11;
	beta11 = new double[N*m];

	REAL *beta12;
	beta12 = new double[N*m];

	REAL *beta21;
	beta21 = new double[N*m];

	REAL *beta22;
	beta22 = new double[N*m];

	
	for(int i = 0; i < xsrcSize; i++){
		beta11[i] = (1./(octpi*Re))*Rp[i/nq]*wc[i]*ux[i/nq]*nx[i/nq];
		//cout << "b11_"<<i<<" = "<<beta11[i]<< endl;
	}

	for(int i= xsrcSize; i < N; i++){
		beta11[i] = 0.0;
		//cout << "b11_"<<i<<" = "<<beta11[i]<< endl;
	}

	
	for(int i = 0; i < xsrcSize; i++){
		beta12[i] = (1./(octpi*Re))*Rp[i/nq]*wc[i]*ux[i/nq]*ny[i/nq];
		//cout << "b12_"<<i<<" = "<<beta12[i]<< endl;
	}

	for(int i= xsrcSize; i < N; i++){
		beta12[i] = 0.0;
		//cout << "b12_"<<i<<" = "<<beta12[i]<< endl;
	}

	
	for(int i = 0; i < xsrcSize; i++){
		beta21[i] = (1./(octpi*Re))*Rp[i/nq]*wc[i]*uy[i/nq]*nx[i/nq];
		//cout << "b21_"<<i<<" = "<<beta21[i]<< endl;
	}

	for(int i= xsrcSize; i < N; i++){
		beta21[i] = 0.0;
		//cout << "b21_"<<i<<" = "<<beta21[i]<< endl;
	}

	for(int i = 0; i < xsrcSize; i++){
		beta22[i] = (1./(octpi*Re))*Rp[i/nq]*wc[i]*uy[i/nq]*ny[i/nq];
		//cout << "b22_"<<i<<" = "<<beta22[i]<< endl;
	}

	for(int i= xsrcSize; i < N; i++){
		beta22[i] = 0.0;
		//cout << "b22_"<<i<<" = "<<beta22[i]<< endl;
	}	

	H2_2D_Tree tree1(nChebNodes, beta11, location, N, m);// Build the fmm tree;
	
	REAL *potentialPIxx;
	potentialPIxx = new double[N*m];

	myKernelPIxx PIxx;
    	PIxx.calculate_Potential(tree1, potentialPIxx);

	
	H2_2D_Tree tree2(nChebNodes, beta12, location, N, m);// Build the fmm tree;
	
	REAL *potentialPIxy; 	
	potentialPIxy = new double[N*m];

	myKernelPIxy PIxy;
    	PIxy.calculate_Potential(tree2,potentialPIxy);


	H2_2D_Tree tree3(nChebNodes, beta21, location, N, m);// Build the fmm tree;
	
	REAL *potentialPIyx;
	potentialPIyx = new double[N*m];

	myKernelPIxy PIyx;
    	PIyx.calculate_Potential(tree3,potentialPIyx);


	H2_2D_Tree tree4(nChebNodes, beta22, location, N, m);// Build the fmm tree;
	
	REAL *potentialPIyy;	
	potentialPIyy = new double[N*m]; 

	myKernelPIyy PIyy;
    	PIyy.calculate_Potential(tree4,potentialPIyy);

	//tree1.~H2_2D_Tree(); tree2.~H2_2D_Tree(); tree3.~H2_2D_Tree(); tree4.~H2_2D_Tree();
	delete[] beta11, beta12, beta21, beta22;

	for(int i= xsrcSize; i<N; i++){
		P_f[i-xsrcSize] += potentialPx[i] + potentialPy[i] + potentialPIxx[i] + potentialPIxy[i] + potentialPIyx[i] + potentialPIyy[i]; 
	}

	
	delete[] potentialPx, potentialPy; 
	delete[] potentialPIxx, potentialPIxy, potentialPIyx, potentialPIyy;
	
}


#endif //__ufield_fmm_hpp__

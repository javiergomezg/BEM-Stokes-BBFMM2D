 #ifndef __SL_potential_hpp__
#define __SL_potential_hpp__

#include"header/BBFMM2D.hpp"
#define REAL double

using namespace std;
using namespace Eigen;

class myKernelGxx: public kernel_Base {
public:
    virtual double kernel_Func(Point r0, Point r1){
        //implement your own kernel here
	double dx = -(r1.x - r0.x);
	double dy = -(r1.y - r0.y);
        double rSquare	= dx*dx + dy*dy;
	double r = pow(rSquare,0.5);
	if(r > 10e-16){
	//if(r0.panel != r1.panel){
	        return -log(r) + dx*dx/rSquare;
	}
    }
};

class myKernelGxy: public kernel_Base { //sirve con Gyx
public:
    virtual double kernel_Func(Point r0, Point r1){
        //implement your own kernel here
	double dx = -(r1.x - r0.x);
	double dy = -(r1.y - r0.y);
        double rSquare	= dx*dx + dy*dy;
	double r = pow(rSquare,0.5);
	if(r > 10e-16){
	//if(r0.panel != r1.panel){
		return dx*dy/rSquare;
	}
    }
};

class myKernelGyy: public kernel_Base { 
public:
    virtual double kernel_Func(Point r0, Point r1){
        //implement your own kernel here
	double dx = -(r1.x - r0.x);
	double dy = -(r1.y - r0.y);
        double rSquare	= dx*dx + dy*dy;
	double r = pow(rSquare,0.5);
	if(r > 10e-16){
	//if(r0.panel != r1.panel){
		return -log(r) + dy*dy/rSquare;
	}
    }
};


void bbfmm_SL_potential_cpp(REAL* xsrc, int xsrcsize,
			REAL *ysrc, int ysrcSize,
			REAL *xtar, int xtarSize,
			REAL *ytar, int ytarSize,
			REAL *Rp, int n_pan,
			REAL *wc, int n_wc,
			REAL *fx, int n_fx,
			REAL *fy, int n_fy,
			REAL Re, 
			REAL *Axfmm, int n_Axfmm,
			int nchebappol);

void bbfmm_SL_potential_cpp(REAL *xsrc, int xsrcSize,
			REAL *ysrc, int ysrcSize,
			REAL *xtar, int xtarSize,
			REAL *ytar, int ytarSize,
			REAL *Rp, int n_pan,
			REAL *wc, int n_wc,
			REAL *fx, int n_fx,
			REAL *fy, int n_fy,
			REAL Re,
			REAL *Axfmm, int n_Axfmm,
			int nchebappol){

	int nq = xsrcSize/n_pan;
	unsigned long N = xsrcSize + xtarSize;
	unsigned long n_src = xsrcSize;
	unsigned m = 1;
	vector<Point> location;
	unsigned short nChebNodes = nchebappol;

	//cout << "Digite los nodos de chebyshev para aproximar el kernel " <<endl;
	//cin >> nChebNodes;

	double pi = acos(-1.), dpi = 2.*pi, ddpi = 4.*pi;

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

	double* gamma1;
	gamma1 = new double[N*m];

	double* gamma2;
	gamma2 = new double[N*m];


	for(unsigned long i = 0; i < n_src; i++){
		gamma1[i] = ((-Re/ddpi)*wc[i]*fx[i/nq]*Rp[i/nq]);
	}

	for(unsigned long i = 0; i < n_src; i++){
		gamma2[i] = ((-Re/ddpi)*wc[i]*fy[i/nq]*Rp[i/nq]);
	}

	for(unsigned long i = n_src; i < N; i++){
		gamma1[i] = 0.0;
	}

	for(unsigned long i = n_src; i < N; i++){
		gamma2[i] = 0.0;
	}

	H2_2D_Tree Atree(nChebNodes, gamma1, location, N, m);
	
	double* potentialgxx; double* potentialgyx;
	potentialgxx = new double[N*m]; potentialgyx = new double[N*m];

	myKernelGxx gxx;
    	gxx.calculate_Potential(Atree, potentialgxx);

	myKernelGxy gyx;
    	gyx.calculate_Potential(Atree, potentialgyx);


	H2_2D_Tree Btree(nChebNodes, gamma2, location, N, m);
	
	double* potentialgyy; double* potentialgxy;
	potentialgyy = new double[N*m]; potentialgxy = new double[N*m];

	myKernelGyy gyy;
    	gyy.calculate_Potential(Btree, potentialgyy);

	myKernelGxy gxy;
    	gxy.calculate_Potential(Btree, potentialgxy);


	double Axxfmm[n_pan], Axyfmm[n_pan];
	
	for(int i = xsrcSize; i < N; i++){
		Axxfmm[i-xsrcSize] = potentialgxx[i] + potentialgxy[i];
		Axfmm[i-xsrcSize] = -Axxfmm[i-xsrcSize];
	}

	for(int i = xsrcSize; i < N; i++){
		Axyfmm[i-xsrcSize] = potentialgyx[i] + potentialgyy[i];
		Axfmm[i-xsrcSize + n_pan] = -Axyfmm[i-xsrcSize];		
	}

	//Atree.~H2_2D_Tree(); Btree.~H2_2D_Tree();
	
	delete[] gamma1, gamma2;
	delete[] potentialgxx, potentialgxy, potentialgyx, potentialgyy;

}


#endif //__DL_potential_hpp__

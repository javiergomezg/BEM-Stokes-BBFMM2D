//#ifndef __DL_potential_dir_hpp__
//#define __DL_potential_dir_hpp__

#include"header/BBFMM2D.hpp"
#define REAL double

using namespace std;
using namespace Eigen;

void direct_DL_potential_cpp(REAL *xsrc, int xsrcsize,
			REAL *ysrc, int ysrcSize,
			REAL *xtar, int xtarSize,
			REAL *ytar, int ytarSize,
			REAL *Rp, int n_pan,
			REAL *wc, int n_wc,
			REAL *ux, int n_ux,
			REAL *uy, int n_uy,
			REAL *nx, int n_nx,
			REAL *ny, int n_ny,
			REAL *b, int n_b);


void direct_DL_potential_cpp(REAL *xsrc, int xsrcSize,
			REAL *ysrc, int ysrcSize,
			REAL *xtar, int xtarSize,
			REAL *ytar, int ytarSize,
			REAL *Rp, int n_pan,
			REAL *wc, int n_wc,
			REAL *ux, int n_ux,
			REAL *uy, int n_uy,
			REAL *nx, int n_nx,
			REAL *ny, int n_ny,
			REAL *b, int n_b){

	clock_t startBuild	=	clock();

	REAL Txxx; REAL Txxy; REAL Tyxx; REAL Tyxy;
	REAL Txyx; REAL Txyy; REAL Tyyx; REAL Tyyy;
	REAL r_aux; REAL dx; REAL dy;
	REAL Tmidp2p;
	REAL pi = acos(-1.),dpi = 2.*pi ,ddpi = 4.*pi;

	REAL *bx; REAL *by; bx = new double[n_pan]; by = new double[n_pan];

	int nc = xsrcSize/n_pan;

	for(int i=0; i < n_pan; i++){
		bx[i] = ux[i];
		by[i] = uy[i];
	}

	for(int q = 0; q < xtarSize; q++){
		for(int i=0; i < xsrcSize; i++){
			dx = xsrc[i] - xtar[q];
			dy = ysrc[i] - ytar[q];
			r_aux = pow(dx*dx + dy*dy, 0.5);

			Txxx = -4.*(dx*dx*dx)/pow(r_aux,4.);
			Txxy = -4.*(dx*dx*dy)/pow(r_aux,4.);
			
			Tyxx = Txxy;
			Tyxy = -4.*(dy*dx*dy)/pow(r_aux,4.);
			
			Txyx = Txxy;
			Txyy = Tyxy;

			Tyyx = Txyy;
			Tyyy = -4.*(dy*dy*dy)/pow(r_aux,4.);

			if(i/nc != q){
				bx[q] += (-1./dpi)*(Rp[i/nc]/2.)*wc[i]*(ux[i/nc]*(Txxx*nx[i/nc] 
					+ Txxy*ny[i/nc]) + uy[i/nc]*(Tyxx*nx[i/nc] + Tyxy*ny[i/nc]));

				by[q] += (-1./dpi)*(Rp[i/nc]/2.)*wc[i]*(ux[i/nc]*(Txyx*nx[i/nc] 
					+ Txyy*ny[i/nc])+ uy[i/nc]*(Tyyx*nx[i/nc] + Tyyy*ny[i/nc]));
			}			
		}
		if(q % 1000 == 0){
			cout<<"El calculo directo va en el punto de colocacion: "<< q << endl;
			clock_t midBuild	=	clock();
			Tmidp2p= double(midBuild-startBuild)/double(CLOCKS_PER_SEC);
			cout<<"han transcurrido: "<< Tmidp2p <<" [s] -> " << Tmidp2p/60 <<" [min]" <<endl;
		}
		if(q == 1000){
			double Tmidp2p_aux = Tmidp2p*(n_pan/1000); 
			cout<< endl << "tiempo estimado: "<< Tmidp2p_aux<< "[s] -> "<< Tmidp2p_aux/60 << "[min] -> "<< Tmidp2p_aux/3600 << "[h]\n"<<endl; 
		}
		//cout << "bx_" << q << " = " << bx[q] << endl;
	}

	for(int i=0; i < n_pan; i++){
		b[i] = bx[i];
		b[i+n_pan] = by[i];
	}

	delete[] bx, by;

}




//#endif //__DL_potential_dir_hpp__

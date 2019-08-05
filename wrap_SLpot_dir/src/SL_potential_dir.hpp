// #ifndef __SL_potential_dir_hpp__
//#define __SL_potential_dir_hpp__

#include"header/BBFMM2D.hpp"
#define REAL double

using namespace std;
using namespace Eigen;

void direct_SL_potential_cpp(REAL *xsrc, int xsrcSize,
			REAL *ysrc, int ysrcSize,
			REAL *xtar, int xtarSize,
			REAL *ytar, int ytarSize,
			REAL *Rp, int n_pan,
			REAL *wc, int n_wc,
			REAL *fx, int n_fx,
			REAL *fy, int n_fy,
			REAL Re,
			REAL *Ax, int n_Ax);

void direct_SL_potential_cpp(REAL *xsrc, int xsrcSize,
			REAL *ysrc, int ysrcSize,
			REAL *xtar, int xtarSize,
			REAL *ytar, int ytarSize,
			REAL *Rp, int n_pan,
			REAL *wc, int n_wc,
			REAL *fx, int n_fx,
			REAL *fy, int n_fy,
			REAL Re,
			REAL *Ax, int n_Ax){

	REAL pi = acos(-1.),dpi = 2.*pi ,ddpi = 4.*pi;

	clock_t startBuild	=	clock();

	REAL Gxx, Gxy;
	REAL Gyx, Gyy;
	REAL r_aux, dx, dy;
	REAL Tmidp2p;

	REAL *gx; REAL *gy;
	gx = new double[n_pan]; gy = new double[n_pan];

	for(int i = 0; i<n_pan;i++){
		gx[i] = 0.0;
		gy[i] = 0.0;
	}
	
	int N_sources = xsrcSize;
	int nc = n_wc/n_pan;
	
	for(int q = 0; q < n_pan; q++){
		for(int i=0; i < N_sources; i++){
			dx = xtar[q] - xsrc[i];
			dy = ytar[q] - ysrc[i];
			r_aux = pow(dx*dx + dy*dy, 0.5);

			Gxx = -log(r_aux) + dx*dx/(r_aux*r_aux);
			Gxy = dx*dy/(r_aux*r_aux);
			
			Gyx = Gxy;
			Gyy = -log(r_aux) + dy*dy/(r_aux*r_aux);
			
			
			if(i/nc != q){
				gx[q] += (Re/(ddpi))*wc[i]*Rp[i/nc]*(Gxx*fx[i/nc] +
					Gxy*fy[i/nc]);

				gy[q] += (Re/(ddpi))*wc[i]*Rp[i/nc]*(Gyx*fx[i/nc] +
					Gyy*fy[i/nc]);
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
	}

	for(int i=0; i < n_pan; i++){
		Ax[i] = gx[i];
		Ax[i+n_pan] = gy[i];
	}

	delete[] gx, gy;


}


//#endif //__DL_potential_hpp__

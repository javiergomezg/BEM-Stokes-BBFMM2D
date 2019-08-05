// #ifndef __SL_potential_dir_hpp__
//#define __SL_potential_dir_hpp__

#include"header/BBFMM2D.hpp"
#define REAL double

using namespace std;
using namespace Eigen;

void direct_pfield_cpp(REAL *xsrc, int xsrcSize,
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
			REAL *p_f, int n_p_f);

void direct_pfield_cpp(REAL *xsrc, int xsrcSize,
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
			REAL *p_f, int n_p_f){

	REAL pi = acos(-1.),dpi = 2.*pi ,ddpi = 4.*pi, octpi = 8.*pi;

	clock_t startBuild	=	clock();

	//REAL Gxx, Gxy;
	//REAL Gyx, Gyy;
	REAL px, py;

	//REAL Txxx; REAL Txxy; REAL Tyxx; REAL Tyxy;
	//REAL Txyx; REAL Txyy; REAL Tyyx; REAL Tyyy;

	REAL PIxx, PIxy, PIyx, PIyy;

	REAL r_aux, dx, dy;
	REAL Tmidp2p;

	int N_sources = xsrcSize;
	int nc = n_wc/n_pan;
	
	for(int q = 0; q < n_p_f; q++){
		for(int i=0; i < N_sources; i++){
			dx = xtar[q] - xsrc[i];
			dy = ytar[q] - ysrc[i];
			r_aux = pow(dx*dx + dy*dy, 0.5);

			//Gxx = -log(L*r_aux) + dx*dx/(r_aux*r_aux);
			//Gxy = dx*dy/(r_aux*r_aux);
			
			//Gyx = Gxy;
			//Gyy = -log(L*r_aux) + dy*dy/(r_aux*r_aux);

			px = 2.*(dx)/(r_aux*r_aux);
			py = 2.*(dy)/(r_aux*r_aux);

			dx = xsrc[i] - xtar[q];
			dy = ysrc[i] - ytar[q];
			r_aux = pow(dx*dx + dy*dy, 0.5);

			//Txxx = -4.*(dx*dx*dx)/pow(r_aux,4.);
			//Txxy = -4.*(dx*dx*dy)/pow(r_aux,4.);
			
			//Tyxx = Txxy;
			//Tyxy = -4.*(dy*dx*dy)/pow(r_aux,4.);
			
			//Txyx = Txxy;
			//Txyy = Tyxy;

			//Tyyx = Txyy;
			//Tyyy = -4.*(dy*dy*dy)/pow(r_aux,4.);
			
			PIxx = 4.*(-1./(r_aux*r_aux) + 2.*(dx*dx)/pow(r_aux,4.));

			PIxy = 8.*(dx*dy)/pow(r_aux,4.);
			PIyx = PIxy;

			PIyy = 4.*(-1./(r_aux*r_aux) + 2.*(dy*dy)/pow(r_aux,4.));
			
			//ux_f[q] += - (Re/(octpi))*wc[i]*Rp[i/nc]*(Gxx*fx[i/nc] +
			//	Gxy*fy[i/nc]) - (1./octpi)*(Rp[i/nc])*wc[i]*(ux[i/nc]*(Txxx*nx[i/nc] 
			//		+ Txxy*ny[i/nc]) + uy[i/nc]*(Tyxx*nx[i/nc] + Tyxy*ny[i/nc]));

			//uy_f[q] += - (Re/(octpi))*wc[i]*Rp[i/nc]*(Gyx*fx[i/nc] +
			//	Gyy*fy[i/nc]) - (1./octpi)*(Rp[i/nc])*wc[i]*(ux[i/nc]*(Txyx*nx[i/nc] 
			//		+ Txyy*ny[i/nc])+ uy[i/nc]*(Tyyx*nx[i/nc] + Tyyy*ny[i/nc]));

			p_f[q] += (-1./octpi)*wc[i]*Rp[i/nc]*(px*fx[i/nc] + 
					py*fy[i/nc]) + (1./(octpi*Re))*wc[i]*(ux[i/nc]*(PIxx*nx[i/nc]
					+ PIxy*ny[i/nc]) + uy[i/nc]*(PIyx*nx[i/nc] + PIyy*ny[i/nc]));

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

}


//#endif //__DL_potential_hpp__

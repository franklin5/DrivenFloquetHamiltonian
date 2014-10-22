/*
 * BdGmatrix.cpp
 *
 *  Created on: Oct 15, 2014
 *      Author: Lin Dong
 */
#include "floquet.h"
void cFloquet::update(int nk){
  if (nk == -1) {
    _bdg_H.setZero(); // This is done only once.
  } else {
    if (string(topo) == "bulk") {
      int nkx = nk % _NKX;     // --> the modulo (because nk = nkx+ nky * NKX )
      int nky = int (nk/_NKX); // --> the floor
      update(_gauss_k[nkx], _gauss_k[nky]);
    } else if (string(topo) == "edge"){
      update(_gauss_k[nk]);
    }
  }
}

void cFloquet::update(double kx, double ky){
  // Bulk property: Chern number
  // Wavefunction is expanded in frequency domain for time-dependent problem; p = 0 for time-independent problem.
  if (_ibdg == 4) {
    double xi = kx*kx + ky*ky - _mu;
    int p,q,m,n;
    int ip, iq, im, in;
    for (int ip1 = 0; ip1 < _pblock; ++ip1) {
      p  = ip1-_PMAX;
      ip = ip1*_ibdg*_NMAX;
      // only lower left part of the matrix is referenced for a self-adjoint matrix in Eigen's notation.
      for (int im1 = 0; im1 < _NMAX; ++im1) {
	m  = im1 + 1;
	im = im1*_ibdg;
	_bdg_H(ip+im,  ip+im)	= complex<double>(xi+_h+pow(m*M_PI/_L,2.0)+2*M_PI*p/_T,0.0);
	_bdg_H(ip+im+1,ip+im)   = complex<double>(_v*kx,_v*ky);
	_bdg_H(ip+im+1,ip+im+1)	= complex<double>(xi-_h+pow(m*M_PI/_L,2.0)+2*M_PI*p/_T,0.0);
	_bdg_H(ip+im+2,ip+im+2) = complex<double>(-(xi+_h+pow(m*M_PI/_L,2.0)+2*M_PI*p/_T),0.0);
	_bdg_H(ip+im+3,ip+im+2) = complex<double>(_v*kx,-_v*ky);
	_bdg_H(ip+im+3,ip+im+3) = complex<double>(-(xi-_h+pow(m*M_PI/_L,2.0)+2*M_PI*p/_T),0.0);
      }
      for (int iq1 = 0; iq1 < ip1; ++iq1){
	q = iq1-_PMAX;
	iq = iq1*_ibdg*_NMAX;
	_bdg_H(ip+im)
      }
    }
    
  } else if (_ibdg == 2){
    
  }
}
 
void cFloquet::update(double kx){
	// Edge property: spectrum under hard wall boundary condition.
	// Wavefunction is expanded in both y-direction and frequency domain.
	if (_ibdg == 4) {

	} else if (_ibdg == 2){

	}

}

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
	  int p, q;
	  // The off-diagonal coupling introduced from time-dependent order parameter should be computed only here.
	  complex<double> Gamma2;
	  double dt = 0.0005;
	  double t = 0.0;
	  VectorXd Delta_t(100000);
	  Delta_t.setZero();
	  int count = 0;
	  while (t<_T){
	    Delta_t(count) = 0.1-0.1*cos(2*M_PI*t/_T);
//		test_output << reD << '\t' << imD << endl;
	  	count++;
		t+=dt;
	  }
	  for (int i = 0; i < _pblock; ++i) {
	  		p = i-_PMAX;
	  		for (int j = 0; j < _pblock; ++j) {
	  			q = j-_PMAX;
	  			Gamma2 = complex<double> (0.0,0.0);
	  			t = 0.0;
	  			for (int ig = 0; ig < count; ++ig) {
					Gamma2 +=  abs(Delta_t(ig)) *
							complex<double> (cos(2*M_PI*(q-p)*t/_T),-sin(2*M_PI*(q-p)*t/_T));
					t += dt;
				}
	  			Gamma2 = Gamma2/_T*dt;
//	  			cout << Gamma2 << endl;
	  			_bdg_H(i+2*_pblock,j+_pblock) = Gamma2;
	  			_bdg_H(i+3*_pblock,j) = -Gamma2;
	  		}
	  }
  } else {
	  int nkx = nk % _NKX;     // --> the modulo (because nk = nkx+ nky * NKX )
	  int nky = int (nk/_NKX); // --> the floor
	  double kx = _gauss_k[nkx], ky = _gauss_k[nky];
	  update_kxky(kx,ky);


  }
}

void cFloquet::update_kxky(double kx, double ky){
	double xi = kx*kx + ky*ky - _mu;
	int p;
	for (int i = 0; i < _pblock; ++i) {
		p = i-_PMAX;
		// only lower left part of the matrix is needed for storage.
		_bdg_H(i,i) 	   = complex<double>(xi+_h+2*M_PI*p/_T,0.0);
		_bdg_H(i+_pblock,i)   = complex<double>(_v*kx,_v*ky);
		_bdg_H(i+_pblock,i+_pblock)= complex<double>(xi-_h+2*M_PI*p/_T,0.0);
		_bdg_H(i+2*_pblock,i+2*_pblock) = complex<double>(-(xi+_h+2*M_PI*p/_T),0.0);
		_bdg_H(i+3*_pblock,i+2*_pblock) = complex<double>(_v*kx,-_v*ky);
		_bdg_H(i+3*_pblock,i+3*_pblock) = complex<double>(-(xi-_h+2*M_PI*p/_T),0.0);
	}
}

/*
 * ioFloquet.cpp
 *
 *  Created on: Oct 15, 2014
 *      Author: ld7
 */
#include "floquet.h"
#include "lgwt.h"
void cFloquet::construction(){
	// read in parameters from file input.txt
	char dummyname[100];
	chernsolver = new char[100];
	double dummyvalue;
	FILE *input_open;
	input_open = fopen("input.txt","r");
	assert(input_open != NULL);
	cout.precision(16);
	if (_rank == _root) cout << "Starting to read in parameters from file input.txt" << endl;
	/* example file of input.txt with comment:
	 *
	 *
		hz      0.9		--> Zeeman field in z direction
		muInf   0.12	--> Chermical potential or the one extracted from quench data
		Tperiod 26.24	--> Oscillation period
		lambda  1.2		--> Spin-orbit coupling strength
		Length  200.0	--> System size when consider hard-wall boundary condition
		kmax    5.0		--> Momentum space integration cutoff for Chern number calculation

		PMAX    21		--> Frequency expansion basis cutoff. If not used, set as 1
		NMAX    1		--> Real space expansion basis cutoff. If not used, set as 1
		NKX     51		--> Good quantum number for spectrum or curvature plot

		chernsolver curvature	--> Solver for Chern number calculation: current available options include
								curvature or connection
								approach.

	 *
	 *
	 */
	fscanf(input_open,"%s %lf", dummyname, &dummyvalue);
	_h = dummyvalue;	if (_rank == _root) cout << dummyname << "=" << _h << endl;
	fscanf(input_open,"%s %lf", dummyname, &dummyvalue);
	_mu = dummyvalue;	if (_rank == _root) cout << dummyname << "=" << _mu << endl;
	fscanf(input_open,"%s %lf", dummyname, &dummyvalue);
	_T = dummyvalue;	if (_rank == _root) cout << dummyname << "=" << _T << endl;
	fscanf(input_open,"%s %lf", dummyname, &dummyvalue);
	_v = dummyvalue;	if (_rank == _root) cout << dummyname << "=" << _v << endl;
	fscanf(input_open,"%s %lf", dummyname, &dummyvalue);
	_L = dummyvalue;	if (_rank == _root) cout << dummyname << "=" << _L << endl;
	fscanf(input_open,"%s %lf", dummyname, &dummyvalue);
	_kmax = dummyvalue;	if (_rank == _root) cout << dummyname << "=" << _kmax << endl;
	fscanf(input_open,"%s %lf", dummyname, &dummyvalue);
	_PMAX = int(dummyvalue);	if (_rank == _root) cout << dummyname << "=" << _PMAX << endl;
	fscanf(input_open,"%s %lf", dummyname, &dummyvalue);
	_NMAX = int(dummyvalue);	if (_rank == _root) cout << dummyname << "=" << _NMAX << endl;
	fscanf(input_open,"%s %lf", dummyname, &dummyvalue);
	_NKX = int(dummyvalue);	if (_rank == _root) cout << dummyname << "=" << _NKX << endl;
	fscanf(input_open,"%s %s", dummyname, chernsolver);
	if (_rank == _root) cout << dummyname << "=" << chernsolver << endl;
	fclose(input_open);

	int ibdg = 4; // size of BdG matrix in operator form: 2 or 4 under consideration.
	_pblock = 2*_PMAX+1;
	_pblock4 = _pblock*ibdg;
	_SMAX = _NMAX*_pblock4;
	_NKX2 = _NKX*_NKX;

	_bdg_E.resize(_SMAX);
	_bdg_V.resize(_SMAX,_SMAX);
	_bdg_H.resize(_SMAX,_SMAX);

	update(-1); // null construction that is only done once for arbitrary momentum state.

	_gauss_k = new double [_NKX];  // integration momentum value
	_gauss_w_k = new double [_NKX];// integration momentum weight
	gauss_lgwt(_NKX,-_kmax,_kmax,_gauss_k,_gauss_w_k);
	// Legendre-Gauss quadrature method of integration:
	// code adapted from a matlab code in FileExchangeMATLAB: (lgwt.m)
	// http://www.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes/content/lgwt.m

	if (_rank == _root) cout << "Construction is completed." << endl;
}

void cFloquet:: destruction(){
	if (_root==_rank) {
		delete []sendbuf;
		delete []sendcounts;
		delete []displs;
		delete []TotalEig;
		delete []recvcounts;
		delete []displs_r;
	}
	delete []recvbuf;
	delete []localEig;
}

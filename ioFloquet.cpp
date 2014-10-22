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
	topo = new char[100];
	double dummyvalue;
	int intdummyvalue;
	FILE *input_open;
	//input_open = fopen("/home/ld7/workspace/DrivenFloquetHamiltonian/input.txt","r");
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
		Delta0  0.1		--> Magnitude of order parameter
		lambda  1.2		--> Spin-orbit coupling strength
		Length  200.0	--> System size when consider hard-wall boundary condition
		kmax    5.0		--> Momentum space integration cutoff for Chern number calculation

		ibdg    4		--> size of BdG matrix in operator form: 2 or 4 under consideration.
		PMAX    21		--> Frequency expansion basis cutoff. If not used, set as 1
		NMAX    1		--> Real space expansion basis cutoff. If not used, set as 1
		NKX     51		--> Good quantum number for spectrum or curvature plot

		chernsolver curvature	--> Solver for Chern number calculation: current available options include
										curvature or connection
										^^^^^^^^^    ^^^^^^^^^^
		topo            bulk	--> topological property that are under consideration: options include
										bulk or edge
										^^^^    ^^^^
		mu              1.0		--> Starting to initialize parameters used in the PRX paper.
		J               1.5
		b               1.5
		a               4.0
		Delta0          1.0
		omega           14.2857
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
	_Delta0 = dummyvalue;	if (_rank == _root) cout << dummyname << "=" << _Delta0 << endl;
	fscanf(input_open,"%s %lf", dummyname, &dummyvalue);
	_v = dummyvalue;	if (_rank == _root) cout << dummyname << "=" << _v << endl;
	fscanf(input_open,"%s %lf", dummyname, &dummyvalue);
	_L = dummyvalue;	if (_rank == _root) cout << dummyname << "=" << _L << endl;
	fscanf(input_open,"%s %lf", dummyname, &dummyvalue);
	_kmax = dummyvalue;	if (_rank == _root) cout << dummyname << "=" << _kmax << endl;

	fscanf(input_open,"%s %lf", dummyname, &intdummyvalue);
	_ibdg = intdummyvalue;	if (_rank == _root) cout << dummyname << "=" << _ibdg << endl;
	if (_ibdg != 2 || _ibdg != 4 ) {
		cerr << "BdG matrix size is not supported. Please check input.txt and readme." << endl;
		exit(1);
	}
	fscanf(input_open,"%s %d", dummyname, &intdummyvalue);
	_PMAX = intdummyvalue;	if (_rank == _root) cout << dummyname << "=" << _PMAX << endl;
	fscanf(input_open,"%s %d", dummyname, &intdummyvalue);
	_NMAX = intdummyvalue;	if (_rank == _root) cout << dummyname << "=" << _NMAX << endl;
	fscanf(input_open,"%s %d", dummyname, &intdummyvalue);
	_NKX = intdummyvalue;	if (_rank == _root) cout << dummyname << "=" << _NKX << endl;

	fscanf(input_open,"%s %s", dummyname, chernsolver);
	if (string(chernsolver) != "curvature" || string(chernsolver) != "connection") {
		cerr << "Solver for Chern number is not supported. Please check input.txt and readme." << endl;
		exit(1);
	}
	if (_rank == _root) cout << dummyname << "=" << chernsolver << endl;
	fscanf(input_open,"%s %s", dummyname, topo);
	if (_rank == _root) cout << dummyname << "=" << topo << endl;
	if (_ibdg == 2) { // PRX parameters
		_mu = dummyvalue;	if (_rank == _root) cout << dummyname << "=" << _mu << endl;
		fscanf(input_open,"%s %lf", dummyname, &dummyvalue);
		_J = dummyvalue;	if (_rank == _root) cout << dummyname << "=" << _J << endl;
		fscanf(input_open,"%s %lf", dummyname, &dummyvalue);
		_b = dummyvalue;	if (_rank == _root) cout << dummyname << "=" << _b << endl;
		fscanf(input_open,"%s %lf", dummyname, &dummyvalue);
		_a = dummyvalue;	if (_rank == _root) cout << dummyname << "=" << _a << endl;
		fscanf(input_open,"%s %lf", dummyname, &dummyvalue);
		_Delta0 = dummyvalue;	if (_rank == _root) cout << dummyname << "=" << _Delta0 << endl;
		fscanf(input_open,"%s %lf", dummyname, &dummyvalue);
		_omega = dummyvalue;	if (_rank == _root) cout << dummyname << "=" << _omega << endl;
		_T = 2*M_PI/_omega;
	}
	fclose(input_open);

	_pblock = 2*_PMAX+1;
	_pblock4 = _pblock*_ibdg;
	_SMAX = _NMAX*_pblock4;
	if (string(topo) == "bulk") {
		_NKX2 = _NKX*_NKX;
	} else if (string(topo) == "edge"){
		_NKX2 = _NKX;
	} else {
		cerr << "Please specify available option: either bulk or edge." << endl;
		exit(1);
	}

	_bdg_E.resize(_SMAX);
	_bdg_V.resize(_SMAX,_SMAX);
	_bdg_H.resize(_SMAX,_SMAX);

	update(-1); // null construction that is only done once for arbitrary momentum state.

	_gauss_k = new double [_NKX];  // integration momentum value
	_gauss_w_k = new double [_NKX];// integration momentum weight
	gauss_lgwt(_NKX,-_kmax,_kmax,_gauss_k,_gauss_w_k); // TODO: option of doing 1D integral to be updated.
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
	delete []_gauss_k;
	delete []_gauss_w_k;
	delete []chernsolver;
	delete []topo;
}

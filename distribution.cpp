/*
 * distribution.cpp
 *
 *  Created on: Oct 15, 2014
 *      Author: ld7
 */

#include "floquet.h"

int cFloquet::compute_count(int rank, int size){
  // compute distribution count: the last rank does the remaineder job while the rest do the most even work.
  // _NKX2 is the total number of work to be distributed,
  // because momentum is good quantum number and each work is indepedent with each other.
  int result;
  if (rank != size-1) {
    result = int(_NKX2/size);
  } else {
    result = int(_NKX2/size) + _NKX2 % size;
 }
  return result;
}

void cFloquet::distribution(){
	  // The paradigm is to send momentum value to different rank
	  // and let each processor work on the big matrix independently.
	  // --> A two-level parallization is to be investigated in the furture
	  //     for a sparse matrix partial eigenvalue spectrum calculation.
  if (_rank == _root){ // send process is only root significant
    sendbuf = new int[_NKX2];
    for(int i = 0; i< _NKX2; ++i){
      sendbuf[i] = i;
    }
    sendcounts = new int[_size];
    displs = new int[_size];
    for(int i=0; i<_size; i++){
      sendcounts[i] = compute_count(i,_size);
      displs[i] = i*int(_NKX2/_size);
    }
  }
  recvcount = compute_count(_rank,_size); // This is a rank dependent variable.
  recvbuf = new int[recvcount]; // So is this array: rank dependent size
  MPI_Scatterv(sendbuf,sendcounts,displs,MPI_INT,recvbuf,recvcount,MPI_INT,_root,COMM_WORLD);
  stride = _SMAX*recvcount;
  localEig = new double[stride];
  SelfAdjointEigenSolver<MatrixXcd> ces;
  for(int ig = 0; ig<_size; ++ig) {
	if (ig ==_rank){
	  for (int i=0; i<recvcount; ++i) {
	update(recvbuf[i]);
	clock_t start = clock();
	ces.compute(_bdg_H,0); // only eigenvalues are computed.
	//	ces.compute(_bdg_H); // eigenvectors are also computed.
	//cout << _bdg_H << endl;
	clock_t end = clock();
	cout << "rank = " << ig << " is doing job recvbuf[" << i << "] out of " << recvcount << "using time " <<  double (end-start)/ (double) CLOCKS_PER_SEC << endl;
	for(int j = 0; j < _SMAX; ++j){
	  localEig[j+i*_SMAX]=ces.eigenvalues()[j];
	  //	  cout << localEig[j+i*pblock4] << " ";
	}
	//	cout << endl;
	  }
	  if (_rank==_root){
	cout << "rank " << _rank << " has finished "<< recvcount << "tasks out of total" << _NKX2 << endl;
	  }
	}
  }
}


void cFloquet:: aggregation(){

	if (_root==_rank) {
		TotalEig = new double [_NKX2*_SMAX];
		recvcounts = new int[_size];
		displs_r = new int[_size];
		offset = 0;
		for(int ig=0;ig<_size;++ig){
		  recvcounts[ig] = compute_count(ig,_size)*_SMAX;
		  displs_r[ig] = offset;
		  offset += recvcounts[ig];
		  //cout << offset << " ";
		}
		//    cout << endl;
	  }
	  MPI_Gatherv(localEig, stride, MPI_DOUBLE, TotalEig, recvcounts, displs_r, MPI_DOUBLE, _root, COMM_WORLD);
	  if (_rank == _root) {
		int itemp;
		ofstream bdg_output;
		bdg_output.open("spectrum.OUT"); // TODO: modify output file name
		assert(bdg_output.is_open());
		for(int i =0;i<_size;++i){
		  itemp = compute_count(i,_size);
		  if( i != _size-1) {
		offset = itemp;
		  } /*else {
		offset = offset; // offset is updated as size-2 position
		  }*/
		  //      cout << itemp << endl;
		  for (int j = 0; j<itemp;++j){
		for (int q = 0; q<_SMAX;++q){
	//	  cout << TotalEig[i*offset*pblock4+j*pblock4+q] << '\t';
		  //	  cout << i*offset*pblock4+j*pblock4+q << '\t';
		  bdg_output << TotalEig[i*offset*_SMAX+j*_SMAX+q] << '\t';
		}
	//	cout << endl;
		bdg_output << endl;
		  }
		  //      bdg_output << endl;
		}
		bdg_output.close();
	  }
}






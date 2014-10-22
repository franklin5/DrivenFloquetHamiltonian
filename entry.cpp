/*
 * main.cpp
 *
 *  Created on: Oct 15, 2014
 *      Author: Lin Dong
 */

#include "stdcpp.h"
#include "floquet.h"
#define root 0
int main(int argc, char** argv){
//	cout << argc << ' ' << argv[0] << argv[1] << endl;
	Init(argc,argv);
	int rank, size;
	rank = COMM_WORLD.Get_rank();
	size = COMM_WORLD.Get_size();
	if (rank == root) {
		cout << "======================================================================\n"
			 << "The purpose of this program is to study\n"
			 << "the Floquet Topological property of a periodically driven Hamiltonian.\n"
			 << "Motivated, proposed, designed, implemented and researched \n"
			 << "by Lin Dong at Rice University. \n"
			 << "at " << __TIME__ << ", on " << __DATE__ << endl
			 << "MPI is initialized and program starts from " << __FILE__ << endl
			 << "======================================================================\n";
	}
	cFloquet Floquet(rank,size,root);
	Floquet.construction();
	Floquet.distribution();
	Floquet.aggregation();
	Floquet.destruction();
	Finalize();
	return 0;
}

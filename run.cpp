#include "pregel_kmeans.h"

int main(int argc, char* argv[]){
	init_workers();
	pregel_hashmin("/kmeansToy", "/kmeansToyOut", 5);
	worker_finalize();
	return 0;
}

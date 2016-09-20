#include <iostream>
#include <ctime>
#include "ga_util.h"


	int main() {

		// init rand()
		srand((unsigned int)time(NULL));

		std::cout << "Start GA Equation Solver" << std::endl;

		GaGenaration generation;
		generation.initGeneration(POP, GTYPE_LEN, GTYPE_MAX, PTYPE_MIN, PTYPE_MAX);

		generation.evaluation();

		for (int i = 0; i < NUM_OF_GENERATION; i++) {
			generation.selection();
			generation.coutResult(i);
			generation.crossover();
			generation.mutation();
			generation.evaluation();
		}

		return 0;
	}
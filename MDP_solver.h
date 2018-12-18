#ifndef MDP_SOLVER
#define MDP_SOLVER

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string.h>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

class MDP
{
	public:
		// action index : 0 - stand, 1 - hit, 2 - dd, 3 - split
		int states;	//num of states
		double* values_sh; //values - stand hit
		double* values_shd; //values - stand hit doubledown
		double** P; //Hit transition matrix (281 * 281)
		double* Qstand; //Q for stand for all states
		double* Qhit; //Q for hit for all states
		double* Qdd; //Q for dd for all states
		int* Pi; //Policy

		MDP(int states, double** P, double* Qstand);
		void value_iterations(double epsilon);

};

#endif

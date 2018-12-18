#ifndef STATE
#define STATE

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string.h>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

class State
{
	public:
		int player_val;	//value of hand
		int dealer_upCard;	//Dealer up card
		bool isSoft;  //A = 11 or not
		double value; //value function of state
		int type;
		int index;

		State();
		State(int type);
		State(int player_val, int dealer_upCard);

};

#endif

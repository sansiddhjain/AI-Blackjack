#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <string>
#include "State.h"
#include "MDP_solver.h"

// double p; // Probability of getting a face card
using namespace std;

double** calculateQValueStand(double p); // Calculates q value for stand for various states
// void printMatrix(double** matrix, int N, int M); // Prints Matrix
double** calculateStateReward(double p); //

State calculateNewState (State current, int cardDrawn);
int getIndex(State state);
State getStateFromIndex(int index);
double** getHitMatrix(double p);
double* getQstand(double** Q);
int** getPairOptimalAction(double* Qstand, double* Qhit, double* Qdd, double* values, double** P, double p, int* Pi);

char convertToLetter(int i);
void printOutput(int* mdp_policy, int** pair_policy);

static int N = 281;

int main(int argc, char const *argv[]) {
  double p = atof(argv[1]);
	double epsilon = 1e-10;
  std::cout << p << '\n';

  double** rewardMatrix = calculateStateReward(p);

	cout << "Printing Reward matrix" <<"\n";
  for (int i = 0; i < 20; i++)
  {
    for (int j = 0; j < 26; j++) {
      cout <<  rewardMatrix[i][j] << ", ";
    }
    cout << "" << '\n';
  }

  rewardMatrix = calculateQValueStand(p);

	double** HitMatrix = getHitMatrix(p);

  // for (int i = 0; i < 281; i++) {
  //   for (int j = 0; j < 281; j++) {
  //     cout <<  HitMatrix[i][j] << ", ";
  //   }
  //   cout << "" << '\n';
  // }

	double* Qstand = getQstand(rewardMatrix);
  std::cout << "Qstand" << '\n';
	for (int i = 0; i < 281; i++) {
		cout <<  Qstand[i] << ", ";
    if (i % 10 == 9) {
      std::cout << '\n';
    }
  }
	cout << endl;

	MDP mdp = MDP(N, HitMatrix, Qstand);
	cout << "Initialised mdp" << "\n";
	mdp.value_iterations(epsilon);

	int** pair_output = getPairOptimalAction(mdp.Qstand, mdp.Qhit, mdp.Qdd, mdp.values_shd, HitMatrix, p, mdp.Pi);


	for(int i = 30; i < mdp.states-1; i++)
	{
		cout << mdp.Pi[i] << ", ";
    if (i%10 == 9)
      std::cout << '\n';
    if (i%190 == 9)
      std::cout << '\n';
	}
	cout << endl;

	cout<<"Soft 13 D 5:"<<"\n";
  std::cout << mdp.Qstand[204] << '\n';
  std::cout << mdp.Qhit[204] << '\n';
  std::cout << mdp.Qdd[204] << '\n';
	std::cout << mdp.values_sh[204] << '\n';
	std::cout << mdp.values_shd[204] << '\n';

	cout<<"Soft 15 D 4:"<<"\n";
  std::cout << mdp.Qstand[223] << '\n';
  std::cout << mdp.Qhit[223] << '\n';
  std::cout << mdp.Qdd[223] << '\n';
	std::cout << mdp.values_sh[223] << '\n';
	std::cout << mdp.values_shd[223] << '\n';

  cout << "Printing value_sh matrix at end" << "\n";
	for (int i = 0; i < 281; i++) {
		cout <<  mdp.values_sh[i] << ", ";
    if (i % 10 == 9) {
      std::cout << '\n';
    }
  }
	cout << endl;

  cout << "Printing value_shd matrix at end" << "\n";
	for (int i = 0; i < 281; i++) {
		cout <<  mdp.values_shd[i] << ", ";
    if (i % 10 == 9) {
      std::cout << '\n';
    }
  }
	cout << endl;

  // calulate the q values for soft 12 (Which is (A, A))

  double* q_hit_s12 = new double[10];
  double* q_dd_s12 = new double[10];
  double* q_stand_s12 = new double[10];
  double* q_split_s12 = new double[10];

  for (int i = 0; i < 10; i++) {
    q_hit_s12[i] = 0;
    q_stand_s12[i] = 0;
    q_dd_s12[i] = 0;
    q_split_s12[i] = 0;
  }

  // Calculate q_stand, q_hit, and q_dd for Soft 12 first
  // i goes from 0-9 : represents dealer card going from A to 10
  for (int i = 0; i < 10; i++) {

    // Updating qStand for soft 12
    q_stand_s12[i] = mdp.Qstand[100+i];

    for (int j = 0; j < 10; j++) {
      if (j < 8) { // Drawn card A-8, we go to states Soft 13-20
        q_hit_s12[i] += ((1-p)/9.0)*mdp.values_sh[200+j*10+i];
        q_dd_s12[i] += 2*((1-p)/9.0)*mdp.Qstand[200+j*10+i];
      }
      if (j == 8) { // Drawn card 9, we go to Hard 21
        q_hit_s12[i] += ((1-p)/9.0)*mdp.values_sh[190+i];
        q_dd_s12[i] += 2*((1-p)/9.0)*mdp.Qstand[190+i];
      }
      if (j == 9) { // Drawn card 10, we go to Hard 12
        q_hit_s12[i] += p*mdp.values_sh[100+i];
        q_dd_s12[i] += 2*p*mdp.Qstand[100+i];
      }
    }
  }

  // Now calculate q_split for Soft 12
  // i goes from 0-9 : represents dealer card going from A to 10
    for (int i = 0; i < 10; i++) {
      // Drawn card varies from A-10, index varies from 0-9
      for (int j = 0; j < 10; j++) {
        if (j == 0) { // If drawn card is A, we go to Soft 12
          double temp = max(q_dd_s12[i], q_dd_s12[i]/2);
          q_split_s12[i] += ((1-p)/9.0)*max(q_stand_s12[i], temp);
        }
        if (j > 0 && j < 9) { // If drawn card is 2-9, we go to Soft 13-20
          int index = 200 + 10*(j-1) + i;
          double temp = max(mdp.Qdd[index], mdp.Qdd[index]/2);
          q_split_s12[i] += ((1-p)/9.0)*max(mdp.Qstand[index], temp);
        }
        if (j == 9) { // If drawn card is 10, we go to Hard 21
          int index = 190 + i;
          double temp = max(mdp.Qdd[index], mdp.Qdd[index]/2);
          q_split_s12[i] += p*max(mdp.Qstand[index], temp);
        }
      }
      q_split_s12[i] = 2*q_split_s12[i];
    }

    for (size_t j = 0; j < 10; j++) {
      double max = q_stand_s12[j];
      pair_output[9][j] = 0;
      if (q_hit_s12[j] > max) {
        max = q_hit_s12[j];
        pair_output[9][j] = 1;
      }
      if (q_dd_s12[j] > max) {
        max = q_dd_s12[j];
        pair_output[9][j] = 2;
      }
      if (q_split_s12[j] > max) {
        max = q_split_s12[j];
        pair_output[9][j] = 3;
      }
    }

    cout<<"Optimal policy for pairs:   "<<"\n";
  	for (int i = 0; i < 10; i++) {
      for (int j = 0; j < 10; j++) {
        cout <<  pair_output[i][j] << ", ";
      }
      cout << "" << '\n';
    }
  	cout<<""<<"\n";

    printOutput(mdp.Pi, pair_output);
}


// Calculates q value of standing, given the player value is x, and the dealer upcard is y
double** calculateQValueStand(double p)
{
  // Rows = 2-21, Soft 13-20
	double** rewardMatrix = calculateStateReward(p);
	double** qValueStand = new double*[20];
	for (int i = 0; i < 20; i++)
		qValueStand[i] = new double[10];

	for (int j = 0; j < 10; j++)
	{
		for (int i = 0; i < 20; i++)
		{
			if( j!= 0 )
				qValueStand[i][j] = rewardMatrix[i][j-1];
			else
				qValueStand[i][j] = 0;
		}
	}

	for (int i = 0; i < 20; i++)
	{
		for (int j = 15; j < 26; j++)
		{
			if (j == 20)
				continue;
			if (j == 19)
				qValueStand[i][0] += p*rewardMatrix[i][j];
			else
				qValueStand[i][0] += ((1-p)/9.0)*rewardMatrix[i][j];
		}
	}

	qValueStand[19][0] -= p; // State (21, A)
	qValueStand[19][9] -= (1-p)/9; // State (21, 10)

	return qValueStand;
}

// Function calculates the reward (value function v) for being in state (x, y)
// where x is the player's total value, and y the dealer's value
// ** This is after the player has chosen to stand, of course **
void printRewardMatrix(double** matrix, int N, int M) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < M; j++) {
      cout <<  matrix[i][j] << ", ";
    }
    cout << "" << '\n';
  }
}

double** calculateStateReward(double p)
{
  int rowSize = 20; // Player value varies from 2-21
  int columnSize = 26; // Dealer total can be 2-21, Bust, Soft 13-16 UPDATE - Soft 12-16
  // For reward matix, index : normal = face value - 2, Bust = 20, Soft 13-16 = 21-25
  double** rewardMatrix = new double*[20];
  for (int i = 0; i < rowSize; i++)
	{
    rewardMatrix[i] = new double[26];
  }
	std::cout << "Initialised Reward Matrix" << '\n';
	// Initialise Reward Matrix
	for (int i = 0; i < 20; i++)
  {
    for (int j = 0; j < 26; j++)
    	rewardMatrix[i][j] = 0;
  }
  std::cout << "Filled it with Zeros" << '\n';
  for (int i = 0; i < 20; i++)
  {
    for (int j = 15; j < 21; j++)
    {
      if (i == j)
        rewardMatrix[i][j] = 0;
      else if (j == 20)
        rewardMatrix[i][j] = 1;
      else if (i < j)
        rewardMatrix[i][j] = -1;
      else
        rewardMatrix[i][j] = 1;
    }
  }

  std::cout << "Filled the initial few columns" << '\n';
  double** probabMatrix = new double*[26]; // Only 2-16, and Soft 13-16
  // (Indices 15-20 correspond to states 17-21 and Bust, which are the states
  // from which you will NEVER TRANSITION FROM. Therefore, the are left empty)
  for (int i = 0; i < 26; i++) {
    probabMatrix[i] = new double[26]; // Dealer states same as before
  }
  std::cout << "Initialised Probability Matrix" << '\n';
	// Initialise probabMatrix
  for (int i = 0; i < 26; i++)
  {
    for (int j = 0; j < 26; j++)
    	probabMatrix[i][j] = 0;
  }
  std::cout << "Filled it with Zeros" << '\n';
	// Fill probability matrix up
  for(int i = 25; i >= 0; i--)
  {
    if (i>=15 && i <=20)
      continue;
		// Next if block - corresponding to probability transition matrix for soft states
		// IFFY - todo : check if unncessary (maybe use Sanyam's getNextState function)
		if (i >= 21) {
		  for(int j = 1; j < 11; j++)
		  {
			  if ((i)+j < 26)
			  {
				  probabMatrix[i][i+j] = (1-p)/9.0;
			  }
			  else if ((i)+j < 31)
			  {
					probabMatrix[i][(i)+j-11] = (1-p)/9.0;
			  }
				else
			  {
					if (j == 10)
						probabMatrix[i][(i)+j-21] = p;
					else
						probabMatrix[i][(i)+j-21] = (1-p)/9.0;
			  }
		  }
      std::cout << "Filled columns corresponding to soft" << '\n';
	  }
		// If block for when dealer total is between 11 and 16
	  else if (i >= 9) {
		  double sum = 0;
		  for(int j = 1; j < 11; j++)
		  {
			  if (i+j < 20)
			  {
				  if(j != 10)
						probabMatrix[i][i+j] = (1-p)/9.0;
				  else
						probabMatrix[i][i+j] = p;
				  sum += probabMatrix[i][i+j];
			  }
			  if (i+j == 20)
			  {
				  probabMatrix[i][i+j] = 1 - sum;
				  break;
			  }
		  }
      std::cout << "Filled columns where dealer total is 11-16" << '\n';
	  }

		// If block for when dealer total is between 6 and 10
	  else if (i >= 4) {
		  for(int j = 2; j <= 11; j++)
		  {
				if(j != 10)
					probabMatrix[i][i+j] = (1-p)/9.0;
				else
					probabMatrix[i][i+j] = p;
		  }
      std::cout << "Filled columns where dealer total is 6-10" << '\n';
	  }

		// If block for when dealer total is between 2 and 5
	  else {
		  for(int j = 2; j <= 11; j++)
		  {
				if(j == 10)
					probabMatrix[i][i+j] = p;
				else if (j == 11)
					probabMatrix[i][i+j+10] = (1-p)/9.0;
				else
					probabMatrix[i][i+j] = (1-p)/9.0;
		  }
      std::cout << "Filled columns where dealer total is 2-15" << '\n';
	  }
  }

  // for (int i = 0; i < 26; i++)
  // {
  //   for (int j = 0; j < 26; j++) {
  //     cout <<  probabMatrix[i][j] << ", ";
  //   }
  //   cout << "" << '\n';
  // }

	// Fill Reward Matrix up

	// First fill Hard 16 to Hard 6 (including), in reverse order
	for(int j = 14; j >= 4; j--)
	{
		for(int i = 0; i < 20; i++)
		{
			for(int k = j+1; k < 21; k++)
				rewardMatrix[i][j] += probabMatrix[j][k]*rewardMatrix[i][k];
		}
	}
  std::cout << "Filled columns of rewardMatrix where dealer total is Hard 6-16" << '\n';

	// Then fill up Soft 12 to Soft 16
	for(int j = 25; j >= 21; j--)
	{
		for(int i = 0; i < 20; i++)
		{
			for(int k = 0; k < 26; k++)
				rewardMatrix[i][j] += probabMatrix[j][k]*rewardMatrix[i][k];
		}
	}
  std::cout << "Filled columns of rewardMatrix where dealer total is soft 13-16" << '\n';

	// At last fill up the first 4 columns - reward for being in state (x, y) where y : 2-5
	for(int j = 3; j >= 0; j--)
	{
		for(int i = 0; i < 20; i++)
		{
			for(int k = j+1; k < 26; k++)
				rewardMatrix[i][j] += probabMatrix[j][k]*rewardMatrix[i][k];
		}
	}
  std::cout << "Filled columns of rewardMatrix where dealer total is Hard 2-5" << '\n';

	return rewardMatrix;
}

State calculateNewState (State current, int cardDrawn)
{
	State newState = State();
  newState.dealer_upCard = current.dealer_upCard;
	newState.player_val = current.player_val + cardDrawn;
	newState.isSoft = current.isSoft;

	if (cardDrawn == 1)
	{
		newState.player_val += 10;
		newState.isSoft = true;
	}
	if (newState.isSoft == true && newState.player_val > 21)
	{
		newState.player_val += -10;
		if (current.isSoft == true && cardDrawn == 1)
		{
			newState.isSoft = true;
		}
		else
		{
			newState.isSoft = false;
		}
	}
	if (newState.player_val > 21)
	{
		newState = State(1);
	}

	return newState;
}

int getIndex(State state)
{
  if (state.type == 1)
	{
		return 280;
	}
  else
  {
		if ((state.isSoft && state.player_val == 21) || !state.isSoft)
  	{
  		int x = state.player_val;
  		int y = state.dealer_upCard;
  		return ((10*(x-2))+y-1);
  	}
  	else
  	{
  		int x = state.player_val;
  		int y = state.dealer_upCard;
      int res = ((20*(x-3)) - (10*(x-13)) +y-1);
      return res;
  	}
  }
}

State getStateFromIndex(int index)
{
	if (index < 200)
	{
		int x = (index/10)+2;
		int y = (index%10)+1;
		State state = State(x, y);
		state.isSoft = false;
		state.type = 0;
		return state;
	}
	else if (index >= 200 && index < 280)
	{
		int x = (index/10) - 7;
		int y = (index%10) + 1;
		State state = State(x, y);
		state.isSoft = true;
		state.type= 0;
		return state;
	}
	else
	{
		State state = State(1);
		return state;
	}
}

double** getHitMatrix(double p)
{
  double** HitMatrix = new double*[N];
  for (int i = 0; i < N; ++i)
  	HitMatrix[i] = new double[N];

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    	HitMatrix[i][j] = 0;
  }

	for (int i = 0; i < N; i++)
	{
		if (i == 280)
    {
			HitMatrix[i][i] = 1.0;
      continue;
    }
		else
		{
			for (int k = 1; k<11; k++)
			{
				State current = getStateFromIndex(i);
				State newState = calculateNewState(current, k);

				int newIndex = getIndex(newState);
				// if (i == 244 && k ==4){
				// 	cout<<"new index is - "<<newIndex<<"\n";
				// 	cout<<newState.player_val<<" "<<newState.dealer_upCard<<"\n";
				// }

				if (k == 10)
					HitMatrix[i][newIndex] = HitMatrix[i][newIndex]+p;
				else
				{
					HitMatrix[i][newIndex] = HitMatrix[i][newIndex]+((1-p)/9);
				}
			}
		}
	}
	return HitMatrix;
}

double* getQstand(double** Q)  //Q is 20*10 dimensional
{
	double* Qstand = new double[N];
	for (int i = 0; i < 20; i++)
	{
		for (int j = 0; j < 10; j++)
		{
			Qstand[i*10+j] = Q[i][j];
		}
	}
	for (int k = 200; k < 280; k++)
	{
		Qstand[k] = Qstand[k-90];
	}
	Qstand[280] = -1;
	return Qstand;
}

int** getPairOptimalAction(double* Qstand, double* Qhit, double* Qdd, double* values, double** P, double p, int* Pi)
{
	int** PairOptimal = new int*[10];
	double** values_pairs = new double*[10];
	double** Q_pairs = new double*[10];

	for (int i = 0; i < 10; i++)
	{
		PairOptimal[i] = new int[10];
		values_pairs[i] = new double[10];
		Q_pairs[i] = new double[10];
	}

	for (int i = 0; i < 9; i++)
	{
		int x = i+2;
		for (int j = 0; j < 10; j++){
			double temp = 0.0;
			int start_index = getIndex(State(x, j+1));
			for (int k = 1; k < 11; k++)
			{
				if (k == i+2){
					continue;
				}
				State newState = calculateNewState(State(x, j+1), k);
				int index = getIndex(newState);
				// if (i == 1 && j == 1){
				// 	cout << "New index is - "<<index<<"\n";
				// }

				if (i == 8 && k == 1) {
					if(0 < j && j < 9){
						temp += 2*P[start_index][index]*1.5;
					}
					else{
						if(j == 0){
							temp += 2*P[start_index][index]*1.5*(1-p);
						}
						if(j == 9){
							temp += 2*P[start_index][index]*1.5*(1-((1-p)/9));
						}
					}
				}
				else{
					temp += 2*P[start_index][index]*values[index];
				}
			}
			if (i == 8)
			{
				// cout<<i<<" "<<j<<" "<<temp<<"\n";
				values_pairs[i][j] = temp/(1-(2*p));
				//Q_pairs[i][j] = temp+(2*p*values_pairs[i][j]);
				Q_pairs[i][j] = values_pairs[i][j];

			}
			else
			{
				values_pairs[i][j] = temp/(1-(2*(1-p)/9.0));
				//Q_pairs[i][j] = temp + ((2*(1-p)/9) * values_pairs[i][j]);
				Q_pairs[i][j] = values_pairs[i][j];
			}

			// cout<<x<<" "<<j<<"\n";
			int index1 = getIndex(State(2*x, j+1));
			// cout<<"Index     ------    "<<index1<<"\n";
			double val = values[index1];
			// cout<<"Value     ------    "<<val<<"\n";
			int action = Pi[getIndex(State(2*x, j+1))];
			double max = 0.0;
			if (action == 0)
			{
				max = Qstand[getIndex(State(2*x, j+1))];
			}
			else if (action == 1)
			{
				max = Qhit[getIndex(State(2*x, j+1))];
			}
			else
			{
				max = Qdd[getIndex(State(2*x, j+1))];
			}

			if((j == 1))
			{
				cout << "Pair" << x << ", " << x << ":, D" << j+1 << endl;
				cout << "Qstand " << Qstand[getIndex(State(2*x, j+1))] << endl;;
				cout << "Qhit " << Qhit[getIndex(State(2*x, j+1))] << endl;;
				cout << "Qdd " << Qdd[getIndex(State(2*x, j+1))]  << endl;
				cout << "Qsplit " << Q_pairs[i][j] << endl;
				cout << endl;
			}

			if (max < Q_pairs[i][j])
			{
				PairOptimal[i][j] = 3;
			}
			else
			{
				PairOptimal[i][j] = action;
			}
		}
	}

    return PairOptimal;
}

char convertToLetter(int i)
{
  switch (i) {
    case 0:
      return 'S';
    case 1:
      return 'H';
    case 2:
      return 'D';
    default:
      return 'P';
  }
}

void printOutput(int* mdp_policy, int** pair_policy)
{
  std::ofstream out("Policy.txt");
  std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(out.rdbuf());
  // Hard 5 - 19
  for (int i = 0; i < 15; i++)
  {
    std::cout << i+5 << "\t";
    for (int j = 1; j < 11; j++) {
      std::cout << convertToLetter(mdp_policy[30+i*10 + (j%10)]) << ' ';
    }
    std::cout << "" << '\n';
  }
  // std::cout << "" << '\n';

  // Soft 13 - 20
  for (int i = 0; i < 8; i++)
  {
    std::cout << "A" << i+2 << "\t";
    for (int j = 1; j < 11; j++) {
      std::cout << convertToLetter(mdp_policy[200+i*10 + (j%10)]) << ' ';
    }
    std::cout << "" << '\n';
  }
  // std::cout << "" << '\n';

  // Pairs (2, 2) - (A, A)
  for (int i = 0; i < 10; i++)
  {
    if (i < 9)
      std::cout << i+2 << i+2 << "\t";
    else
      std::cout << "A" << "A" << "\t";
    for (int j = 1; j < 11; j++) {
      std::cout << convertToLetter(pair_policy[i][j%10]) << ' ';
    }
    if (i < 9)
      std::cout << "" << '\n';
  }
   std::cout.rdbuf(coutbuf);
}

#include "MDP_solver.h"
#include "State.h"
#include <ctime>
#include <cstdlib>


MDP::MDP (int states, double** P, double* Qstand)
{
    this->states = states;
    this->P = P;
	this->Qstand = Qstand;
	this->Qhit = new double[states-1];
	this->Qdd = new double[states-1];
	this->values_sh = new double[states];
	this->values_shd = new double[states];
	this->Pi = new int[states-1];
	srand(time(NULL));

	for(int i = 0; i < states-1; i++)
	{
		this->Pi[i] = -1;
		// this->values_sh[i] = 2*((double) rand() / RAND_MAX) - 1;
		// this->values_shd[i] = 2*((double) rand() / RAND_MAX) - 1;
    this->values_sh[i] = 0.0;
		this->values_shd[i] = 0.0;
	}
	this->values_sh[states-1] = -1.0;
	this->values_shd[states-1] = -1.0;

	// for(int i = 0; i < states; i++)
	// {
	// 	cout << values[i] << ", ";
	// }
	// cout << endl;
}


void MDP::value_iterations(double epsilon)
{
	bool converged = false;
	double max_diff;
	double maximum;
	int iter = 0;
	while (!converged)
	{
		max_diff = -1.0;

		for (int s = 0; s < states-1; s++)
		{
			maximum = -999999.0;

			double temp = 0.0;
			double temp1 = 0.0;
			for (int s1 = 0; s1 < states; s1++)
			{
				temp = temp + (P[s][s1]*values_sh[s1]);
				temp1 = temp1 + (2*P[s][s1]*Qstand[s1]);
				// if (s == 72 && P[s][s1]!=0.0){
				// 	cout<<"Value summed --"<<(2*P[s][s1]*Qstand[s1])<<"\n";
				// }
			}

			Qhit[s] = temp;
			Qdd[s] = temp1;
		}

		// cout << "Qhit" <<  endl;
		// for(int i = 0; i < states-1; i++)
		// {
		// 	cout << Qhit[i] << ", ";
		// }
		// cout <<  endl;
    //
		// cout << "Qdd" <<  endl;
		// for(int i = 0; i < states-1; i++)
		// {
		// 	cout << Qdd[i] << ", ";
		// }
		// cout <<  endl;


		cout << "Calculated Qhit and Qdd" << endl;
		// cout<<"Prevvalue - "<<values[244]<<"\n";
		// cout<<"254 - "<<values[254]<<P[244][254]<<"\n";
		// cout<<"264 - "<<values[264]<<P[244][264]<<"\n";
		// cout<<"274 - "<<values[274]<<P[244][274]<<"\n";
		// cout<<"104 - "<<values[104]<<P[244][104]<<"\n";
		// cout<<"114 - "<<values[114]<<P[244][114]<<"\n";
		// cout<<"124 - "<<values[124]<<P[244][124]<<"\n";
		// cout<<"134 - "<<values[134]<<P[244][134]<<"\n";
		// cout<<"144 - "<<values[114]<<P[244][144]<<"\n";
		// cout<<"154 - "<<values[154]<<P[244][154]<<"\n";
		// cout<<"194 - "<<values[194]<<P[244][194]<<"\n";

		// cout<<"Hit - "<<Qhit[244]<<"\n";
		// cout<<"DD - "<<Qdd[244]<<"\n";
		// cout<<"Stand - "<<Qstand[244]<<"\n";

		int count = 0;
		for (int s = 0; s < states-1; s++)
		{
			double max = Qstand[s];
			Pi[s] = 0;
			double prev_sh = values_sh[s];
			double prev_shd = values_shd[s];

			if (Qhit[s] > max)
			{
				max = Qhit[s];
				Pi[s] = 1;
			}
			if (Qdd[s] > max)
			{
				max = Qdd[s];
				Pi[s] = 2;
			}

			values_shd[s] = max;

			double max1 = Qstand[s];
			if (Qhit[s] > max1)
				max1 = Qhit[s];

			values_sh[s] = max1;

			double diff = abs(values_sh[s] - prev_sh);
			if (diff > epsilon)
				count++;

			if (diff > max_diff)
			{
				max_diff = diff;
			}

			diff = abs(values_shd[s] - prev_shd);
			if (diff > epsilon)
				count++;

			if (diff > max_diff)
			{
				max_diff = diff;
			}
		}
		// cout<<values[244]<<"\n";

		// for(int i = 0; i < states; i++)
		// {
		// 	cout << values[i] << ", ";
		// }
		// cout << endl;

		iter++;
		cout<<"iteration :"<<iter<<endl;
		cout<<"Count>epsilon :"<<count<<endl;

		if (max_diff < epsilon) {
			converged = true;
		}
	}
}

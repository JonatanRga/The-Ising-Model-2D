#include <iostream>
#include <math.h>
#include <cmath>
#include <iomanip>		// for proper show Ising Matrix
#include <cstdlib>
#include <fstream>
#include <random>		

using namespace std;

const int LatSize = 10;

mt19937 generator(125);
uniform_real_distribution<double> dis(0.00, 1.00);

//gsl_rng_uniform(r)     decided to use in build to random libary MersneTwister generator instead GSL LIB

class cSpin
{
public:
	int iMatrix[LatSize][LatSize];

	cSpin() : M(0)			// 2 dim main randomly choosen array of Ising model  IsingMatrix
	{
		for (int i = 0; i < LatSize; i++)			 // loop for table hight
			for (int j = 0; j < LatSize; j++)        // loop for the table lenght
			{
				int U = (int)(dis(generator) + 0.5);	//choose radomly 0 or 1
				if (U < 1)
					U = -1;								//then 0 -> -1
				else
					U = 1;								// 1 stays 1 :)

				iMatrix[i][j] = U;
				M += U;
			}
	}


	int& operator()(int x, int y)	//for periodic bonduary conditions
	{
		if (x < 0)
			x += LatSize;

		else if (x >= LatSize)
			x -= LatSize;

		else if (y < 0)
			y += LatSize;

		else if (y >= LatSize)
			y -= LatSize;

		return iMatrix[x][y];
	}

	void display(int x, int y)
	{
		for (int i = 0; i < LatSize; i++)  // loop for table hight
		{
			for (int j = 0; j < LatSize; j++)  // loop for the table lenght
			{

				if ((i == x) && (j == y))
				{
					//system("Color A4");
					std::cout << setw(1) << "("<< setw(1) <<iMatrix[i][j] << ")" << setw(1);  // display the current element out of the array
					//system("Color A1");
				}
				else
					std::cout << setw(3) << iMatrix[i][j] << " ";  // display the current element out of the array
			}
			std::cout << endl;  // when the inner loop is done, go to a new line
		}
		
	}

	void LogMatrix(int logFlag1, float Temp, float Energy)
	{

		std::ofstream logM("logMatrix.txt");

		if (logFlag1 == 0)		// to be opened only first time
		{
			int counter1 = 0;
			logM << "DATA SPLITED USING FOLLOWING: ;" << endl;
			logM << "Spins Matrixes:" << endl;
			logM << "Number      Temperature     Energy" << endl;
			logFlag1 = 1;
		}

		for (int i = 0; i < LatSize; i++)
		{
			logM << endl << i << Temp << Energy << endl;
			for (int j = 0; j < LatSize; j++)
				logM << setw(2) << iMatrix[j][i] << " ";  // save the current matrix
		}

		logM.close();
	}

	int SpinSum()
	{
		return M;
	}

private:

	int M;

};





int main()
{
	const int J = 1;
	const int B = 1;

	//cSpin	spin;

	//srand(static_cast<unsigned> (time(0)));

	//spin.display();


	const int N = 20000;

	int M1;
	int LogFlag1 = 0;		//using for log creation
	int neighbour = 4;    // nearest neighbour count
	int flipCounter = 1; // number of flips count when energy is lower
	int flipCounter1 = 1;	//numer of flips when energy is else
	int currConfig = 1; //number of sweeps count
	const int Nflips = 2000;  // number of sweeps	520
	double avgM;          // average magnetic moment
	int M[N] = {};	  //will store current sum Spin - Z
	int sumM = 0;   //Cumulative Magnetic Moment starts from 0
	double sumE = 0.0;
	double E[N] = {};		//Energy
	double T[N] = {};		//Temperature
	double kB = 0.01 ;//8.6173303e-5;	//Boltzman constant in eV/K


	int x;		// coordinates on the ising matrix
	int y;

	float pFlip;	//probability of flip

	double delta_M = 0.0;
	double delta_neighbour = 0.0;
	double delta_E = 0.0;

	int bezpiecznik = 0;

	cSpin	spin;

	/*

	int E = -J * neighbour - B * M;   //Total Energy

	int delta_M = -2 * spin(x, y);	//change in Magnetic Moment

	int delta_neighbour = delta_M * (spin(x - 1, y) + spin(x + 1, y) + spin(x, y - 1) + spin(x, y + 1));

	int delta_E = -J * delta_neighbour - B * delta_M; //Change in total energy

	int M_Sum = 0;   //Cumulative Magnetic Moment starts from 0


	

	do       //	we are choosing state answering to temperature range in while condition
	{
		cSpin	spin;									// we are choosing initial state
		//M[0] = spin.SpinSum();								//statistic sum in Ising model called also Z
		int sSum = 0;

		for (int i = 0; i < LatSize; i++)
		{
			for (int j = 0; j < LatSize; j++)
			{
				sSum += spin(i, j);
				
				//delta_M += -2 * spin(x, y);    //change in Magnetic Moment	2->0.5
				//delta_neighbour = delta_M * (spin(x - 1, y) + spin(x + 1, y) + spin(x, y - 1) + spin(x, y + 1));

			}

		}

		std::cout << endl << sSum << endl;

		delta_E = exp(sSum / kB);// - B * M[0]; //Change in total energy
		//E[0] = delta_E;
		//T[0] = (2 * E[0]) / (neighbour *kB);
		//E = (-1) * ((double)neighbour) / (2 - kB * M);   //Total Energy

		std::cout << "    E: " << delta_E << "        T: " << T[0] << endl;

		bezpiecznik++;		// in case temerature won't be found


	} while (!(((T[0] > 200) && (T[0] < 300)) || (bezpiecznik > 10)));
	*/


	std::cout << endl << "-------------------------------------" << endl;
	std::cout << endl << "Choosen state to following Energy and Temperature:" << endl;
	std::cout << "    E: " << E[0] << "        T: " << T[0] << endl;
	std::cout << endl << "corresponding matrix to this state:" << endl;
	spin.display(11, 11);
	std::cout << endl << "cooling down in progress..." << endl;
	M[0] += delta_M;
	std::cout << endl << "-------------------------------------" << endl;


	//E = (-1) * ((double)neighbour) / (2 - kB * M);   //Total Energy
	//T = (2 * delta_E) / (delta_neighbour *kB);

	for (int s = 0; s < Nflips;)    //Loop for Nflips of sweeps
	{										//22 flips??!!
														//generate_random_coordinate_location
			y = (int)(dis(generator) * (LatSize + 1));       //Randomly choose x-coordinate add1
			x = (int)(dis(generator) * (LatSize + 1));        //Randomly choose y-coordinate add1
			spin(x, y) *= (-1);
			delta_M += -2.0 * spin(x, y);    //change in Magnetic Moment

			cout << "x: " << x << "    y: " << y << endl;
			std::cout << "s0: " << spin(x, y) << " sL: " << spin(x - 1, y) << " sR: " << spin(x + 1, y) << " sU: " << spin(x, y - 1) << " sD: " << spin(x, y + 1) << endl;


			M[flipCounter] = spin.SpinSum();
			delta_neighbour = delta_M * (spin(x - 1, y) + spin(x + 1, y) + spin(x, y - 1) + spin(x, y + 1));
			delta_E = exp(delta_neighbour / kB);														//Change in total energy
			
		
			//delta_E = (-1.0) * delta_neighbour /( 2.0 - kB * delta_M); //Change in total energy
																 //E = (-1) * neighbour / (2 - kB * M);   //Total Energy

			if ( (dis(generator) < (double)exp(-kB * delta_E)) ) //delta_E < 0.0) ||
			{
				//spin(x, y) *= -1;

				M[flipCounter] = spin.SpinSum() + delta_M;											// New Magnetization energy of flipped configuration
				sumE += E[flipCounter] = -delta_E;												// New total energy of flipped configuration
				T[flipCounter] = (2.0 * (-sumE)) / (neighbour *kB);
				sumM += M[flipCounter];												    //Summing up M to find sumM


																						//currConfig = (double)flipCounter / (LatSize * LatSize);	//Total number of sweeps till this point
																						//currConfig = (double)flipCounter / (LatSize * LatSize);	    //Total number of sweeps till this point

				std::cout << (s + 1) << ".-----------Success!------------------" << endl;
				std::cout << "E: " << sumE <<"    delta E: "<< delta_E<< "        T: " << T[flipCounter] << "        M: " << M[flipCounter] << endl;
				std::cout << "x: " << y << "    y: " << x << endl;
				spin.display(x, y);
				std::cout << "----------------------------------------" << endl;

				//if (currConfig == 0)
				//	currConfig = 1;
				//avgM = (sumM) / currConfig;									//Average M
				// Total number of sites chosen up till this point
				flipCounter++;
				delta_E = 0.0;
				s++;
			}
			else
			{
				spin(x, y) *= (-1);
				//cSpin	spin;
				//currConfig = (double)flipCounter / (LatSize * LatSize);     //Number of Sweeps
				flipCounter1++;   //Number of Flips
			}
		
		//spin.LogMatrix(logFlag1,  T, E);


		//std::cout << "The Magnetic Moment after " << currConfig << " sweeps = " << M << endl;
		//std::cout << "The Cumulative Magnetic Moment after " << currConfig << " sweeps = " << sumM << endl;
		//std::cout << "The Average Magnetic Moment after " << currConfig << "sweeps = " << avgM << endl;

	}

	/*
	std::cout << "The site chosen is at (" << x << ", " << y << ") = " << spin(x, y) << endl;
	std::cout << "The left neighbour has spin = " << spin(x-1, y) << endl;
	std::cout << "The right neighbour has spin = " << spin(x+1, y) << endl;
	std::cout << "The above neighbour has spin = " << spin(x, y-1) << endl;
	std::cout << "The below neighbour has spin = " << spin(x, y+1) << endl;
	std::cout << "The change in neighbour energy = " << delta_neighbour << endl;
	std::cout << "The change in total energy dE =" << delta_E << endl;
	*/
	std::cout << "flips summary: " << (flipCounter + flipCounter1) << endl;
	std::cout << "The number of flips E<0: " << flipCounter << "   The number of flips E>=0: " << flipCounter1 << endl;
	//std::cout << "The number of sweeps = " << currConfig << endl;
	//std::cout << "The final Magnetic Moment = " << M << endl;
	//std::cout << "The final average Magnetic Moment = " << avgM << endl;
	
	system("pause");
	return 0;
}




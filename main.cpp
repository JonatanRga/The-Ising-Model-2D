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

	int size = LatSize * LatSize;

	int M1;
	int LogFlag1 = 0;		//using for log creation
	int neighbour = 4;    // nearest neighbour count
	int flipCounter = 1; // number of flips count when energy is lower
	int magnetCounter = 0;	//numer of flips when energy is else
	int currConfig = 1; //number of sweeps count
	const int N = 200000;  // number of sweeps	520
	double avgM;          // average magnetic moment
	double M = 0.;	  //will store current sum Spin - Z
	int sumM = 0;   //Cumulative Magnetic Moment starts from 0
	double sumE = 0.0;
	
	double kB = 0.01 ;//8.6173303e-5;	//Boltzman constant in eV/K
	

	int x;		// coordinates on the ising matrix
	int y;

	float pFlip;	//probability of flip

	double delta_neighbour = 0.0;
	double delta_E = 0.0;
	double E = 0.;

	cSpin	spin;					// set all spins radomly

	double T = 1.5;		// Temperature - critical condtion 

	do
	{					// This will be BIG loop			
		
		for (int k = 0; k < N; k++)
		{
			for (int i = 0; i < LatSize; i++)
				for (int j = 0; j < LatSize; j++)
				{
					// 1. Generate random coordinate location
					y = (int)(dis(generator) * (LatSize + 1));       //Randomly choose x-coordinate
					x = (int)(dis(generator) * (LatSize + 1));       //Randomly choose y-coordinate
			// 2. Than flip it
					spin(x, y) *= (-1);
					// 3. Than we do this counting:
					delta_E = -2.0 * spin(x, y)* (spin(x - 1, y) + spin(x + 1, y) + spin(x, y - 1) + spin(x, y + 1));

					//cout << "x: " << x << "    y: " << y << endl;
					//std::cout << "s0: " << spin(x, y) << " sL: " << spin(x - 1, y) << " sR: " << spin(x + 1, y) << " sU: " << spin(x, y - 1) << " sD: " << spin(x, y + 1) << endl;
					//delta_M += ;    //change in Magnetic Moment
					//M[flipCounter] = spin.SpinSum();

			// 4. Than we chceck condition for acceptance, first: Is it smaller than previous energy?, second: Is it smaller than some random? - just to not stuck overflow :)
					if (     (delta_E < 0) || (dis(generator) < exp(delta_E / T))       )
					{
						//4a. Condtion was succefully accepted, so we accept this state to be our new state
						spin(x, y) *= -1;							//get back previous state
						flipCounter++;								//counts successfoul events
					}
				}
			// 5. If we did enough step to let the Ising matrix normalize a bit, than we can messure magnetization also


			if ((flipCounter % 1000 == 0) && (flipCounter > 30000))
			{
				M += (double)(abs(spin.SpinSum()) / size);
			}
		}
		cout <<"Energy = " << delta_E << "     Temperature = " << T << "    <Magnetization> = " << M << endl;
		M = 0.;
		T += 0.1;	// Temperature step
	} while (T < 3.55);



			
	//spin.LogMatrix(logFlag1,  T, E);


	//std::cout << "The Magnetic Moment after " << currConfig << " sweeps = " << M << endl;
	//std::cout << "The Cumulative Magnetic Moment after " << currConfig << " sweeps = " << sumM << endl;
	//std::cout << "The Average Magnetic Moment after " << currConfig << "sweeps = " << avgM << endl;

	

	/*
	std::cout << "The site chosen is at (" << x << ", " << y << ") = " << spin(x, y) << endl;
	std::cout << "The left neighbour has spin = " << spin(x-1, y) << endl;
	std::cout << "The right neighbour has spin = " << spin(x+1, y) << endl;
	std::cout << "The above neighbour has spin = " << spin(x, y-1) << endl;
	std::cout << "The below neighbour has spin = " << spin(x, y+1) << endl;
	std::cout << "The change in neighbour energy = " << delta_neighbour << endl;
	std::cout << "The change in total energy dE =" << delta_E << endl;
	*/
	
	//std::cout << "The number of sweeps = " << currConfig << endl;
	//std::cout << "The final Magnetic Moment = " << M << endl;
	//std::cout << "The final average Magnetic Moment = " << avgM << endl;
	
	system("pause");
	return 0;
}




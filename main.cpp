#include <iostream>
#include <math.h>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <random>		



using namespace std;

const int LatSize = 10;

std::mt19937 generator(124);
std::uniform_real_distribution<double> dis(0.00, 1.00);


//gsl_rng_uniform(r)
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

				iMatrix[j][i] = U;
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

	void display()
	{
		for (int i = 0; i < LatSize; i++)  // loop for table hight
		{
			for (int j = 0; j < LatSize; j++)  // loop for the table lenght
				cout << setw(2) << iMatrix[j][i] << " ";  // display the current element out of the array

			cout << endl;  // when the inner loop is done, go to a new line
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

	cSpin	spin;

	//srand(static_cast<unsigned> (time(0)));

	spin.display();



	int aMatrix[LatSize][LatSize];
	int M1;

	int neighbour = 4;    // nearest neighbour count
	int flipCounter = 1; // number of flips count when energy is lower
	int flipCounter1 = 1;	//numer of flips when energy is else
	int currConfig = 0; //number of sweeps count
	int Nflips = 20;  // number of sweeps	520
	double avgM;          // average magnetic moment
	unsigned int M = 0;	  //will store current sum Spin - Z
	int sumM = 0;   //Cumulative Magnetic Moment starts from 0

	double E;		//Energy
	double T = 200;
	double kB = 1.0;		//8.6173303e-5;	//Boltzman constant in eV/K

	int x;		// coordinates on the ising matrix
	int y;

	float pFlip;	//probability of flip


	int LogFlag1 = 0;


	for (int s = 0; s < Nflips; ++s)    //Loop for Nflips of sweeps
	{

		for (int t = 0; t < LatSize * LatSize; ++t)   //Loop for 1 Sweep
		{
			M = spin.SpinSum();

			x = (int)(dis(generator) * (LatSize + 1));       //Randomly choose x-coordinate add1
			y = (int)(dis(generator) * (LatSize + 1));        //Randomly choose y-coordinate add1

			int delta_M = -2 * spin(x, y);    //change in Magnetic Moment
			int delta_neighbour = delta_M * (spin(x - 1, y) + spin(x + 1, y) + spin(x, y - 1) + spin(x, y + 1));
			int delta_E = (-1) * delta_neighbour / 2 - kB * delta_M; //Change in total energy

			E = (-1) * neighbour / (2 - kB * M);   //Total Energy
			T = (2 * delta_E) / (delta_neighbour *kB);

			float pFlip = exp((-1)*delta_E);

			if ((delta_E < 0) || (dis(generator) < (double)exp(-delta_E)))
			{
				//spin(x, y) *= -1;

				M += delta_M;												// New Magnetization energy of flipped configuration
				E += delta_E;												// New total energy of flipped configuration

				flipCounter++;												// Total number of sites chosen up till this point
			}
			else
			{
				cSpin flip;
				//currConfig = (double)flipCounter / (LatSize * LatSize);     //Number of Sweeps

				flipCounter1++;   //Number of Flips

			}

			currConfig = (double)flipCounter / (LatSize * LatSize);	//Total number of sweeps till this point
			if (currConfig == 0)
				currConfig = 1;
		}


		sumM += M;                       //Summing up M to find sumM
		avgM = (sumM) / currConfig;    //Average M

		cout << s << "    E: " << E << "        T: " << T << endl;

		spin.display();

		//spin.LogMatrix(logFlag1,  Temp, E);


		//cout << "The Magnetic Moment after " << currConfig << " sweeps = " << M << endl;
		//cout << "The Cumulative Magnetic Moment after " << currConfig << " sweeps = " << sumM << endl;
		//cout << "The Average Magnetic Moment after " << currConfig << "sweeps = " << avgM << endl;

	}

	/*
	cout << "The site chosen is at (" << x << ", " << y << ") = " << spin(x, y) << endl;
	cout << "The left neighbour has spin = " << spin(x-1, y) << endl;
	cout << "The right neighbour has spin = " << spin(x+1, y) << endl;
	cout << "The above neighbour has spin = " << spin(x, y-1) << endl;
	cout << "The below neighbour has spin = " << spin(x, y+1) << endl;
	cout << "The change in neighbour energy = " << delta_neighbour << endl;
	cout << "The change in total energy dE =" << delta_E << endl;
	*/
	cout << "The number of flips E<0: " << flipCounter << "   The number of flips E>=0: " << flipCounter1 << endl;
	cout << "The number of sweeps = " << currConfig << endl;
	cout << "The final Magnetic Moment = " << M << endl;
	cout << "The final average Magnetic Moment = " << avgM << endl;

	system("pause");
	return 0;
}
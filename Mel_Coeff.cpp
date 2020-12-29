#include <math.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include<vector> 
#include <complex.h>
#include "Mel_Coeff.h"

#include<iostream>
#include<fstream>

#include <chrono>
using namespace std::chrono; 

//#define ARMA_DONT_USE_CXX11
#include <armadillo>

#ifdef _DEBUG
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
// Replace _NORMAL_BLOCK with _CLIENT_BLOCK if you want the
// allocations to be of _CLIENT_BLOCK type
#else
#define DBG_NEW new
#endif

#define PI 3.141592653589793

int main()
{

	MEL_C::Mel_coeff mc;

	/*
	// open file stream
	std::ofstream file;

	file.open("Sin_10_Hz.txt");

	std::vector<double> audio_i(16000);
	for (int kk = 0; kk < 16000; kk++)
	{

		audio_i.at(kk) = sin(2*PI*10*(double)kk/16000);

		file << audio_i.at(kk);

	}

	file.close();
	*/

	int length_sine = 16000;
	std::vector<double> audio_i(length_sine);
	for (int kk = 0; kk < length_sine; kk++)
	{

		audio_i.at(kk) = sin(2 * PI * 10 * (double)kk / 16000);

	}

	int num_coeff = 13;
	int sf = 16000;
	int size_win = 400;
	int size_overlap = 320;
	bool log_energy = true;

	//auto start = high_resolution_clock::now(); 

	std::vector<std::vector<double>> matrix_out = mc.mel_coeff_output(audio_i, num_coeff, sf, size_win, size_overlap, log_energy);

	//auto stop = high_resolution_clock::now(); 

	//auto duration = duration_cast<microseconds>(stop - start); 

	//std::cout << "Execution time: " << duration.count()/1000 << " ms" << std::endl; 

/*
std::ofstream file;
file.open("Mel_Coeff.txt");
std::vector<double> audio_o(matrix_out.size());

for (int ll = 0; ll < matrix_out[0].size(); ll++)
{

        for (int kk = 0; kk < matrix_out.size(); kk++)
        {

        audio_o.at(kk) = matrix_out[kk][ll];

        file << audio_o.at(kk) << "\t";

        }

        file << "\n";

}
        file.close();

*/
}


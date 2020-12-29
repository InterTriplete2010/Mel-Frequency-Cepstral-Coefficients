#include <iostream>
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include<vector> 
#include <complex.h>
#include <math.h>
#include <numeric>
#include <algorithm>
#include<iostream>
#include<fstream>
#include "Mel_Coeff.h"

//#define ARMA_DONT_USE_CXX11
#include <armadillo>

#ifdef _DEBUG
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
// Replace _NORMAL_BLOCK with _CLIENT_BLOCK if you want the
// allocations to be of _CLIENT_BLOCK type
#else
#define DBG_NEW new
#endif

int num_iter;	//Number of columns of the MEL-Coeff column

//Number of the default band edges
int num_bands = 42;

//Array for the band edges
std::vector<double> band_edge(num_bands);

//Global variable for the FFT
int n_length;
std::complex<double> complex_imag_neg(0.0, -1.0);

#define PI 3.141592653589793

using namespace MEL_C;

//Initalize the 42 default values of the band edge
void Mel_coeff::initialize_band_edge()
{

		band_edge.at(0) = 133.333333333333;
		band_edge.at(1) = 200;
		band_edge.at(2) = 266.666666666667;
		band_edge.at(3) = 333.333333333333;
		band_edge.at(4) = 400; 
		band_edge.at(5) = 466.666666666667;
		band_edge.at(6) = 533.333333333333;
		band_edge.at(7) = 600;
		band_edge.at(8) = 666.666666666667;
		band_edge.at(9) = 733.333333333333;
		band_edge.at(10) = 800;
		band_edge.at(11) = 866.666666666667;
		band_edge.at(12) = 933.333333333333;
		band_edge.at(13) = 999.758946666667;
		band_edge.at(14) = 1070.91209082862;
		band_edge.at(15) = 1147.12922560652;
		band_edge.at(16) = 1228.77075673170;
		band_edge.at(17) = 1316.22274011952;
		band_edge.at(18) = 1409.89870740065;
		band_edge.at(19) = 1510.24162137597;
		band_edge.at(20) = 1617.72597064178;
		band_edge.at(21) = 1732.86001329015;
		band_edge.at(22) = 1856.18818029401;
		band_edge.at(23) = 1988.29364994199;
		band_edge.at(24) = 2129.80110549646;
		band_edge.at(25) = 2281.37968911497;
		band_edge.at(26) = 2443.74616600319; 
		band_edge.at(27) = 2617.66831376149;
		band_edge.at(28) = 2803.96855295239;
		band_edge.at(29) = 3003.52783605657;
		band_edge.at(30) = 3217.28981320707;
		band_edge.at(31) = 3446.26529439996;
		band_edge.at(32) = 3691.53702928199;
		band_edge.at(33) = 3954.26482711710;
		band_edge.at(34) = 4235.69104114247;
		band_edge.at(35) = 4537.14644324789;
		band_edge.at(36) = 4860.05651675778;
		band_edge.at(37) = 5205.94819707239;
		band_edge.at(38) = 5576.45709204249;
		band_edge.at(39) = 5973.33521622028;
		band_edge.at(40) = 6398.45927555924;
		band_edge.at(41) = 6853.83954173857;
}

//Convert the audio signal into a matrix
std::vector<std::vector<double>> Mel_coeff::auto_internal_buffer(std::vector<double> audio_i,int wl, int hl)
{

	int track_index;
	int track_cycle;

	num_iter = std::floor((audio_i.size() - wl) / hl) + 1;	//Number of columns

	//The new matrix will be (wl x num_iter)
	std::vector<double> temp_v;

	for (int ff = 0; ff < num_iter; ff++)
	{

		temp_v.push_back(0);

	}

	std::vector<std::vector<double>> matrix_audio;
	for (int kk = 0; kk < wl; kk++)
	{

		matrix_audio.push_back(temp_v);

	}

	//Save the data into the new matrix
	for (int kk = 0; kk < num_iter; kk++)
	{
		
		track_cycle = hl * (kk);
		track_index = 0;
		
		for (int ll = track_cycle; ll < wl + track_cycle; ll++)
		{

			matrix_audio[track_index][kk] = audio_i.at(ll);
			track_index++;

		}

	}

	return matrix_audio;

}

//Calcuate the ln(energy)
std::vector<double> Mel_coeff::calc_log_energy(std::vector<std::vector<double>> audio_i, int window_length)
{
	//Create the vector where to save the energy
	std::vector<double> store_log_energy(audio_i[0].size());

	//Calculate the energy 
	for (int kk = 0; kk < audio_i[0].size(); kk++)
	{

		for (int ll = 0; ll < window_length; ll++)
		{

			store_log_energy[kk] += std::pow(audio_i[ll][kk],2);

		}

			store_log_energy[kk] = std::log(store_log_energy[kk]);

	}

	return store_log_energy;

}

//Create a Hamming window to be used with the FFT. The length corresponds to the number of rows of the audio matrix
std::vector<double> Mel_coeff::hamm_w(int length_w)
{

	double coeff_I = 0.54;
	double coeff_II = 0.46;

	std::vector<double> hamm_window(length_w);

	for (int kk = 0; kk < length_w; kk++)
	{

		hamm_window[kk] = coeff_I - coeff_II * cos(2 * PI * kk / length_w);

	}

	return hamm_window;

}

//Transpose the matrix
std::vector<std::vector<double>> Mel_coeff::transpose_matrix(std::vector<std::vector<double>> matrix_in)
{

	std::vector<std::vector<double>> matrix_out;

	std::vector<double> temp_v;

	for (int ff = 0; ff < matrix_in.size(); ff++)
	{

		temp_v.push_back(0);

	}

	for (int kk = 0; kk < matrix_in[0].size(); kk++)
	{

		matrix_out.push_back(temp_v);

	}

	for (int kk = 0; kk < matrix_in[0].size(); kk++)
	{

		for (int ll = 0; ll < matrix_in.size(); ll++)
		{

			matrix_out[kk][ll] = matrix_in[ll][kk];

		}

	}

	return matrix_out;

}

//Multiply the input by the Hamming Window
arma::mat Mel_coeff::Ham_Win(std::vector<std::vector<double>> audio_i, std::vector<double> hamm_window)
{

	arma::mat output_HW(audio_i.size(), audio_i[0].size());

	//Multiply the input signal by the Hamming window
	for (int ll = 0; ll < audio_i[0].size(); ll++)
	{

		for (int kk = 0; kk < audio_i.size(); kk++)
		{

			output_HW(kk,ll) = audio_i[kk][ll] * hamm_window[kk];
			
		}

	}

	return output_HW;

}

//Calculate the FFT
arma::mat Mel_coeff::FFT_Win(arma::mat audio_i)
{
	
	arma::mat fft_output = arma::abs(arma::fft(audio_i));
	
	return fft_output;

}

//Design the filter banks
std::vector<std::vector<double>> Mel_coeff::filter_bank(int sf, std::vector<double> band_edges_in, int window_length)
{

	//	Determine the number of valid bands
	int	validNumBands = band_edges_in.size() - 2;

	//Build the vector where to save the filter banks
	std::vector<std::vector<double>> filter_bank_out;

	std::vector<double> temp_v;

	for (int ff = 0; ff < validNumBands; ff++)
	{

		temp_v.push_back(0);

	}

	for (int kk = 0; kk < window_length; kk++)
	{

		filter_bank_out.push_back(temp_v);

	}

	//Center frequencies
	std::vector<double> center_freq(band_edges_in.size() - 2);
	for (int kk = 1; kk < band_edges_in.size() - 1; kk++)
	{

		center_freq.at(kk - 1) = band_edges_in.at(kk);

	}

	//	The following alorithm is specified by the documentation of Slaney's Auditory Toolbox
	std::vector<double> lin_freq(window_length);
	for (int kk = 0; kk < window_length; kk++)
	{

		lin_freq.at(kk) = ((double)kk*(double)sf) / ((double)window_length);

	}

	//Code lines 77 to 87 in Matlab
	std::vector<double> p_vector(band_edges_in.size());

	for (int kk = 0; kk < band_edges_in.size(); kk++)
	{

		for (int ll = 0; ll < lin_freq.size(); ll++)
		{

			if (lin_freq.at(ll) > band_edges_in.at(kk))
			{
			
				p_vector.at(kk) = ll;
				break;

			}

		}


	}

	//	Code 89 to 106 in Matlab 
	//	Create triangular filters for each band
	std::vector<double> bw_vector(band_edges_in.size() - 1);
	for (int kk = 0; kk < bw_vector.size(); kk++)
	{

		bw_vector.at(kk) = band_edges_in.at(kk + 1) - band_edges_in.at(kk);

	}

	for (int kk = 0; kk < validNumBands; kk++)
	{

		//	Rising side of triangle
		for (int ll = p_vector.at(kk); ll < p_vector.at(kk + 1); ll++)
		{

			filter_bank_out[ll][kk] = (lin_freq.at(ll) - band_edges_in.at(kk)) / bw_vector.at(kk);

		}

		// Falling side of triangle
		for (int ll = p_vector.at(kk + 1); ll < p_vector.at(kk + 2); ll++)
		{

			filter_bank_out[ll][kk] = (band_edges_in.at(kk + 2) - lin_freq.at(ll)) / bw_vector.at(kk + 1);

		}

	}

	//	Apply normalization
	//	Weight by bandwidth
	std::vector<double> filter_band_width(validNumBands);
	for (int kk = 2; kk < band_edges_in.size(); kk++)
	{

		filter_band_width.at(kk - 2) = band_edges_in.at(kk) - band_edges_in.at(kk - 2);

	}

	std::vector<double> weight_band(validNumBands);
	for (int kk = 0; kk < validNumBands; kk++)
	{

		weight_band.at(kk) = 2/ filter_band_width.at(kk);

	}

	for (int kk = 0; kk < validNumBands; kk++)
	{

		for (int ll = 0; ll < filter_bank_out.size(); ll++)
		{

			filter_bank_out[ll][kk] = filter_bank_out[ll][kk]* weight_band.at(kk);

		}
	}

	return filter_bank_out;
}

//Calculate the cepstral coefficients
std::vector<std::vector<double>> Mel_coeff::cepstral_coefficients(arma::mat abs_fft, std::vector<std::vector<double>> filter_b, int num_coeff, int colms_mel_coeff, bool log_energy)
{

	//Matrix used to save the filter banks*abs(FFT)
	std::vector<std::vector<double>> filter_fft;

	std::vector<double> temp_v;

	for (int ff = 0; ff < colms_mel_coeff; ff++)
	{

		temp_v.push_back(0);

	}

	for (int kk = 0; kk < filter_b[0].size(); kk++)
	{

		filter_fft.push_back(temp_v);

	}

	//Multiply the filter banks by the audio
	for (int kk = 0; kk < filter_b[0].size(); kk++)
	{

		for (int hh = 0; hh < abs_fft.n_rows; hh++)
		{

			for (int ll = 0; ll < filter_b.size(); ll++)
			{

				filter_fft[kk][hh] += filter_b[ll][kk] * abs_fft(hh,ll);
			
			}

		}

		//Take the log10
		for (int hh = 0; hh < abs_fft.n_rows; hh++)
		{

			filter_fft[kk][hh] = std::log10(filter_fft[kk][hh]);

		}

	}

	//Create a DCTmatrix
	std::vector<std::vector<double>> DCTmatrix;

	std::vector<double> temp_vv;

	for (int ff = 0; ff < filter_fft.size(); ff++)
	{

		temp_vv.push_back(0);

	}

	for (int kk = 0; kk < num_coeff; kk++)
	{

		DCTmatrix.push_back(temp_vv);

	}

	double A = std::sqrt(1 / (double)filter_fft.size());
	double B = std::sqrt(2 / (double)filter_fft.size());
	double C = 2 * (double)filter_fft.size();
	
	for (double kk = 0; kk < filter_fft.size(); kk++)
	{
		for (double nn = 0; nn < num_coeff; nn++)
		{
			if (nn == 0)
			{

				DCTmatrix[nn][kk] = A * std::cos(PI * nn * (2 * (kk + 1) - 1) / C);

			}

			else
			{

				DCTmatrix[nn][kk] = B * std::cos(PI * nn * (2 * (kk + 1) - 1) / C);
			
			}
		}
	}

	//Calculate the coefficients
	std::vector<std::vector<double>> mel_coeff_out;

	std::vector<double> temp_vvv;

	int mel_coeff_start = 0;

	if (log_energy)
	{

		mel_coeff_start++;

	}

	for (int ff = 0; ff < colms_mel_coeff; ff++)
	{

		temp_vvv.push_back(0);

	}

	for (int kk = 0; kk < num_coeff + mel_coeff_start; kk++)
	{

		mel_coeff_out.push_back(temp_vvv);

	}

	//Final step to get the mel coefficients
	for (int kk = 0; kk < num_coeff; kk++)
	{

		for (int hh = 0; hh < colms_mel_coeff; hh++)
		{

			for (int ll = 0; ll < DCTmatrix[0].size(); ll++)
			{

				mel_coeff_out[kk + mel_coeff_start][hh] += DCTmatrix[kk][ll] * filter_fft[ll][hh];

			}

		}

	}

	return mel_coeff_out;

}

//Calcuate the MEL-COEFF
std::vector<std::vector<double>> Mel_coeff::mel_coeff_output(std::vector<double> audio_i, int num_coeff, int sf, int window_length, int overlap_length, bool log_energy)
{

	int hop_length = window_length - overlap_length;
	
	//Convert the audio input into a matrix
	std::vector<std::vector<double>> matrix_audio = auto_internal_buffer(audio_i, window_length, hop_length);

	//Calculate the logarithmic energy, if the flag is true
	std::vector<double> energy_audio_i;
	if (log_energy)
	{

		energy_audio_i = calc_log_energy(matrix_audio, window_length);

	}

	//Calculate the Hamming window
	std::vector<double> hamm_window = hamm_w(window_length);

	//Multiply the input by the Hamming Window
	arma::mat matrix_audio_HW = Ham_Win(matrix_audio, hamm_window);

	//Calculate the FFT
	arma::mat abs_fft_output = FFT_Win(matrix_audio_HW);//FFT_Win(matrix_audio_HW[kk]);
	
	//Initialize the default Band Edges
	initialize_band_edge();

	//Design the filter banks
	std::vector<std::vector<double>> filter_banks = filter_bank(sf, band_edge, window_length);

	abs_fft_output = abs_fft_output.t();

	//Calculate the Mel Coefficients
	std::vector<std::vector<double>> cep_coeff = cepstral_coefficients(abs_fft_output,filter_banks,num_coeff, matrix_audio[0].size(), log_energy);

	//Check if the energy of the signal needs to be added to the coefficients
	if (log_energy)
	{

		for (int kk = 0; kk < cep_coeff[0].size(); kk++)
		{

			cep_coeff[0][kk] = energy_audio_i[kk];

		}
	}

	//Transpose the mel coeff matrix to  be consistent with Matlab's output
	cep_coeff = transpose_matrix(cep_coeff);

	return cep_coeff;

}

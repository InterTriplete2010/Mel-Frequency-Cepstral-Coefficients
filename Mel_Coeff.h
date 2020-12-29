#pragma once
#include <stdio.h> 
#include <complex.h>

//#define ARMA_DONT_USE_CXX11
#define ARMA_DONT_USE_CXX11_MUTEX
#include <armadillo>

#ifndef MEL_COEFF_H
#define MEL_COEFF_H

#ifdef __cplusplus
extern "C" {  // only need to export C interface if
              // used by C++ source code
#endif

    namespace MEL_C
    {

        class Mel_coeff

        {


        private:

            //Populate the array with Band Edge frequencies
            void initialize_band_edge();

            //Convert the audio signal into a matrix
            std::vector<std::vector<double>> auto_internal_buffer(std::vector<double>, int, int);

            //Calcuate the ln(energy)
            std::vector<double> calc_log_energy(std::vector<std::vector<double>>, int);

            //Create a Hamming window to be used with the FFT. The length corresponds to the number of rows of the audio matrix
            std::vector<double> hamm_w(int);

            //Transpose the matrix
            std::vector<std::vector<double>> transpose_matrix(std::vector<std::vector<double>>);

            //Multiply the input by the Hamming Window
            arma::mat Ham_Win(std::vector<std::vector<double>>, std::vector<double>);

            //Calculate the FFT
            arma::mat FFT_Win(arma::mat);

            //Design the filter banks
            std::vector<std::vector<double>> filter_bank(int, std::vector<double>, int);
             
            //Calculate the cepstral coefficients
            std::vector<std::vector<double>> cepstral_coefficients(arma::mat, std::vector<std::vector<double>>,int,int, bool);

        public:

            //calcuate the mel-coefficients
            std::vector<std::vector<double>> mel_coeff_output(std::vector<double>, int, int, int, int, bool);

        };

    }

#endif

#ifdef __cplusplus

}

#endif

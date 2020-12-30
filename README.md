# Mel-Frequency-Cepstral-Coefficients
Mel Frequency Cepstral Coefficients

This code calculates the Mel-Frequency-Cepstral-Coefficients by following the same steps as in Matlab (function: mfcc). The code uses the default 40-band filter bank that spans approximately 133 Hz to 6864 Hz, as reported in Matlab.

The function mel_coeff_output takes 6 arguments:

    A) The input vector
    B) The number of coefficients that need to be calculated
    C) The sampling frequency
    D) The WindowLength 
    E) The OverlapLength
    F) A boolean variable that allows the user to decided whether or not to calcuate and save the Log(Energy) of the signal
    
In the example reported in Mel_coeff, we have a sine waveform of length 16000, oscillating at 10 Hz and with sampling frequency of 16000 Hz. The number of coefficients chosen is 13, the WindowLength is 400 samples, the OverlapLength is 320 and the boolean variable is set to true => the Log(Energy) of the signal will be calculated and saved in the output matrix "matrix_out". The "matrix_out" will have a final dimension of 196 x 14, where the first column represents the Log(Energy) of the signal and the remaining columns are the Mel Frequency Cepstral Coefficients.

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
        

Please, keep in mind that this code needs the library Armadillo to be compiled: http://arma.sourceforge.net/docs.html. 
You can compile it in Linux by typing the following command line: g++ -ggdb Mel_Coeff.cpp Mel_Coeff_H.o -larmadillo -o Mel_exe

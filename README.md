# Perceptual_Cycles_Travel
The behavioral raw data and codes provided here allow performing the main analyses (Phase Effect on Detection Performance; Behavioral Optimal Phase Shift Between Target Positions; Figure 4) from Fakche C. &amp; Dugué L., Perceptual cycles travel across retinotopic space (2023), Journal of Cognitive Neuroscience, with Matlab R2014b 64-bit. 

The Data folder contains 17 subfolders, one for each participant. In each participant folder, there are four matlab files:
CodeParticipant_result_matrix_4Hz.mat, CodeParticipant_result_matrix_6Hz.mat, CodeParticipant_result_matrix_8Hz.mat
and CodeParticipant_result_matrix_10Hz.mat, for the peripheral disk oscillating at 4, 6, 8 and 10 Hz, respectively. 
In each file, there is a Matlab variable named result_matrix_all, that is a 2-D matrix, with each row corresponding to one trial. 
The first column indicates the position of the target (1, 2, or 3; note that 1 corresponds to Position 3, and 3 corresponds to Position 1), 
the second column, the response of the participant (1: perceived, 0: not perceived), and the third column the phase bin of the oscillating disk 
at the moment the target is displayed (from phase bin 1 to 7).

WavesLocal_run1_DataFitting
- Step 1: Computation of the hit rate for each phase bin.
- Step 2: Average hit rate across participants, and 95% CI Computation. 
- Step 3: Fit data to a complex sine function composed of the induced frequency and the first harmonic.
- Step 4: Plot the real data and the fit. 

WavesLocal_run2_MonteCarlo
- Step 1: Create 50,000 datasets of hit rate by randomly assigning a performance value for each trial according to the average performance. 
- Step 2: Average the 50,000 datasets across participants. 
- Step 3: Fit 50,000 datasets to the complex sine function. 
- Step 4: Compute p-value on the fitted amplitude.  

WavesLocal_run3_PhaseShift
- Step 1: Compute the optimal phase on data averaged across participants. 
- Step 2: Compute the optimal phase for individual data. 
- Step 3: Rose plot. 
- Step 4: Compute Harrison-Kanji test on optimal phase in radians.
- Step 5: Compute absolute phase difference in radians for each participants.
- Step 6: Compute Watson-Williams test on absolute phase difference against zero.  

WavesLocal_run4_PropagationSpeed
- Step 1: Compute the phase shift in degrees. 
- Step 2: Estimate the propagation speed according to the phase shift (in degrees) and the cortical distance between targets. 

This code simulates one frame of a MIMO radar system with tunable number of targets. 

The system uses an analog beamformed phased array transmitter and fully digital array receiver. Each frame consists of N_beacon subframes, and each subframe contains N_chirp chirps. On each subframe, the transmitter modulates its spatial domain signal with a different compressive beacon. The beacon responses in each range-doppler bin are used to estimate angle of departure compressively. Angle of arrival is directly estimated since the receiver is fully digital (for now, receiver is single element).

In the current script, you can either run a single test which plots the range-doppler heatmaps after each target identification, or you may run in monte carlo mode to get histograms of range, doppler, and angle estimates. These are reported as normalizes to the nominal FFT grid size of the fast FFT, slow FFT, and spatial FFT, respectively.

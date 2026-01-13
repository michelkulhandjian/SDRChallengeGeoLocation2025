# How to runProcesses 
1. Open runProcesses.m
2. Enable logging in line 6 if needed.
3. It has four modes (e.g., GENDATASET, CROSSCORRELATE, LOCANALYSIS, and LOCALIZE)
4. The mode GENDATASET generates IQ dataset or load IQ dataset given by the parameters with input/output folders (e.g, input_path_dir and output_path_dir)
5. The mode LOCANALYSIS will provide geometric studies for the given parameters such the positions of receiving nodes and the transmitter. Currently it is not working
6. The Mode CROSSCORRELATE will load IQ from the input directory or generate dataset provided params.m file  (e.g., fs, position of receivers, Tx positions, modulations etc.) it will compare the cross-correlation TDOA vs true TDOA computed from the geometry. 
6. The mode LOCALIZE 
      a. load from IQ for the given parameter settings and correlate to obtain TDOAs - geolocate the transmitter
	  b. generate IQ and use correlation to obtain TDOAs  - geolocate the transmitter
	  c. compute true TDOAs  - geolocate the transmitter
	  d. read TDOAs from tdoa.txt  - geolocate the transmitter
	  e. read from cvs file provided in the input folder  - geolocate the transmitter
	  f. Can take LAT/LONG or Euclidean coordinates 
	  g. Provides plots for different SNR provided at params.m 
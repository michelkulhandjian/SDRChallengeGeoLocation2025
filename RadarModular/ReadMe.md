# Create SEC Output
```
runProcesses.m
```
Lines 9,10  -> input & out directory
Input formats -> BIN1-A : IQ samples (9.6K) .dat file.  Fixed time duration
->Specify the input and output directories
->Produces Sparse and Full sec output files (currently N = 6000)

- params.m contains some settable parameters

# Matlab files
1. Place matlab files in a folder (e.g., /space/radar/ , etc.). 

# To run runProcesses.m as GENDATASET
1. You need to specify outputPathDir path directory. If not specified then will take hard coded outputPathDir. 
2. You need to specify the runType = "GENDATASET". 
3. You need to specify the mode type. Please refer to generate_radar_dataset.m for mode explanation.
4. Example of Linux command:  /usr/local/MATLAB/R2023a/bin/matlab -nodisplay -nosplash -nodesktop -r "inputPathDir=''; outputPathDir='/space/dataset';runType='GENDATASET'; mode=3; run /space/radar/runProcesses(inputPathDir, outputPathDir, runType, mode);exit;"
5. It will generate dataset of .dat files given the parameters in params.m file and save it in outputPathDir directory.

# To run runProcesses.m as DETECTION
1. You need to specify the runType = "DETECTION".
2. This can run in different modes. Please refer the mode types in processRadarDetection.m for mode explanation.
3. In case of mode = 2, point the IQ captured waveform dataset .dat  (e.g., /dataset/BIN1-A , etc.)
4. You need to specify your inputPathDir and outputPathDir path directories. If not specified, it will take 'pwd'. 
5. Example of Linux command:  /usr/local/MATLAB/R2023a/bin/matlab -nodisplay -nosplash -nodesktop -r "inputPathDir='/dataset/BIN1-A'; outputPathDir='/space';runType='DETECTION'; mode=2; run /space/radar/runProcesses(inputPathDir, outputPathDir, runType, mode);exit;"
6. The detection results will be saved in results.txt file in outputPathDir/results/. 

# To run runProcesses.m as ROC
1. You need to specify the runType = "ROC".
2. This can run in different modes. Please refer the mode types in roc_analysis.m for mode explanation.
3. In case of mode = 2, point the IQ captured waveform dataset .dat  (e.g., /dataset/BIN1-A , etc.)
4. You need to specify your inputPathDir and outputPathDir path directories. If not specified, it will take 'pwd'. 
5. Example of Linux command:  /usr/local/MATLAB/R2023a/bin/matlab -nodisplay -nosplash -nodesktop -r "inputPathDir='/dataset/BIN1-A'; outputPathDir='/space';runType='ROC'; mode=2; run /space/radar/runProcesses(inputPathDir, outputPathDir, runType, mode);exit;"
6. The ROC analysis results will be saved in res.dat file in outputPathDir/results/. 

# To run runProcesses.m as ESTIMATION
1. You need to specify the runType = "ESTIMATION".
2. This can run in different modes. Please refer the mode types in processRadarEstimation.m for mode explanation.
3. In case of mode = 4, point the IQ captured waveform dataset .dat  (e.g., /dataset/BIN1-A , etc.)
4. You need to specify your inputPathDir and outputPathDir path directories. If not specified, it will take 'pwd'. 
5. Example of Linux command:  /usr/local/MATLAB/R2023a/bin/matlab -nodisplay -nosplash -nodesktop -r "inputPathDir='/dataset/BIN1-A'; outputPathDir='/space';runType='ESTIMATION'; mode=2; run /space/radar/runProcesses(inputPathDir, outputPathDir, runType, mode);exit;"
6. The estimation results will be saved in results.txt file in outputPathDir/results/. 
6. If the radar is detected, then the .json file will be saved in the output folder specified in outputPathDir/results/ under the name of csi_wf_as_pilot.json

# To run runProcesses.m as ESTIMATION_ANALYSIS
1. Not implemented

# To run runSanityCheck.m
First you need to clone project-nsc-rice and checkout develop branch.
git clone https://gitlab.ad.sklk.us/projectrice/proj-nsc-rice
cd proj-nsc-rice
git checkout develop
cd code/skyscan/native
mkdir build && cd build
cmake ..
make -j

# To run runSanityCheck.m as SEC
1. You need to specify the runType = "SEC".
2. You need to specify your inputPathDir, outputPathDir, and sourcePathDir path directories. If not specified, it will take default hardcoded ones. Note, sourcePathDir is folder of clone project. 
3. Example of Linux command:  /usr/local/MATLAB/R2023a/bin/matlab -nodisplay -nosplash -nodesktop -r "inputPathDir='/dataset/BIN1-A'; outputPathDir='/space';runType='SEC'; sourcePathDir='/home/proj-nsc-rice'; run /space/radar/runSanityCheck(inputPathDir, outputPathDir, sourcePathDir, runType);exit;"

# To run runSanityCheck.m as DETECTION
1. You need to specify the runType = "DETECTION".
2. You need to specify your inputPathDir, outputPathDir, and sourcePathDir path directories. If not specified, it will take default hardcoded ones. Note, sourcePathDir is folder of clone project. 
3. Example of Linux command:  /usr/local/MATLAB/R2023a/bin/matlab -nodisplay -nosplash -nodesktop -r "inputPathDir='/dataset/BIN1-A'; outputPathDir='/space';runType='DETECTION'; sourcePathDir='/home/proj-nsc-rice'; run /space/radar/runSanityCheck(inputPathDir, outputPathDir, sourcePathDir, runType);exit;"
4. The max values of absolute errors will be reported. Those include Full and Sparse SEC and MF outputs values. 

# To run runProcesses.m as ESTIMATION
1. Not implemented

# To run analyzeData.m
1. You need to specify the modes = [6]; % modes 1-signal visual, 2-MF output, 3- ROC plots, 4-compute SNR, 5-signal level 6-Estimation 7-Detection
2. You need to specify the inputPathDir, movedPathDir, outputPathDir
3. Example of Linux command: /usr/local/MATLAB/R2023a/bin/matlab -nodisplay -nosplash -nodesktop -r "run('analyzeData.m');exit;"

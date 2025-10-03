
% Set the parameters
function [rP] = setRadarRange( rP)

switch rP.mode
    case 'BIN1-A'
        rP.gRange = [ 500, 800,  1000, 1200, 1500]*1e9;  % Chip Rate (GHz/s)
        rP.pwRange = [1, 2, 3, 4]*1e-6;  % pulsewidth
        rP.fRange = [0:4]*1e6 ; % Initial frequency f0 (MHz)
        rP.pw = 1*1e-6;

    case 'BIN1-B'
        rP.gRange = [ 3500, 3800,  4000, 4200, 4500]*1e9;  % Chip Rate (GHz/s)
        rP.pwRange = [1, 2, 3, 4]*1e-6;  % pulsewidth
        rP.fRange = [0:4]*1e6 ; % Initial frequency f0 (MHz)
        rP.pw = 2*1e-6;

    case 'BIN2-A'
        rP.gRange = [ 100, 300,  400, 500, 700]*1e9;  % Chip Rate (GHz/s)
        rP.pwRange = [45, 50, 55, 60]*1e-6;  % pulsewidth
        rP.fRange = [0:4]*1e6 ; % Initial frequency f0 (MHz)
        rP.pw = 50*1e-6;

    case 'BIN2-B'
        %rP.gRange = [ 50, 80,  100, 120, 150]*1e9;  % Chip Rate (GHz/s)
        % new one
        rP.gRange = [ 100, 200,  300, 400, 600]*1e9;  % Chip Rate (GHz/s)
        rP.pwRange = [75, 100, 125, 150]*1e-6;  % pulsewidth
        rP.fRange = [0:4]*1e6 ; % Initial frequency f0 (MHz)
        rP.pw = 100*1e-6;

    case 'BIN3-A'
        rP.gRange = [ 200, 400,  500, 600, 800]*1e9;  % Chip Rate (GHz/s)
        rP.pwRange = [5, 10, 15, 20]*1e-6;  % pulsewidth
        rP.fRange = [0:4]*1e6 ; % Initial frequency f0 (MHz)
        rP.pw = 10*1e-6;

    case 'BIN3-B'
        rP.gRange = [ 100, 200,  300, 400, 600]*1e9;  % Chip Rate (GHz/s)
        rP.pwRange = [25, 50, 75, 100]*1e-6;  % pulsewidth
        rP.fRange = [0:4]*1e6 ; % Initial frequency f0 (MHz)
        rP.pw = 50*1e-6;

    case 'OTHER-A'
        rP.gRange = [ 100, 300,  400, 500, 700]*1e9;  % Chip Rate (GHz/s)
        rP.pwRange = [40, 50, 60, 70]*1e-6;  % pulsewidth
        rP.fRange = [0:4]*1e6 ; % Initial frequency f0 (MHz)
        rP.pw = 50*1e-6;

    case 'OTHER-B'
        rP.gRange = [ 200, 400,  500, 600, 800]*1e9;  % Chip Rate (GHz/s)
        rP.pwRange = [5, 10, 15, 20]*1e-6;  % pulsewidth
        rP.fRange = [0:4]*1e6 ; % Initial frequency f0 (MHz)
        rP.pw = 10*1e-6;

    case {'DETECT','ESTIMATE'} 
        rP.gRange = [ 100, 200,  300, 400, 600]*1e9; %
		rP.gRange = [ 1e-9*0.05, 0.000005, 0.001, 0.05, 1, 50, 300,  400,  500, 1000, 4000]*1e9; %working[ 300, 400,  4000]*1e9; %[ 3500, 3800,  4000, 4200, 4500]*1e9; %[ 100 150 200 300 350 400 450 500 600 900 1000 2000 3000 4000 ]*1e9;  % Chip Rate (GHz/s)
        rP.pwRange = [1, 2, 3, 9, 10, 11, 49, 50,  99, 100 ]*1e-6;  % pulsewidth
        rP.fRange = [0:4]*1e6 ; % Initial frequency f0 (MHz)
        rP.pw = 10*1e-6;
end
rP.fRange1 = [-2:2]*1e6 ;
rP.fRange1 = [-1:1]*1e12*1e-6 ;
end
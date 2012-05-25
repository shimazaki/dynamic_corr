% This is a tutorial m-file for analysis of dynamic spike interactions
%
% Shimazaki H., Amari S., Brown E. N., and Gruen S.
% State-space Analysis of Time-varying Higher-order Spike Correlation 
% for Multiple Neural Spike Train Data. 
% PLoS Computational Biology 8(3): e1002385. 
% http://dx.doi.org/10.1371/journal.pcbi.1002385
%
% Dec 12, 2011 Author Hideaki Shimazaki
% http://2000.jukuin.keio.ac.jp/shimazaki

% Note: 
% Currently, the method works for up to N~10 neurons for pairwise analysis. 
% N~5 for triple-wise analysis. Try starting your analysis from small 
% number of subset neurons.
%
% Basic terminolgy:
% N: The number of neurons.
% R: The order of interactions in the model.
%   R = [1]: Rate analysis
%   R = [1,2]: Pairwise analysis 
%   R = [1,2,3]: Triple-wise analysis

%%%%%%%%%%%%%%%%%% Let's start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By running the script, you will obtain the result for analysis of N=3
% neurons. 
% Figure 1 displays spike timing of each neuron from 
% top to bottom. In each panel, the ordinate represents trials. 
% Figure 2 displays occurence rates of joint spike event of r (=1,2,3 from 
% top to bottom) neuron(s) within the bin size that will be specified below.
% Figure 3 displays time-varying log-linear parameters of the r-th order. 
% Figure 4 displays optimized hyper-parameters in a state equation. 
% Figure 5 displays snapshots of the time-varying log-linear parameters. 

%%%%%%%%%%%%%%%%%% How to use the results %%%%%%%%%%%%%%%%%%%%%%
% After running this script, you will obtain a structure array 'data'. 
% 'data.raw' contains an input spike timing data used for analysis.
% 'data.model' contains estimated dynamic interactions from the data.
% Each values of dynamic interactions are written in 
% data.model.param.smoother.theta
% 95% credible interval is obtained using
% 1.95996 * sqrt( data.model.param.smoother.diagW ) 
% Each row in the 'theta' vector indicates interactions among subset neurons.
% The id for the subset neuros (which neurons' interaction?) is in a matrix
% data.model.binary.id 
% '1' indicates an active neuron, '0' indicates an inactive neuron.
% The order of interactions for each element (the number of '1's) is in 
% data.model.binary.order
% Hence, to extract the dynamics of pairwise interactions only, use
% data.model.param.smoother.theta(data.model.binary.order==2,:)

clear all; close all;

%%%%%%%%%%%%%%%%%%% Matlab settings %%%%%%%%%%%%%%%%%%%%%%
% Please add the paths to 'ssloglin' and 'share' directories.
% 'ssloglin' contains m-files for state-space analysis of ineteractions.
% 'share' directory contains m-files for the log-linear analysis in general.
addpath(genpath('../ssloglin')); 
addpath(genpath('../share')); 

%%%%%%%%%%%%%%%%%%% Data %%%%%%%%%%%%%%%%%%%%%%
load Example_3Neurons.mat raw; 
data.raw = raw;
% The file Example_3Neurons.mat contains 3 neurons spike timing data
% and additional information (Number of neurons (N=3), and trials (n=300), 
% sampling resolution, and a period of observation(period)).
% The spike timing data is in a cell array, `raw.xs'.  
% Here, raw.xs is {{1x300 cell}  {1x300 cell}  {1x300 cell}}. 
% Each sub cell array contains data from a single neuron.
% In raw.xs{1}, spike timing data of 300 trials from the first neuron
% is stored in vectors with different lengths.
% 
% Complete structure of spike timing data is as follows:
%data.raw.xs: {{1x300 cell}  {1x300 cell}  {1x300 cell}}; 
%data.raw.N = 3;            %Number of neurons.
%data.raw.n = 300;          %Number of repeated trials.
%data.raw.D = 0.001;        %Sampling resolution [s] of spike timing data.
%data.raw.period = [0 0.5]; %Period of observation.
%data.raw.miss = [];        %Optional
%data.raw.ext = [];         %Optional

%%%%%%%%%%%%%%%%%%% Model Settings %%%%%%%%%%%%%%%%%%%%%%
% 
% The following settings specify the log-linear model used for analysis.
data.model.struct = struct;     %Please do not uncomment this line.
data.model.struct.N = 3;        %Number of neurons.
data.model.struct.D = 0.001;    %Bin size [s] to make binary data.
data.model.struct.n = 100;      %Number of repeated trials for analysis.
data.model.struct.R = [1 2 3];  %The order of interactions in the model.
                                %Use R = [1 2] for pairwise analysis.

%%%%%%%%%%%%%%%%%%% Optional Settings %%%%%%%%%%%%%%%%%%%%%%
% There are variables you can tune according to your data. 
% The default values are shown. Uncomment and change the variables.
% The only free parameter in the analysis is data.model.init.Sig (and 
% the bin size, data.model.struct.D).
data.model.init = struct;       %Please do not comment out this line.

%data.model.init.state = struct('Q',1,'F',0,'G',0); 
% init.state determins the structure of the state process. 
% Q: 1: Q=\gamma I, a single smoothing parameter, \gamma, is optimized.
% Q: 2: a diagonal matrix Q is optimzied.
% Q: 3: a full matrix Q is optimzied.
% Q: 0: No noise in the state process => stationary log-linear analysis.
%
% F: 1: the first order autoregressive paramter, F, is optimized.
% F: 0: F=I, random gaussian walk model
%
% G: Either spike history or stimulus effect. Optional.
%
% Initial values of the hyper-parameters. 
% Except for 'Sig', they will be updated and optimized. 
%data.model.init.Sig = diag( 1*ones(1,d) ); 
%data.model.init.Q = diag(.05*ones(1,d)); 
%
% Termination criterion
%data.model.init.em_max_iteration = 200; %Maximum number of iterations
%data.model.init.em_min_logLdiff = 0.1; %Iteration stops when increase 
% in log-likelihood is less than 0.1.

%%%%%%%%%%%%%%%%%%% State-space analysis %%%%%%%%%%%%%%%%%%%%%
data.model = GenModel(data.raw, data.model.struct, data.model.init);

%%%%%%%%%%%%%%%%%%% Display results %%%%%%%%%%%%%%%%%%%%%%%%%%
% If your monitor display does not match with default settings, 
% edit 'Disp.m'. 
Disp(data.raw,data.model);

%%%%%%%%%%%%%%%%%%% Save your data %%%%%%%%%%%%%%%%%%%%%%
state = data.model.init.state;
str = strcat('data/example_','Q',num2str(state.Q),...
    '_F',num2str(state.F),'_n',num2str(data.model.struct.n));
%save(str,'data');      %remove comment to save the data.


% Revision history
% ver.0.11  2012/04/07 Added more explanations.
% ver.0.1   2012/03/16 Extended.  
% ver.0.0   2011/12/12 The first draft of tutorial was created. 
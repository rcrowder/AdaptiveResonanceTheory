%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an implementation of an ARTMAP variant,
% as described in:
% Gaddam, C. S. (2007).
% Feature Selection via transmitter depletion in ARTMAP. Online Document, xx(x) xxx-xxx.
% Boston, MA: Boston University.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programmed by Chaitanya Sai (August 2007-08)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The notation follows the article's notation,
% as does the headers for each step of the
% algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% close all;

function [a,b,c] = fuzzyARTMAPTester(dataStructTemp,lambda_Attention)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


inputNet = struct('M', {[]}, ... % Number of features
    'TestOnly', {0}, ... % Train & Test or Test only? unless specified otherwise
    'test_fast', {1}, ... % Use fast tester unless specified otherwise
    'compute_overlap', {0}, ... % Indirectly measure overlap by looking at how many nodes overlap with a test_point
    'C', {0}, ... % Number of committed coding nodes (paper says set to 0, but harder to Matlab)
    'plotSteps',{0}, ... % % Iterate through individual learning steps
    'NodeList', {[]}, ... % List of node associated with each input
    'maxNumCategories', {100}, ...
    'numClasses', {[]}, ... % The total number of classes
    'w', {[]}, ... % coding node weight vector
    'W', {[]}, ... % output class weight vector
    'rho_bar', {0}, ... % rho_bar, the baseline vigilance. baseline=0 maximizes code compression
    'rho', {0.0}, ... % vigilance
    'alpha', {.1}, ... % Signal Rule Parameter
    'alphaTrain',{[]},... % Store when reverting alpha value in TESTING
    'beta', {1.0}, ...  % Learning fraction
    'gamma',{.000001},... % fraction additive in numerator and denominator of matching signal
    'epsilon', {-.00001}, ... % Match Tracking (codes inconsistent cases)
    'base_vigilance', {0.0}, ... % Baseline vigilance
    'p', {1.0}, ... % CAM rule power
    'dataSubsets', {1}, ... % 4-fold cross-validation
    'votes', {5}, ... % Number of voting systems
    'learningRate', {1.0}, ...
    'win_sequence',{[]},...
    'act_sequence',{[]},...
    'search_cycles', {0},... % Total number of search cycles during training
    'lambda_Attention',{[]},...
    'dataStore',{[]},... % Store data if using structs as inputs
    'numEpochs',{1},...%Number of times to present input
    'model_num',{[]},... % Depletion Model Number
    'type_depletion',{4},... %Options are 0 1 2 3 4
    'deplete_fast',{[]},...
    'learn_trail',{[]},...
    'learn_description',{''},...
    'short_learn_descript',{''},...
    'uncommitted_fail',{0},...
    'NI_fail',{[]},...
    'NI_fail_addNode',{[]},...
    'node_created',{[]},...
    'start_conservative',{[]},...
    'ovlp_net',{[]},...
    'A_dep_dist',{[]},...
    'e_mod_dist',{[]},...
    'e_realLow',{0},...
    'e_zero',{0},...
    'diffE_count',{[]},...
    'perc_error',{[]},... % percentage error added to labels in training
    'var_error_features',{[]},... % variance of error added to features in training
    'testAccuracy',{[]}); % Attention gain parameter



% inputNet.lambda_Attention=0.5;

inputNet.lambda_Attention=lambda_Attention;
inputNet.deplete_fast = 0;  %TURN OFF (=0) IF DEPLETING USING ITERATIVE VERSION
inputNet.start_conservative =0;
inputNet.DEFAULT_RETEST=1;
inputNet.alpha=10^(-8);
inputNet.gamma=0.0000001;
inputNet.plotSteps=0;
inputNet.model_num=32;
dataSetNum=2;

%Train for number of epochs
inputNet.numEpochs=1;

%%%Inject label errors into dataset | e.g perc_error = 10 implies 10% of
%%%training data have erroneous class labels during training
perc_error=0;
inputNet.perc_error=perc_error;

%%%Inject feature noise into dataset | e.g var_error_features = 10 implies 10% of
%%%training data have noise added to input features
var_error_features=0;
inputNet.var_error_features=var_error_features;

inputNet.compute_overlap=1;

if (dataSetNum==13)
    %Matrix manipulation becomes too memory intensive for movie genres
    %benchmark
    inputNuet.test_fast=0;
end


if (inputNet.perc_error~=0)
    disp('Adding error in training set!');
elseif (inputNet.var_error_features~=0)
    disp('Adding feature error in training set!');
end



if ~exist('dataStructTemp','var')
    disp('Training and Test data must be provided in the form of a MATLAB struct dataStructTemp.');
    error('Use createDataStruct.m to create an input MATLAB struct');
else


    [a,b,c]=Default_ARTMAP(dataStructTemp,inputNet);

    disp('Output class predictions are stored in variable ''a''.');
    disp('Distributed output predictions are stored in ''b''.');
    disp('The biased ARTMAP network details are stored in ''c''');

end
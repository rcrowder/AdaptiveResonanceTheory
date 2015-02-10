%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an implementation of an ARTMAP variant,
% as described in:
% Gaddam, C. S. (2007).
% Feature Selection via transmitter depletion in ARTMAP. Online Document, xx(x) xxx-xxx.
% Boston, MA: Boston University.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programmed by Sai Chaitanya Gaddam (August 2007-08)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The notation follows the article's notation,
% as does the headers for each step of the
% algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function artmap_network = DEFARTMAP_Create_Network(M, numClasses)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an implementation of Default ARTMAP,
% as described in CAS/CNS Technical Report TR-2003-008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programmed by Ogi Ogas
%%%%%%%%%%%%%%%%%%%%%%%%%
% The notation follows the article's notation,
% as does the headers for each step of the
% algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: numFeatures = 2*M
% 02/21/05  Added search cycle to structure
% 08/31/05  Added list of nodes (for each training input)
%
% This is all B.2 Initialization

% Check the ranges of the input parameters.
if(M < 1)
    error('The number of features must be a positive integer.');
end
if(numClasses < 2)
    error('There must be more than one class.');
end

% Create and initialize the weights to F2 (the coding node weight vector).
% w(i,j)= weight(F0 Node, F2 Node)
w = ones(2*M, 1);

% Create and initialize the weights to Fab (the output class weight vector/the map field).
% W(j,k) = weight (F2 Node, Fab Node/Class)
W = [];

% Create the structure and return.
artmap_network = struct('M', {M}, ... % Number of features
    'plotSteps',{0}, ... % % Iterate through individual learning steps
    'C', {0}, ... % Number of committed coding nodes (paper says set to 0, but harder to Matlab)
    'NodeList', {[]}, ... % List of node associated with each input
    'maxNumCategories', {100}, ...
    'numClasses', {numClasses}, ... % The total number of classes
    'w', {w}, ... % coding node weight vector
    'W', {W}, ... % output class weight vector
    'rho_bar', {0.0}, ... % rho_bar, the baseline vigilance. baseline=0 maximizes code compression
    'rho', {0.0}, ... % vigilance
    'alpha', {.01}, ... % Signal Rule Parameter
    'alphaTrain',{[]},... % Store when reverting alpha value in TESTING
    'beta', {1.0}, ...  % Learning fraction
    'gamma',{[]},... % fraction additive in numerator and denominator of matching signal
    'epsilon', {-.001}, ... % Match Tracking (codes inconsistent cases)
    'base_vigilance', {0.0}, ... % Baseline vigilance
    'p', {1.0}, ... % CAM rule power
    'dataSubsets', {1}, ... % 4-fold cross-validation
    'votes', {5}, ... % Number of voting systems
    'learningRate', {1.0}, ...
    'win_sequence',{[]},...
    'act_sequence',{[]},...
    'search_cycles', {0},... % Total number of search cycles during training
    'lambda_Attention',{0},...
    'dataStore',{[]},... % Store data if using structs as inputs
    'model_num',{[]},... % Depletion Model Number
    'type_depletion',{4},... %Options are 0 1 2 3 4
    'deplete_in_mismatch',{0},...
    'deplete_in_X_denom',{0},...
    'use_type_three_tracking',{0},...
    'start_conservative',{0},...
    'num_type_three_used',{0},...
    'learn_trail',{[]},...
    'learn_description',{''},...
    'short_learn_descript',{''},...
    'NI_fail',{[]},...
    'NI_fail_addNode',{[]},...
    'node_created',{[]},...
    'uncommitted_fail',{0},...
    'e_realLow',{0},...
    'e_zero',{0},...
    'diffE_count',{[]},...
    'testAccuracy',{[]}); % Attention gain parameter
return


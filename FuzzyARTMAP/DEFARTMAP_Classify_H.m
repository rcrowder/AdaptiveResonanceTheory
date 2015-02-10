function [classes, sigmaS, yMatrix]= DEFARTMAP_Classify_H(artmap, testdata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an implementation of Default ARTMAP,
% as described in CAS/CNS Technical Report TR-2003-008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programmed by Ogi Ogas, Modified By (Chaitanya Sai)
%%%%%%%%%%%%%%%%%%%%%%%%%
% The notation follows the article's notation,
% as does the headers for each step of the
% algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DISTRIBUTED TESTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 02/21/05: fixed a bug where y was in column format instead of row format
% during IGCam rule
% 09/11/05: checks to make sure test data has same features as training
%           data
[numFeatures, n] = size(testdata);
M_test = numFeatures/2;
if (artmap.M~=M_test)
    error('The test data does not have the same number of features as the training data!')
end

% Set up the return variables.
classes = zeros(1, n);
      
% Classify and learn on each sample.
for sampleNumber = 1:n
        
        % C.2: SELECT NEXT INPUT VECTOR A
        A = testdata(:, sampleNumber);
        artmap.rho = artmap.rho_bar;
        
        % C.3: RESET THE CODE
        y = zeros(artmap.C,1); % the code at F2. At each node j, activity reset to 0.
        
        % C.4 CALCULATE SIGNALS TO COMMITTED CODING NODES
        T = DEFARTMAP_choice_function(artmap,A);
        
        % C.5: CALCULATE LAMBDA
        Lambda=[]; Lambda_prime=[];     
        Lambda = find(T>artmap.alpha.*artmap.M);
        Lambda_prime = find(T==artmap.M);
            

        % C.6: INCREASED GRADIENT (IG) CAM RULE
        if ~isempty(Lambda_prime)  % Point Box Case -- w_j = A 
              y(Lambda_prime)=(1./length(Lambda_prime));
              y=y';
        else
              sumL = sum(1./(artmap.M-T(Lambda)));
              y(Lambda)=(1./(artmap.M-T(Lambda)))./sumL;
              y=y';             
        end
        % Store the distributed y values for all test points. This is used
        % in the pseudo-inverse part
        yMatrix(sampleNumber,:)= y;
        % C.7 DISTRIBUTED OUTPUT CLASS PREDICTIONS
        sigma = artmap.W'*y';
        [dummy, class]=max(sigma);
        sigmaS(:,sampleNumber)=sigma;
        classes(sampleNumber)=class;
   
end      % for sampleNumber = 1:numSamples

return
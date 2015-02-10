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



function [classes, yMatrix, varargout]= DEFARTMAP_Classify(artmap, testdata)


[numFeatures, n] = size(testdata);
yMatrix=[];
M_test = numFeatures/2;
if (artmap.M~=M_test)
    error('The test data does not have the same number of features as the training data!')
end

% Set up the return variables.
classes = zeros(1, n);
top_node_num=round(artmap.C*.02);
node_num_penalty=(sum(artmap.W,1)').^(2/3);

%ovlp_net = zeros(1, n);
% Classify and learn on each sample.
for sampleNumber = 1:n

    if (mod(sampleNumber,400)==0)
        disp('Sample Number:')
        disp(sampleNumber')
    end

    % C.2: SELECT NEXT INPUT VECTOR A
    A = testdata(:, sampleNumber);
    artmap.rho = artmap.rho_bar;

    % C.3: RESET THE CODE
    y = zeros(artmap.C,1); % the code at F2. At each node j, activity reset to 0.

    % C.4 CALCULATE SIGNALS TO COMMITTED CODING NODES
    T = DEFARTMAP_choice_function(artmap,A);
    %[T,ovlp_net(sampleNumber)] = DEFARTMAP_choice_function(artmap,A);

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
    % sampleNumber;
    % C.7 DISTRIBUTED OUTPUT CLASS PREDICTIONS
            dum=sort(y,'descend');
            y_rect=y.*(y>mean(y)).*(y>dum(top_node_num));
            
            sigma = artmap.W'*y_rect';
            [dummy, class]=max(sigma./node_num_penalty);
            classes(sampleNumber)=class;
    % C.7 DEFAULT WTA OUTPUT CLASS PREDICTIONS

    %WHAT TO DO IN A CLASH
    %classes(sampleNumber) = min(find(artmap.W(y==max(y),:)));

end      % for sampleNumber = 1:numSamples

if (nargout==3)
    varargout(1)={ovlp_net};
end
% yMatrix=0;
return
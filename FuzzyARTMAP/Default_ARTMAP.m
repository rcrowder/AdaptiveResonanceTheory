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
function [test_classes, yMatrix, trainingNet] = Default_ARTMAP(varargin)

if nargin==1
    %Input is a datastruct containing training and test data%
    trainingdata=varargin{1}.training_input;
    training_classes=varargin{1}.training_output;
    testdata=varargin{1}.test_input;
    test_classes=varargin{1}.test_output;
elseif nargin==2
    %The parameters for initializing biased ARTMAP are provided as a struct%
    initNet=varargin{2};
    trainingdata=repmat(varargin{1}.training_input,1,initNet.numEpochs);
    training_classes=repmat(varargin{1}.training_output,1,initNet.numEpochs);
    testdata=varargin{1}.test_input;
    test_classes=varargin{1}.test_output;
elseif nargin==3
    %Training, test input data and training classes are provided as separate matrices%
    trainingdata=varargin{1};
    training_classes=varargin{2};
    testdata=varargin{3};
elseif nargin==4
    %Training, test input data and training classes are provided as
    %separate matrices%
    %The parameters for initializing biased ARTMAP are provided as a
    %struct%
    trainingdata=varargin{1};
    training_classes=varargin{2};
    testdata=varargin{3};
    initNet=varargin{4};
end




ON=1; OFF=0;
VISUALS = OFF; % Only for 2-D visualization, 2-d input


% Make sure the user specifies the input parameters.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Class Vector: class(k) from class(1) to class(n)
%   ** Classes should be integers greater than 0
if (~isempty(find(training_classes==0)))
    error('The classes must be consecutive integers greater than 0.');
end
numClasses = max(training_classes); % Calculate the number of classes
% if numClasses<2
%     error('Needs to be more than one training class.')
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feature Vector: a, a(1) to a(M)
% M is the number of features (pre-complement coding!)
%
% Data should be array of feature vectors, M x n,
% where n is the number of vectors (data points)
%

[M,n]=size(trainingdata);
% READ THE TEST DATA
[test_M,test_n] = size(testdata);
if (test_M ~= M)
    error ('The test data has a different number of features than the training data.')
end

if (VISUALS==ON)
    GRAPH_DATA(trainingdata, training_classes)
    title('training data')
    input('hit to close')
    close all
end

% Normalize the feature vectors.
alldata = [trainingdata testdata];
%alldata = DEFARTMAP_Normalize(alldata);
alldata = DEFARTMAP_Complement_Code(alldata);

trainingVectors = alldata(:,1:n);
testVectors = alldata(:,n+1:n+test_n);

% B.1 COMPLEMENT CODE THE INPUT
%%% THIS IS F0 %%%
% A = complement coded feature vector, A = [a a^c]
% (actually, in this implementation, a and a^c are interwoven in A, but
% this has no effect on the algorithm)
% A is 2*M x n
% Mx2 = 2*M
% but remember, |A|=M -- this is also a kind of normalization.


% B.2: INITIALIZE THE NETWORK

if nargin==1
    initNet = DEFARTMAP_Create_Network(M, numClasses);
elseif (nargin==2) && (initNet.TestOnly==0)
    w = ones(2*M, 1);
    W = [];
    initNet.M=M;
    initNet.numClasses=numClasses;
    initNet.w=w;
    initNet.W=W;
    initNet.dataStore=varargin{1};
    clear w W;
elseif nargin==3
    initNet = DEFARTMAP_Create_Network(M, numClasses);
elseif (nargin==4) && (initNet.TestOnly==0)
    w = ones(2*M, 1);
    W = [];
    initNet.M=M;
    initNet.numClasses=numClasses;
    initNet.w=w;
    initNet.W=W;
    clear w W
end

% B.3-B.11: TRAIN THE NETWORK ON ALL THE INPUT VECTORS
%display('___________NOW TRAINING________________')



if (initNet.TestOnly==0)
    if initNet.plotSteps==0;
        trainingNet = DEFARTMAP_biasedLearn(initNet, trainingVectors, training_classes);
    elseif initNet.plotSteps==1;
        trainingNet = DEFARTMAP_biasedLearnPlot(initNet, trainingVectors, training_classes);
    end

else
    trainingNet=initNet;
end



if (initNet.compute_overlap==0)
    if (initNet.test_fast==0)
        [test_classes, yMatrix]= DEFARTMAP_Classify(trainingNet, testVectors);
    elseif (initNet.test_fast==1)
        %%%%%%%%%%%%%%%%%Faster Classification using matrix comparisons%%%%%%%%%
         try
            [test_classes, yMatrix]= DEFARTMAP_Classify_Fast(trainingNet, testVectors);
        catch
            [lastmsg,lastid]=lasterr;
            if isempty(strfind(lastid,'nomem'))
                error(lastid,lastmsg)
            else
                warning(lastid,'Not enough memory; reverting to slow testing; may take significantly longer');
                [test_classes, yMatrix]= DEFARTMAP_Classify(trainingNet, testVectors);
            end

        end
    end
    
else

    if (initNet.test_fast==0)
        [test_classes, yMatrix,trainingNet.ovlp_net]= DEFARTMAP_Classify(trainingNet, testVectors);
    elseif (initNet.test_fast==1)
        %%%%%%%%%%%%%%%%%Faster Classification using matrix comparisons%%%%%%%%%
         try
            [test_classes, yMatrix]= DEFARTMAP_Classify_Fast(trainingNet, testVectors);
        catch
            [lastmsg,lastid]=lasterr;
            if isempty(strfind(lastid,'nomem'))
                error(lastid,lastmsg)
            else
                warning(lastid,'Not enough memory; reverting to slow testing; may take significantly longer');
                [test_classes, yMatrix]= DEFARTMAP_Classify(trainingNet, testVectors);
            end

        end
    end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   [test_classes, yMatrix,ovlp_net]= DEFARTMAP_Classify(trainingNet, testVectors);
%   trainingNet.ovlp_net=ovlp_net;
trainingNet.dataStore=[];

try (sum(trainingdata(:)~=testdata(:))==0)
    'YESSSS'
    trainingNet.testAccuracy=sum(test_classes(:)==training_classes(:))/length(training_classes)*100;
catch
    'Not Testing on Training Set'
end


if (VISUALS==ON)
    GRAPH_DATA(testdata, test_classes)
    title('test data, as classified by Default ARTMAP')
    input('hit to close')
    close all
end


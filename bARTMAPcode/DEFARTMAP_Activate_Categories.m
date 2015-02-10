function categoryActivation = ART_Activate_Categories(input, weight, bias, net)
% ART_Activate_Categories    Activates the categories in an ART network.
%    CATEGORYACTIVATION = ART_Activate_Categories(INPUT, WEIGHT, BIAS)
%    This function returns a vector of category activations for the given
%    input vector, weight matrix, and bias value.
% 
%    The input parameters are as follows:
%    The INPUT is a vector of size NumFeatures that contains the input
%    signal into the network. The WEIGHT is a matrix of size 
%    NumFeatures-by-NumCategories which holds the weights of the network.
%    The BIAS is the constant that is used to differentiate between very
%    similar category activation values. The length of the INPUT vector
%    must equal the number of rows in the WEIGHT matrix, and the BIAS
%    value must be within the range [0, 1] (although values very near
%    0 are best).
%
%    The return parameter is as follows:
%    The CATEGORYACTIVATION is a vector of size NumCategories that
%    holds the activation value for each category.


% Make sure the user supplied the required parameters.
if(nargin ~= 3)
    error('You must specify the 3 input parameters.');
end

% Check the size and range of the parameters.
[numFeatures, numCategories] = size(weight);
if(length(input) ~= numFeatures)
    error('The length of the input and rows of the weights do not match.');
end
if((bias < 0) | (bias > 1))
    error('The bias must be within the range [0, 1].');
end

% Set up the return variable.
categoryActivation = ones(1, numCategories);

% B.7
% Calculate signals to committed coding nodes
% j=1...C
%        Activation(j) = |Input^Weight(j)| + (1-alpha)*(M-|weight(j)|)
for j = 1:numCategories    
    matchVector = min(input, weight(:, j));
    weightLength = sum(weight(:, j));
    T(j) = sum(matchVector) + (1 - net.alpha).*(numFeatures - weightLength);
end


return
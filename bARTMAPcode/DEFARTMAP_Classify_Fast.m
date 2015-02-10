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

function [classes, yMatrix, varargout]= DEFARTMAP_Classify_Fast(artmap, testdata)

%Repeat data in 3rd dimension
testdata_r=repmat(testdata,[1 1 artmap.C]);
[numFeatures, n] = size(testdata);

%Repeat data in 3rd dimension and permute
w_r=permute(repmat(artmap.w,[1 1 n]),[1 3 2]);

%This finds the minimum of testdata_r and w_r
min_vals=w_r.*sign(max(testdata_r-w_r+eps,0))+testdata_r.*sign(max(w_r-testdata_r,0));

%Sum up along the feature dimension
min_vals=permute(sum(min_vals,1),[2 3 1]);

%Provide sum of weight vectors
w_sums=repmat(sum(artmap.w,1),[n 1]);

%The choice values
T_vals=min_vals + (1-artmap.alpha).*(artmap.M - w_sums);
yMatrix=T_vals;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DISTRIBUTED TESTING
% Lambda_prime=(T_vals==artmap.M);
% Lambda_prime_inv=repmat(1./sum(Lambda_prime,2),[1 artmap.C]).*Lambda_prime;
% Lambda_prime_inv(isnan(Lambda_prime_inv))=0;
% Lambda=T_vals.*(T_vals>artmap.alpha*artmap.M);
% Lambda_inv=1./(artmap.M-Lambda);
% Lambda_inv_norm=Lambda_inv./repmat(sum(Lambda_inv,2),[1 artmap.C]);
% Choose_val=repmat(sign(sum(Lambda_prime_inv,2)),[1 artmap.C]);
% yMatrix=Lambda_prime_inv.*Choose_val+Lambda_inv_norm.*(~Choose_val);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           WTA Class Mapping
[dum,MaxY]=max(yMatrix,[],2);
Class_map=artmap.W*[1:1:size(artmap.W,2)]';
classes=Class_map(MaxY);

if nargout==3
    varargout{1}=sum(w_sums==min_vals,2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function complementCodedData = DEFARTMAP_Complement_Code(data)
% ART_Complement_Code    Complement codes the data for use with an ART network.
% vertical is M, the FEATURES 
% Horizontal is N, the NUMBER of data points.
%complementCodedData = [data; 1-data];

%data=(data+4)/9;
complementCodedData = [data; 1-data];
% Do the complement coding for each sample.



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

%load from flixLoadData.m


if (~exist('moviMAT','var'))
    fid=fopen('E:\netflix_noprobe\download\netflixresultsgenerator\algorithms\svd\release\moviFile_300_m.bin','rb');
    moviMAT=fread(fid,[17770 64],'float32');
    fclose(fid)
end


if (~exist('movie_genres_numeric','var'))
    load MovieGenreData
end
num_in_genre=zeros(1,21);
for i=1:21;num_in_genre(i)=sum(movie_genres_numeric==i);end
[dum,dum_a]=ismember(movie_genres_numeric,find(num_in_genre>400));
movie_genres_numeric_reduced=movie_genres_numeric(dum);
movie_labels_reduced=movie_genres_mapping(unique(movie_genres_numeric_reduced));
rev_map(num_in_genre>400)=(1:1:length(find(num_in_genre>400)));
movie_genres_numeric_reduced=rev_map(movie_genres_numeric_reduced);


% 
% dataStructTemp=createDataStruct([moviMAT movie_genres_numeric],'rows','end',0,'movie_genres',(1:1:17770));
% dataStructTemp.labels=movie_genres_mapping;



dataStructTemp=createDataStruct([moviMAT(dum,:) movie_genres_numeric_reduced'],'rows','end',85,'movie_genres',setdiff((1:1:17770).*(dum'),0));
dataStructTemp.labels=movie_labels_reduced;
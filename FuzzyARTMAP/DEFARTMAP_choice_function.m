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

function [T,varargout] = DEFARTMAP_choice_function(artmap,A)


if (nargout<2)

    T = zeros(1,artmap.C);
    for j=1:artmap.C
        size_minAw=sum(min(A,artmap.w(:,j)));
        size_w=sum(artmap.w(:,j));
        T(j) = size_minAw + (1-artmap.alpha).*(artmap.M - size_w);
    end

elseif (nargout==2)
%Also compute number of coding boxes with which the presented point A
%overlaps
    T = zeros(1,artmap.C);
    Ovlp=zeros(1,artmap.C);
    for j=1:artmap.C
        size_minAw=sum(min(A,artmap.w(:,j)));
        size_w=sum(artmap.w(:,j));
        T(j) = size_minAw + (1-artmap.alpha).*(artmap.M - size_w);
        Ovlp(j) = (size_minAw==size_w);
    end
    % The number of coding boxes with which the presented point A overlaps
    varargout(1)={sum(Ovlp)};

end



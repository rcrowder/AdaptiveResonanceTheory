%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of an implementation of the biased ARTMAP Matlab code package,
% as described in CAS/CNS Technical Report TR-2009-xxx
% Biased AR: A neural architecture that shifts attention towards previously disregarded features 
% following an incorrect prediction.
% Boston, MA: Boston University.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programmed by Sai Chaitanya Gaddam(May 2007-09)
%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The notation follows the article's notation,
% as does the headers for each step of the
% algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function dataStructTemp  = createDataStruct(in_mat,row_col_string,labels_end_begin,perc_train,data_description,varargin)

if strcmp(row_col_string,'rows')
    rand_order=randperm(size(in_mat,1));
    %rand_order=[1:1:size(in_mat,1)];
    in_mat= in_mat(rand_order,:);
elseif strcmp(row_col_string,'cols')
    in_mat=in_mat';
    rand_order=randperm(size(in_mat,1));
    %rand_order=[1:1:size(in_mat,1)];
    in_mat= in_mat(rand_order,:);
end

if strcmp(labels_end_begin,'end')
    label_vals=in_mat(:,end);
    in_mat=in_mat(:,1:end-1);

elseif strcmp(labels_end_begin,'begin')
    label_vals=in_mat(:,1);
    in_mat=in_mat(:,2:end);
end

%%%%Make sure all features are 0-1 range

max_vals=max(in_mat,[],1);
min_vals=min(in_mat,[],1);
mean_vals=mean(in_mat,1);
std_vals=std(in_mat,1);

if size(in_mat,2)<50
disp('max')
disp(max_vals)
disp('min')
disp(min_vals);
disp('mean')
disp(mean_vals)
disp('std')
disp(std_vals)
else
    figure(1)
    subplot(4,1,1)
    bar(max_vals);
    title('max')
    subplot(4,1,2)
    bar(min_vals);
    title('min');
    subplot(4,1,3)
    bar(mean_vals);
    title('mean');
    subplot(4,1,4)
    bar(std_vals);
    title('std');

    
end


%%%Shift into 0+%%%%%
in_mat=in_mat-repmat(min(min_vals,0),size(in_mat,1),1);

disp('After shift into 0+')

max_vals=max(in_mat,[],1);
min_vals=min(in_mat,[],1);
mean_vals=mean(in_mat,1);
std_vals=std(in_mat,1);

if size(in_mat,2)<50
disp('max')
disp(max_vals)
disp('min')
disp(min_vals);
disp('mean')
disp(mean_vals)
disp('std')
disp(std_vals)
else
    figure(2)
    subplot(4,1,1)
    bar(max_vals);
    title('max')
    subplot(4,1,2)
    bar(min_vals);
    title('min');
    subplot(4,1,3)
    bar(mean_vals);
    title('mean');
    subplot(4,1,4)
    bar(std_vals);
    title('std');

    
end

%%%Compress into 0-1
in_mat=in_mat./repmat(max(max_vals,1),size(in_mat,1),1);

disp('After compression into [0-1]')

max_vals=max(in_mat,[],1);
min_vals=min(in_mat,[],1);
mean_vals=mean(in_mat,1);
std_vals=std(in_mat,1);

if size(in_mat,2)<50
disp('max')
disp(max_vals)
disp('min')
disp(min_vals);
disp('mean')
disp(mean_vals)
disp('std')
disp(std_vals)
else
    figure(3)
    subplot(4,1,1)
    bar(max_vals);
    title('max')
    subplot(4,1,2)
    bar(min_vals);
    title('min');
    subplot(4,1,3)
    bar(mean_vals);
    title('mean');
    subplot(4,1,4)
    bar(std_vals);
    title('std');

    
end
%%%%%%%%%%%%%

label_vals=label_vals(:);
if (min(label_vals)==0)
    disp('Shifting labels from [0 n-1]to [1 n]')
    label_vals=label_vals+1;
elseif (min(label_vals)==0)
    error('Weird: Labels are negative');
end


num_train=round(perc_train*size(in_mat,1)/100);
dataStructTemp.training_input=in_mat(1:num_train,:)';
dataStructTemp.training_output=label_vals(1:num_train);
dataStructTemp.test_input=in_mat(num_train+1:end,:)';
dataStructTemp.test_output=label_vals(num_train+1:end);
dataStructTemp.description=data_description;
dataStructTemp.descriptionVerbose=data_description;
if (nargin==5)
dataStructTemp.randOrderTrain=rand_order(1:num_train);
dataStructTemp.randOrderTest=rand_order(num_train+1:end);
else
dataStructTemp.randOrderTrain=varargin{1}(rand_order(1:num_train));
dataStructTemp.randOrderTest=varargin{1}(rand_order(num_train+1:end));
end



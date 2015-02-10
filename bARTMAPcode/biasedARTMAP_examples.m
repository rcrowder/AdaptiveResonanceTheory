%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an implementation of an ARTMAP variant,
% as described in:
% Gaddam, C. S. (2007).
% Feature Selection via transmitter depletion in ARTMAP. Online Document, xx(x) xxx-xxx.
% Boston, MA: Boston University.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programmed by Chaitanya Sai (August 2007-08)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The notation follows the article's notation,
% as does the headers for each step of the
% algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% close all;


if ~exist('dataStruct_biasedART','var')
    load biasedARTbenchmarks
    load dataset_descriptions
    disp('Finished loading datasets')
end

close all;

%%%%%%%%DATASETS%%%%%%%%%%%%%
%
%     1'six_point_test'
%     2'Stripes_Sparse'
%     3'Stripes_Dense'
%     4'CIS_Sparse'
%     5'CIS_Dense'
%     6'Checkerboard_Sparse'
%     7'Checkerboard_Dense'
%     8'Binary_five'
%     9'Test on Stripe 3 Capped'
%     10'Test on Stripe 2 Capped'
%     11'Test on Stripe 1 Capped'
%     12'Test on Stripe 0 Capped'
%     13'Genre prediction'
%%%%%%%%DATASETS%%%%%%%%%%%%%


disp({dataset_strings.string_val}')
dataSetNum=input('Input dataset?');

lambda_Attention=input('lambda Value (0 -> fuzzy ARTMAP)?');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    dataStruct_biasedART(dataSetNum).test_output=dataStruct_biasedART(dataSetNum).test_output(:);
catch
    disp('Test Output Not Available')
end


inputNet = struct('M', {[]}, ... % Number of features
    'TestOnly', {0}, ... % Train & Test or Test only? unless specified otherwise
    'test_fast', {1}, ... % Use fast tester unless specified otherwise
    'compute_overlap', {0}, ... % Indirectly measure overlap by looking at how many nodes overlap with a test_point
    'C', {0}, ... % Number of committed coding nodes (paper says set to 0, but harder to Matlab)
    'plotSteps',{0}, ... % % Iterate through individual learning steps
    'NodeList', {[]}, ... % List of node associated with each input
    'maxNumCategories', {100}, ...
    'numClasses', {[]}, ... % The total number of classes
    'w', {[]}, ... % coding node weight vector
    'W', {[]}, ... % output class weight vector
    'rho_bar', {0}, ... % rho_bar, the baseline vigilance. baseline=0 maximizes code compression
    'rho', {0.0}, ... % vigilance
    'alpha', {.1}, ... % Signal Rule Parameter
    'alphaTrain',{[]},... % Store when reverting alpha value in TESTING
    'beta', {1.0}, ...  % Learning fraction
    'gamma',{.000001},... % fraction additive in numerator and denominator of matching signal
    'epsilon', {-.00001}, ... % Match Tracking (codes inconsistent cases)
    'base_vigilance', {0.0}, ... % Baseline vigilance
    'p', {1.0}, ... % CAM rule power
    'dataSubsets', {1}, ... % 4-fold cross-validation
    'votes', {5}, ... % Number of voting systems
    'learningRate', {1.0}, ...
    'win_sequence',{[]},...
    'act_sequence',{[]},...
    'search_cycles', {0},... % Total number of search cycles during training
    'lambda_Attention',{[]},...
    'dataStore',{[]},... % Store data if using structs as inputs
    'numEpochs',{1},...%Number of times to present input
    'model_num',{[]},... % Depletion Model Number
    'type_depletion',{4},... %Options are 0 1 2 3 4
    'deplete_fast',{[]},...
    'learn_trail',{[]},...
    'learn_description',{''},...
    'short_learn_descript',{''},...
    'uncommitted_fail',{0},...
    'NI_fail',{[]},...
    'NI_fail_addNode',{[]},...
    'node_created',{[]},...
    'start_conservative',{[]},...
    'ovlp_net',{[]},...
    'A_dep_dist',{[]},...
    'e_mod_dist',{[]},...
    'e_realLow',{0},...
    'e_zero',{0},...
    'diffE_count',{[]},...
    'perc_error',{[]},... % percentage error added to labels in training
    'var_error_features',{[]},... % variance of error added to features in training
    'testAccuracy',{[]}); % Attention gain parameter



% inputNet.lambda_Attention=0.5;

inputNet.lambda_Attention=lambda_Attention;
inputNet.deplete_fast = 0;  %TURN OFF (=0) IF DEPLETING USING ITERATIVE VERSION
inputNet.start_conservative =0;
inputNet.DEFAULT_RETEST=1;
inputNet.alpha=10^(-8);
inputNet.gamma=0.0000001;
inputNet.plotSteps=0;
inputNet.model_num=32;

%Train for number of epochs
inputNet.numEpochs=1;

%%%Inject label errors into dataset | e.g perc_error = 10 implies 10% of
%%%training data have erroneous class labels during training
perc_error=0;
inputNet.perc_error=perc_error;

%%%Inject feature noise into dataset | e.g var_error_features = 10 implies 10% of
%%%training data have noise added to input features
var_error_features=0;
inputNet.var_error_features=var_error_features;

inputNet.compute_overlap=1;

if (dataSetNum==13)
    %Matrix manipulation becomes too memory intensive for movie genres
    %benchmark
    inputNuet.test_fast=0;
end


if (inputNet.perc_error~=0)
    disp('Adding error in training set!');
elseif (inputNet.var_error_features~=0)
    disp('Adding feature error in training set!');
end
%%%%%%%%%%%%%%%%%%%
%Create uniform plot grid%%
%%%%%%%%%%%%%%%%%%%



if ismember(dataSetNum,[1 2 3 4 5 6 7])

    if (perc_error>0)
        orig_dataset_labels=dataStruct_biasedART(dataSetNum).training_output;
        dataStruct_biasedART(dataSetNum).training_output=addClassError(dataStruct_biasedART(dataSetNum).training_output,inputNet.perc_error);
        temp_labels_noise=dataStruct_biasedART(dataSetNum).training_output;
    elseif (var_error_features>0)
        orig_dataset_features=dataStruct_biasedART(dataSetNum).training_input;
        dataStruct_biasedART(dataSetNum).training_input=max(min(dataStruct_biasedART(dataSetNum).training_input+randn(size(dataStruct_biasedART(dataSetNum).training_input))*var_error_features,1),0);
        temp_features_noise=dataStruct_biasedART(dataSetNum).training_input;
    end
    if ~exist('dataTestUnif','var')
        load ARTMAP_data_test dataTestUnif dataTestOutUnif test_x test_y
    end
    dataStruct_biasedART(dataSetNum).test_input=dataTestUnif;
    if dataSetNum>1
        dataStruct_biasedART(dataSetNum).test_output=dataTestOutUnif{floor(dataSetNum/2)}';
    end
end






[a,b,c]=Default_ARTMAP(dataStruct_biasedART(dataSetNum),inputNet);
if ismember(dataSetNum,[1 2 3 4 5 6 7])
    if (perc_error>0)
        inputNet.lambda_Attention=0;
        [a_zero,b_zero,c_zero]=Default_ARTMAP(dataStruct_biasedART(dataSetNum),inputNet);
        disp('Same training order lambda = 0')
        disp(sum(a_zero==dataStruct_biasedART(dataSetNum).test_output)*100/length(a_zero));
        labels_switch_positions=find(orig_dataset_labels~=temp_labels_noise);
        disp(intersect(labels_switch_positions,c.NI_fail))
        dataStruct_biasedART(dataSetNum).training_output=orig_dataset_labels;
    elseif (var_error_features>0)
        dataStruct_biasedART(dataSetNum).training_input=orig_dataset_features;
    end
end
if (ismember(inputNet.model_num,[31 2^c.M]) && (inputNet.deplete_fast==1) && (inputNet.lambda_Attention==1))
    c.lambda_Attention='\infty';
end



if (c.M==2) && (inputNet.plotSteps==0)

    a=a(:);
    figure(1)
    if dataSetNum~=1
         colormap([1 .4 .4; 1 .8 .8;.8 .8 1; .4 .4 1])

        imagesc(test_x(:),test_y(:),(reshape(a+(a~=dataStruct_biasedART(dataSetNum).test_output).*((a==2)-(a==1)),length(test_x),length(test_x))));
          set(gca,'ydir','normal')
        hold on

    else
       
        colormap([1 .7 .7;.7 .7 1])
        imagesc(test_x(:),test_y(:),(reshape(a,length(test_x),length(test_x))));
        set(gca,'ydir','normal')
        hold on

    end
    drawRectsArtmap;

    axis equal
    axis tight
    axis([-.1 1.1 -.1 1.1]);
    title(['biased ARTMAP: feature selection in ART search cycle; Type (One Epoch)'],'fontsize',12);
    ylabel(['\lambda = ' num2str(c.lambda_Attention,3)],'fontsize',12)
    if dataSetNum~=1
        xlabel(['Accuracy% : ', num2str(sum(dataStruct_biasedART(dataSetNum).test_output==a)*100/length(dataStruct_biasedART(dataSetNum).test_output),3),', C Nodes: ',num2str(c.C), ...
            ', \alpha: ',num2str(c.alpha)......
            ', Training Pts : ',num2str(length(dataStruct_biasedART(dataSetNum).training_input))...
            ', U Fail : ',num2str(c.uncommitted_fail) ', NIT fail : ',num2str(length(c.NI_fail))],'fontsize',12);
    end


    if (perc_error>0)
        figure(2)
        if dataSetNum~=1
            
            colormap([1 .4 .4; 1 .8 .8;.8 .8 1; .4 .4 1])

           
            a_test=dataStruct_biasedART(dataSetNum).test_output;
            colormap([1 .8 .8;.8 .8 1])
            imagesc(test_x(:),test_y(:),(reshape(a_test+(a_test~=dataStruct_biasedART(dataSetNum).test_output).*((a_test==2)-(a_test==1)),length(test_x),length(test_x))));
            set(gca,'ydir','normal')
            hold on

        else
          
            colormap([1 .7 .7;.7 .7 1])
            imagesc(test_x(:),test_y(:),(reshape(a,length(test_x),length(test_x))));
            set(gca,'ydir','normal')
            hold on

        end
        for i=1:size(c_zero.w,2)
            if (sum(c_zero.w(:,i))==2)
                rect_created=[(min(c_zero.w(1,i),1-c_zero.w(3,i))) (min(c_zero.w(2,i),1-c_zero.w(4,i))) abs(1-c_zero.w(3,i)-c_zero.w(1,i))+.01 abs(1-c_zero.w(4,i)-c_zero.w(2,i))+.01];

                %IF RESCALING
                %rect_created=[(min((9*c.w(1,i)-4),1-(9*c.w(3,i)-4))) (min((9*c.w(2,i)-4),1-(9*c.w(4,i)-4))) abs(1-(9*c.w(3,i)-4)-(9*c.w(1,i)-4))+.01 abs(1-(9*c.w(4,i)-4)-(9*c.w(2,i)-4))+.01];

                if c_zero.W(i,1)==0
                    rectangle('Position',rect_created,'EdgeColor',[0 0 1],'LineWidth',1.5)
                    hold on
                    %text((rect_created(1)+rect_created(3)/2),(rect_created(2)+rect_created(4)/2),'\beta')
                else
                    rectangle('Position',rect_created,'EdgeColor',[1 0 0],'LineWidth',1.5)
                    hold on
                    %text((rect_created(1)+rect_created(3)/2),(rect_created(2)+rect_created(4)/2),'\alpha')
                end
            end
        end

        axis equal
        axis tight
        axis([-.1 1.1 -.1 1.1]);
        title(['biased ARTMAP: feature selection in ART search cycle; Type ' num2str(c_zero.model_num) ' (One Epoch)'],'fontsize',12);
        ylabel(['\lambda = ' num2str(c_zero.lambda_Attention,3)],'fontsize',12)
        if dataSetNum~=1
            xlabel(['Accuracy% : ', num2str(sum(dataStruct_biasedART(dataSetNum).test_output==a_zero)*100/length(dataStruct_biasedART(dataSetNum).test_output),3),', C Nodes: ',num2str(c_zero.C), ...
                ', \alpha: ',num2str(c_zero.alpha)......
                ', Training Pts : ',num2str(length(dataStruct_biasedART(dataSetNum).training_input))...
                ', U Fail : ',num2str(c_zero.uncommitted_fail) ', NIT fail : ',num2str(length(c_zero.NI_fail))],'fontsize',12);
        end
    end


    try
        depleted_count=0;
        for i=1:length(c.learn_trail)
            if ~isempty(c.learn_trail(dataSetNum).MT_nodes)
                depleted_count=depleted_count+1;
                depleted_sample(depleted_count)=i;
                depleted_num(depleted_count)=length(c.learn_trail(dataSetNum).MT_nodes);
                depleted_A(depleted_count,:)=c.learn_trail(dataSetNum).MT_Adash_after{end};
            end
        end
        depleted_sample=depleted_sample(:);
        depleted_num=depleted_num(:);
        min(sum(depleted_A,2))
    catch
        disp('Trail Not Available');
    end


elseif ismember(c.M,6)

    inputNet.lambda_Attention=0;

    num_mapping=c.W*[1 2]';
    max_vals=(b==repmat(max(b,[],2),1,3));
    ties=(mod(max_vals*(2*num_mapping+1),2)==0);


    inputNet.gamma=0.01;

    [a_nodep,b_nodep,c_nodep]=Default_ARTMAP(dataStruct_biasedART(dataSetNum),inputNet);
    test_out=num2str(dataStruct_biasedART(dataSetNum).test_input');
    test_out=test_out(:,1:3:size(test_out,2));
    num_mapping=c_nodep.W*[1 2]';
    max_vals=(b_nodep==repmat(max(b_nodep,[],2),1,3));
    ties_nodep=(mod(max_vals*(2*num_mapping+1),2)==0);
    show_val={'+','-','T'};
    a_with_ties=a.*(~ties)+3*ties;
    a_nodep_with_ties=a_nodep.*(~ties_nodep)+3*ties_nodep;

    [strvcat(num2str(round(b_nodep*1000)/1000,'%1.3f\t')) ...
        repmat(' ',2^(c.M),1)...
        repmat(' ',2^(c.M),1)...
        repmat(' ',2^(c.M),1)...
        repmat(' ',2^(c.M),1)...
        repmat(' ',2^(c.M),1)...
        repmat(' ',2^(c.M),1)...
        [show_val{a_nodep_with_ties}]'...
        repmat(' ',2^(c.M),1)...
        repmat(' ',2^(c.M),1)...
        repmat(' ',2^(c.M),1)...
        repmat(' ',2^(c.M),1)...
        repmat(' ',2^(c.M),1)...
        repmat(' ',2^(c.M),1)...
        [show_val{a_with_ties}]'...
        repmat(' ',2^(c.M),1)...
        repmat(' ',2^(c.M),1)...
        repmat(' ',2^(c.M),1)...
        repmat(' ',2^(c.M),1)...
        repmat(' ',2^(c.M),1)...
        repmat(' ',2^(c.M),1)...
        strvcat(num2str(round(b*1000)/1000,'%1.3f\t')) ...
        repmat(' ',2^(c.M),1)...
        repmat(' ',2^(c.M),1)...
        repmat(' ',2^(c.M),1)...
        repmat(' ',2^(c.M),1)...
        repmat(' ',2^(c.M),1)...
        repmat(' ',2^(c.M),1)...
        strvcat(test_out)]


    disp('Not Displaying Outputs: More than 2 input dimensions')

    disp('Total number of differences: Depletion versus no depletion')

    sum(a~=a_nodep)

elseif c.M==41

    disp('Finished running Boston dataset simulation')

    if dataSetNum<13
        title_boston=['Boston dataset ' num2str(dataSetNum-8) ' (uncapped)  \lambda=' num2str(c.lambda_Attention,2)];
    end

    [conf_nonnorm,conf_norm,avg_acc]= confMatrixPlot(dataStruct_biasedART(dataSetNum).test_output(:),a(:),c,dataStruct_biasedART(dataSetNum).labels,title_boston);


elseif strcmp(dataStruct_biasedART(dataSetNum).description,'movie_genres')
    b_mean=mean(b,1);
    b_std=std(b,1);
    b_new=(b-repmat(b_mean,length(dataStruct_biasedART(dataSetNum).test_output),1))./repmat(b_std,length(dataStruct_biasedART(dataSetNum).test_output),1);
    top_node_num=round(c.C*.02);
    dum=sort(b_new,2,'descend');
    y_rect=b_new.*(b_new>repmat(mean(b_new,2),1,c.C)).*(b_new>repmat(dum(:,top_node_num),1,c.C));
    sigma = c.W'*y_rect';
    sigma_mean=mean(sigma,2);
    sigma_std=std(sigma',1)';
    sigma_new=(sigma-repmat(sigma_mean,1,length(dataStruct_biasedART(dataSetNum).test_output)))./repmat(sigma_std,1,length(dataStruct_biasedART(dataSetNum).test_output));
    [dum,classes]=max(sigma,[],1);    
    
     [conf_nonnorm,conf_norm,avg_acc]= confMatrixPlot(dataStruct_biasedART(dataSetNum).test_output(:),classes,c,dataStruct_biasedART(dataSetNum).labels,'');
end

disp('Output class predictions are stored in variable ''a''.');
disp('Distributed output predictions are stored in ''b''.');
disp('The biased ARTMAP network details are stored in ''c''');

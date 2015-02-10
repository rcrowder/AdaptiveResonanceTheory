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


function updated_artmap = DEFARTMAP_LearnTwinV10dynamicPlot(artmap, inputVectors, classes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an implementation of Default ARTMAP,
% as described in CAS/CNS Technical Report TR-2003-008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programmed by Ogi Ogas
%%%%%%%%%%%%%%%%%%%%%%%%%
% The notation follows the article's notation,
% as does the headers for each step of the
% algorithm
% 02/21/05 Added search cycle count to artmap structure
% 08.31.05 Added list of nodes activated for each training input (This is
% only used for certain priming simulations; otherwise, this is ignored.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         WINNER-TAKE-ALL TRAINING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ON=1; OFF=0;

LEARN_TRAIL=ON;
PLOT_PARTITION=OFF;
if (artmap.deplete_fast==1)
    filesavedir=['ARTMAP_Attention_Trail_' artmap.dataStore.description '_' num2str(round(artmap.model_num)) '_' num2str(round(artmap.lambda_Attention*1000)) '_fast\' ];
elseif (artmap.deplete_fast==0)
    filesavedir=['ARTMAP_Attention_Trail_' artmap.dataStore.description '_' num2str(round(artmap.model_num)) '_' num2str(round(artmap.lambda_Attention*1000)) '_slow\' ];
else
    filesavedir=['ARTMAP_Attention_Trail_' artmap.dataStore.description '_' num2str(round(artmap.model_num)) '_' num2str(round(artmap.lambda_Attention*1000)) '\' ];
end
%filesavedir=['Attention_Feb_19_08'];
filename_prefix='output';

%GRAPH_DATA(artmap.dataStore.test_input,artmap.dataStore.test_output,artmap.dataStore.test_output);
axis equal
axis tight
axis([-.1 1.1 -.1 1.1]);
k_step=1;
save k_stepTrack k_step
MT_nodes=0;
space_amt=' ';


disp('Using Twin-Trail Catcher Version 16 (dyn) of DEFARTMAP_Learn!')

% Make sure that the data is appropriate for the given network.
[numFeatures, n] = size(inputVectors);
if (numFeatures ~= 2*artmap.M)
    error('The input vector does not contain the same number of features as the network.');
end



[artmap.learn_trail(1:n).A]=deal([]);
[artmap.learn_trail(1:n).K]=deal([]);
[artmap.learn_trail(1:n).winning_node]=deal([]);
[artmap.learn_trail(1:n).weight_before]=deal([]);
[artmap.learn_trail(1:n).weight_after]=deal([]);
[artmap.learn_trail(1:n).MT_nodes]=deal([]);
[artmap.learn_trail(1:n).MT_mismatch]=deal([]);
[artmap.learn_trail(1:n).MT_x_before]=deal([]);
[artmap.learn_trail(1:n).MT_e_before]=deal([]);
[artmap.learn_trail(1:n).MT_e_after]=deal([]);
[artmap.learn_trail(1:n).MT_Adash_before]=deal([]);
[artmap.learn_trail(1:n).MT_Adash_after]=deal([]);
artmap.compare_rho=zeros(n,3);
comp_rho=1;
% Initialize the return artmap network

%%%%%%%%%%%%for version x
if ismember(artmap.model_num,[0 1 19 20 21 22 23 24 25 26 28 29 30 31 32])
    e_AttentionDefault=zeros(2*artmap.M,1);
elseif ismember(artmap.model_num,[17 18])
    e_AttentionDefault=ones(2*artmap.M,1);
end
search_cycles = 0;

% Classify and learn on each sample.
for sampleNumber = 1:n

    % Reset the activated node (NOTE: this is Moshe Bar addition)
    activated_node=[];

    % Get the next input vector.
    % B.3 / B.6
    A = inputVectors(:, sampleNumber);

    MT_nodes=0;
    MT_ALL_nodes=0;

    if LEARN_TRAIL==ON
        artmap.learn_trail(sampleNumber).A=A;
    end

    % Get the next class
    % B.3 / B.6
    K = classes(sampleNumber);
    if LEARN_TRAIL==ON
        artmap.learn_trail(sampleNumber).K=K;
    end

    %     if sampleNumber==5
    %         disp('Begin Debug')
    %     end


    title(['Input Number: ' num2str(sampleNumber)]);
    %     if (ismember(artmap.model_num,[31 32]) && (artmap.deplete_fast==1) && (artmap.lambda_Attention==1))
    %     xlabel(['\lambda = \infty A: [',regexprep(num2str(round(A*100)'/100),'\s+',' '),'] '])
    %     else
    %     xlabel(['\lambda = ',num2str(round(artmap.lambda_Attention*1000)/1000),' A: [',regexprep(num2str(round(A*100)'/100),'\s+',' '),'] '])
    %     end
    plotInputPoint(A,filesavedir,filename_prefix,sampleNumber,K)

    % B.5
    artmap.rho = artmap.rho_bar;
    artmap.rho_dep = artmap.rho_bar;

    if (sampleNumber<15) && (artmap.start_conservative==1)
        artmap.rho = .98;
        artmap.rho_dep = .98;
    end

    if (artmap.C==0)
        artmap.C=1;
        artmap.w(:,1) = A;
        artmap.W = zeros(1,artmap.numClasses);
        artmap.W(1,K)=1;
        if PLOT_PARTITION==ON
            colormap([1 .7 .7;.7 .7 1])
            [labels_predict, dum]= DEFARTMAP_Classify_Fast(artmap, artmap.testVectors);
            h_back=imagesc(artmap.test_x(:),artmap.test_y(:),(reshape(labels_predict,length(artmap.test_x),length(artmap.test_x))));
            set(gca,'ydir','normal')
            hold on
        end

        coding_handle{artmap.C}=plotCodingNode(A,0,filesavedir,filename_prefix,sampleNumber,K);

        if LEARN_TRAIL==ON
            artmap.learn_trail(sampleNumber).winning_node=artmap.C;
            artmap.learn_trail(sampleNumber).weight_before = A';
            artmap.learn_trail(sampleNumber).weight_after = A';
        end
        activated_node=1; % The first and only node is activated.
    else
        T = DEFARTMAP_choice_function(artmap,A);
        [T_sorted, T_j] = sort(-T);
        T_j = T_j(logical(-T_sorted>artmap.alpha.*artmap.M));
        MATCH=0;

        e_AttentionBefore=e_AttentionDefault;
        A_dashBefore=A;

        for j=T_j


            x_sum=sum(min(max(artmap.w(:,j)-artmap.lambda_Attention*e_AttentionBefore,0),...
                max(A-artmap.lambda_Attention*e_AttentionBefore,0)));
            z_sum=sum(min(max(1-artmap.lambda_Attention*e_AttentionBefore,0),...
                max(A-artmap.lambda_Attention*e_AttentionBefore,0)));

            MATCH_TRACK_TEST=(((x_sum+artmap.gamma)/(sum(max(A-artmap.lambda_Attention*e_AttentionBefore,0))+eps+artmap.gamma)) >= artmap.rho);



            if LEARN_TRAIL==ON
                MT_ALL_nodes=MT_ALL_nodes+1;
                artmap.learn_trail(sampleNumber).rho_calc{MT_ALL_nodes}=[(x_sum+artmap.gamma)./(sum(A_dashBefore)+eps+artmap.gamma) ...
                    x_sum./artmap.M ];
                artmap.learn_trail(sampleNumber).rho{MT_ALL_nodes} = [artmap.rho_dep artmap.rho];
            end

            if exist('coding_handle','var')
                if (length(coding_handle)>=j && ~isempty(coding_handle{j}))
                    delete(coding_handle{j})
                end
            end

            coding_handle{j}=plotCodingNode(artmap.w(:,j),artmap.lambda_Attention*e_AttentionBefore,filesavedir,filename_prefix,sampleNumber,find(artmap.W(j,:)),'');


            if MATCH_TRACK_TEST
                search_cycles = search_cycles+1;
                J=find(artmap.W(j,:));
                if J==K
                    MATCH=1;
                    artmap.win_sequence(length(artmap.win_sequence)+1)=j;
                    if LEARN_TRAIL==ON
                        artmap.learn_trail(sampleNumber).winning_node=j;
                        artmap.learn_trail(sampleNumber).weight_before = artmap.w(:,j);
                    end

                    artmap.w=DEFARTMAP_Update_Weights(artmap,j,A);

                   

                    if LEARN_TRAIL==ON
                        artmap.learn_trail(sampleNumber).weight_after = artmap.w(:,j);
                    end

                    activated_node=j;
                    if exist('coding_handle','var')
                        if (length(coding_handle)>=j && ~isempty(coding_handle{j}))
                            delete(coding_handle{j})
                        end
                    end

                    if PLOT_PARTITION==ON
                        colormap([1 .7 .7;.7 .7 1])
                        [labels_predict, dum]= DEFARTMAP_Classify_Fast(artmap, artmap.testVectors);

                        set(h_back,'Cdata',reshape(labels_predict,length(artmap.test_x),length(artmap.test_x)));
                        set(gca,'ydir','normal')
                        hold on
                    end
                    coding_handle{j}=plotCodingNode(artmap.w(:,j),artmap.lambda_Attention*e_AttentionBefore,filesavedir,filename_prefix,sampleNumber,K);
                    break; % Don't check any more T_j
                else % if J~=K

                     artmap.rho = ((x_sum+artmap.gamma)./(sum(max(A-artmap.lambda_Attention*e_AttentionBefore,0))+eps+artmap.gamma)) + artmap.epsilon;
                    artmap.rho_uc = ((z_sum+artmap.gamma)./(sum(max(A-artmap.lambda_Attention*e_AttentionBefore,0))+eps+artmap.gamma));


                    [e_AttentionAfter,A_dashAfter,x_Attention] =...
                        ART_depleter_new(sampleNumber,artmap.w(:,j),A,e_AttentionBefore,artmap.lambda_Attention,artmap.M,artmap.model_num,artmap.deplete_fast);

                    %tmp=plotCodingNode(A_dashAfter,artmap.lambda_Attention*e_AttentionAfter,filesavedir,filename_prefix,sampleNumber);
                    %delete(tmp)

                    if sum(e_AttentionAfter)<.0001
                        artmap.e_realLow=artmap.e_realLow+1;
                    end
                    if sum(e_AttentionAfter)<2*eps
                        artmap.e_zero=artmap.e_zero+1;
                    end

                    if LEARN_TRAIL==ON
                        MT_nodes=MT_nodes+1;
                        artmap.learn_trail(sampleNumber).MT_nodes{MT_nodes}=j;
                        artmap.learn_trail(sampleNumber).MT_rho{MT_nodes} =artmap.rho;
                        artmap.learn_trail(sampleNumber).MT_x_before{MT_nodes} =x_Attention';
                        artmap.learn_trail(sampleNumber).MT_e_before{MT_nodes} =e_AttentionBefore';
                        artmap.learn_trail(sampleNumber).MT_e_after{MT_nodes} =e_AttentionAfter';
                        artmap.learn_trail(sampleNumber).MT_Adash_before{MT_nodes} =A_dashBefore';
                        artmap.learn_trail(sampleNumber).MT_Adash_after{MT_nodes} =A_dashAfter';
                    end

                    A_dashBefore=A_dashAfter;
                    e_AttentionBefore=e_AttentionAfter;


                end % if J==K b
            else % if (sum(min(A,artmap.w(:,j)))./artmap.M >= artmap.rho)

                

            end % if (sum(min(A,artmap.w(:,j)))./artmap.M >= artmap.rho)
        end % for J=T_j

        if (MATCH==0) %
            %sampleNumber
            artmap.C=artmap.C+1;
            artmap.w(:,artmap.C)=A;
            artmap.W=[artmap.W; zeros(1,artmap.numClasses)];
            artmap.W(artmap.C,K)=1;
            if PLOT_PARTITION==ON
                colormap([1 .7 .7;.7 .7 1])
                [labels_predict, dum]= DEFARTMAP_Classify_Fast(artmap, artmap.testVectors);
                set(h_back,'Cdata',reshape(labels_predict,length(artmap.test_x),length(artmap.test_x)));
                set(gca,'ydir','normal')
                hold on
            end
            coding_handle{artmap.C}=plotCodingNode(A,artmap.lambda_Attention*e_AttentionBefore,filesavedir,filename_prefix,sampleNumber,K);


            z_sum=sum(min(max(1-artmap.lambda_Attention*e_AttentionBefore,0),...
                max(A-artmap.lambda_Attention*e_AttentionBefore,0)));
            rho_uc = ((z_sum+artmap.gamma)./(sum(max(A-artmap.lambda_Attention*e_AttentionBefore,0))+eps+artmap.gamma));


            if artmap.rho>rho_uc
                artmap.uncommitted_fail=artmap.uncommitted_fail+1;
            end


            if LEARN_TRAIL==ON
                artmap.learn_trail(sampleNumber).winning_node=artmap.C;
                artmap.learn_trail(sampleNumber).weight_before = A';
                artmap.learn_trail(sampleNumber).weight_after = A';
            end



            activated_node=artmap.C; % Activated node is new node
        end % (MATCH==0)

    end      %  if (artmap.C==0)
    artmap.NodeList = [artmap.NodeList, activated_node];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%THIS LOOP CHECKS IF THE NEXT_INPUT TEST IS MET
    y = zeros(artmap.C,1);
    T = DEFARTMAP_choice_function(artmap,A);
    Lambda = find(T>artmap.alpha.*artmap.M);
    Lambda_prime = find(T==artmap.M);

    if ~isempty(Lambda_prime)  % Point Box Case -- w_j = A
        y(Lambda_prime)=(1./length(Lambda_prime));
    else
        sumL = sum(1./(artmap.M-T(Lambda)));
        y(Lambda)=(1./(artmap.M-T(Lambda)))./sumL;
    end

    class_val = find(artmap.W(y==max(y),:),1,'first');

    if K~=class_val
        artmap.NI_fail(end+1)=sampleNumber;

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end      % for sampleNumber = 1:numSamples

artmap.search_cycles = search_cycles;
updated_artmap = artmap;

return
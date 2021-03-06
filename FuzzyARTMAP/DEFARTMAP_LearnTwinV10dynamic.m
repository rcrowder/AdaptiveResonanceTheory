%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an implementation of an ARTMAP variant,
% as described in:
% Gaddam, C. S. (2007).
% Feature Selection via transmitter depletion in ARTMAP. Online Document,
% xx(x) xxx-xxx.
% Boston, MA: Boston University.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programmed by Chaitanya Sai (August 2007-08)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The notation follows the article's notation,
% as does the headers for each step of the
% algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function updated_artmap = DEFARTMAP_LearnTwinV10dynamic(artmap, inputVectors, classes)
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

LEARN_TRAIL=OFF;

disp('Using Twin-Trail Catcher Version 28 (dyn) of DEFARTMAP_Learn!')

% Make sure that the data is appropriate for the given network.
[numFeatures, n] = size(inputVectors);
if (numFeatures ~= 2*artmap.M)
    error('The input vector does not contain the same number of features as the network.');
end

if (LEARN_TRAIL==ON)
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
end
artmap.compare_rho=zeros(n,3);
comp_rho=1;
% Initialize the return artmap network

%%%%%%%%%%%%Initializing e before every round e should start at 1
%%%%%%%%%%%%for version x
if ismember(artmap.model_num,[0 1 19 20 21 22 23 24 25 26 28 29 30 31 32])
    e_AttentionDefault=zeros(2*artmap.M,1);
elseif ismember(artmap.model_num,[17 18])
    e_AttentionDefault=ones(2*artmap.M,1);
end

search_cycles = 0;
artmap.e_mod_dist=zeros(n,1);
artmap.A_dep_dist=zeros(n,1);
% Classify and learn on each sample.
for sampleNumber = 1:n


    if (mod(sampleNumber,400)==0)
        disp('Sample Number:')
        disp(sampleNumber')
    end

    %     if sampleNumber==6
    %         disp('Start Debugging Here')
    %     end

    %%%%%%%%%%%%%%%%%

    %     if sampleNumber==12268
    %         artmap.lambda_Attention=0;
    %     end

    %%%%%%%%%%%%%%%

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

    % B.5
    artmap.rho = artmap.rho_bar;
    artmap.rho_dep = artmap.rho_bar;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Try dynamically changing alpha

    %     artmap.alpha=10/(artmap.C+10);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if  (artmap.start_conservative==1)
        %artmap.rho = .98;
        %artmap.rho_dep = .98;
        artmap.alpha=10/(artmap.C+10);
    end

    if (artmap.C==0)
        artmap.C=1;
        artmap.w(:,1) = A;
        artmap.W = zeros(1,artmap.numClasses);
        artmap.W(1,K)=1;

        if LEARN_TRAIL==ON
            artmap.learn_trail(sampleNumber).winning_node=artmap.C;
            artmap.learn_trail(sampleNumber).weight_before = A';
            artmap.learn_trail(sampleNumber).weight_after = A';
        end

        y = zeros(1,artmap.C); % the code at F2. At each node j, activity reset to 0.
        activated_node=1; % The first and only node is activated.
    else
        T = DEFARTMAP_choice_function(artmap,A);
        [T_sorted, T_j] = sort(-T);
        T_j = T_j(logical(-T_sorted>artmap.alpha.*artmap.M*(3.5)));
        MATCH=0;

        e_AttentionBefore=e_AttentionDefault;
        A_dashBefore=A;

        for j=T_j

            x_sum=sum(min(A_dashBefore,artmap.w(:,j)));
            z_sum=sum(min(A_dashBefore,1));
            if sum(artmap.model_num==[24 25])==1
                x_sum = sum(max(min(A_dashBefore,artmap.w(:,j))-artmap.lambda_Attention*e_AttentionBefore,0));
                z_sum = sum(max(min(A_dashBefore,1)-artmap.lambda_Attention*e_AttentionBefore,0));
            elseif sum(artmap.model_num==26)==1
                x_sum=max(sum(min(A_dashBefore,artmap.w(:,j)))-(artmap.M-sum(A_dashBefore)),0);
                z_sum=max(sum(min(A_dashBefore,1))-(artmap.M-sum(A_dashBefore)),0);

            elseif sum(artmap.model_num==27)==1
                x_sum =  sum(max(min(A_dashBefore,artmap.w(:,j))-((2*artmap.M-sum(artmap.w(:,j)))/(2*artmap.M))*artmap.lambda_Attention*e_AttentionBefore,0));
                z_sum =  sum(max(min(A_dashBefore,1)-((2*artmap.M-2*artmap.M)/(2*artmap.M))*artmap.lambda_Attention*e_AttentionBefore,0));
            elseif sum(artmap.model_num==28)==1
                x_sum=sum(max(min(A_dashBefore,artmap.w(:,j))-(A-A_dashBefore),0));
                z_sum=sum(max(min(A_dashBefore,1)-(A-A_dashBefore),0));
            elseif sum(artmap.model_num==29)==1
                x_sum=sum(min(max(artmap.w(:,j)-artmap.lambda_Attention*e_AttentionBefore,0),...
                    max(A_dashBefore-artmap.lambda_Attention*e_AttentionBefore,0)));
                z_sum=sum(min(max(1-artmap.lambda_Attention*e_AttentionBefore,0),...
                    max(A_dashBefore-artmap.lambda_Attention*e_AttentionBefore,0)));
            elseif sum(artmap.model_num==[30 31 32])==1
                x_sum=sum(min(max(artmap.w(:,j)-artmap.lambda_Attention*e_AttentionBefore,0),...
                    max(A-artmap.lambda_Attention*e_AttentionBefore,0)));
                z_sum=sum(min(max(1-artmap.lambda_Attention*e_AttentionBefore,0),...
                    max(A-artmap.lambda_Attention*e_AttentionBefore,0)));
            end




            if (artmap.use_type_three_tracking==1)


                if (sum(artmap.model_num==28)==1)
                    MATCH_TRACK_TEST=(((x_sum+artmap.gamma)/(sum(max(2*A_dashBefore-A,0))+eps+artmap.gamma)) >= artmap.rho);
                elseif (sum(artmap.model_num==29)==1)
                    MATCH_TRACK_TEST=(((x_sum+artmap.gamma)/(sum(max(A_dashBefore-artmap.lambda_Attention*e_AttentionBefore,0))+eps+artmap.gamma)) >= artmap.rho);
                elseif (sum(artmap.model_num==[30 31 32])==1)
                    MATCH_TRACK_TEST=(((x_sum+artmap.gamma)/(sum(max(A-artmap.lambda_Attention*e_AttentionBefore,0))+eps+artmap.gamma)) >= artmap.rho);
                else
                    MATCH_TRACK_TEST=(((x_sum+artmap.gamma)/(sum(A_dashBefore)+eps+artmap.gamma)) >= artmap.rho);
                end

                %MATCH_FACTOR=(sum(A_dashBefore)-sum(min(A_dashBefore,artmap.w(:,j))))/((sum(A)-sum(min(A,artmap.w(:,j))))+eps);
                %MATCH_TRACK_TEST=((sum(min(A,artmap.w(:,j)))*MATCH_FACTOR)/artmap.M >= artmap.rho);

            elseif (artmap.use_type_three_tracking==0)
                MATCH_TRACK_TEST=(x_sum/artmap.M >= artmap.rho);
            elseif (artmap.use_type_three_tracking==2)


                MATCH_TRACK_1= (((x_sum+artmap.gamma)/(sum(A_dashBefore)+eps+artmap.gamma)) >= artmap.rho_dep);
                MATCH_TRACK_2= (x_sum/artmap.M >= artmap.rho);
                artmap.compare_rho(comp_rho,:)=[MATCH_TRACK_1 MATCH_TRACK_2 sampleNumber];
                comp_rho=comp_rho+1;

                MATCH_TRACK_TEST = ((((x_sum+artmap.gamma)/(sum(A_dashBefore)+eps+artmap.gamma)) >= artmap.rho_dep) ...
                    && (x_sum/artmap.M >= artmap.rho));

            elseif (artmap.use_type_three_tracking==3)
                MATCH_TRACK_TEST=((x_sum+artmap.gamma)/...
                    ((1-sum(artmap.w(:,j))/(2*artmap.M))*artmap.M+sum(artmap.w(:,j))*sum(A_dashBefore)/(2*artmap.M)+artmap.gamma))>=artmap.rho;
            end



            if LEARN_TRAIL==ON
                MT_ALL_nodes=MT_ALL_nodes+1;
                artmap.learn_trail(sampleNumber).rho_calc{MT_ALL_nodes}=[(x_sum+artmap.gamma)./(sum(A_dashBefore)+eps+artmap.gamma) ...
                    x_sum./artmap.M ];
                artmap.learn_trail(sampleNumber).rho{MT_ALL_nodes} = [artmap.rho_dep artmap.rho];
            end


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

                    storeWeights=artmap.w;
                    artmap.w=DEFARTMAP_Update_Weights(artmap,j,A);
                    artmap.e_mod_dist(sampleNumber)=artmap.lambda_Attention*sum(e_AttentionBefore)/2;
                    artmap.A_dep_dist(sampleNumber)=artmap.lambda_Attention*sum(A_dashBefore);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %                     %%THIS LOOP CHECKS IF THE NEXT_INPUT TEST IS MET
                    %                     y = zeros(artmap.C,1);
                    %                     T = DEFARTMAP_choice_function(artmap,A);
                    %                     Lambda = find(T>artmap.alpha.*artmap.M);
                    %                     Lambda_prime = find(T==artmap.M);
                    %
                    %                     if ~isempty(Lambda_prime)  % Point Box Case -- w_j = A
                    %                         y(Lambda_prime)=(1./length(Lambda_prime));
                    %                     else
                    %                         sumL = sum(1./(artmap.M-T(Lambda)));
                    %                         y(Lambda)=(1./(artmap.M-T(Lambda)))./sumL;
                    %                     end
                    %
                    %                     %                     %%%%%%%%%%%%%%%%DEFAULT TEST%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %                                         if artmap.DEFAULT_RETEST
                    %                                             class_val = find(artmap.W(y==max(y),:),1,'first');
                    %                                         else
                    %                                             %%%%%%%%%%%DISTRIBUTED TEST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %                                             [dum,class_val] = max(y'*artmap.W);
                    %                                         end
                    %
                    %                                         if K~=class_val
                    %                                             artmap.NI_fail_addNode(end+1)=sampleNumber;
                    %                     %
                    %                     %                         %%%%%%%%%%%%UNDO RECENT LEARNING AND ADD NODE%%%%%%%%%%%
                    %                     %                         artmap.w=storeWeights;
                    %
                    %
                    %                                             artmap.C=artmap.C+1;
                    %                                             artmap.w(:,artmap.C)=A;
                    %
                    %                                             if LEARN_TRAIL==ON
                    %                                                 artmap.learn_trail(sampleNumber).winning_node=artmap.C;
                    %                                                 artmap.learn_trail(sampleNumber).weight_before = A';
                    %                                                 artmap.learn_trail(sampleNumber).weight_after = A';
                    %                                             end
                    %
                    %
                    %                                             artmap.W=[artmap.W; zeros(1,artmap.numClasses)];
                    %                                             artmap.W(artmap.C,K)=1;
                    %
                    %                                             % Activated node is new node
                    %
                    %                                         end


                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




                    if LEARN_TRAIL==ON
                        artmap.learn_trail(sampleNumber).weight_after = artmap.w(:,j);
                    end

                    activated_node=j;
                    break; % Don't check any more T_j
                else % if J~=K

                    if (artmap.use_type_three_tracking==1)

                        if sum(artmap.model_num==28)==1
                            %%%If Model Number 28
                            artmap.rho = ((x_sum+artmap.gamma)./(sum(max(2*A_dashBefore-A,0))+eps+artmap.gamma)) + artmap.epsilon;
                            artmap.rho_uc = ((z_sum+artmap.gamma)./(sum(max(2*A_dashBefore-A,0))+eps+artmap.gamma));
                        elseif sum(artmap.model_num==29)==1
                            artmap.rho = ((x_sum+artmap.gamma)./(sum(max(A_dashBefore-artmap.lambda_Attention*e_AttentionBefore,0))+eps+artmap.gamma)) + artmap.epsilon;
                            artmap.rho_uc = ((z_sum+artmap.gamma)./(sum(max(A_dashBefore-artmap.lambda_Attention*e_AttentionBefore,0))+eps+artmap.gamma));
                        elseif sum(artmap.model_num==[30 31 32])==1
                            artmap.rho = ((x_sum+artmap.gamma)./(sum(max(A-artmap.lambda_Attention*e_AttentionBefore,0))+eps+artmap.gamma)) + artmap.epsilon;
                            artmap.rho_uc = ((z_sum+artmap.gamma)./(sum(max(A-artmap.lambda_Attention*e_AttentionBefore,0))+eps+artmap.gamma));
                        else
                            %%%If Model Number is not 28, 29, or 30
                            artmap.rho = ((x_sum+artmap.gamma)./(sum(A_dashBefore)+eps+artmap.gamma)) + artmap.epsilon;
                            artmap.rho_uc = ((z_sum+artmap.gamma)./(sum(A_dashBefore)+eps+artmap.gamma));
                        end


                        %MATCH_FACTOR=(sum(A_dashBefore)-sum(min(A_dashBefore,artmap.w(:,j))))/((sum(A)-sum(min(A,artmap.w(:,j))))+eps);
                        %artmap.rho=(sum(min(A,artmap.w(:,j)))*MATCH_FACTOR)/artmap.M;

                    elseif (artmap.use_type_three_tracking==0)
                        artmap.rho = x_sum./artmap.M + artmap.epsilon;
                        artmap.rho_uc = z_sum./artmap.M ;
                    elseif (artmap.use_type_three_tracking==2)
                        artmap.rho = x_sum./artmap.M + artmap.epsilon;
                        artmap.rho_uc = z_sum./artmap.M;
                        artmap.rho_dep = ((x_sum+artmap.gamma)./(sum(A_dashBefore)+eps+artmap.gamma)) + artmap.epsilon;
                        artmap.rho_dep_uc = ((z_sum+artmap.gamma)./(sum(A_dashBefore)+eps+artmap.gamma));
                        artmap.rho=min(artmap.rho_dep,artmap.rho);
                        artmap.rho_uc=min(artmap.rho_dep_uc,artmap.rho_uc);
                        artmap.rho_dep=artmap.rho;

                    elseif (artmap.use_type_three_tracking==3)
                        %THIS SEEMS TO BE WRONG%
                        artmap.rho = ((x_sum+artmap.gamma)./...
                            ((1-sum(artmap.w(:,j))/(2*artmap.M))*artmap.M+sum(artmap.w(:,j))*sum(A_dashBefore)/(2*artmap.M)+artmap.gamma));
                    end



                    [e_AttentionAfter,A_dashAfter,x_Attention] =...
                        ART_depleter_new(sampleNumber,artmap.w(:,j),A,e_AttentionBefore,artmap.lambda_Attention,artmap.M,artmap.model_num,artmap.deplete_fast);


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
            else % if (sum(min(A,artmap.w(:,j)))./artmap.M >= artmap.rho)  %RESET, NO MATCH TRACKING

                if (artmap.deplete_in_mismatch==1)


                    [e_AttentionAfter,A_dashAfter,x_Attention] =...
                        ART_depleter_new(sampleNumber,artmap.w(:,j),A,e_AttentionBefore,artmap.lambda_Attention,artmap.M,artmap.model_num,artmap.deplete_fast);


                    if sum(e_AttentionAfter)<.0001
                        artmap.e_realLow=artmap.e_realLow+1;
                    end
                    if sum(e_AttentionAfter)<2*eps
                        artmap.e_zero=artmap.e_zero+1;
                    end


                    if LEARN_TRAIL==ON
                        MT_nodes=MT_nodes+1;
                        artmap.learn_trail(sampleNumber).MT_mismatch{MT_nodes}=1;
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


                end

            end % if (sum(min(A,artmap.w(:,j)))./artmap.M >= artmap.rho)
        end % for J=T_j

        if (MATCH==0) %
            %sampleNumber
            artmap.C=artmap.C+1;
            artmap.w(:,artmap.C)=A;
            artmap.node_created=[artmap.node_created sampleNumber];
            artmap.e_mod_dist(sampleNumber)=artmap.lambda_Attention*sum(e_AttentionBefore)/2;
            artmap.A_dep_dist(sampleNumber)=artmap.lambda_Attention*sum(A_dashBefore);

            if (artmap.use_type_three_tracking==1)

                if sum(artmap.model_num==28)==1
                    z_sum=sum(max(min(A_dashBefore,1)-(A-A_dashBefore),0));
                    rho_uc = ((z_sum+artmap.gamma)./(sum(max(2*A_dashBefore-A,0))+eps+artmap.gamma))+0*artmap.epsilon;
                elseif sum(artmap.model_num==29)==1
                    z_sum=sum(min(max(1-artmap.lambda_Attention*e_AttentionBefore,0),...
                        max(A_dashBefore-artmap.lambda_Attention*e_AttentionBefore,0)));
                    rho_uc = ((z_sum+artmap.gamma)./(sum(max(A_dashBefore-artmap.lambda_Attention*e_AttentionBefore,0))+eps+artmap.gamma));
                    %rho_uc = ((z_sum+artmap.gamma)./(sum(max(2*A_dashBefore-A,0))+eps+artmap.gamma))+0*artmap.epsilon;
                elseif sum(artmap.model_num==[30 31 32])==1
                    z_sum=sum(min(max(1-artmap.lambda_Attention*e_AttentionBefore,0),...
                        max(A-artmap.lambda_Attention*e_AttentionBefore,0)));
                    rho_uc = ((z_sum+artmap.gamma)./(sum(max(A-artmap.lambda_Attention*e_AttentionBefore,0))+eps+artmap.gamma));
                    %rho_uc =
                    %((z_sum+artmap.gamma)./(sum(max(2*A_dashBefore-A,0))+eps+artmap.gamma))+0*artmap.epsilon;
                elseif sum(artmap.model_num==[24 25])==1
                    z_sum = sum(max(min(A_dashBefore,1)-artmap.lambda_Attention*e_AttentionBefore,0));
                    rho_uc = ((z_sum+artmap.gamma)./(sum(A_dashBefore)+eps+artmap.gamma));
                elseif sum(artmap.model_num==26)==1
                    z_sum=max(sum(min(A_dashBefore,1))-(artmap.M-sum(A_dashBefore)),0);
                    rho_uc = ((z_sum+artmap.gamma)./(sum(A_dashBefore)+eps+artmap.gamma));

                elseif sum(artmap.model_num==27)==1
                    z_sum =  sum(max(min(A_dashBefore,1)-((2*artmap.M-2*artmap.M)/(2*artmap.M))*artmap.lambda_Attention*e_AttentionBefore,0));
                    rho_uc = ((z_sum+artmap.gamma)./(sum(A_dashBefore)+eps+artmap.gamma));
                end



            end


            if artmap.rho>(rho_uc+.00001)
                artmap.uncommitted_fail=artmap.uncommitted_fail+1;
                disp('It is happening!');
            end



            if LEARN_TRAIL==ON
                artmap.learn_trail(sampleNumber).winning_node=artmap.C;
                artmap.learn_trail(sampleNumber).weight_before = A';
                artmap.learn_trail(sampleNumber).weight_after = A';
            end


            artmap.W=[artmap.W; zeros(1,artmap.numClasses)];
            artmap.W(artmap.C,K)=1;

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
% if (sampleNumber<=20)    
% drawRectsARTMAP_internal
% end
end      % for sampleNumber = 1:numSamples

artmap.search_cycles = search_cycles;
updated_artmap = artmap;

return
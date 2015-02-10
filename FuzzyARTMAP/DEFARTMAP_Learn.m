function updated_artmap = DEFARTMAP_Learn(artmap, inputVectors, classes)
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
DISPLAY_WEIGHTBOXES=OFF;
DISPLAY_ONLYFINALWEIGHTBOX=OFF;
VISUALS=OFF;
DESCRIPT=OFF;
SHORT_DESCRIPT=OFF;
LEARN_TRAIL=ON;

disp('Using Original Bloated Version!!')

MT_nodes=0;
space_amt=' ';


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
[artmap.learn_trail(1:n).MT_x_after]=deal([]);
[artmap.learn_trail(1:n).MT_e_before]=deal([]);
[artmap.learn_trail(1:n).MT_e_after]=deal([]);
[artmap.learn_trail(1:n).MT_Adash_before]=deal([]);
[artmap.learn_trail(1:n).MT_Adash_after]=deal([]);
% Initialize the return artmap network
updated_artmap = {};

A=[]; K=[];
search_cycles = 0;

% Classify and learn on each sample.
for sampleNumber = 1:n
    % Reset the activated node (NOTE: this is Moshe Bar addition)
    activated_node=[];
    
    % Get the next input vector.
    % B.3 / B.6
    A = inputVectors(:, sampleNumber);
    MT_nodes=0;

    
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
    
    if (DESCRIPT==ON) 
        artmap.learn_description=strvcat(artmap.learn_description,['Sample number: ',num2str(sampleNumber) ' =======================================================']);
        artmap.learn_description=strvcat(artmap.learn_description,['Input presented: ',regexprep(num2str(A(:)',3),'\s*',' ')]);
        artmap.learn_description=strvcat(artmap.learn_description,['Category of new input: ',num2str(K)]);
    end
    
    
    
    % Create a new category if this supervisory signal
    % has never been seen before.
    if (artmap.C==0)
        % B.4 Set initial weights for newly committed coding node j=C=1
        artmap.C=1;
        j=1;
        artmap.w(:,1) = A;
        artmap.W = zeros(1,artmap.numClasses);
        artmap.W(1,K)=1;
        
        if LEARN_TRAIL==ON
            artmap.learn_trail(sampleNumber).winning_node=artmap.C;
            artmap.learn_trail(sampleNumber).weight_before = A';
            artmap.learn_trail(sampleNumber).weight_after = A';
        end
        
        
        if (DESCRIPT==ON)
            artmap.learn_description=strvcat(artmap.learn_description,['This is the first category box. Weights are set to: ', regexprep(num2str(artmap.w(:,1)',3) ,'\s*',' ')]);
            artmap.learn_description=strvcat(artmap.learn_description,'First Point--no category boxes!'); 
            artmap.learn_description=strvcat(artmap.learn_description,['Adding New Node']);
            artmap.learn_description=strvcat(artmap.learn_description,['Input: ', regexprep(num2str(A',3),'\s*',space_amt)]);
            artmap.learn_description=strvcat(artmap.learn_description,['Weight Node Number: ', num2str(artmap.C)]);
            artmap.learn_description=strvcat(artmap.learn_description,['Newly Added Weight: ', regexprep(num2str(artmap.w(:,artmap.C)',3),'\s*',space_amt)]);
        end
        
        if (SHORT_DESCRIPT==ON)
            artmap.short_learn_descript=strvcat(artmap.short_learn_descript,['Adding New Node']);
            artmap.short_learn_descript=strvcat(artmap.short_learn_descript,['Input: ', regexprep(num2str(A',3),'\s*',space_amt)]);
            artmap.short_learn_descript=strvcat(artmap.short_learn_descript,['Weight Node Number: ', num2str(artmap.C)]);
            artmap.short_learn_descript=strvcat(artmap.short_learn_descript,['Newly Added Weight: ', regexprep(num2str(artmap.w(:,artmap.C)',3),'\s*',space_amt)]);
        end
        
        
        y = zeros(1,artmap.C); % the code at F2. At each node j, activity reset to 0.
        activated_node=1; % The first and only node is activated.
    else
        % B.7 Calculate signals to committed coding nodes
        T = DEFARTMAP_choice_function(artmap,A);
        
        
        % B.8 Sort committed coding node signals (A passing through
        % w_ij)
        % This next line sorts Tj in decreasing order 
        [T_sorted, T_j] = sort(-T);
        % This next line orders the coding nodes j with Tj>alpha*M, still in decreasing
        % order. We will need j itself (remember, j is an integer
        % from 1 to C) for next step.
        T_j = T_j(find(-T_sorted>artmap.alpha.*artmap.M));
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % B.9 Find appropriate coding node (you're only going to pick one
        % coding node j, the best j for the input vector A).
        % Go down coding nodes, as ordered by the sorted T_j
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % B.9a Take next coding node j=J (from sorted T_j), and see if
        % there is a suitable match.
        MATCH=0;
        
        %These additions are specific to ARTMAP Attention
        %=========ARTMAP Attention=======
        %attention parameter resets to zero at every input presentation
        e_Attention=zeros(size(A));
        A_dash=A;
        %================================
        
        for j=T_j
            
            if (artmap.use_type_three_tracking==1)
                %%%%%%%%%%Uncommitted nodes will always pass this test%%%%%%%%%%%%%% 
                comparing_to_rho = sum(min(A_dash,artmap.w(:,j)))./(sum(A_dash)+eps);
                artmap.num_type_three_used=artmap.num_type_three_used+(sum(A_dash)~=sum(A));
            else    
                comparing_to_rho = sum(min(A_dash,artmap.w(:,j)))./artmap.M;
            end
            
            
            if (DESCRIPT==ON) 
                artmap.learn_description=strvcat(artmap.learn_description,['Computing match between depleted input A and weight ' num2str(j) ': ' regexprep( num2str(artmap.w(:,j)',3) ,'\s*',space_amt)]);    
                artmap.learn_description=strvcat(artmap.learn_description,['Original input A: ', regexprep(num2str(A',3),'\s*',space_amt)]);
                artmap.learn_description=strvcat(artmap.learn_description,['Depleted input A: ', regexprep(num2str(A_dash',3),'\s*',space_amt)]);
                artmap.learn_description=strvcat(artmap.learn_description,['Match value sum(min(A_depleted , w_j)/M) rho: ',num2str(comparing_to_rho)]);
            end
            
            
            
            %%%%%%%%%%%Match Tracking: Change for Type III%%%%%%%%%%%%%%
            if (artmap.use_type_three_tracking==1)
                MATCH_TRACK_TEST=(sum(min(A_dash,artmap.w(:,j)))./(sum(A_dash)+eps) >= artmap.rho);
            else
                MATCH_TRACK_TEST=(sum(min(A_dash,artmap.w(:,j)))./artmap.M >= artmap.rho);
            end
            
            
            if LEARN_TRAIL==ON
                artmap.learn_trail(sampleNumber).rho_calc=sum(min(A_dash,artmap.w(:,j)))./(sum(A_dash)+eps);
                artmap.learn_trail(sampleNumber).rho = artmap.rho;
            end
            
            
            if MATCH_TRACK_TEST 
                search_cycles = search_cycles+1;
                % Then set Yj=1: Winner-take-all
                % B.9b (Not necessary with WTA learning)
                % B.9c CORRECT PREDICTION?
                J=find(artmap.W(j,:));
                if J==K
                    % if sigma_k = W_JK = 1
                    MATCH=1; 
                    % B.11 LEARNING: UPDATING WEIGHTS OF EXISTING
                    % CODING NODE IN F2
                    artmap.win_sequence(length(artmap.win_sequence)+1)=j;
                    if (DESCRIPT==ON)                
                        artmap.learn_description=strvcat(artmap.learn_description,'Correct match between input and the (next) nearest category box.');
                        artmap.learn_description=strvcat(artmap.learn_description,['Changing weight Node Number: ', num2str(j)]);
                        artmap.learn_description=strvcat(artmap.learn_description,['Weight for this category box before modification: ', regexprep(num2str( artmap.w(:,j)',3),'\s*',' ')]);
                    end
                    
                    
                    
                    if (SHORT_DESCRIPT==ON)
                        artmap.short_learn_descript=strvcat(artmap.short_learn_descript,['Changing weight']);
                        artmap.short_learn_descript=strvcat(artmap.short_learn_descript,['Input: ', regexprep(num2str(A',3),'\s*',space_amt)]);
                        artmap.short_learn_descript=strvcat(artmap.short_learn_descript,['Weight node number: ', num2str(j)]);
                        artmap.short_learn_descript=strvcat(artmap.short_learn_descript,['Weight before modification: ', regexprep( num2str(artmap.w(:,j)',3) ,'\s*',space_amt)]);
                    end
                    
                    if LEARN_TRAIL==ON
                        artmap.learn_trail(sampleNumber).winning_node=j;
                        artmap.learn_trail(sampleNumber).weight_before = artmap.w(:,j);
                    end
                    
                    artmap.w=DEFARTMAP_Update_Weights(artmap,j,A);
                    
                    if LEARN_TRAIL==ON
                        artmap.learn_trail(sampleNumber).weight_after = artmap.w(:,j);
                    end
                    
                    if (SHORT_DESCRIPT==ON)
                        artmap.short_learn_descript=strvcat(artmap.short_learn_descript,['Weight after modification: ', regexprep( num2str(artmap.w(:,j)',3) ,'\s*',space_amt)]);
                    end
                    
                    if (DESCRIPT==ON)
                        artmap.learn_description=strvcat(artmap.learn_description,['Weight after modification: ', regexprep( num2str(artmap.w(:,j)',3) ,'\s*',space_amt)]);
                    end
                    
                    
                    if (DESCRIPT==ON)
                        
                        %disp('the updated weights for this category box are now:')
                        %artmap.w(:,j)
                        %disp('Weight boxes: one should be stretched, or holds new point')
                        
                        artmap.learn_description=strvcat(artmap.learn_description,'Finished Learning on this input; moving to next input');
                    end
                    if (DISPLAY_WEIGHTBOXES==ON)
                        GRAPH_WEIGHTS(artmap,A,K)
                    end
                    
                    
                    
                    activated_node=j;
                    break; % Don't check any more T_j
                else % if J~=K
                    % if sigma_k = 0
                    % B.9d MATCH TRACKING: RAISE VIGILANCE
                    %%%%%%%%%%%Match Tracking: Change for Type III%%%%%%%%%%%%%%
                    
                    if LEARN_TRAIL==ON
                        MT_nodes=MT_nodes+1;
                        artmap.learn_trail(sampleNumber).MT_mismatch{MT_nodes}=1;
                        artmap.learn_trail(sampleNumber).MT_nodes{MT_nodes}=j;
                        artmap.learn_trail(sampleNumber).MT_rho_before{MT_nodes} =artmap.rho;
                    end
                    
                    
                    
                    if (artmap.use_type_three_tracking==1)
                        artmap.rho = sum(min(A_dash,artmap.w(:,j)))./(sum(A_dash)+eps) + artmap.epsilon; 
                    else
                        artmap.rho = sum(min(A_dash,artmap.w(:,j)))./artmap.M + artmap.epsilon;
                    end
                    
                    
                    if LEARN_TRAIL==ON
                        artmap.learn_trail(sampleNumber).MT_rho_after{MT_nodes} =artmap.rho;
                    end
                    
                    
                    if (DESCRIPT==ON)
                        %disp('INCORRECT match. Raise vigilance. Go to next closest category box (This is called match tracking, not really vigilance!')
                        artmap.learn_description=strvcat(artmap.learn_description,['INCORRECT match; match tracking and updating rho: ', num2str(artmap.rho)]);
                    end
                    
                    %                     %=========ARTMAP Attention=======
                    
                    %This is without normalization of any sort
                    x_Attention=min(A_dash,artmap.w(:,j));
                    
                    if LEARN_TRAIL==ON
                        artmap.learn_trail(sampleNumber).MT_x_before{MT_nodes} =x_Attention';
                    end
                    
                    
                    
                    %This is with normalization
                    %======================================================
                    %=============
                    if artmap.type_depletion==4
                        
                        if (DESCRIPT==ON)
                            artmap.learn_description=strvcat(artmap.learn_description,['Match Tracking triggers depletion; dyad before normalization: ', regexprep(num2str(x_Attention(:)',3),'\s*',' ')]);
                        end
                        
                        if (artmap.deplete_in_X_denom==1)
                            x_Attention=max(x_Attention+[x_Attention(artmap.M+1:end,1) ; x_Attention(1:artmap.M,1)]-sum(x_Attention)/(sum(A_dash)+eps),0);
                        else
                            x_Attention=max(x_Attention+[x_Attention(artmap.M+1:end,1) ; x_Attention(1:artmap.M,1)]-sum(x_Attention)/artmap.M,0);
                        end
                        
                        
                        if (DESCRIPT==ON)
                            
                            artmap.learn_description=strvcat(artmap.learn_description,['Match Tracking triggers depletion; dyad after normalization: ', regexprep(num2str(x_Attention(:)',3),'\s*',' ')]);
                        end
                    end
                    
                    if LEARN_TRAIL==ON
                        artmap.learn_trail(sampleNumber).MT_x_after{MT_nodes} =x_Attention';
                    end
                    
                    if LEARN_TRAIL==ON
                        artmap.learn_trail(sampleNumber).MT_e_before{MT_nodes} =e_Attention';
                    end
                    
                    
                    e_Attention=max(e_Attention,artmap.lambda_Attention*x_Attention);
                    
                    if LEARN_TRAIL==ON
                        artmap.learn_trail(sampleNumber).MT_e_after{MT_nodes} =e_Attention';
                    end
                    
                    
                    
                    %Attention can be too harsh if it is depleting along
                    %many dimension. To offset this, we could scale
                    %according to
                    
                    if artmap.type_depletion==1
                        %A) the number of dimensions that are non-zero
                        %==========
                        e_Attention=e_Attention/max(sum(sign02(e_Attention)),1);
                        %B) The sum total of attention, implying that attention
                        %is a constant
                    elseif artmap.type_depletion==2
                        e_Attention=e_Attention/(sum(e_Attention)+eps)*artmap.lambda_Attention;
                        %C) The attention graph is used to comparatively scale things implying that attention
                        %is a constant
                    elseif artmap.type_depletion==3
                        e_Attention=e_Attention/(sum(e_Attention)+eps).*e_Attention;
                    end
                    
                    if LEARN_TRAIL==ON
                        artmap.learn_trail(sampleNumber).MT_Adash_before{MT_nodes} =A_dash';
                    end
                    
                    
                    
                    A_dash=rect(A-e_Attention,0);
                    %                     %================================
                    
                    if LEARN_TRAIL==ON
                        artmap.learn_trail(sampleNumber).MT_Adash_after{MT_nodes} =A_dash';
                    end
                    
                    
                    
                    if (DESCRIPT==ON)
                        
                        artmap.learn_description=strvcat(artmap.learn_description,['Match Tracking triggers depletion. e_Attention: ', regexprep(num2str(e_Attention(:)',3),'\s*',' ')]);
                        artmap.learn_description=strvcat(artmap.learn_description,['Depleted input: ', regexprep(num2str(A_dash(:)',3),'\s*',' ')]);
                        
                    end
                    
                    
                    
                end % if J==K b
            else % if (sum(min(A,artmap.w(:,j)))./artmap.M >= artmap.rho)
                
                
                
                
                
                if (DESCRIPT==ON)
                    %disp('Did not pass vigilance test. Do not change vigilance. Moving to next category box.')
                    artmap.learn_description=strvcat(artmap.learn_description,'Did not pass vigilance test. Do not change vigilance. Moving to next category box.');
                end
                
                if (artmap.deplete_in_mismatch==1)
                    
                    if LEARN_TRAIL==ON
                        MT_nodes=MT_nodes+1;   
                        artmap.learn_trail(sampleNumber).MT_nodes{MT_nodes}=j;
                    end
                    
                    
                    
                    if (DESCRIPT==ON)
                        %disp('INCORRECT match. Raise vigilance. Go to next closest category box (This is called match tracking, not really vigilance!')
                        artmap.learn_description=strvcat(artmap.learn_description,['INCORRECT match; match tracking and updating rho: ', num2str(artmap.rho)]);
                    end
                    
                    %                     %=========ARTMAP Attention=======
                    
                    %This is without normalization of any sort
                    x_Attention=min(A_dash,artmap.w(:,j));
                    
                    if LEARN_TRAIL==ON
                        artmap.learn_trail(sampleNumber).MT_x_before{MT_nodes} =x_Attention';
                    end
                    
                    
                    
                    
                    
                    
                    
                    %This is with normalization
                    %======================================================
                    %=============
                    if artmap.type_depletion==4
                        
                        if (DESCRIPT==ON)
                            artmap.learn_description=strvcat(artmap.learn_description,['Match Tracking triggers depletion; dyad before normalization: ', regexprep(num2str(x_Attention(:)',3),'\s*',' ')]);
                        end
                        
                        if (artmap.deplete_in_X_denom==1)
                            x_Attention=max(x_Attention+[x_Attention(artmap.M+1:end,1) ; x_Attention(1:artmap.M,1)]-sum(x_Attention)/(sum(A_dash)+eps),0);
                        else
                            x_Attention=max(x_Attention+[x_Attention(artmap.M+1:end,1) ; x_Attention(1:artmap.M,1)]-sum(x_Attention)/artmap.M,0);
                        end
                        
                        if LEARN_TRAIL==ON
                            artmap.learn_trail(sampleNumber).MT_x_after{MT_nodes} =x_Attention';
                        end
                        
                        
                        if (DESCRIPT==ON)
                            
                            artmap.learn_description=strvcat(artmap.learn_description,['Match Tracking triggers depletion; dyad after normalization: ', regexprep(num2str(x_Attention(:)',3),'\s*',' ')]);
                        end
                    end
                    
                    if LEARN_TRAIL==ON
                        artmap.learn_trail(sampleNumber).MT_e_before{MT_nodes} =e_Attention';
                    end
                    
                    
                    
                    e_Attention=max(e_Attention,artmap.lambda_Attention*x_Attention);
                    %Attention can be too harsh if it is depleting along
                    %many dimension. To offset this, we could scale
                    %according to
                    
                    
                    if LEARN_TRAIL==ON
                        artmap.learn_trail(sampleNumber).MT_e_after{MT_nodes} =e_Attention';
                    end
                    
                    if artmap.type_depletion==1
                        %A) the number of dimensions that are non-zero
                        %==========
                        e_Attention=e_Attention/max(sum(sign02(e_Attention)),1);
                        %B) The sum total of attention, implying that attention
                        %is a constant
                    elseif artmap.type_depletion==2
                        e_Attention=e_Attention/(sum(e_Attention)+eps)*artmap.lambda_Attention;
                        %C) The attention graph is used to comparatively scale things implying that attention
                        %is a constant
                    elseif artmap.type_depletion==3
                        e_Attention=e_Attention/(sum(e_Attention)+eps).*e_Attention;
                    end
                    
                    if LEARN_TRAIL==ON
                        artmap.learn_trail(sampleNumber).MT_Adash_before{MT_nodes} =A_dash';
                    end
                    
                    A_dash=rect(A-e_Attention,0);
                    %                     %================================
                    
                    
                    if LEARN_TRAIL==ON
                        artmap.learn_trail(sampleNumber).MT_Adash_after{MT_nodes} =A_dash';
                    end
                    
                    if (DESCRIPT==ON)
                        
                        artmap.learn_description=strvcat(artmap.learn_description,['Match Tracking triggers depletion. e_Attention: ', regexprep(num2str(e_Attention(:)',3),'\s*',' ')]);
                        artmap.learn_description=strvcat(artmap.learn_description,['Depleted input: ', regexprep(num2str(A_dash(:)',3),'\s*',' ')]);
                        
                    end
                end
                
                
                %=========ARTMAP Attention=======
                %Two Options: a) Reset e_Attention
                %              b)continue with depletion
                %================================
                
                %=========ARTMAP Attention=======
                %                  e_Attention=max(e_Attention,artmap.lamda_Attention*min(A_dash,artmap.w(:,j)));
                %                  A_dash=rect(A-e_Attention,0);
                %================================
                
                
            end % if (sum(min(A,artmap.w(:,j)))./artmap.M >= artmap.rho)
        end % for J=T_j
        
        % B.10: ADD A NEW CODING NODE IN F2
        % w(i,C) = A
        if (MATCH==0) % 
            %sampleNumber
            artmap.C=artmap.C+1;
            artmap.w(:,artmap.C)=A;
            
            if LEARN_TRAIL==ON
                artmap.learn_trail(sampleNumber).winning_node=artmap.C;
                artmap.learn_trail(sampleNumber).weight_before = A';
                artmap.learn_trail(sampleNumber).weight_after = A';
            end
            
            
            if (DESCRIPT==ON)
                %sprintf('None of the closest boxes are the right match. So create a new category box %d, with weights:,',artmap.C)  
                artmap.learn_description=strvcat(artmap.learn_description,['None of the closest boxes are the right match; create a new category box ',...
                        num2str(artmap.C), ' with weights: ', regexprep(num2str(artmap.w(:,artmap.C)',3),'\s*',space_amt)]);
                %artmap.w(:,artmap.C)
                %input('We will graph current weights, which should have a new point box, and then we can move to the next input.')
                artmap.learn_description=strvcat(artmap.learn_description,'Finished Learning on this input; moving to next input');
            end
            
            
            if (SHORT_DESCRIPT==ON)
                artmap.short_learn_descript=strvcat(artmap.short_learn_descript,['Adding New Node']);
                artmap.short_learn_descript=strvcat(artmap.short_learn_descript,['Input: ', regexprep(num2str(A',3),'\s*',space_amt)]);
                artmap.short_learn_descript=strvcat(artmap.short_learn_descript,['Weight Node Number: ', num2str(artmap.C)]);
                artmap.short_learn_descript=strvcat(artmap.short_learn_descript,['Newly Added Weight: ', regexprep(num2str(artmap.w(:,artmap.C)',3),'\s*',space_amt)]);
            end
            
            artmap.W=[artmap.W; zeros(1,artmap.numClasses)];
            artmap.W(artmap.C,K)=1;
            if (DISPLAY_WEIGHTBOXES==ON)
                GRAPH_WEIGHTS(artmap,A,K)
            end
            activated_node=artmap.C; % Activated node is new node
        end % (MATCH==0)
        
    end      %  if (artmap.C==0)
    artmap.NodeList = [artmap.NodeList, activated_node];
    
end      % for sampleNumber = 1:numSamples
if (DISPLAY_ONLYFINALWEIGHTBOX==ON)
    GRAPH_WEIGHTS(artmap,A,K)
end

% GRAPH_WEIGHTS(artmap,A,1)


% Fill the new network with the appropriate values.
artmap.search_cycles = search_cycles;
updated_artmap = artmap;

return
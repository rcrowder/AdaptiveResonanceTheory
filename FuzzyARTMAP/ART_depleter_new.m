function [out_e_Attention,out_A_Attention,out_x_Attention] = ART_depleter_new(sampleNumber,weight_j,A,e_Attention,lambda_Attention,M,model_num,varargin)
ON=1;OFF=0;

if nargin==8
    deplete_fast=varargin{1};
else
    deplete_fast=-1;
end


if (((model_num==26) ||  (model_num==28)) && (deplete_fast==0)) %(same as 22, but the mismatch comparison is different)
    %Same as 4, but no deplete in mismatch
    keep_dynamic_depleting=ON;
    e_AttentionTemp=e_Attention;
    dyn_count=1;
    while keep_dynamic_depleting==ON;
        %X_dash=max(weight_j-min(A,weight_j)-lambda_Attention*e_Attention,0);
        X_dash=min(max(A-lambda_Attention*e_AttentionTemp,0),weight_j);
        X_dashDoubled=X_dash+[X_dash(M+1:2*M,1) ; X_dash(1:M,1)];
        e_inc=.01*max(max(X_dashDoubled-sum(X_dash)/M,0)-e_AttentionTemp,0);
        e_AttentionTemp=e_AttentionTemp+e_inc;

        dyn_count=dyn_count+1;
        if (dyn_count==5000)
            disp('reached second threshold');
        end

        if (sum(e_inc)<.0001)
            keep_dynamic_depleting=OFF;
        end

    end
    out_e_Attention=e_AttentionTemp;
    out_A_Attention=max(A-lambda_Attention*e_AttentionTemp,0);
    out_x_Attention=min(out_A_Attention,weight_j);

end



if ((model_num==28) && (deplete_fast==1))  %FAST VERSION: SLOW VERSION COMMENTED OUT BELOW
    %Same as 4, but no deplete in mismatch
    keep_dynamic_depleting=ON;
    e_AttentionTemp=e_Attention;
    dyn_count=1;
    while keep_dynamic_depleting==ON;
        %X_dash=max(weight_j-min(A,weight_j)-lambda_Attention*e_Attention,0);
        X_dash=min(max(A-lambda_Attention*e_AttentionTemp,0),weight_j);
        X_dashDoubled=X_dash+[X_dash(M+1:2*M,1) ; X_dash(1:M,1)];
        e_AttentionPrev=e_AttentionTemp;
        A_prev=max(A-lambda_Attention*e_AttentionTemp,0)';
        e_AttentionTemp=max(max(X_dashDoubled-sum(X_dash)/M,0),e_AttentionTemp);
        A_curr=max(A-lambda_Attention*e_AttentionTemp,0)';
        out_e_Attention=e_AttentionTemp;
        out_A_Attention=max(A-lambda_Attention*e_AttentionTemp,0);
        out_x_Attention=min(max(A-lambda_Attention*e_AttentionTemp,0),weight_j);

        dyn_count=dyn_count+1;
        if (dyn_count==5000)
            disp('reached second threshold');
        end

        %if (sum(abs(e_AttentionPrev-e_AttentionTemp))<10^(-7))
        if (sum(abs(A_prev-A_curr))<10^(-7))
            keep_dynamic_depleting=OFF;
        end

    end

end



if (model_num==29) %(same (or almost the sam) as 28, but the MT formulation is different
    keep_dynamic_depleting=ON;
    e_AttentionTemp=e_Attention;
    dyn_count=1;
    while keep_dynamic_depleting==ON;
        %X_dash=max(weight_j-min(A,weight_j)-lambda_Attention*e_Attention,0);
        A_dep=max(A-lambda_Attention*e_AttentionTemp,0);
        X_dash=min(max(weight_j-lambda_Attention*e_AttentionTemp,0),max(A_dep-lambda_Attention*e_AttentionTemp,0));
        X_dashDoubled=X_dash+[X_dash(M+1:2*M,1) ; X_dash(1:M,1)];
        e_inc=.01*max(max(X_dashDoubled-sum(X_dash)/M,0)-e_AttentionTemp,0);
        e_AttentionTemp=e_AttentionTemp+e_inc;

        dyn_count=dyn_count+1;
        if (dyn_count==5000)
            disp('reached second threshold');
        end

        if (sum(e_inc)<.0001)
            keep_dynamic_depleting=OFF;
        end

    end
    out_e_Attention=e_AttentionTemp;
    out_A_Attention=max(A-lambda_Attention*e_AttentionTemp,0);
    out_x_Attention=min(out_A_Attention,weight_j);

end


if ((model_num==30) && (deplete_fast==0)) %SLOW VERSION (same (or almost the sam) as 28, but the MT formulation is different
    keep_dynamic_depleting=ON;
    e_AttentionTemp=e_Attention;
    dyn_count=1;
    while keep_dynamic_depleting==ON;
        %X_dash=max(weight_j-min(A,weight_j)-lambda_Attention*e_Attention,0
        %);
        X_dash=min(max(weight_j-lambda_Attention*e_AttentionTemp,0),max(A-lambda_Attention*e_AttentionTemp,0));
        X_dashDoubled=X_dash+[X_dash(M+1:2*M,1) ; X_dash(1:M,1)];
        e_inc=.01*max(max(X_dashDoubled-sum(X_dash)/M,0)-e_AttentionTemp,0);
        e_AttentionTemp=e_AttentionTemp+e_inc;

        dyn_count=dyn_count+1;
        if (dyn_count==5000)
            disp('reached second threshold');
        end

        if (sum(e_inc)<.0001)
            keep_dynamic_depleting=OFF;
        end

    end
    out_e_Attention=e_AttentionTemp;
    out_A_Attention=max(A-lambda_Attention*e_AttentionTemp,0);
    out_x_Attention=min(out_A_Attention,max(weight_j-lambda_Attention*e_AttentionTemp,0));
end


if ((model_num==30) && (deplete_fast==1))  %FAST VERSION:

    e_AttentionTemp=e_Attention;


    %X_dash=max(weight_j-min(A,weight_j)-lambda_Attention*e_Attention,0);
    X_dash=min(max(weight_j-lambda_Attention*e_AttentionTemp,0),max(A-lambda_Attention*e_AttentionTemp,0));
    X_dashDoubled=X_dash+[X_dash(M+1:2*M,1) ; X_dash(1:M,1)];
    e_AttentionTemp=max(max(X_dashDoubled-sum(X_dash)/M,0),e_AttentionTemp);
    out_e_Attention=e_AttentionTemp;
    out_A_Attention=max(A-lambda_Attention*e_AttentionTemp,0);
    out_x_Attention=min(max(A-lambda_Attention*e_AttentionTemp,0),max(weight_j-lambda_Attention*e_AttentionTemp,0));



end


if ((model_num==31) && (deplete_fast==0)) %SLOW VERSION (same (or almost the sam) as 28, but the MT formulation is different
    keep_dynamic_depleting=ON;
    e_AttentionTemp=e_Attention;
    dyn_count=1;
    if (lambda_Attention>0)
    while keep_dynamic_depleting==ON;
        %X_dash=max(weight_j-min(A,weight_j)-lambda_Attention*e_Attention,0
        %);
        X_dash=min(max(weight_j-lambda_Attention*e_AttentionTemp,0),max(A-lambda_Attention*e_AttentionTemp,0));
        X_dashDoubled=X_dash+[X_dash(M+1:2*M,1) ; X_dash(1:M,1)];
        e_inc=.0001*max(max(X_dashDoubled-sum(min(A,weight_j))/M,0)-e_AttentionTemp,0);
        e_AttentionTemp=e_AttentionTemp+e_inc;

        dyn_count=dyn_count+1;
        if (dyn_count==50000)
            disp('reached second threshold');
        end

        if (sum(e_inc)<.000001)
            keep_dynamic_depleting=OFF;
        end

    end
    out_e_Attention=e_AttentionTemp;
    out_A_Attention=max(A-lambda_Attention*e_AttentionTemp,0);
    out_x_Attention=min(out_A_Attention,max(weight_j-lambda_Attention*e_AttentionTemp,0));
    else
    out_e_Attention=e_AttentionTemp;
    out_A_Attention=A;
    out_x_Attention=min(A,weight_j);
        
    end
end




if ((model_num==31) && (deplete_fast==1))  %FAST VERSION:

    e_AttentionTemp=e_Attention;


    %X_dash=max(weight_j-min(A,weight_j)-lambda_Attention*e_Attention,0);
    X_dash=min(max(weight_j-e_AttentionTemp,0),max(A-e_AttentionTemp,0));
    X=min(weight_j,A);
    X_dashDoubled=X_dash+[X_dash(M+1:2*M,1) ; X_dash(1:M,1)];
    e_increasers=sign02(X_dashDoubled-sum(min(A,weight_j))/M);
    e_C_case=sign02(sum(X)/M-abs(X(1:M)-X(M+1:2*M)));
    e_C=.5*(X(1:M)+X(M+1:2*M)-sum(X)/M);
    e_A=(X(1:M)-sum(X)/M).*sign02(max(e_AttentionTemp(1:M),e_C)-X(M+1:2*M)).*(~e_C_case);
    e_B=(X(M+1:2*M)-sum(X)/M).*sign02(max(e_AttentionTemp(1:M),e_C)-X(1:M)).*(~e_C_case);
    e_C=e_C.*e_C_case;
    %        X_Doubled=X+[X(M+1:2*M,1) ; X(1:M,1)];
    e_new_possible=[e_A+e_B+e_C;e_A+e_B+e_C];
    if (sum((e_new_possible<e_AttentionTemp).*(e_new_possible>0).*(e_increasers))>0)
        error('Weirdness!');
    end
    out_e_Attention=e_AttentionTemp.*(~e_increasers)+e_increasers.*e_new_possible;
    out_A_Attention=max(A-out_e_Attention,0);
    out_x_Attention=min(max(A-out_e_Attention,0),max(weight_j-out_e_Attention,0));



end



if ((model_num==32) && (deplete_fast==0))  %SLOW VERSION NOT ACTUALLY SLOW:

    if lambda_Attention>0
        e_AttentionTemp=e_Attention;
        f_AttentionTemp=e_Attention.*lambda_Attention;

        %X_dash=max(weight_j-min(A,weight_j)-lambda_Attention*e_Attention,0);
        X_dash=min(max(weight_j-f_AttentionTemp,0),max(A-f_AttentionTemp,0));
        X=min(weight_j,A);
        e_increasers=sign02(X_dash-sum(min(A,weight_j))/(2*M)-e_AttentionTemp);
        f_new_possible=(X-sum(X)/(2*M))/(1+(1/lambda_Attention));
        e_new_possible=f_new_possible/lambda_Attention;
        f_new=f_AttentionTemp.*(~e_increasers)+e_increasers.*f_new_possible;
        if (sum((e_new_possible<e_AttentionTemp).*(e_new_possible>0).*(e_increasers))>0)
            error('Weirdness!');
        end
        out_e_Attention=f_new/lambda_Attention;
        out_A_Attention=max(A-lambda_Attention.*out_e_Attention,0);
        out_x_Attention=min(max(A-lambda_Attention.*out_e_Attention,0),max(weight_j-lambda_Attention.*out_e_Attention,0));
    else
        out_e_Attention=e_Attention;
        out_A_Attention=A;
        out_x_Attention=min(weight_j,A);
    end



end



if ((model_num==32) && (deplete_fast==1))  %FAST VERSION:

     e_AttentionTemp=e_Attention;
     f_AttentionTemp=e_Attention;

        %X_dash=max(weight_j-min(A,weight_j)-lambda_Attention*e_Attention,0);
        X_dash=min(max(weight_j-f_AttentionTemp,0),max(A-f_AttentionTemp,0));
        X=min(weight_j,A);
        e_increasers=sign02(X_dash-sum(min(A,weight_j))/(2*M)-e_AttentionTemp);
        f_new_possible=(X-sum(X)/(2*M));
        e_new_possible=f_new_possible;
        f_new=f_AttentionTemp.*(~e_increasers)+e_increasers.*f_new_possible;
        if (sum((e_new_possible<e_AttentionTemp).*(e_new_possible>0).*(e_increasers))>0)
            error('Weirdness!');
        end
        out_e_Attention=f_new;
        out_A_Attention=max(A-out_e_Attention,0);
        out_x_Attention=min(max(A-out_e_Attention,0),max(weight_j-out_e_Attention,0));


end


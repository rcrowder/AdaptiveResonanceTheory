
lambdaGrid=[(.1:.1:1) (1.25:.25:10) (15:15:100)];
store_acc=zeros(52,1);

for i=1:52

    close all;numEpochs=1;lambda_Attention=lambdaGrid(i);test_fast=0;test_only=0;biasedARTMAPTester
    title_genres=['Genre Prediction  \lambda=' num2str(c.lambda_Attention,2)];
    [conf_nonnorm,conf_norm,avg_acc]= confMatrixGail(dataStructTemp.test_output(:),classes,c,dataStructTemp.labels,title_genres);
    disp(avg_acc)
    store_acc(i)=avg_acc;
    save store_acc_biased_ARTMAP_movieGenres store_acc
end
    
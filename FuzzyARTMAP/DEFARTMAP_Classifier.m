%This code performs the same as DEFARTMAP_Classify_H, but probably runs
%faster. I say probably because hope should never wilt.
M_test=trainingNetRevised.M;
C_test=trainingNetRevised.C;
w_test=trainingNetRevised.w;
W_test=trainingNetRevised.W;
alpha_test=trainingNetRevised.alpha;
test_inp_vect=(sampleNumber)


T_chosen= sum(min(repmat(test_inp_vect,1,C_test),w_test),1)+(1-alpha_test).*(M-sum(w_test,1));
Lambda=find(T_chosen>alpha_test*M_test);
Lambda_prime=find(T_chosen==M_test);

y=zeros(1,C_test);
if ~isempty(Lambda_prime)
    y(Lambda_prime)=1./length(Lambda_prime);
else
    sumL = sum(1./(M_test-T_chosen(Lambda)))';
    y(Lambda)=(1./(artmap.M-T(Lambda)))./sumL;
    
end

%yMatrix(sampleNumber,:)= y;
%[dummy, class]=max(sigma);
%classes(sampleNumber)=class;


sigma=W_test'*y';
sigmaS(:,sampleNumber)=sigma
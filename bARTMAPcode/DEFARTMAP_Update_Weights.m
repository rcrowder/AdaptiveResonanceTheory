function updated_weights = DEFART_Update_Weights(artmap, J, A)
% B.11 LEARNING
artmap.w(:,J) = artmap.beta.*(min(A,artmap.w(:,J))) + (1-artmap.beta).*artmap.w(:,J);
updated_weights = artmap.w;

return
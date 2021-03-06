function out_labels=addClassError(input_labels,perc_error)
num_points=length(input_labels);
class_vals=1:1:max(input_labels);
[dum,revs_map]=sort(class_vals);
num_errors=round(perc_error*num_points/100);
dum=randperm(num_points);
change_inds=dum(1:num_errors);
rand_vals=ceil(rand(1,num_errors)*(length(class_vals)-1));
new_labels=rand_vals+revs_map(input_labels(change_inds));
new_labels(new_labels>length(class_vals))=new_labels(new_labels>length(class_vals))-length(class_vals);
new_labels=class_vals(new_labels);
out_labels=input_labels;
out_labels(change_inds)=new_labels;
out_labels=reshape(out_labels,size(input_labels));

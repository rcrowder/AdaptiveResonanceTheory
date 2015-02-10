
dataStruct(1).training_input=dep_test(:,1:2)';
dataStruct(1).training_output=dep_test(:,3)'+1;
dataStruct(1).test_input=dep_test(:,1:2)';
dataStruct(1).test_output=dep_test(:,3)'+1;
dataStruct(1).description='six_point';
dataStruct(1).descriptionVerbose='Canonical Depletion Test';

check_rand=randperm(2000);
check_train_in=checkerboard_train_denser(check_rand,1:2)';
check_train_out=checkerboard_train_denser(check_rand,3)'+1;
check_train_in_sparse=checkerboard_train_denser(check_rand(1:500),1:2)';
check_train_out_sparse=checkerboard_train_denser(check_rand(1:500),3)'+1;


load ARTMAP_data dataStruct
load interpatch_leave_3
inter_rand=randperm(1000);
inter_train_in=interpatch_train_denser(inter_rand,1:2)';
inter_train_out=interpatch_train_denser(inter_rand,3)'+1;
inter_train_in_sparse=interpatch_train_denser(inter_rand(1:200),1:2)';
inter_train_out_sparse=interpatch_train_denser(inter_rand(1:200),3)'+1;



dataStruct(2).training_input=inter_train_in_sparse;
dataStruct(2).training_output=inter_train_out_sparse;
dataStruct(2).test_input=interpatch_test_denser(:,1:2)';
dataStruct(2).test_output=interpatch_test_denser(:,3)'+1;
dataStruct(2).description='stripes_sparse';
dataStruct(2).descriptionVerbose='stripes_sparse 6X1 200 training points';

dataStruct(3).training_input=inter_train_in;
dataStruct(3).training_output=inter_train_out;
dataStruct(3).test_input=interpatch_test_denser(:,1:2)';
dataStruct(3).test_output=interpatch_test_denser(:,3)'+1;
dataStruct(3).description='stripes_dense';
dataStruct(3).descriptionVerbose='stripes_dense 6X1 1000 training points';



dataStruct(4).training_input=cis_test(1:100,1:2)';
dataStruct(4).training_output=cis_train2(1:100,3)'+1;
dataStruct(4).test_input=cis_train2(:,1:2)';
dataStruct(4).test_output=cis_test(:,3)'+1;
dataStruct(4).description='cis_sparse';
dataStruct(4).descriptionVerbose='CIS 100 training points';

dataStruct(5).training_input=cis_test(:,1:2)';
dataStruct(5).training_output=cis_train2(:,3)'+1;
dataStruct(5).test_input=cis_train2(:,1:2)';
dataStruct(5).test_output=cis_test(:,3)'+1;
dataStruct(5).description='cis_dense';
dataStruct(5).descriptionVerbose='CIS 1000 training points';

dataStruct(6).training_input=check_train_in_sparse;
dataStruct(6).training_output=check_train_out_sparse;
dataStruct(6).test_input=checkerboard_test_denser(:,1:2)';
dataStruct(6).test_output=checkerboard_test_denser(:,3)'+1;
dataStruct(6).description='checkboard_sparse';
dataStruct(6).descriptionVerbose='checkboard_sparse 6X6 500 training points';

dataStruct(7).training_input=check_train_in;
dataStruct(7).training_output=check_train_out;
dataStruct(7).test_input=checkerboard_test_denser(:,1:2)';
dataStruct(7).test_output=checkerboard_test_denser(:,3)'+1;
dataStruct(7).description='checkboard_dense';
dataStruct(7).descriptionVerbose='checkboard_dense 6X6 2000 training points';


dataStruct(8).training_input=bin_six_check(:,1:5)';
dataStruct(8).training_output=bin_six_check(:,6)'+1;
dataStruct(8).test_input=bin_five_test';
dataStruct(8).test_output=[];
dataStruct(8).description='binary_test';
dataStruct(8).descriptionVerbose='Synthetic 5-D Binary Dataset';


frac_red=.55/3;
frac_blue=(1-frac_red*3)/3;
frac_net=frac_red+frac_blue;

training_pos=rand(2,1000);
test_pos=rand(2,5000);
training_label=ceil(mod(training_pos(2,:),frac_net)/frac_red);
test_label=ceil(mod(test_pos(2,:),frac_net)/frac_red);


dataStruct(9).training_input=training_pos(:,1:200);
dataStruct(9).training_output=training_label(1:200);
dataStruct(9).test_input=test_pos;
dataStruct(9).test_output=test_label;
dataStruct(9).description='stripes_sparse_unequal';
dataStruct(9).descriptionVerbose='stripes_sparse_unequal 6X1 200 training points';

dataStruct(10).training_input=training_pos;
dataStruct(10).training_output=training_label;
dataStruct(10).test_input=test_pos;
dataStruct(10).test_output=test_label;
dataStruct(10).description='stripes_dense_unequal';
dataStruct(10).descriptionVerbose='stripes_dense_unequal 6X1 1000 training points';

for i=1:13
    min_val=min(dataStruct(11).training_input(i,:));
    min_val=min(min_val,min(dataStruct(11).test_input(i,:)));
    max_val=max(dataStruct(11).training_input(i,:));
    max_val=max(max_val,max(dataStruct(11).test_input(i,:)));
    dataStruct(11).training_input(i,:)=(dataStruct(11).training_input(i,:)-min_val+.05)/max_val;
    dataStruct(11).test_input(i,:)=(dataStruct(11).test_input(i,:)-min_val+.05)/max_val;
    
end


save ARTMAP_data dataStruct
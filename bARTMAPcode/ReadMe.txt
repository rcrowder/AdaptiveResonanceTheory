Author:

Sai Chaitanya Gaddam
v 1.0 04/10/2009




Description:

This folder contains the code and data required to run all the examples in the biased ART paper [2]


The following examples can be run by executing biasedARTMAP_examples:

    '1) Six Point Dataset'
    '2) Stripes benchmark (sparse)'
    '3) Stripes benchmark (dense)'
    '4) Circle in Square benchmark (sparse)'
    '5) Circle in Square benchmark (dense)'
    '6) Checkerboard benchmark (sparse)'
    '7) Checkerboard benchmark (dense)'
    '8) 6-D binary dataset'
    '9) Boston Benchmark: test on stripe 1'
    '10) Boston Benchmark: test on stripe 2'
    '11) Boston Benchmark: test on stripe 3'
    '12) Boston Benchmark: test on stripe 4'
    '13) Movie Genre Benchmark'


To provide your own dataset, run biasedARTMAPTester

Usage: [a,b,c] = biasedARTMAPTester(dataStruct,lambda_value)

lambda_value = 0 -> biased ARTMAP is equivalent to fuzzy ARTMAP
lambda_value = 10 -> found to be optimal for many benchmarks 




The MATLAB struct dataStruct should have the following format.

The datastruct fields are:


        training_input: [f features x m records]
       training_output: [m labels x1]
            test_input: [f features x n records]
           test_output: [n labels x1]
           description: 'dataset_title'
    descriptionVerbose: 'A more verbose description of the dataset'



References:

[1]Bennet & Lanning. (2007). The Netflix Prize
	http://www.cs.uic.edu/~liub/KDD-cup-2007/NetflixPrize-description.pdf

[2]Carpenter & Gaddam, (2009). Biased ART: A neural architecture that shifts attention toward previously disregarded features following an incorrect prediction. Technical Report CAS/CNS TR-2009-003, Boston, MA: Boston University.
        http://cns.bu.edu/~gsc/biasedART/TR-2009-003_CarpenterGaddam_.pdf


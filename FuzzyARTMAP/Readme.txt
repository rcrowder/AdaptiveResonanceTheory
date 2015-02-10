Description:

This folder contains the code and data required to run all the examples in the description below.

To run GUI, run ARTMAPgui

The following examples can be run by executing the gui ARTMAPgui, or from the command line using fuzzyARTMAP_examplesWrapper:

    '1) Circle in Square benchmark (sparse)'
    '2) Circle in Square benchmark (dense)'
    '3) Stripes benchmark (sparse)'
    '4) Stripes benchmark (dense)'
    '5) Checkerboard benchmark (sparse)'
    '6) Checkerboard benchmark (dense)'
    '7) Boston Benchmark: test on stripe 1'
    '8) Boston Benchmark: test on stripe 2'
    '9) Boston Benchmark: test on stripe 3'
    '10) Boston Benchmark: test on stripe 4'
    '11) Movie Genre Benchmark'


To provide your own dataset, run fuzzyARTMAPTester

Usage: [a,b,c] = fuzzyARTMAPTester(dataStruct)




The MATLAB struct dataStruct should have the following format.

The datastruct fields are:


        training_input: [f features x m records]
       training_output: [m labels x1]
            test_input: [f features x n records]
           test_output: [n labels x1]
           description: 'dataset_title'
    descriptionVerbose: 'A more verbose description of the dataset'



References:



[1](Carpenter, G. A., Grossberg, S., Markuzon, N., Reynolds, J.H., & Rosen, D.B. (1992) 
Fuzzy ARTMAP: A neural network architecture for incremental supervised learning of analog multidimensional maps. IEEE Transactions on Neural Networks, 3, 698-713.) 
        


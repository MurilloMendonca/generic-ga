# Generic GA

This project is a genetic implmentation of a Genetic Algorithm (GA) in C++. The classes use template parameters and virtual function to allow for a wide range of applications. The user needs to define: the type of data a Individual will hold (by template parameter); The features encoding e decoding to a string of '0's and '1's strategy (by a static objetc to the Individual class); The fitness scoring function (by a virtual function of the GA class); The generation of the initial population (by a virtual function of the GA class); And a function to verify if a Individual is valid (by a virtual function of the GA class).

A simple example is provided in the test.cpp file.

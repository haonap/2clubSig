# Code for computing 2-club signatures
This repository contains C++ code used for computing 2-club signatures in the article "Graph Signatures: Identification and Optimization" which has been submitted to European Journal of Operational Research. If you wish to use or cite this code, please cite:
      
    @article{BBJBHP2021g-sign,
        author = {Balabhaskar Balasundaram and Juan S. Borrero and Hao Pan},
        journal = {European Journal of Operational Research},
        note = {Under Review},
        title = {Graph Signatures: {I}dentification and Optimization},
        year = {2021}}

# Understanding and using the code
The code should be straightforward if you start to read from file main.cpp. Necessary comments have been added in the code for easiness of understanding. As stated previously, the code is used for computing 2-club signatures. We present two methods for computing 2-club signatures in the code, GSIP-F2(named as GSIP_F2 in the code) and MW. 

Let's start from file main.cpp. Line 7-13 is reading in three input parameters(instance name, tau, and method) from file instance.txt. For instance name, you can go into folder graphSequences in this repository and pick any one of them. Each instance is actually a graph sequence. For tau, you can pick any positive integer less than or equal to the length of the graph sequence. For method, you can either pick 1 or 2, 1 for GSIP-F2 and 2 for MW. With these settings, the code is going to compute a tau-persistent 2-club signature of the graph sequence. For example, the default instance.txt file contains "karate_10_0.8 3 2" which is telling the code to compute a 3-persistent 2-club signature of the graph sequence karate_10_0.8 using MW method. Line 15 is reading in the graph sequence specified in file instance.txt. Line 18-22 triggers computation. Functions ReadIn, GSIP_F2, MW are defined in files functions.cpp and functions.h. There are other functions defined in files functions.cpp and functions.h, and they would be called by those three foregoing functions during execution. 



# Compilation and execution in Linux environment
1. Download or clone the repository to your machine. 
2. From terminal, go to the repository. 
3. Type "make" and hit enter to compile. 
4. Open instance.txt file, enter instance name, tau, and method sequentially as input parameters. 
5. In terminal, type "./main" and hit enter to execute. 




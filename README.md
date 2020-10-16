# Code for computing 2-club signatures
This repository contains C++ code used for computing 2-club signatures in the article "Graph Signatures: Identification and Optimization" which has been submitted to European Journal of Operational Research. If you wish to use or cite this code, please cite:
        
    @misc{MW2CLB-GithubRepo,
      Author = {Hao Pan and Balabhaskar Balasundaram and Juan S. Borrero},
      Date-Added = {2020-10-15 17:38:24 -0500},
      Date-Modified = {2020-10-15 17:38:59 -0500},
      Howpublished = {Codes and instances online at: \url{https://github.com/haonap/2clubSig}},
      Month = {October },
      Title = {Implementation of the moving window method for the maximum 2-club signature problem.},
      Year = {2020}}
      
# Compilaton and execution in Linux environment
1. Download or clone the repository. 
2. From terminal, go to the repository. 
3. Type "make" and hit enter to compile. 
4. Open instance.txt file, which contains input information. There are 3 entries in this file: instance name(directory copy one from the "graphSequences" folder), tau, method(1 for GSIP-F2, 2 for MW). With these settings, it is going to compute a tau-persistent 2-club signature of the instance. 
5. In terminal, type "./main" and hit enter to execute. 
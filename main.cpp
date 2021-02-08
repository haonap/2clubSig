#include <iostream>
#include "functions.h"

using namespace std;

int main(){
    string seqName;
    int tau, method;
    fstream fin("instance.txt", fstream::in); // read in input parameters
    fin >> seqName;
    fin >> tau;
    fin >> method;
    fin.close();

    ReadIn(seqName); // read in graph sequence

    if(method == 1){ // pick a method, 1 for GSIP-F2, 2 for MW-2CLB, 3 for MW-F2
        GSIP_F2(tau);
    }else{
        MW(tau, method);
    }


    return 0;
}

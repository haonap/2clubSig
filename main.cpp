#include <iostream>
#include "functions.h"

using namespace std;

int main(){
    string seqName;
    int tau, method;
    fstream fin("instance.txt", fstream::in);
    fin >> seqName;
    fin >> tau;
    fin >> method;
    fin.close();

    ReadIn(seqName);

    if(method == 1){
        GSIP_F2(tau);
    }else if(method == 2){
        MW(tau);
    }

    return 0;
}

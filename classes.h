//
// Created by Hao Pan on 9/6/20.
//

#ifndef INC_2CLUBSIG_CLASSES_H
#define INC_2CLUBSIG_CLASSES_H


#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include "gurobi_c++.h"
using namespace std;

class node{
public:
    int name;
    int degree = 0;
    vector<int> neighbors;
    vector<int> kNeighbors;
    vector<int> adjacency; // used in CommunityPeel
    vector<int> F; // used in CommunityPeel
    vector<int> numCommonNeighbors; // used in CommunityPeel
    vector<vector<int>> commonNeighbors; // used in CommunityPeel

    //constructor
    node(int);
};

class graph {
public:

    string name;
    int n, m;
    vector<node> nodeList;
    int maxDeg;
    int maxDegNode;
    vector<int> kCore;

    // constructor
    graph(string);
    graph(){};

    // functions
    void GetKCore(int); // this function gets kCore of a single static graph without reordering vertices
    bool IsAdj(int, int); // this function checks if two nodes are adjacent
    bool IsKAdj(int, int);
    vector<int> FindCommonNeighbors(int, int);
    void FindKNeighbors(int k);
    vector<int> BFS(int);
    vector<int> KBFS(int, int);
};

class commonNeighborLazyCut: public GRBCallback{
public:
    GRBVar* xvar;
    vector<int>* nodesP;
    vector<int>* nodeMapP;
    vector<graph>* graphSeqP;
    int windowHead;
    int tau;

    commonNeighborLazyCut(GRBVar*, vector<int>*, vector<int>*, vector<graph>*, int, int);

protected:
    //callback function
    void callback();
};


#endif //INC_2CLUBSIG_CLASSES_H

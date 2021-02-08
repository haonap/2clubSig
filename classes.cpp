#include <iostream>
#include "classes.h"
using namespace std;


node::node(int nameN): name(nameN){}


graph::graph(string nameG){

    name = nameG;
    ifstream fin(name);
    string buf;
    int u, v;
    fin >> buf;
    fin >> n;
    fin >> buf;
    fin >> m;

    for(int i = 0; i < n; i++){
        node pNode(i);
        nodeList.push_back(pNode);
    }

    while(fin >> buf){
        if(buf == "e"){
            fin >> u;
            fin >> v;
            nodeList[u - 1].neighbors.push_back(v - 1);
            nodeList[v - 1].neighbors.push_back(u - 1);
        }
    }

    maxDeg = 0;
    for(int i = 0; i < n; i++){
        nodeList[i].degree = (int)nodeList[i].neighbors.size();
        if(nodeList[i].neighbors.size() > maxDeg){
            maxDeg = (int)nodeList[i].neighbors.size();
            maxDegNode = i;
        }
    }
}


// get k-core of a graph
void graph::GetKCore(int k){
    queue<int> Q;
    vector<int> F, D;
    kCore.clear();

    for(int i = 0; i < n; i++){
        D.push_back((int)nodeList[i].neighbors.size());
        if(D[i] < k){
            F.push_back(1);
            Q.push(i);
        }else{
            F.push_back(0);
        }
    }

    while (Q.size()) {
        int v = Q.front();
        Q.pop();
        if(D[v] > 0){
            D[v] = 0;
            for(int u = 0; u < nodeList[v].neighbors.size(); u++){
                int pNeighbor = nodeList[v].neighbors[u];
                if(D[pNeighbor] > 0){
                    D[pNeighbor]--;
                    if(D[pNeighbor] < k){
                        if(F[pNeighbor] == 0){
                            Q.push(pNeighbor);
                            F[pNeighbor] = 1;
                        }
                    }
                }
            }
        }
    }

    for(int i = 0; i < n; i++){
        if(D[i] > 0){
            kCore.push_back(i);
        }
    }
}


// check if u is adjacent to v
bool graph::IsAdj(int u, int v){
    if(binary_search(nodeList[u].neighbors.begin(), nodeList[u].neighbors.end(), v)){
        return 1;
    }else{
        return 0;
    }
}


// check if u is in distance at most k to v
bool graph::IsKAdj(int u, int v){
    if(binary_search(nodeList[u].kNeighbors.begin(), nodeList[u].kNeighbors.end(), v)){
        return 1;
    }else{
        return 0;
    }
}


// find common neighbors of u and v
vector<int> graph::FindCommonNeighbors(int u, int v){
    vector<int> comNeighbors;
    int i = 0, j = 0;
    while(i < nodeList[u].neighbors.size() && j < nodeList[v].neighbors.size()){
        if(nodeList[u].neighbors[i] == nodeList[v].neighbors[j]){
            comNeighbors.push_back(nodeList[u].neighbors[i]);
            i++;
            j++;
        }else if (nodeList[u].neighbors[i] > nodeList[v].neighbors[j]){
            j++;
        }else{
            i++;
        }
    }
    return comNeighbors;
}


// find distance-k neighbors of all vertices
void graph::FindKNeighbors(int k) {
    for(int i = 0; i < n; i++){
        nodeList[i].kNeighbors.clear();
        nodeList[i].kNeighbors = KBFS(i, k);
    }
}


// find all vertices connected to a vertex
vector<int> graph::BFS(int s){
    vector<int> tempVec;
    vector<int> distance(n, n);
    queue<int> Q;
    distance[s] = 0;
    Q.push(s);
    int u, v;
    int goFlag = 1;

    while(Q.size() && goFlag){
        u = Q.front();
        Q.pop();
        for(int i = 0; i < nodeList[u].neighbors.size(); i++){
            v = nodeList[u].neighbors[i];
            if(distance[v] > distance[u] + 1){
                distance[v] = distance[u] + 1;
                Q.push(v);
            }
        }
    }

    for(int i = 0; i < distance.size(); i++){
        if(distance[i] < n && i != s){
            tempVec.push_back(i);
        }
    }
    return tempVec;
}


// find distance-k neighbors of a vertex
vector<int> graph::KBFS(int s, int k){
    vector<int> tempVec;
    vector<int> distance(n, n);
    queue<int> Q;
    distance[s] = 0;
    Q.push(s);
    int u, v;
    int goFlag = 1;

    while(Q.size() && goFlag){
        u = Q.front();
        Q.pop();
        for(int i = 0; i < nodeList[u].neighbors.size(); i++){
            v = nodeList[u].neighbors[i];
            if(distance[v] > distance[u] + 1){
                distance[v] = distance[u] + 1;
                if(distance[v] > k){ // KBFS will stop at level k
                    goFlag = 0;
                    break;
                }else{
                    Q.push(v);
                }

            }
        }
    }

    for(int i = 0; i < n; i++){
        if(distance[i] <= k && i != s){
            tempVec.push_back(i);
        }
    }
    return tempVec;
}


// Gurobi callback, adding common neighbor constraints as lazy cuts
commonNeighborLazyCut::commonNeighborLazyCut(GRBVar* xvar, vector<int>* nodesP, vector<int>* nodeMapP, vector<graph>* graphSeqP, int windowHead, int tau): xvar(xvar), nodesP(nodesP), nodeMapP(nodeMapP), graphSeqP(graphSeqP), windowHead(windowHead), tau(tau){}
void commonNeighborLazyCut::callback(){
    try{
        if(where==GRB_CB_MIPSOL){ //if get an integer solution
            double* xval = getSolution(xvar, (int)nodesP->size());
            vector<int> curSol;

            for(int i = 0; i < (int)nodesP->size(); i++){
                if(xval[i] > 0.5){
                    curSol.push_back((*nodesP)[i]);
                }
            }

            vector<int> commonNeighbors;
            int count;
            for(int u = 0; u < curSol.size(); u++){
                int flag3 = 0;
                int uNode = curSol[u];
                for(int v = u + 1; v < curSol.size(); v++){
                    int flag2 = 0;
                    int vNode = curSol[v];
                    for(int t = windowHead; t < windowHead + tau; t++){
                        int flag1 = 0;
                        if((*graphSeqP)[t].IsAdj(uNode, vNode) == 0){
                            commonNeighbors = (*graphSeqP)[t].FindCommonNeighbors(uNode, vNode);
                            count = 0;
                            for(int w = 0; w < commonNeighbors.size(); w++){
                                int pNeighbor = commonNeighbors[w];

                                if((*nodeMapP)[pNeighbor] > -1){
                                    if(xval[(*nodeMapP)[pNeighbor]] > 0.5){
                                        break;
                                    }else{
                                        count++;
                                    }
                                }else{
                                    count++;
                                }
                            }
                            if(count == (int)commonNeighbors.size()){
                                flag1 = 1;
                                flag2 = 1;
                                flag3 = 1;

                                int mCount = 0;
                                GRBLinExpr expr = 0;
                                for(int p = 0; p < commonNeighbors.size(); p++){
                                    int pNeighbor = commonNeighbors[p];
                                    if((*nodeMapP)[pNeighbor] > -1){
                                        expr += xvar[(*nodeMapP)[pNeighbor]];
                                    }
                                }
                                addLazy(xvar[(*nodeMapP)[uNode]] + xvar[(*nodeMapP)[vNode]] - expr <= 1);
                            }
                        }
                        if(flag1 == 1){
                            break;
                        }
                    }
                    if(flag2 == 1){
                        break;
                    }
                }
                if(flag3 == 1){
                    break;
                }
            }
        }
    }catch(GRBException e){
        cout << "Error number during: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }catch(...){
        cout << "Error during callback" << endl;
    }
}

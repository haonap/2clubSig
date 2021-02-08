//
// Created by Hao Pan on 9/6/20.
//

#ifndef INC_2CLUBSIG_FUNCTIONS_H
#define INC_2CLUBSIG_FUNCTIONS_H


#include <iostream>
#include "classes.h"
#include <fstream>
#include <vector>
#include <sys/time.h>
#include <iomanip>
#include <queue>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <cstdlib>
#include <cstring>
#include "gurobi_c++.h"
#include <queue>

/* for getDir() */
#include <stdio.h>  // defines FILENAME_MAX
// #define WINDOWS  // uncomment this line to use it for windows
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

using namespace std;

void ReadIn(string);
void GSIP_F2(int);
void MW(int, int);
vector<int> GetPersistent2ClubWindow(int, int, int*, vector<int>*);
vector<int> GetPersistent2ClubWindowGeneric(int, int, vector<int>*);
graph GetIntersectionGraphOfPowerGraphs(int, int);

graph GetIntersectionGraph(int, int);
void CorePeel(graph*, int*);
void CommunityPeel(graph*, int*);
vector<vector<int>> GetComponents(graph*, int*);
bool cmp(const vector<int> &,const vector<int> &);

string itos_c(int);
int stoi_c(string);
double stod_c(string);
string dtos_c(double);

void makeDir(string);
void emptyDir(string);
void goToDir(string);
string getDir();

double get_wall_time();
double get_cpu_time();


#endif //INC_2CLUBSIG_FUNCTIONS_H

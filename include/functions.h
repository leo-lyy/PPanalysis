#pragma once
#include<vector>
#include<cmath>
#include<iostream>
#include"system.h"

using namespace std;

void unwrap(System& system);
void wrap(System& system);
double norm(const axis& rr);
void cross(axis& zz, axis& xx, axis yy);
double dot(axis& A, axis& B);
axis unitVec(axis& vec);
void vecInit(vector<vector<axis>>& v, long int molNum, const long int maxhelix);
double degVec(axis& a, axis& b);
void printProgressBar(int k, int kmax);
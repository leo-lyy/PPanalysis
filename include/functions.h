#pragma once
#include<vector>
#include<cmath>
#include"system.h"

using namespace std;

void unwrap(System& system);
void wrap(System& system);
double norm(const axis& rr);
void cross(axis& zz, axis& xx, axis yy);
double dot(axis& A, axis& B);
axis unitVec(axis& vec);
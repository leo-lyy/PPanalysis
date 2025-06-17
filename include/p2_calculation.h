#pragma once
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include<chrono>
#include<cmath>
#include<algorithm>
#include<stdexcept>
#include"system.h"
#include"functions.h"

using namespace std;

void printProgressP2(int k, int kmax, double p2);
double p2calculation(System& system, long int f);
void dumpIO_p2(System& system, ifstream& dumpFilein);

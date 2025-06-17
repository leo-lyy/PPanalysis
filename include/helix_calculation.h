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

void helixCalculation(System& system);
void dumpIO_helix(System& system, ifstream& dumpFilein, ofstream& dumpFileout);
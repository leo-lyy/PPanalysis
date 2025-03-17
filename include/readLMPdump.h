#pragma once
#include <iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include<chrono>
#include<cmath>
#include"system.h"
#include"functions.h"
#include"p2calculation.h"
using namespace std;

long int countDumpFrame(const string& dumpFileName);

void dumpIO_p2(System& system, ifstream& dumpFilein);

void Profile_velocity(System& system, ifstream& dumpFilein, double zlol, double zhil, int nslice);
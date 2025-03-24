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

using namespace std;

void Profile_velocity(System& system, ifstream& dumpFilein, double zlol, double zhil, int nslice);
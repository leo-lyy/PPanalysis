#pragma once
#include <iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include<chrono>
#include<cmath>
#include <iomanip>
#include"system.h"
#include"functions.h"

using namespace std;
void dumpIO_Ree(System& system, ifstream& dumpFilein);
void reeCalculation(System& system, long int f, ofstream& Reeout);

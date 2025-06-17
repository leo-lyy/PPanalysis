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


void dumpIO_Rg(System& system, ifstream& dumpFilein, int tarchain, int endchain);
#pragma once
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <iostream>
#include <sstream>
#include <algorithm>

#include "system.h"
#include "data_types.h"
#include "functions.h"
using namespace std;

bool compareAtomID(const Atom& a, const Atom& b);
void readGeneralInfo(ifstream& file, System& system);
void readBox(ifstream& file, System& system);
void readMasses(ifstream& file, System& system);
void readAtoms(ifstream& file, System& system);
void readVelocities(ifstream& file, System& system);
void readBonds(ifstream& file, System& system);
void readAngles(ifstream& file, System& system);
void readDihedrals(ifstream& file, System& system);
void readImpropers(ifstream& file, System& system);
void readLammpsData(const string& dataFileName, System& system);


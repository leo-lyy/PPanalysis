
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include<chrono>
#include<cmath>
#include<algorithm>
#include<stdexcept>
#include"../include/system.h"
#include"../include/readLMPdata.h"
#include"../include/readLMPdump.h"
#include"../include/p2calculation.h"
#include"../include/functions.h"

using namespace std;
void showMenu()
{
    // ╝╚╔ ╗ ∥ =
    cout << "======= LAMMPS Trajectory Analysis Program for Polypropylene ========" << endl;
    cout << "|  Make sure the data file are written in the LAMMPS data format    |" << endl;
    cout << "|  through **write_data xxx.data nocoeff**, and the dump file is    |" << endl;
    cout << "|  written in the format of 'id mol type x y z ix iy iz vx vy vz'   |" << endl;
    cout << "=====================================================================" << endl;
    cout << "1. Calculate the orientarion order parameter P2;" << endl;
    cout << "2. Calculate the number of helix structure (crystall degree);" << endl;
    cout << "3. Calculate the velocity profile in z direction;" << endl;
    cout << "0. Exit" << endl;
    cout << "Please choose the analysis type: ";
}


// 清理输入缓冲区
void clearInputBuffer() {
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}
int main()
{
    System system;
    int choice;
    string dataFileName,dumpFileName;

    while (true)
    {
        showMenu();
        cin >> choice;

        // 检查输入是否有效
        if (cin.fail()) 
        {
            cout << "Invalid input! Please enter a number." << endl;
            clearInputBuffer();
            continue;
        }

        switch (choice)
        {
            case 1: // Calculate the orientarion order parameter P2
            {
                cout << "Please input the data file name: ";
                cin >> dataFileName;
                cout << "Please input the dump file name: ";
                cin >> dumpFileName;
                cout << "Reading files ..." << endl;
                readLammpsData(dataFileName, system);
                cout << "Counting the number of frames in the dump file ..." << endl;
                system.frames = countDumpFrame(dumpFileName);
                cout << "The dump file contains " << system.frames << " frames." << endl;
                cout << "Calculating the P2 ..." << endl;
                ifstream dumpFilein(dumpFileName);
                dumpIO_p2(system, dumpFilein);
                cout << "Done!" << endl;
                return 0;
            }
            case 2: // Calculate the number of helix structure (crstall degree)
            {
                cout << "Developing ..." << endl;
                return 0;
            }
            case 3: // Calculate the velocity profile in z direction
            {
                cout << "Please input the data file name: ";
                cin >> dataFileName;
                cout << "Please input the dump file name: ";
                cin >> dumpFileName;
                cout << "Reading files ..." << endl;
                readLammpsData(dataFileName, system);
                
                double zlol, zhil;
                int nslice;
                cout << endl << "Please enter Z_low_limit, Z_high_limit and n_slice to analyse the velocity:"<<endl;
                cin >> zlol >> zhil >> nslice;

                cout << "Counting the number of frames in the dump file ..." << endl;
                system.frames = countDumpFrame(dumpFileName);
                cout << "The dump file contains " << system.frames << " frames." << endl;
                cout << "Calculating the P2 ..." << endl;
                ifstream dumpFilein(dumpFileName);
                Profile_velocity(system, dumpFilein, zlol, zhil, nslice);
                cout << "Done!" << endl;


                return 0;
            }
            
            case 0: // Exit
            {
                cout << "Goodbye!" << endl;
                return 0;
            }
            default:
            {
                cout << "Invalid input!" << endl;
                break;
            }

        }

    
    }
    return 0;
}

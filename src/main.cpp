
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
#include"../include/p2_calculation.h"
#include"../include/functions.h"
#include"../include/profileV_calculation.h"
#include"../include/helix_calculation.h"
#include"../include/rg_calculation.h"
#include"../include/endtoend_calculation.h"

using namespace std;
void showMenu()
{
/*
+----------------------------------------------------------------+
|             _                                                  |
|       __   | |           __                                    |
|      | _|  | |          |_ |                                   |
|      | |   | |           | |                                   |
|      | |   | |           | |                                   |
|      | |   | |           | | __                                |
|      | |   | |           | |/ /                                |
|      | |  /   \          | / /                                 |
|      | | / / \ \         |/ /    LAMMPS Trajectory Analysis    |
|      | |/ /   \ \        / /                                   |
|      | / /     \ \      / /|     Program for Polypropylene     |
|      |/ /       \ \    / / |                                   |
|      / /         \ \  / /| |   _ ___                           |
|     / /|          \ \/ / | |  | `___ \                         |
|    / / |           \__/  | |  | |   | |                        |
|   /_/| |                 | |  | |   | |                        |
|      |__|               |__|  | |   |__|                       |
|                                                                |
|================================================================|
| Make sure the data file are written in the LAMMPS data format  |
| through **write_data xxx.data nocoeff**, and the dump file is  |
| written in the format of 'id mol type x y z ix iy iz vx vy vz' |
+----------------------------------------------------------------+
 */ 
    cout << "+----------------------------------------------------------------+" << endl;
    cout << "|             _                                                  |" << endl;
    cout << "|       __   | |           __                                    |" << endl;
    cout << "|      | _|  | |          |_ |                                   |" << endl;
    cout << "|      | |   | |           | |                                   |" << endl;
    // cout << "|      | |   | |           | |                                   |" << endl;
    cout << "|      | |   | |           | | __                                |" << endl;
    cout << "|      | |   | |           | |/ /                                |" << endl;
    cout << "|      | |  /   \\          | / /                                 |" << endl;
    cout << "|      | | / / \\ \\         |/ /    LAMMPS Trajectory Analysis    |" << endl;
    cout << "|      | |/ /   \\ \\        / /                                   |" << endl;
    cout << "|      | / /     \\ \\      / /|     Program for Polypropylene     |" << endl;
    cout << "|      |/ /       \\ \\    / / |                                   |" << endl;
    cout << "|      / /         \\ \\  / /| |   _ ___                           |" << endl;
    cout << "|     / /|          \\ \\/ / | |  | `___ \\                         |" << endl;
    cout << "|    / / |           \\__/  | |  | |   | |                        |" << endl;
    cout << "|   /_/| |                 | |  | |   | |                        |" << endl;
    cout << "|      |__|               |__|  | |   |__|                       |" << endl;
    cout << "|                                                                |" << endl;
    cout << "|================================================================|" << endl;
    cout << "| Make sure the data file are written in the LAMMPS data format  |" << endl;
    cout << "| through **write_data xxx.data nocoeff**, and the dump file is  |" << endl;
    cout << "| written in the format of 'id mol type x y z ix iy iz vx vy vz' |" << endl;
    cout << "+----------------------------------------------------------------+" << endl;

    // Print the menu options
    cout << "  [1] Calculate the orientarion order parameter P2;                " << endl;
    cout << "  [2] Calculate the number of helix structure (crystall degree);   " << endl;
    cout << "  [3] Calculate the velocity profile in z direction(velocity);     " << endl;
    cout << "  [4] Calculate the velocity profile in z direction(displacement); " << endl;
    cout << "  [5] Calculate the radius of gyration;                            " << endl;
    cout << "  [6] Calculate the <end to end distance>;                           " << endl;
    cout << "                                                                   " << endl;
    cout << "  [0] Exit                                                         " << endl;
    cout << "------------------------------------------------------------------" << endl;
    cout << "Please choose the analysis type:                                   " << endl;
}


// 清理输入缓冲区
void clearInputBuffer() {
    cin.clear();
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
}
class FileError : public std::runtime_error {
    public:
        explicit FileError(const std::string& message) 
            : std::runtime_error(message) {}
    };
// Function to process file with error checking
void processFile(const string& filename) {
    ifstream inputFile(filename);
    
    try {
        // Check if file opened successfully
        if (!inputFile.is_open()) {
            throw FileError("Failed to open file: " + filename);
        }
        
        // Check if file is empty
        if (inputFile.peek() == ifstream::traits_type::eof()) {
            throw FileError("File is empty: " + filename);
        }
        // cout << "File processed successfully!" << endl;
    }
    catch (const FileError& e) {
        cerr << "Error: " << e.what() << endl;
        throw; // Re-throw to let caller handle if needed
    }
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
                processFile(dataFileName);
                cout << "Please input the dump file name: ";
                cin >> dumpFileName;
                processFile(dumpFileName);
                cout << "Reading data files ..." << endl;
                readLammpsData(dataFileName, system);
                cout << "Counting the number of frames in the dump file ..." << endl;
                system.frames = countDumpFrame(dumpFileName);
                cout << "The dump file contains " << system.frames << " frames." << endl;
                cout << "Calculating the P2 ..." << endl;
                ifstream dumpFilein(dumpFileName);
                dumpIO_p2(system, dumpFilein);
                cout << endl << "Done!" << endl;
                return 0;
            }
            case 2: // Calculate the number of helix structure (crstall degree)
            {
                cout << "Please input the data file name: ";
                cin >> dataFileName;
                processFile(dataFileName);
                cout << "Please input the dump file name: ";
                cin >> dumpFileName;
                processFile(dumpFileName);
                cout << "Reading data files ..." << endl;
                readLammpsData(dataFileName, system);
                cout << "Counting the number of frames in the dump file ..." << endl;
                system.frames = countDumpFrame(dumpFileName);
                cout << "The dump file contains " << system.frames << " frames." << endl;
                cout << "Calculating the number of helix beads..." << endl;
                ifstream dumpFilein(dumpFileName);
                ofstream dumpFileout(dumpFileName.substr(0, dumpFileName.find_last_of(".")) + "_crystal.dump");
                dumpIO_helix(system, dumpFilein, dumpFileout);
                cout << endl << "Done!" << endl;
                return 0;
            }
            case 3: // Calculate the velocity profile in z direction(read the real instant velocity)
            {
                int mode = 0; // 0 for velocity, 1 for displacement
                cout << "Please input the data file name: ";
                cin >> dataFileName;
                processFile(dataFileName);
                cout << "Please input the dump file name: ";
                cin >> dumpFileName;
                processFile(dumpFileName);
                cout << "Reading data files ..." << endl;
                readLammpsData(dataFileName, system);
                
                double zlol, zhil;
                int nslice;
                cout << endl << "Please enter Z_low_limit, Z_high_limit and n_slice to analyse the velocity:"<<endl;
                cin >> zlol >> zhil >> nslice;

                cout << "Counting the number of frames in the dump file ..." << endl;
                system.frames = countDumpFrame(dumpFileName);
                cout << "The dump file contains " << system.frames << " frames." << endl;
                // cout << "Calculating the P2 ..." << endl;
                ifstream dumpFilein(dumpFileName);
                Profile_velocity(system, dumpFilein, zlol, zhil, nslice, mode);
                cout << endl << "Done!" << endl;


                return 0;
            }
            case 4: // Calculate the velocity profile in z direction(calculate the velocity from displacement)
            {
                int mode = 1; // 0 for velocity, 1 for displacement
                cout << "Please input the data file name: ";
                cin >> dataFileName;
                processFile(dataFileName);
                cout << "Please input the dump file name: ";
                cin >> dumpFileName;
                processFile(dumpFileName);
                cout << "Reading data files ..." << endl;
                readLammpsData(dataFileName, system);
                
                double zlol, zhil;
                int nslice;
                cout << endl << "Please enter Z_low_limit, Z_high_limit and n_slice to analyse the velocity:"<<endl;
                cin >> zlol >> zhil >> nslice;

                cout << "Counting the number of frames in the dump file ..." << endl;
                system.frames = countDumpFrame(dumpFileName);
                cout << "The dump file contains " << system.frames << " frames." << endl;
                // cout << "Calculating the P2 ..." << endl;
                ifstream dumpFilein(dumpFileName);
                Profile_velocity(system, dumpFilein, zlol, zhil, nslice, mode);
                cout << endl << "Done!" << endl;

                return 0;
            }
            case 5: // Calculate the radius of gyration
            {
                cout << "Please input the data file name: ";
                cin >> dataFileName;
                cout << "Please input the dump file name: ";
                cin >> dumpFileName;
                int tarchain, endchain;
                cout << "Please enter the start and the end chain number to calculate the Mean Square Radious of Gyration:"<<endl;
                cin >> tarchain >> endchain;
                cout << "Reading files ..." << endl;
                readLammpsData(dataFileName, system);
                cout << "Counting the number of frames in the dump file ..." << endl;
                system.frames = countDumpFrame(dumpFileName);
                cout << "The dump file contains " << system.frames << " frames." << endl;
                cout << "Calculating the radius of gyration ..." << endl;
                ifstream dumpFilein(dumpFileName);
                dumpIO_Rg(system, dumpFilein, tarchain, endchain);
                // dumpFilein.close();
                cout << endl << "Done!" << endl;
                return 0;
            }
            case 6: // Calculate the end to end distance
            {
                cout << "Please input the data file name: ";
                cin >> dataFileName;
                cout << "Please input the dump file name: ";
                cin >> dumpFileName;
                // int tarchain, endchain;
                // cout << "Please enter the start and the end chain number to calculate End to End Distance:"<<endl;
                // cin >> tarchain >> endchain;
                cout << "Reading files ..." << endl;
                readLammpsData(dataFileName, system);
                cout << "Counting the number of frames in the dump file ..." << endl;
                system.frames = countDumpFrame(dumpFileName);
                cout << "The dump file contains " << system.frames << " frames." << endl;
                cout << "Calculating the end to end distance ..." << endl;
                ifstream dumpFilein(dumpFileName);
                // dumpIO_Ree(system, dumpFilein, tarchain, endchain);
                dumpIO_Ree(system, dumpFilein);
                // dumpFilein.close();
                cout << endl << "Done!" << endl;
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

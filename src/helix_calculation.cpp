#include "../include/helix_calculation.h"

using namespace std;

const double p = 2.1666667; // one third of the iPP helical repeat unit
const double sigma = 3.95;

void helixCalculation(System& system)
{
    vector<vector<axis>> r1, r3, r5, rh, bh; //  rh: postion of helix   bh: orientation of helix;
    vector<long int> helixNum(system.molnumMAX + 1, 0);// calculate the helix number in each chain
    cout << "Initializing the helix calculation ..." << endl;
    for (long int i = 1; i <= system.molnumMAX; i++)
    {
        // the number of the helix in each chain should be (the number of monomers - 2) * 2
        helixNum[i] = (system.molen[i] / 3 - 2) * 2;
    }
    long int maxhelix = *max_element(helixNum.begin(), helixNum.end());
    // cout << "The max helix number is: " << maxhelix << endl;
    // calculate the intermolecular oder parameter: Crystallinity chic
    vecInit(r1, system.molnumMAX, maxhelix + 1);
    vecInit(r3, system.molnumMAX, maxhelix + 1);
    vecInit(r5, system.molnumMAX, maxhelix + 1);
    vecInit(rh, system.molnumMAX, maxhelix + 1);
    vecInit(bh, system.molnumMAX, maxhelix + 1);

    vector<vector<long int>> helixNeighbor(system.molnumMAX + 1, vector<long int>(maxhelix + 1, 0));
    axis rminus, XX, YY, ZZ;
    long int plen = 0;
    long int k = 0;
    cout << "Calculating the helix position and orientation ..." << endl;
    for(long int i = 0; i < system.num_atoms; i++)
    {
        if(i % 100 == 0) printProgressBar(i + 100, system.num_atoms);
        long int molid = system.atoms[i].mol;
        if(i + 1 < system.num_atoms && system.atoms[i + 1].mol != molid) // next chain
        {
            plen += system.molen[molid];
            k = 0;
        }
        if (molid == 0) continue;  // ignore the wall atoms
        long int atomid = i + 1; // the atom id in the whole box, unique
        long int parid = atomid - plen; // particle id in the chain
        if ((i) % 3 != 2 && parid <= system.molen[molid] - 6 ) // - 6 : ignore the last methyl group , - 5 : consider the last methyl group, but bugs may occur
        {
            k++;
            long int p1 = i;
            long int p3 = i + 3;
            long int p5 = i + 6;
            // if (p5 >= system.molen[molid] ) p5--; // ignore the last methyl group
            // fout << molid <<" : "<< p1 + 1 <<" "<< p3 + 1 <<" "<< p5 + 1 <<endl;
            r1[molid][k].x = system.atoms[p1].x;
            r1[molid][k].y = system.atoms[p1].y;
            r1[molid][k].z = system.atoms[p1].z;

            r3[molid][k].x = system.atoms[p3].x;
            r3[molid][k].y = system.atoms[p3].y;
            r3[molid][k].z = system.atoms[p3].z;

            r5[molid][k].x = system.atoms[p5].x;
            r5[molid][k].y = system.atoms[p5].y;
            r5[molid][k].z = system.atoms[p5].z;
            // fout << i + 1 << endl;
            // fout << "r1: " << r1[molid][k].x <<" "<< r1[molid][k].y <<" "<< r1[molid][k].z << endl;
            // fout << "r3: " << r3[molid][k].x <<" "<< r3[molid][k].y <<" "<< r3[molid][k].z << endl;
            // fout << "r5: " << r5[molid][k].x <<" "<< r5[molid][k].y <<" "<< r5[molid][k].z << endl;
            rh[molid][k].x = (r1[molid][k].x + r3[molid][k].x + r5[molid][k].x) / 3;
            rh[molid][k].y = (r1[molid][k].y + r3[molid][k].y + r5[molid][k].y) / 3;
            rh[molid][k].z = (r1[molid][k].z + r3[molid][k].z + r5[molid][k].z) / 3;

            // wrap the rh into the box, wrap the helix center postions into the box
            rh[molid][k].x = fmod(rh[molid][k].x, system.bx);
            rh[molid][k].y = fmod(rh[molid][k].y, system.by);
            rh[molid][k].z = fmod(rh[molid][k].z, system.bz);
            rh[molid][k].id = system.atoms[p1].id;

            rminus.x = r5[molid][k].x - r1[molid][k].x;
            rminus.y = r5[molid][k].y - r1[molid][k].y;
            rminus.z = r5[molid][k].z - r1[molid][k].z;
            // fout <<"|r5 - r1| = "<< norm(rminus)<<endl;
            
            double cosphi = 2 * p / norm(rminus);
            // fout << "cosphi = " << cosphi <<endl;
            
            if (cosphi > 1)
            {
                bh[molid][k] = {-1, -1, -1};
                continue;
            }
            
            bh[molid][k].x = cosphi;
            bh[molid][k].y = sin(acos(cosphi));
            bh[molid][k].z = 0;
            // fout << "bh(XYZ): " << bh[molid][k].x <<" "<< bh[molid][k].y <<" "<< bh[molid][k].z <<endl;

            XX.x = r5[molid][k].x - r1[molid][k].x;
            XX.y = r5[molid][k].y - r1[molid][k].y;
            XX.z = r5[molid][k].z - r1[molid][k].z;
            ZZ.x = r3[molid][k].x - rh[molid][k].x;
            ZZ.y = r3[molid][k].y - rh[molid][k].y;
            ZZ.z = r3[molid][k].z - rh[molid][k].z;
            cross(YY, XX, ZZ);  // YY = XX cdot ZZ 
            // fout << "XX: " << XX.x <<" "<< XX.y <<" "<< XX.z <<endl;
            // fout << "YY: " << YY.x <<" "<< YY.y <<" "<< YY.z <<endl;
            // fout << "ZZ: " << ZZ.x <<" "<< ZZ.y <<" "<< ZZ.z <<endl;
            // fout << "X * Z = " << XX.x * ZZ.x + XX.y * ZZ.y + XX.z * ZZ.z <<endl;
            // fout << YY.x <<" "<< YY.y <<" "<< YY.z <<endl;

            XX = unitVec(XX);
            YY = unitVec(YY);
            ZZ = unitVec(ZZ);
            // fout << "u_XX: " << XX.x <<" "<< XX.y <<" "<< XX.z <<endl;
            // fout << "u_YY: " << YY.x <<" "<< YY.y <<" "<< YY.z <<endl;
            // fout << "u_ZZ: " << ZZ.x <<" "<< ZZ.y <<" "<< ZZ.z <<endl;

            double original_x = bh[molid][k].x;
            double original_y = bh[molid][k].y;
            double original_z = bh[molid][k].z;
            bh[molid][k].x = original_x * XX.x + original_y * YY.x + original_z * ZZ.x;
            bh[molid][k].y = original_x * XX.y + original_y * YY.y + original_z * ZZ.y;
            bh[molid][k].z = original_x * XX.z + original_y * YY.z + original_z * ZZ.z;
            // fout << "bh: " << bh[molid][k].x <<" "<< bh[molid][k].y <<" "<< bh[molid][k].z <<endl;
            bh[molid][k].id = system.atoms[p1].id;
            // fout << "======================" <<endl;
        }

        
    }
    // cout << endl << "Done!" << endl;
    //////////// Searching the neighbour helix/////////////
    cout <<endl<< "Search the neighbouring helix..." << endl;
    // 四层嵌套循环，外两层循环链A及其螺旋，内两层循环链B及其螺旋，搜索AB之间的相邻螺旋
    for(long int i = 1; i <= system.molnumMAX; i++)
    {
        long int helixlen = helixNum[i];
        printProgressBar(i, system.molnumMAX);
        // cout<< "Processing chain " << i << " ..." << endl;
        for (long int j = 1; j <= helixlen ; j++)   ////////////////////////////////////////
        {
            // if(calMode == 2) cout <<"i j : "<< i <<"  "<< j <<endl;
            if (bh[i][j].x == -1 || bh[i][j].y == -1 || bh[i][j].z == -1 || bh[i][j].x * bh[i][j].y * bh[i][j].z == 0) continue;
            for (long int ii = i; ii <= system.molnumMAX; ii++)
            {
                long int helixlen2 = helixNum[ii];
                for(long int jj = j + 1; jj <= helixlen2 ; jj++)    /////////////////////////////
                {
                    // if(calMode == 2 && i==50) cout <<"bh(i,j) = "<< bh[i][j].x <<" "<< bh[i][j].y <<" "<< bh[i][j].z <<"  //  "<< "bh(ii,jj) = "<< bh[ii][jj].x <<" "<< bh[ii][jj].y <<" "<< bh[ii][jj].z <<endl;
                    // if (i == ii && j == jj ) continue;
                    if (bh[ii][jj].x == -1 || bh[ii][jj].y == -1 || bh[ii][jj].z == -1 || bh[ii][jj].x * bh[ii][jj].y * bh[ii][jj].z == 0) continue;
                    double deghelix = degVec(bh[i][j], bh[ii][jj]);
                    if (deghelix < 30)
                    {
                        double drx = rh[i][j].x - rh[ii][jj].x;
                        double dry = rh[i][j].y - rh[ii][jj].y;
                        double drz = rh[i][j].z - rh[ii][jj].z;
                        double rxz = sqrt(drx * drx + drz * drz);
                        if ((1.4 * sigma < rxz && rxz < 1.9 * sigma) && fabs(dry) < p / 2)
                        {
                            // This is a crystal
                            helixNeighbor[i][j]++;
                            helixNeighbor[ii][jj]++;
                        }
                    }
                }
            }
        }
    }
    for(long int i = 1; i <= system.molnumMAX; i++)
    {
        long int helixlen = helixNum[i];
        for (long int j = 0; j < helixlen; j++)
        {
            // hout << helixNeighbor[i][j] <<"    ";
            if(helixNeighbor[i][j] >= 2)
            {
                // system.atoms[rh[i][j].id - 1].crystal = true;
                // system.atoms[rh[i][j].id + 3 - 1].crystal = true;
                // system.atoms[rh[i][j].id + 6 - 1].crystal = true;
                for(int k = -1; k < 6; k++) // 7 beads in a helix
                {
                    system.atoms[rh[i][j].id + k].crystal = true;
                }
            }
        }
    }
    cout << endl << "Done!" << endl;

    
}

void dumpIO_helix(System& system, ifstream& dumpFilein, ofstream& dumpFileout)
{
    string line;
    istringstream ss(line);
    ofstream crydeg("crystal_degree.txt");
    // ofstream cryzpos("crystal_position.txt");
    ofstream outHistCryzPos("hist_crystal_zpos.txt");
    ofstream outHistzPos("hist_zpos.txt");
    long int totalFrame = system.frames;
    long int f = 0;
    while (f < totalFrame)        // read the dump file frame by frame
    {
        auto fstart = std::chrono::high_resolution_clock::now();// 记录开始的时间
        f++;
        cout << endl << "Calculating frame " << f << "/" << totalFrame << " ..." << endl;
        long int id, flable;
        int mol, type;
        double x, y, z, vx, vy, vz;
        long int ix, iy, iz;
        vector<double> zPos;    //store the z position of atoms
        vector<double> cryzPos; //store the z position of crystal atoms
        for (long int k = 0; k < 9; k++)
        {
            getline(dumpFilein, line);
            dumpFileout << line << endl;
            // cout << line << endl;
            if (line == "ITEM: TIMESTEP")
            {
                getline(dumpFilein, line);
                dumpFileout << line << endl;
                ss >> flable;
                k++;
            }
            else if (line == "ITEM: BOX BOUNDS pp pp pp")
            {
                getline(dumpFilein, line);
                dumpFileout << line << endl;
                ss >> system.xlo >> system.xhi;
                getline(dumpFilein, line);
                dumpFileout << line << endl;
                ss >> system.ylo >> system.yhi;
                getline(dumpFilein, line);
                dumpFileout << line << endl;
                ss >> system.zlo >> system.zhi;
                system.bx = system.xhi - system.xlo;
                system.by = system.yhi - system.ylo;
                system.bz = system.zhi - system.zlo;
                k+=3;
            }
            else if (line == "ITEM: BOX BOUNDS xy xz yz pp pp pp")
            {
                system.box_type = 1; // triclinic box
                getline(dumpFilein, line);
                ss >> system.xlo >> system.xhi >> system.xtilt;
                getline(dumpFilein, line);
                ss >> system.ylo >> system.yhi >> system.ytilt;
                getline(dumpFilein, line);
                ss >> system.zlo >> system.zhi >> system.ztilt;
                system.bx = system.xhi - system.xlo;
                system.by = system.yhi - system.ylo;
                system.bz = system.zhi - system.zlo;
                k += 3;
            }
            else if (line == "ITEM: ATOMS id mol type x y z ix iy iz")
            {
                for (long int i = 0; i < system.num_atoms; i++)
                {
                    getline(dumpFilein, line);
                    istringstream ss(line);
                    ss >> id >> mol >> type >> x >> y >> z >> ix >> iy >> iz;
                    id--;
                    system.atoms[id].x = x;
                    system.atoms[id].y = y;
                    system.atoms[id].z = z;
                    system.atoms[id].mol = mol;
                    system.atoms[id].type = type;
                    system.atoms[id].ix = ix;
                    system.atoms[id].iy = iy;
                    system.atoms[id].iz = iz;
                    system.atoms[id].crystal = false;
                }
                unwrap(system);
            }
            else if (line == "ITEM: ATOMS id mol type x y z ix iy iz vx vy vz")
            {
                for (long int i = 0; i < system.num_atoms; i++)
                {
                    getline(dumpFilein, line);
                    istringstream ss(line);
                    ss >> id >> mol >> type >> x >> y >> z >> ix >> iy >> iz >> vx >> vy >> vz;
                    id--;
                    system.atoms[id].x = x;
                    system.atoms[id].y = y;
                    system.atoms[id].z = z;
                    system.atoms[id].mol = mol;
                    system.atoms[id].type = type;
                    system.atoms[id].ix = ix;
                    system.atoms[id].iy = iy;
                    system.atoms[id].iz = iz;
                    system.atoms[id].vx = vx;
                    system.atoms[id].vy = vy;
                    system.atoms[id].vz = vz;
                    system.atoms[id].crystal = false;
                }
                unwrap(system);
            }
        }
        helixCalculation(system);
        cout << "Writing frame " << f << " ..." << endl;

        dumpFileout.flush();  // 强制刷新输出流
        wrap(system);
        long int crystalCount = 0;
        for (long int i = 0; i < system.num_atoms; i++)
        {
            // store all iPP atom z positions 
            if (system.atoms[i].type != 4) zPos.push_back(system.atoms[i].z);
            if (system.atoms[i].crystal) 
            {
                crystalCount++;
                dumpFileout << system.atoms[i].id << " " << system.atoms[i].mol << " " << 5 << " " 
                << system.atoms[i].x << " " << system.atoms[i].y << " " << system.atoms[i].z << " " 
                << system.atoms[i].ix << " " << system.atoms[i].iy << " " << system.atoms[i].iz <<" "
                << system.atoms[i].vx << " " << system.atoms[i].vy << " " << system.atoms[i].vz << endl;
                // cryzpos << system.atoms[i].z << endl;
                cryzPos.push_back(system.atoms[i].z);   // store crystal iPP atoms z positions
                
            }
            else dumpFileout << system.atoms[i].id << " " << system.atoms[i].mol << " " << system.atoms[i].type << " " 
            << system.atoms[i].x << " " << system.atoms[i].y << " " << system.atoms[i].z << " " 
            << system.atoms[i].ix << " " << system.atoms[i].iy << " " << system.atoms[i].iz << " "
            << system.atoms[i].vx << " " << system.atoms[i].vy << " " << system.atoms[i].vz << endl;
        }
        crydeg << f <<" "<< crystalCount << endl;
        // make a histogram of zPos and cryzPos, with each bin width close to 1.0 (the number of bins = box length in z direction)
        int numBins = static_cast<int>(system.bz) + 1;
        vector<int> histZPos(numBins, 0);
        vector<int> histCryZPos(numBins, 0);
        for (double zval : zPos)
        {
            // consider the box limit is below 0, so the zpos should be shifted by the zlo of the box
            int binIndex = static_cast<int>(zval - system.zlo);
            if (binIndex >= 0 && binIndex < numBins)
            {
                histZPos[binIndex]++;
            }
        }
        for (double cryzval : cryzPos)
        {
            int binIndex = static_cast<int>(cryzval - system.zlo);
            if (binIndex >= 0 && binIndex < numBins)
            {
                histCryZPos[binIndex]++;
            }
        }
        // output the histogram data to files, the first line is bin centers, the rest lines are counts for each frame
        if (f == 1)
        {
            for (int b = 0; b < numBins; b++)
            {
                double binCenter = b + 0.5;
                outHistzPos << binCenter << " ";
                outHistCryzPos << binCenter << " ";
            }
            outHistzPos << endl;
            outHistCryzPos << endl;
        }
        for (int b = 0; b < numBins; b++)
        {
            outHistzPos << histZPos[b] << " ";
            outHistCryzPos << histCryZPos[b] << " ";
        }
        outHistzPos << endl;
        outHistCryzPos << endl;

        crydeg.flush();
        dumpFileout.flush();  // 强制刷新输出流

        auto fend = std::chrono::high_resolution_clock::now();// 记录读取结束的时间
        std::chrono::duration<double> elapsed = fend - fstart;
        double remainTime = double(totalFrame - f) * elapsed.count();
        cout << "Elapsed time for frame " << f <<" : " << elapsed.count() << " s" << ", approximate remaining time: " << remainTime <<" s."<< endl;
    }
}
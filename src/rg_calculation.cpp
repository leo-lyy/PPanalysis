# include "../include/rg_calculation.h"
using namespace std;

void dumpIO_Rg(System& system, ifstream& dumpFilein, int tarchain, int endchain)
{
    vector<double> molMass(system.molnumMAX+1,0);
    vector<double> msrg(system.molnumMAX+1,0);
    vector<vector<double>> com(system.molnumMAX+1, vector<double>(4,0));
    string line;
    ofstream msrgout("msrg.txt");
    ofstream comXout("comX.txt");
    ofstream comYout("comY.txt");
    ofstream comZout("comZ.txt");

    long int f = 0;
    while (f < system.frames)
    {
        f++;
        long long timestep = -1; bool gotAtoms=false; bool gotBox=false; bool gotTimestep=false; 
        // Read frame header + atom section
        while (std::getline(dumpFilein, line)) {
            if (line.rfind("ITEM: TIMESTEP",0)==0) {
                if(!std::getline(dumpFilein,line)) return; { istringstream iss(line); iss >> timestep; }
                gotTimestep=true;
            } else if (line.rfind("ITEM: NUMBER OF ATOMS",0)==0) {
                if(!std::getline(dumpFilein,line)) return; /* optional parse */
            } else if (line.rfind("ITEM: BOX BOUNDS",0)==0) {
                bool triclinic = (line.find("xy")!=string::npos);
                system.box_type = triclinic ? 1 : 0;
                if(!std::getline(dumpFilein,line)) return; { istringstream iss(line); if(triclinic) iss >> system.xlo >> system.xhi >> system.xtilt; else iss >> system.xlo >> system.xhi; }
                if(!std::getline(dumpFilein,line)) return; { istringstream iss(line); if(triclinic) iss >> system.ylo >> system.yhi >> system.ytilt; else iss >> system.ylo >> system.yhi; }
                if(!std::getline(dumpFilein,line)) return; { istringstream iss(line); if(triclinic) iss >> system.zlo >> system.zhi >> system.ztilt; else iss >> system.zlo >> system.zhi; }
                system.bx = system.xhi - system.xlo; system.by = system.yhi - system.ylo; system.bz = system.zhi - system.zlo;
                gotBox=true;
            } else if (line.rfind("ITEM: ATOMS",0)==0) {
                bool hasVel = (line.find("vx")!=string::npos);
                for (long long c=0;c<system.num_atoms;c++) {
                    long long id; int mol,type; double x,y,z; int ix,iy,iz; double vx,vy,vz; 
                    if(hasVel) { if(!(dumpFilein>>id>>mol>>type>>x>>y>>z>>ix>>iy>>iz>>vx>>vy>>vz)) return; }
                    else { if(!(dumpFilein>>id>>mol>>type>>x>>y>>z>>ix>>iy>>iz)) return; }
                    if(id<=0 || id>system.num_atoms) return; 
                    size_t idx = static_cast<size_t>(id-1);
                    auto &a = system.atoms[idx];
                    a.mol=mol; a.type=type; a.x=x; a.y=y; a.z=z; a.ix=ix; a.iy=iy; a.iz=iz; 
                    if (f==1) molMass[mol] += system.mass[type];
                }
                gotAtoms=true;
                break;
            }
            if(gotAtoms) break; 
        }
        if(!gotAtoms) break; 
        if(!(gotBox && gotTimestep)) { break; }
        if(system.bx<=0||system.by<=0||system.bz<=0) { return; }

        unwrap(system);

        for (int i = tarchain; i <= endchain; i++) { com[i][1]=com[i][2]=com[i][3]=0.0; msrg[i]=0.0; }
        for (long int i = 0; i < system.num_atoms; i++) {
            int mol = system.atoms[i].mol;
            if (mol >= tarchain && mol <= endchain) {
                double m = system.mass[ system.atoms[i].type ];
                com[mol][1] += system.atoms[i].x * m;
                com[mol][2] += system.atoms[i].y * m;
                com[mol][3] += system.atoms[i].z * m;
            }
        }
        for (int i = tarchain; i <= endchain; i++) if(molMass[i]>0){ com[i][1]/=molMass[i]; com[i][2]/=molMass[i]; com[i][3]/=molMass[i]; }
        for (long int j=0;j<system.num_atoms;j++) {
            int mol = system.atoms[j].mol; if(mol < tarchain || mol> endchain) continue; 
            double m = system.mass[ system.atoms[j].type ];
            double dx=system.atoms[j].x-com[mol][1]; double dy=system.atoms[j].y-com[mol][2]; double dz=system.atoms[j].z-com[mol][3];
            msrg[mol] += m*(dx*dx+dy*dy+dz*dz);
        }
        for (int i=tarchain;i<=endchain;i++) if(molMass[i]>0) msrg[i]=sqrt(msrg[i]/molMass[i]); else msrg[i]=0.0;
        for (int i=tarchain;i<=endchain;i++) {
            if(system.bx>0){ while(com[i][1]<system.xlo) com[i][1]+=system.bx; while(com[i][1]>=system.xhi) com[i][1]-=system.bx; }
            if(system.by>0){ while(com[i][2]<system.ylo) com[i][2]+=system.by; while(com[i][2]>=system.yhi) com[i][2]-=system.by; }
            if(system.bz>0){ while(com[i][3]<system.zlo) com[i][3]+=system.bz; while(com[i][3]>=system.zhi) com[i][3]-=system.bz; }
        }
        msrgout<<fixed<<setprecision(0)<<f<<"        ";
        comXout<<fixed<<setprecision(0)<<f<<"        ";
        comYout<<fixed<<setprecision(0)<<f<<"        ";
        comZout<<fixed<<setprecision(0)<<f<<"        ";
        msrgout<<fixed<<setprecision(5); comXout<<fixed<<setprecision(5); comYout<<fixed<<setprecision(5); comZout<<fixed<<setprecision(5);
        for(int i=tarchain;i<=endchain;i++){ msrgout<<msrg[i]<<"      "; comXout<<com[i][1]<<"      "; comYout<<com[i][2]<<"      "; comZout<<com[i][3]<<"      "; }
        msrgout<<endl; comXout<<endl; comYout<<endl; comZout<<endl;
        printProgressBar(f, system.frames);
    }
    msrgout.close(); comXout.close(); comYout.close(); comZout.close();
}

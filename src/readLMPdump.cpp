
#include"../include/readLMPdump.h"

using namespace std;


long int countDumpFrame(const string& dumpFileName)
{
    ifstream dumpFile(dumpFileName);
    string line;
    long int frameCount = 0;
    while (getline(dumpFile, line))
    {
        if (line == "ITEM: TIMESTEP")
        {
            frameCount++;
        }
    }
    dumpFile.close();
    return frameCount;
}


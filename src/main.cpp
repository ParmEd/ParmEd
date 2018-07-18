#include <iostream>
#include "readparm.h"

using namespace std;

int main() {

    for (int i = 0; i < 10; i++) {
        ParmDataMap parmData;
        ParmStringMap parmComments, unkParmData;
        ParmFormatMap parmFormats;
        string version;
        string fname = "trx.prmtop";
        vector<string> flagList;

        ExitStatus retval = readparm(fname, flagList, parmData, parmComments,
                                     unkParmData, parmFormats, version);
    }

    /*
    cout << "My flags are:" << endl;
    for (ParmDataMap::const_iterator it=parmData.begin(); it!=parmData.end(); it++) {
        cout << "\t" << it->first << " size = " << (it->second).size() << endl;
    }
    */
    return 0;
}

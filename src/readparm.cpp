/* readparm.cpp
 *
 * This is optimized code written in C++ designed to speed up topology file
 * reading compared to what is possible in pure Python.
 */
#include <cstdlib>    //< atof and atoi
#include <cstdio>     //< sscanf
#include <cstring>    //< strncpy
#include <fstream>
#include "readparm.h"

using namespace std;

/// Determine how the data is formatted based on the Fortran format statement
ParmDataType parseFormat(const string &fmt, int &ncols, int &width) {
    string up = upper(fmt);
    int i, j, k;
    if (sscanf(fmt.c_str(), "%dA%d", &i, &j) == 2) {
        // Must be characters, but in that case our data type only supports 4
        // characters, so mark it as unknown if the strings are longer
        if (j > 4) return UNKNOWN;
        // OK, we have a string format
        ncols = i;
        width = j;
        return HOLLERITH;
    } else if (sscanf(fmt.c_str(), "%dE%d.%d", &i, &j, &k) == 3) {
        // Floating point numbers
        ncols = i;
        width = j;
        return FLOAT;
    } else if (sscanf(fmt.c_str(), "%dI%d", &i, &j) == 2) {
        ncols = i;
        width = j;
        return INTEGER;
    } else if (sscanf(fmt.c_str(), "%d(F%d.%d)", &i, &j, &k) == 3) {
        ncols = i;
        width = j;
        return FLOAT;
    } else if (sscanf(fmt.c_str(), "%dF%d.%d", &i, &j, &k) == 3) {
        ncols = i;
        width = j;
        return FLOAT;
    }
    return UNKNOWN;
}

/// Parse the actual topology file and store all of the data in hash tables
ExitStatus readparm(const string &fname, vector<string> &flagList,
                    ParmDataMap &parmData, ParmStringMap &parmComments,
                    ParmStringMap &unkParmData, ParmFormatMap &parmFormats,
                    string &version) {

    string line;

    // First open the file, and make sure we can
    ifstream parm(fname.c_str());

    if (!parm)
        return NOOPEN;

    const string VERSIONFLAG = "%VERSION";
    const string COMMENTFLAG = "%COMMENT";
    const string FORMATFLAG = "%FORMAT";
    const string DATAFLAG = "%FLAG";

    const size_t VERSIONLEN = VERSIONFLAG.size();
    const size_t COMMENTLEN = COMMENTFLAG.size();
    const size_t FORMATLEN = FORMATFLAG.size();
    const size_t FLAGLEN = DATAFLAG.size();

    // The first line needs to start with %VERSION
    if (getline(parm, line).eof()) {
        parm.close();
        return EMPTY;
    }

    if (line.substr(0,VERSIONLEN) != VERSIONFLAG) {
        parm.close();
        return NOVERSION;
    }
    // Now get our version
    version = strip(line.substr(VERSIONLEN));

    // We have successfully parsed our version. Now parse the rest of the file
    string curflag = "";
    int ncols = -1;
    int width = -1;
    size_t linewidth = -1;
    ParmDataType curtype = UNKNOWN;
    ParmData d;
    string word;

    while (!getline(parm, line).eof()) {
        // RESIDUE_ICODE can have blank entries, so don't strip it...
        if (curflag != "RESIDUE_ICODE")
            line = rstrip(line);
        if (line.substr(0, FLAGLEN) == DATAFLAG) {
            // This is a new flag -- push the data back to 
            curflag = strip(line.substr(FLAGLEN));
            flagList.push_back(curflag);
            continue;
        }
        if (line.substr(0, COMMENTLEN) == COMMENTFLAG) {
            parmComments[curflag].push_back(strip(line.substr(COMMENTLEN)));
            continue;
        }
        if (line.substr(0, FORMATLEN) == FORMATFLAG) {
            size_t start = line.find_first_of('(') + 1;
            size_t end = line.find_last_of(')');
            string fmt = line.substr(start, end-start);
            curtype = parseFormat(fmt, ncols, width);
            ParmFormatType typ;
            typ.dataType = curtype;
            typ.fmt = fmt;
            parmFormats[curflag] = typ;
            linewidth = ncols * width;
            continue;
        }

        // This is where the actual data processing occurs
        switch (curtype) {
            case UNKNOWN:
                unkParmData[curflag].push_back(line);
                break;
            case FLOAT:
                for (size_t i = 0; i < line.size(); i += width) {
                    d.f = atof(line.substr(i, width).c_str());
                    parmData[curflag].push_back(d);
                }
                break;
            case INTEGER:
                for (size_t i = 0; i < line.size(); i += width) {
                    d.i = atoi(line.substr(i, width).c_str());
                    parmData[curflag].push_back(d);
                }
                break;
            case HOLLERITH:
                for (size_t i = 0; i < line.size(); i += width) {
                    word = strip(line.substr(i, width));
                    strncpy(d.c, word.c_str(), MAX_HOLLERITH_SIZE);
                    parmData[curflag].push_back(d);
                }
                break;
            default:
                // Should not reach here
                parm.close();
                return ERR;
                break;
        }
    }

    parm.close();

    return OK;
}

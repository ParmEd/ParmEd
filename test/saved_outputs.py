
"""
This module simply contains a repository of saved outputs so it doesn't flood
the actual test scripts
"""

PRINT_DETAILS = """
The mask @1 matches 1 atoms:

   ATOM    RES  RESNAME  NAME  TYPE   At.#   LJ Radius    LJ Depth      Mass    Charge GB Radius GB Screen
      1      1      SER     N    N3      7      1.8240      0.1700   14.0100    0.1849    1.5500    0.7900
"""

PRINT_BONDS = """\
             Atom 1              Atom 2       R eq   Frc Cnst
      1    N (  N3)       5   CA (  CT)     1.4710   367.0000
      1    N (  N3)       2   H1 (   H)     1.0100   434.0000
      1    N (  N3)       3   H2 (   H)     1.0100   434.0000
      1    N (  N3)       4   H3 (   H)     1.0100   434.0000
"""

PRINT_ANGLES = """\
             Atom 1               Atom 2               Atom 3   Frc Cnst   Theta eq
      1    N (  N3)        5   CA (  CT)        7   CB (  CT)    80.0000   111.2000
      1    N (  N3)        5   CA (  CT)       12    C (   C)    80.0000   111.2000
      4   H3 (   H)        1    N (  N3)        5   CA (  CT)    50.0000   109.5000
      3   H2 (   H)        1    N (  N3)        4   H3 (   H)    35.0000   109.5000
      3   H2 (   H)        1    N (  N3)        5   CA (  CT)    50.0000   109.5000
      2   H1 (   H)        1    N (  N3)        3   H2 (   H)    35.0000   109.5000
      2   H1 (   H)        1    N (  N3)        4   H3 (   H)    35.0000   109.5000
      2   H1 (   H)        1    N (  N3)        5   CA (  CT)    50.0000   109.5000
      1    N (  N3)        5   CA (  CT)        6   HA (  HP)    50.0000   109.5000
"""

PRINT_DIHEDRALS = """\
               Atom 1               Atom 2               Atom 3               Atom 4     Height  Periodic.      Phase  EEL Scale  VDW Scale
        1    N (  N3)        5   CA (  CT)        7   CB (  CT)       10   OG (  OH)     0.1556     3.0000     0.0000     1.2000     2.0000
        1    N (  N3)        5   CA (  CT)       12    C (   C)       13    O (   O)     0.0000     2.0000     0.0000     1.2000     2.0000
        1    N (  N3)        5   CA (  CT)       12    C (   C)       14    N (   N)     0.0000     2.0000     0.0000     1.2000     2.0000
        4   H3 (   H)        1    N (  N3)        5   CA (  CT)        6   HA (  HP)     0.1556     3.0000     0.0000     1.2000     2.0000
        4   H3 (   H)        1    N (  N3)        5   CA (  CT)        7   CB (  CT)     0.1556     3.0000     0.0000     1.2000     2.0000
        4   H3 (   H)        1    N (  N3)        5   CA (  CT)       12    C (   C)     0.1556     3.0000     0.0000     1.2000     2.0000
        3   H2 (   H)        1    N (  N3)        5   CA (  CT)        6   HA (  HP)     0.1556     3.0000     0.0000     1.2000     2.0000
        3   H2 (   H)        1    N (  N3)        5   CA (  CT)        7   CB (  CT)     0.1556     3.0000     0.0000     1.2000     2.0000
        3   H2 (   H)        1    N (  N3)        5   CA (  CT)       12    C (   C)     0.1556     3.0000     0.0000     1.2000     2.0000
        2   H1 (   H)        1    N (  N3)        5   CA (  CT)        6   HA (  HP)     0.1556     3.0000     0.0000     1.2000     2.0000
        2   H1 (   H)        1    N (  N3)        5   CA (  CT)        7   CB (  CT)     0.1556     3.0000     0.0000     1.2000     2.0000
        2   H1 (   H)        1    N (  N3)        5   CA (  CT)       12    C (   C)     0.1556     3.0000     0.0000     1.2000     2.0000
        1    N (  N3)        5   CA (  CT)        7   CB (  CT)        8  HB2 (  H1)     0.1556     3.0000     0.0000     1.2000     2.0000
        1    N (  N3)        5   CA (  CT)        7   CB (  CT)        9  HB3 (  H1)     0.1556     3.0000     0.0000     1.2000     2.0000
"""

SET_BOND = """\
             Atom 1              Atom 2       R eq   Frc Cnst
    288   CA (  CT)     290   CB (  CT)     1.5000   300.0000
    288   CA (  CT)     294    C (   C)     1.5220   317.0000
    286    N (   N)     288   CA (  CT)     1.4490   337.0000
    317   CA (  CT)     319   CB (  CT)     1.5000   300.0000
    317   CA (  CT)     323    C (   C)     1.5220   317.0000
    315    N (   N)     317   CA (  CT)     1.4490   337.0000
    438   CA (  CT)     440   CB (  CT)     1.5000   300.0000
    438   CA (  CT)     444    C (   C)     1.5220   317.0000
    436    N (   N)     438   CA (  CT)     1.4490   337.0000
    586   CA (  CT)     588   CB (  CT)     1.5000   300.0000
    586   CA (  CT)     592    C (   C)     1.5220   317.0000
    584    N (   N)     586   CA (  CT)     1.4490   337.0000
    694   CA (  CT)     696   CB (  CT)     1.5000   300.0000
    694   CA (  CT)     700    C (   C)     1.5220   317.0000
    692    N (   N)     694   CA (  CT)     1.4490   337.0000
    847   CA (  CT)     849   CB (  CT)     1.5000   300.0000
    847   CA (  CT)     853    C (   C)     1.5220   317.0000
    845    N (   N)     847   CA (  CT)     1.4490   337.0000
   1009   CA (  CT)    1011   CB (  CT)     1.5000   300.0000
   1009   CA (  CT)    1015    C (   C)     1.5220   317.0000
   1007    N (   N)    1009   CA (  CT)     1.4490   337.0000
   1331   CA (  CT)    1333   CB (  CT)     1.5000   300.0000
   1331   CA (  CT)    1337    C (   C)     1.5220   317.0000
   1329    N (   N)    1331   CA (  CT)     1.4490   337.0000
   1341   CA (  CT)    1343   CB (  CT)     1.5000   300.0000
   1341   CA (  CT)    1347    C (   C)     1.5220   317.0000
   1339    N (   N)    1341   CA (  CT)     1.4490   337.0000
   1410   CA (  CT)    1412   CB (  CT)     1.5000   300.0000
   1410   CA (  CT)    1416    C (   C)     1.5220   317.0000
   1408    N (   N)    1410   CA (  CT)     1.4490   337.0000
   1603   CA (  CT)    1605   CB (  CT)     1.5000   300.0000
   1603   CA (  CT)    1609    C (   C)     1.5220   317.0000
   1601    N (   N)    1603   CA (  CT)     1.4490   337.0000
   1646   CA (  CT)    1648   CB (  CT)     1.5000   300.0000
   1646   CA (  CT)    1652    C (   C)     1.5220   317.0000
   1644    N (   N)    1646   CA (  CT)     1.4490   337.0000
    288   CA (  CT)     289   HA (  H1)     1.0900   340.0000
    317   CA (  CT)     318   HA (  H1)     1.0900   340.0000
    438   CA (  CT)     439   HA (  H1)     1.0900   340.0000
    586   CA (  CT)     587   HA (  H1)     1.0900   340.0000
    694   CA (  CT)     695   HA (  H1)     1.0900   340.0000
    847   CA (  CT)     848   HA (  H1)     1.0900   340.0000
   1009   CA (  CT)    1010   HA (  H1)     1.0900   340.0000
   1331   CA (  CT)    1332   HA (  H1)     1.0900   340.0000
   1341   CA (  CT)    1342   HA (  H1)     1.0900   340.0000
   1410   CA (  CT)    1411   HA (  H1)     1.0900   340.0000
   1603   CA (  CT)    1604   HA (  H1)     1.0900   340.0000
   1646   CA (  CT)    1647   HA (  H1)     1.0900   340.0000
"""

SET_ANGLE = """\
             Atom 1               Atom 2               Atom 3   Frc Cnst   Theta eq
    290   CB (  CT)      288   CA (  CT)      294    C (   C)    63.0000   111.1000
    286    N (   N)      288   CA (  CT)      290   CB (  CT)    80.0000   109.7000
    319   CB (  CT)      317   CA (  CT)      323    C (   C)    63.0000   111.1000
    315    N (   N)      317   CA (  CT)      319   CB (  CT)    80.0000   109.7000
    440   CB (  CT)      438   CA (  CT)      444    C (   C)    63.0000   111.1000
    436    N (   N)      438   CA (  CT)      440   CB (  CT)    80.0000   109.7000
    588   CB (  CT)      586   CA (  CT)      592    C (   C)    63.0000   111.1000
    584    N (   N)      586   CA (  CT)      588   CB (  CT)    80.0000   109.7000
    696   CB (  CT)      694   CA (  CT)      700    C (   C)    63.0000   111.1000
    692    N (   N)      694   CA (  CT)      696   CB (  CT)    80.0000   109.7000
    849   CB (  CT)      847   CA (  CT)      853    C (   C)    63.0000   111.1000
    845    N (   N)      847   CA (  CT)      849   CB (  CT)    80.0000   109.7000
   1011   CB (  CT)     1009   CA (  CT)     1015    C (   C)    63.0000   111.1000
   1007    N (   N)     1009   CA (  CT)     1011   CB (  CT)    80.0000   109.7000
   1333   CB (  CT)     1331   CA (  CT)     1337    C (   C)    63.0000   111.1000
   1329    N (   N)     1331   CA (  CT)     1333   CB (  CT)    80.0000   109.7000
   1343   CB (  CT)     1341   CA (  CT)     1347    C (   C)    63.0000   111.1000
   1339    N (   N)     1341   CA (  CT)     1343   CB (  CT)    80.0000   109.7000
   1412   CB (  CT)     1410   CA (  CT)     1416    C (   C)    63.0000   111.1000
   1408    N (   N)     1410   CA (  CT)     1412   CB (  CT)    80.0000   109.7000
   1605   CB (  CT)     1603   CA (  CT)     1609    C (   C)    63.0000   111.1000
   1601    N (   N)     1603   CA (  CT)     1605   CB (  CT)    80.0000   109.7000
   1648   CB (  CT)     1646   CA (  CT)     1652    C (   C)    63.0000   111.1000
   1644    N (   N)     1646   CA (  CT)     1648   CB (  CT)    80.0000   109.7000
    292  HB2 (  HC)      290   CB (  CT)      293  HB3 (  HC)    35.0000   109.5000
    291  HB1 (  HC)      290   CB (  CT)      292  HB2 (  HC)    35.0000   109.5000
    291  HB1 (  HC)      290   CB (  CT)      293  HB3 (  HC)    35.0000   109.5000
    289   HA (  H1)      288   CA (  CT)      290   CB (  CT)    50.0000   109.5000
    288   CA (  CT)      290   CB (  CT)      291  HB1 (  HC)    40.0000   100.0000
    288   CA (  CT)      290   CB (  CT)      292  HB2 (  HC)    50.0000   109.5000
    288   CA (  CT)      290   CB (  CT)      293  HB3 (  HC)    50.0000   109.5000
    321  HB2 (  HC)      319   CB (  CT)      322  HB3 (  HC)    35.0000   109.5000
    320  HB1 (  HC)      319   CB (  CT)      321  HB2 (  HC)    35.0000   109.5000
    320  HB1 (  HC)      319   CB (  CT)      322  HB3 (  HC)    35.0000   109.5000
    318   HA (  H1)      317   CA (  CT)      319   CB (  CT)    50.0000   109.5000
    317   CA (  CT)      319   CB (  CT)      320  HB1 (  HC)    40.0000   100.0000
    317   CA (  CT)      319   CB (  CT)      321  HB2 (  HC)    50.0000   109.5000
    317   CA (  CT)      319   CB (  CT)      322  HB3 (  HC)    50.0000   109.5000
    442  HB2 (  HC)      440   CB (  CT)      443  HB3 (  HC)    35.0000   109.5000
    441  HB1 (  HC)      440   CB (  CT)      442  HB2 (  HC)    35.0000   109.5000
    441  HB1 (  HC)      440   CB (  CT)      443  HB3 (  HC)    35.0000   109.5000
    439   HA (  H1)      438   CA (  CT)      440   CB (  CT)    50.0000   109.5000
    438   CA (  CT)      440   CB (  CT)      441  HB1 (  HC)    40.0000   100.0000
    438   CA (  CT)      440   CB (  CT)      442  HB2 (  HC)    50.0000   109.5000
    438   CA (  CT)      440   CB (  CT)      443  HB3 (  HC)    50.0000   109.5000
    590  HB2 (  HC)      588   CB (  CT)      591  HB3 (  HC)    35.0000   109.5000
    589  HB1 (  HC)      588   CB (  CT)      590  HB2 (  HC)    35.0000   109.5000
    589  HB1 (  HC)      588   CB (  CT)      591  HB3 (  HC)    35.0000   109.5000
    587   HA (  H1)      586   CA (  CT)      588   CB (  CT)    50.0000   109.5000
    586   CA (  CT)      588   CB (  CT)      589  HB1 (  HC)    40.0000   100.0000
    586   CA (  CT)      588   CB (  CT)      590  HB2 (  HC)    50.0000   109.5000
    586   CA (  CT)      588   CB (  CT)      591  HB3 (  HC)    50.0000   109.5000
    698  HB2 (  HC)      696   CB (  CT)      699  HB3 (  HC)    35.0000   109.5000
    697  HB1 (  HC)      696   CB (  CT)      698  HB2 (  HC)    35.0000   109.5000
    697  HB1 (  HC)      696   CB (  CT)      699  HB3 (  HC)    35.0000   109.5000
    695   HA (  H1)      694   CA (  CT)      696   CB (  CT)    50.0000   109.5000
    694   CA (  CT)      696   CB (  CT)      697  HB1 (  HC)    40.0000   100.0000
    694   CA (  CT)      696   CB (  CT)      698  HB2 (  HC)    50.0000   109.5000
    694   CA (  CT)      696   CB (  CT)      699  HB3 (  HC)    50.0000   109.5000
    851  HB2 (  HC)      849   CB (  CT)      852  HB3 (  HC)    35.0000   109.5000
    850  HB1 (  HC)      849   CB (  CT)      851  HB2 (  HC)    35.0000   109.5000
    850  HB1 (  HC)      849   CB (  CT)      852  HB3 (  HC)    35.0000   109.5000
    848   HA (  H1)      847   CA (  CT)      849   CB (  CT)    50.0000   109.5000
    847   CA (  CT)      849   CB (  CT)      850  HB1 (  HC)    40.0000   100.0000
    847   CA (  CT)      849   CB (  CT)      851  HB2 (  HC)    50.0000   109.5000
    847   CA (  CT)      849   CB (  CT)      852  HB3 (  HC)    50.0000   109.5000
   1013  HB2 (  HC)     1011   CB (  CT)     1014  HB3 (  HC)    35.0000   109.5000
   1012  HB1 (  HC)     1011   CB (  CT)     1013  HB2 (  HC)    35.0000   109.5000
   1012  HB1 (  HC)     1011   CB (  CT)     1014  HB3 (  HC)    35.0000   109.5000
   1010   HA (  H1)     1009   CA (  CT)     1011   CB (  CT)    50.0000   109.5000
   1009   CA (  CT)     1011   CB (  CT)     1012  HB1 (  HC)    40.0000   100.0000
   1009   CA (  CT)     1011   CB (  CT)     1013  HB2 (  HC)    50.0000   109.5000
   1009   CA (  CT)     1011   CB (  CT)     1014  HB3 (  HC)    50.0000   109.5000
   1335  HB2 (  HC)     1333   CB (  CT)     1336  HB3 (  HC)    35.0000   109.5000
   1334  HB1 (  HC)     1333   CB (  CT)     1335  HB2 (  HC)    35.0000   109.5000
   1334  HB1 (  HC)     1333   CB (  CT)     1336  HB3 (  HC)    35.0000   109.5000
   1332   HA (  H1)     1331   CA (  CT)     1333   CB (  CT)    50.0000   109.5000
   1331   CA (  CT)     1333   CB (  CT)     1334  HB1 (  HC)    40.0000   100.0000
   1331   CA (  CT)     1333   CB (  CT)     1335  HB2 (  HC)    50.0000   109.5000
   1331   CA (  CT)     1333   CB (  CT)     1336  HB3 (  HC)    50.0000   109.5000
   1345  HB2 (  HC)     1343   CB (  CT)     1346  HB3 (  HC)    35.0000   109.5000
   1344  HB1 (  HC)     1343   CB (  CT)     1345  HB2 (  HC)    35.0000   109.5000
   1344  HB1 (  HC)     1343   CB (  CT)     1346  HB3 (  HC)    35.0000   109.5000
   1342   HA (  H1)     1341   CA (  CT)     1343   CB (  CT)    50.0000   109.5000
   1341   CA (  CT)     1343   CB (  CT)     1344  HB1 (  HC)    40.0000   100.0000
   1341   CA (  CT)     1343   CB (  CT)     1345  HB2 (  HC)    50.0000   109.5000
   1341   CA (  CT)     1343   CB (  CT)     1346  HB3 (  HC)    50.0000   109.5000
   1414  HB2 (  HC)     1412   CB (  CT)     1415  HB3 (  HC)    35.0000   109.5000
   1413  HB1 (  HC)     1412   CB (  CT)     1414  HB2 (  HC)    35.0000   109.5000
   1413  HB1 (  HC)     1412   CB (  CT)     1415  HB3 (  HC)    35.0000   109.5000
   1411   HA (  H1)     1410   CA (  CT)     1412   CB (  CT)    50.0000   109.5000
   1410   CA (  CT)     1412   CB (  CT)     1413  HB1 (  HC)    40.0000   100.0000
   1410   CA (  CT)     1412   CB (  CT)     1414  HB2 (  HC)    50.0000   109.5000
   1410   CA (  CT)     1412   CB (  CT)     1415  HB3 (  HC)    50.0000   109.5000
   1607  HB2 (  HC)     1605   CB (  CT)     1608  HB3 (  HC)    35.0000   109.5000
   1606  HB1 (  HC)     1605   CB (  CT)     1607  HB2 (  HC)    35.0000   109.5000
   1606  HB1 (  HC)     1605   CB (  CT)     1608  HB3 (  HC)    35.0000   109.5000
   1604   HA (  H1)     1603   CA (  CT)     1605   CB (  CT)    50.0000   109.5000
   1603   CA (  CT)     1605   CB (  CT)     1606  HB1 (  HC)    40.0000   100.0000
   1603   CA (  CT)     1605   CB (  CT)     1607  HB2 (  HC)    50.0000   109.5000
   1603   CA (  CT)     1605   CB (  CT)     1608  HB3 (  HC)    50.0000   109.5000
   1650  HB2 (  HC)     1648   CB (  CT)     1651  HB3 (  HC)    35.0000   109.5000
   1649  HB1 (  HC)     1648   CB (  CT)     1650  HB2 (  HC)    35.0000   109.5000
   1649  HB1 (  HC)     1648   CB (  CT)     1651  HB3 (  HC)    35.0000   109.5000
   1647   HA (  H1)     1646   CA (  CT)     1648   CB (  CT)    50.0000   109.5000
   1646   CA (  CT)     1648   CB (  CT)     1649  HB1 (  HC)    40.0000   100.0000
   1646   CA (  CT)     1648   CB (  CT)     1650  HB2 (  HC)    50.0000   109.5000
   1646   CA (  CT)     1648   CB (  CT)     1651  HB3 (  HC)    50.0000   109.5000
"""

PRINT_LJMATRIX = """
                  Atom Type 1                   Atom Type 2   A coefficient   B coefficient      R i,j    Eps i,j
------------------------------------------------------------------------------------------------------------------
            N,N2,N3,NA,NB [1]             N,N2,N3,NA,NB [1]   944293.233000      801.323529   3.648000   0.170000
            N,N2,N3,NA,NB [1]                         H [2]     2126.011810       20.960420   2.424000   0.051662
            N,N2,N3,NA,NB [1]                        CT [3]   995480.466000      736.907417   3.732000   0.136374
            N,N2,N3,NA,NB [1]                        HP [4]    20179.142500       64.575606   2.924000   0.051662
            N,N2,N3,NA,NB [1]                        H1 [5]    62066.599700      113.252061   3.211000   0.051662
            N,N2,N3,NA,NB [1]                        OH [6]   744975.864000      750.714425   3.545000   0.189124
            N,N2,N3,NA,NB [1]                        HO [7]        0.000000        0.000000   0.000000   0.000000
            N,N2,N3,NA,NB [1] C,C*,CA,CB,CC,CN,CR,CV,CW [8]   882619.071000      653.361429   3.732000   0.120913
            N,N2,N3,NA,NB [1]                      O,O2 [9]   606829.342000      677.220874   3.485200   0.188944
            N,N2,N3,NA,NB [1]                       HC [10]    89677.698900      136.131731   3.311000   0.051662
            N,N2,N3,NA,NB [1]                       H5 [11]    54614.725300      105.031585   3.183000   0.050498
            N,N2,N3,NA,NB [1]                       H4 [12]    65847.387000      115.327881   3.233000   0.050498
            N,N2,N3,NA,NB [1]                       HA [13]    79162.715400      126.451907   3.283000   0.050498
            N,N2,N3,NA,NB [1]                        S [14]  2015621.900000     1289.234040   3.824000   0.206155
"""

SUMMARY = """\
Amino Acid Residues:   3
Nucleic Acid Residues: 0
Number of cations:     10
Number of anions:      10
Num. of solvent mols:  696
Num. of unknown res:   6
Total charge (e-):     -0.0000
Total mass (amu):      13618.1680
Number of atoms:       2174
Number of residues:    725
Residue set:           ACE, ALA, CYX, Cl-, NME, Na+, WAT
Residue count:         ACE: 3, ALA: 1, CYX: 2, Cl-: 10, NME: 3, Na+: 10, WAT: 696
System volume (ang^3): 26461.48
System density (g/mL): 0.854596
"""

PRINT_DETAILSC = """
The mask @1 matches 1 atoms:

   ATOM    RES  RESNAME  NAME  TYPE   At.#   LJ Radius    LJ Depth      Mass    Charge GB Radius GB Screen
      1      1      ALA     N   NH3      7      1.8500      0.2000   14.0070   -0.3000    1.5500    0.7900
"""

PRINT_BONDSC = """\
             Atom 1              Atom 2       R eq   Frc Cnst
      1    N ( NH3)       5   CA ( CT1)     1.4800   200.0000
      2  HT1 (  HC)       1    N ( NH3)     1.0400   403.0000
      3  HT2 (  HC)       1    N ( NH3)     1.0400   403.0000
      4  HT3 (  HC)       1    N ( NH3)     1.0400   403.0000
"""

PRINT_ANGLESC = """\
             Atom 1               Atom 2               Atom 3   Frc Cnst   Theta eq
      1    N ( NH3)        5   CA ( CT1)        7   CB ( CT3)    67.7000   110.0000
      1    N ( NH3)        5   CA ( CT1)       11    C (   C)    43.7000   110.0000
      2  HT1 (  HC)        1    N ( NH3)        3  HT2 (  HC)    44.0000   109.5000
      2  HT1 (  HC)        1    N ( NH3)        4  HT3 (  HC)    44.0000   109.5000
      2  HT1 (  HC)        1    N ( NH3)        5   CA ( CT1)    30.0000   109.5000
      3  HT2 (  HC)        1    N ( NH3)        4  HT3 (  HC)    44.0000   109.5000
      3  HT2 (  HC)        1    N ( NH3)        5   CA ( CT1)    30.0000   109.5000
      4  HT3 (  HC)        1    N ( NH3)        5   CA ( CT1)    30.0000   109.5000
      1    N ( NH3)        5   CA ( CT1)        6   HA (  HB)    51.5000   107.5000
"""

PRINT_DIHEDRALSC = """\
               Atom 1               Atom 2               Atom 3               Atom 4     Height  Periodic.      Phase  EEL Scale  VDW Scale
        1    N ( NH3)        5   CA ( CT1)       11    C (   C)       12    O (   O)     0.0000     1.0000     0.0000     1.0000     1.0000
        1    N ( NH3)        5   CA ( CT1)       11    C (   C)       13    N ( NH1)     0.6000     1.0000     0.0000     1.0000     1.0000
        1    N ( NH3)        5   CA ( CT1)        7   CB ( CT3)        8  HB1 (  HA)     0.2000     3.0000     0.0000     1.0000     1.0000
        1    N ( NH3)        5   CA ( CT1)        7   CB ( CT3)        9  HB2 (  HA)     0.2000     3.0000     0.0000     1.0000     1.0000
        1    N ( NH3)        5   CA ( CT1)        7   CB ( CT3)       10  HB3 (  HA)     0.2000     3.0000     0.0000     1.0000     1.0000
        2  HT1 (  HC)        1    N ( NH3)        5   CA ( CT1)        6   HA (  HB)     0.1000     3.0000     0.0000     1.0000     1.0000
        2  HT1 (  HC)        1    N ( NH3)        5   CA ( CT1)        7   CB ( CT3)     0.1000     3.0000     0.0000     1.0000     1.0000
        2  HT1 (  HC)        1    N ( NH3)        5   CA ( CT1)       11    C (   C)     0.1000     3.0000     0.0000     1.0000     1.0000
        3  HT2 (  HC)        1    N ( NH3)        5   CA ( CT1)        6   HA (  HB)     0.1000     3.0000     0.0000     1.0000     1.0000
        3  HT2 (  HC)        1    N ( NH3)        5   CA ( CT1)        7   CB ( CT3)     0.1000     3.0000     0.0000     1.0000     1.0000
        3  HT2 (  HC)        1    N ( NH3)        5   CA ( CT1)       11    C (   C)     0.1000     3.0000     0.0000     1.0000     1.0000
        4  HT3 (  HC)        1    N ( NH3)        5   CA ( CT1)        6   HA (  HB)     0.1000     3.0000     0.0000     1.0000     1.0000
        4  HT3 (  HC)        1    N ( NH3)        5   CA ( CT1)        7   CB ( CT3)     0.1000     3.0000     0.0000     1.0000     1.0000
        4  HT3 (  HC)        1    N ( NH3)        5   CA ( CT1)       11    C (   C)     0.1000     3.0000     0.0000     1.0000     1.0000
"""

SET_BONDC = """\
             Atom 1              Atom 2       R eq   Frc Cnst
      7   CB ( CT3)       5   CA ( CT1)     1.5000   300.0000
      1    N ( NH3)       5   CA ( CT1)     1.4800   200.0000
     11    C (   C)       5   CA ( CT1)     1.4900   250.0000
     17   CB ( CT3)      15   CA ( CT1)     1.5000   300.0000
     13    N ( NH1)      15   CA ( CT1)     1.4300   320.0000
     21    C (   C)      15   CA ( CT1)     1.4900   250.0000
     27   CB ( CT3)      25   CA ( CT1)     1.5000   300.0000
     23    N ( NH1)      25   CA ( CT1)     1.4300   320.0000
     31    C (  CC)      25   CA ( CT1)     1.5220   200.0000
      5   CA ( CT1)       6   HA (  HB)     1.0800   330.0000
     15   CA ( CT1)      16   HA (  HB)     1.0800   330.0000
     25   CA ( CT1)      26   HA (  HB)     1.0800   330.0000
"""

SET_ANGLEC = """\
             Atom 1               Atom 2               Atom 3   Frc Cnst   Theta eq
      1    N ( NH3)        5   CA ( CT1)        7   CB ( CT3)    67.7000   110.0000
      7   CB ( CT3)        5   CA ( CT1)       11    C (   C)    52.0000   108.0000
     13    N ( NH1)       15   CA ( CT1)       17   CB ( CT3)    70.0000   113.5000
     17   CB ( CT3)       15   CA ( CT1)       21    C (   C)    52.0000   108.0000
     23    N ( NH1)       25   CA ( CT1)       27   CB ( CT3)    70.0000   113.5000
     27   CB ( CT3)       25   CA ( CT1)       31    C (  CC)    52.0000   108.0000
      6   HA (  HB)        5   CA ( CT1)        7   CB ( CT3)    35.0000   111.0000
      5   CA ( CT1)        7   CB ( CT3)        8  HB1 (  HA)    40.0000   100.0000
      5   CA ( CT1)        7   CB ( CT3)        9  HB2 (  HA)    33.4300   110.1000
      5   CA ( CT1)        7   CB ( CT3)       10  HB3 (  HA)    33.4300   110.1000
      8  HB1 (  HA)        7   CB ( CT3)        9  HB2 (  HA)    35.5000   108.4000
      8  HB1 (  HA)        7   CB ( CT3)       10  HB3 (  HA)    35.5000   108.4000
      9  HB2 (  HA)        7   CB ( CT3)       10  HB3 (  HA)    35.5000   108.4000
     16   HA (  HB)       15   CA ( CT1)       17   CB ( CT3)    35.0000   111.0000
     15   CA ( CT1)       17   CB ( CT3)       18  HB1 (  HA)    40.0000   100.0000
     15   CA ( CT1)       17   CB ( CT3)       19  HB2 (  HA)    33.4300   110.1000
     15   CA ( CT1)       17   CB ( CT3)       20  HB3 (  HA)    33.4300   110.1000
     18  HB1 (  HA)       17   CB ( CT3)       19  HB2 (  HA)    35.5000   108.4000
     18  HB1 (  HA)       17   CB ( CT3)       20  HB3 (  HA)    35.5000   108.4000
     19  HB2 (  HA)       17   CB ( CT3)       20  HB3 (  HA)    35.5000   108.4000
     26   HA (  HB)       25   CA ( CT1)       27   CB ( CT3)    35.0000   111.0000
     25   CA ( CT1)       27   CB ( CT3)       28  HB1 (  HA)    40.0000   100.0000
     25   CA ( CT1)       27   CB ( CT3)       29  HB2 (  HA)    33.4300   110.1000
     25   CA ( CT1)       27   CB ( CT3)       30  HB3 (  HA)    33.4300   110.1000
     28  HB1 (  HA)       27   CB ( CT3)       29  HB2 (  HA)    35.5000   108.4000
     28  HB1 (  HA)       27   CB ( CT3)       30  HB3 (  HA)    35.5000   108.4000
     29  HB2 (  HA)       27   CB ( CT3)       30  HB3 (  HA)    35.5000   108.4000
"""

PRINT_LJMATRIXC = """
Atom Type 1 Atom Type 2   A coefficient   B coefficient      R i,j    Eps i,j
------------------------------------------------------------------------------
  NH3 [1]   NH3 [1]  1316590.401168     1026.290564   3.700000   0.200000
  NH3 [1]  H,HC [2]      609.333680       15.289896   2.074500   0.095917
  NH3 [1]   CT1 [3]  1535031.988974      623.165940   4.125000   0.063246
  NH3 [1] HA,HB [4]    68302.639169      134.620719   3.170000   0.066332
  NH3 [1]   CT3 [5]  1615031.807669      903.962743   3.910000   0.126491
  NH3 [1]     C [6]  1573041.189142      966.063587   3.850000   0.148324
  NH3 [1]     O [7]   620648.710655      620.162833   3.550000   0.154919
  NH3 [1]   NH1 [8]  1316590.401168     1026.290564   3.700000   0.200000
  NH3 [1]    CC [9]  1254852.764965      770.652143   3.850000   0.118322
  NH3 [1]   OC [10]   620648.710655      620.162833   3.550000   0.154919
"""

SUMMARYC1 = """\
Amino Acid Residues:   154
Nucleic Acid Residues: 0
Number of cations:     0
Number of anions:      0
Num. of solvent mols:  17856
Num. of unknown res:   5
Total charge (e-):     -11.0000
Total mass (amu):      339672.2974
Number of atoms:       56057
Number of residues:    18015
Residue set:           ALA, ARG, ASN, ASP, CYS, GLN, GLU
                       GLY, HSD, ILE, LEU, LYS, MET, PHE
                       PRO, SER, THR, TIP3, TRP, TYR, VAL
Residue count:         ALA: 13, ARG: 9, ASN: 5, ASP: 14, CYS: 2, GLN: 4, GLU: 12
                       GLY: 10, HSD: 5, ILE: 12, LEU: 11, LYS: 6, MET: 5, PHE: 6
                       PRO: 10, SER: 9, THR: 6, TIP3: 17856, TRP: 5, TYR: 4, VAL: 11
System volume (ang^3): 615109.34
System density (g/mL): 0.916989
"""

PRINT_DETAILSA = """
The mask :1-2 matches 12 atoms:

   ATOM    RES  RESNAME  NAME  TYPE   At.#      Mass
      1      1      ACE   CH3     C      6   12.0110
      2      1      ACE     C     C      6   12.0110
      3      1      ACE     O     O      8   15.9990
      4      1      ACE  H31H     H      1    1.0080
      5      1      ACE  H32H     H      1    1.0080
      6      1      ACE  H33H     H      1    1.0080
      7      2      NME     N     N      7   14.0070
      8      2      NME   CH3     C      6   12.0110
      9      2      NME     H    HN      1    1.0080
     10      2      NME  H31H     H      1    1.0080
     11      2      NME  H32H     H      1    1.0080
     12      2      NME  H33H     H      1    1.0080
"""

PRINT_BONDSA = """\
             Atom 1              Atom 2       R eq   Frc Cnst
      1  CH3 (   C)       2    C (   C)     1.5090   345.0000
      1  CH3 (   C)       4 H31H (   H)     1.1120   341.0000
      1  CH3 (   C)       5 H32H (   H)     1.1120   341.0000
      1  CH3 (   C)       6 H33H (   H)     1.1120   341.0000
"""

PRINT_ANGLESA = """\
             Atom 1               Atom 2               Atom 3   Frc Cnst   Theta eq
      2    C (   C)        1  CH3 (   C)        4 H31H (   H)    39.0000   109.5000
      2    C (   C)        1  CH3 (   C)        5 H32H (   H)    39.0000   109.5000
      2    C (   C)        1  CH3 (   C)        6 H33H (   H)    39.0000   109.5000
      4 H31H (   H)        1  CH3 (   C)        5 H32H (   H)    40.0000   107.8000
      4 H31H (   H)        1  CH3 (   C)        6 H33H (   H)    40.0000   107.8000
      5 H32H (   H)        1  CH3 (   C)        6 H33H (   H)    40.0000   107.8000
"""

PRINT_DIHEDRALSA = """\
               Atom 1               Atom 2               Atom 3               Atom 4     Height  Periodic.      Phase  EEL Scale  VDW Scale
        4 H31H (   H)        1  CH3 (   C)        2    C (   C)        3    O (   O)     0.2350     3.0000     0.0000        N/A        N/A
        4 H31H (   H)        1  CH3 (   C)        2    C (   C)        7    N (   N)    -0.0100     3.0000     0.0000        N/A        N/A
        5 H32H (   H)        1  CH3 (   C)        2    C (   C)        3    O (   O)     0.2350     3.0000     0.0000        N/A        N/A
        5 H32H (   H)        1  CH3 (   C)        2    C (   C)        7    N (   N)    -0.0100     3.0000     0.0000        N/A        N/A
        6 H33H (   H)        1  CH3 (   C)        2    C (   C)        3    O (   O)     0.2350     3.0000     0.0000        N/A        N/A
        6 H33H (   H)        1  CH3 (   C)        2    C (   C)        7    N (   N)    -0.0100     3.0000     0.0000        N/A        N/A
        1  CH3 (   C)        2    C (   C)        7    N (   N)        8  CH3 (   C)    -1.0000     1.0000     0.0000        N/A        N/A
        1  CH3 (   C)        2    C (   C)        7    N (   N)        8  CH3 (   C)     2.0000     2.0000   180.0000        N/A        N/A
        1  CH3 (   C)        2    C (   C)        7    N (   N)        8  CH3 (   C)     2.0000     3.0000     0.0000        N/A        N/A
        1  CH3 (   C)        2    C (   C)        7    N (   N)        9    H (  HN)     1.0000     2.0000   180.0000        N/A        N/A
        1  CH3 (   C)        2    C (   C)        7    N (   N)        9    H (  HN)     0.8000     3.0000     0.0000        N/A        N/A
"""

SUMMARYA1 = """\
Amino Acid Residues:   0
Nucleic Acid Residues: 0
Number of cations:     0
Number of anions:      0
Num. of solvent mols:  818
Num. of unknown res:   2
Total charge (e-):     0.0000
Total mass (amu):      14809.3650
Number of atoms:       2466
Number of residues:    820
Residue set:           ACE, NME, WAT
Residue count:         ACE: 1, NME: 1, WAT: 818
System volume (ang^3): 1.00
System density (g/mL): 24591.940605
"""

SUMMARYA2 = """\
Amino Acid Residues:   0
Nucleic Acid Residues: 0
Number of cations:     0
Number of anions:      0
Num. of solvent mols:  818
Num. of unknown res:   2
Total charge (e-):     0.0000
Total mass (amu):      14809.3650
Number of atoms:       2466
Number of residues:    820
Residue set:           ACE, NME, WAT
Residue count:         ACE: 1, NME: 1, WAT: 818
System volume (ang^3): 27031.52
System density (g/mL): 0.909751
"""

PRINT_BONDS_2MASKS = """\
             Atom 1              Atom 2       R eq   Frc Cnst
      1    N (  N3)       3   H2 (   H)     1.0100   434.0000
"""

PRINT_ANGLES_2MASKS_1 = """\
             Atom 1               Atom 2               Atom 3   Frc Cnst   Theta eq
      1    N (  N3)        5   CA (  CT)        7   CB (  CT)    80.0000   111.2000
      1    N (  N3)        5   CA (  CT)       12    C (   C)    80.0000   111.2000
      1    N (  N3)        5   CA (  CT)        6   HA (  HP)    50.0000   109.5000
"""

PRINT_ANGLES_2MASKS_2 = """\
             Atom 1               Atom 2               Atom 3   Frc Cnst   Theta eq
      4   H3 (   H)        1    N (  N3)        5   CA (  CT)    50.0000   109.5000
      3   H2 (   H)        1    N (  N3)        4   H3 (   H)    35.0000   109.5000
      3   H2 (   H)        1    N (  N3)        5   CA (  CT)    50.0000   109.5000
      2   H1 (   H)        1    N (  N3)        3   H2 (   H)    35.0000   109.5000
      2   H1 (   H)        1    N (  N3)        4   H3 (   H)    35.0000   109.5000
      2   H1 (   H)        1    N (  N3)        5   CA (  CT)    50.0000   109.5000
"""

PRINT_ANGLES_2MASKS_3 = """\
             Atom 1               Atom 2               Atom 3   Frc Cnst   Theta eq
      1    N (  N3)        5   CA (  CT)        7   CB (  CT)    80.0000   111.2000
      1    N (  N3)        5   CA (  CT)       12    C (   C)    80.0000   111.2000
      1    N (  N3)        5   CA (  CT)        6   HA (  HP)    50.0000   109.5000
"""

PRINT_ANGLES_3MASKS = """\
             Atom 1               Atom 2               Atom 3   Frc Cnst   Theta eq
      1    N (  N3)        5   CA (  CT)        7   CB (  CT)    80.0000   111.2000
"""

PRINT_DIHEDRALS_2MASKS = """\
               Atom 1               Atom 2               Atom 3               Atom 4     Height  Periodic.      Phase  EEL Scale  VDW Scale
        1    N (  N3)        5   CA (  CT)        7   CB (  CT)       10   OG (  OH)     0.1556     3.0000     0.0000     1.2000     2.0000
        1    N (  N3)        5   CA (  CT)       12    C (   C)       13    O (   O)     0.0000     2.0000     0.0000     1.2000     2.0000
        1    N (  N3)        5   CA (  CT)       12    C (   C)       14    N (   N)     0.0000     2.0000     0.0000     1.2000     2.0000
        1    N (  N3)        5   CA (  CT)        7   CB (  CT)        8  HB2 (  H1)     0.1556     3.0000     0.0000     1.2000     2.0000
        1    N (  N3)        5   CA (  CT)        7   CB (  CT)        9  HB3 (  H1)     0.1556     3.0000     0.0000     1.2000     2.0000
"""

PRINT_DIHEDRALS_3MASKS = """\
               Atom 1               Atom 2               Atom 3               Atom 4     Height  Periodic.      Phase  EEL Scale  VDW Scale
        4   H3 (   H)        1    N (  N3)        5   CA (  CT)        6   HA (  HP)     0.1556     3.0000     0.0000     1.2000     2.0000
        4   H3 (   H)        1    N (  N3)        5   CA (  CT)        7   CB (  CT)     0.1556     3.0000     0.0000     1.2000     2.0000
        4   H3 (   H)        1    N (  N3)        5   CA (  CT)       12    C (   C)     0.1556     3.0000     0.0000     1.2000     2.0000
"""

PRINT_DIHEDRALS_4MASKS = """\
               Atom 1               Atom 2               Atom 3               Atom 4     Height  Periodic.      Phase  EEL Scale  VDW Scale
        7   CB (  CT)        5   CA (  CT)       12    C (   C)       14    N (   N)     0.0700     2.0000     0.0000     1.2000     2.0000
M       7   CB (  CT)        5   CA (  CT)       12    C (   C)       14    N (   N)     0.1000     4.0000     0.0000     1.2000     2.0000
"""

PRINTLJMATRIX_MGNACL = """
Atom Type 1 Atom Type 2   A coefficient   B coefficient   C coefficient      R i,j    Eps i,j
----------------------------------------------------------------------------------------------
   HC [1]  Mg2+ [8]     7354.219030       23.534357       35.617936   2.924000   0.018828
CT,CX [2]  Mg2+ [8]    97526.006400      139.243106       97.650208   3.345000   0.049701
    C [3]  Mg2+ [8]    86469.113200      123.456587      124.432687   3.345000   0.044066
    O [4]  Mg2+ [8]    53861.861300      121.802066       52.368490   3.098200   0.068860
    N [5]  Mg2+ [8]    89596.771500      149.010737      100.319252   3.261000   0.061956
    H [6]  Mg2+ [8]       96.094933        2.690198       35.617936   2.037000   0.018828
   H1 [7]  Mg2+ [8]     4843.780440       19.099689       35.617936   2.824000   0.018828
 Mg2+ [8]  Mg2+ [8]     7170.637930       25.448794        4.417729   2.874000   0.022580
 Mg2+ [8]   Na+ [9]     9783.484750       32.222975       22.962985   2.910000   0.026532
 Mg2+ [8]  Cl- [10]   492351.740000      462.292889      296.473338   3.587000   0.108517
 Mg2+ [8]   OW [11]    68897.634400      127.063907      132.900000   3.205300   0.058584
 Mg2+ [8]   HW [12]        0.000000        0.000000        0.000000   0.000000   0.000000
"""

PRINTLJMATRIX_NACLMG = """
Atom Type 1 Atom Type 2   A coefficient   B coefficient   C coefficient      R i,j    Eps i,j
----------------------------------------------------------------------------------------------
   HC [1] Mg2+ [10]     7354.219030       23.534357       35.617936   2.924000   0.018828
CT,CX [2] Mg2+ [10]    97526.006400      139.243106       97.650208   3.345000   0.049701
    C [3] Mg2+ [10]    86469.113200      123.456587      124.432687   3.345000   0.044066
    O [4] Mg2+ [10]    53861.861300      121.802066       52.368490   3.098200   0.068860
    N [5] Mg2+ [10]    89596.771500      149.010737      100.319252   3.261000   0.061956
    H [6] Mg2+ [10]       96.094933        2.690198       35.617936   2.037000   0.018828
   H1 [7] Mg2+ [10]     4843.780440       19.099689       35.617936   2.824000   0.018828
  Na+ [8] Mg2+ [10]     9783.484750       32.222975       22.962985   2.910000   0.026532
  Cl- [9] Mg2+ [10]   492351.740000      462.292889      296.473338   3.587000   0.108517
Mg2+ [10] Mg2+ [10]     7170.637930       25.448794        4.417729   2.874000   0.022580
Mg2+ [10]   OW [11]    68897.634400      127.063907      132.900000   3.205300   0.058584
Mg2+ [10]   HW [12]        0.000000        0.000000        0.000000   0.000000   0.000000
"""

PDB_SUMMARY = """\
Amino Acid Residues:   129
Nucleic Acid Residues: 0
Number of cations:     0
Number of anions:      0
Num. of solvent mols:  139
Num. of unknown res:   6
Total charge (e-):     0.0000
Total mass (amu):      15942.3372
Number of atoms:       1164
Number of residues:    274
Residue set:           ALA, ARG, ASN, ASP, CYS, GLN, GLU
                       GLY, HIS, HOH, ILE, LEU, LYS, MET
                       NO3, PHE, PRO, SER, THR, TRP, TYR
                       VAL
Residue count:         ALA: 12, ARG: 11, ASN: 14, ASP: 7, CYS: 8, GLN: 3, GLU: 2
                       GLY: 12, HIS: 1, HOH: 139, ILE: 6, LEU: 8, LYS: 6, MET: 2
                       NO3: 6, PHE: 3, PRO: 2, SER: 10, THR: 7, TRP: 6, TYR: 3
                       VAL: 6
System volume (ang^3): 25998.98
System density (g/mL): 1.018244
"""

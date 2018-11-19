
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
System volume (ang^3): 27031.52
System density (g/mL): 0.909751
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

AMOEBA_SMALL_MDIN = """\
Input file for AMOEBA simulations.
 &cntrl
     ! Add whatever variables you need here
     ntb=1, ntt=1, ntp=0, ! PBC, thermostat, barostat
     irest=0, ntx=1,      ! restart flags
 /
 &amoeba
     ! Some basic potential parameters. For better
     ! energy conservation you need to adjust these
     ! defaults
     beeman_integrator=1,   ! Use Beeman integrator
     dipole_scf_tol=0.01,   ! 10e-6 gives good NVE

     ! You should not generally modify these variables:
     do_valence=1, do_bond=1, do_ureyb=1,
     do_reg_angle=1, do_trig_angle=1, do_opbend=1,
     do_torsion=1, do_pi_torsion=1, do_strbend=1,
     do_torsion_torsion=0,
 /
"""

AS4_TITR_OUTPUT = """\
AS4 	pKa =   4.0
    ATOM     STATE 0     STATE 1     STATE 2     STATE 3     STATE 4
       N     -0.4157     -0.4157     -0.4157     -0.4157     -0.4157
       H      0.2719      0.2719      0.2719      0.2719      0.2719
      CA      0.0341      0.0341      0.0341      0.0341      0.0341
      HA      0.0864      0.0864      0.0864      0.0864      0.0864
      CB     -0.1783     -0.0316     -0.0316     -0.0316     -0.0316
     HB2     -0.0122      0.0488      0.0488      0.0488      0.0488
     HB3     -0.0122      0.0488      0.0488      0.0488      0.0488
      CG      0.7994      0.6462      0.6462      0.6462      0.6462
     OD1     -0.8014     -0.5554     -0.5554     -0.6376     -0.6376
     OD2     -0.8014     -0.6376     -0.6376     -0.5554     -0.5554
    HD21      0.0000      0.4747      0.0000      0.0000      0.0000
       C      0.5973      0.5973      0.5973      0.5973      0.5973
       O     -0.5679     -0.5679     -0.5679     -0.5679     -0.5679
    HD22      0.0000      0.0000      0.4747      0.0000      0.0000
    HD11      0.0000      0.0000      0.0000      0.4747      0.0000
    HD12      0.0000      0.0000      0.0000      0.0000      0.4747
--------------------------------------------------------------------
Prot Cnt           0           1           1           1           1
--------------------------------------------------------------------
pKa Corr      0.0000      4.0000      4.0000      4.0000      4.0000
--------------------------------------------------------------------
Reference Energies (ES = Explicit solvent, IS = Implicit solvent)

igb=1 IS     0.00000    21.42980    21.42980    21.42980    21.42980
igb=2 IS     0.00000    26.88946    26.88946    26.88946    26.88946
igb=5 IS     0.00000    26.59805    26.59805    26.59805    26.59805
igb=7 IS     0.00000    23.41811    23.41811    23.41811    23.41811
igb=8 IS     0.00000    26.34489    26.34489    26.34489    26.34489
igb=1 ES     0.00000     Not Set     Not Set     Not Set     Not Set
igb=2 ES     0.00000    33.26130    33.26130    33.26130    33.26130
igb=5 ES     0.00000    32.06435    32.06435    32.06435    32.06435
igb=7 ES     0.00000    28.26235    28.26235    28.26235    28.26235
igb=8 ES     0.00000    31.28604    31.28604    31.28604    31.28604
--------------------------------------------------------------------
Reference Energies for Internal Dielectric of 2.0

igb=1 IS     0.00000     Not Set     Not Set     Not Set     Not Set
igb=2 IS     0.00000    12.67691    12.67691    12.67691    12.67691
igb=5 IS     0.00000    13.08491    13.08491    13.08491    13.08491
igb=7 IS     0.00000     Not Set     Not Set     Not Set     Not Set
igb=8 IS     0.00000     Not Set     Not Set     Not Set     Not Set
igb=1 ES     0.00000     Not Set     Not Set     Not Set     Not Set
igb=2 ES     0.00000     Not Set     Not Set     Not Set     Not Set
igb=5 ES     0.00000     Not Set     Not Set     Not Set     Not Set
igb=7 ES     0.00000     Not Set     Not Set     Not Set     Not Set
igb=8 ES     0.00000     Not Set     Not Set     Not Set     Not Set
"""

PRINT_BONDS_MEASURE = """\
             Atom 1              Atom 2       R eq   Frc Cnst   Distance     Energy
      5    C (   C)       6    O (   O)     1.2290   570.0000     1.2290     0.0000
      5    C (   C)       7    N (   N)     1.3350   490.0000     1.3350     0.0000
      2  CH3 (  CT)       5    C (   C)     1.5220   317.0000     1.5300     0.0203
     18    C (   C)      19    O (   O)     1.2290   570.0000     1.2290     0.0000
     18    C (   C)      20    N (   N)     1.3350   490.0000     1.3350     0.0000
     14   CG (   C)      15  OD1 (   O)     1.2290   570.0000     1.2600     0.5478
     14   CG (   C)      16  OD2 (  OH)     1.3640   450.0000     1.2600     4.8672
     11   CB (  2C)      14   CG (   C)     1.5220   317.0000     1.5270     0.0079
      9   CA (  CX)      11   CB (  2C)     1.5260   310.0000     1.5250     0.0003
      9   CA (  CX)      18    C (   C)     1.5220   317.0000     1.5220     0.0000
      7    N (   N)       9   CA (  CX)     1.4490   337.0000     1.4490     0.0000
     20    N (   N)      22  CH3 (  CT)     1.4490   337.0000     1.4490     0.0000
      2  CH3 (  CT)       3 HH32 (  HC)     1.0900   340.0000     1.0900     0.0000
      2  CH3 (  CT)       4 HH33 (  HC)     1.0900   340.0000     1.0900     0.0000
      1 HH31 (  HC)       2  CH3 (  CT)     1.0900   340.0000     1.0900     0.0000
     16  OD2 (  OH)      17  HD2 (  HO)     0.9600   553.0000     0.9600     0.0000
     11   CB (  2C)      12  HB2 (  HC)     1.0900   340.0000     1.0900     0.0000
     11   CB (  2C)      13  HB3 (  HC)     1.0900   340.0000     1.0900     0.0000
      9   CA (  CX)      10   HA (  H1)     1.0900   340.0000     1.0900     0.0000
      7    N (   N)       8    H (   H)     1.0100   434.0000     1.0100     0.0000
     22  CH3 (  CT)      23 HH31 (  H1)     1.0900   340.0000     1.0900     0.0000
     22  CH3 (  CT)      24 HH32 (  H1)     1.0900   340.0000     1.0900     0.0000
     22  CH3 (  CT)      25 HH33 (  H1)     1.0900   340.0000     1.0900     0.0000
     20    N (   N)      21    H (   H)     1.0100   434.0000     1.0100     0.0000
"""

PRINT_ANGLES_MEASURE = """\
             Atom 1               Atom 2               Atom 3   Frc Cnst   Theta eq      Angle     Energy
      6    O (   O)        5    C (   C)        7    N (   N)    80.0000   122.9001   122.8999     0.0000
      5    C (   C)        7    N (   N)        9   CA (  CX)    50.0000   121.9001   121.9001     0.0000
      2  CH3 (  CT)        5    C (   C)        6    O (   O)    80.0000   120.4001   120.5001     0.0002
      2  CH3 (  CT)        5    C (   C)        7    N (   N)    70.0000   116.6000   116.6001     0.0000
     19    O (   O)       18    C (   C)       20    N (   N)    80.0000   122.9001   122.8999     0.0000
     18    C (   C)       20    N (   N)       22  CH3 (  CT)    50.0000   121.9001   121.9001     0.0000
     15  OD1 (   O)       14   CG (   C)       16  OD2 (  OH)    80.0000   120.0001   125.5999     0.7642
     11   CB (  2C)        9   CA (  CX)       18    C (   C)    63.0000   111.1000   107.7945     0.2097
     11   CB (  2C)       14   CG (   C)       15  OD1 (   O)    80.0000   120.4001   117.2001     0.2495
     11   CB (  2C)       14   CG (   C)       16  OD2 (  OH)    80.0000   110.0000   117.2000     1.2633
      9   CA (  CX)       11   CB (  2C)       14   CG (   C)    63.0000   111.1000   109.4700     0.0510
      9   CA (  CX)       18    C (   C)       19    O (   O)    80.0000   120.4001   120.5000     0.0002
      9   CA (  CX)       18    C (   C)       20    N (   N)    70.0000   116.6000   116.6001     0.0000
      7    N (   N)        9   CA (  CX)       11   CB (  2C)    80.0000   109.7000   111.1000     0.0478
      7    N (   N)        9   CA (  CX)       18    C (   C)    63.0000   110.1000   111.1001     0.0192
      5    C (   C)        7    N (   N)        8    H (   H)    50.0000   120.0001   120.0000     0.0000
      4 HH33 (  HC)        2  CH3 (  CT)        5    C (   C)    50.0000   109.5000   108.6358     0.0114
      3 HH32 (  HC)        2  CH3 (  CT)        4 HH33 (  HC)    35.0000   109.5000   109.4423     0.0000
      3 HH32 (  HC)        2  CH3 (  CT)        5    C (   C)    50.0000   109.5000   108.6358     0.0114
      1 HH31 (  HC)        2  CH3 (  CT)        3 HH32 (  HC)    35.0000   109.5000   109.5001     0.0000
      1 HH31 (  HC)        2  CH3 (  CT)        4 HH33 (  HC)    35.0000   109.5000   109.5001     0.0000
      1 HH31 (  HC)        2  CH3 (  CT)        5    C (   C)    50.0000   109.5000   111.1000     0.0390
     18    C (   C)       20    N (   N)       21    H (   H)    50.0000   120.0001   120.0001     0.0000
     14   CG (   C)       16  OD2 (  OH)       17  HD2 (  HO)    50.0000   113.0000   109.5000     0.1866
     13  HB3 (  HC)       11   CB (  2C)       14   CG (   C)    50.0000   109.5000   109.4574     0.0000
     12  HB2 (  HC)       11   CB (  2C)       13  HB3 (  HC)    35.0000   109.5000   109.4423     0.0000
     12  HB2 (  HC)       11   CB (  2C)       14   CG (   C)    50.0000   109.5000   109.4575     0.0000
     10   HA (  H1)        9   CA (  CX)       11   CB (  2C)    50.0000   109.5000   108.6356     0.0114
     10   HA (  H1)        9   CA (  CX)       18    C (   C)    50.0000   109.5000   108.6358     0.0114
      9   CA (  CX)       11   CB (  2C)       12  HB2 (  HC)    50.0000   109.5000   109.5001     0.0000
      9   CA (  CX)       11   CB (  2C)       13  HB3 (  HC)    50.0000   109.5000   109.5001     0.0000
      8    H (   H)        7    N (   N)        9   CA (  CX)    50.0000   118.0401   118.0999     0.0001
      7    N (   N)        9   CA (  CX)       10   HA (  H1)    50.0000   109.5000   109.5000     0.0000
     24 HH32 (  H1)       22  CH3 (  CT)       25 HH33 (  H1)    35.0000   109.5000   109.4425     0.0000
     23 HH31 (  H1)       22  CH3 (  CT)       24 HH32 (  H1)    35.0000   109.5000   109.4424     0.0000
     23 HH31 (  H1)       22  CH3 (  CT)       25 HH33 (  H1)    35.0000   109.5000   109.4423     0.0000
     21    H (   H)       20    N (   N)       22  CH3 (  CT)    50.0000   118.0401   118.0999     0.0001
     20    N (   N)       22  CH3 (  CT)       23 HH31 (  H1)    50.0000   109.5000   109.5000     0.0000
     20    N (   N)       22  CH3 (  CT)       24 HH32 (  H1)    50.0000   109.5000   109.5001     0.0000
     20    N (   N)       22  CH3 (  CT)       25 HH33 (  H1)    50.0000   109.5000   109.5000     0.0000
"""

PRINT_DIHEDRALS_MEASURE = """\
               Atom 1               Atom 2               Atom 3               Atom 4     Height  Periodic.      Phase  EEL Scale  VDW Scale   Dihedral     Energy
        6    O (   O)        5    C (   C)        7    N (   N)        9   CA (  CX)     2.5000     2.0000   180.0001     1.2000     2.0000    -0.0001     0.0000
        5    C (   C)        7    N (   N)        9   CA (  CX)       11   CB (  2C)     2.0000     1.0000     0.0000     1.2000     2.0000    60.0002     3.0000
M       5    C (   C)        7    N (   N)        9   CA (  CX)       11   CB (  2C)     1.8000     2.0000     0.0000     1.2000     2.0000    60.0002     0.9000
M       5    C (   C)        7    N (   N)        9   CA (  CX)       11   CB (  2C)     0.8000     3.0000     0.0000     1.2000     2.0000    60.0002     0.0000
M       5    C (   C)        7    N (   N)        9   CA (  CX)       11   CB (  2C)     0.0000     4.0000     0.0000     1.2000     2.0000    60.0002     0.0000
        5    C (   C)        7    N (   N)        9   CA (  CX)       18    C (   C)     0.0000     1.0000     0.0000     1.2000     2.0000  -179.9998     0.0000
M       5    C (   C)        7    N (   N)        9   CA (  CX)       18    C (   C)     0.2700     2.0000     0.0000     1.2000     2.0000  -179.9998     0.5400
M       5    C (   C)        7    N (   N)        9   CA (  CX)       18    C (   C)     0.4200     3.0000     0.0000     1.2000     2.0000  -179.9998     0.0000
M       5    C (   C)        7    N (   N)        9   CA (  CX)       18    C (   C)     0.0000     4.0000     0.0000     1.2000     2.0000  -179.9998     0.0000
        2  CH3 (  CT)        5    C (   C)        7    N (   N)        9   CA (  CX)     2.5000     2.0000   180.0001     1.2000     2.0000   179.9999     0.0000
       19    O (   O)       18    C (   C)       20    N (   N)       22  CH3 (  CT)     2.5000     2.0000   180.0001     1.2000     2.0000    -0.0002     0.0000
       14   CG (   C)       11   CB (  2C)        9   CA (  CX)       18    C (   C)     0.6680     1.0000     0.0000     1.2000     2.0000    58.0551     1.0214
M      14   CG (   C)       11   CB (  2C)        9   CA (  CX)       18    C (   C)     0.2050     2.0000   180.0001     1.2000     2.0000    58.0551     0.2952
M      14   CG (   C)       11   CB (  2C)        9   CA (  CX)       18    C (   C)     1.0280     3.0000   180.0001     1.2000     2.0000    58.0551     2.0507
M      14   CG (   C)       11   CB (  2C)        9   CA (  CX)       18    C (   C)     0.0030     4.0000     0.0000     1.2000     2.0000    58.0551     0.0012
       11   CB (  2C)        9   CA (  CX)       18    C (   C)       19    O (   O)     0.0000     2.0000     0.0000     1.2000     2.0000   121.9450     0.0000
       11   CB (  2C)        9   CA (  CX)       18    C (   C)       20    N (   N)     0.2000     1.0000     0.0000     1.2000     2.0000   -58.0551     0.3058
M      11   CB (  2C)        9   CA (  CX)       18    C (   C)       20    N (   N)     0.2000     2.0000     0.0000     1.2000     2.0000   -58.0551     0.1120
M      11   CB (  2C)        9   CA (  CX)       18    C (   C)       20    N (   N)     0.4000     3.0000     0.0000     1.2000     2.0000   -58.0551     0.0021
M      11   CB (  2C)        9   CA (  CX)       18    C (   C)       20    N (   N)     0.0000     4.0000     0.0000     1.2000     2.0000   -58.0551     0.0000
        9   CA (  CX)       11   CB (  2C)       14   CG (   C)       15  OD1 (   O)     1.2830     1.0000     0.0000     1.2000     2.0000    90.0001     1.2830
M       9   CA (  CX)       11   CB (  2C)       14   CG (   C)       15  OD1 (   O)     0.4310     2.0000     0.0000     1.2000     2.0000    90.0001     0.0000
M       9   CA (  CX)       11   CB (  2C)       14   CG (   C)       15  OD1 (   O)     0.0570     3.0000   180.0001     1.2000     2.0000    90.0001     0.0570
M       9   CA (  CX)       11   CB (  2C)       14   CG (   C)       15  OD1 (   O)     0.1180     4.0000     0.0000     1.2000     2.0000    90.0001     0.2360
        9   CA (  CX)       11   CB (  2C)       14   CG (   C)       16  OD2 (  OH)     0.1190     1.0000     0.0000     1.2000     2.0000   -89.9999     0.1190
M       9   CA (  CX)       11   CB (  2C)       14   CG (   C)       16  OD2 (  OH)     1.2550     2.0000   180.0001     1.2000     2.0000   -89.9999     2.5100
M       9   CA (  CX)       11   CB (  2C)       14   CG (   C)       16  OD2 (  OH)     0.1730     3.0000   180.0001     1.2000     2.0000   -89.9999     0.1730
M       9   CA (  CX)       11   CB (  2C)       14   CG (   C)       16  OD2 (  OH)     0.1570     4.0000   180.0001     1.2000     2.0000   -89.9999     0.0000
        9   CA (  CX)       18    C (   C)       20    N (   N)       22  CH3 (  CT)     2.5000     2.0000   180.0001     1.2000     2.0000   179.9999     0.0000
        7    N (   N)        9   CA (  CX)       11   CB (  2C)       14   CG (   C)     0.7020     1.0000   180.0001     1.2000     2.0000  -179.9999     1.4040
M       7    N (   N)        9   CA (  CX)       11   CB (  2C)       14   CG (   C)     0.1990     2.0000   180.0001     1.2000     2.0000  -179.9999     0.0000
M       7    N (   N)        9   CA (  CX)       11   CB (  2C)       14   CG (   C)     0.8560     3.0000     0.0000     1.2000     2.0000  -179.9999     0.0000
M       7    N (   N)        9   CA (  CX)       11   CB (  2C)       14   CG (   C)     0.0080     4.0000     0.0000     1.2000     2.0000  -179.9999     0.0160
        7    N (   N)        9   CA (  CX)       18    C (   C)       19    O (   O)     0.0000     2.0000     0.0000     1.2000     2.0000    -0.0000     0.0000
        7    N (   N)        9   CA (  CX)       18    C (   C)       20    N (   N)     0.4500     1.0000   180.0001     1.2000     2.0000   179.9999     0.9000
M       7    N (   N)        9   CA (  CX)       18    C (   C)       20    N (   N)     1.5800     2.0000   180.0001     1.2000     2.0000   179.9999     0.0000
M       7    N (   N)        9   CA (  CX)       18    C (   C)       20    N (   N)     0.5500     3.0000   180.0001     1.2000     2.0000   179.9999     1.1000
M       7    N (   N)        9   CA (  CX)       18    C (   C)       20    N (   N)     0.0000     4.0000     0.0000     1.2000     2.0000   179.9999     0.0000
I       6    O (   O)        2  CH3 (  CT)        5    C (   C)        7    N (   N)    10.5000     2.0000   180.0001     0.0000     0.0000   179.9999     0.0000
I      19    O (   O)        9   CA (  CX)       18    C (   C)       20    N (   N)    10.5000     2.0000   180.0001     0.0000     0.0000   179.9999     0.0000
I      16  OD2 (  OH)       11   CB (  2C)       14   CG (   C)       15  OD1 (   O)    10.5000     2.0000   180.0001     0.0000     0.0000   180.0000     0.0000
        6    O (   O)        5    C (   C)        7    N (   N)        8    H (   H)     2.0000     1.0000     0.0000     1.2000     2.0000   179.9999     0.0000
M       6    O (   O)        5    C (   C)        7    N (   N)        8    H (   H)     2.5000     2.0000   180.0001     1.2000     2.0000   179.9999     0.0000
        5    C (   C)        7    N (   N)        9   CA (  CX)       10   HA (  H1)     0.0000     2.0000     0.0000     1.2000     2.0000   -59.9997     0.0000
        4 HH33 (  HC)        2  CH3 (  CT)        5    C (   C)        6    O (   O)     0.8000     1.0000     0.0000     1.2000     2.0000  -120.5122     0.3938
M       4 HH33 (  HC)        2  CH3 (  CT)        5    C (   C)        6    O (   O)     0.0000     2.0000     0.0000     1.2000     2.0000  -120.5122     0.0000
M       4 HH33 (  HC)        2  CH3 (  CT)        5    C (   C)        6    O (   O)     0.0800     3.0000   180.0001     1.2000     2.0000  -120.5122     0.0000
        4 HH33 (  HC)        2  CH3 (  CT)        5    C (   C)        7    N (   N)     0.0000     2.0000     0.0000     1.2000     2.0000    59.4877     0.0000
        3 HH32 (  HC)        2  CH3 (  CT)        5    C (   C)        6    O (   O)     0.8000     1.0000     0.0000     1.2000     2.0000   120.5122     0.3938
M       3 HH32 (  HC)        2  CH3 (  CT)        5    C (   C)        6    O (   O)     0.0000     2.0000     0.0000     1.2000     2.0000   120.5122     0.0000
M       3 HH32 (  HC)        2  CH3 (  CT)        5    C (   C)        6    O (   O)     0.0800     3.0000   180.0001     1.2000     2.0000   120.5122     0.0000
        3 HH32 (  HC)        2  CH3 (  CT)        5    C (   C)        7    N (   N)     0.0000     2.0000     0.0000     1.2000     2.0000   -59.4879     0.0000
        2  CH3 (  CT)        5    C (   C)        7    N (   N)        8    H (   H)     2.5000     2.0000   180.0001     1.2000     2.0000    -0.0000     0.0000
        1 HH31 (  HC)        2  CH3 (  CT)        5    C (   C)        6    O (   O)     0.8000     1.0000     0.0000     1.2000     2.0000     0.0000     1.6000
M       1 HH31 (  HC)        2  CH3 (  CT)        5    C (   C)        6    O (   O)     0.0000     2.0000     0.0000     1.2000     2.0000     0.0000     0.0000
M       1 HH31 (  HC)        2  CH3 (  CT)        5    C (   C)        6    O (   O)     0.0800     3.0000   180.0001     1.2000     2.0000     0.0000     0.0000
        1 HH31 (  HC)        2  CH3 (  CT)        5    C (   C)        7    N (   N)     0.0000     2.0000     0.0000     1.2000     2.0000   179.9999     0.0000
       19    O (   O)       18    C (   C)       20    N (   N)       21    H (   H)     2.0000     1.0000     0.0000     1.2000     2.0000   179.9999     0.0000
M      19    O (   O)       18    C (   C)       20    N (   N)       21    H (   H)     2.5000     2.0000   180.0001     1.2000     2.0000   179.9999     0.0000
       18    C (   C)       20    N (   N)       22  CH3 (  CT)       23 HH31 (  H1)     0.0000     2.0000     0.0000     1.2000     2.0000     0.0001     0.0000
       18    C (   C)       20    N (   N)       22  CH3 (  CT)       24 HH32 (  H1)     0.0000     2.0000     0.0000     1.2000     2.0000   120.0002     0.0000
       18    C (   C)       20    N (   N)       22  CH3 (  CT)       25 HH33 (  H1)     0.0000     2.0000     0.0000     1.2000     2.0000  -119.9997     0.0000
       15  OD1 (   O)       14   CG (   C)       16  OD2 (  OH)       17  HD2 (  HO)     1.9000     1.0000     0.0000     1.2000     2.0000     0.0001     3.8000
M      15  OD1 (   O)       14   CG (   C)       16  OD2 (  OH)       17  HD2 (  HO)     2.3000     2.0000   180.0001     1.2000     2.0000     0.0001     0.0000
       13  HB3 (  HC)       11   CB (  2C)        9   CA (  CX)       18    C (   C)     0.1556     3.0000     0.0000     1.2000     2.0000   -61.9449     0.0008
       13  HB3 (  HC)       11   CB (  2C)       14   CG (   C)       15  OD1 (   O)     0.8000     1.0000     0.0000     1.2000     2.0000  -149.9738     0.1074
M      13  HB3 (  HC)       11   CB (  2C)       14   CG (   C)       15  OD1 (   O)     0.0000     2.0000     0.0000     1.2000     2.0000  -149.9738     0.0000
M      13  HB3 (  HC)       11   CB (  2C)       14   CG (   C)       15  OD1 (   O)     0.0800     3.0000   180.0001     1.2000     2.0000  -149.9738     0.0799
       13  HB3 (  HC)       11   CB (  2C)       14   CG (   C)       16  OD2 (  OH)     0.0000     2.0000     0.0000     1.2000     2.0000    30.0262     0.0000
       12  HB2 (  HC)       11   CB (  2C)        9   CA (  CX)       18    C (   C)     0.1556     3.0000     0.0000     1.2000     2.0000   178.0552     0.0008
       12  HB2 (  HC)       11   CB (  2C)       14   CG (   C)       15  OD1 (   O)     0.8000     1.0000     0.0000     1.2000     2.0000   -30.0261     1.4926
M      12  HB2 (  HC)       11   CB (  2C)       14   CG (   C)       15  OD1 (   O)     0.0000     2.0000     0.0000     1.2000     2.0000   -30.0261     0.0000
M      12  HB2 (  HC)       11   CB (  2C)       14   CG (   C)       15  OD1 (   O)     0.0800     3.0000   180.0001     1.2000     2.0000   -30.0261     0.0801
       12  HB2 (  HC)       11   CB (  2C)       14   CG (   C)       16  OD2 (  OH)     0.0000     2.0000     0.0000     1.2000     2.0000   149.9739     0.0000
       11   CB (  2C)       14   CG (   C)       16  OD2 (  OH)       17  HD2 (  HO)     2.3000     2.0000   180.0001     1.2000     2.0000  -179.9999     0.0000
       10   HA (  H1)        9   CA (  CX)       11   CB (  2C)       12  HB2 (  HC)     0.1556     3.0000     0.0000     1.2000     2.0000    60.5122     0.0001
       10   HA (  H1)        9   CA (  CX)       11   CB (  2C)       13  HB3 (  HC)     0.1556     3.0000     0.0000     1.2000     2.0000  -179.4879     0.0001
       10   HA (  H1)        9   CA (  CX)       11   CB (  2C)       14   CG (   C)     0.1556     3.0000     0.0000     1.2000     2.0000   -59.4879     0.0001
       10   HA (  H1)        9   CA (  CX)       18    C (   C)       19    O (   O)     0.8000     1.0000     0.0000     1.2000     2.0000  -120.5122     0.3938
M      10   HA (  H1)        9   CA (  CX)       18    C (   C)       19    O (   O)     0.0000     2.0000     0.0000     1.2000     2.0000  -120.5122     0.0000
M      10   HA (  H1)        9   CA (  CX)       18    C (   C)       19    O (   O)     0.0800     3.0000   180.0001     1.2000     2.0000  -120.5122     0.0000
       10   HA (  H1)        9   CA (  CX)       18    C (   C)       20    N (   N)     0.0000     2.0000     0.0000     1.2000     2.0000    59.4877     0.0000
        9   CA (  CX)       18    C (   C)       20    N (   N)       21    H (   H)     2.5000     2.0000   180.0001     1.2000     2.0000     0.0000     0.0000
        8    H (   H)        7    N (   N)        9   CA (  CX)       10   HA (  H1)     0.0000     2.0000     0.0000     1.2000     2.0000   120.0002     0.0000
        8    H (   H)        7    N (   N)        9   CA (  CX)       11   CB (  2C)     0.0000     2.0000     0.0000     1.2000     2.0000  -119.9999     0.0000
        8    H (   H)        7    N (   N)        9   CA (  CX)       18    C (   C)     0.0000     2.0000     0.0000     1.2000     2.0000     0.0002     0.0000
        7    N (   N)        9   CA (  CX)       11   CB (  2C)       12  HB2 (  HC)     0.1556     3.0000     0.0000     1.2000     2.0000   -59.9998     0.0000
        7    N (   N)        9   CA (  CX)       11   CB (  2C)       13  HB3 (  HC)     0.1556     3.0000     0.0000     1.2000     2.0000    60.0001     0.0000
       21    H (   H)       20    N (   N)       22  CH3 (  CT)       23 HH31 (  H1)     0.0000     2.0000     0.0000     1.2000     2.0000  -179.9999     0.0000
       21    H (   H)       20    N (   N)       22  CH3 (  CT)       24 HH32 (  H1)     0.0000     2.0000     0.0000     1.2000     2.0000   -59.9999     0.0000
       21    H (   H)       20    N (   N)       22  CH3 (  CT)       25 HH33 (  H1)     0.0000     2.0000     0.0000     1.2000     2.0000    60.0002     0.0000
I       8    H (   H)        5    C (   C)        7    N (   N)        9   CA (  CX)     1.1000     2.0000   180.0001     0.0000     0.0000   179.9999     0.0000
I      21    H (   H)       18    C (   C)       20    N (   N)       22  CH3 (  CT)     1.1000     2.0000   180.0001     0.0000     0.0000   179.9999     0.0000
"""

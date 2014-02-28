""" Creates a leaprc file for use in amberParm-created OFF and frcmod files """

from sys import stderr

def MakeLeaprc(frcmod, off, leaprc, radius_set = 'modified Bondi radii (mbondi)', 
               logfile = 'leap.log'):
   """ Writes out the leaprc file from the given parm, OFF files """

   file = open(leaprc,'w',0)

   if radius_set.strip() == "Bondi radii (bondi)": radii = "bondi"
   elif radius_set.strip() == "amber6 modified Bondi radii (amber6)": radii = "amber6"
   elif radius_set.strip() == "modified Bondi radii (mbondi)": radii = 'mbondi'
   elif radius_set.strip() == "H(N)-modified Bondi radii (mbondi2)": radii = 'mbondi2'
   else:
      print >> stderr, 'Error: Unknown GB Radius set! [ %s ]' % radius_set
      print >> stderr, '       GB radii will not be reliable unless they are default mbondi'
      radii = 'mbondi'


   file.write("""logFile %s
#
# ----- leaprc for loading the prmtop-based force field
# ----- NOTE: this is designed for PDB format 3!
#
#  Load the main parameter set.
#
parm = loadamberparams %s

#
#  Load main library file
#
loadOFF %s

#
#  This needs to be kept to properly recognize termini
#
addPdbResMap {
  { 0 "ALA" "NALA" } { 1 "ALA" "CALA" }
  { 0 "ARG" "NARG" } { 1 "ARG" "CARG" }
  { 0 "ASN" "NASN" } { 1 "ASN" "CASN" }
  { 0 "ASP" "NASP" } { 1 "ASP" "CASP" }
  { 0 "CYS" "NCYS" } { 1 "CYS" "CCYS" }
  { 0 "CYX" "NCYX" } { 1 "CYX" "CCYX" }
  { 0 "GLN" "NGLN" } { 1 "GLN" "CGLN" }
  { 0 "GLU" "NGLU" } { 1 "GLU" "CGLU" }
  { 0 "GLY" "NGLY" } { 1 "GLY" "CGLY" }
  { 0 "HID" "NHID" } { 1 "HID" "CHID" }
  { 0 "HIE" "NHIE" } { 1 "HIE" "CHIE" }
  { 0 "HIP" "NHIP" } { 1 "HIP" "CHIP" }
  { 0 "ILE" "NILE" } { 1 "ILE" "CILE" }
  { 0 "LEU" "NLEU" } { 1 "LEU" "CLEU" }
  { 0 "LYS" "NLYS" } { 1 "LYS" "CLYS" }
  { 0 "MET" "NMET" } { 1 "MET" "CMET" }
  { 0 "PHE" "NPHE" } { 1 "PHE" "CPHE" }
  { 0 "PRO" "NPRO" } { 1 "PRO" "CPRO" }
  { 0 "SER" "NSER" } { 1 "SER" "CSER" }
  { 0 "THR" "NTHR" } { 1 "THR" "CTHR" }
  { 0 "TRP" "NTRP" } { 1 "TRP" "CTRP" }
  { 0 "TYR" "NTYR" } { 1 "TYR" "CTYR" }
  { 0 "VAL" "NVAL" } { 1 "VAL" "CVAL" }
  { 0 "HIS" "NHIS" } { 1 "HIS" "CHIS" }
  { 0 "G" "G5"  } { 1 "G" "G3"  } 
  { 0 "A" "A5"  } { 1 "A" "A3"  } 
  { 0 "C" "C5"  } { 1 "C" "C3"  } 
  { 0 "U" "U5"  } { 1 "U" "U3"  } 
  { 0 "DG" "DG5"  } { 1 "DG" "DG3"  }  
  { 0 "DA" "DA5"  } { 1 "DA" "DA3"  }  
  { 0 "DC" "DC5"  } { 1 "DC" "DC3"  }  
  { 0 "DT" "DT5"  } { 1 "DT" "DT3"  }  
#  some old Amber residue names for RNA:
  { 0  "RA" "RA5" } { 1 "RA" "RA3"} 
  { 0  "RC" "RC5" } { 1 "RC" "RC3"} 
  { 0  "RG" "RG5" } { 1 "RG" "RG3"} 
  { 0  "RU" "RU5" } { 1 "RU" "RU3"} 
}

# Set the proper radius set
set default PBRadii %s
""" % (logfile, frcmod, off, radii))

   file.close()

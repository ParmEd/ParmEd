"""
Holds an Amber residue class for easy writing of mol2 and
OFF files
"""

from chemistry import periodic_table

AMINO = [ "ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH", "GLN", "GLU", "GLY", 
          "HID", "HIE", "HIP", "ILE", "LEU", "LYN", "LYS", "MET", "PHE", "PRO", "SER", "THR", 
          "TRP", "TYR", "VAL" ]

NUCLEIC = [ "RA", "A", "DA", "RC", "C", "DC", "RG", "G", "DG", "U", "RU", "DT" ]

def ToResidue(molecule, resnum):
   """ Returns a residue object from a specified resnum in a molecule object """

   begin = molecule.residue_pointers[resnum]
   try:
      end = molecule.residue_pointers[resnum+1] # 1 after the end, as it needs to be for array slices
   except IndexError:
      end = len(molecule.atoms) # catch an error if we pick the last residue

   # Gather the bonding information for the given residue. Create a new bond array that is specific
   # to this residue (and whose partners are indexed as though this residue starts from atom 0). Also
   # find the "connect" atoms that bond to adjacent residues.
   bonds = []
   connects = [0,0]
   head = tail = 0
   for i in range(begin, end):
      bonds.append([])
      for j in molecule.bonds[i]:
         if j < begin and head == 0:
            head = i - begin + 1
            connects[0] = head
         elif j < begin:
            connects.append(i - begin + 1)
         elif j >= end:
            if molecule.residue_container[j] == resnum + 1:
               # we keep track of index in connects and the atom number so we can
               # remove it from here and replace it in connects[1] at the end.
               # we do this so that we only count the last "tail" residue as connect1
               tail = len(connects)
               connects.append(i - begin + 1)
            else:
               connects.append(i - begin + 1)
         else: # otherwise it's an intra-residue bond. adjust j for new numbering and add it to bonds
            bonds[i-begin].append(j - begin + 1)

   # Now find the tail (if there is one) and put it in the 2nd location
   if tail > 1:
      connects[1] = connects.pop(tail)
   tail = connects[1]

   return Residue(atoms = molecule.atoms[begin:end], atom_types = molecule.atom_types[begin:end], 
                  charges = molecule.charges[begin:end], bonds = bonds, coords = molecule.coords[3*begin:3*end], 
                  head = head, tail = tail, connects = connects, resname = molecule.residues[resnum],
                  elements = molecule.elements[begin:end])

class Residue:
   """ This is a defined residue in a biomolecular system """

   def __init__(self, atoms=[], atom_types=[], charges=[], bonds=[], coords=[],
                head=0, tail=0, connects=[0,0], resname = '', elements=[]):
      """ initializes the residue """

      global NUCLEIC, AMINO

      self.atoms = atoms
      self.atom_types = atom_types
      self.charges = charges
      self.bonds = bonds
      self.coords = coords
      self.head = head
      self.tail = tail
      self.natom = len(self.atoms)
      self.is_nterm = head == 0
      self.is_cterm = tail == 0
      self.name = self.orig_name = resname
      self.elements = elements
      self.connects = connects

      # Rename termini to distinguish from non-termini
      if self.is_nterm and not self.is_cterm and self.name in AMINO:
         self.name = "N" + self.name
      elif self.is_cterm and not self.is_nterm and self.name in AMINO:
         self.name = "C" + self.name
      elif self.is_nterm and not self.is_cterm and self.name in NUCLEIC:
         self.name += '5'
      elif self.is_cterm and not self.is_nterm and self.name in NUCLEIC:
         self.name += '3'
   
   def OFF(self):
      """ Returns a string that is this residue/unit's entry in OFF file format. For OFF file format
          specification, see http://ambermd.org/doc/OFF_file_format.txt """

      # Write out the atoms table section
      ret_string = ("!entry.%s.unit.atoms table  str name  str type  int typex " +
                   "int resx  int flags  int seq  int elmnt  dbl chg\n") % self.name
      for i in range(self.natom):
         ret_string += ' "%s" "%s" 0 1 %d %d %s %.6f\n' % (self.atoms[i], self.atom_types[i], 0x20000,
                                                           i+1, periodic_table.AtomicNum[self.elements[i]], 
                                                           self.charges[i])

      # Write out perturbed atoms information table
      ret_string += ("!entry.%s.unit.atomspertinfo table  str pname  str ptype  " +
                     "int ptypex  int pelmnt  dbl pchg\n") % self.name
      for i in range(self.natom):
         ret_string += ' "%s" "%s" 0 -1 0.0\n' % (self.atoms[i], self.atom_types[i])

      # Write out boundbox array
      ret_string += "!entry.%s.unit.boundbox array dbl\n" % self.name
      ret_string += " -1.000000\n 0.0\n 0.0\n 0.0\n 0.0\n"

      # Write out childsequence section
      ret_string += "!entry.%s.unit.childsequence single int\n 2\n" % self.name

      # Write out connect array (connect0, connect1)
      ret_string += "!entry.%s.unit.connect array int\n %d\n %d\n" % (self.name, self.head, self.tail)

      # Write out connectivity table
      ret_string += "!entry.%s.unit.connectivity table  int atom1x  int atom2x  int flags\n" % self.name
      for i in range(len(self.bonds)):
         for partner in self.bonds[i]:
            if partner > i+1:
               ret_string += " %d %d 1\n" % (i+1, partner)

      # Write out hierarchy table
      ret_string += ('!entry.%s.unit.hierarchy table  str abovetype  int abovex  str belowtype  ' +
                     'int belowx\n "U" 0 "R" 1\n') % self.name
      for i in range(self.natom):
         ret_string += ' "R" 1 "A" %d\n' % (i+1)
         
      # Write out unit name
      ret_string += '!entry.%s.unit.name single str\n "%s"\n' % (self.name, self.name)
      
      # Write out positions table
      ret_string += '!entry.%s.unit.positions table  dbl x  dbl y  dbl z\n' % self.name
      for i in range(self.natom):
         index = i * 3
         ret_string += " %.6g %.6g %.6g\n" % (self.coords[index], self.coords[index+1], self.coords[index+2])

      # Write out residueconnect array
      ret_string += ("!entry.%s.unit.residueconnect table  int c1x  int c2x  " + 
                     "int c3x  int c4x  int c5x  int c6x\n") % self.name
      for i in range(6):
         try:
            ret_string += ' %d' % self.connects[i]
         except IndexError:
            ret_string += ' 0'
      ret_string += '\n'
      
      # Write out residues table
      if self.orig_name in AMINO: attype = 'p'
      elif self.orig_name in NUCLEIC: attype = 'n'
      elif self.orig_name == 'WAT': attype = 'w'
      else: attype = '?'
      ret_string += ("!entry.%s.unit.residues table  str name  int seq  int childseq  " +
                     "int startatomx  str restype  int imagingx\n") % self.name
      ret_string += ' "%s" 1 %d 1 "%s" 0\n' % (self.name, self.natom+1, attype)

      # Write out PDB Sequence Number array
      ret_string += "!entry.%s.unit.residuesPdbSequenceNumber array int\n 1\n" % self.name

      # Write out solventcap array
      ret_string += "!entry.%s.unit.solventcap array dbl\n -1.000000\n 0.0\n 0.0\n 0.0\n 0.0\n" % self.name

      # Write out velocities table
      ret_string += "!entry.%s.unit.velocities table  dbl x  dbl y  dbl z\n" % self.name
      for i in range(self.natom):
         ret_string += " 0.0 0.0 0.0\n"

      # Done writing this entry

      return ret_string

   def __eq__(self, other):
      """ 2 Residue classes are equivalent if they have the same number of atoms,
          bonds, same name, both head and tail are either both 0 or both non-zero,
          and the connects array is the same length """
      if self.name != other.name:
         return False
      
      if len(self.atoms) != len(other.atoms) or len(self.bonds) != len(other.bonds):
         return False

      if (self.head == 0 and other.head != 0) or (self.head != 0 and other.head == 0):
         return False

      if (self.tail == 0 and other.tail != 0) or (self.tail != 0 and other.tail == 0):
         return False

      if len(self.connects) != len(other.connects):
         return False

      return True

   def __ne__(self, other):
      """ 2 Residue classes are not equivalent if they are not equivalent :) """
      return not self.__eq__(other)

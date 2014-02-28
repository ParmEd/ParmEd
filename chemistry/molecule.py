""" 
Molecule class for manipulating biomolecular structures and storing
molecular data for use in various transformations.
"""

from chemistry import exceptions

class Molecule:
   """ 
   Molecule class that contains biomolecular data:
      o  atoms       : sequential list of atoms in an array
      o  atom_types  : atom types associated with each atom
      o  elements    : element of each atom
      o  charges     : partial charges on each atom
      o  residues    : array of residue names
      o  bonds       : array in which each element lists the atom indices
                       of the other atoms that atom is bonded to
      o  residue_pointers : index of first atom of each residue
      o  residue_container: which residue each atom belongs to
      o  coords      : cartesian coordinates (x1,y1,z1,x2,y2,z2,...) of
                       each atom in the molecule
   """

   def __init__(self, atoms=[], atom_types=[], charges=[], residues=[], bonds=[], 
                residue_pointers=[], coords=[], elements=[], title='',radii=[]):
      """ Initializing and checking the molecular data """
      self.atoms = atoms
      self.residues = residues
      self.residue_pointers = residue_pointers
      self.coords = coords
      self.bonds = bonds
      self.atom_types = atom_types
      self.charges = charges
      self.elements = elements
      self.title = title
      self.radii = radii
      
      self._check()
      if self.valid:
         self._fillcontainer()

   def DeleteAtom(self, atomno):
      """ Deletes an atom from the molecule, removing all bonds it's involved with
          and translating all bonds by one """
      # Remove atomno atom from all of the arrays
      self.atoms.pop(atomno)
      self.coords.pop(atomno*3 + 2)
      self.coords.pop(atomno*3 + 1)
      self.coords.pop(atomno*3    )
      self.atom_types.pop(atomno)
      self.charges.pop(atomno)
      self.elements.pop(atomno)
      if len(self.radii) > 0:
         self.radii.pop(atomno)

      # Remove atomno from all of the bonds
      self.bonds.pop(atomno)
      for i in range(len(self.bonds)):
         for j in range(len(self.bonds[i])):
            k = len(self.bonds[i]) - 1 - j # go in reverse order
            if self.bonds[i][k] == atomno:
               self.bonds[i].pop(k)
            elif self.bonds[i][k] > atomno:
               self.bonds[i][k] -= 1
      
      # Adjust residue_pointers
      for i in range(self.residue_container[atomno] + 1, len(self.residue_pointers)):
         self.residue_pointers[i] -= 1

      # Get rid of now-vacant residues if we removed the last atom from that residue
      if self.residue_pointers[self.residue_container[atomno]+1] - \
         self.residue_pointers[self.residue_container[atomno]+1] == 0:
         self.residue_pointers.pop(self.residue_container[atomno])
         self.residues.pop(self.residue_container[atomno])

      # Now get rid of atomno from residue_container
      self.residue_container.pop(atomno)


   def _check(self):
      """ Checks for consistency in the molecule """
      self.valid = False
      if len(self.atoms) != len(self.atom_types):
         raise(exceptions.MoleculeError('len(atoms) != len(atom_types)'))

      if len(self.charges) == 0:
         for i in range(len(self.atoms)):
            self.charges.append(0.0)

      if len(self.charges) != len(self.atoms):
         raise(exceptions.MoleculeError('len(atoms) != len(charges)'))

      if len(self.coords) != 3 * len(self.atoms):
         raise(exceptions.MoleculeError('len(atoms) != len(coords) * 3'))

      if len(self.residue_pointers) != len(self.residues):
         raise(exceptions.MoleculeError('len(residue_pointers) != len(residues)'))

      if len(self.elements) != len(self.atoms):
         raise(exceptions.MoleculeError('len(elements) != len(atoms)'))

      if len(self.radii) != 0 and len(self.radii) != len(self.atoms):
         raise(exceptions.MoleculeError('len(radii) != len(atoms)'))

      self.valid = True

   def _fillcontainer(self):
      """ Fills residue_container so we know what residue each atom is in """
      self.residue_container = []
      # Fill residues 1 - nres-1
      for i in range(len(self.residues)-1):
         curres = self.residue_pointers[i]
         nextres = self.residue_pointers[i+1]
         for j in range(curres, nextres):
            self.residue_container.append(i)
      # Fill the last residue
      for i in range(self.residue_pointers[len(self.residue_pointers)-1], len(self.atoms)):
         self.residue_container.append(len(self.residues)-1)


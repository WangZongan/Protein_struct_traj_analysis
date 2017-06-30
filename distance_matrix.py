#!/usr/bin/env python
import sys
import Bio.PDB 

def make_dij(structure, activeatomlist=['CA'], cutoff=1000, colsep=" ", defbfactor=1.00, onlylikeatoms=True):
    """
    This function builds a pairwise distance matrix from a PDB file or open file handle pointing to a PDB file.
            activeatomlist  (List of Strings) contains a list of atoms to be considered for the distance matrix.
            cutoff          (Float or Int)    contains a distance in angstroms beyond which pairs are not included
            colsep          (String)          contains a separator to use between columns
            defbfactor      (Float)           contains a weight to include in the final column. 
            onlylikeatoms   (Bool)            if set as true, only like atom distances will be included. 
                                              e.g. CA-CA and CB-CB but no CA-CB
    """   
    aa_def_types={
            'N': 1,'CA': 2,'C': 3,'O': 4,'HN': 5,'HA1': 6,'CAP': 7,'CB': 8,
            'HA2': 9,'CG': 10,'CG1': 11,'CG2': 12,'OG': 13,'OG1': 14,'SG': 15,'CD': 16,
            'CD1': 17,'CD2': 18,'OD': 19,'OD1': 20,'OD2': 21,'ND': 22,'ND1': 23,'ND2': 24,
            'SD': 25,'CE': 26,'CE1': 27,'CE2': 28,'CE3': 29,'OE': 30,'OE1': 31,'OE2': 32,
            'NE': 33,'NE1': 34,'NE2': 35,'CZ': 36,'CZ2': 37,'CZ3': 38,'NZ': 39,'CH2': 40,
            'OH': 41,'NH': 42,'NH1': 43,'NH2': 44,'H': 45,'X': 46}

    protein = Bio.PDB.PDBParser(QUIET=True).get_structure('protein', structure)

    atomlist=[]
    for atom in protein.get_atoms():
        if (atom.get_name() in activeatomlist):             
            atomlist.append(atom) 

    #building residue ref
    residueid={}
    currentresidue=0
    for res in protein.get_residues():
        het = res.get_id()[0]
        resnum = res.get_id()[1]
        if het == ' ':  # check to make sure not a water or other heteroatom
            residueid[resnum] = currentresidue 
            currentresidue += 1

    distlist=[]
    for i in range(len(atomlist)):
        atom1 = atomlist[i]
        res1 = residueid[atom1.get_parent().get_id()[1]]
        atom1code = aa_def_types[atom1.get_name()]
        for j in range(len(atomlist))[i+1:]:
            atom2 = atomlist[j]
            res2 = residueid[atom2.get_parent().get_id()[1]]
            atom2code = aa_def_types[atom2.get_name()]

            dist = atom2-atom1
            disttup = (res1,atom1code,res2,atom2code,dist,defbfactor)
            diststring = colsep.join( [str(x) for x in disttup] )

            if dist < cutoff:
                if ((onlylikeatoms and atom1code==atom2code) or (not onlylikeatoms)):
                    distlist.append( diststring ) 

    distlist.insert(0, str(len(distlist)))
    return "\n".join(distlist)

def main():
    filename = sys.argv[1][:-4]+'.dij'
    f = open(filename,"w")
    print >> f, make_dij(sys.argv[1]) 
    f.close()

if __name__=="__main__":
    main()

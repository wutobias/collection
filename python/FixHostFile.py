from rdkit import Chem
from rdkit.Chem import AllChem
import sys

### This script renumbers and reorders the atoms in a mol2/pdb file
### of a cyclodextrin molecule and writes it out in pdb format. Note
### that it won't write a connection table at the end of the file.

if len(sys.argv)!=2:
    print('Usage: FixHostFile.py \"host.pdb\" > \"host-fixed.pdb\"')
    exit(0)

input = sys.argv[1]
if input.endswith(".pdb"):
    mol = Chem.MolFromPDBFile(input, removeHs=False)
elif input.endswith(".mol2"):
    mol = Chem.MolFromMol2File(input, removeHs=False)
else:
    print("File type not understood")
    exit(1)

gluc_smi = '[#6]1(-[#8])-[#6](-[#8,#7])-[#6](-[#8,#7])-[#6](-[#8])-[#6](-[#6]-[#8,#7])-[#8]-1'
gluc     = Chem.MolFromSmarts(gluc_smi)
NameList = [0]*mol.GetNumAtoms()
ResNames = [0]*mol.GetNumAtoms()
CoordsList = [0]*mol.GetNumAtoms()
EleNames = [0]*mol.GetNumAtoms()
ResIdx   = 1
Res2Atom = dict()
Conf     = mol.GetConformer()
for match in mol.GetSubstructMatches(gluc):
    Res2Atom[ResIdx] = list()
    for idx in range(gluc.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(match[idx])
        if NameList[match[idx]]==0:
            NameList[match[idx]] = atom.GetSymbol()+str(idx+1)
            ResNames[match[idx]] = "GL"+str(ResIdx)
            Res2Atom[ResIdx].append(match[idx])
            pos                    = Conf.GetAtomPosition(match[idx])
            CoordsList[match[idx]] = [pos.x, pos.y, pos.z]
            EleNames[match[idx]]   = atom.GetSymbol()
        nidx_shift = 0
        for natom in atom.GetNeighbors():
            nidx  = natom.GetIdx()
            if NameList[nidx]==0 and natom.GetSymbol()=="H":
                nidx  = natom.GetIdx()
                NameList[nidx] = "H"+str(idx+1+nidx_shift)
                ResNames[nidx] = "GL"+str(ResIdx)
                Res2Atom[ResIdx].append(nidx)
                EleNames[nidx]   = "H"
                pos              = Conf.GetAtomPosition(nidx)
                CoordsList[nidx] = [pos.x, pos.y, pos.z]
                nidx_shift    += 1
    ResIdx += 1

Res2Atom[ResIdx] = list()
for atom in mol.GetAtoms():
    idx   = atom.GetIdx()
    if NameList[idx]==0:
        NameList[idx] = atom.GetSymbol()+str(idx+gluc.GetNumAtoms()+1)
        ResNames[idx] = "SDC"
        Res2Atom[ResIdx].append(idx)
        EleNames[idx]   = atom.GetSymbol()
        pos             = Conf.GetAtomPosition(idx)
        CoordsList[idx] = [pos.x, pos.y, pos.z]
        
AtomCount = 1
for residx in range(1,ResIdx+1):
    for idx in Res2Atom[residx]:
        line    = list("HETATM")
        line   += [" "]*74
        
        name    = NameList[idx].ljust(4," ")
        idxname = str(AtomCount).rjust(5, " ")
        resname = ResNames[idx].ljust(3," ")
        residx  = str(residx).rjust(4," ")
        ele     = str(EleNames[idx]).rjust(2, " ")
        x       = ("%.3f" %CoordsList[idx][0]).rjust(8, " ")
        y       = ("%.3f" %CoordsList[idx][1]).rjust(8, " ")
        z       = ("%.3f" %CoordsList[idx][2]).rjust(8, " ")
        occ     = "1.00".rjust(6, " ")
        bfac    = "0.00".rjust(6, " ")
        ### Atom id
        for i in range(5):
            line[6+i]  = idxname[i]
        ### Edit atom name
        for i in range(4):
            line[13+i] = name[i]
        ### Edit residue name
        for i in range(3):
            line[17+i] = resname[i]
        ### Edit residue idx
        for i in range(4):
            line[23+i] = residx[i]
        ### x coords
        for i in range(8):
            line[30+i] = x[i]
        ### y coords
        for i in range(8):
            line[38+i] = y[i]
        ### z coords
        for i in range(8):
            line[46+i] = z[i]
        ### occupancy
        for i in range(6):
            line[54+i] = occ[i]
        ### b factor
        for i in range(6):
            line[60+i] = bfac[i]
        ### Element
        for i in range(2):
            line[76+i] = ele[i]
        print("".join(line))
        AtomCount += 1
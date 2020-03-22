This is a collection of bash scripts. Many of them depend on paths to programs or libraries in the computer system in my previous lab. Most of them are pretty standard and should be easy to adapt.

#### resp_it_mol2/resp_it_mol2.sh
Script for calculating RESP charges. Requires AmberTools, MOE and fconv (http://pc1664.pharmazie.uni-marburg.de/download/).

#### smi2mol.sh
Convert SMILES to mol2 file.

#### moeminimize.sh
quick minimization with MOE.

#### moepreparePDB.sh
Prepare PDB file using MOE. Splits apart protein and ligand.

#### obabel_conformers.sh
Generate conformers using openbabel.

#### compare_halfs.sh
Quickly analyze trajectories by comparing RMSD of first and second half.

#### pca_cpptraj.sh
Perform PCA with cpptraj.

#### count_wat.sh
Count water molecules in PDB file.

#### resetscratch.sh
Reset gaussian scratch directory.
import mdtraj as md
import numpy as np
import argparse
from scipy.spatial import cKDTree


def centrar_proteina(prot,nom_prot):
    centro = prot.xyz[0].mean(axis=0)
    prot_centrada = prot.atom_slice(range(prot.n_atoms))  # Crea copia segura
    prot_centrada.xyz[0] -= centro
    prot_centrada.save_pdb(f"{nom_prot}_center.pdb")
    return prot_centrada


def Unir_proteina_box(prot_name,prot,solv):
    #### Selec de atomos ####
    protein_atoms = prot.topology.select("protein")
    solvent_sel_str = " or ".join(f"resname {res}" for res in args.solventes)
    solvent_atoms = solv.topology.select(solvent_sel_str)

    #### Coordenadas ####
    coords_prot = prot.xyz[0][protein_atoms]
    coords_solv = solv.xyz[0]

    #### Buscar colisiones ####
    dist_min_nm = args.distancia / 10.0 # Distancia minima en nm (mdtraj usa nm)
    tree = cKDTree(coords_prot)
    dists, _ = tree.query(coords_solv, distance_upper_bound=dist_min_nm)

    mask = dists > dist_min_nm
    keep_atom_indices = np.where(mask)[0]
    residues_to_keep = set(solv.topology.atom(i).residue.index for i in keep_atom_indices)
    atom_indices_to_keep = [atom.index for atom in solv.topology.atoms if atom.residue.index in residues_to_keep]

    
    solv_filtrado = solv.atom_slice(atom_indices_to_keep)

    # Combinar la proteína con el solvente filtrado
    sistema_completo = prot.stack(solv_filtrado)
    sistema_completo.save_pdb(f"{prot_name}_SV.pdb")



if __name__ == '__main__':

    #### Cargar Argumentos ####

    parser = argparse.ArgumentParser(description="Eliminar solapamientos entre proteína y caja de solventes")
    parser.add_argument("-p", "--proteina", required=True, help="Archivo PDB de la proteína")
    parser.add_argument("-b", "--box", required=True, help="Archivo PDB de la caja de solventes")
    parser.add_argument("-d", "--distancia", type=float, default=4, help="Distancia mínima (en Å) para eliminar solventes")
    parser.add_argument("-s", "--solventes", nargs="+", default=["WAT", "ETH" , "PHE"], help="Lista de resnames de solventes")

    args = parser.parse_args()

    prot = md.load_pdb(args.proteina)
    solv = md.load_pdb(args.box)

    #### Centrar Proteina ####
    prot_name = args.proteina.split(".")[0]
    prot = centrar_proteina(prot , prot_name)

    #### Solvatar #####
 
    Unir_proteina_box(prot_name,prot,solv)







# This app is build by Subodh Niroula. This converts SMILES string to 2d and 3D molecular structure

import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
import py3Dmol
from stmol import showmol


#Title and subtitle
st.title("Molecular Structure Drawing Tool")
st.write("Draw molecules using SMILES notation or choose from examples")

# Example molecules with verified SMILES 
example_molecules = {
    "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Paracetamol": "CC(=O)NC1=CC=C(O)C=C1",
    "Water": "O",
    "Ethanol": "CCO",
    "Benzene": "c1ccccc1",
    "Methane": "C",
    "Glucose": "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O"
}

# Input method selection
input_method = st.radio(
    "Choose input method:",
    ["Select from examples", "Draw using SMILES"]
)

# Initialize SMILES variable
smiles = None

if input_method == "Draw using SMILES":
    smiles = st.text_input("Enter SMILES notation:", placeholder="e.g., CCO for ethanol")
else:
    # Display example molecules 
    selected_molecule = st.selectbox(
        "Select a molecule:",
        list(example_molecules.keys())
    )
    if selected_molecule:
        smiles = example_molecules[selected_molecule]
        st.text_input("SMILES notation:", value=smiles, disabled=True)

#Function that converts smiles string to 3d structure
def smiles_to_3d(smiles):
    try:
        # Generate a molecule object from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
                # Generate 2D coordinates
                AllChem.Compute2DCoords(mol)

                # Add hydrogens to the molecule
                mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol)  # UFF optimization
        
        return mol
    except Exception as e:
        st.error(f"An error occurred: {e}")
        return None

# Convert SMILES to 3D if input is provided
mol_3d = None
if smiles:
    mol_3d = smiles_to_3d(smiles)

# Display 2D and 3D structures
if mol_3d:
    # Display the 2D structure
    st.subheader("2D Structure")
    img = Draw.MolToImage(mol_3d)
    st.image(img)

    # Display the 3D structure
    st.subheader("3D Structure")
    mol_block = Chem.MolToMolBlock(mol_3d)
    viewer = py3Dmol.view(width=500, height=400)
    viewer.addModel(mol_block, "mol")
    viewer.setStyle({"stick": {}})
    viewer.zoomTo()
    viewer.spin(True)
    showmol(viewer, height=500, width=800)


# Download button 
if mol_3d:
    pdb_block = Chem.MolToPDBBlock(mol_3d)
    st.download_button("Download 3D Structure as PDB", data=pdb_block, file_name="molecule.pdb", mime="chemical/x-pdb")

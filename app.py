# This app is built by Subodh Niroula.
# It converts a SMILES string into 2D and 3D molecular structures.

import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import py3Dmol
from stmol import showmol

# Title and subtitle
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
    "Glucose": "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O",
}

# Input method selection
input_method = st.radio(
    "Choose input method:",
    ["Select from examples", "Draw using SMILES"],
)

# Initialize SMILES variable
smiles = None
if input_method == "Draw using SMILES":
    smiles = st.text_input("Enter SMILES notation:", placeholder="e.g., CCO for ethanol")
else:
    selected_molecule = st.selectbox("Select a molecule:", list(example_molecules.keys()))
    if selected_molecule:
        smiles = example_molecules[selected_molecule]
        st.text_input("SMILES notation:", value=smiles, disabled=True)


def smiles_to_3d(smiles_string):
    """Convert a SMILES string into an RDKit molecule with 3D coordinates."""
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            st.error("Invalid SMILES notation. Please check your input and try again.")
            return None

        # Generate 2D coordinates before adding hydrogens, so the 2D depiction
        # (drawn later from this same molecule) stays clean and readable.
        AllChem.Compute2DCoords(mol)

        # Add explicit hydrogens, then embed and optimize a 3D conformer.
        mol = Chem.AddHs(mol)
        embed_result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if embed_result != 0:
            st.error("Could not generate 3D coordinates for this molecule.")
            return None

        AllChem.UFFOptimizeMolecule(mol)
        return mol
    except Exception as e:
        st.error(f"An error occurred: {e}")
        return None


# Convert SMILES to a 3D-ready molecule if input is provided
mol_3d = None
if smiles:
    mol_3d = smiles_to_3d(smiles)

# Display 2D and 3D structures
if mol_3d:
    st.subheader("2D Structure")
    img = Draw.MolToImage(mol_3d)
    st.image(img)

    st.subheader("3D Structure")
    mol_block = Chem.MolToMolBlock(mol_3d)
    viewer = py3Dmol.view(width=500, height=400)
    viewer.addModel(mol_block, "mol")
    viewer.setStyle({"stick": {}})
    viewer.zoomTo()
    viewer.spin(True)
    showmol(viewer, height=500, width=800)

    # Download button
    pdb_block = Chem.MolToPDBBlock(mol_3d)
    st.download_button(
        "Download 3D Structure as PDB",
        data=pdb_block,
        file_name="molecule.pdb",
        mime="chemical/x-pdb",
    )

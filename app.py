# This app is built by Subodh Niroula.
# It converts a molecule name or SMILES string into 2D and 3D molecular structures.

import requests
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import py3Dmol
from stmol import showmol

st.set_page_config(page_title="Molecular Structure Drawing Tool", layout="centered")

# ---------------------------------------------------------------- styling ---
# Same palette and type as index.html
st.markdown(
    """
    <style>
      :root{
        --paper:#F5F3EE; --panel:#FFFFFF; --ink:#1C1E1A; --muted:#6B6A61;
        --line:#DAD5C8; --accent:#2B4F49; --accent-soft:#E4ECE9; --error:#A3392A;
      }
      .stApp{ background:var(--paper); color:var(--ink); }
      header[data-testid="stHeader"]{ background:transparent; }
      .block-container{ max-width:760px; padding-top:3rem; padding-bottom:5rem; }

      .kicker{ font-family:'SFMono-Regular',Consolas,Menlo,monospace; font-size:13px;
               color:var(--accent); letter-spacing:.02em; margin:0 0 10px; }
      h1.title{ font-family:Georgia,'Times New Roman',serif; font-weight:600;
                font-size:38px; line-height:1.15; margin:0 0 14px; color:var(--ink); }
      .lede{ color:var(--muted); font-size:16px; max-width:56ch; margin:0 0 28px; }
      h2.sect{ font-family:Georgia,'Times New Roman',serif; font-size:20px;
               font-weight:600; margin:0 0 12px; color:var(--ink); }

      /* panels */
      div[data-testid="stVerticalBlockBorderWrapper"]{
        background:var(--panel); border:1px solid var(--line) !important;
        border-radius:10px; }

      /* inputs */
      div[data-baseweb="input"] > div, div[data-baseweb="select"] > div{
        background:var(--paper) !important; border:1px solid var(--line) !important;
        border-radius:7px !important; }
      input{ font-family:'SFMono-Regular',Consolas,Menlo,monospace !important; }
      label, .stRadio label p{ color:var(--muted); font-size:13px; }

      /* mode selector as tab-style pills */
      div[role="radiogroup"]{ gap:8px; flex-wrap:wrap; }
      div[role="radiogroup"] > label{
        border:1px solid var(--line); border-radius:7px; padding:6px 14px;
        background:transparent; cursor:pointer; }
      div[role="radiogroup"] > label > div:first-child{ display:none; }
      div[role="radiogroup"] > label p{ font-size:14px; color:var(--muted); margin:0; }
      div[role="radiogroup"] > label:has(input:checked){
        background:var(--accent); border-color:var(--accent); }
      div[role="radiogroup"] > label:has(input:checked) p{ color:#fff; }

      /* buttons */
      .stDownloadButton button{
        background:transparent; color:var(--accent);
        border:1px solid var(--accent); border-radius:7px; font-size:13px; }
      .stDownloadButton button:hover{ background:var(--accent-soft); color:var(--accent);
        border-color:var(--accent); }

      /* 2D drawing stage */
      .stage{ background:var(--paper); border:1px solid var(--line); border-radius:8px;
              min-height:280px; display:flex; align-items:center;
              justify-content:center; overflow:hidden; }
      .stage svg{ max-width:100%; height:auto; }

      .note{ background:var(--accent-soft); color:var(--accent); border-radius:7px;
             padding:10px 12px; font-size:13.5px; margin-top:12px; }
      footer{ visibility:hidden; }
    </style>
    <p class="kicker">Name / SMILES → Structure</p>
    <h1 class="title">Molecular Structure Drawing Tool</h1>
    <p class="lede">Search a molecule by name or SMILES and get a 2D diagram and an
    interactive 3D model, with downloads.</p>
    """,
    unsafe_allow_html=True,
)

# --------------------------------------------------------------- examples ---
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


# ---------------------------------------------------------------- helpers ---
@st.cache_data(show_spinner=False)
def smiles_from_name(name):
    """Look up a molecule name on PubChem. Returns (smiles, cid) or (None, None)."""
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{requests.utils.quote(name)}/property/SMILES,ConnectivitySMILES/JSON"
    )
    try:
        r = requests.get(url, timeout=15)
        if r.status_code != 200:
            return None, None
        p = r.json()["PropertyTable"]["Properties"][0]
        smiles = p.get("SMILES") or p.get("ConnectivitySMILES")
        return smiles, p.get("CID")
    except Exception:
        return None, None


@st.cache_data(show_spinner=False)
def build_molecule(smiles):
    """Return (2D SVG, 3D mol block) for a SMILES string, or (None, None) if invalid."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None

    # 2D depiction from the molecule without explicit hydrogens
    flat = Chem.Mol(mol)
    AllChem.Compute2DCoords(flat)
    drawer = rdMolDraw2D.MolDraw2DSVG(480, 380)
    drawer.drawOptions().setBackgroundColour((0.96, 0.95, 0.93, 1))
    drawer.DrawMolecule(flat)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    # 3D conformer with hydrogens, then force-field optimization
    mol3d = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol3d, randomSeed=42) != 0:
        if AllChem.EmbedMolecule(mol3d, randomSeed=42, useRandomCoords=True) != 0:
            return svg, None
    try:
        AllChem.UFFOptimizeMolecule(mol3d)
    except Exception:
        pass
    return svg, Chem.MolToMolBlock(mol3d)


def file_base(label):
    keep = "".join(c.lower() if c.isalnum() else "-" for c in label).strip("-")
    return keep or "molecule"


# ------------------------------------------------------------------ input ---
with st.container(border=True):
    mode = st.radio(
        "Input method",
        ["Examples", "Search by name", "Draw using SMILES"],
        horizontal=True,
        label_visibility="collapsed",
    )

    smiles, label, info = None, "molecule", None

    if mode == "Examples":
        label = st.selectbox("Select a molecule", list(example_molecules.keys()))
        smiles = example_molecules[label]

    elif mode == "Search by name":
        name = st.text_input(
            "Enter a molecule name", placeholder="e.g., ibuprofen, dopamine, cholesterol"
        ).strip()
        if name:
            with st.spinner("Searching PubChem…"):
                smiles, cid = smiles_from_name(name)
            if smiles:
                label = name
                info = f"Found <b>{name}</b> (PubChem CID {cid})."
            else:
                st.error(f'No molecule found for "{name}". Try a common or IUPAC name.')

    else:
        smiles = st.text_input(
            "Enter SMILES notation", placeholder="e.g., CCO for ethanol"
        ).strip() or None
        label = "custom-smiles"

    if info:
        st.markdown(f'<div class="note">{info}</div>', unsafe_allow_html=True)

# ---------------------------------------------------------------- results ---
if smiles:
    svg, mol_block = build_molecule(smiles)

    if svg is None:
        st.error("Invalid SMILES notation. Please check your input and try again.")
    else:
        base = file_base(label)
        col2d, col3d = st.columns(2)

        with col2d, st.container(border=True):
            st.markdown('<h2 class="sect">2D structure</h2>', unsafe_allow_html=True)
            st.markdown(f'<div class="stage">{svg}</div>', unsafe_allow_html=True)
            st.download_button(
                "Download SVG",
                data=svg,
                file_name=f"{base}-2d.svg",
                mime="image/svg+xml",
            )

        with col3d, st.container(border=True):
            st.markdown('<h2 class="sect">3D structure</h2>', unsafe_allow_html=True)
            if mol_block is None:
                st.error("Could not generate 3D coordinates for this molecule.")
            else:
                viewer = py3Dmol.view(width=340, height=280)
                viewer.addModel(mol_block, "mol")
                viewer.setStyle({}, {"stick": {"radius": 0.15}, "sphere": {"scale": 0.25}})
                viewer.setBackgroundColor("0xF5F3EE")
                viewer.zoomTo()
                viewer.spin(True)
                showmol(viewer, height=280, width=340)

                mol3d = Chem.MolFromMolBlock(mol_block, removeHs=False)
                pdb_block = Chem.MolToPDBBlock(mol3d)
                sdf_block = mol_block.rstrip() + "\n$$$$\n"
                d1, d2 = st.columns(2)
                d1.download_button(
                    "Download SDF",
                    data=sdf_block,
                    file_name=f"{base}-3d.sdf",
                    mime="chemical/x-mdl-sdfile",
                )
                d2.download_button(
                    "Download PDB",
                    data=pdb_block,
                    file_name=f"{base}-3d.pdb",
                    mime="chemical/x-pdb",
                )

st.markdown(
    '<p style="margin-top:36px;font-size:13px;color:#6B6A61;">Built by Subodh Niroula. '
    "Name search uses PubChem. 2D and 3D structures are generated locally with RDKit.</p>",
    unsafe_allow_html=True,
)
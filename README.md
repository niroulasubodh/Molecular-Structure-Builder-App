# Molecular Structure Drawing Tool

A small web app that turns a SMILES string into a 2D diagram and an interactive, rotatable 3D molecular model. Type or select a molecule and see both representations side by side.

Built by **Subodh Niroula**.

## What it does

- Accepts a [SMILES](https://en.wikipedia.org/wiki/Simplified_Molecular_Input_Line_Entry_System) string, either typed in directly or picked from a built-in list of example molecules (aspirin, caffeine, paracetamol, glucose, and more).
- Renders a static 2D structural diagram of the molecule.
- Generates a 3D conformer and displays it as an interactive model you can rotate and zoom.
- Lets you download the generated 3D structure as a `.pdb` file.

## Two ways to run this

This repo includes two versions of the tool, for two different situations.

### 1. `app.py` — the full Python app

This is the original tool, built with [Streamlit](https://streamlit.io/), [RDKit](https://www.rdkit.org/), and [py3Dmol](https://pypi.org/project/py3Dmol/). It generates 3D coordinates on demand for *any* valid SMILES string using RDKit's conformer embedding and force-field optimization.

**Run it locally:**

```bash
pip install streamlit rdkit py3Dmol stmol
streamlit run app.py
```

Then open the local URL Streamlit prints in your terminal (usually `http://localhost:8501`).

**Deploy it so others can use it from a link:** Streamlit apps need a running server, so GitHub alone can't host it — push this repo to [Streamlit Community Cloud](https://share.streamlit.io) (free) or [Hugging Face Spaces](https://huggingface.co/spaces), then link the live URL from this README.

### 2. `index.html` — a static, browser-only version

This is a self-contained HTML page with no server and no Python required. It runs entirely in the browser using two WebAssembly/JS libraries:

- [RDKit.js](https://github.com/rdkit/rdkit-js) for parsing SMILES and drawing the 2D structure
- [3Dmol.js](https://3dmol.org/) for the interactive 3D viewer

Because it's just static files, GitHub Pages can host it directly, and it will actually run and respond to input for anyone who opens the link — unlike a plain code preview.

**View it locally:** just open `index.html` in a browser.

**Publish it via GitHub Pages:**

1. Push `index.html` to this repository (repo root, or a `/docs` folder).
2. In the repo, go to **Settings → Pages**.
3. Under "Build and deployment," set **Source** to "Deploy from a branch," pick your branch, and the root (or `/docs`) folder.
4. Save. GitHub will publish it at `https://<your-username>.github.io/<repo-name>/`.
5. Add that link to the top of this README so visitors can try it straight from GitHub.

**One limitation to know about:** the browser-only version has full, reliable 3D generation for the built-in example molecules (their 3D coordinates are precomputed with RDKit ahead of time and bundled into the page). For a custom SMILES string you type in, it tries a live lookup against PubChem's database for a matching 3D structure, and falls back to showing just the 2D diagram if no match is found. The Python app (`app.py`) doesn't have this limitation — it generates 3D coordinates for any valid molecule, since it runs RDKit's real conformer-generation algorithm rather than looking one up.

## Tech stack

| Piece | Python app | Static page |
|---|---|---|
| SMILES parsing & 2D depiction | RDKit | RDKit.js |
| 3D conformer generation | RDKit (`AllChem.EmbedMolecule` + UFF) | Precomputed (examples) / PubChem lookup (custom) |
| 3D viewer | py3Dmol / stmol | 3Dmol.js |
| Interface | Streamlit | Plain HTML/CSS/JS |

## Example molecules included

Aspirin, caffeine, paracetamol, water, ethanol, benzene, methane, and glucose.


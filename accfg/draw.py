from .main import AccFG
from .compare import compare_mols, get_RascalMCES
from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D
import matplotlib as mpl
import networkx as nx

from IPython.display import Image
from PIL import Image as pilImage
from PIL import ImageDraw, ImageFont
import io

def draw_RascalMCES(smiles1, smiles2, legends=None, subImgSize=(500, 500)):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    res = get_RascalMCES(smiles1, smiles2)
    highlightAtomLists = [[x[0] for x in res[0].atomMatches()], 
                          [x[1] for x in res[0].atomMatches()]]
    highlightBondLists = [[x[0] for x in res[0].bondMatches()], 
                          [x[1] for x in res[0].bondMatches()]]
    img = Draw.MolsToGridImage([mol1, mol2], 
                               legends=legends,
                               highlightAtomLists=highlightAtomLists, 
                               highlightBondLists=highlightBondLists,
                               useSVG=False,
                               subImgSize=subImgSize)
    return img

def set_alpha(color, alpha):
    return tuple(list(color[:3]) + [alpha])

def show_atom_idx(smi, label = 'molAtomMapNumber'):
    #https://chemicbook.com/2021/03/01/how-to-show-atom-numbers-in-rdkit-molecule.html
    if isinstance(smi, str):
        mol  = Chem.MolFromSmiles(smi)
    else:
        mol = smi
    for atom in mol.GetAtoms():
        atom.SetProp(label,str(atom.GetIdx()))
    return mol

def draw_mol_with_fgs_dict(smi, fgs_dict, with_legend = True, with_atom_idx = True, alpha = 1, cmp = mpl.colormaps['Pastel1'], img_size = (500, 400)):
    mol = Chem.MolFromSmiles(smi)

    highlight_atoms = {}
    highlight_bonds = {}
    if with_legend:
        legend = ''.join(f"'{key}': {value}\n" for key, value in fgs_dict.items())
    else:
        legend = ''
    for i, (fg, atom_list) in enumerate(fgs_dict.items()):
        for atom_set in atom_list:
            for atom_idx in atom_set:
                highlight_atoms.setdefault(atom_idx, []).append(cmp(i))   
            if len(atom_set) > 1:
                for bond in mol.GetBonds():
                    at1_idx = bond.GetBeginAtomIdx()
                    at2_idx = bond.GetEndAtomIdx()
                    if at1_idx in atom_set and at2_idx in atom_set:
                        highlight_bonds.setdefault(bond.GetIdx(), []).append(set_alpha(cmp(i),alpha))
    if with_atom_idx:
        mol = show_atom_idx(mol)
    # Draw the molecule with highlighted atoms
    d2d = rdMolDraw2D.MolDraw2DCairo(img_size[0], img_size[1])
    dopts = d2d.drawOptions()
    dopts.atomHighlightsAreCircles = True
    dopts.legendFraction = 0.3
    
    d2d.DrawMoleculeWithHighlights(mol, legend, highlight_atoms,highlight_bonds, None, None)
    d2d.GetDrawingText()
    d2d.FinishDrawing()
    return d2d.GetDrawingText()

def draw_mol_with_fgs(smi, afg = AccFG(), canonical=True, with_legend = True, with_atom_idx = True, alpha = 1, cmp = mpl.colormaps['Pastel1'], img_size = (500, 400)):
    
    if canonical:
        smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
    
    fgs = afg.run(smi, canonical=canonical)
    
    mol = Chem.MolFromSmiles(smi)

    highlight_atoms = {}
    highlight_bonds = {}
    if with_legend:
        legend = ''.join(f"'{key}': {value}\n" for key, value in fgs.items())
    else:
        legend = ''
    for i, (fg, atom_list) in enumerate(fgs.items()):
        for atom_set in atom_list:
            for atom_idx in atom_set:
                highlight_atoms.setdefault(atom_idx, []).append(cmp(i))   
            if len(atom_set) > 1:
                for bond in mol.GetBonds():
                    at1_idx = bond.GetBeginAtomIdx()
                    at2_idx = bond.GetEndAtomIdx()
                    if at1_idx in atom_set and at2_idx in atom_set:
                        highlight_bonds.setdefault(bond.GetIdx(), []).append(set_alpha(cmp(i),alpha))
    if with_atom_idx:
        mol = show_atom_idx(mol)
    # Draw the molecule with highlighted atoms
    d2d = rdMolDraw2D.MolDraw2DCairo(img_size[0], img_size[1])
    dopts = d2d.drawOptions()
    dopts.atomHighlightsAreCircles = True
    dopts.legendFraction = 0.3
    
    d2d.DrawMoleculeWithHighlights(mol, legend, highlight_atoms,highlight_bonds, None, None)
    d2d.GetDrawingText()
    d2d.FinishDrawing()
    return d2d.GetDrawingText()

def canonical_smiles(smi):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi))

def draw_compare_mols(smi_1, smi_2, afg = AccFG(), similarityThreshold=0.7, canonical=True, img_size = (500, 400)):
    if canonical:
        smi_1 = canonical_smiles(smi_1)
        smi_2 = canonical_smiles(smi_2)
    
    (target_fgs,target_alkanes),(ref_fgs,ref_alkanes) = compare_mols(smi_1, smi_2, afg, similarityThreshold)
    target_fg_dict = [(fg, atom_list) for fg,_,atom_list in target_fgs]
    target_alkane_dict = [(fg, atom_list) for fg,_,atom_list in target_alkanes]
    target_dict = dict(target_fg_dict + target_alkane_dict)
    
    target_img = draw_mol_with_fgs_dict(smi_1, target_dict, img_size=img_size)
    
    ref_fg_dict = [(fg, atom_list) for fg,_,atom_list in ref_fgs]
    ref_alkane_dict = [(fg, atom_list) for fg,_,atom_list in ref_alkanes]
    ref_dict = dict(ref_fg_dict + ref_alkane_dict)
    ref_img = draw_mol_with_fgs_dict(smi_2, ref_dict, img_size=img_size)
    
    img_list = [molimg(target_img), molimg(ref_img)]
    return img_list

def molimg(x):
    # https://github.com/bbu-imdea/efgs/blob/main/try_efgs.py
    # Utility function to draw image of labeled mol from PNG binary string
    try:
        return pilImage.open(io.BytesIO(x))
    except:
        return None

def img_grid(images, titles = None, num_columns = 5,font = "arial.ttf", font_size = 18, title_height = 30, cell_width = 500, cell_height = 400, bg_color = (255, 255, 255)):
    # https://github.com/bbu-imdea/efgs/blob/main/try_efgs.py
    """
    Creates a grid of images with optional titles under each image.

    Parameters:
    - images (list): List of PIL Image objects.
    - titles (list, optional): List of titles as strings for each image. If None, no titles will be added.
    - num_columns (int): Number of columns in the grid.
    - font_sizes: Font size of titles (if provided).
    - title_height: Height of titles (if provieded).
    - cell_width (int): Width of each image cell.
    - cell_height (int): Height of each image cell (not including title space).
    - bg_color (tuple): Background color as an RGB tuple.

    Returns:
    - PIL Image object containing the image grid.
    """
    # Set title height if titles are provided; otherwise, set to zero
    title_height = title_height if titles else 0
    cell_height_with_title = cell_height + title_height + 2

    # Calculate grid dimensions
    num_images = len(images)
    num_rows = (num_images + num_columns - 1) // num_columns  # Round up for last row

    # Create a blank background image for the grid
    grid_width = num_columns * cell_width
    grid_height = num_rows * cell_height_with_title
    grid_image = pilImage.new("RGB", (grid_width, grid_height), bg_color)

    # Set up font for titles if needed
    if titles:
        try:
            font = ImageFont.truetype(font, font_size)  # Change font and size as desired
        except IOError:
            font = ImageFont.load_default()

    # Iterate through each cell in the grid
    for i in range(num_rows * num_columns):
        row = i // num_columns
        col = i % num_columns
        x = col * cell_width
        y = row * cell_height_with_title

        if i < num_images:
            # Resize the image to fit in the cell
            img_resized = images[i].resize((cell_width, cell_height))
            grid_image.paste(img_resized, (x, y))

            # Draw title if provided
            if titles and i < len(titles):
                draw = ImageDraw.Draw(grid_image)
                title = titles[i]

                # Calculate title position and center it
                text_bbox = draw.textbbox((0, 0), title, font=font)
                text_width = text_bbox[2] - text_bbox[0]
                title_x = x + (cell_width - text_width) // 2  # Center title horizontally
                title_y = y + cell_height  # Position directly beneath the image

                # Draw title text
                draw.text((title_x, title_y), title, fill="black", font=font)
        else:
            # Leave blank space if no image for this cell
            blank_cell = pilImage.new("RGB", (cell_width, cell_height), bg_color)
            grid_image.paste(blank_cell, (x, y))

    return grid_image


def print_directed_graph_as_tree(graph, roots, header='',last=False, show_atom_idx=False):
    """Prints a directed graph as a tree, starting from a given root node."""
    elbow = "└──"
    pipe = "│  "
    tee = "├──"
    blank = "   "
    
    for j,root in enumerate(roots):
        if len(roots) > 1 and j == len(roots) - 1:
            last = True
        if show_atom_idx:
            print(header + (elbow if last else tee) + root+ ": "+str(graph.nodes[root]['mapped_atoms']))
        else:
            print(header + (elbow if last else tee) + root)
        children = list(graph.successors(root))
        if children != []:
            for i,child in enumerate(children):
                # print(child)
                print_directed_graph_as_tree(graph, [child], header=header + (blank if last else pipe), last=i == len(children) - 1, show_atom_idx=show_atom_idx)
                
def print_fg_tree(graph, roots, show_atom_idx=False):
    """Prints a functional group tree."""
    for u, v in list(graph.edges()):
    # Temporarily remove edge (u,v)
        graph.remove_edge(u, v)
        if not nx.has_path(graph, u, v):
            graph.add_edge(u, v)
    print_directed_graph_as_tree(graph, roots, header='',last=False, show_atom_idx=show_atom_idx)
    
    
def SmilesMCStoGridImage(smiles: list[str] or dict[str, str], align_substructure: bool = True, verbose: bool = False, **kwargs):
     # https://bertiewooster.github.io/2022/10/09/RDKit-find-and-highlight-the-maximum-common-substructure-between-molecules.html
     """
     Convert a list (or dictionary) of SMILES strings to an RDKit grid image of the maximum common substructure (MCS) match between them

     :returns: RDKit grid image, and (if verbose=True) MCS SMARTS string and molecule, and list of molecules for input SMILES strings
     :rtype: RDKit grid image, and (if verbose=True) string, molecule, and list of molecules
     :param molecules: The SMARTS molecules to be compared and drawn
     :type molecules: List of (SMARTS) strings, or dictionary of (SMARTS) string: (legend) string pairs
     :param align_substructure: Whether to align the MCS substructures when plotting the molecules; default is True
     :type align_substructure: boolean
     :param verbose: Whether to return verbose output (MCS SMARTS string and molecule, and list of molecules for input SMILES strings); default is False so calling this function will present a grid image automatically
     :type verbose: boolean
     """
     mols = [Chem.MolFromSmiles(smile) for smile in smiles]
     res = rdFMCS.FindMCS(mols, **kwargs)
     mcs_smarts = res.smartsString
     mcs_mol = Chem.MolFromSmarts(res.smartsString)
     smarts = res.smartsString
     smart_mol = Chem.MolFromSmarts(smarts)
     smarts_and_mols = [smart_mol] + mols

     smarts_legend = "Max. substructure match"

     # If user supplies a dictionary, use the values as legend entries for molecules
     if isinstance(smiles, dict):
          mol_legends = [smiles[molecule] for molecule in smiles]
     else:
          mol_legends = ["" for mol in mols]

     legends =  [smarts_legend] + mol_legends
    
     matches = [""] + [mol.GetSubstructMatch(mcs_mol) for mol in mols]

     subms = [x for x in smarts_and_mols if x.HasSubstructMatch(mcs_mol)]

     Chem.Compute2DCoords(mcs_mol)

     if align_substructure:
          for m in subms:
               _ = Chem.GenerateDepictionMatching2DStructure(m, mcs_mol)

     drawing = Draw.MolsToGridImage(smarts_and_mols, highlightAtomLists=matches, legends=legends)

     if verbose:
          return drawing, mcs_smarts, mcs_mol, mols
     else:
          return drawing
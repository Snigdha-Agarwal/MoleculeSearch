from rdkit import Chem
from rdkit.Chem import Draw
from jinja2 import Environment, select_autoescape

env = Environment(
    autoescape=select_autoescape(['html', 'xml'])
)
def printStruct(value):
    return Draw.MolsToGridImage([Chem.MolFromSmiles(value)])

env.filters['printStruct'] = printStruct

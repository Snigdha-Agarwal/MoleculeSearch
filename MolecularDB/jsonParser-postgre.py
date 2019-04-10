import io
import json
import os
import json
import pgdb
from rdkit import Chem
from rdkit.Chem import Draw


def floatOrZero(val):
    if val == "" or val == None:
        return None
    else:
        return float(val)

def main():
    cnx = pgdb.connect(database='molecularDB')
    dbcursor = cnx.cursor()
    mol_info_sql = "INSERT INTO molecular_information (structure) VALUES(%s) where inchi_key = (%s)"
    mol_ener_sql = "INSERT INTO molecular_energies (inchi_key, elec_state, solv_state, total_elec_energy," \
                   " g_energy, h_energy, s_energy, zpve, homo, lumo, geom) " \
                   "VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s::jsonb)"

    path = "../Molecules/"
    infovals = []
    enerVals = []
    for filename in os.listdir(path):
        with io.open(path + filename, 'r', encoding='windows-1252') as json_data:
            print(filename)
            d = json.load(json_data)

            # #adding structure
            # values = dbcursor.execute("SELECT * from molecular_information")
            # for val in values:
            #     inchi_key = val[0]
            #     smile = val[2]
            #     Draw.MolToFile(Chem.MolFromSmiles(smile), "img.png")
            #     with open("img.png", "rb") as image:
            #         f = image.read()
            #         b = pgdb.Binary(f)
            #         dbcursor.execute(mol_info_sql,(b,inchi_key))

            #energy data
            for estate in ['s0','s1','t1','cat-rad']:
                if estate in d:
                    for solv in ['vac','solv']:
                        if solv in d[estate]:
                            geom = json.dumps(d[estate][solv]['geom'])

                            # formatted_geom = []
                            # for item in geom:
                            #     formatted_geom.append(item.replace(',','|'))
                            energies = d[estate][solv]['energies']
                            print(energies['homo'], energies['lumo'], energies['S'], energies['zpve'])
                            enerVals.append((d['inchi-key'],estate,solv,floatOrZero(energies['total_electronic_energy']),
                                             floatOrZero(energies['G']),floatOrZero(energies['H']),
                                             floatOrZero(energies['S']),floatOrZero(energies['zpve']),
                                             floatOrZero(energies['homo']),floatOrZero(energies['lumo']), geom))

    dbcursor.executemany(mol_ener_sql,enerVals)
    cnx.commit()
    
    cnx.close()

if __name__=='__main__':
    main()
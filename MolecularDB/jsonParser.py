import json
import os
import mysql.connector

def floatOrZero(val):
    if val == "":
        return 0
    else:
        return float(val)

def main():
    cnx = mysql.connector.connect(user='root',
                                  host='127.0.0.1',
                                  database='molecularDB')
    dbcursor = cnx.cursor()
    mol_info_sql = "INSERT INTO molecular_information (inchi_key, inchi, smile, formula, molecular_weight, ion_pot, " \
                   "red_pot, ver_en, 0_0_S1, 0_0_T1) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
    mol_ener_sql = "INSERT INTO molecularDB.molecular_energies (inchi_key, elec_state, solv_state, total_elec_energy," \
                   " g_energy, h_energy, s_energy, zpve, homo, lumo, geom) " \
                   "VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"

    path = "../Molecules/"
    infovals = []
    enerVals = []
    for filename in os.listdir(path):
        with open(path +filename) as json_data:
            d = json.load(json_data)

            #create data
            properties = d['properties']
            infovals.append((d['inchi-key'],d['inchi'],d['smiles'],d['formula'],floatOrZero(properties['mw']),
                             floatOrZero(properties['ip']), floatOrZero(properties['rp']),floatOrZero(properties['ve']),
                             floatOrZero(properties['0-0']), floatOrZero(properties['0-0'])))

            #energy data
            for estate in ['s0','s1','t1','cat-rad']:
                if estate in d:
                    for solv in ['vac','solv']:
                        if solv in d[estate]:
                            geom = d[estate][solv]['geom']
                            formatted_geom = []
                            for item in geom:
                                formatted_geom.append(item.replace(',','|'))
                            energies = d[estate][solv]['energies']
                            enerVals.append((d['inchi-key'],estate,solv,floatOrZero(energies['total_electronic_energy']),
                                             floatOrZero(energies['G']),floatOrZero(energies['H']),
                                             floatOrZero(energies['S']),floatOrZero(energies['zpve']),
                                             floatOrZero(energies['homo']),floatOrZero(energies['lumo']),
                                             ','.join(formatted_geom)))

    dbcursor.executemany(mol_info_sql,infovals)
    dbcursor.executemany(mol_ener_sql,enerVals)
    cnx.commit()
    
    cnx.close()

if __name__=='__main__':
    main()
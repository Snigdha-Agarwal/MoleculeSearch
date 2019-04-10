import json
import sys

from flask import Flask, render_template, request
from flask_sqlalchemy import SQLAlchemy

from config import Config
from searchForm import searchMoleculeForm
from flask_bootstrap import Bootstrap

app = Flask(__name__)
app.config.from_object(Config)
bootstrap = Bootstrap(app)
POSTGRES = {
    'user': '',
    'pw': '',
    'db': 'molecularDB',
    'host': 'localhost',
    'port': '5432',
}
app.config['SQLALCHEMY_DATABASE_URI'] = 'postgresql://%(user)s:\
%(pw)s@%(host)s:%(port)s/%(db)s' % POSTGRES
app.config['SQLALCHEMY_ECHO'] = True
db = SQLAlchemy(app)


class MolecularInformation(db.Model):
    __tablename__ = 'molecular_information'
    # enerR = db.relationship("MolecularEnergies", backref='molecular_information', lazy="dynamic")
    inchi_key = db.Column(db.String(50), primary_key=True)
    inchi = db.Column(db.String(100))
    smile = db.Column(db.String(100))
    formula = db.Column(db.String(20))
    molecular_weight = db.Column(db.Float)
    red_pot = db.Column(db.Float)
    ion_pot = db.Column(db.Float)
    ver_en = db.Column(db.Float)
    zero_S1 = db.Column("zero_zero_s1", db.Float)
    zero_T1 = db.Column("zero_zero_t1", db.Float)

    @property
    def serialize(self):
        return {
            "smile": self.smile,
            "inchi_key": self.inchi_key,
            "inchi": self.inchi,
            "formula": self.formula,
            "molecular_weight": self.molecular_weight,
            "red_pot": self.red_pot,
            "ion_pot": self.ion_pot,
            "ver_en": self.ver_en,
            "zero_S1": self.zero_S1,
            "zero_T1": self.zero_T1,
            "energies": self.serialize_many2many
        }

    @property
    def serialize_many2many(self):
        """
        Return object's relations in easily serializable format.
        NB! Calls many2many's serialize property.
        """
        return [item.serialize for item in self.energies]


class MolecularEnergies(db.Model):
    __tablename__ = 'molecular_energies'
    inchi_key = db.Column(db.String(50), db.ForeignKey('molecular_information.inchi_key'), primary_key=True)
    enerR = db.relationship('MolecularInformation', backref=db.backref('energies', lazy='dynamic'))
    elec_state = db.Column(db.String(10), primary_key=True)
    solv_state = db.Column(db.String(10), primary_key=True)
    total_elec_energy = db.Column(db.Float)
    g_energy = db.Column(db.Float)
    h_energy = db.Column(db.Float)
    s_energy = db.Column(db.Float)
    zpve = db.Column(db.Float)
    homo = db.Column(db.Float)
    lumo = db.Column(db.Float)
    geom = db.Column(db.String)

    @property
    def serialize(self):
        """Return object data in easily serializable format"""

        return {
            "elec_state": self.elec_state,
            "solv_state": self.solv_state,
            "total_elec_energy": self.total_elec_energy,
            "g_energy": self.g_energy,
            "h_energy": self.h_energy,
            "s_energy": self.s_energy,
            "zpve": self.zpve,
            "homo": self.homo,
            "lumo": self.lumo,
            "geom": self.geom
        }


@app.route('/', methods=['GET', 'POST'])
def main():
    formResponse = searchMoleculeForm()
    if formResponse.validate_on_submit():
        return search_results(formResponse)
    return render_template('index.html', title='Search', form=formResponse)


@app.route('/results')
def search_results(formResponse):
    lowerRange = sys.float_info.min
    upperRange = sys.float_info.max
    if formResponse.data['lowerRange'] != '':
        lowerRange = float(formResponse.data['lowerRange'])
    if formResponse.data['upperRange'] != '':
        upperRange = float(formResponse.data['upperRange'])
    searchField = formResponse.data['searchField']
    if lowerRange == sys.float_info.min and upperRange == sys.float_info.max:
        result = MolecularInformation.query.join(MolecularEnergies, MolecularEnergies.inchi_key ==
                                                 MolecularInformation.inchi_key).all()
    else:
        result = MolecularInformation.query.join(MolecularEnergies, MolecularEnergies.inchi_key ==
                                                 MolecularInformation.inchi_key).filter(
            get_field(searchField).between(lowerRange, upperRange)).all()
    return render_template('results.html', result=[i.serialize for i in result])


def get_field(searchField):
    if searchField == "Homo":
        return MolecularEnergies.homo
    elif searchField == "Lumo":
        return MolecularEnergies.lumo
    elif searchField == "Redox Potential":
        return MolecularInformation.red_pot


@app.route('/detailed/<molecule>')
def detailed_result(molecule):
    return render_template('details.html', title=molecule, result=json.loads(request.args['object'].replace('\'', "\"")
                                                                             .replace('None', '\"\"')))


if __name__ == "__main__":
    # jinja2.filters.FILTERS['printStruct'] = printStruct
    app.run()

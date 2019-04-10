from flask_wtf import FlaskForm
from wtforms import StringField, SelectField, SubmitField


class searchMoleculeForm(FlaskForm):
    choices = [('Homo', 'Homo'),
               ('Lumo', 'Lumo'),
               ('Redox Potential', 'Redox Potential')]
    searchField = SelectField('Search for molecule with:', choices=choices)
    lowerRange = StringField('Start Range','')
    upperRange = StringField('End Range','')
    submit = SubmitField('Search')

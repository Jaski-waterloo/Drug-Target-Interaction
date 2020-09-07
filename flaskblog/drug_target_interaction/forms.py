from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed
from wtforms import StringField, PasswordField, SubmitField, BooleanField, FieldList, SelectField, RadioField
from wtforms.validators import DataRequired, Length, Email, EqualTo, ValidationError
from flask import g, request, session
from flaskblog import pockets
from flask_login import login_user, current_user, logout_user, login_required

class PDBForm(FlaskForm):
    pdb = FileField('Upload PDB File', validators=[FileAllowed(['pdb']), DataRequired()])
    fasta = FileField('Upload .fasta file', validators=[FileAllowed(['fasta']), DataRequired()])
    ligand = FileField('Upload ligand file(Supported: .smiles, .mol, .mol2, .pdb, .sdf)', validators=
            [FileAllowed(['smiles', 'mol', 'mol2', 'pdb','sdf']), DataRequired()])
    
    pdb_name = StringField("Enter PDB Name", validators=[DataRequired()])
    ligand_name = StringField("Enter Ligand Name",validators=[DataRequired()])

    submit = SubmitField('Submit')  


def generate_form(choices = [(1,1),(2,2)]):
    class BForm(FlaskForm):
        radio = RadioField('Select Pocket', choices = choices, validators=[DataRequired()])
        submit = SubmitField('Submit')
    return BForm()

class BSForm(FlaskForm):
	#temp = [(i[1], i[2], i[3]) for i in pockets]
	temp = [0,1]
	choices = []
	for i in range(len(temp)):
		choices.append((temp[i], temp[i]))
	radio = RadioField('Select Pocket', choices = choices, validators=[DataRequired()])
	submit = SubmitField('Submit')


class DrugBankApproved(FlaskForm):
    pdb = FileField('Upload PDB File', validators=[FileAllowed(['pdb']), DataRequired()])
    fasta = FileField('Upload .fasta file', validators=[FileAllowed(['fasta']), DataRequired()])

    pdb_name = StringField("Enter PDB Name", validators=[DataRequired()])

    submit = SubmitField('Submit')


class DrugBankAll(FlaskForm):
    pdb = FileField('Upload PDB File', validators=[FileAllowed(['pdb']), DataRequired()])
    fasta = FileField('Upload .fasta file', validators=[FileAllowed(['fasta']), DataRequired()])

    pdb_name = StringField("Enter PDB Name", validators=[DataRequired()])

    submit = SubmitField('Submit')

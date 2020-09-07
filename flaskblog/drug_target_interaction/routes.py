from flask import render_template, url_for, flash, redirect, request, Blueprint, g, session, send_file
from flask_login import login_user, current_user, logout_user, login_required
from flaskblog import db, bcrypt
from flaskblog.models import User, Post
from flaskblog.drug_target_interaction.forms import (PDBForm, generate_form, BSForm, DrugBankAll, DrugBankApproved)
from flaskblog.users.utils import save_picture, send_reset_email
from flaskblog.drug_target_interaction.utils import get_pockets, get_pocket_with_id, ConverToLigand, GetAllSequenceFeatures, GetAllLigandFeatures, GetAllBindingSiteFeatures, calcVolume
from rdkit import Chem
from flaskblog import pockets

import pandas as pd
import numpy as np
import pickle

drug_target_interaction = Blueprint('drug_target_interaction', __name__)


@drug_target_interaction.route("/input_pdb", methods=['GET', 'POST'])
@login_required
def input_pdb():
    form = PDBForm()
    if form.validate_on_submit():
        pdb_text = request.files['pdb'].read().decode('utf-8')
        fasta_text = request.files['fasta'].read().decode('utf-8')
        ligand_text = request.files['ligand'].read().decode('utf-8')
        ligand_type = str(request.files['ligand'])
        session['pdb_text'] = pdb_text
        session['fasta_text'] = fasta_text
        session['ligand_text'] = ligand_text
        session['ligand_type'] = ligand_type
        session['pdb_name'] = form.pdb_name.data
        session['ligand_name'] = form.ligand_name.data
        #g.pdb_text = pdb_text
        #g.fasta_text = fasta_text
        #g.ligand_text = ligand_text
        #g.ligand_type = g.ligand_type
        #session['pockets'] = [1,2,3]
        try:
            mol = Chem.MolFromPDBBlock(pdb_text)
            Chem.SanitizeMol(mol)
        except:
            flash('PDB is invalid or Rdkit could not sanitize PDB /nPlease enter another PDB', 'danger')
            return render_template('input_pdb.html', title="Predict", form=form)

        try:
            if 'smiles' in ligand_type:
                ligand_type = 'smiles'
            elif 'mol' in ligand_type:
                ligand_type = 'mol'
            elif 'mol2' in ligand_type:
                ligand_type = 'mol2'
            elif 'pdb' in ligand_type:
                ligand_type = 'pdb'
            elif 'sdf' in ligand_type:
                ligand_type = 'sdf'
            ligand = ConverToLigand(ligand_text, ligand_type)
            session['ligand_type'] = ligand_type
            session['ligand_name'] = "temp"
            #g.ligand_type = ligand_type
        except Exception as e:
            flash(e, 'danger')
            return render_template('input_pdb.html', title='Predict', form=form)

        try:
            #global pockets
            print("getting pockets")
            pockets = get_pockets(pdb_text)
            #user_id = current_user.user_id
            session["pockets"] = pockets
            print("done pockets")

        except Exception as e:
            flash(e, 'danger')
            return render_template('input_pdb.html', title='Predict', form=form)

        flash('This can take up to 5 minuites; If it is taking longer than that, please consider making another prediction', 'danger')

        #return redirect(url_for('drug_target_interaction.test', pockets=pockets))
        #return render_template('test.html', title='Test', data=pockets)
        return redirect(url_for('drug_target_interaction.select_site'))
    return render_template('input_pdb.html', title='Predict', form=form)


@drug_target_interaction.route("/select_site", methods=['GET', 'POST'])
@login_required
def select_site():
    features_flag = [1,1,1,1]
    pwd = '/home/jaskiratsinghbhatia4/12-Error-Pages-1/flaskblog/drug_target_interaction/'
    #temp = session['ligand_type']
    pockets = session['pockets']
    #pockets = g.get('pockets')
    pockets = [(i,i) for i in pockets]
    form = generate_form(pockets)
    if form.validate_on_submit():
        #flash('This might take some time', 'danger')
        selected_option = form.radio.data
        session['selected_id'] = selected_option
        fasta_text = session['fasta_text']
        pdb_text = session['pdb_text']
        ligand_text = session['ligand_text']
        ligand_type = session['ligand_type']
        ligand_name = session['ligand_name']
        ligand = ConverToLigand(ligand_text, ligand_type)

        
        selected_option = selected_option.strip('()').split(', ')
        selected_option[0] = int(selected_option[0])
        selected_option[1] = str(selected_option[1][1:-1])
        selected_option[2] = str(selected_option[2][1:-1])
        selected_option = tuple(selected_option)

        pocket_output = get_pocket_with_id(session['pdb_text'], selected_option)

        #return render_template('test.html' ,data=pocket_output)

        binding_site_features = GetAllBindingSiteFeatures(pocket_output)

        chain_ = binding_site_features[0][1]

        binding_site_features = binding_site_features[1]
        if None in binding_site_features.values() or np.nan in binding_site_features.values():
            features_flag[0] = 0


        import tempfile
        fo = tempfile.NamedTemporaryFile(mode='w+')
        fo.write(fasta_text)
        file_name = fo.name
        fo.seek(0)

        sequence_features = GetAllSequenceFeatures(fo)

        fo.close()

        ligand_features = GetAllLigandFeatures(ligand, ligand_name)

        if ligand.GetNumAtoms() < 17:
            flash('ligand size is small, results may be inaccurate')
        
        if None in ligand_features.values() or np.nan in ligand_features.values():
            features_flag[1] = 0

        #chain_ = binding_site_features[0][1]

        #return render_template('test.html' ,data=chain_)

        flag = 0

        #from itertools import chain
        #sequence_features = list(chain(*sequence_features))
        temp = []

        for features in sequence_features:
            temp.append(features['sequence_Chain'])
            if chain_ in features['sequence_Chain']:
                sequence_features = features
                flag = 1
                break
        #return render_template('test.html' ,data=flag)

        if flag != 1:
            flash('Chain for Pocket not found', 'danger')
            return redirect(url_for('drug_target_interaction.input_pdb'))

        if None in sequence_features.values() or np.nan in sequence_features.values():
            features_flag[2] = 0


        pdb_name = session['pdb_name'].upper()
        selected_option = list(selected_option)
        selected_option = [pdb_name] + selected_option
        selected_option = tuple(selected_option)
        volume = calcVolume(selected_option)[1]

        #return render_template('test.html', data=selected_option)
        
        if volume is None:
            features_flag[3] = 0




        if sum(features_flag) < 4:
            if features_flag[0] == 0:
                flash('Unable to Calculate Binding Site Features','danger')
            if features_flag[1] == 0:
                flash('Unable to Calculate Ligand Features', 'danger')
            if features_flag[2] == 0:
                flash('Unable to calculate Sequence Features', 'danger')
            if features_flag[3] == 0:
                flash('Unable to calculate volume of binding site', 'danger')

        #return render_template('test.html', data=binding_site_features)        

        final_dict_of_features = {}
        for k,v in binding_site_features.items():
            final_dict_of_features[k] = v
        for k,v in sequence_features.items():
            final_dict_of_features[k] = v
        for k,v in ligand_features.items():
            final_dict_of_features[k] = v
        final_dict_of_features['binding_site_Volume'] = volume

        df = pd.DataFrame(final_dict_of_features)
        with open(pwd+'features_to_take.p', 'rb') as p:
            to_take = pickle.load(p)
        to_take.remove('Class')
        df = df[to_take]

        df = df.replace([np.inf, -np.inf], np.nan)

        with open(pwd+"clf_xgb_final.p", "rb") as p:
            clf = pickle.load(p)

        preds1 = clf.predict(df)

        with open(pwd+"clf_svm_final.p", "rb") as p:
            clf = pickle.load(p)

        preds2 = clf.predict(df)

        with open(pwd+"clf_rf_final.p", "rb") as p:
            clf = pickle.load(p)

        preds3 = clf.predict(df)



        return render_template('see_predictions.html', data = [preds1[0], preds2[0], preds3[0]], title="Predictions")
    return render_template('select_site.html', title='Predict', form=form)



#@drug_target_interaction.route("calc_binding_site_features"






@drug_target_interaction.route("/test", methods=['GET', 'POST'])
def test():
    #pockets = session['pockets']
    #pockets = [i for i in g]
    return render_template('test.html', data=pockets)




@drug_target_interaction.route("/input_pdb_approved_drugbank", methods=['GET', 'POST'])
@login_required
def input_pdb_approved_drugbank():
    form = DrugBankApproved()
    if form.validate_on_submit():
        pdb_text = request.files['pdb'].read().decode('utf-8')
        fasta_text = request.files['fasta'].read().decode('utf-8')
        session['pdb_text'] = pdb_text
        session['fasta_text'] = fasta_text
        session['pdb_name'] = form.pdb_name.data

        try:
            mol = Chem.MolFromPDBBlock(pdb_text)
            Chem.SanitizeMol(mol)
        except:
            flash('PDB is invalid or Rdkit could not sanitize PDB /nPlease enter another PDB', 'danger')
            return render_template('input_pdb_approved_drugbank.html', title="Predict", form=form)

        try:
            #global pockets
            print("getting pockets")
            pockets = get_pockets(pdb_text)
            #user_id = current_user.user_id
            session["pockets"] = pockets
            print("done pockets")

        except Exception as e:
            flash(e, 'danger')
            return render_template('input_pdb_approved_drugbank.html', title='Predict', form=form)

        #return redirect(url_for('drug_target_interaction.test', pockets=pockets))
        #return render_template('test.html', title='Test', data=pockets)
        return redirect(url_for('drug_target_interaction.select_site_approved'))
    return render_template('input_pdb_approved_drugbank.html', title='Predict', form=form)


@drug_target_interaction.route("/input_pdb_all_drugbank", methods=['GET', 'POST'])
@login_required
def input_pdb_all_drugbank():
    form = DrugBankAll()
    if form.validate_on_submit():
        pdb_text = request.files['pdb'].read().decode('utf-8')
        fasta_text = request.files['fasta'].read().decode('utf-8')
        session['pdb_text'] = pdb_text
        session['fasta_text'] = fasta_text
        session['pdb_name'] = form.pdb_name.data

        try:
            mol = Chem.MolFromPDBBlock(pdb_text)
            Chem.SanitizeMol(mol)
        except:
            flash('PDB is invalid or Rdkit could not sanitize PDB /nPlease enter another PDB', 'danger')
            return render_template('input_pdb_all_drugbank.html', title="Predict", form=form)

        try:
            #global pockets
            print("getting pockets")
            pockets = get_pockets(pdb_text)
            #user_id = current_user.user_id
            session["pockets"] = pockets
            print("done pockets")

        except Exception as e:
            flash(e, 'danger')
            return render_template('input_pdb_all_drugbank.html', title='Predict', form=form)

        #return redirect(url_for('drug_target_interaction.test', pockets=pockets))
        #return render_template('test.html', title='Test', data=pockets)
        return redirect(url_for('drug_target_interaction.select_site_all'))
    return render_template('input_pdb_all_drugbank.html', title='Predict', form=form)


@drug_target_interaction.route("/select_site_all", methods=['GET', 'POST'])
@login_required
def select_site_all():
    features_flag = [1,1,1]
    pwd = '/home/jaskiratsinghbhatia4/12-Error-Pages-1/flaskblog/drug_target_interaction/'
    #temp = session['ligand_type']
    pockets = session['pockets']
    #pockets = g.get('pockets')
    pockets = [(i,i) for i in pockets]
    form = generate_form(pockets)
    if form.validate_on_submit():
        #flash('This might take some time', 'danger')
        selected_option = form.radio.data
        session['selected_id'] = selected_option
        fasta_text = session['fasta_text']
        pdb_text = session['pdb_text']

        selected_option = selected_option.strip('()').split(', ')
        selected_option[0] = int(selected_option[0])
        selected_option[1] = str(selected_option[1][1:-1])
        selected_option[2] = str(selected_option[2][1:-1])
        selected_option = tuple(selected_option)

        pocket_output = get_pocket_with_id(session['pdb_text'], selected_option)

        #return render_template('test.html' ,data=pocket_output)

        binding_site_features = GetAllBindingSiteFeatures(pocket_output)

        chain_ = binding_site_features[0][1]

        binding_site_features = binding_site_features[1]
        if None in binding_site_features.values() or np.nan in binding_site_features.values():
            features_flag[0] = 0


        import tempfile
        fo = tempfile.NamedTemporaryFile(mode='w+')
        fo.write(fasta_text)
        file_name = fo.name
        fo.seek(0)

        sequence_features = GetAllSequenceFeatures(fo)

        fo.close()



        #chain_ = binding_site_features[0][1]

        #return render_template('test.html' ,data=chain_)

        flag = 0

        #from itertools import chain
        #sequence_features = list(chain(*sequence_features))
        temp = []

        for features in sequence_features:
            temp.append(features['sequence_Chain'])
            if chain_ in features['sequence_Chain']:
                sequence_features = features
                flag = 1
                break
        #return render_template('test.html' ,data=flag)

        if flag != 1:
            flash('Chain for Pocket not found', 'danger')
            return redirect(url_for('drug_target_interaction.input_pdb'))

        if None in sequence_features.values() or np.nan in sequence_features.values():
            features_flag[2] = 0


        pdb_name = session['pdb_name'].upper()
        selected_option = list(selected_option)
        selected_option = [pdb_name] + selected_option
        selected_option = tuple(selected_option)
        volume = calcVolume(selected_option)[1]

        #return render_template('test.html', data=selected_option)

        if volume is None:
            features_flag[1] = 0




        if sum(features_flag) < 3:
            if features_flag[0] == 0:
                flash('Unable to Calculate Binding Site Features','danger')
            if features_flag[1] == 0:
                flash('Unable to calculate volume of binding site', 'danger')
            if features_flag[2] == 0:
                flash('Unable to calculate Sequence Features', 'danger')



        with open(pwd+"LigandFeatures_drugbank_all_with_volume.p", "rb") as p:
            ligands_to_test = pickle.load(p)
        for key in ligands_to_test.keys():
            temp = ligands_to_test[key]
            ligands_to_test[key] = {"ligand_"+k: v for k, v in temp.items()}

        with open(pwd+"final_dict_of_features.p", "rb") as p:
            final_dict_of_features = pickle.load(p)


        #return render_template('test.html', data=binding_site_features)

        binding_site_features['binding_site_Volume'] = volume

        for j in range(len(ligands_to_test)):
            try:
                ligand_f = ligands_to_test[j]
            except:
                continue

            for k,v in binding_site_features.items():
                try:
                    final_dict_of_features[k].append(v)
                except:
                    pass

            for k,v in ligand_f.items():
                final_dict_of_features[k].append(v)

            for k,v in sequence_features.items():
                final_dict_of_features[k].append(v)

        #return render_template('test.html', data = [len(i) for i in final_dict_of_features.values()])

        df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in final_dict_of_features.items() ]))


        #df = pd.DataFrame(final_dict_of_features)
        with open(pwd+'features_to_take.p', 'rb') as p:
            to_take = pickle.load(p)
        to_take.remove('Class')
        to_take.append('ligand_Name')
        #ligand_names = df.pop('ligand_Name')
        df = df[to_take]

        df = df.replace([np.inf, -np.inf], np.nan)
        df.dropna(inplace=True)

        ligand_names = df.pop('ligand_Name')

        with open(pwd+"clf_xgb_final.p", "rb") as p:
            clf = pickle.load(p)

        preds1 = clf.predict(df)

        with open(pwd+"clf_svm_final.p", "rb") as p:
            clf = pickle.load(p)

        preds2 = clf.predict(df)

        with open(pwd+"clf_rf_final.p", "rb") as p:
            clf = pickle.load(p)

        preds3 = clf.predict(df)

        pred_df = pd.DataFrame()

        pred_df["ligand_Name"] = ligand_names
        pred_df["xgb"] = preds1
        pred_df["svc"] = preds2
        pred_df["rf"] = preds3

        from rdkit import Chem

        ligand = Chem.SDMolSupplier(pwd + '3D structures.sdf')

        ligands = [i for i in ligand]

        ligands = [[j,i] for i,j in zip(ligands, range(len(ligands)))]

        ligands = dict(ligands)

        #return render_template('test.html', data = len(pred_df))


        #df['kn'] = [1-i for i in df['kn']]
        #df_gt = df.loc[(df['xgb'] == 1) | (df['rf'] == 1) | (df['ets'] == 1)]
        #pred_df['sum'] = pred_df['xgb'] + pred_df['rf'] + pred_df['svc']
        #df_gt = pred_df.loc[df['sum'] >= 1]
        #df_gt = df.loc[df['kn']==1]
        final = []
        for i in pred_df['ligand_Name']:
            try:
                temp = ligand[i]
                #if temp.GetNumAtoms() > 20:
                final.append((temp.GetProp('DATABASE_ID'), temp.GetProp("GENERIC_NAME"), temp))
            #     if len(final) == 50:
            #         break
            except:
                final.append(('NA', 'NA', 'NA'))


        final_df = {}

        final_df['DrugBank ID'] = [i[0] for i in final]

        #final_df['Smiles'] = [i[2].GetProp("SMILES") for i in final]

        final_df["Generic Name"] = [i[1] for i in final]

        final_df['xgb'] = pred_df['xgb']

        final_df['rf'] = pred_df['rf']

        final_df['svm'] = pred_df['svc']

        session['preds'] = final_df



        return redirect(url_for('drug_target_interaction.download_file_helper'))

    return render_template('select_site.html', title='Predict', form=form)




@drug_target_interaction.route("/select_site_approved", methods=['GET', 'POST'])
@login_required
def select_site_approved():
    features_flag = [1,1,1]
    pwd = '/home/jaskiratsinghbhatia4/12-Error-Pages-1/flaskblog/drug_target_interaction/'
    #temp = session['ligand_type']
    pockets = session['pockets']
    #pockets = g.get('pockets')
    pockets = [(i,i) for i in pockets]
    form = generate_form(pockets)
    if form.validate_on_submit():
        #flash('This might take some time', 'danger')
        selected_option = form.radio.data
        session['selected_id'] = selected_option
        fasta_text = session['fasta_text']
        pdb_text = session['pdb_text']

        selected_option = selected_option.strip('()').split(', ')
        selected_option[0] = int(selected_option[0])
        selected_option[1] = str(selected_option[1][1:-1])
        selected_option[2] = str(selected_option[2][1:-1])
        selected_option = tuple(selected_option)

        pocket_output = get_pocket_with_id(session['pdb_text'], selected_option)

        #return render_template('test.html' ,data=pocket_output)

        binding_site_features = GetAllBindingSiteFeatures(pocket_output)

        chain_ = binding_site_features[0][1]

        binding_site_features = binding_site_features[1]
        if None in binding_site_features.values() or np.nan in binding_site_features.values():
            features_flag[0] = 0


        import tempfile
        fo = tempfile.NamedTemporaryFile(mode='w+')
        fo.write(fasta_text)
        file_name = fo.name
        fo.seek(0)

        sequence_features = GetAllSequenceFeatures(fo)

        fo.close()



        #chain_ = binding_site_features[0][1]

        #return render_template('test.html' ,data=chain_)

        flag = 0

        #from itertools import chain
        #sequence_features = list(chain(*sequence_features))
        temp = []

        for features in sequence_features:
            temp.append(features['sequence_Chain'])
            if chain_ in features['sequence_Chain']:
                sequence_features = features
                flag = 1
                break
        #return render_template('test.html' ,data=flag)

        if flag != 1:
            flash('Chain for Pocket not found', 'danger')
            return redirect(url_for('drug_target_interaction.input_pdb'))

        if None in sequence_features.values() or np.nan in sequence_features.values():
            features_flag[2] = 0


        pdb_name = session['pdb_name'].upper()
        selected_option = list(selected_option)
        selected_option = [pdb_name] + selected_option
        selected_option = tuple(selected_option)
        volume = calcVolume(selected_option)[1]

        #return render_template('test.html', data=selected_option)

        if volume is None:
            features_flag[1] = 0




        if sum(features_flag) < 3:
            if features_flag[0] == 0:
                flash('Unable to Calculate Binding Site Features','danger')
            if features_flag[1] == 0:
                flash('Unable to calculate volume of binding site', 'danger')
            if features_flag[2] == 0:
                flash('Unable to calculate Sequence Features', 'danger')



        with open(pwd+"LigandFeatures_drugbank_approved_with_volume.p", "rb") as p:
            ligands_to_test = pickle.load(p)
        for key in ligands_to_test.keys():
            temp = ligands_to_test[key]
            ligands_to_test[key] = {"ligand_"+k: v for k, v in temp.items()}

        with open(pwd+"final_dict_of_features.p", "rb") as p:
            final_dict_of_features = pickle.load(p)


        #return render_template('test.html', data=binding_site_features)

        binding_site_features['binding_site_Volume'] = volume

        for j in range(len(ligands_to_test)):
            try:
                ligand_f = ligands_to_test[j]
            except:
                continue

            for k,v in binding_site_features.items():
                try:
                    final_dict_of_features[k].append(v)
                except:
                    pass

            for k,v in ligand_f.items():
                final_dict_of_features[k].append(v)

            for k,v in sequence_features.items():
                final_dict_of_features[k].append(v)

        #return render_template('test.html', data = [len(i) for i in final_dict_of_features.values()])

        df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in final_dict_of_features.items() ]))


        #df = pd.DataFrame(final_dict_of_features)
        with open(pwd+'features_to_take.p', 'rb') as p:
            to_take = pickle.load(p)
        to_take.remove('Class')
        to_take.append('ligand_Name')
        #ligand_names = df.pop('ligand_Name')
        df = df[to_take]

        df = df.replace([np.inf, -np.inf], np.nan)
        df.dropna(inplace=True)

        ligand_names = df.pop('ligand_Name')

        with open(pwd+"clf_xgb_final.p", "rb") as p:
            clf = pickle.load(p)

        preds1 = clf.predict(df)

        with open(pwd+"clf_svm_final.p", "rb") as p:
            clf = pickle.load(p)

        preds2 = clf.predict(df)

        with open(pwd+"clf_rf_final.p", "rb") as p:
            clf = pickle.load(p)

        preds3 = clf.predict(df)

        pred_df = pd.DataFrame()

        pred_df["ligand_Name"] = ligand_names
        pred_df["xgb"] = preds1
        pred_df["svc"] = preds2
        pred_df["rf"] = preds3

        from rdkit import Chem

        ligand = Chem.SDMolSupplier(pwd + 'structures.sdf')

        ligands = [i for i in ligand]

        ligands = [[j,i] for i,j in zip(ligands, range(len(ligands)))]

        ligands = dict(ligands)

        #return render_template('test.html', data = len(pred_df))


        #df['kn'] = [1-i for i in df['kn']]
        #df_gt = df.loc[(df['xgb'] == 1) | (df['rf'] == 1) | (df['ets'] == 1)]
        #pred_df['sum'] = pred_df['xgb'] + pred_df['rf'] + pred_df['svc']
        #df_gt = pred_df.loc[df['sum'] >= 1]
        #df_gt = df.loc[df['kn']==1]
        final = []
        for i in pred_df['ligand_Name']:
            try:
                temp = ligand[i]
                #if temp.GetNumAtoms() > 20:
                final.append((temp.GetProp('DATABASE_ID'), temp.GetProp("GENERIC_NAME"), temp))
            #     if len(final) == 50:
            #         break
            except:
                final.append(('NA', 'NA', 'NA'))


        final_df = {}

        final_df['DrugBank ID'] = [i[0] for i in final]

        #final_df['Smiles'] = [i[2].GetProp("SMILES") for i in final]

        final_df["Generic Name"] = [i[1] for i in final]

        final_df['xgb'] = pred_df['xgb']

        final_df['rf'] = pred_df['rf']

        final_df['svm'] = pred_df['svc']

        session['preds'] = final_df



        return redirect(url_for('drug_target_interaction.download_file_helper'))

    return render_template('select_site.html', title='Predict', form=form)




@drug_target_interaction.route("/download_file_helper", methods=['GET', 'POST'])
@login_required
def download_file_helper():
    return render_template('download_file_helper.html', title="Download")



@drug_target_interaction.route("/download_file", methods=['GET', 'POST'])
@login_required
def download_file():
    import tempfile
    dic = session['preds']
    import pandas as pd
    df = pd.DataFrame(dic)
    with tempfile.NamedTemporaryFile(mode='w+') as fo:
        df.to_csv(fo.name, index=False)

        return send_file(fo.name, as_attachment=True)

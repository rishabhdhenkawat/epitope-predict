from __future__ import absolute_import, print_function
import sys, os, math
import numpy as np
import pandas as pd
from collections import OrderedDict
from epitopepredict import peptutils, sequtils
from flask import Flask, render_template, request, url_for, redirect, session
from flask import Flask, request, url_for,Response, render_template_string
#from flask_mail import Mail, Message
from flask import Flask, stream_with_context, request, Response, flash     

app = Flask(__name__)

nlf = pd.read_csv(('data/NLF.csv'),index_col=0)
blosum62 = pd.read_csv('data/blosum62.csv')
blosum50 = pd.read_csv('data/blosum50.csv')

codes = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
         'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


def check_valid_prot_seq(seq):
    for i in seq:
        if i not in mhclearn.codes:
            return False

def aff2log50k(a):
    return 1 - (math.log(a) / math.log(50000))

def log50k2aff(a):
    return np.power(50000,1-a)

def random_encode(p):
    """Random encoding of a peptide for testing"""

    return [np.random.randint(20) for i in pep]

def one_hot_encode(seq):
    """One hot encoding of peptide"""

    o = list(set(codes) - set(seq))
    s = pd.DataFrame(list(seq))
    x = pd.DataFrame(np.zeros((len(seq),len(o)),dtype=int),columns=o)
    a = s[0].str.get_dummies(sep=',')
    a = a.join(x)
    a = a.sort_index(axis=1)
    #a = a.set_index([a.index,list(seq)])
    #show_matrix(a)
    e = a.values.flatten()
    return e

def blosum_encode(seq):
    """Blosum62 encoding of peptide"""

    s=list(seq)
    x = pd.DataFrame([blosum62[i] for i in seq]).reset_index(drop=True)
    x = x.iloc[:,:-4]
    e = x.values.flatten()
    return e

def nlf_encode(seq):
    """NLF encoding of a peptide """

    x = pd.DataFrame([nlf[i] for i in seq]).reset_index(drop=True)
    e = x.values.flatten()
    return e


def auc_score(true, sc, cutoff=None):

    from sklearn import metrics
    if cutoff!=None:
        true = (true<=cutoff).astype(int)
        sc = (sc<=cutoff).astype(int)
    fpr, tpr, thresholds = metrics.roc_curve(true, sc, pos_label=1)
    r = metrics.auc(fpr, tpr)
    return  r


def get_protein_set():
    syf = os.path.join('data/SYF_set.fasta')
    return sequtils.fasta_to_dataframe(syf)

def get_training_set(allele=None, length=None):
    """Get training set for MHC-I data."""

    b = pd.read_csv('data/curated_training_data.no_mass_spec.csv')
    eval1 = get_evaluation_set1()
    df = b.loc[~b.peptide.isin(eval1.peptide)].copy()
    if allele is not None:
        df = b.loc[b.allele==allele].copy()

    df['log50k'] = df.ic50.apply(lambda x: aff2log50k(x))
    df['length'] = df.peptide.str.len()
    if length != None:
        df = df[(df.length==length)]
    df = df[df.ic50<50000]
    df = df[df.measurement_type=='quantitative']
    #df['binder'] = df.loc[df.ic50<500].astype(int)
    return df

def get_evaluation_set1(allele=None, length=None):
    """Get eval set of peptides"""

    e = pd.read_csv('data/binding_data_2013.csv',comment='#')
    if allele is not None:
        e = e[e.allele==allele]
    if length != None:
        e = e[(e.length==length) ]
    e['log50k'] = e.ic50.apply(lambda x: aff2log50k(x)).round(2)
    return e

def prepare_data( df, name, allele):
        """Put raw prediction data into DataFrame and rank,
           override for custom processing. Can be overriden for
           custom data."""

        df['name'] = name
        df['pos'] = df.index+1
        df['allele'] = allele
        get_ranking(df)
        return df

def get_ranking( df):
    """Add a ranking column according to scorekey"""

    s='score'
    df['rank'] = df[s].rank(method='min',ascending=1)
    df.sort_values(by=['rank','name','allele'], ascending=True, inplace=True)
    return df



def get_allele_names():
    """Get allele names from training set"""
    b = get_training_set()
    a = b.allele.value_counts()
    a = a[a>100]
    return list(a.index)

def get_model(allele, length):
    """Get a regression model."""
    try:
        import sklearn
    except:
        print ('you need scikit-learn to use this predictor')
    import joblib
    #os.makedirs(models_path, exist_ok=True)
    allele=allele.replace(":","_")
    allele=allele.replace("*","_")
    fname = os.path.join("weights/"+allele+'_'+str(length)+'.joblib')
    print(fname)
    
    reg = joblib.load(fname)
    return reg

def clean_sequence(seq):
    """clean a sequence of invalid characters before prediction"""
    import re
    new = re.sub('[-*_#X]', '', seq)
    return new

def clean_sequence(seq):
    """clean a sequence of invalid characters before prediction"""
    import re
    new = re.sub('[-*_#X]', '', seq)
    return new
def main(alleleList,length,seq):
    result=[]
    for allele in alleleList:
        for l in length:
            try:
                reg=get_model(allele, l)
                peptides,k=peptutils.create_fragments(seq=seq, length=l, overlap=1)
                s=pd.Series(peptides)
                X = s.apply(lambda x: pd.Series(blosum_encode(x)),1)
                l=reg.predict(X)
                sc = l
                in_=pd.DataFrame()
                in_["peptide"]=peptides
                res = pd.DataFrame(np.column_stack([in_["peptide"],sc]),columns=['peptide','log50k'])
                res['log50k'] = res.log50k.astype('float')
                res['score'] = res.log50k.apply(lambda x: log50k2aff(x))
                #print(res.dtypes)
                allele=allele
                df = prepare_data(res,"temp",allele)
                result.append(df)
            except Exception as e:
                print(e)
                f=pd.DataFrame()
                result.append(f)
    finalresult = pd.concat(result)
    return finalresult       


@app.route("/", methods=['post', 'get'])
def index1():
    headings=get_allele_names()
    headings.append("All")
    return render_template('index.html',headings=headings)

@app.route("/process", methods=['POST', 'get'])
def process():
    if request.method == "POST":
        seq=request.form.get("seq")
        allele=request.form.getlist("allele")
        length=request.form.getlist("length")
        if allele[0]=="All":
            alleleList=get_allele_names()
        else:
            alleleList=allele
        print(seq)
        print(allele)
        print(length)
        test_list = [int(i) for i in length]
        print(test_list)
        seq=seq.replace("\n","").replace(" ","").replace("\t","")
        print(seq)
        df=main(alleleList,test_list,seq)
        print(df)
        df.to_csv("result.csv")
        return render_template('index2.html', tables=[df.to_html(classes='data')], titles=df.columns.values)
from flask import send_file, send_from_directory, safe_join, abort
app.config['UPLOAD_FOLDER']="./"
@app.route('/database_download/result.csv')
def database_download():
    # uploads = os.path.join(app.root_path, app.config['UPLOAD_FOLDER'])
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename="result.csv")



if __name__ == "__main__":
    app.run(host='0.0.0.0', port="8000",debug=True, threaded=True)

from pdfminer.pdfinterp import PDFResourceManager, PDFPageInterpreter
from pdfminer.converter import TextConverter
from pdfminer.layout import LAParams
from pdfminer.pdfpage import PDFPage
import requests
from io import StringIO
import pandas as pd 
import numpy as np 
import os
from spacy.matcher import Matcher
import spacy
import json
from ast import literal_eval

import spacy
from spacy.matcher import Matcher

def email_phone(st):
    print('email new method')
    st = convertPDFToText(st)
    
    nlp = spacy.load("en_core_web_sm")
    phone_matcher = Matcher(nlp.vocab)
    phone_pattern = [{'TEXT': {'REGEX': '^(?:(?:-{0,}|:{0,}))(?:(?:\+|0{0,2})91(\s*[\-]\s*)?|[0]?)?[789]\d{9}$'}}]
    phone_matcher.add('phonematcher', None, phone_pattern)

    doc_st = nlp(st)
    matches = phone_matcher(doc_st)
    ph_numbers = list(set([doc_st[s:e].text for (_, s, e) in matches]))
    if len(ph_numbers) == 1:
        number = ph_numbers[0]
    else:
        number = ', '.join(ph_numbers)
    email_pattern = [{'TEXT': {'REGEX': "^[a-zA-Z0-9+_.-]+@[a-zA-Z0-9.-]+$"}}]
    email_matcher = Matcher(nlp.vocab)
    email_matcher.add('emailmatcher', None, email_pattern)
    matches = email_matcher(doc_st)
    emails = list(set([doc_st[s:e].text for (_, s, e) in matches]))
    if len(emails) == 1:
        email = emails[0]
    else:
        email = ', '.join(emails)

    return email, number




def libreoffice_exec():
            # TODO: Provide support for more platforms
            #if sys.platform == 'darwin':
              #  return '/Applications/LibreOffice.app/Contents/MacOS/soffice'
            return 'libreoffice'

def convert_to(folder, source, timeout=None):
    args = ['soffice', '--headless', '--convert-to', 'pdf', '--outdir', folder, source]
    import subprocess
    import re
    process = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout)
    filename = re.search('-> (.*?) using filter', process.stdout.decode())

    if filename is None:
        raise LibreOfficeError(process.stdout.decode())
    else:
        return filename.group(1)
    


def convertPDFToText_(remote_file):
    ext = remote_file.split('.')[-1].lower()
    if ext=='PDF':
        remote_file=remote_file.split('.')[-2]
        remote_file=remote_file+'.'+'pdf'     
    ext = remote_file.split('.')[-1].lower()
    if ext=='pdf':
        print('Extracting data from: {}'.format(remote_file))
        url = remote_file
        r = requests.get(url, allow_redirects=True)
        path = remote_file.split('/')[-1]
        n=path
        open(path, 'wb').write(r.content)
        rsrcmgr = PDFResourceManager()
        retstr = StringIO()
        codec = 'utf-8'
        laparams = LAParams()
        device = TextConverter(rsrcmgr, retstr, laparams=laparams)
        fp = open(path, 'rb')
        interpreter = PDFPageInterpreter(rsrcmgr, device)
        password = ""
        maxpages = 0
        caching = True
        pagenos = set()
        for page in PDFPage.get_pages(fp, pagenos, maxpages=maxpages, password=password, caching=caching,
                                      check_extractable=True):
            interpreter.process_page(page)
        fp.close()
        device.close()
        string = retstr.getvalue()
        retstr.close()
        #os.remove(n)
        return string
    else :
        print('Extracting data from: {}'.format(remote_file))
        url = remote_file
        r = requests.get(url, allow_redirects=True)
        path = remote_file.split('/')[-1]
        n = path
        open(path, 'wb').write(r.content)
        import fitz
        convert_to('./' ,n)
        pdf_name=n
        pdf_name=pdf_name.split('.')[-2]+'.pdf'
        doc = fitz.open(pdf_name)
        string=''
        for i in range(0,doc.pageCount):
            string=string+doc.getPageText(pno=i)
        
        return string
"""    
# don't use this function - just for reference
def get_keywords(jd):
    url = "https://jobfit.perspectico.com/parser"
    input_dict = {"description": {"0": jd}, "skills": {"0": ""}}
    j = json.dumps(input_dict)
    data = {"inputFile": j}
    try:
        response = requests.post(url, data)
        result = response.json()
        # print(result)
        if result['Success']:
            return result['Keywords']
        else:
            print("ERROR: Keywords not fetched")
            result = []
            return result
    except:
        print("ERROR")
        result = []
        return result
"""    
def find_keywords(text):
    bad_chars = ["!", ";", "$", "#", "&", "%", "\n", "\uf0b7", ""]
    try:
        text = text.lower().strip()
    except:
        pass
    for i in bad_chars:
        text = text.replace(i,"")
    resume_words = set(text.split())
    return resume_words



from resume_parser import *
#result={}
major_list, minor_list=[],[]
def show_result_thread(url): # , major_list, minor_list):
    jd_words = ""
    with open('skills_dict.json', 'r') as fp:
        skills_dict = json.load(fp)
    with open('skills_set.txt','r') as f:
        skills_set = literal_eval(f.read())
    with open('traits_set.txt','r') as f2:
        traits_set = literal_eval(f2.read())
    #global result
    result = {}
    global major_list
    global minor_list
    st_list = []
    
    try:
                resume_name = url.split("/")[-1]
                result[resume_name] = {}
                st = convertPDFToText_(url)   
                result[resume_name]["Name"]=temp[0]['name']
                result[resume_name]["Phone Number"]=temp[0]['mobile_number']
                result[resume_name]["Email ID"]=temp[0]['email']
                
                
    except Exception as e:
                print(e)
    
                #print(e)
                resume_name = url.split("/")[-1]
                print("fail: {}".format(resume_name))
                result[resume_name]["text"]='NA'
                result[resume_name]["Name"]='NA'
                result[resume_name]["Phone Number"]='NA'
                result[resume_name]["Email ID"]='NA'
                pass
        #print(result)        
        #return result
            
    return result
            
            
import requests
import time
import concurrent.futures     
def show_result(resume_urls, major, minor, jd_words=""):
            #global result
            #result={}
            global major_list
            major_list=major
            global minor_list
            minor_list= minor
            result={}
            with concurrent.futures.ProcessPoolExecutor() as executor:
                        #executor.map(show_result_thread, resume_urls)
                        for d in executor.map(show_result_thread, resume_urls):
#                                     print(d)
#                                     print(type(d))
                                    result.update(d)
                        
#             finalresult={}
#             for i in resume_urls:
#                         resume_name = i.split("/")[-1]
#                         if result[resume_name]["Name"]!='NA':
#                                     finalresult[resume_name]=result[resume_name]
                                    
                        
                        
            #print("from main",result)
            return result
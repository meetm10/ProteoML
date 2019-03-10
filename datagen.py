import csv
import pandas as pd
import scipy.stats as st
import math
import statistics

def get_seq():
    
    normal_file = open('human_proteome.txt','r').read()
    normal_file = normal_file.split('>')

    normal_data = []

    for entry in normal_file:
        if entry == "":
            continue
        entry = entry.replace('\n','')
        normal_data.append(entry.split('|'))

    normal_dict = {}

    for entry in normal_data:
        if entry[-1] in normal_dict:
            continue
        normal_dict[entry[-1]] = [[x for x in entry[:-1]]]
        
    normal_seq = [key for key in normal_dict.keys()]

    disprot_idp_file = open('disprot.csv','r',encoding="utf8")
    disprot_idp = csv.reader(disprot_idp_file)
    disprot_idp_get_index = [1,5,8,9,11,12,14]

    data = []

    for entry in disprot_idp:
        data.append(entry)

    disprot_idp_cols = []

    temp_data = []
    for x in disprot_idp_get_index:
        temp_data.append(data[0][x])
    disprot_idp_cols.append(temp_data)

    disprot_idp_data = []

    for i in range(1,len(data)):
        temp_data = []
        for x in disprot_idp_get_index:
            temp_data.append(data[i][x])
        disprot_idp_data.append(temp_data)

    # print(disprot_idp_cols)
    
    disprot_idp_seq = [disprot_idp_data[x][4] for x in range(len(disprot_idp_data))]

    disprot_idr_file = open('disprot_idr.csv','r',encoding="utf8")
    disprot_idr = csv.reader(disprot_idr_file)
    disprot_idr_get_index = [1,2,3,8,10,12]

    data = []

    for entry in disprot_idr:
        data.append(entry)

    disprot_idr_cols = []

    temp_data = []

    for x in disprot_idr_get_index:
        temp_data.append(data[0][x])
    disprot_idr_cols.append(temp_data)    

    disprot_idr_data = []

    for i in range(1,len(data)):
        temp_data = []
        for x in disprot_idr_get_index:
            temp_data.append(data[i][x])
        disprot_idr_data.append(temp_data)

    # print(disprot_idr_cols)
    
    disprot_idr_seq = [disprot_idr_data[x][5] for x in range(len(disprot_idr_data))]
    
    return disprot_idp_seq, disprot_idr_seq, normal_seq

def std_dev(A,mean):
    if(len(A)> 0):
        A[:] = [(x - mean)**2 for x in A]
        return math.sqrt(sum(A)/len(A))
    else:
        return math.inf

def build_eveness_data(seq,amino_acids):
    A = []
    for x in amino_acids:
        A_ind = [i for i, ltr in enumerate(seq) if ltr == x]
        A.extend(A_ind)
    if(len(A) > 0):
        mean = len(seq)/2
        std_dev_A = std_dev(A, mean)
        median = statistics.median(A)
        if(std_dev_A > 0):
            skewness = ((mean-median)*3)/std_dev_A
        else:
            skewness = 0
    else:
        skewness = 0
    return skewness

def neg_charge_density(s):
    A=[]
    i=0
    for x in ['D','E']:
        A_ind = []
        for i, ltr in enumerate(s):
            if i == 0:
                if ltr == x and (s[i+1] == 'D' or s[i+1] == 'E'):
                    A_ind.append(i)
            elif i == len(s)-1:
                if ltr == x and (s[i-1] == 'D' or s[i-1] == 'E'):
                    A_ind.append(i)
            else:
                if ltr == x and (s[i-1] == 'D' or s[i-1] == 'E' or s[i+1] == 'D' or s[i+1] == 'E'):
                    A_ind.append(i)                
        A.extend(A_ind)   
    return len(set(A))/len(s)

def build_feature_dict(seq_list,skew,charge_density):
    
    hydrophobic = {'A': 1, 'V': 1, 'L': 2, 'I': 1, 'M': 2, 'F': 3, 'Y': 3, 'W': 3}
    hydrophobic_acids = ['A','V','L','I','M','F','Y','W']

    uncharged = {'G': 1, 'S': 1, 'T': 1, 'N': 1, 'Q': 1}
    uncharged_acids = ['G','S','T','N','Q']

    pos_charged = {'H': 1, 'K': 1, 'R': 1}
    pos_charged_acids = ['H','K','R']

    neg_charged = {'D': 1, 'E': 1}
    neg_charged_acids = ['D','E']

    both_charged = {'H': 1, 'K': 1, 'R': 1, 'D': 1, 'E': 1}
    both_charged_acids = ['H','K','R','D','E']

    special = {'G': 1, 'P': 1, 'C': 1}
    special_acids = ['G','P','C']

    aromatic = {'F': 1, 'Y': 2, 'W': 2, 'H': 1, 'P': 1}
    aromatic_acids = ['F','Y','W','H','P']

    aliphatic = {'A': 1,'V': 1,'L': 1,'I': 1,'E': 1,'M': 1,'S': 1,'T': 1,'N': 1,'Q': 1,'G': 1,'C': 1,'K': 1,'R': 1,'D': 1}
    aliphatic_acids = ['A','V','L','I','E','M','S','T','N','Q','G','C','K','R','D']
    
    target_dict = {}
    
    for s in seq_list:
        if len(s)!=0:
            target_dict[s] = {}

            count = 0
            for x in hydrophobic_acids:
                count += s.count(x)*hydrophobic[x]
            target_dict[s]['hydrophobic'] = count/len(s)

            count = 0
            for x in uncharged_acids:
                count += s.count(x)*uncharged[x]
            target_dict[s]['uncharged'] = count/len(s)

            count = 0
            for x in pos_charged_acids:
                count += s.count(x)*pos_charged[x]
            target_dict[s]['pos_charged'] = count/len(s)

            count = 0
            for x in neg_charged_acids:
                count += s.count(x)*neg_charged[x]
            target_dict[s]['neg_charged'] = count/len(s)

            count = 0
            for x in both_charged_acids:
                count += s.count(x)*both_charged[x]
            target_dict[s]['both_charged'] = count/len(s)

            count = 0
            for x in special:
                count += s.count(x)*special[x]
            target_dict[s]['special'] = count/len(s)

            count = 0
            for x in aromatic_acids:
                count += s.count(x)*aromatic[x]
            target_dict[s]['aromatic'] = count/len(s)

            count = 0
            for x in aliphatic_acids:
                count += s.count(x)*aliphatic[x]
            target_dict[s]['aliphatic'] = count/len(s)
            
            if skew:
                target_dict[s]['pos_skew'] = build_eveness_data(s,pos_charged_acids)
                target_dict[s]['neg_skew'] = build_eveness_data(s,neg_charged_acids)
            
            if charge_density:
                target_dict[s]['neg_charge_density'] = neg_charge_density(s)
        
    return target_dict


def fetch_train_test_data():
    
    disprot_idp_seq, disprot_idr_seq, normal_seq = get_seq()

    skew = True
    charge_density = True
    
    disprot_idp_feature_dict = build_feature_dict(disprot_idp_seq,skew,charge_density)
    disprot_idr_feature_dict = build_feature_dict(disprot_idr_seq,skew,charge_density)
    normal_feature_dict = build_feature_dict(normal_seq,skew,charge_density)
    
    normalData = pd.DataFrame.from_dict(normal_feature_dict[key] for key in normal_feature_dict.keys())
    normalData['label'] = 1
    idpData = pd.DataFrame.from_dict(disprot_idp_feature_dict[key] for key in disprot_idp_feature_dict.keys())
    idrData = pd.DataFrame.from_dict(disprot_idr_feature_dict[key] for key in disprot_idr_feature_dict.keys())
    abnormalData = pd.concat([idpData, idrData])
    abnormalData['label'] = 0
    X_normal = normalData.iloc[:, :-1].values
    X_abnormal = abnormalData.iloc[:, :-1].values
    Y_normal = normalData.iloc[:, -1].values
    Y_abnormal = abnormalData.iloc[:, -1].values
    
    return X_normal, X_abnormal, Y_normal, Y_abnormal

def validate():
    disprot_idp_seq, disprot_idr_seq, normal_seq = get_seq()
    
    skew = True
    charge_density = True

    disprot_idp_feature_dict = build_feature_dict(disprot_idp_seq,skew,charge_density)
    disprot_idr_feature_dict = build_feature_dict(disprot_idr_seq,skew,charge_density)
    normal_feature_dict = build_feature_dict(normal_seq,skew,charge_density)
    
    before_result = []
    after_result = []
    for entry in normal_seq:
        regexp = re.compile(r'P[SGYQ]{0,4}G')
        if len(entry) > 100 and (entry.count('S')+entry.count('Y')+entry.count('G')+entry.count('Q'))/len(entry) > 0.4 and normal_feature_dict[entry]['hydrophobic']>0.16 and normal_feature_dict[entry]['aromatic']>0.08 and normal_feature_dict[entry]['neg_charged']>=0.007:
            before_result.append(entry)
            if regexp.search(entry): 
                after_result.append(entry)
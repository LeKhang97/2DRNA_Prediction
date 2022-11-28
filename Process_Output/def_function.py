import numpy as np
import pandas as pd
import json
import math
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import product

def pairing_prob(seq, threshold = 0.05):
    pairing_dict = {}
    (propensity,ensemble_energy) = RNA.pf_fold(seq)  
    for i,j in itertools.combinations(range(len(seq)+1),2):
        prob = RNA.get_pr(i, j)
        if prob != 0:
            if -math.log10(prob) < -math.log10(threshold):
                pairing_dict[(j-1,i-1)] = round(prob,5)
    
    return pairing_dict

def flatten(lol):
    l = []
    for i in lol:
        l += i
    
    return l

def mattews_corr_coeff(l):
	tp = l[0]
	tn = l[1]
	fp = l[2]
	fn = l[3]
	if (tp+fp == 0):
		print("We have an issue : no positives detected ! (linear structure)")
	if (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn) == 0:
		return None
	else:
		x = (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)


	return (tp*tn-fp*fn) / math.sqrt(x)

def accuracy(tp, tn, fp, fn):
	return (tp+tn)/(tp+fp+tn+fn)

def recall_sensitivity(tp, tn, fp, fn):
	return tp/(tp+fn)

def specificity(tp, tn, fp, fn):
	return tn/(tn+fp)

def precision_ppv(tp, tn, fp, fn):
	return tp/(tp+fp)

def npv(tp, tn, fp, fn):
	return tn/(tn+fn)

def f1_score(tp, tn, fp, fn):
	return 2*tp/(2*tp+fp+fn)

def dbn_to_basepairs(structure):
	parenthesis = []
	brackets = []
	braces = []
	rafters = []
	basepairs = []
	As = []
	Bs = []
	for i, c in enumerate(structure):
		if c == '(':
			parenthesis.append(i)
		if c == '[':
			brackets.append(i)
		if c == '{':
			braces.append(i)
		if c == '<':
			rafters.append(i)
		if c == 'A':
			As.append(i)
		if c == 'B':
			Bs.append(i)
		if c == '.':
			continue
		if c == ')':
			basepairs.append((i, parenthesis.pop()))
		if c == ']':
			basepairs.append((i, brackets.pop()))
		if c == '}':
			basepairs.append((i, braces.pop()))
		if c == '>':
			basepairs.append((i, rafters.pop()))
		if c == 'a':
			basepairs.append((i, As.pop()))
		if c == 'b':
			basepairs.append((i, Bs.pop()))
	return basepairs

def compare_two_structures(true2d, prediction, pred_struct = True):
	true_basepairs = dbn_to_basepairs(true2d)
	if pred_struct == True:
		pred_basepairs = dbn_to_basepairs(prediction)
	else:
		pred_basepairs = prediction
        
	tp = 0
	fp = 0
	tn = 0
	fn = 0
	for bp in true_basepairs:
		if bp in pred_basepairs:
			tp += 1
		else:
			fn += 1
	for bp in pred_basepairs:
		if bp not in true_basepairs:
			fp += 1
	tn = len(true2d) * (len(true2d) - 1) * 0.5 - fp - fn - tp
	return [tp, tn, fp, fn]

def compare_two_contacts(true2d, prediction):
	true_contacts = [p for p,v in enumerate(true2d) if v == "*"]
	pred_contacts = [p for p,v in enumerate(prediction) if v == "*"]
	tp = 0
	fp = 0
	tn = 0
	fn = 0
	for bp in true_contacts:
		if bp in pred_contacts:
			tp += 1
		else:
			fn += 1
	for bp in pred_contacts:
		if bp not in true_contacts:
			fp += 1
	tn = len(true2d) - fp - fn - tp
	return [tp, tn, fp, fn]

def get_contact(file):
    file = [x for x in file if x != '']
    contact_list = []
    for p,v in enumerate(file):
        if (p%2 == 1) & (p >= 3):
            contact_list += [v] 

    return contact_list

def get_struct(file, contact = True):
    file = [x for x in file if x != '']
    struct_list = []
    for p,v in enumerate(file):
        if contact == True:
            if (p%2 == 0) & (p >= 2):
                struct_list += [v.split("\t")[0].split("+")[0].replace(" ","")] 
        else:
            if p >= 2:
                struct_list += [v.split("\t")[0].split("+")[0].replace(" ","")] 

    return struct_list

def get_score(file, contact = True):
    file = [x for x in file if x != '']
    score_list = []
    for p,v in enumerate(file):
        if contact == True:
            if (p%2 == 0) & (p >= 2):
                score_list += [[float(i) for i in v.split("\t")[1:]]]
        else:
            if p >= 2:
                score_list += [[float(i) for i in v.split("\t")[1:]]]

    return score_list

def get_motif(file, contact = True):
    file = [x for x in file if x != '']
    motif_list = []
    for p,v in enumerate(file):
        if contact == True:
            if (p%2 == 0) & (p >= 2):
                motif_list += [[i.replace(" ", "") for i in v.split("\t")[0].split("+")[1:]]]
        else:
            if p >= 2:
                motif_list += [[i.replace(" ", "") for i in v.split("\t")[0].split("+")[1:]]]

    return motif_list

def get_filtered_sol(get_motif_file, filter_motif):
    keep_list = []
    for seq_name in filter_motif:
        if seq_name in file[0]:
            for p,v in enumerate(get_motif_file):
                if all([not(motif in filter_motif[seq_name]) for motif in v]):
                    keep_list += [p]
            break
    return keep_list

def normalize(get_score_file):
    normalized_score = []
    max_score = []
    num_obj = len(get_score_file[0])
    for i in range(num_obj):
        max_score += [max([scores[i] for scores in get_score_file])]
    
    for scores in get_score_file:
        normalized_score += [[scores[i]/max_score[i] if max_score[i] != 0 else 0 for i in range(num_obj)]]
        
    return normalized_score

def n_align(s,c, identity = 0):
    pos_dict = {'value': [], 'pos': []}
    for i in range(len(s)-len(c) + 1):
        mm = 0
        for j in range(len(c)):
            if s[i:(i+len(c))][j] != c[j]:
                mm += 1
        if (len(c)-mm)/len(c) >= identity:
            pos_dict['pos'] += [[i+j for j in range(len(c))]]
            #pos_dict['pos'] += [i]
            pos_dict['value'] += [(len(c)-mm)/len(c)]
  
    return pos_dict

def multi_n_align(s,c, global_identity = 0.5, local_identity = 0):
    c = c.replace('&&','&')
    c_list = c.split('&')
    pos = []
    value = []
    pos_dict = {'value': [], 'pos': []}
    gap_pos = [len(c_list[i]) for i in range(len(c_list))]

    for i in range(len(c_list)):
        local_align = n_align(s,c_list[i],local_identity)
        pos += [local_align['pos']]
        value += [local_align['value']]

  
    for i,t in zip(product(*pos),product(*value)):
        gap_align = [i[j+1][0] - i[j][0] for j in range(len(i) - 1)]
        threshold = sum(np.multiply(t,gap_pos))/sum([len(u) for u in c_list])
        if all([gap_align[j] >= gap_pos[j] for j in range(len(gap_pos) - 1)]) and threshold >= global_identity:
            pos_dict['pos'] += [i]
            pos_dict['value'] += [sum(np.multiply(t,gap_pos))/sum([len(u) for u in c_list])]

    return pos_dict

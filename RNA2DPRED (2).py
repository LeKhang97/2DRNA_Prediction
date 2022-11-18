import sys
import os
import RNA
import json
import math
import itertools
from itertools import product
from docplex.mp.model import Model
from motif_alignment import *
from operator import itemgetter

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

def checking_pseudo(l):
    for i in itertools.combinations(l,2):
        if (((min(i[0]) in range(min(i[1]), max(i[1]) + 1)) and not(max(i[0]) in range(min(i[1]), max(i[1]) + 1))) or
            ((min(i[1]) in range(min(i[0]), max(i[0]) + 1)) and not(max(i[1]) in range(min(i[0]), max(i[0]) + 1)))):
            return True
    return False

def visual_structure(l, length):
    pseudo_depth = [('(',')'), ('[',']'), ('{','}'), ('<','>'),
                ('a','a'), ('b','b'), ('c','c'), ('d','d')]
    sorted_l = sorted(l, key=itemgetter(1))
    depth = 1
    dict_depth = {}
    seq = '.'* length
    for p,v in enumerate(sorted_l):
        flag = 0
        if p >= 1:
            for i in range(1,depth + 1):
                if checking_pseudo(list(dict_depth[str(i)]) + [v]) == False:
                    dict_depth[str(i)] += [v]
                    flag = 1
                    break
            if flag == 0:
                depth += 1
                dict_depth[str(depth)] = [v]
        else:
            dict_depth[str(depth)] = [v]

    for depth in dict_depth:
        for pair in dict_depth[depth]:
            seq = seq[:min(pair)] + pseudo_depth[int(depth)-1][0] + seq[min(pair)+1: max(pair)] + pseudo_depth[int(depth)-1][1] + seq[max(pair)+1:]
            
    return seq

def valid_motif(struct, pseudoknot = False):
    if pseudoknot == False:
        if any((i[0] in struct) and (i[1] in struct) for i in itertools.combinations(['(','[','<'],2)):
            return False
    
    if all(i == '.' for i in struct.replace('&','').replace('_','')):
        return False

    if struct[0] == '&' or struct[len(struct)-1] == '&':
        return False

    if len(struct.replace('&','').replace('_','')) <= 3:
        return False

    if '&' in struct or '_' in struct:
        if len(struct.replace('&','').replace('_',''))/(struct.count('&') + struct.count('_')) < 3:
            return False

    if any(struct.count(i) != struct.count(j) for i,j in zip(['(', '[', '<'], [')', ']', '>'])) or (struct.count('(') + struct.count('[') + struct.count('<') == 1):
        return False

    i1 = 0; i2 = 0; i3 = 0; j1 = 0; j2 = 0; j3 = 0
    for u,v in zip(range(len(struct)), range(len(struct)-1,-1,-1)):
        if struct[u] == '(':
            i1 += 1
        elif struct[u] == '[':
            i2 += 1
        elif struct[u] == '<':
            i3 += 1
        if struct[u] == ')':
            i1 -= 1
        elif struct[u] == ']':
            i2 -= 1
        elif struct[u] == '>':
            i3 -= 1
        
        if struct[v] == '(':
            j1 -= 1
        elif struct[v] == '[':
            j2 -= 1
        elif struct[v] == '<':
            j3 -= 1
        if struct[v] == ')':
            j1 += 1
        elif struct[v] == ']':
            j2 += 1
        elif struct[v] == '>':
            j3 += 1
        
        if any(i < 0 for i in [i1, i2, i3, j1, j2, j3]):
            return False

    return True

def main_argument():
    parser = argparse.ArgumentParser(description ="This tool is used to predict RNA 2D structure! :))")

    parser.add_argument('-v', '--verbose',
        action ='store_true', 
        help ='verbose mode.')

    parser.add_argument('-i',
        '--identity',
		type = float,
		default = 1,
        help ='identity threshold of the alignment. Default = 1 (no mismatching).')

    parser.add_argument('-t',
        '--threshold',
		type = float,
		default = 0.00001,
        help ='Pairing probability threshold. Default = 0.00001.')

    parser.add_argument('-s',
        '--sequence',
		help ='target sequence. Can be input from terminal or from file.')

    parser.add_argument('-m',
        '--motif',
    	help ="motif to be aligned. The input is motif filename (json format). " + 
		"Components of motif can be separated by '&' or '_'.")

    parser.add_argument('-p',
        '--pseudoknot',
        action ='store_true', 
    	help ="allowing the formation of pseudoknots")
    
    parser.add_argument('-l',
        '--lone_basepair',
        action ='store_true', 
    	help ="allowing the formation of lone basepairing.")

    parser.add_argument('-f',
        '--filter',
    	help ="filter file (json format). This file contain motifs that need to be removed for the target RNAs.")

    parser.add_argument('-o', '--outfile',
        	        #default = None,
                	action ='store',
                	help ='output file.')

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    list_nu = ['A','U','T','G','C', 'a','u','t','g','c']

    args = main_argument()

    print('Predicting seq/file:{}'.format(args.sequence))

    if args.pseudoknot:
        pkn = True
    else:
        pkn = False 

    # Read sequence or sequence file
    if all([i in list_nu for i in args.sequence]):
        seq = args.sequence
    else:
        if args.sequence in os.listdir():
            with open(args.sequence,'r') as infile:
                file = infile.read().split('\n')
                if '>' in file[0]:
                    if len(file) < 2:
                        sys.exit('This format is not supported!')
                    else:
                        if any([not(i in list_nu) for i in file[1]]):
                            sys.exit('Cannot detect sequence!')
                        else:
                            seq = file[1]
                            if 'test_' in file[0]:
                                seq_name = file[0][6:]
                            else:
                                seq_name = file[0][1:]
                else:
                    if all([i in list_nu for i in file[0]]):
                        seq = file[0]
                        seq_name = ''
                    else:
                        sys.exit('Cannot detect sequence!')
        else:
            sys.exit('filename has not been found!')
    
    # Read motif file and create filtered motif file from motif file and filter file (if any)
    if args.verbose:
        print("Aligning motifs to target sequence...")

    filtered_motif_dict = {}
    motif_file = json.load(open(args.motif))   
    for motif in motif_file:
        if args.filter:
            filter_file = json.load(open(args.filter))
            if seq_name in filter_file.keys():
                state = not(motif in [i[4:] if 'JSON' in i else i for i in filter_file[seq_name]])
            else:
                state = True
        else:
            state = True

        if state and (len(motif_file[motif]['sequence']) < len(seq)) and valid_motif(motif_file[motif]['struct2d'], pkn):
            if any(i in motif_file[motif]['struct2d'] for i in ['&', '_']):
                x = select_comp_pos(single_align_comp(seq, motif_file[motif]['sequence'], args.identity))
            else:
                x = single_align(seq, motif_file[motif]['sequence'], args.identity)
            if bool(x['pos']):
                pair = dbn_to_basepairs(motif_file[motif]['struct2d'].replace('&','')) 
                for i in range(len(x['pos'])):
                    pair2 = [(flatten(x['pos'][i])[j[0]], flatten(x['pos'][i])[j[1]]) for j in pair]
                    if any(i in motif_file[motif]['struct2d'] for i in ['&', '_']):
                        filtered_motif_dict[f'JSON{motif}_{i}'] = {'pos': x['pos'][i], 'value': x['value'][i],
                                            'pairing_pos': pair2}
                    else:
                        filtered_motif_dict[f'JSON{motif}_{i}'] = {'pos': [x['pos'][i]], 'value': x['value'][i],
                                            'pairing_pos': pair2}

    # Calculate pairing probabilities
    if args.verbose:
        print("Calculating pairing probabilities...")

    print(seq)
    pairing_dict = {}
    (propensity,ensemble_energy) = RNA.pf_fold(seq)
    for i,j in itertools.combinations(range(1,len(seq)+1),2): #Have to start from 1, otherwise it will cause error sometimes
        prob = RNA.get_pr(i, j) 
        if prob >= 10**(-10):
            if -math.log10(prob) < -math.log10(args.threshold):
                pairing_dict[(j-1,i-1)] = round(prob,int(-math.log10(args.threshold)) + 1)

    print(pairing_dict)
    if args.verbose:
        print("Creating forbidden combination of motifs...")
    list_not_combine = []
    for motif1, motif2 in itertools.combinations(filtered_motif_dict,2):
        for i in flatten(filtered_motif_dict[motif1]['pos']):
            if i in flatten(filtered_motif_dict[motif2]['pos']):
                if args.verbose:
                    print("adding {} to list...".format((motif1, motif2)))
                list_not_combine += [(motif1, motif2)]
                break
        if pkn == False:
            if not((motif1, motif2) in list_not_combine):
                for i,j in itertools.product(filtered_motif_dict[motif1]['pairing_pos'], 
                filtered_motif_dict[motif2]['pairing_pos']):
                    if ((j[1] >= i[1] and j[1] <= i[0] and j[0] >= i[0]) or 
                    (i[1] >= j[1] and i[1] <= j[0] and i[0] >= j[0])):
                        if args.verbose:
                            print("adding {} to list...".format((motif1, motif2)))
                        list_not_combine += [(motif1, motif2)]
                        break
        if len(list_not_combine) > 100000:
            sys.exit("Combinatorial explosion!")

    print(filtered_motif_dict)
    
    #sys.exit()
    # Integer Programming
    num = 1
    superpose = True
    if args.verbose:
        print("\nBegin to solve integer programming problem...\n")

    sol_val_list = []
    sol_var_list = []
    while True:
        m = Model(name='RNA2DPRED')
        x = m.binary_var_dict(filtered_motif_dict, name = '')
        y = m.binary_var_dict(pairing_dict, name = '')

        # Avoid overlapping motifs
        if args.verbose:
            if pkn == True:
                print("Adding constraints to avoid overlapping motifs...")
            else:
                print("Adding constraints to avoid overlapping motifs and pseudoknots from motifs...")
        for i in list_not_combine:
            m.add_constraint(x[i[0]] + x[i[1]] <= 1)

        if args.verbose:
            print("Adding constraints to unexpected pairings in motif regions...")
        # Avoid probabilistic basepairs in motif region:
        for motif in filtered_motif_dict:   
            for pair in pairing_dict:
                if ((pair[0] in flatten(filtered_motif_dict[motif]['pos'])) or
                (pair[1] in flatten(filtered_motif_dict[motif]['pos']))):
                    m.add_constraint(x[motif] + y[pair] <= 1)

        # Avoid pseudoknot (if selected)
        if pkn == False:
            if args.verbose:
                print("Adding constraints to avoid pseudoknots from probabilistic model...")
            #Avoid pseudoknot from probabilistic model
            for i,j in itertools.combinations(pairing_dict,2):
                if ((j[1] >= i[1] and j[1] <= i[0] and j[0] >= i[0]) or 
                (i[1] >= j[1] and i[1] <= j[0] and i[0] >= j[0])):
                    m.add_constraint(y[i] + y[j] <= 1)
            if args.verbose:
                print("Adding constraints to avoid pseudoknots between motifs and probabilistic model...")
            #Avoid pseudoknot between probabilistic model and motifs
            for motif in filtered_motif_dict:
                for i,j in itertools.product(pairing_dict, filtered_motif_dict[motif]['pairing_pos']):
                    if ((j[1] >= i[1] and j[1] <= i[0] and j[0] >= i[0]) 
                    or (i[1] >= j[1] and i[1] <= j[0] and i[0] >= j[0])):
                        m.add_constraint(y[i] + x[motif] <= 1)

        else:
            for i,j in itertools.combinations(pairing_dict,2):
                if (((j[1] in [i[1],i[0]]) or (j[0] in [i[1],i[0]])) or 
                ((i[1] in [j[1],j[0]]) or (i[0] in [j[1],j[0]]))):
                    m.add_constraint(y[i] + y[j] <= 1)

        # Avoid lone basepair (Attemp)
        if not(args.lone_basepair):
            if args.verbose:
                print("Adding constraints to avoid lone base pairings...")

            for i in pairing_dict:
                if ((not(i[0]+1 in [pair[0] for pair in pairing_dict]) and
                not(i[0]-1 in [pair[0] for pair in pairing_dict])) or 
                (not(i[1]+1 in [pair[1] for pair in pairing_dict]) and
                not(i[1]-1 in [pair[1] for pair in pairing_dict]))):
                    m.add_constraint(y[i] <= 0)
                
                else: 
                    for pair1, pair2 in itertools.combinations(pairing_dict,2):
                        if ((((pair1[0] == i[0]+1) and (pair2[0] == i[0]+1)) or
                        ((pair2[0] == i[0]+1) and (pair1[0] == i[0]+1))) or
                        (((pair1[1] == i[1]+1) and (pair2[1] == i[1]+1)) or
                        ((pair2[1] == i[0]+1) and (pair1[1] == i[1]+1)))):
                            m.add_constraint(y[pair1] - y[i] + y[pair2] >= 0)



        # Avoid base pairing too closed:
        if args.verbose:
            print("Adding constraints to avoid close base pairings..")

        for i in pairing_dict:
            if i[0] <= (i[1] + 2):
                m.add_constraint(y[i] <= 0)

        for motif in filtered_motif_dict:
            for i in filtered_motif_dict[motif]['pairing_pos']:
                if i[0] <= (i[1] + 2):
                    m.add_constraint(x[motif] <= 0)

        # Objective function
        prob_from_motifs = 0
        for motif in filtered_motif_dict:
            for pair in filtered_motif_dict[motif]['pairing_pos']:
                if pair in pairing_dict:
                    prob_from_motifs += x[motif]*pairing_dict[pair]

        if num == 1:
            n = m
            n.set_objective('max',sum([x[i]*sum([len(comp)**2 
            for comp in filtered_motif_dict[i]['pos']]) for i in filtered_motif_dict]))
            t = n.solve()

        # Function A and MEA (Prioritize MEA) (the greater number, the more priority)
        m.set_multi_objective('max', [sum([x[j]*sum([len(comp)**2 for comp in filtered_motif_dict[j]['pos']]) 
                            for j in filtered_motif_dict]),
                            sum([y[i]*pairing_dict[i] for i in pairing_dict]) + prob_from_motifs],
                            priorities = [0,1])
        

        # Get pareto set
        if len(sol_val_list) > 0:
            if superpose == True:
                # Find superposed solutions:
                m.add_constraints(sum(1-x[i] for i in filtered_motif_dict if i in pre_sol) + 
                sum(1 - y[i] for i in pairing_dict if i in pre_sol) + 
                sum(x[i] for i in filtered_motif_dict if not(i in pre_sol)) + 
                sum(y[i] for i in pairing_dict if not(i in pre_sol)) >= 1 
                for pre_sol in sol_var_list)

                m.add_constraint(sum([x[i]*sum(len(comp)**2 for comp in filtered_motif_dict[i]['pos']) for i in filtered_motif_dict]) >= sol_val_list[-1][0])
            else:
                #if  sol_val_list[-1][0] <= t.objective_value - args.threshold:
                if  sol_val_list[-1][0] <= t.objective_value and t.objective_value != 0:
                    m.add_constraint(sum([x[i]*sum(len(comp)**2 for comp in filtered_motif_dict[i]['pos']) for i in filtered_motif_dict]) >= sol_val_list[-1][0] + args.threshold)
                else:
                    print("\nNo more non-dorminated solutions!")
                    break

            #m.add_constraint(sum([y[i]*pairing_dict[i] for i in pairing_dict]) + prob_from_motifs >= sol_val_list[-1][1] + 0.005)

        sol = m.solve()

        # This code is not optimal!
        if bool(sol) and num <= 100: 
            print("The solution number {0} is: {1}".format(num, sol.multi_objective_values))
            if bool(sol_val_list):
                if (sol.multi_objective_values[1] < sol_val_list[len(sol_val_list) - 1][1] and 
                    sol.multi_objective_values[0] <= sol_val_list[len(sol_val_list) - 1][0]):
                    print("But it's dorminated!\n")
                    print("Finding on top of {0}:".format(sol.multi_objective_values[0]))
                    superpose = False
                else:
                    superpose = True
                    sol_val_list.append(sol.multi_objective_values)
                    sol_var_list.append([str(var_name)[1:] if 'JSON' in str(var_name)
                    else (int(str(var_name)[1:].split('_')[0]), int(str(var_name)[1:].split('_')[1])) 
                    for var_name in sol._var_value_map])
                    print("\nKeep finding superposed solutions at {0}:".format(sol.multi_objective_values[0]))
            else:
                superpose = True
                sol_val_list.append(sol.multi_objective_values)
                sol_var_list.append([str(var_name)[1:] if 'JSON' in str(var_name)
                else (int(str(var_name)[1:].split('_')[0]), int(str(var_name)[1:].split('_')[1])) 
                for var_name in sol._var_value_map])
                print("Keep finding superposed solutions at {0}:\n".format(sol.multi_objective_values[0]))

        else:
 #           print("\nNo more non-dorminated solutions!")
            break

        num += 1

    print("List of solutions:")
    # Show the result
    sol_pair_list = []
    out_align = ''
    sol_dict = {}
    s = 1
    for sol, score in zip(sol_var_list,sol_val_list):
        x = []
        align = ''
        for i in sol:
            if type(i) == tuple:
                x += [i] 
            else:
                x += filtered_motif_dict[i]['pairing_pos']
                align_site = ' '*len(seq)
                for u in flatten(filtered_motif_dict[i]['pos']):
                    align_site = align_site[:u] + '-' + align_site[(u + 1):]

                align += align_site + '\t' + i.split('_')[0] + '\n'
        struct2d = visual_structure(x, len(seq))
        
        print(struct2d + '\t' + str(score[0]) + '\t' + str(score[1]))
        print(align + '\n')
        out_align += f'Solution_{s}:\n' + struct2d + '\t' + str(score[0]) + '\t' + str(score[1]) + '\n' + align + '\n'
        sol_pair_list += [x]
        sol_dict[f'Solution_{s}'] = {'Score': score,
                                    'Motif':[i.split('_')[0] for i in sol_var_list[s-1] if type(i) == str],
                                    'Base_pair': x}
        s += 1

    #print("\nThe whole Pareto set is: {0} \n".format(sol_val_list))
    #print("The whole variable set is: {0} ".format(sol_var_list))
    print("The result dict is: {0} ".format(sol_dict))
    #print("\nThe whole pair set is: {0} \n".format(sol_pair_list))
    
    if args.outfile:
        with open('aligned_' + args.outfile,'w') as outfile1:
                outfile1.write(out_align)

        with open(args.outfile,'w') as outfile2:
                outfile2.write(json.dumps(sol_dict, indent= 2))
    sys.exit()
    
    
import sys
import os
import RNA
import json
import math
import itertools
from itertools import product
from docplex.mp.model import Model
from motif_alignment import *

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
    	help ="motif to be aligned. The input is motif filename (json format). " + 
		"Components of motif can be separated by '&' or '_'.")

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
            state = not(motif in [i[4:] if 'JSON' in i else i for i in filter_file[seq_name]])
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

    pairing_dict = {}
    (propensity,ensemble_energy) = RNA.pf_fold(seq)
    for i,j in itertools.combinations(range(len(seq)+1),2):
        prob = RNA.get_pr(i, j) 
        if prob != 0:
            if -math.log10(prob) < -math.log10(args.threshold):
                pairing_dict[(j-1,i-1)] = round(prob,int(-math.log10(args.threshold)) + 1)

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

        # Avoid lone basepair (under development)

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
                if  sol_val_list[-1][0] <= t.objective_value:
                    m.add_constraint(sum([x[i]*sum(len(comp)**2 for comp in filtered_motif_dict[i]['pos']) for i in filtered_motif_dict]) >= sol_val_list[-1][0] + args.threshold)
                else:
                    print("\nNo more non-dorminated solutions!")
                    break

            #m.add_constraint(sum([y[i]*pairing_dict[i] for i in pairing_dict]) + prob_from_motifs >= sol_val_list[-1][1] + 0.005)

        sol = m.solve()

        # This code is not optimal!
        if bool(sol): 
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
        struct2d = '.'*len(seq)
        for i in x:
            struct2d = struct2d[:i[1]] + '(' + struct2d[(i[1] + 1):i[0]] + ')' + struct2d[(i[0] + 1):]
        
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
    
    
import os
import sys
import json
from def_function import *
import argparse

def main_argument():
    parser = argparse.ArgumentParser(description ="This tool is used to process output of RNA2DPRED! :))")
    
    parser.add_argument('-v', '--verbose',
        action ='store_true', 
        help ='verbose mode.')
    
    parser.add_argument('-a',
        '--all_solutions',
        action ='store_false', 
        help ="Get all solutions for each benchmark. Otherwise only takes the best (default).")

    parser.add_argument('-b',
        '--benchmark',
        help ="The benchmark filename (json format).")

    parser.add_argument('-o', '--outfile',
        #default = None,
        action ='store',
        help ='output file.')
    
    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = main_argument()

    benchmark = json.load(open(args.benchmark))
    dict_result = {}
    for filename in os.listdir():
        if '.json'in filename and not('benchmark' in filename):
            file = json.load(open(filename))
            seq_name = "{0}_{1}".format(filename.split('_')[1], filename.split('_')[2])
            if seq_name.split('_')[0] in benchmark:
                dict_result[seq_name] = {}
                for i in range(len(file)):
                    z = [(t[0], t[1]) for t in file[f'Solution_{i+1}']['Base_pair']]
                    dict_result[seq_name][f'Solution_{i+1}'] = mattews_corr_coeff(compare_two_structures(
                        benchmark[seq_name.split('_')[0]][seq_name.split('_')[1]]['struct2d'], z, False))
    
    if args.all_solutions == False:
        if args.verbose:
            print('All solutions are:', dict_result, sep = '\n')
        if args.outfile:
            with open(args.outfile,'w') as outfile:
                outfile.write(json.dumps(dict_result, indent = 4))
        sys.exit()
    
    else:
        max_dict_result = {}
        for seq_name in dict_result:
            max_dict_result[seq_name] = {}
            max_dict_result[seq_name]['MCC'] = max(dict_result[seq_name][i] for i in dict_result[seq_name])
            max_dict_result[seq_name]['Solution'] = max(i for i in dict_result[seq_name] if dict_result[seq_name][i] == max_dict_result[seq_name]['MCC'])
        
        if args.verbose:
            print('The best solutions are:', max_dict_result, sep = '\n')
        if args.outfile:
            with open(args.outfile,'w') as outfile:
                outfile.write(json.dumps(max_dict_result, indent = 4))
        sys.exit()              
    

from optparse import Option
import timeit
import sys
import os
from itertools import product
import threading
import json
import math
import queue
import argparse

def flatten(lol):
	result = []
	for l in lol:
		result += l

	return result

def visual_alignment(s, pos_list, components = False):
	string = ''
	sol_name_list = [f"sol{i+1}" for i in range(len(pos_list))]
	if components == False:
		for sol_name, pos in zip(sol_name_list, pos_list):
			align_string = ''.join('-' if i in pos else ' ' for i in range(len(s)))
			string += f"{sol_name}\n\t{s}\n\t"
			string += align_string + '\n'

	else:
		flatten_post_list = [flatten(i) for i in pos_list]
		for sol_name, pos in zip(sol_name_list, flatten_post_list):
			align_string = ''.join('-' if i in pos else ' ' for i in range(len(s)))
			string += f"{sol_name}\n\t{s}\n\t"
			string += align_string + '\n'
	
	return string



def single_align(s,c, identity = 1, shift = 0):
	c = c.replace('&&','&').replace('_','&')
	pos_dict = {'pos': [], 'value': []}
	if identity == 1:
		for i in range(len(s) - len(c) + 1):
			if s[i:(i+len(c))].upper() == c.upper():
				pos_dict['pos'] += [[i+j + shift for j in range(len(c))]]
				pos_dict['value'] += [1]
	else:
		mm_allow = (1 - identity)*len(c)
		for i in range(len(s) - len(c) + 1):
			mm = 0
			flag = 0
			for j in range(len(c)):
				if s[i:(i+len(c))][j].upper() != c[j].upper():
					mm += 1
					if mm > mm_allow:
						flag = 1
						break
					#The subtitute code below is shorter but slower:
					#mm = sum(s[i:(i+len(c))][j] != c[j] for j in range(len(c)))
			#if (len(c)-mm)/len(c) >= identity:
			if flag == 0:
				pos_dict['pos'] += [[i + j + shift for j in range(len(c))]]
				pos_dict['value'] += [(len(c)-mm)/len(c)]

	return pos_dict

def single_align_comp(s,c, identity = 1, shift = 0):
	pos_dict = {'pos': [], 'value': []}
	if '&' in c or '_' in c:
		c = c.replace('&&','&').replace('_','&')
		list_c = c.split('&')
		for component in list_c:
			x = single_align(s,component, identity, shift)
			pos_dict['pos'] += [x['pos']]
			pos_dict['value'] += [x['value']]

		return  pos_dict
	else:
	        return single_align(s,c, identity, shift)


def select_comp_pos(post_dict, thread = 0, threads = 1):
        selected_post_dict = {'pos': [], 'value': []}
        y = post_dict['pos'][0]
        z = post_dict['value'][0]
        l_post = [y[int(thread*len(y)/threads):int((thread + 1)*len(y)/threads)]] + list(post_dict['pos'][1:])
        l_val = [z[int(thread*len(z)/threads):int((thread + 1)*len(z)/threads)]] + list(post_dict['value'][1:])

        for u,v in zip(product(*l_post), product(*l_val)):
                if all([u[j][len(u[j]) - 1] < u[j+1][0] for j in range(len(u)-1)]):
                        selected_post_dict['pos'] += [u]
                        selected_post_dict['value'] += [sum(len(u[i])*v[i] for i in range(len(v)))/len(flatten(u))]

        return selected_post_dict

def argument():
	parser = argparse.ArgumentParser(description ="This tool is used to align motif to the target sequence. It's mainly for the RNA2DPRED tool :))")

	parser.add_argument('-v', '--verbose',
                    action ='store_true', help ='verbose mode.')

	parser.add_argument('-i',
		        '--identity',
			type = float,
			default = 1,
                	help ='identity threshold. Default = 1 (no mismatching).')

	parser.add_argument('-t',
        	        '--threads',
			type = int,
			default = 1,
	                help ='number of threads.')

	parser.add_argument('-s',
        	        '--sequence',
			 help ='target sequence. Can be input from terminal or rom file.')

	parser.add_argument('-m',
        	        '--motif',
                	help ="motif to be aligned. Can be input from terminal or from motif file. For motif which contains " + 
			"many components, components are separated by '&' or '_' (note that '&' cannot be used for input from terminal).")

	parser.add_argument('-o', '--outfile',
        	        #default = None,
                	action ='store',
                	help ='output file.')

	parser.add_argument('-f', 
					'--format',
					choices= ['0', '1'],
					default = '0',
					help="format of the output. Either dict of solutions, each solution contains 'pos' and 'value' field (0, default)" + 
					" or dict with 'pos' and 'value' field, each field is a list of solution (1).")

	args = parser.parse_args()

	return args

if __name__ == "__main__":
	list_nu = ['A','U','T','G','C', 'a','u','t','g','c']

	args = argument()

	# Read input sequence
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
				else:
					if all([i in list_nu for i in file[0]]):
						seq = file[0]
					else:
						sys.exit('Cannot detect sequence!')
		else:
			sys.exit('filename has not been found!')

	# Read motif sequence
	if all([i in list_nu + ['_'] for i in args.motif]):
		motif = args.motif
		motif = motif.replace('_','&')
	else:
		if args.motif in os.listdir():
			with open(args.motif,'r') as infile:
				file = infile.read().split('\n')
				if any([not(i in list_nu + ['_', '&']) for i in file[0]]):
					sys.exit('Cannot detect sequence!')
				else:
					motif = file[0].replace('_','&')
		else:
			sys.exit('filename has not been found!')

	que = queue.Queue()

	list_threads = []
	if '&' in motif:
		dict_motif = single_align_comp(seq,motif,args.identity)

	for num_thread in range(args.threads):
		if '&' in motif:
			list_threads.append(threading.Thread(target = lambda q, x, y, z: q.put(select_comp_pos(x,y,z)),
			args = (que, dict_motif, num_thread, args.threads)))
		else:
			if num_thread != args.threads - 1:
				list_threads.append(threading.Thread(target = lambda q, x, y, z, s: q.put(single_align(x,y,z,s)),
				args = (que,seq[int(len(seq)*num_thread/args.threads):(int(len(seq)*(num_thread+1)/args.threads) + len(motif) - 1)],motif,args.identity,int(len(seq)*num_thread/args.threads))))
			else:
				list_threads.append(threading.Thread(target = lambda q, x, y, z, s: q.put(single_align(x,y,z,s)),
        			args = (que,seq[int(len(seq)*num_thread/args.threads):int(len(seq)*(num_thread+1)/args.threads)],motif,args.identity,int(len(seq)*num_thread/args.threads))))

	start = timeit.default_timer()

	for key in list_threads:
		key.start()

	for key in list_threads:
		key.join()

	stop = timeit.default_timer()

	class Result(object):
		def __init__(self, pos, value):
			self.pos = pos
			self.value = value

		def first_pos(self):
			l = []
			for i in self.pos:
				l += [i[0]]

			return l
		def last_pos(self):
			l = []
			for i in self.pos:
				l += [i[len(i)-1]]

			return l

	result = {'pos': [], 'value': []}
	while not que.empty():
		subresult = que.get()
		result['pos'] += subresult['pos']
		result['value'] += subresult['value']

	if not bool(result['pos']):
		sys.exit('No solutions found!')

	if args.format == '0':
		result2 = {}
		for i in range(len(result['pos'])):
			result2[f'sol{i+1}'] = {'pos': result['pos'][i], 'value': result['value'][i]}

	if args.verbose:
		if args.format == '1':
			print(result, end = '\n\n')
		else: 
			print(result2, end = '\n\n')
			
		print("Proccessing time: ", stop - start, sep = '')
		print("Number of threads used: ", args.threads , sep = '', end = '\n\n')

	
	print(visual_alignment(seq, result['pos'], '&' in motif))

	if args.outfile:
		with open(args.outfile, 'w') as out:
			if args.format == '1':
				out.write(json.dumps(result, indent = 2))	
			else:
				out.write(json.dumps(result2, indent = 2))		
		
		with open('align_' + args.outfile, 'w') as out:
			out.write(visual_alignment(seq, result['pos'], '&' in motif))


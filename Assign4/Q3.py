from __future__ import division
import os
import glob
import numpy as np
import csv

SPIDER_PATH = 'SPIDER2/'
TEST_RAPTOR_PATH = 'RAPTORX/4ympA/303301.ss3_simp.txt'
RAPTOR_PATH = 'RAPTORX/'
TEST_DSSP_PATH = 'DSSP/4ymp_DSSP.txt'
DSSP_PATH = 'DSSP/'
truth1 = "HHHCCCEEE"
test1 = "HHHHCCEEE"
OUTFILE = 'Q3.csv'


def get_Q3(truth, test):
	'''
	given two secondary structures as strings
	returns the Q3 score
	'''
	#if (len(truth) != len(test)):
		#print("Mismatched Lengths")
	#else:
	tot = 0
	score = 0
	for i in range(len(truth)):
		if (truth[i] != ' '):
			tot = tot + 1
		if (truth[i] == test[i]):
			score = score + 1
	Q3 = score/tot
	return Q3

def read_DSSP(dssp_file, chain):
	'''
	given a DSSP file
	returns the ground truth SS3
	'''
	f = open(dssp_file, "r")
	SS = ""
	for i in range(28):
		f.readline()
	for line in f:
		fields = line.split()
		if (fields[2] == chain):
			s = fields[4]
			if (s == 'E' or s == 'B'):
				state = 'E'
			elif (s == 'G' or s == 'H'):
				state = 'H'
			elif (s.isalpha() or s == ' '):
				state = 'C'
			else:
				state = ' '
			SS = SS + state
	return SS

def read_RAPTORX(raptorx_file):
	'''
	given ss3_simp file from raptor
	takes the SS3 FASTA string and returns
	'''
	f = open(raptorx_file, "r")
	f.readline()
	f.readline()
	SS = f.readline()
	return SS

def read_SPIDER2(spider2_file):
	'''
	given spider2 file
	produces SS3 FASTA string
	'''
	SS = ""
	f = open(spider2_file, "r")
	f.readline()
	for line in f:
		fields = line.split("	")
		SS = SS + fields[2]
	return SS


def main():
	spider_dict = {}
	raptor_dict = {}
	dssp_dict = {}
	pro_list = []

	for file in os.listdir(SPIDER_PATH):
		#print(SPIDER_PATH + file)
		SS = read_SPIDER2(SPIDER_PATH + file)
		fields = file.split("_")
		spider_dict[fields[0]] = SS 
		pro_list.append(fields[0])

	for folder in os.listdir(RAPTOR_PATH):
		if (folder == '.DS_Store'):
			continue
		path = RAPTOR_PATH + folder + '/'
		os.chdir(path)
		file = glob.glob('*ss3_simp.txt')
		os.chdir('..')
		os.chdir('..')
		SS = read_RAPTORX(path + file[0])
		SS = SS.strip('\n')
		raptor_dict[folder] = SS

	for file in os.listdir(DSSP_PATH):
		code = file.split('_')[0]
		if (code == '5a7d'):
			SS = read_DSSP(DSSP_PATH + file, 'B')
			dssp_dict[code + 'B'] = SS
			SS = read_DSSP(DSSP_PATH + file, 'L')
			dssp_dict[code + 'L'] = SS
		elif (code == '5j4a'):
			SS = read_DSSP(DSSP_PATH + file, 'A')
			dssp_dict[code + 'A'] = SS
			SS = read_DSSP(DSSP_PATH + file, 'B')
			dssp_dict[code + 'B'] = SS
		elif (code == '5j5v'):
			SS = read_DSSP(DSSP_PATH + file, 'A')
			dssp_dict[code + 'A'] = SS
			SS = read_DSSP(DSSP_PATH + file, 'B')
			dssp_dict[code + 'B'] = SS
			SS = read_DSSP(DSSP_PATH + file, 'C')
			dssp_dict[code + 'C'] = SS
		else:
			SS = read_DSSP(DSSP_PATH + file, 'A')
			dssp_dict[code + 'A'] = SS

	scores = np.zeros([2,15])
	#[0][1:16] = pro_list

	for i, pro in enumerate(pro_list):
		check1 = dssp_dict[pro]
		check2 = spider_dict[pro]
		check3 = raptor_dict[pro]
		scores[0][i] = get_Q3(dssp_dict[pro], spider_dict[pro])
		scores[1][i] = get_Q3(dssp_dict[pro], raptor_dict[pro])

	np.savetxt(OUTFILE, scores, delimiter = ',',fmt = '%1.2f')

	with open('NAMES.csv', 'wb') as f:
		wr = csv.writer(f, quoting=csv.QUOTE_ALL)
		wr.writerow(pro_list)


if __name__ == '__main__':
	main()
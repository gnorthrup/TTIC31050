from __future__ import division
import os
import glob


SPIDER_PATH = 'SPIDER2/'
TEST_RAPTOR_PATH = 'RAPTORX/4ympA/303301.ss3_simp.txt'
RAPTOR_PATH = 'RAPTORX/'
truth1 = "HHHCCCEEE"
test1 = "HHHHCCEEE"


def get_Q3(truth, test):
	'''
	given two secondary structures as strings
	returns the Q3 score
	'''
	if (len(truth) != len(test)):
		print("Mismatched Lengths")
	else:
		tot = 0
		score = 0
		for i in range(len(truth)):
			if (truth[i] != ' '):
				tot = tot + 1
			if (truth[i] == test[i]):
				score = score + 1
		Q3 = score/tot
		return Q3

# def read_DSSP(dssp_file):

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
	Q3 = get_Q3(truth1, test1)
	print(Q3)
	spider_dict = {}
	raptor_dict = {}
	pro_list = []
	for file in os.listdir(SPIDER_PATH):
		print(SPIDER_PATH + file)
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

if __name__ == '__main__':
	main()
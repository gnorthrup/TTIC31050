#This file produces the code for the first assignment of TTIC 301
#Graham Northrup

import numpy as np

GAP_PEN = -1
MATCH_SCORE = 1
NONMATCH_SCORE = -1
OPEN_GAP = 0 #-5,-10
KEEP_GAP = -1 #0
GAP = False
AFFINE = True

OXBENCH1 = 'ATGAAATTTACCGTAGAACGTGAGCATTTATTAAAACCGCTACAACAGGTGAGCGGTCCGTTAGGTGGTCGTCCTACGCTACCGATTCTCGGTAATCTGCTGTTACAGGTTGCTGACGGTACGTTGTCGCTGACCGGTACTGATCTCGAGATGGAAATGGTGGCACGTGTTGCGCTGGTTCAGCCACACGAGCCAGGAGCGACGACCGTTCCGGCGCGCAAATTCTTTGATATCTGCCGTGGTCTGCCTGAAGGCGCGGAAATTGCCGTGCAGCTGGAAGGTGAACGGATGCTGGTACGCTCCGGGCGTAGCCGTTTTTCGCTGTCTACCCTGCCAGCGGCGGATTTCCCGAACCTCGATGACTGGCAG'
OXBENCH2 = 'AGTGAAGTCGAATTTACCCTGCCGCAGGCAACGATGAAGCGTCTGATTGAAGCGACCCAGTTTTCTATGGCGCATCAGGACGTTCGCTATTACTTAAATGGTATGCTGTTTGAAACCGAAGGTGAAGAACTGCGCACCGTGGCAACCGACGGCCACCGTCTGGCGGTCTGTTCAATGCCAATTGGTCAATCTTTGCCAAGCCATTCGGTGATCGTACCGCGTAAAGGCGTGATTGAACTGATGCGTATGCTCGACGGCGGCGACAATCCGCTGCGCGTACAGATTGGCAGCAACAACATTCGCGCCCACGTTGGCGACTTTATCTTCACCTCCAAACTGGTGGATGGTCGCTTCCCGGATTATCGCCGCGTT'

#>2pola-1-AS.dna|gb|U00096.2|.revcom
OXBENCH1GT = '------ATGAAATTTACCGTAGAACGTGAGCATTTATTAAAACCGCTACAACAGGTGAGCGGTCCGTTAGGTGGTCGTCCTACGCTACCGATTCTCGGTAATCTGCTGTTACAGGTTGCTGACGGTACGTTGTCGCTGACCGGTACTGATCTCGAGATGGAAATGGTGGCACGTGTTGCGCTGGTTCAGCCACACGAGCCAGGAGCGACGACCGTTCCGGCGCGCAAATTCTTTGATATCTGCCGTGGTCTGCCT---GAAGGCGCGGAAATTGCCGTGCAGCTGGAAGGTGAACGGATGCTGGTACGCTCCGGGCGTAGCCGTTTTTCGCTGTCTACCCTGCCAGCGGCGGATTTCCCG------AACCTCGATGACTGGCAG'
#>2pola-2-AS.dna|gb|U00096.2|.revcom
OXBENCH2GT = 'AGTGAAGTCGAATTTACCCTGCCGCAGGCAACGATGAAGCGTCTGATTGAAGCGACCCAGTTTTCTATGGCGCATCAGGACGTTCGCTATTACTTAAATGGTATGCTGTTTGAAACCGAAGGTGAAGAACTGCGCACCGTGGCAACCGACGGCCACCGTCTGGCGGTCTGTTCAATGCCAATTGGTCAATCTTTGCCAAGCCATTCGGTGATCGTACCGCGTAAAGGCGTGATTGAACTGATGCGTATGCTCGACGGCGGCGACAATCCGCTGCGCGTACAGATTGGCAGCAACAACATTCGCGCCCACGTTGGCGACTTTATCTTCACCTCCAAACTGGTGGAT---GGTCGCTTCCCGGATTATCGCCGCGTT---------'

A = 'ATGAAATTT'
B = 'TGAAGTCGAATTT'
AGT = '----ATGAAATTT'
BGT = 'TGAAGTCGAATTT'

def get_SW_mat(seq1, seq2):
	'''
	given two sequences as strings
	populates and returns a Smith-Waterman matrix of alignment scores
	'''
	N = len(seq1)
	M = len(seq2)
	score = np.zeros((N+1,M+1))
	for i in range(1,N+1):
		for j in range(1, M+1):
			if (seq1[i-1] == seq2[j-1]):
				align = MATCH_SCORE
			else:
				align = NONMATCH_SCORE
			match = score[i-1][j-1] + align
			gap1 = score[i][j-1] + GAP_PEN
			gap2 = score[i-1][j] + GAP_PEN
			score[i][j] = max(match, gap1, gap2, 0)
	return score


def walk_back_SW(mat, seq1, seq2):
	'''
	given the sequences as strings and the scoring matrix as
	numpy array
	returns the optimal local alignment
	'''
	N = len(seq1)
	M = len(seq2)
	local_pairs = []
	best_i = 0
	best_j = 0
	best_score = 0
	for i in range(1,N+1):
		for j in range(1, M+1):
			if (mat[i][j] > best_score):
				best_i = i
				best_j = j
				best_score = mat[i][j]
	#print(best_score)
	
	N_seq = [seq1[best_i - 1]]
	M_seq = [seq2[best_j - 1]]
	local_pairs.append((best_i, best_j))
	score = best_score
	row = best_i
	column = best_j
	while (score > 0):
		match = mat[row - 1][column - 1]
		gap_l = mat[row][column - 1]
		gap_u = mat[row - 1][column]
		step = max(match, gap_l, gap_u)
		check = min(match, gap_l, gap_u)
		if (match == 0):
			score = 0
			continue
		elif (step == match):
			N_seq.append(seq1[row-2])
			M_seq.append(seq2[column-2])
			local_pairs.append((row-1,column-1))
			score = step
			row = row - 1
			column = column - 1
		elif (step == gap_l):
			N_seq.append('-')
			M_seq.append(seq2[column-2])
			score = step
			column = column - 1
		elif (step == gap_u):
			N_seq.append(seq1[row-2])
			M_seq.append('-')
			score = step
			row = row - 1
	N_seq.reverse()
	M_seq.reverse()
	N_seq = "".join(N_seq)
	M_seq = "".join(M_seq)
	return (N_seq, M_seq, local_pairs)

def get_NW_mat(seq1, seq2):
	'''
	given two sequences as strings,
	produces global alignment matrix
	'''
	N = len(seq1)
	M = len(seq2)
	score = np.zeros((N+1,M+1))
	for i in range(0,N+1):
		score[i][0] = GAP_PEN*i
	for j in range(0,M+1):
		score[0][j] = GAP_PEN*j
	for i in range(1,N+1):
		for j in range(1, M+1):
			if (seq1[i-1] == seq2[j-1]):
				align = MATCH_SCORE
			else:
				align = NONMATCH_SCORE
			match = score[i-1][j-1] + align
			gap1 = score[i][j-1] + GAP_PEN
			gap2 = score[i-1][j] + GAP_PEN
			score[i][j] = max(match, gap1, gap2)
	return score

def get_NW_mat_AFFINE(seq1, seq2):
	'''
	given two sequences as strings,
	produces global alignment matrix that uses an affine gap penalty
	'''
	N = len(seq1)
	M = len(seq2)
	score = np.zeros((N+1,M+1))
	scoreS = np.zeros((N+1,M+1))
	scoreT = np.zeros((N+1,M+1))
	for i in range(0,N+1):
		score[i][0] = KEEP_GAP*i + OPEN_GAP
		scoreS[i][0] = KEEP_GAP*i + OPEN_GAP
		scoreT[i][0] = KEEP_GAP*i + OPEN_GAP
	for j in range(0,M+1):
		score[0][j] = KEEP_GAP*j + OPEN_GAP
		scoreS[0][j] = KEEP_GAP*j + OPEN_GAP
		scoreT[0][j] = KEEP_GAP*j + OPEN_GAP
	for i in range(1,N+1):
		for j in range(1, M+1):
			if (seq1[i-1] == seq2[j-1]):
				align = MATCH_SCORE
			else:
				align = NONMATCH_SCORE
			a = score[i-1][j-1] + align
			b = scoreS[i-1][j-1] + align
			c = scoreT[i-1][j-1] + align
			score[i][j] = max(a,b,c)
			d = score[i][j-1] - KEEP_GAP - OPEN_GAP
			e = scoreS[i][j-1] - KEEP_GAP
			f = scoreT[i][j-1] - KEEP_GAP - OPEN_GAP
			scoreS[i][j] = max(d,e,f)
			x = score[i-1][j] - KEEP_GAP - OPEN_GAP
			y = scoreS[i-1][j] - KEEP_GAP - OPEN_GAP
			z = scoreT[i-1][j] - KEEP_GAP 
			scoreT[i][j] = max(x,y,z)
	return score

 
def walk_back_NW(mat, seq1, seq2):
	'''
	given the global alignment matrix, and two sequences
	returns alignment
	'''
	N = len(seq1)
	M = len(seq2)
	N_seq = [seq1[N- 1]]
	M_seq = [seq2[M - 1]]
	row = N
	column = M
	while (row > 1 or column > 1):
		if (row == 0):
			match = -1
			gap_l = 1
			gap_u = -1
		elif(column == 0):
			match = -1
			gap_l = -1
			gap_u = 1
		else:
			match = mat[row - 1][column - 1]
			gap_l = mat[row][column - 1]
			gap_u = mat[row - 1][column]
		step = max(match, gap_l, gap_u)
		if (step == match):
			N_seq.append(seq1[row-2])
			M_seq.append(seq2[column-2])
			score = step
			row = row - 1
			column = column - 1
		elif (step == gap_l):
			N_seq.append('-')
			M_seq.append(seq2[column-2])
			score = step
			column = column - 1
		elif (step == gap_u):
			N_seq.append(seq1[row-2])
			M_seq.append('-')
			score = step
			row = row - 1
	N_seq.reverse()
	M_seq.reverse()
	N_seq = "".join(N_seq)
	M_seq = "".join(M_seq)
	return (N_seq, M_seq)

def check_alignment_acc(test1, test2, gt1, gt2):
	gt_pairs = []
	gt1_index = 0
	gt2_index = 0
	test1_index = 0
	test2_index = 0
	num_correct = 0

	for i in range(0, len(gt1) - 1):
		if (gt1[i] != '-'):
			gt1_index += 1
		if (gt2[i] != '-'):
			gt2_index += 1
		if (gt1[i] != '-' and gt2[i] != '-'):
			gt_pairs.append((gt1_index, gt2_index))

	for i in range(0, len(test1) - 1):
		if (test1[i] != '-'):
			gt1_index += 1
		if (test2[i] != '-'):
			gt2_index += 1
		if (test1[i] != '-' and test2[i] != '-'):
			if ((test1_index, test2_index) in gt_pairs):
				num_correct += 1
	score = num_correct/len(gt_pairs)
	return score

def check_local_align(local_pairs, gt1, gt2):
	gt_pairs = []
	gt1_index = 0
	gt2_index = 0
	test1_index = 0
	test2_index = 0
	num_correct = 0

	for i in range(0, len(gt1) - 1):
		if (gt1[i] != '-'):
			gt1_index += 1
		if (gt2[i] != '-'):
			gt2_index += 1
		if (gt1[i] != '-' and gt2[i] != '-'):
			gt_pairs.append((gt1_index, gt2_index))
	for j in range(0, len(local_pairs) - 1):
		if (local_pairs[j] in gt_pairs):
			num_correct += 1
	score = num_correct/len(gt_pairs)
	return score


def main():
	local_ox_mat = get_SW_mat(OXBENCH1, OXBENCH2)
	(my_local1, my_local2, pairs) = walk_back_SW(local_ox_mat, OXBENCH1, OXBENCH2)
	global_ox_mat = get_NW_mat(OXBENCH1, OXBENCH2)
	(my_global1, my_global2) = walk_back_NW(global_ox_mat, OXBENCH1, OXBENCH2)
	local_score = check_local_align(pairs, OXBENCH1GT, OXBENCH2GT)
	global_score = check_alignment_acc(my_global1, my_global2, OXBENCH1GT, OXBENCH2GT)
	print("Local Alignment Score:")
	print(local_score)
	print("Global Alignment Score:")
	print(global_score)


	affine_ox_mat = get_NW_mat_AFFINE(OXBENCH1, OXBENCH2)
	(affine_global1, affine_global2) = walk_back_NW(affine_ox_mat, OXBENCH1, OXBENCH2)
	affine_score = check_alignment_acc(affine_global1, affine_global2, OXBENCH1GT, OXBENCH2GT)
	print("Affine Gap Score:")
	print(affine_score)


if __name__ == '__main__':
	main()

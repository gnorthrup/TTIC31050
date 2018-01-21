#This file produces the code for the first assignment of TTIC 301
#Graham Northrup

import numpy as np

GAP_PEN = -2
MATCH_SCORE = 1
NONMATCH_SCORE = -1
OPEN_GAP = -1
KEEP_GAP = -1
GAP = False

A = 'AAAGTCGGTC'
B = 'GTCAGTC'

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

def walk_back_NW(mat, seq1, seq2):
	'''
	given the global alignment matrix, and two sequences
	returns alignment
	'''
	N = len(seq1)
	M = len(seq2)
	N_seq = [seq1[N- 1]]
	M_seq = [seq2[M - 1]]
	row = N - 1
	column = M - 1
	while (row > 0 or column > 0):
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


def main():
	score = get_SW_mat(A,B)
	print(score)
	(A_loc, B_loc) = walk_back_SW(score, A, B)
	print(A_loc)
	print(B_loc)
	mat = get_NW_mat(A,B)
	print(mat)
	(A_glob, B_glob) = walk_back_NW(mat, A, B)
	print(A_glob)
	print(B_glob)

if __name__ == '__main__':
	main()

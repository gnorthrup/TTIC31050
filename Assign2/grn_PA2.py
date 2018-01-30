import numpy as np
import math

SUP_MAT_FILENAME = "superfamily_matrix.csv"
TWI_MAT_FILENAME = "twilight_matrix.csv"

sup1 = "VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKF"
sup2 = "VLSSQNKKAIEELGNLIKANAEAWGADALARLFELHPQTKTYF"
sup3 = "MLTQKTKDIVKATAPVLAEHGYDIIKCFYQRMFEAHPELKNVF"
sup4 = "SFSEEQEALVLKSWAILKKDSANIALRFFLKIFEVAPSASQMF"
sup5 = "AFTACEKQTIGKIAQVLAKSPEAYGAECLARLFVTHPGSKSYF"
sup6 = "MLDAQTIATVKATIPLLVETGPKLTAHFYDRMFTHNPELKEIF"
sup7 = "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYF"
sup8 = "TLTDGDKKAINKIWPKIYKEYEQYSLNILLRFLKCFPQAQASF"
sup9 = "SLSDKDKAAVRALWSKIGKSSDAIGNDALSRMIVVYPQTKIYF"
sup10 = "SLSAAEADLAGKSWAPVFANKNANGLDFLVALFEKFPDSANFF"
sup11 = "ALTESQAALVKSSWEEFNANIPKHTHRFFILVLEIAPAAKDLF"
sup12 = "PLSAAEKTKIRSAWAPVYSTYETSGVDILVKFFTSTPAAQEFF"

twi1 = "EDNCIAEDYGKCTWGGTKCCRGRPCRCSMIGTNCECTPRLIME"
twi2 = "VKDGYIVDDVNCTYFCGRNAYCNEECTKLKGESGYCQWASPYG"
twi3 = "MSSEHRCIDTNVPENAACYRYLDGTEEWRCLLYFKEDAGKCVP"
twi4 = "NMTCKDKNGGCAPEAECKMNDKNEIVCKCTKEGSEPLFEGVFC"
twi5 = "GNTCGGETCSAAQVCLKGKCVCNEVHCRIRCKYGLKKDENGCE"
twi6 = "KRPWKCCDEAVCTRSIPPICTCMDEVFECPKTCKSCGPSMGDP"
twi7 = "YSRCQLQGFNCVVRSYGLPTIPCCRGLTCRSYFPGSTYGRCQR"
twi8 = "VEPVDPCFRANCEYQCQPLDQTSYLCVCAEGFAPIPHEPHRCQ"
twi9 = "IRFGMGKVPCPDGEVGYTCDCGEKICLYGQSCNDGQCSGDPKP"
twi10 = "TKPGSCPIILIRCAMLNPPNRCLKDTDCPGIKKCCEGSCGMAC"
twi11 = "LLACLFGNGRCSSNRDCCELTPVCKRGSCVSSGPGLVGGILGG"
twi12 = "KICRRRSAGFKGPCMSNKNCAQVCQQEGWGGGNCDGPFRRCKC"
twi13 = "SALAEGQSCGVYTERCAQGLRCLPRQDEEKPLHALLHGRGVCL"
twi14 = "CVRLHESCLGQQVPCCDPCATCYCRFFNAFCYCRKLGTAMNPC"
twi15 = "DKLIGSCVWGAVNYTSDCNGECKRRGYKGGHCGSFANVNCWCE"

AA_ORDER = {'A':0,
			'R':1,
			'N':2,
			'D':3,
			'C':4,
			'Q':5,
			'E':6,
			'G':7,
			'H':8,
			'I':9,
			'L':10,
			'K':11,
			'M':12,
			'F':13,
			'P':14,
			'S':15,
			'T':16,
			'W':17,
			'Y':18,
			'V':19}

sl = [sup1,sup2,sup3,sup4,sup5,sup6,sup7,sup8,sup9,sup10,sup11,sup12]
bl = 43
tl = [twi1,twi2,twi3,twi4,twi5,twi6,twi7,twi8,twi9,twi10,twi11,twi12,twi13,twi14,twi15]
tbl = 43
def make_mat(seq_list, block_len):
	'''
	Given a block of sequences, makes an alignment matrix
	'''
	counts = np.ones((20,20))

	for (i,seq) in enumerate(seq_list):
		if (i != len(seq_list)):
			for seq2 in seq_list[i+1:]:
				for j in range(block_len):
					counts[AA_ORDER[seq[j]]][AA_ORDER[seq2[j]]] += 1
					counts[AA_ORDER[seq2[j]]][AA_ORDER[seq[j]]] += 1


	mut_prob = counts/sum(sum(counts))
	abundance = sum(counts)/sum(sum(counts))

	out_mat = np.zeros((20,20))

	for i in range(20):
		for j in range(20):
			v = mut_prob[i][j]/(abundance[i]*abundance[j])
			out_mat[i][j] = "{0:.2f}".format(2*math.log(v,2))
	print(out_mat)

	return out_mat

def main():
	sup_mat = make_mat(sl, bl)
	np.savetxt(SUP_MAT_FILENAME, sup_mat, delimiter=",", fmt = '%1.2f')
	twi_mat = make_mat(tl, tbl)
	np.savetxt(TWI_MAT_FILENAME, twi_mat, delimiter=",", fmt = '%1.2f')


if __name__ == '__main__':
	main()




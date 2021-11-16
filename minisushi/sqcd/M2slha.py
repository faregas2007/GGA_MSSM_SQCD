# M2slha.py
#To convert the input from Micheal format into slha format.
#To include the coversion from At to Ab in diagrams with squark Higgs coupling

import sys
import os.path
import numpy as np    

name = sys.argv[1]
output = name+".in"

def read_to_slha(name, output):
	if(name[1] == 't'):
		squark = 1
	elif(name[1] == 'b'):
		squark = 2

	if squark == 1:
		minputs = {'M_A':{}, 'M_Q':{}, 'M_G':{}, 'TG(BETA)':{}, 'MU':{}, 'QQ':{}, 'M_Q(QQ)':{}, 'A_t':{}, 'M_S1':{}, 'M_S2':{}, 'G_Q^A':{}, 'ALPHAS(QS)':{}, 'SIN(THETA)':{}, 'COS(THETA)':{}, 'GAQQ(1,1)':{}, 'GAQQ(1,2)':{}, 'GAQQ(2,1)':{}, 'GAQQ(2,2)':{}}
	elif squark == 2:
		minputs = {'M_A':{}, 'M_Q':{}, 'M_G':{}, 'TG(BETA)':{}, 'MU':{}, 'QQ':{}, 'M_Q(QQ)':{}, 'A_b':{}, 'M_S1':{}, 'M_S2':{}, 'G_Q^A':{}, 'ALPHAS(QS)':{}, 'SIN(THETA)':{}, 'COS(THETA)':{}, 'GAQQ(1,1)':{}, 'GAQQ(1,2)':{}, 'GAQQ(2,1)':{}, 'GAQQ(2,2)':{}}

	lines = []
	if not os.path.isfile(name):
		print ("file does not exist")
	else:
		with open(name) as fin:
			for line in fin:
				lines.append(line.split())

	for index in range(len(lines)-1):
		for key, value in minputs.items():
			if lines[index][0] == key:
				minputs[key] = lines[index][2]

	mA = float(minputs.get('M_A'))
	mQ = float(minputs.get('M_Q'))
	mG = float(minputs.get('M_G'))
	tanb = float(minputs.get('TG(BETA)'))
	mu = float(minputs.get('MU'))
	mQQ = float(minputs.get('QQ'))
	mQQQ = float(minputs.get('M_Q(QQ)'))
	if squark == 1:
		Aq = float(minputs.get('A_t'))
	elif squark == 2:
		Aq = float(minputs.get('A_b'))
	msq1 = float(minputs.get('M_S1'))
	msq2 = float(minputs.get('M_S2'))
	gqA = float(minputs.get('G_Q^A'))
	alphaqs = float(minputs.get('ALPHAS(QS)'))
	sintheta = float(minputs.get('SIN(THETA)'))
	costheta = float(minputs.get('COS(THETA)'))
	gq11 = float(minputs.get('GAQQ(1,1)'))
	gq12 = float(minputs.get('GAQQ(1,2)'))
	gq21 = float(minputs.get('GAQQ(2,1)'))

	#corrections At to Ab in diagram with squark Higgs couplings	
	factorb = -(1. + 1./tanb**2)
	
	if squark == 1:
		Aqtree = tanb*(2*gq12 - mQQQ*mu)/mQQQ
		totalcorrfactor = 1
	elif squark == 2:
		Aqtree = (2*gq12 - mQQQ*mu)/(mQQQ*tanb)
		totalcorrfactor = (mu + Aqtree*tanb)/((mu*tanb + Aqtree)*tanb)
	
	if(squark == 1):
		text="""Block QCD
	1	%(alphaqs).12f	#alphaqs
Block REAL
	1	%(mu).12f	#mu
	2	%(Aqtree).12f	#Aq
	3	%(beta).12f	#beta
	4	0	#squark mixing
	5	%(sintheta).12f	#sintheta
	6	%(costheta).12f	#costheta
	7	%(totalcorr).12f	#correction_to _squark_Higgs_coupling
	8	%(gq12).12f	#gq12_coupling
Block INTEGER
	1	1	#top=1_bottom=2
Block MASS
	1	%(mA).12f	#mA
	2	%(mQ).12f	#mQ
	3	%(mQQQ).12f	#mQQQ
	4	%(mG).12f	#mG
	5	%(msq1).12f	#msq1
	6	%(msq2).12f	#msq2
Block VEGAS
	1	30	#iter
	2	10000	#startdiv
	3	3000000	#enddiv
	4	10000	#startfin
	5	3000000	#endfin
	6	0.4	#alpha
	7	1000	#increase
	8	Cuhre	#integrator
Block IBPPARA
	1	0.0001	#lambda
	2	0.00000001	#epsilon
"""
	elif(squark == 2):
		text = """Block QCD
	1	%(alphaqs).12f	#alphaqs
Block REAL
	1	%(mu).12f	#mu
	2	%(Aqtree).12f	#Aq
	3	%(beta).12f	#beta
	4	0	#squark mixing
	5	%(sintheta).12f	#sintheta
	6	%(costheta).12f	#costheta
	7	%(totalcorr).12f	#correction_to_A_squark_coupling
	8	%(gq12).12f	#gq12_coupling
Block INTEGER
	1	2	#stop_1_sbot_2
Block MASS
	1	%(mA).12f	#mA
	2	%(mQ).12f	#mQ
	3	%(mQQQ).12f	#mQQQ
	4	%(mG).12f	#mG
	5	%(msq1).12f	#msq1
	6	%(msq2).12f	#msq2
Block VEGAS
	1	30	#iter
	2	10000	#stardiv
	3	3000000	#enddiv
	4	10000	#starfin
	5	3000000	#endfin
	6	0.4	#alpha
	7	1000	#increase
	8	Cuhre	#integrator
Block IBPPARA
	1	0.0001	#lambda
	2	0.00000001	#epsilon
"""
	dico = {'alphaqs':alphaqs, 'mu':mu, 'Aqtree':Aqtree, 'beta':np.arctan(tanb), 'sintheta':sintheta, 'costheta':costheta, 'totalcorr':totalcorrfactor, 'gq12':gq12, 'squark':squark, 'mA':mA, 'mQ':mQ, 'mQQQ':mQQQ, 'mG':mG, 'msq1':msq1, 'msq2':msq2}
	return open(output, 'w').write(text%dico)

read_to_slha(name, output)
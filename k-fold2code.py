#!/usr/bin/env python
#This tool allow users to plot SVM-prob ROC curve and get origin data
#This script should be run under the ../libsvm/python
from svmutil import *
from sys import argv, platform
from os import path, popen
from random import randrange , seed
from operator import itemgetter
from time import sleep

#transfer train_file,test_file into symbol file, those are 'trainfile' and 'testfile'
def score(train_file,test_file):
	n = 0
	m = 0
	p_symbol = {}
	n_symbol = {}
	list1 = ['A','Q','W','E','R','T','Y','U','I','O','P','S','D','F','G','H','J','K','L','Z','X','C','V','B','N','M']
	for i in list1:
		for j in range(31):
			p_symbol[i + str(j)] = 0
			n_symbol[i + str(j)] = 0
	train = open("trainscore","w")
	test1 = open('testscore',"w")
	with open(train_file) as t:
		for line in t.readlines():
			line =line.strip('\n')
			data = line.split('\t')
			peptide = data[0]
			mi = len(peptide)
			if data[1] == '1':
				n += 1
				for i in range(mi):
					p_symbol[peptide[i] + str(i)] += 1
			if data[1] == '-1':
				m += 1
				for i in range(mi):
					n_symbol[peptide[i] + str(i)] += 1
	with open(train_file) as t:
		for line in t.readlines():
			line =line.strip('\n')
			temp = line.split('\t')
			peptide = temp[0]
			data = []
			for i in range(63):
				data.insert(i,0)
			data[0] = temp[1]
			for i in range(mi):
				j = i + 1
				k = i + 32
				data[j] = p_symbol[peptide[i] + str(i)] / float(n)
				data[k] = n_symbol[peptide[i] + str(i)] / float(m)
				data[j] = str(j) + ":" + str(data[j])
				data[k] = str(k) + ":" + str(data[k])
			result = " ".join(data)
			train.write(result + "\n")
	train.close()
	with open(test_file) as t:
		for line in t.readlines():
			line =line.strip('\n')
			temp = line.split('\t')
			peptide = temp[0]
			data = []
			for i in range(63):
				data.insert(i,0)
			data[0] = temp[1]
			for i in range(mi):
				j = i + 1
				k = i + 32
				data[j] = p_symbol[peptide[i] + str(i)] / float(n)
				data[k] = n_symbol[peptide[i] + str(i)] / float(m)
				data[j] = str(j) + ":" + str(data[j])
				data[k] = str(k) + ":" + str(data[k])
			result = " ".join(data)
			test1.write(result + "\n")
	test1.close()


#svm_train svm_predict and it's output.deci is the score and model is traing by svm_train
def get_pos_deci(train_y, train_x, test_y, test_x, param):
	model = svm_train(train_y, train_x, param)
	#predict and grab decision value, assure deci>0 for label+,
	#the positive descision value = val[0]*labels[0]
	labels = model.get_labels()
	py, evals, deci = svm_predict(test_y, test_x, model)
	deci = [labels[0]*val[0] for val in deci]
	return deci,model

def get_cv_deci(train_file, param, nr_fold):
	if nr_fold == 1 or nr_fold == 0:
		score(train_file,train_file)
		trainfile = "trainscore"
		prob_y, prob_x = svm_read_problem(trainfile)
		if set(prob_y) != set([1,-1]):
			print("ROC is only applicable to binary classes with labels 1, -1")
			raise SystemExit
		deci,model = get_pos_deci(prob_y, prob_x, prob_y, prob_x, param)
		return deci
	deci, model, fileline = [], [], []
	with open(train_file) as f:
		for line in f.readlines():
			fileline.append(line)
	file_l = len(fileline)
	
	#random permutation by swapping i and j instance
	for i in range(file_l):
		j = randrange(i,file_l)
		fileline[i], fileline[j] = fileline[j], fileline[i]
	
	#cross training : folding
	for i in range(nr_fold):
		begin = i * file_l // nr_fold
		end = (i + 1) * file_l // nr_fold
		train_tem = fileline[:begin] + fileline[end:]
		test_tem = fileline[begin:end]
		train_te = open('trainfile_t.txt',"w")
		test_te = open('testfile_t.txt',"w")
		train_te.write("".join(train_tem))
		test_te.write("".join(test_tem))
		train_te.close()
		test_te.close()
		train_t = 'trainfile_t.txt'
		test_t = 'testfile_t.txt'
		score(train_t, test_t)
		trainfile = "trainscore"
		train_y, train_x = svm_read_problem(trainfile)
		testfile = "testscore"
		test_y, test_x = svm_read_problem(testfile)
		subdeci, submdel = get_pos_deci(train_y, train_x, test_y, test_x, param)
		deci += subdeci
	return deci

#deal with the command line
def proc_argv(argv = argv):
	#print("Usage: %s " % argv[0])
	#The command line : ./plotroc.py [-v cv_fold | -T testing_file] [libsvm-options] training_file
	train_file = argv[-1]
	test_file = None
	fold = 5
	options = []
	i = 1
	while i < len(argv)-1:
		if argv[i] == '-T': 
			test_file = argv[i+1]          #here my testfile is none
			i += 1
		elif argv[i] == '-v':
			fold = int(argv[i+1])
			i += 1
		else :
			options += [argv[i]]
		i += 1
	return ' '.join(options), fold, train_file, test_file    #other param are added to options


def plot_roc(deci, label, output, title):
	#count of postive and negative labels
	db = []
	pos, neg = 0, 0 		
	for i in range(len(label)):
		if label[i]>0:
			pos+=1
		else:	
			neg+=1
		db.append([deci[i], label[i]])

	#sorting by decision value
	db = sorted(db, key=itemgetter(0), reverse=True)

	#calculate ROC 
	xy_arr = []
	roc_data = []
	tp, fp = 0., 0.			#assure float division
	for i in range(len(db)):
		if db[i][1]>0:		#positive
			tp+=1
		else:
			fp+=1
		xy_arr.append([fp/neg,tp/pos])
		roc_data.append([db[i][0],(neg-fp+tp)/(neg+pos),tp/pos,1-fp/neg])
	#area under curve
	aoc = 0.			
	prev_x = 0
	for x,y in xy_arr:
		if x != prev_x:
			aoc += (x - prev_x) * y
			prev_x = x	
	return roc_data,aoc
	
#the main 
def main():
	if len(argv) <= 1:
		print("Usage: %s [-v cv_fold | -T testing_file] [libsvm-options] training_file" % argv[0])
		raise SystemExit
	param,fold,train_file,test_file = proc_argv()
	output_file = path.split(train_file)[1] + '-roc.png'
	score(train_file,train_file)
	trainfile = "trainscore"
	train_y, train_x = svm_read_problem(trainfile)
	if set(train_y) != set([1,-1]):
		print("ROC is only applicable to binary classes with labels 1, -1")
		raise SystemExit
	seed(0) #reset random seed
	output_title = path.split(train_file)[1]
	deci = get_cv_deci(train_file, param, fold)
	haha, auc = plot_roc(deci, train_y, output_file, output_title)
	lala=open('demo.txt','w')
	for i in haha:
		k='\t'.join([str(j) for j in i])
		lala.write(k+"\n")
	lala.close()
	print(auc)
	
if __name__ == '__main__':
	main()


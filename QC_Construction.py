##################################
#Author: Arjun Bhagoji
#Description: Program to find two girth 6 LDPC codes whose duals are contained in one another using the method given in Hagiwara and Imai (2007)
#Usage: The program takes in a prime number P and its perfect fulfillment given by the program Perfume.py and returns two compressed parity-check matrices using cPickle
##################################

import sys
from matplotlib.pyplot import *
import numpy as np
import random as random
import networkx as nx
import time as time
from scipy.sparse import *
import cPickle

#time1=time.time()
########################################
#Obtaining parameters P,sigma and tau for the code
print "Size of permutation matrix to be used:" 
P=int(sys.stdin.readline())
print "Fulfillment(sigma) to be used:"
sigma=int(sys.stdin.readline())
print "Tau co-prime to sigma:"
tau=int(sys.stdin.readline())
########################################

########################################
#Obtaining order of sigma mod P
count=1
while (pow(sigma,count)%P)!=1:
	count=count+1

order=count
L=2*order
########################################

########################################
#Creating model matrices
Ta=np.zeros((L/2,L/2))
Tb=np.zeros((L/2,L/2))

row_Ta=[]
row_Tb=[]
for i in range(L/2):
	elem_a=(pow(sigma,i))%P
	elem_b=(tau*elem_a)%P
	row_Ta.append(elem_a)
	row_Tb.append(elem_b)

row_Ta=np.array(row_Ta)
row_Tb=np.array(row_Tb)	
for i in range(L/2):
	Ta[i,:]=np.roll(row_Ta,i)
	Tb[i,:]=np.roll(row_Tb,i)
	
model_C=np.hstack([Ta,Tb])	
########################################

########################################
#Expanding Ta as a sparse matrix
rows=np.arange(P)
data=np.ones(P)
zero_zero_roll=Ta[0,0]
columns_zero=np.roll(rows,-int(zero_zero_roll))

zero_zero=coo_matrix((data,(rows,columns_zero)), shape=(P,P))

for k in range(1,L/2):
		current_roll=Ta[0,k]
		columns=np.roll(rows,-int(current_roll))
		curr_mat=coo_matrix((data,(rows,columns)), shape=(P,P))
		if k==1:
			row0=hstack([zero_zero,curr_mat])
		if k>1:
			row0=hstack([row0,curr_mat])

one_zero_roll=Ta[1,0]
columns_one=np.roll(rows,-int(one_zero_roll))

one_zero=coo_matrix((data,(rows,columns_one)), shape=(P,P))

for k in range(1,L/2):
		current_roll=Ta[1,k]
		columns=np.roll(rows,-int(current_roll))
		curr_mat=coo_matrix((data,(rows,columns)), shape=(P,P))
		if k==1:
			row1=hstack([one_zero,curr_mat])
		if k>1:
			row1=hstack([row1,curr_mat])
HC_part_A=vstack([row0,row1])

count=2

while count<L/2:
	rolli=Ta[count,0]
	columns_count=np.roll(rows,-int(rolli))
	count_zero=coo_matrix((data,(rows,columns_count)), shape=(P,P))
	for k in range(1,L/2):
		current_roll=Ta[count,k]
		columns=np.roll(rows,-int(current_roll))
		curr_mat=coo_matrix((data,(rows,columns)), shape=(P,P))
		if k==1:
			row_count=hstack([count_zero,curr_mat])
		if k>1:
			row_count=hstack([row_count,curr_mat])
	HC_part_A=vstack([HC_part_A,row_count])
	count=count+1
########################################

########################################
#Exapnding Tb as a sparse matrix
zero_zero_roll=Tb[0,0]
columns_zero=np.roll(rows,-int(zero_zero_roll))

zero_zero=coo_matrix((data,(rows,columns_zero)), shape=(P,P))

for k in range(1,L/2):
		current_roll=Tb[0,k]
		columns=np.roll(rows,-int(current_roll))
		curr_mat=coo_matrix((data,(rows,columns)), shape=(P,P))
		if k==1:
			row0=hstack([zero_zero,curr_mat])
		if k>1:
			row0=hstack([row0,curr_mat])

one_zero_roll=Tb[1,0]
columns_one=np.roll(rows,-int(one_zero_roll))

one_zero=coo_matrix((data,(rows,columns_one)), shape=(P,P))

for k in range(1,L/2):
		current_roll=Tb[1,k]
		columns=np.roll(rows,-int(current_roll))
		curr_mat=coo_matrix((data,(rows,columns)), shape=(P,P))
		if k==1:
			row1=hstack([one_zero,curr_mat])
		if k>1:
			row1=hstack([row1,curr_mat])
HC_part_B=vstack([row0,row1])

count=2

while count<L/2:
	rolli=Tb[count,0]
	columns_count=np.roll(rows,-int(rolli))
	count_zero=coo_matrix((data,(rows,columns_count)), shape=(P,P))
	for k in range(1,L/2):
		current_roll=Tb[count,k]
		columns=np.roll(rows,-int(current_roll))
		curr_mat=coo_matrix((data,(rows,columns)), shape=(P,P))
		if k==1:
			row_count=hstack([count_zero,curr_mat])
		if k>1:
			row_count=hstack([row_count,curr_mat])
	HC_part_B=vstack([HC_part_B,row_count])
	count=count+1
########################################

########################################
#Creating the parity check matrix for C
HC=hstack([HC_part_A,HC_part_B])
#time2=time.time()
row_indices_C=HC.row
column_indices_C=HC.col

########################################
#Creating the parity check matrix for D
HD_part_A=coo_matrix.transpose(HC_part_B)
HD_part_B=coo_matrix.transpose(HC_part_A)
HD=hstack([HD_part_A,HD_part_B])
row_indices_D=HD.row
column_indices_D=HD.col
########################################

########################################
#Storing the parity check matrices 
f = open('qcC.dat','wb')
cPickle.dump(HC,f,0)
f.close()

f = open('qcD.dat','wb')
cPickle.dump(HD,f,0)
f.close()
########################################

########################################
#Creating Tanner graphs from both C and D to find out girth
n=L*P

tannerC=nx.Graph()
tannerD=nx.Graph()


for k in range(HC.nnz):
	tannerC.add_edge(column_indices_C[k],n+row_indices_C[k],weight=0)
	tannerD.add_edge(column_indices_D[k],n+row_indices_D[k],weight=0)
#time3=time.time()	
cycles_C=nx.cycle_basis(tannerC)
girth_C=min(map(len,cycles_C))
cycles_D=nx.cycle_basis(tannerD)
girth_D=min(map(len,cycles_D))
########################################
print "C is a (" + str(order) + "," + str(L) + ") LDPC code with block length " + str(n)
print "The girth of C is %d" %girth_C
print "D is a (" + str(order) + "," + str(L) + ") LDPC code with block length " + str(n)
print "The girth of D is %d" %girth_D

#print time2-time1
#print time3-time1

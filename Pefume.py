##################################
#Author: Arjun Bhagoji
#Description: Program to find a Perfect Fulfillment as defined in Hagiwara and Imai (2007)
#Limited: We only consider order 3
#Usage: The program takes in a prime number P and outputs a perfect fulfillment if it exists
##################################

import numpy as np
from fractions import gcd
import sys

##################################
def numpy_gcd(a, b):
    a, b = np.broadcast_arrays(a, b)
    a = a.copy()
    b = b.copy()
    pos = np.nonzero(b)[0]
    while len(pos) > 0:
        b2 = b[pos]
        a[pos], b[pos] = b2, a[pos] % b2
        pos = pos[b[pos]!=0]
    return a
##################################

##################################
def order_find(x):
	count=1
	while (pow(x,count)%P)!=1:
		count=count+1
	orde=count
	return orde
##################################


print "Please input the prime number for which the perfect fulfillment is to be found \n"

P=int(sys.stdin.readline())
order=3
sigma=[]

for i in range(P-2):
	curr=i+2
	#print curr
	co_prime=False
	correct_order=False
	flags=np.zeros(order-1)
	if order_find(curr)==order:
		correct_order=True
	if gcd(curr,P)==1:
		co_prime=True
	for j in range(order-1):
		if gcd((curr**(j+1))-1,P)==1:
			flags[j]=1
	if correct_order and co_prime and np.array_equal(flags,np.ones(order-1))==True:
		sigma.append(curr)
	#if curr==92:
	#	sys.exit()

print "The perfect fulfillment is (sigma,tau): \n" 
print sigma

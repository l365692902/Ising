#fin_T1 = open("Ising_N10_T1.txt", 'r')
#fin_T4 = open("Ising_N10_T4.txt", 'r')
#fout=open("Ising_N10.txt", 'w')

#for x in range(10000):
#    fout.write("L: "+fin_T1.readline())
#    fout.write("H: "+fin_T4.readline())

#fin_T1.close()
#fin_T4.close()
#fout.close()

import numpy as np
import os
print(os.getcwd())

l_images=np.genfromtxt('Ising_N10_T1.txt',dtype=np.float32)
print(l_images.shape)
h_images=np.genfromtxt('Ising_N10_T4.txt',dtype=np.float32)
print(h_images.shape)
print(h_images.dtype)
images=np.concatenate((l_images,h_images),axis=0)
print(images.shape)

label_one=np.ones((10000,1),dtype=np.float32)
label_zero=np.zeros((10000,1),dtype=np.float32)

l_label=np.concatenate((label_zero,label_one),axis=1)
h_label=np.concatenate((label_one,label_zero),axis=1)
print(l_label)
print(h_label)
labels=np.concatenate((l_label,h_label),axis=0)
print(labels)

perm=np.arange(20000)
np.random.shuffle(perm)
images=images[perm]
labels=labels[perm]

print(perm)
print(labels)

np.savez('./Ising_N10_20000',images=images,labels=labels)

#data=np.load('./savetest.npz')
#y=np.random.random_sample((3,4))
#print(data['y'])
#x=np.random.random_sample((3,2))
#print(data['x'])
#np.savez('./savetest',y=y,x=x)
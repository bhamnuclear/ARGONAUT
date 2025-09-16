import numpy as np
import matplotlib.pyplot as plt

from RmTools import Argonaut

ex_list = [0.0,2.345,2.75]
scales = [5e-2,10,5]

xlist = np.append(np.arange(-0.5,-0.01,0.001),np.arange(-0.01,0.01,0.00005))
xlist = np.append(xlist,np.arange(0.01,3.5,0.001))


stdev = np.zeros_like(xlist) + 0.01


m1_dch1 = 8 
m2_dch1 = 1
z1_dch1 = 4
z2_dch1 = 1
qval_8bep = 0.1859
l_dch1 = [1,3,2]

m1_dch2 = 5 
m2_dch2 = 4
z1_dch2 = 3
z2_dch2 = 2
qval_5lia = -1.6873
l_dch2 = [0,2,1]


gamma_uni = [0.54e-3,0.081,0.614]
branches = [[1,0.01,1.],[0.,0.99,0.]] 



ls = [l_dch1,l_dch2]
qvals = [qval_8bep,qval_5lia]
masses = [[m1_dch1,m2_dch1],[m1_dch2,m2_dch2]]
zs = [[z1_dch1,z2_dch1],[z1_dch2,z2_dch2]]



lnshp = Argonaut(xlist,ex_list,qvals,gamma_uni,branches,ls,masses,zs,stdev,scales)

state1_contrb = Argonaut(xlist,[ex_list[0]],qvals,[gamma_uni[0]],
                         [[branches[0][0]],[branches[1][0]]],[[ls[0][0]],[ls[1][0]]],
                         masses,zs,stdev,[scales[0]])

state2_contrb = Argonaut(xlist,[ex_list[1]],qvals,[gamma_uni[1]],
                         [[branches[0][1]],[branches[1][1]]],[[ls[0][1]],[ls[1][1]]],
                         masses,zs,stdev,[scales[1]])

state3_contrb = Argonaut(xlist,[ex_list[2]],qvals,[gamma_uni[2]],
                         [[branches[0][2]],[branches[1][2]]],[[ls[0][2]],[ls[1][2]]],
                         masses,zs,stdev,[scales[2]])



fig = plt.figure(figsize=(8,6),layout="tight")
ax = fig.add_subplot(1,1,1)

ax.plot(xlist,lnshp,c='r',label='Combined')
ax.plot(xlist,state1_contrb,ls='-.',label='State 1')
ax.plot(xlist,state2_contrb,ls='-.',label='State 2')
ax.plot(xlist,state3_contrb,ls='-.',label='State 3')

ax.legend(frameon=False,draggable=True)
plt.show()
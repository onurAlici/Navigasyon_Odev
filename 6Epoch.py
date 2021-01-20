import math
import numpy as np

set = ["G02","G06","G12","G13","G15","G17","G19","G24","G25","G29","G32"]
y_alici = [4208830.7150, 2334850.1024, 4171267.0925]
c1 = [23123177.528, 22460793.836, 20757135.281, 25563871.306, 24260511.875, 23606970.487, 21664607.314, 20261240.595,
      23771970.989, 25536690.367,24591756.095]
p2 = [23123172.790, 22460793.738, 20757139.357, 25563871.021, 24260509.776, 23606966.978, 21664601.548, 20261239.568,
      23771969.962, 25536688.228,24591755.414]
koord = np.array([
    [14667.409870, 22167.121056, -1169.070842, -249.908174],
    [4506.444368,  23846.938435,  10790.948740,    162.582573],
    [17532.597734,  -1991.791699,  19580.957159,    228.121742],
    [18537.275434 , 12832.891056 ,-14189.872859 ,   -56.633532],
    [25723.866065  , 3499.082198,  -6779.156975,   -303.707829],
    [-5477.497451,  14281.723371,  22112.429239 ,    80.947866],
    [2058.622839  ,14847.957342,  21614.589026,   -299.803023],
    [19111.230112,   8293.614740,  16504.184606,    -68.193145],
    [18520.851896 ,-14926.535021 , 11340.009006 ,  -770.048186],
    [23974.391372,  -4625.447593, -10508.249800  ,  162.834479],
    [4117.007230, -15085.485562,  21519.570219  ,  -21.087791]])


ısık = 299792458
cd = np.add(c1,ısık*(koord[:,3]/10**6)   )
pd = np.add(p2,ısık*(koord[:,3]/10**6)   )
koord_m = koord[:,:3] * 1000





def q():
    x = y_alici[0]
    y = y_alici[1]
    z = y_alici[2]
    xyz = np.array([[math.sqrt((k[0]-x)**2+(k[1]-y)**2+(k[2]-z)**2) for k in koord_m]])
    return xyz

q = q()

l = np.subtract(q,pd).T



def A():
    coord = np.array([[y_alici[0] for i in range(11)],[y_alici[1] for i in range(11)],[y_alici[2] for i in range(11)]]).T
    ara = np.subtract(coord,koord_m)
    sonuc = np.divide(ara,np.tile(q,(3,1)).T)
    c_ekle = np.append(sonuc, np.array([[ısık for i in range(11)]]).T,axis=1)
    return c_ekle

A = A()

N = np.matmul(A.T,A)

Q = np.linalg.inv(N)

n = np.matmul(A.T,l)

x = np.matmul(Q,n)

v = np.subtract(np.matmul(A,x),l)

vtv = np.matmul(v.T,v)[0][0]

s0 = math.sqrt(vtv/len(set))

qxx = Q[0][0]
qyy = Q[1][1]
qzz = Q[2][2]
sx = s0 + math.sqrt(qxx)
sy = s0 + math.sqrt(qyy)
sz = s0 + math.sqrt(qzz)

konumsapma = math.sqrt(sx**2 + sy**2 + sz**2 )

PDOP = math.sqrt(qxx + qyy + qzz)

def q2(xx):
    x = y_alici[0] + xx[0][0]
    y = y_alici[1] + xx[1][0]
    z = y_alici[2] + xx[2][0]
    xyz2 = np.array([[math.sqrt((k[0]-x)**2+(k[1]-y)**2+(k[2]-z)**2)+ ısık * xx[3][0] for k in koord_m]])

    return xyz2

dq = q2(x)

cson = np.add(cd,v.T)

pson = np.add(pd,v.T)

fark = np.subtract(dq,pson)
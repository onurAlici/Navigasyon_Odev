import math
import numpy as np

set = ["G01", "G03", "G08", "G11", "G14", "G17", "G18", "G22", "G23", "G31", "G32"]
y_alici = [4208830.7150, 2334850.1024, 4171267.0925]
c1 = [20242409.795, 21337105.198, 24110870.179, 20522456.564, 21879220.337, 24413507.774, 20309310.027, 20940223.040,
      23236714.929, 22401231.486, 23162819.400]
p2 = [20242409.227, 21337104.289, 24110871.195, 20522451.404, 21879216.863, 24413504.368, 20309307.249, 20940216.621,
      23236710.040, 22401227.562, 23162818.269]
koord = np.array([
    [17734.832987, 6362.324408, 18659.390158, -55.647160],
    [18198.966069, -4947.272407, 18623.424103, 183.306038],
    [24955.644260, 5479.242121, -7612.538702, -149.536099],
    [23028.750566, 6143.837990, 10914.388347, -571.550886],
    [3656.291316, 15558.887256, 21565.447508, -73.517531],
    [5798.905820, -14251.888542, 22047.445371, 81.233159],
    [18665.325180, 14336.635141, 11933.032450, 66.460078],
    [16218.788887, 2931.001319, 21048.493423, -725.100787],
    [26763.957993, -252.708243, -521.068660, -179.227400],
    [4834.835267, 24103.964152, 9452.833703, 27.275273],
    [-4458.733232, 15071.523501, 21459.576879, -20.411587]])
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

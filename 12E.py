import math
import numpy as np

set = ["G10", "G13", "G15", "G16", "G20", "G21", "G25", "G26", "G27", "G29", "G31"]
y_alici = [4208830.7150, 2334850.1024, 4171267.0925]

c1 = [23260507.595, 25404373.085, 24309138.782, 21574916.122, 21234197.986, 20946075.730, 25598803.247, 20806390.236,
      23862798.574, 22508656.198, 25510184.570]
p2 = [23260508.389, 25404371.141, 24309137.425, 21574912.310, 21234194.824, 20946071.327, 25598805.199, 20806390.550,
      23862797.833, 22508653.322, 25510184.293]
koord = np.array([
    [23933.005055, 10653.682785, -4927.134635, -1.489137],
    [-12687.219408, 18454.832890, 14130.587360, -56.596360],
    [-3135.341771, 25303.669533, 6410.045923, -303.656028],
    [15582.628127, -5595.339298, 20680.172478, -52.546332],
    [22069.999028, 14088.276871, 4837.048737, 525.727336],
    [14949.402182, 8458.238108, 21018.792942, -162.851994],
    [14865.569199, 18816.779447, -11900.863906, -770.111247],
    [22303.873428, 992.537587, 14499.786512, 240.528010],
    [12005.006209, -16010.943105, 17194.774063, -172.722768],
    [4551.061112, 24055.224504, 10249.666956, 162.627349],
    [24629.461063, -5234.977000, -9128.458256, 27.332636]])
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

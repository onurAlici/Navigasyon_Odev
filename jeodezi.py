import math
import numpy as np

y_alici = [4208830.7150,  2334850.1024,  4171267.0925]
set = [
"G02",
"G05",
"G06",
"G07",
"G08",
"G09",
"G13",
"G23",
"G27",
"G28",
"G30" ]


c1 = np.array([[24833572.746,22006938.883,25521548.613,20875657.148,24767366.920,21329534.125,24650865.997,24293787.311,24955140.621,22719337.867,20503335.870]])
p2 = np.array([[24833567.012,22006935.228,25521548.179,20875652.645,24767366.519,21329532.679,24650861.741,24293781.788,24955138.766,22719332.701,20503333.926]])

koord =np.array([
[21324.613932,	-15515.371491,	2982.570564,	-249.718929],
[15019.722616,	-7988.797021,	20322.237779,	0.302647],
[23945.035781,	-4660.331125,	-10504.958406,	162.819348],
[10554.770834,	11793.945778,	21645.770940,	-55.156770],
[-5533.370036,	24853.668757,	7203.784396,	-149.438247],
[8954.453201,	21162.659242,	13272.023023,	380.623615],
[12695.543065,	-18677.922092,	13828.514853,	-56.670101],
[996.973499,	26265.528951,	2029.230276,	-179.342026],
[-11992.598049,	16282.758274,	16944.395465,	-172.270788],
[22405.558236,	13120.115741,	-4729.476710,	767.000723],
[19786.452815,	4524.922698,	17258.625081,	-134.189812]])
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

###
# Analysis computes an overview of the evolution of the massive bodies
###

from numpy import pi, sqrt
import numpy as np
import matplotlib.pyplot as plt
from SketchOloc import SketchOloc
import os

G = 4*(pi**2)

#########

def orbel_xv2el(xo,yo,zo,vxo,vyo,vzo,gm):

    Cx = yo*vzo - zo*vyo
    Cy = zo*vxo - xo*vzo
    Cz = xo*vyo - yo*vxo

    C = sqrt(Cx*Cx + Cy*Cy + Cz*Cz)
    v2 = vxo**2 + vyo**2 + vzo**2
    vr = xo*vxo + yo*vyo + zo*vzo
    r = sqrt(xo**2 + yo**2 + zo**2)
    E = 0.5*v2 - gm/r

    i = np.arccos(Cz/C)

    if sqrt(Cx**2+Cy**2)/C < 0.001:
        O = 0.
        aux = np.arctan2(yo,xo)
        if (abs(i-pi)<0.01):
            aux = -aux
    else:
        O = np.arctan2(Cx,-Cy)
        aux = np.arctan2(zo/np.sin(i), xo*np.cos(O)+yo*np.sin(O))

    if E < 0:
        a = - 0.5*gm/E
        e = sqrt(1.-C**2/(gm*a))
        if (e != 0.):
            if ((a-r)/(a*e)>1):
                u = 0.
            elif ((a-r)/(a*e)<-1):
                u = pi
            else:
                u = np.arccos((a-r)/(a*e))
            if (vr<0):
                u = 2*pi-u
            w = np.arctan2(sqrt(1-e**2)*np.sin(u)/(1-e*np.cos(u)), (np.cos(u)-e)/(1-e*np.cos(u)))
            w = aux-w
            M = u-e*np.sin(u)
        else:
            w = 0.
            M = aux

    elif E > 0:
        a = 0.5*gm/E
        e = sqrt(1.+C**2/(gm*a))
        u = np.log((a+r)/(a*e)+sqrt(((a+r)/(a*e))**2-1))
        if (vr<0):
            u = -u
        w = np.arctan2(sqrt(e**2-1)*np.sinh(u)/(e*np.cosh(u)-1), (e-np.cosh(u))/(e*np.cosh(u)-1))
        w = aux-w
        M = e*np.sinh(u)-u

    return(a,e,i,O,w,M)

def oloc2mat(N,oloc,mass):

   eta = np.zeros(N)
   mu = np.zeros(N)

   for k in range(1, N):
      for j in range(N):
         if oloc[k,j] == 1:
            mu[k] += mass[j]
         elif oloc[k,j] == -1:
            eta[k] += mass[j]

   mat = np.zeros((N, N))
   umat = np.zeros((N, N))

   mtot = 0
   for j in range(N):
      mtot += mass[j]

   mat[0,:] = mass/mtot

   for k in range(1, N):
      for j in range(N):
         if oloc[k,j] == 1:
            mat[k, j] = mass[j]/mu[k]
         elif oloc[k,j] == -1:
            mat[k, j] = -mass[j]/eta[k]


   umat[:,0] = np.full(N, 1.)

   for k in range(1, N):
      for j in range(N):
         if oloc[k, j] == 1:
            umat[j, k] = eta[k]/(eta[k]+mu[k])
         elif oloc[k, j] == -1:
            umat[j, k] = -mu[k]/(eta[k]+mu[k])

   return(mu, eta, mat, umat)

#########

dr = pi/180.

inparfile = open("paramhjs.in",'r')

for j in range(4):
    inparfile.readline()
dirs = inparfile.readline().strip()
gname = inparfile.readline().strip()
diro = dirs+'/'+gname
print diro

inparfile.close()

os.chdir(diro)

elbodies = open("elbodies.dat",'r')

nbod = int(elbodies.readline().strip().split(" ")[0])
Norbits = nbod-1
print "Number of bodies:",nbod
data = np.loadtxt(elbodies,dtype='float')

elbodies.close()

Nframe = data.shape[0]/nbod

t = np.zeros(Nframe)
a = np.zeros((Norbits,Nframe))
e = np.zeros((Norbits,Nframe))
i = np.zeros((Norbits,Nframe))
O = np.zeros((Norbits,Nframe))
w = np.zeros((Norbits,Nframe))
M = np.zeros((Norbits,Nframe))
for j in range(Nframe):
    t[j] = data[nbod*j][0]
    for k in range(Norbits):
        num = int(-data[nbod*j+1+k][0]-2)
        a[num][j] = data[nbod*j+1+k][1]
        e[num][j] = data[nbod*j+1+k][2]
        i[num][j] = data[nbod*j+1+k][3]/dr
        O[num][j] = data[nbod*j+1+k][4]/dr
        w[num][j] = data[nbod*j+1+k][5]/dr
        M[num][j] = data[nbod*j+1+k][6]/dr

fig = plt.figure(figsize=(25, 15))

for j in range(Norbits):
    for k in range(6):
        plt.subplot(Norbits, 6, j*6+k+1)
        plt.xlabel("Time (yr)")
        if k == 0:
            plt.plot(t, a[j], '.', color='C%s'%j)
            if j == 0:
                plt.title(r"$a$ (au)")
        if k == 1:
            plt.plot(t, e[j], '.', color='C%s'%j)
            if j == 0:
                plt.title(r"$e$")
        if k == 2:
            plt.plot(t, i[j], '.', color='C%s'%j)
            if j == 0:
                plt.title(r"$i$ (deg)")
        if k == 3:
            plt.plot(t, O[j], '.', color='C%s'%j)
            if j == 0:
                plt.title(r"$\Omega$ (deg)")
        if k == 4:
            plt.plot(t, w[j], '.', color='C%s'%j)
            if j == 0:
                plt.title(r"$\omega$ (deg)")
        if k == 5:
            plt.plot(t, M[j], '.', color='C%s'%j)
            if j == 0:
                plt.title(r"$M$ (deg)")

fig.subplots_adjust(left=0.05, right=0.98)
plt.savefig("Evolution.png")
plt.close(fig)

try:
    olocfile = open("oloc.out","r")

    data = np.loadtxt(olocfile)
    tchange = SketchOloc(data, Nframe, nbod)

    olocfile.close()
except:
    tchange = []

ce = []
try:
    cefile = open("ce.out", "r")
except:
    print "No close encounter"
    noce = True
else:
    noce = False
    for line in cefile.readlines():
        if line != '':
            line = line.split()
            k = int(line[0])
            tcebeg = float(line[1])
            if (len(line)>2):
                tceend = float(line[2])
            else:
                tceend = t[-1]
            ce += [[k, tcebeg, tceend]]
    cefile.close()

    ce = np.array(ce)

    tcebeg = [ ce[ce[:,0]==k+2][:,1] for k in range(Norbits)]
    tceend = [ ce[ce[:,0]==k+2][:,2] for k in range(Norbits)]

    for k in range(Norbits):
        print len(tcebeg[k]),"close encounters for orbit",k

print len(tchange), "changes of hierarchy"

if len(tchange) > 0:

    plhjs = open("plhjs.in","r")

    if int(plhjs.readline().strip()) != nbod:
        print "Problem with nbod"

    mass = np.zeros(nbod)
    for j in range(nbod):
        mass[j] = float(plhjs.readline().strip())/G
        plhjs.readline(); plhjs.readline()

    plhjs.close()

    oloc = np.zeros((nbod, nbod, Nframe))
    for n in range(Nframe):
        for j in range(nbod):
            oloc[1:nbod, j, n] = data[n][1+(nbod-1)*j:nbod+(nbod-1)*j]

    mu = np.zeros((nbod, Nframe)); eta = np.zeros((nbod, Nframe))
    mat = np.zeros((nbod, nbod, Nframe)); umat = np.zeros((nbod, nbod, Nframe))
    for n in range(Nframe):
        mu[:, n], eta[:, n], mat[:, :, n], umat[:, :, n] = oloc2mat(nbod, oloc[:, :, n], mass)

    xvbodies = open("xvbodies.dat",'r')

    nbod = int(xvbodies.readline().strip().split(" ")[0])
    data = np.loadtxt(xvbodies)

    xvbodies.close()

    x = np.zeros((nbod, Nframe))
    y = np.zeros((nbod, Nframe))
    z = np.zeros((nbod, Nframe))
    vx = np.zeros((nbod, Nframe))
    vy = np.zeros((nbod, Nframe))
    vz = np.zeros((nbod, Nframe))
    for j in range(Nframe):
        for k in range(nbod):
            num = int(data[(nbod+1)*j+1+k][0]-1)
            x[num][j] =  data[(nbod+1)*j+1+k][1]
            y[num][j] =  data[(nbod+1)*j+1+k][2]
            z[num][j] =  data[(nbod+1)*j+1+k][3]
            vx[num][j] =  data[(nbod+1)*j+1+k][4]
            vy[num][j] =  data[(nbod+1)*j+1+k][5]
            vz[num][j] =  data[(nbod+1)*j+1+k][6]

    xj = [sum(mat[k]*x) for k in range(nbod)]
    yj = [sum(mat[k]*y) for k in range(nbod)]
    zj = [sum(mat[k]*z) for k in range(nbod)]
    vxj = [sum(mat[k]*vx) for k in range(nbod)]
    vyj = [sum(mat[k]*vy) for k in range(nbod)]
    vzj = [sum(mat[k]*vz) for k in range(nbod)]

    for j in range(Nframe):
        for k in range(Norbits):
            a[k][j], e[k][j], i[k][j], O[k][j], w[k][j], M[k][j] = orbel_xv2el(xj[k+1][j], yj[k+1][j], zj[k+1][j], vxj[k+1][j], vyj[k+1][j], vzj[k+1][j], G*(mu[k+1][j]+eta[k+1][j]))

            i[k][j] = (i[k][j]/dr) % 360
            w[k][j] = (w[k][j]/dr) % 360
            O[k][j] = (O[k][j]/dr) % 360
            if e[k][j] < 1:
                M[k][j] = (M[k][j]/dr) % 360
            else:
                M[k][j] /= dr

    fig = plt.figure(figsize=(25, 15))

    for j in range(Norbits):
        for k in range(6):
            plt.subplot(Norbits, 6, j*6+k+1)
            plt.xlabel(r"Time (yr)")
            if k == 0:
                plt.plot(t, a[j], '.', color='C%s'%j)
                if j == 0:
                    plt.title(r"$a$ (au)")
            if k == 1:
                plt.plot(t, e[j], '.', color='C%s'%j)
                if j == 0:
                    plt.title(r"$e$")
            if k == 2:
                plt.plot(t, i[j], '.', color='C%s'%j)
                if j == 0:
                    plt.title(r"$i$ (deg)")
            if k == 3:
                plt.plot(t, O[j], '.', color='C%s'%j)
                if j == 0:
                    plt.title(r"$\Omega$ (deg)")
            if k == 4:
                plt.plot(t, w[j], '.', color='C%s'%j)
                if j == 0:
                    plt.title(r"$\omega$ (deg)")
            if k == 5:
                plt.plot(t, M[j], '.', color='C%s'%j)
                if j == 0:
                    plt.title(r"$M$ (deg)")
            for n in range(len(tchange)):
                plt.axvline(x=tchange[n], color='r')
            if not noce:
                for tbeg, tend in zip(tcebeg[j], tceend[j]):
                    plt.axvspan(tbeg, tend, color=(1,0.8,0.8), alpha=0.5)

    fig.subplots_adjust(left = 0.05, right = 0.98)
    plt.savefig("Evolution_withchange.png")
    plt.close(fig)

    del x, y, z, vx, vy, vz
    del xj, yj, zj, vxj, vyj, vzj

del t, a, e, i, O, w, M

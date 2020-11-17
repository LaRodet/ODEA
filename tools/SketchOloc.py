###
# SketchOloc represents the evolution of the hierarchy in a txt file
###

import os
import numpy as np

def SketchOloc(data, Nframe, N):

    t = data[:,0]
    oloc = data[:,1:]

    Nframe = len(t)
    l = len(oloc[0])
    N = int((1+int(np.sqrt(1+4*l)))/2)

    tchange = []; nchange = [0]
    for n in range(1,Nframe):
        if ((oloc[n]!=oloc[n-1]).any()):
            tchange += [0.5*(t[n]+t[n-1])]
            nchange += [n]

    Nchange = len(tchange)

    hierarchy = open("hierarchy.dat","w")

    for nch in range(Nchange+1):

        if (nch>0):
            hierarchy.write("%s yr\n" % tchange[nch-1])

        ### Reorganize
        olocmat = (np.reshape(oloc[nchange[nch]],(N,N-1))).astype('int')

        r = range(1, N+1)
        order = [*r]
        for n in range(N-1):
            okaux = np.ones(N,dtype=bool)
            for j in range(N-1):
                if (olocmat[n,j]!=0):
                    aux = (olocmat[n+1:,j] == olocmat[n,j])
                    if (sum(aux)>0):
                        okaux[n+1:] = okaux[n+1:] & aux
                    else:
                        aux = (olocmat[n+1:,j] == -olocmat[n,j])
                        if (sum(aux)>0):
                            okaux[n+1:] = okaux[n+1:] & aux

            ok = False; m = n
            while (not ok and m<N-1):
                m += 1
                if okaux[m]:
                    ok = True
            if (m>n+1):
                #        print "switch",n+1,m
                aux = np.copy(olocmat[n+1,:])
                olocmat[n+1,:] = olocmat[m,:]
                olocmat[m,:] = aux
                aux = order[n+1]
                order[n+1] = order[m]
                order[m] = aux

        text = [ "" for n in range(2*N)]

        for n in order:
            text[0] += "%s " % n
        text[0] += "\n"

        orbits = [[] for k in range(N)] # Listes des orbites contenant k elements
        for o in range(N-1):
            orbits[len(olocmat[abs(olocmat[:,o])!=0,o])-1] += [o]

        for k in range(1,N):

            for n in range(N-1):
                if any([olocmat[n,o]*olocmat[n+1,o]==-1 for o in orbits[k]]):
                    text[k] += " _"
                else:
                    text[k] += "  "
            text[k] += "\n"

        for k in range(2*N):
            hierarchy.write(text[k])

        hierarchy.write("\n")
        hierarchy.write("********************")
        hierarchy.write("\n")

    hierarchy.close()

    return(tchange)

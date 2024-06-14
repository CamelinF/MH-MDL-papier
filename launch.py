
import os
import sys
import time
import math
import random
import skmine
import time
import bisect
import numpy as np



from skmine.datasets.utils import describe
from skmine.datasets.fimi import fetch_file
from skmine.itemsets import LCM
from random import sample,random,choices,randint

dataTypeKrimp = {
    "connect" : "bai32",
    "hepatitis" :"bm128",
    "mushroom" : "bm128",
    "splice1" : "bai32",
    "pumsb" : "uint16",
    "pumsb_star" : "uint16",
    "eisen" : "uint16",
    "retail" :"uint16",
    "chess" :"bm128",
}
def get_Transaction_Weight(transaction):
    x=len(transaction)
    return 2**x-x-1

def get_Sample(data,info,n,f,choixP):
    indexT=[]
    sampleT=[]
    ind=sample(choixP,n)

    for i in ind:
        indexT.append(i)
        sampleT.append(data[i])
    #freq={}
    # x=set()
    # xP=set()
    # for i in sampleT:
    #     for j in i:
    #         if j in freq:
    #             freq[j]=freq[j]+1
    #             if j in x:
    #                 if freq[j]>=n*f:
    #                     x.remove(j)
    #                     xP.add(j)
    #         else:
    #             freq[j]=1
    #             x.add(j)
    xP=[]
    test=sampleT[0]

    for j in test:
        acc=True
        for i in range(1,n):
            if acc and not(j in sampleT[i]):
                acc=False
        if acc:
            xP.append(j)
    # print(f"             taille de l'intersection {len(xP)}")
    poids= get_Transaction_Weight(xP)
    potentiel=1
    ensPot=set()

    if info["isRect"]:
        potentiel=get_Transaction_Weight(sampleT[0])
    else:
        for i in sampleT:
        # for j in i:
            # ensPot=ensPot.union([j])
            potentiel=potentiel*(get_Transaction_Weight(i)**(1/n))
            #print(potentiel)
    #potentiel=get_Transaction_Weight(ensPot)
    #print(f"poids={poids},potentiel={potentiel}, valeur=, log/log ="+str(math.log2(poids)/math.log2(potentiel))+" log(valeur) ="+str(math.log2(poids/potentiel)))
    #print(f"potentiel={potentiel } , poids={poids}")
    if poids==0:
        res=0
        poids=1
    else:
        res=poids/potentiel
    #res=poids/potentiel
    #print(f"res={res}")
    return sampleT,res,indexT#[poids,potentiel]
    #return sampleT,len(xP)/len(ensPot)

def cftp_Sampler(data,info,n,f,choixP):
    currValue=1
    indexT=[]
    currState=[]
    seenSample=[]
    seenValue=[]
    seenIndex=[]
    seenRand=[]
    i=0
    while currState==[]:
        i=i+1
        k=2**i
        #print(f"    {k}")
        maxT=0
        if i > 10:
            newState,newValue,newIndex=get_Sample(data,info,n,f,choixP)
            currState=newState
            if newValue==0:
                currValue=0.1
            else:
                currValue=newValue
            indexT=newIndex
        while k>0:

            newState,newValue,newIndex=get_Sample(data,info,n,f,choixP)
            test=random()
            #print(math.log(test))
            x=newValue/currValue
            #x=math.log(currValue[1])+math.log(newValue[0])-math.log(currValue[0])-math.log(newValue[1])
            #if math.log(test)> x:
            #print(f"valeur de k ={newValue}, valeur curr={currValue}, valeur du test ={x} et valeur à battre={test}")
            if test<x:
                currState=newState
                currValue=newValue
                indexT=newIndex
            seenSample.insert(0,newState)
            seenValue.insert(0,newValue)
            seenRand.insert(0,test)
            seenIndex.insert(0,newIndex)
            k-=1
        if currState!=[]:
            #print("ok")
            for k in range(2*i,len(seenSample)):
                #x=math.log(currValue[1])+math.log(seenValue[0])-math.log(currValue[0])-math.log(seenValue[1])
                x=seenValue[k]/currValue

                if seenRand[k]< x:
                    #print(f"valeur de k ={seenValue[k]}, valeur curr={currValue}, valeur du test ={x} et valeur à battre={seenRand[k]}")
                    currState=seenSample[k]
                    currValue=seenValue[k]
                    indexT=seenIndex[k]
    # print(f" valeur choisie={currValue}, valeur max {max(seenValue)},valeur moy={sum(seenValue)/len(seenValue)}")
    # print(f" on s'est arreter au bout de {i}epoch")
    return currState,indexT

def run_cftp(dataFile,tailleSample):
    os.system(f"java -jar cftp-sampling/target/CouplingFromThePastRunner-jar-with-dependencies.jar data/{dataFile}.dat result.txt 1 -m1{tailleSample} > CFTP_results.txt")
    f=open("CFTP_results.txt","r")
    s=f.readline()
    temp=""
    ens=s
    while temp!="Ensemble" and s!="":
        temp=s[0]
        if len(s)>7:
            for x in range(1,8):
                temp=temp+s[x]
            ens=s
        s=f.readline()
    f.close()
    res=[]
    actu=""
    for i in range(10,len(ens)):
        if ens[i]!=' ':
            actu=actu+ens[i]
        else:
            res.append(actu)
            actu=""
    ensTrans=[]
    for i in range(tailleSample):
        ensTrans.append(int(res[i]))
    return ensTrans

def getFreq(d,m):
    couv=d[m[0]]
    for i in m:
        couv=couv.intersection(d[i])
        #print(f"{m} à uen couv de {len(couv)}")
    return couv,len(couv)

# OBTENTION DU DECODAGE
# sens ==FALSE => original vers Krimp, sens True => Krimp vers original
def getDecodage(dataFile,sens):
    f=open(f"data/{dataFile}.db")
    s=f.readline()
    while s[0]!="a" and s[1] !="b":
        s=f.readline()
    alphabet=[]
    splitS=s.split()
    for i in range(1,len(splitS)):
        alphabet.append(int(splitS[i]))

    while s[0]!="i" and s[1] !="t":
        s=f.readline()
    items=[]
    splitS=s.split()
    for i in range(1,len(splitS)):
        items.append(int(splitS[i]))
    f.close()
    dico={}
    for i in range(len(items)):
        if sens :
            dico[alphabet[i]]=items[i]
        else:
            dico[items[i]]=alphabet[i]
    return dico


def run_Extract(dataFile,dataT,ensTrans,ensMotif,type,seuilF):
    ensTrans.sort()
    dicoK= getDecodage(dataFile,False)
   # printDico(dicoK)

    r=open("data/sampleTrans.dat","w")
    x=0
    for i in range(len(ensTrans)):
        for j in range(len(dataT[ensTrans[i]])):
            r.write(str(dataT[ensTrans[i]][j])+" ")
        r.write("\n")
    r.close()
    os.system(f"./lcm53 F data/sampleTrans.dat {len(ensTrans)} resLCM.txt >log.txt")

    f=open("resLCM.txt")
    s=f.readline()
    compteurMotif=0
    while s!="":
        compteurMotif=compteurMotif+1
        i=1
        splitS=s.split()
        motifM=[]
        motifK=[]
        while splitS[i]!="]":
            motifM.append(int(splitS[i]))
            motifK.append(dicoK[int(splitS[i])])
            i=i+1
        motifK.sort()
        #print(f"motif={motif}")
        present=False
        for i in range(len(ensMotif)):

            if motifK ==ensMotif[i][2]:

                present=True
        if not(present):
            freq=getFreq(dataT,motifM)
            ensMotif.append([freq,len(motifK),motifK])
        s=f.readline()
    f.close()
    ensMotif.sort(key=lambda x: (-x[0],-x[1],x[2]))


    r=open(f"data/candidates/{dataFile}-{type}-"+str(seuilF)+"d.isc","w")
    r.write("ficfis-1.3\n")
    r.write("mi: numSets="+str(len(ensMotif)))
    maxL=[]
    minS=[]
    for i in ensMotif:
        maxL.append(i[1])
        minS.append(i[0])
    r.write(" minSup="+str(min(minS))+" maxLen="+str(max(maxL)))
    r.write(f" sepRows=0 iscOrder=d patType={type} dbName={dataFile}\n")
    for i in range(len(ensMotif)):
        freqM=ensMotif[i][0]
        longM=ensMotif[i][1]
        motif=ensMotif[i][2]
        r.write(str(longM)+":")
        for j in range(longM):
            r.write(" "+str(motif[j]))
        r.write(" (" + str(freqM)+")\n")
    r.close()
    return ensMotif

def extract(dataFile,dataT,dataI,ensTrans,ensMotif,typeM,minSupp):
    lcm=LCM(min_supp=minSupp)
    #print("extraction des motifs")
    temps=time.time()
    newCandidats=lcm.fit_transform(ensTrans)
    dicoK= getDecodage(dataFile,False)
    nbS=0
    longEns=len(ensMotif)
    for i in newCandidats.index:
        itemset=newCandidats.iat[i,0]
        motifK=[]
        for i in itemset:
            motifK.append(dicoK[int(i)])
        motifK.sort()
        nbI=len(itemset)
        if nbI>1:
            present=False
            for i in range(len(ensMotif)):
                if motifK ==ensMotif[i][2]:
                    present=True
            if not(present):
                couv,freq=getFreq(dataI,itemset)
                #print(f" on a fait le motif {itemset} avec un freq de {freq}")
                ensMotif.append([freq,nbI,motifK,couv])
        else:
            nbS=nbS+1
    r=open(f"Res/{dataFile}/trouve.txt","a")
    r.write(str(len(newCandidats)-nbS))
    r.write(",")
    r.close()
    #print(" on a extraits "+str(len(newCandidats)-nbS))
    r=open(f"Res/{dataFile}/ajout.txt","a")
    r.write(str(len(ensMotif)-longEns))
    r.write(",")
    r.close()
    #print(f" Il y a {len(ensMotif)-longEns} nouveaux motifs")
    ensMotif.sort(key=lambda x: (-x[0],-x[1],x[2]))
    #print(" on encode pour Krimp")
    r=open(f"data/candidates/{dataFile}-{typeM}-"+str(minSupp)+"d.isc","w")
    r.write("ficfis-1.3\n")
    r.write("mi: numSets="+str(len(ensMotif)))
    maxL=[1]
    for i in ensMotif:
        maxL.append(i[1])
    r.write(" minSup="+str(minSupp)+" maxLen="+str(max(maxL)))
    r.write(f" sepRows=0 iscOrder=d patType={typeM} dbName={dataFile}\n")
    for i in range(len(ensMotif)):
        freqM=ensMotif[i][0]
        longM=ensMotif[i][1]
        motif=ensMotif[i][2]
        r.write(str(longM)+":")
        for j in range(longM):
            r.write(" "+str(motif[j]))
        r.write(" (" + str(freqM)+")\n")
    r.close()
    #print(f" ça a pris {time.time()-temps} secondes")

    return ensMotif


def isRectangular(data):
    ref=len(data[0])
    rectangular=True
    for i in data:
        if len(i)!=ref:
            rectangular=False
    return rectangular


def motifDecode(m,decode):
    motifD=[]
    for i in m:
        motifD.append(decode[i])
    return motifD

def motifStrToInt(m):
    motifInt=[]
    tmp=m.split(',')
    for i in tmp:
        motifInt.append(int(i))
    return motifInt
def motifIntToStr(m):
    motifStr=""
    for i in range(len(m)-2):
        motifStr+=f"{m[i]},"
    motifStr+=f"{m[-2]}"
    return motifStr

# def extractCouvKrimp(ensCT,dicoK):
#     f=open("OccurKrimp.txt")
#     s=f.readline()
#     #print(ensCT)
#     couvK=dict()
#     while s!="":
#         splitS=s.split(":")
#         if len(splitS)>1 and splitS[1] !='\n':
#             itemset=splitS[0]
#             couv= splitS[1]
#             tmp=itemset.split()
#             motif=""
#             for i in range(len(tmp)-2):
#                 if tmp[i][0]=="(":
#                     motif=tmp[i].split(")")[1]
#                     motif+=","
#                 else:
#                     motif+=f"{tmp[i]},"
#             motif+=f"{tmp[-2]}"

#             motifInt=motifStrToInt(motif)
#             motifD=motifDecode(motifInt,dicoK)
#             if  motifD in ensCT:
#                 ind=ensCT.index(motifD)
#                 couvInt=[]
#                 tmp=couv.split(",")
#                 for i in tmp:
#                     if i !="\n":
#                         x=int(i)
#                         couvInt.append(x)
#                 couvK[ind]=couvInt
#         s=f.readline()

#     couvI=dict()
#     for i in couvK:
#         for j in couvK[i]:
#             if j in couvI:
#                 couvI[j].append(i)
#             else:
#                 couvI[j]=[i]
#     f.close()
#     #print(couvI)

#     return couvI

def extractCT(dataFile,dataI,indRun):
    ensMotif=[]
    ensCT=[]
    ensCTFC=[]
    sumTot=0
    poids=[]
    dicoK= getDecodage(dataFile,True)
    os.system(f"rm experiments/{dataFile}{indRun}/compress/{dataFile}*/*-0.ct")
    os.system(f"mv experiments/{dataFile}{indRun}/compress/{dataFile}*/*.ct resKrimp.ct")
    f=open(f"resKrimp.ct")
    s=f.readline()
    s=f.readline()
    s=f.readline()
    splitS=s.split()
    while s!="":
        itemset=[]
        motif=[]
        for indIt in range(len(splitS)-1):
            item=splitS[indIt]
            itemset.append(int(item))
            motif.append(dicoK[int(item)])
        itemset.sort()
        tmp=splitS[-1].split(",")[0]
        tmpp=int(tmp.split("(")[1])
        tmp=splitS[-1].split(",")[1]
        freq=int(tmp.split(")")[0])
        #couv,freq=getFreq(dataI,motif)
        if len(motif)>1:
            ensMotif.append([freq,len(itemset),itemset])


        sumTot+=tmpp
        poids.append(tmpp)
        ensCTFC.append(motif)
        ensCT.append([motif,freq])
        s=f.readline()
        splitS=s.split()

    tailleComp=[]
    # tailleDb=0
    for i in range(len(poids)):
        x=poids[i]/sumTot
        tailleComp.append(x)
    #     tailleDb+=x*poids[i]
    # print(f"tailleDb =>{int(round(tailleDb))}")
    os.system(f"mv resKrimp.ct Res/{dataFile}/Codetables/run{indRun}.ct")

    #couvKrimp=extractCouvKrimp(ensCTFC,dicoK)
    #os.system(f"mv OccurKrimp.txt Res/{dataFile}/Codetables/OccurRun{indRun}.txt")
    #return ensMotif,ensCT,couvKrimp,tailleComp,ensCTFC

    return ensMotif,ensCT,tailleComp

def newSampler(data,poidsD,c):
    ensTrans=[]
    #print(sum(poidsD))
    # ensIndex=choices(setIndex,weights=poidsD,k=c)
    #ensIndex=np.random.choice(len(data),c,replace=False,p=poidsD)
    # for i in range(c):
    #     for j in range(i+1,c):
    #         if ensIndex[i]==ensIndex[j]:
    #             print(" c'est la merde")
    #             print(ensIndex)
    wMax=max(poidsD)
    tailleData=len(data)
    ensIndex=[]
    while len(ensIndex)<c:
        accept=False
        while not(accept):
            x=randint(1,tailleData)
            while x in ensIndex:
                x=randint(1,tailleData)
            accept= random()< (poidsD[x-1]/wMax)
        ensIndex.append(x)
    for i in ensIndex:
        tmp=[]
        for j in range(1,len(data[i-1])): ################################# MODIFIER ICI POUR CLASS
            tmp.append(data[i-1][j])
        #ensTrans.append(data[i-1])
        ensTrans.append(tmp)
    return ensTrans,ensIndex


def newCalculPoidsD(data,tailleI,ensCT,tailleD):
    poidsD=[]
    comp=[]
    f=open("OccurKrimp.txt")
    transCover=dict()
    s=f.readline()
    while s!="":
        splitS=s.split(":")
        trans=int(splitS[0])
        itemsets=splitS[1].split()[0].split(",")
        for i in itemsets:
            if trans in transCover:
                transCover[trans].append(int(i))
            else:
                transCover[trans]=[int(i)]
        s=f.readline()
    #print(f"ensCT={ensCT}")
    for i in range(len(data)):
        tailleC=0
        if i in transCover:
            d=[]
            for j in data[i]:
                d.append(j)
            #print(f"d actu={d}")
            #print(f"motif a enlever={couvK[i]}")
            for j in transCover[i]:
                if tailleD[j]>0:
                    tailleC+=-1*math.log2(tailleD[j])
            p=round(tailleC)/round(tailleI[i])
            if p >1:
                comp.append(1)#len(data[i]))
            else:
                comp.append(round(p,4))#*len(data[i]))
        else:
            comp.append(1)
    somC=sum(comp)
    for i in comp:
        poidsD.append((i/somC))

    return poidsD

def calculPoidsD(data,tailleI,ensCT,tailleD):
    poidsD=[]
    comp=[]
    timeCalcul=0
    for i in range(len(data)):
        d=[]
        for j in data[i]:
            d.append(j)
        tailleC=0
        # for j in range(len(ensCT)):
        j=0

        while j <len(ensCT):
            if not(len(ensCT[j][0])>len(d)):
                dedans=True
                for k in ensCT[j][0]:
                    if not(k in d):
                        dedans=False
                if dedans:
                    tmpTime=time.time()
                    tailleC+=-1*math.log2(tailleD[j])
                    for k in ensCT[j][0]:
                        d.remove(k)
                    timeCalcul+=time.time()-tmpTime
            if len(ensCT[j][0])<2  and len(d)==len(data[i]):
                j=len(ensCT)
                tailleC=tailleI[i]
            j+=1
        tmpTime=time.time()
        p=round(tailleC)/round(tailleI[i])
        if p >1:
            comp.append(1)#len(data[i]))
        else:
            comp.append(round(p,4))#*len(data[i]))
        timeCalcul+=time.time()-tmpTime
    tmpTime=time.time()
    somC=sum(comp)
    for i in comp:
        poidsD.append((i/somC))
    timeCalcul+=time.time()-tmpTime
    return poidsD,timeCalcul


def main(dataFile,nbR,nbT,typeM,seuilF,iterK):
    currDir=os.getcwd()
    tempsTotal=time.time()
    data=fetch_file(f'data/{dataFile}.dat',int_values=True)
    poidsT=[1/len(data)]*len(data)
    choixP=[]
    #dict_keys(['n_items', 'avg_transaction_size', 'n_transactions', 'density'])
    infoD=describe(data)
    infoD["isRect"]=isRectangular(data)
    dataI={}

    for t in range(len(data)):
        choixP.append(t)
        for i in data[t]:
            if i in dataI:
                dataI[i].add(t)
            else:
                dataI[i]={t}
    sumInit=0
    for i in dataI:
        sumInit+=len(dataI[i])
    tailleI=[]
    for i in data:
        tailleTmp=0
        for j in i:
            tailleTmp+=-1*math.log2(len(dataI[j])/sumInit)
        tailleI.append(tailleTmp)
   # for i in dataI:
     #   print(len(dataI[i]))
    minSupp=int(math.ceil(nbT*seuilF))

    temps=[]
    tempsCP=[]
    #print(f" on mine dans {dataFile} les motifs {typeM} avec un seuil de {minSupp}")
    i=1
    ensMotif=[]
    ensRemove=[]
    tmpTaille=[]
    tmpTemps=[]
    while i<nbR+1:
       # print(i)
        #print(f" ChoixP={choixP}")
        #print(f" len CoixP={len(choixP)}")
        times=time.time()
        # if len(choixP)>nbT:
        #     sampleT,indexT=cftp_Sampler(data,infoD,nbT,seuilF,choixP)
        #    # print(f"{indexT}")
        #     for indT in indexT:
        #         #print(f" on tente de remove {indT}")
        #         choixP.remove(indT)
        # else:
        #     i=nbR+2
        #     sampleT=choixP
        #     if len(choixP)!=nbT:
        #         nbRestant=nbT-len(choixP)
        #         ajout=sample(len(data),nbRestant)
        #         ok=True
        #         while ok:
        #             ajout=sample(len(data),nbRestant)
        #             ok=False
        #             for ajt in ajout:
        #                 if ajt in choixP:
        #                     ok=True
        #         for ajt in ajout:
        #             sampleT.append(ajt)

        # sampleT,indexT=cftp_Sampler(data,infoD,nbT,seuilF,choixP)
        indexTransaction=run_cftp(dataFile,nbT)
        sampleT=[]
        for k in indexTransaction:
            sampleT.append(data[k])

        # sampleT,indexTransaction=newSampler(data,poidsT,nbT)

        ensMotif=extract(dataFile,data,dataI,sampleT,ensMotif,typeM,minSupp)
        while ensMotif==[]:
            # sampleT,indexTransaction=newSampler(data,poidsT,nbT)
            # sampleT,indexT=cftp_Sampler(data,infoD,nbT,seuilF,choixP)
            indexTransaction=run_cftp(dataFile,nbT)
            sampleT=[]
            for k in indexTransaction:
                sampleT.append(data[k])
            ensMotif=extract(dataFile,data,dataI,sampleT,ensMotif,typeM,minSupp)
        #print(" fin sample")
        temps.append(time.time()-times)
        f=open(f"Res/{dataFile}/Candidates/runT{i}.txt","a")
        f.write(f"{indexTransaction},{sampleT}")
        f.write("\n")
        f.close()
        #ensMotif=run_Extract(dataFile,data,indT,ensMotif,typeM,seuilF)
        if i%iterK==0:# and (i-9)%10==0:
            
            f=open("Krimp/bin/datadir.conf","w")
            f.write(f"dataDir ={currDir}/data/\n")
            f.write(f"expDir ={currDir}/experiments/{dataFile}{i}/\n")
            f.close()
            os.system("Krimp/bin/krimp > result.txt")
            times=time.time()
            #ensMotif,ensCT,couvKrimp,tailleComp,ensCTFC=extractCT(dataFile,dataI,i)
            ensMotif,ensCT,tailleComp=extractCT(dataFile,dataI,i)
            #tempsCP+=time.time()-times
            # if len(ensMotif) >99:
            #     print(f"finsi à {i}")
            #     i=100
           # print(couvKrimp)
            #poidsT=calculPoidsD(data,tailleI,ensCT,tailleComp)
            # print(ensCTFC)
            # for i in ensCT:
            #     ind=ensCTFC.index(i[0])
            #     print(f"{i},{ensCTFC[ind]}")
            # if i<nbR-1:
            #     poidsT=newCalculPoidsD(data,tailleI,ensCT,tailleComp)
                #poidsT,timeCalcul=calculPoidsD(data,tailleI,ensCT,tailleComp)
                #tempsCP+=timeCalcul
            # tempsCP.append(time.time()-times)
            # if i==nbR-1:
            #     os.system(f"mv OccurKrimp.txt Res/{dataFile}/Occur/OccurRun{i}.txt")
            # else:
            #     os.system(f"rm OccurKrimp.txt")
        #print(ensMotif)
            f=open("result.txt")
            stop=True
            while stop :
                s=f.readline()
                splitS=s.split()
                if len(splitS)>1:
                    if splitS[1]=='Time:':
                        stop=False
            tmpTemps.append(float(splitS[6]))
            stop=True
            while stop:
                s=f.readline()
                splitS=s.split()
                if len(splitS)>1:
                    if splitS[1]=='Result:':
                        stop=False
            tmpTaille.append(int(splitS[2].split(',')[5].split(')')[0]))
            f.close()
        os.system(f"mv data/candidates/{dataFile}-{typeM}-{minSupp}d.isc Res/{dataFile}/Candidates/run{i}.isc")
        i=i+1
    tempsTotalF=time.time()-tempsTotal
    #sampleT=cftp_Sampler(data,infoD,nbT,seuilF)
    #ensMotif=extract(dataFile,data,dataI,sampleT,ensMotif,typeM,minSupp)
    # os.system("Krimp/bin/krimp > result.txt")
    # os.system(f"mv data/candidates/{dataFile}-{typeM}-{minSupp}d.isc Res/{dataFile}/Candidates/run{i}.isc")
    # f=open(f"Res/{dataFile}/tempsSample.txt","a")
    # f.write(str(temps))
    # f.write("\n")
    # f.close()
    f=open(f"Res/{dataFile}/tempsTotal.txt","a")
    f.write(str(tempsTotalF))
    f.write("\n")
    f.close()
    # f=open(f"Res/{dataFile}/tempsCP.txt","a")
    # f.write(str(tempsCP))
    # f.write("\n")
    # f.close()
    f=open(f"Res/{dataFile}/tempsKrimp.txt","a")
    f.write(str(tmpTemps))
    f.write("\n")
    f.close()
    f=open(f"Res/{dataFile}/taille.txt","a")
    f.write(str(tmpTaille))
    f.write("\n")
    f.close()




if __name__ == "__main__":
    data=str(sys.argv[1])
    os.system(f"mkdir Res/{data}")	
    os.system(f"mkdir Res/{data}/Candidates")
    os.system(f"mkdir Res/{data}/Codetables")
    os.system(f"touch Res/{data}/ajout.txt")
    os.system(f"touch Res/{data}/trouve.txt")
    nbRun=int(sys.argv[2])
    tailleS=int(sys.argv[3])
    typeMotif="closed"
    iterationK=float(sys.argv[4])
    seuilFreq=0.5
    
    f=open("compress.conf","r")
    r=open("compress.txt","w")
    s=f.readline()
    while(s!=""):
        if s[0]=="i" and s[1] =="s" and s[2]=="c" and s[3]=="N":
            r.write(f"iscName={data}-{typeMotif}-"+str(math.ceil(tailleS*seuilFreq))+"d\n")
        elif s[0]=="d" and s[1] =="a" and s[2]=="t" and s[3]=="a" and s[3]=="T" and s[3]=="y":
            r.write(f"dataType= {dataTypeKrimp[data]}")
        else:
            r.write(s)
        s=f.readline()
    r.close()
    f.close()
    os.system("rm compress.conf")
    f=open("compress.txt","r")
    r=open("compress.conf","w")
    s=f.readline()
    while(s!=""):
       r.write(s)
       s=f.readline()
    r.close()
    f.close()
    os.system("rm compress.txt")
    
    main(data,nbRun,tailleS,typeMotif,seuilFreq,iterationK)

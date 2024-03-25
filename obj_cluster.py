import sys,itertools,os,fileinput
from decimal import *
import sys,fileinput,math,itertools,os
bps = ['A','T','G','C','U','N','M','R','W','S','Y','K','V','H','D','P','B']

def change_ext_gaps(sequence):
	start_pos, end_pos = 0, 0
	for i,bp in enumerate(sequence):
		if bp in bps:
			start_pos = i - 1
			break
	for i,bp in enumerate(sequence[:: - 1]):
		if bp in bps:
			end_pos = len(sequence) - i
			break
	new_sequence = "?" * (start_pos + 1) + sequence[start_pos + 1 : end_pos] + "?" * (len(sequence) - end_pos)
	return new_sequence

def checkmatrix(indict):
	flag = 0
	exp_len = len(indict.keys()[0])
	inchars_exp =  bps + ["?"] + ["-"]
	for v in indict.keys():
		if len(v) != exp_len:
			flag = 1
			break
		for bp in v:
			if bp not in inchars_exp:
				flag = 2
				break
	return flag

def measurepair_distancem(seq1, seq2):
	num_d = 0
	num_s = 0
	diff=[]

	for ele_x, ele_y in zip(seq1,seq2):
		if ele_x!=ele_y:
			if ele_x in ["A","T","G","C"] and ele_y in ["A","T","G","C"]:
					num_d += 1
					num_s+=1
		if ele_x==ele_y:
			if ele_x in ["A","T","G","C","N"]:
					num_s += 1
	perc_d = Decimal(num_d) / Decimal( num_s)
	return perc_d


def measurepair_distance5(seq1, seq2):
	num_d = 0
	num_s = 0
	for i, each in enumerate(seq1):
		if each != seq2[i]:
			if "?" not in [each, seq2[i]]:
				if "N" not in [each, seq2[i]]:
					num_d += 1
		if each == seq2[i]:
			if "?" not in [each,seq2[i]]:
				num_s += 1
	perc_d = Decimal(num_d) / Decimal( num_s)
	return perc_d


def loadseqs(infilename):
	seqdict = {}
	flag = 0
	for line in fileinput.input([infilename]):
		if ">" in line:
			try:
				seqdict[seqID] = change_ext_gaps(sequence.upper())
			except UnboundLocalError:
				pass
			seqID = line.strip().replace(">", "")
			sequence = ''
		else:
			try:
				sequence += line.strip()
			except:
				print "odd"
				flag = 1
				break
	try:
		seqdict[seqID] = change_ext_gaps(sequence.upper())
	except:
		if flag == 0:
			print "odd"
	fileinput.close()
	return seqdict

def rungene(infile,seqdict):
	def countchildren(k,n,listn,masteriddict,terminals):
		if [k]!=masteriddict[k]:
			for item in masteriddict[k]:
				if item not in terminals:
					n,listn=countchildren(item,n,listn,masteriddict,terminals)
				else:
					n=n+1
					listn.append(item)
		else:
			n+=1
		return n,listn
	def buildtree(tree):
		new={}
		for k,v in tree.iteritems():
			newlist=[]
			if len(v)>1:
				v.sort(lambda x,y: cmp(masterdictcount[x], masterdictcount[y]))
				for item in v:
					newlist.append(buildtree({item:masteriddict[item]}))
			new[k]=newlist
		return new

	def buildtree2(tree,master):
		newlist=[]
		for i,v in enumerate(tree):
			if v not in terminals:
				pathdict[v]=master
				clustoutfile=open(os.path.join(master,str(v)+".fasta"),"w")
				for each in masterdictterminals[v]:
					for seq in seqdict.keys():
						if namesdict[each] in seqdict[seq]:
							clustoutfile.write(">"+"Cluster "+str(v)+"; "+namesdict[each]+'\n'+seq+'\n')
				clustoutfile.close()
				newlist.append(buildtree2(masteriddict[v],master))
			else:
				pathdict[namesdict[v]]=master
				clustoutfile=open(os.path.join(master,namesdict[v]+".fasta"),"w")
				for each in masterdictterminals[v]:
					for seq in seqdict.keys():
						if namesdict[each] in seqdict[seq]:
							clustoutfile.write(">"+namesdict[each]+'\n'+seqdict[namesdict[each]]+'\n')
				newlist.append(namesdict[v])
		return newlist,master


	def keychange(tree,n,terminals):
		newlist=[]
		for c in tree:
			new={}
			if isinstance(c, dict):
				for k,v in c.iteritems():
					if k in terminals:
							new={"terminal": str(namesdict[k]), "node_size":str(5), "children": keychange(v,n,terminals), "yv": str(n+masterdictcount[k]/2*10), "xv":str(int(((max(fusepoints.values())-fusepoints[k])*10))) }
					else:
							new={"name": str(k), "node_size":str(5), "fusepoint":str(fusepoints[k]), "curl":"file:///"+os.path.join(infile+"_clusterfastaouts",str(k)+".fasta"), "maindir":infile, "children": keychange(v,n,terminals), "yv": str(n+masterdictcount[k]/2*10), "xv":str(int(((max(fusepoints.values())-fusepoints[k])*10))) }
					n+=masterdictcount[k]*10
			else:
				new={"name": str(c)}
				n+=10
			newlist.append (new)
		return newlist

	pdists=[]
	IDs={}
	dists=[]
	for each in fileinput.input([infile]):
		m=each.split('\t')
		IDs[m[0]]=''
		IDs[m[1]]=''
		dist=float("{0:.1f}".format(float(m[2].strip())*100))
		pdists.append([m[0],m[1],dist])
		if dist not in dists:
			dists.append(dist)
	fileinput.close()
	IDlist=IDs.keys()
	pdists.sort(key=lambda x: x[-1])
	print pdists[0:10]
	dists.sort()
	clusters_by_thresh={-1:{}}

	for dist in dists:
		clusters_by_thresh[dist]={}
	IDdict={}
	namesdict={}
	fusepoints={}
	masteriddict={}
	id_cluster={}
	for i,each in enumerate(IDlist):
		IDdict[each]=int(i+1)
		clusters_by_thresh[-1][i+1]=[i+1]
		namesdict[int(i+1)]=each
		masteriddict[i+1]=[i+1]
		id_cluster[i+1] = i+1
		fusepoints[i+1]=0
	counts={}
	tempfile=open("tempfile.txt",'w')
	tempfile.write(str(namesdict))
	for each in IDlist:
		counts[IDdict[each]]=1
	y=[0]
	masterdicttree=[]

	terminals=range(1,len(masteriddict)+1)
	pos=0
	cutoff=min(dists)
	for i,cutoff in enumerate(dists):
		curr_id_positions=id_cluster
		currmasterdict={}
		for k in masteriddict.keys():
			currmasterdict[k]=masteriddict[k]
		currdict={}
		if cutoff==0:
			for k in clusters_by_thresh[-1].keys():
				currdict[k]=clusters_by_thresh[-1][k]
		else:
			for k in clusters_by_thresh[dists[i-1]].keys():
				currdict[k]=clusters_by_thresh[dists[i-1]][k]
		if len(currdict)==0:
			break
		maxindex=max(id_cluster.values())
		newindices={}
		while pos<len(pdists):
			if pdists[pos][-1]<=cutoff:
				if id_cluster[IDdict[pdists[pos][0]]]!=id_cluster[IDdict[pdists[pos][1]]]:
					newcluster=list(set(currdict[id_cluster[IDdict[pdists[pos][0]]]]) | set(currdict[id_cluster[IDdict[pdists[pos][1]]]]))
					if len(newcluster)>len(currdict[id_cluster[IDdict[pdists[pos][0]]]]):
						maxindex += 1
						if id_cluster[IDdict[pdists[pos][0]]] not in newindices.keys() and id_cluster[IDdict[pdists[pos][1]]] not in newindices.keys():
							newindices[maxindex]=[id_cluster[IDdict[pdists[pos][0]]],id_cluster[IDdict[pdists[pos][1]]]]
							currmasterdict[maxindex]=[id_cluster[IDdict[pdists[pos][0]]],id_cluster[IDdict[pdists[pos][1]]]]
						elif id_cluster[IDdict[pdists[pos][0]]] in newindices.keys() and id_cluster[IDdict[pdists[pos][1]]] in newindices.keys():
							newindices[maxindex]=list(set(newindices[id_cluster[IDdict[pdists[pos][0]]]]) | set(newindices[id_cluster[IDdict[pdists[pos][1]]]]))
							currmasterdict[maxindex]=list(set(newindices[id_cluster[IDdict[pdists[pos][0]]]]) | set(newindices[id_cluster[IDdict[pdists[pos][1]]]]))
						elif id_cluster[IDdict[pdists[pos][0]]] in newindices.keys() and id_cluster[IDdict[pdists[pos][1]]] not in newindices.keys():
							templist=newindices[id_cluster[IDdict[pdists[pos][0]]]]
							templist.append(id_cluster[IDdict[pdists[pos][1]]])
							currmasterdict[maxindex]=templist
							newindices[maxindex]=templist
						elif id_cluster[IDdict[pdists[pos][0]]] not in newindices.keys() and id_cluster[IDdict[pdists[pos][1]]] in newindices.keys():
							templist=newindices[id_cluster[IDdict[pdists[pos][1]]]]
							templist.append(id_cluster[IDdict[pdists[pos][0]]])
							currmasterdict[maxindex]=templist
							newindices[maxindex]=templist
						currdict[maxindex]=newcluster
						currdict.pop(id_cluster[IDdict[pdists[pos][0]]],None)
						currdict.pop(id_cluster[IDdict[pdists[pos][1]]],None)
						for item in  newcluster:
							id_cluster[item] = maxindex
				pos+=1
			else:
				break
		if cutoff==0:
			for id in list(set(currdict.keys())-set(clusters_by_thresh[-1].keys())):
				masteriddict[id]=currmasterdict[id]
				fusepoints[id]=cutoff
			clusters_by_thresh[cutoff]=currdict
			if len(list(set(currdict.keys())-set(clusters_by_thresh[-1].keys())))>0:
				y.append(cutoff)
		else:
			for id in list(set(currdict.keys())-set(clusters_by_thresh[dists[i-1]].keys())):
				masteriddict[id]=currmasterdict[id]
				fusepoints[id]=cutoff
			clusters_by_thresh[cutoff]=currdict
			if len(list(set(currdict.keys())-set(clusters_by_thresh[dists[i-1]].keys())))>0:
				y.append(cutoff)
		if len(currdict)==0:
			break
	y.sort()
	y=y[::-1]
	print y
#	print clusters_by_thresh[0]
	masterdicttree={max(masteriddict.keys()):masteriddict[max(masteriddict.keys())]}

	masterdictcount={}
	for each in masteriddict.keys():
		masterdictcount[each]=0

	masterdictterminals={}
	for each in masteriddict.keys():
		masterdictterminals[each]=0
	print "masterdict is",masteriddict
	for each in masteriddict.keys():
		cchildren,listchildren=countchildren(each,0,[],masteriddict,terminals)
		masterdictcount[each]= cchildren
		masterdictterminals[each]=listchildren

	clustfile=open(infile+"_clusters",'w')

	for each in masteriddict.keys():
		keylist=masterdictterminals[each]
		keylist.sort()
		clustfile.write(str(each)+"\t"+str(masterdictcount[each])+"\t"+str(fusepoints[each])+'\t'+str(keylist)+'\n')
	masterdicttree1=buildtree(masterdicttree)
	pathdict={}
	if os.path.isdir(infile+"_clusterfastaouts")==False:
		os.mkdir(infile+"_clusterfastaouts")
	masterdicttree2=buildtree2([max(masteriddict.keys())],infile+"_clusterfastaouts")
	if os.path.isdir(infile+"_threshfastaouts")==False:
		os.mkdir(infile+"_threshfastaouts")
	print terminals
	for thresh in clusters_by_thresh.keys():
		if thresh in y:
			clustoutfile=open(os.path.join(infile+"_threshfastaouts",str(thresh)+".fasta"),'w')
			for k in clusters_by_thresh[thresh].keys():
				if k in terminals:
					for seq in seqdict.keys():
						if namesdict[k] in seqdict[seq]:
							clustoutfile.write(">"+"Singleton "+str(k)+" "+namesdict[k]+'\n'+seqdict[namesdict[k]]+'\n')
				else:
					for each in masterdictterminals[k]:
						for seq in seqdict.keys():
							if namesdict[k] in seqdict[seq]:
								clustoutfile.write(">"+"Cluster "+str(k)+"; "+namesdict[each]+'\n'+seqdict[namesdict[each]]+'\n')
			clustoutfile.close()


	testfile=open(infile+"_test",'w')
	testfile.write(str(masterdicttree2))
	n=0
	masterdicttree2=keychange([masterdicttree1],n,terminals)

	outfile=open(infile+"_clusterlist",'w')
	clusters_by_thresh_namedict={}
	for each in clusters_by_thresh.keys():
		clusters_by_thresh_namedict[each]={}

	clusters_by_thresh2={}
	for each in clusters_by_thresh.keys():
		clusters_by_thresh2[each]={}

	for each in clusters_by_thresh.keys():
		for ind in clusters_by_thresh[each]:
			newlist=[]
			newlist2=[]
			for item in clusters_by_thresh[each][ind]:
				if len(masterdictterminals[item])>0:
					newlist+=masterdictterminals[item]
					for v in masterdictterminals[item]:
						newlist2.append(namesdict[v])
				else:
					newlist.append(item)
					newlist2.append(namesdict[item])
			clusters_by_thresh2[each][ind]=newlist
			clusters_by_thresh_namedict[each][ind]=newlist2
	for each in dists:
		outfile.write(str(each)+ " : "+ str(len(clusters_by_thresh[each])) +" : "+str(clusters_by_thresh_namedict[each])+"\n")
	return fusepoints, y,masterdicttree2

class measuredist:
    def __init__(self, seqdict, pdistfile):
        self.seqdict = seqdict
        self.pdistfile = pdistfile
        #		dict2=self.seqdict.keys()[len(self.seqdict.keys()):]
        self.n = 1
        lines=[]
        for seq in self.seqdict.keys():
            seqlist = self.seqdict[seq]
            del self.seqdict[seq]
            for remseq in self.seqdict.keys():
                dist = measurepair_distancem(seq, remseq)
                incombs = list(itertools.product(seqlist, self.seqdict[remseq]))
                for incomb in incombs:
                    lines.append(incomb[0] + '\t' + incomb[1] + '\t' + str(dist) + '\n')
                self.pdistfile.write("".join(lines))
                lines = []
        self.pdistfile.close()
infilename=sys.argv[1]
seqdict = {}
flag = 0
seqcounts = 0 
for line in fileinput.input([infilename]):
	if ">" in line:
		seqcounts += 1
		try:
			seqdict[change_ext_gaps(sequence.upper())] = []
		except NameError:
			pass
		sequence = ''
	else:
		try:
			sequence += line.strip()
		except:
			flag = 1
			break
try:
	seqdict[change_ext_gaps(sequence.upper())] = []
except:
	if flag == 0:
		print "Invalid Fasta"
fileinput.close()
for line in fileinput.input([infilename]):
	if ">" in line:
		try:
			seqdict[change_ext_gaps(sequence.upper())].append(seqID)
		except NameError:
			pass
		seqID = line.strip().replace(">", "")
		sequence = ''
	else:
		try:
			sequence += line.strip()
		except:
			print "Invalid Fasta"
			flag = 1
			break
try:
	seqdict[change_ext_gaps(sequence.upper())].append(seqID)
except:
	if flag == 0:
		print "Invalid Fasta"
fileinput.close()


combs = list(itertools.combinations(seqdict.keys(), 2))
pdistfile = open(infilename + "_pmatrix", 'w')
n = 1
#       n2 = 1
lines = []
for seq in seqdict.keys():
	if len(seqdict[seq]) > 1:
		incomb = list(itertools.combinations(seqdict[seq], 2))
		for comb in incomb:
			lines.append(comb[0] + '\t' + comb[1] + '\t' + "0.0" + '\n')
			if len(lines) % 100000 == 0:
				pdistfile.write("".join(lines))
				lines = []
pdistfile.write("".join(lines))
lines = []

mymeasuredist = measuredist(seqdict, pdistfile)
rungene(infilename + "_pmatrix", seqdict)




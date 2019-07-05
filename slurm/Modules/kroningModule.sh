#!/bin/bash 
#SBATCH -J kroningModule
#SBATCH -o kroningModule.out
#SBATCH -e kroningModule.err
#SBATCH -t KRONINGTIME
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=USERMAIL
#SBATCH --partition=KRONINGPARTITIONS
#SBATCH -N 1

if [[ "$@" =~ "--debug" ]]; then
        set -ex
else
        set -e
fi  

source CHAINREACTION


function prokkagff2bedFunction {
	infile=$1

	if [ "$infile" == "" ] ; then
   		echo "Usage: prokkagff2bed.sh <PROKKA gff file>"
   		exit 0
	fi

	grep -v "#" $infile | grep "ID=" | cut -f1 -d ';' | sed 's/ID=//g' |awk -F"\t" '{if($3=="CDS")print}' |awk -v OFS='\t' '{print $1, $4-1, $5,$9}' |awk '{if(NF==4)print}'
}
function getcovFunction {
	echo '#!/usr/bin/env python
###################################
import gzip as gz, csv, sys

def GetFilelist(i):
    import os
    ## If user supplies an existing hist file, use it
    if os.path.exists(i) and not "/dev/" in i: return [i]
    ## Otherwise, read the files
    files = []
    hin = open(i, '"'r'"')
    for line in hin: files.append(line.rstrip())
    hin.close()
    return files

def WriteCov(samples, genes, houtcsv):
    cover = {}
    houtcsv.writerow(["Gene"]+sorted(samples.keys()))    
    for gene in genes.keys():
        out = [gene]
        for sample in sorted(samples.keys()):
            try: c = samples[sample][gene]
            except KeyError: c = 0.0
            try: cover[sample].append(c)
            except KeyError: cover[sample] = [c]
            out.append(c)
        houtcsv.writerow(out)
    return cover

def WriteStat(cover, houtcsv):
    import numpy as np
    houtcsv.writerow(["Sample","MeanCoverage","MedianCoverage","MaxCoverage","MinCoverage"])
    for sample,l in cover.iteritems():
        median_cov = np.median(l)
        max_cov = np.max(l)
        min_cov = np.min(l)
        mean_cov = np.mean(l)
        houtcsv.writerow([sample,mean_cov,median_cov,max_cov,min_cov])

def CalcCov(files):
    samples = {}
    genes = {}
    cover = {}
    hout = sys.stdout
    houtcsv = csv.writer(hout, delimiter = '"'\t'"')
    for f in files:
        if ".gz" in f: hin = gz.open(f)
        else: hin = open(f)
        sample=(f.split("/")[-1]).split(".")[0]
        samples[sample] = {}
        cover[sample] = []
        hincsv = csv.reader(hin, delimiter = '"'\t'"')
        gene = ""
        cov = 0.0
        for i,row in enumerate(hincsv):
            #if len(row)<8: break
            if gene != row[3]:
                if len(files)==1:
                    if i==0:houtcsv.writerow(["Gene",sample])
                    else: houtcsv.writerow([gene,cov])
                else: samples[sample][gene] = cov
                gene = row[3]
                cov = 0.0
            genes[gene] = ""
            depth = float(row[-4])
            depth_f = float(row[-1])
            cov += depth*depth_f
        hin.close()
        ## Handle last gene in the file
        if len(files)==1: houtcsv.writerow([gene,cov])
        else: samples[sample][gene] = cov
    if len(files)==1: sys.exit()
    cover = WriteCov(samples, genes, houtcsv)
    hout.close()
    hout = sys.stderr
    houtcsv = csv.writer(hout, delimiter = '"'\t'"')
    WriteStat(cover, houtcsv)
    hout.close()

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infiles", type=str,
            help="Infile histogram(s) from bedtools coverage")
    args = parser.parse_args()
    if not args.infiles: sys.exit(parser.print_help())
    files = GetFilelist(args.infiles)
    CalcCov(files)

if __name__ == "__main__": main()' > get_coverage_for_genes.py
}

function minpathFunction {

	echo '#!/usr/bin/env python
# MinPath
# developed by Yuzhen Ye (yye@indiana.edu)
# Indiana University, Bloomington

# whats new in version 1.2 (released on Oct 19, 2010)
#  MinPath1.2 works on any pathway system
#  what you need: a pathway-function mapping file (e.g, data/ec2path), and your input file

# whats new in version 1.1 (released on July 31, 2010) 
#  * numpy independent
#  * add detailed output of ko assignments to each pathway

import sys
import os
import re
import operator
import math

minpath = os.environ.get('"'MinPath'"')
path0 = "/home/yye/Pathways/MinPath"
if minpath or os.path.exists(path0):
        if os.path.exists(path0):
                minpath = path0
else:
        sys.exit("Environment variable MinPath not set")

keggPath0, seedPath0, mapPath0, glpsol0 = minpath + "/data", minpath + "/data", minpath + "/data", minpath + "/glpk-4.6/examples/glpsol"

def intmatrix(dim1, dim2):
	mat = []
	for i in range(dim1):
		tmp = [0] * dim2		
		mat.append(tmp)

	return mat
	
class MinPath:

	def __init__(self, whichdb = "KEGG", pathdir="", getgene = False, givenspe = "", mapfile = ""):
		self.whichDB = whichdb
		self.speID = givenspe
		global keggPath0
		global seedPath0
		global mapPath0
		if pathdir == "":
			if whichdb == "KEGG": 
				self.dataDir = keggPath0
			elif whichdb == "SEED":
				self.dataDir = seedPath0
			else:
				self.dataDir = mapPath0
		else:
			self.dataDir = pathdir 
		self.famTot = 0
		self.pathwayTot = 0
		self.famList = []
		self.famID = []
		self.famName = []
		self.famCount = []
		self.fam2Path = []
		self.pathList = []
		self.pathName = []
		self.path2Fam = []
		self.path2FamUni = []
		self.orgList = []
		self.orgGene2Fam = []
		self.orgGeneList = []
		self.famMapped = []
		self.pathMapped = []
	
		if whichdb == "SEED":
			print "now get SEED"
			fig2ssfile = self.dataDir + "/figfam_subsystem.dat"
			self.ReadFigSubsytem(fig2ssfile)
		elif whichdb == "KEGG":
			pathfile = self.dataDir + "/map_title.tab"	
			self.ReadKEGGPath(pathfile)

			kofile = self.dataDir + "/ko" 
			self.ReadKO(kofile, getgene, givenspe) #KO: KEGG family
		else:
			if os.path.exists(mapfile):
				print "mapfile", mapfile
				self.ReadAnyMap(mapfile)
			elif os.path.exists(self.dataDir + "/" + os.path.basename(mapfile)):
				print "file: ", self.dataDir + "/" + os.path.basename(mapfile)
				self.ReadAnyMap(self.dataDir + "/" + os.path.basename(mapfile))
			else:
				sys.exit("file " + mapfile + " not found")

		self.CheckUniqueFam()

        def GetPathList(self):
		return self.pathList;

        def GetPathName(self):
		return self.pathName;

	def ReadAnyMap(self, mapfile):
		try:
			file = open(mapfile, "r")
		except IOError:
			print "open file %s error" % mapfile

		for aline in file:
			if aline[0] == '"'#'"':
				continue
			subs = aline.strip().split()
			if len(subs) < 2:
				continue
			path, fam = subs[0], subs[1]
			if path in self.pathName:
				pathidx = self.pathName.index(path)
			else:
				pathidx = len(self.pathList)
				self.pathList.append(str(pathidx + 1))
				self.pathName.append(path)
				self.path2Fam.append([])
			if fam in self.famID:
				famidx = self.famID.index(fam)
			else:
				famidx = len(self.famList)
				self.famList.append(str(famidx + 1))
				self.famName.append(fam)
				self.famID.append(fam)
				self.fam2Path.append([])
			if famidx not in self.path2Fam[pathidx]:
				self.path2Fam[pathidx].append(famidx)
				self.fam2Path[famidx].append(pathidx)
				#same path & fun could be listed multiple times

		self.famTot = len(self.famList)
		self.pathTot = len(self.pathList)		
		print "total family", self.famTot, " pathway", self.pathTot
		#for idx in range(self.pathTot):
		#	print "path", self.pathName[idx], " include fam", len(self.path2Fam[idx])

	#read SEED figfam to subsytem mapping
	def ReadFigSubsytem(self, fig2ssfile):
		print "fig2ssfile=%s" % fig2ssfile
		try:
			file = open(fig2ssfile, "r")
		except IOError:
			print "open file %s error" % fig2ssfile
		for aline in file:
			aline = aline.strip()
			#note: in seed, a "family" (or "function") can have multiple FIG "sub"families
			#the families and the subsytems have names, but not codes
			#(subfam, fam, path) = aline.split('"'\t'"')
			cols = aline.split('"'\t'"')
			if len(cols) < 3: 
				continue  # skip the families that have NO subsystems assignment
			fam = cols[1]
			if fam in self.famName:
				famidx = self.famName.index(fam)	
			else:
				famidx = len(self.famName)
				self.famList.append("F" + str(famidx + 1))
				self.famID.append(cols[0])
				self.famName.append(fam)
				row = []
				self.fam2Path.append(row)
			for path in cols[2:]: 
				# note: a line could have multiple subsytems??
				if path[:4] == "CBSS":
					continue  # skip the clustering-based subsystems
				if path in self.pathName:
					pathidx = self.pathName.index(path)
				else:
					pathidx = len(self.pathName)
					self.pathList.append("S" + str(pathidx + 1))
					self.pathName.append(path)
					row = []
					self.path2Fam.append(row)

				if pathidx not in self.fam2Path[famidx]:
					self.fam2Path[famidx].append(pathidx)
				if famidx not in self.path2Fam[pathidx]:
					self.path2Fam[pathidx].append(famidx)
		file.close()
		self.pathTot = len(self.pathList)
		self.famTot = len(self.famList)
		print "total SEED subsystem=%d" % self.pathTot
		print "total SEED functions(families)=%d" % self.famTot

	#read KEGG pathways from ~/pathway/map-title.tab
	def ReadKEGGPath(self, pathfile):
		try:
			file = open(pathfile, "r")	
		except IOError:
			print "open file %s error" % pathfile
			sys.exit()

		for aline in file:
			aline = aline[:-1]
			row = aline.split("\t")
			self.pathList.append(row[0])
			self.pathName.append(row[1])
			row = []
			self.path2Fam.append(row)
		file.close()

		self.pathTot = len(self.pathList)
		print "total KEGG pathway=%d" % self.pathTot

	#read the KEGG families (KO) from ~/genes/ko
	def ReadKO(self, kofile = "", ifreadgene = False, ifgivenspe = ""):
		try:
			file = open(kofile, "r")	
		except IOError:
			print "open file %s error" % kofile
		print "read kofile=%s" % kofile
		for aline in file:
			aline = aline.strip()
			m = re.search(r'"'^ENTRY\s+(?P<ko>[\S]+)'"', aline)
			if m:
				self.famList.append(m.group('"'ko'"'))
				koidx = len(self.famList) - 1
				row = []
				self.fam2Path.append(row)
				for aline in file:
					name = aline[12:-1]
					self.famName.append(name)
					break
				#print "ko=%d %s name=%s" % (koidx, m.group('"'ko'"'), name)
				continue
			m = re.search(r'"'\[PATH:ko(?P<path>[^\]]+)'"', aline)
			if m:
				thispath = m.group('"'path'"')
				if thispath in self.pathList:
					pathidx = self.pathList.index(thispath)
					#print "thispath=%s %d idx=%d" % (thispath, pathidx, len(self.famList)) 
					self.fam2Path[koidx].append(pathidx)
					self.path2Fam[pathidx].append(koidx)
				continue

			#this information is only needed for KEGG pathway inference (blast-result to ko assignment)
			#reading this information is time-consuming
			if not ifreadgene:
				continue
       	 		m = re.search(r'"'^GENES'"', aline)
			if not m:
				continue
			#raw_input(" continue with reading gene...")
			lines = [aline]
			for aline in file:
				if aline[:3] == "///":
					break
				lines.append(aline)
			ifvalid = False
			for aline in lines:
                                #print '"'now check line'"', aline
				if aline[15] == '"':'"':
					org = aline[12:15].lower()
                                        #print "org", org, "ifgivenspe", ifgivenspe
					if ifgivenspe != "" and org != ifgivenspe:
						ifvalid = False 
					elif (ifgivenspe != "" and org == ifgivenspe) or (ifgivenspe == ""):
						ifvalid = True 
						if org in self.orgList:
							orgidx = self.orgList.index(org)
							else:
							orgidx = len(self.orgList)
							self.orgList.append(org)
							row = []
							self.orgGeneList.append(row)
							row = []
							self.orgGene2Fam.append(row)
				if not ifvalid:
					continue
				#else: the org is the same
                       		info = aline[17:]
                        	info2 = re.sub(r'"'\([^\)]+\)'"', "", info)
	                        items = info2.split()
        	                for aquery in items:
					#thisgene = org + ":" + aquery
					thisgene = aquery
					#print "thisgene = ", thisgene
					if thisgene in self.orgGeneList[orgidx]:
						geneidx = self.orgGeneList[orgidx].index(thisgene)
					else:
						geneidx = len(self.orgGeneList[orgidx])
						self.orgGeneList[orgidx].append(thisgene)
						row = []
						self.orgGene2Fam[orgidx].append(row)
					self.orgGene2Fam[orgidx][geneidx].append(koidx)

		file.close()
		self.famTot = len(self.famList)
		self.famID = self.famList
		self.orgTot = len(self.orgList)
		print "total KEGG fam=%d" % (self.famTot)
		if ifreadgene:
			print "total organisms involved =%d" % (self.orgTot)
			one2one = 0
			one2mul = 0
			for i in range(self.orgTot):
				#print " org %s has %d genes assigned to KO" % (self.orgList[i], len(self.orgGeneList[i]))
				for j in range(len(self.orgGeneList[i])):
					if len(self.orgGene2Fam[i][j]) > 1:
						one2mul += 1
					else:
						one2one += 1
			print "total %d genes matched to multiple KO; %d matched to single KO" % (one2mul, one2one)

	#calcualte the uniqueness of Fam (KEGG families or SEED families)
	def CheckUniqueFam(self):
		for i in range(self.famTot):
			if len(self.fam2Path[i]) > 10:
				print "fam=%d %s [%s] map-to-path=%d" % (i, self.famList[i], self.famName[i], len(self.fam2Path[i]))			
				for p in self.fam2Path[i]:
					print "   path-%d[%s %s]" % (p, self.pathList[p], self.pathName[p])
		ubifam = 0
		for i in range(self.famTot):
			if len(self.fam2Path[i]) > 1:
				ubifam += 1
		print "total families that mapped to more than one pathway = %d" % ubifam
					
		for p in range(self.pathTot):
			row = []
			totfam = len(self.path2Fam[p])
			uniquefam = 0 
			print "pathway-%d[%s; %s] fam=%d" % (p, self.pathList[p], self.pathName[p], totfam)
			for k in range(totfam):
				fam = self.path2Fam[p][k]
				if len(self.fam2Path[fam]) == 1:
					uniquefam = uniquefam + 1 
					row.append(fam)
					print "   unique fam-%d %s %s" % (fam, self.famList[fam], self.famName[fam])
			self.path2FamUni.append(row)
			print ">>>pathway-%d[%s; %s] fam=%d unique-fam=%d" % (p, self.pathList[p], self.pathName[p], totfam, uniquefam)
	def GetPath2FamUniMapped(self, apath, what):
		if apath in self.pathList:
			p = self.pathList.index(apath)	
			famlist = []
			for fam in self.path2FamUni[p]:
				if fam in self.famMapped:
					if what == "name":
						famlist.append(self.famList[fam])	
					else:
						famlist.append(fam)	
			return famlist
		else:
			return []

	def GetPath2FamMapped(self, apath, what):
		if apath in self.pathList:
			p = self.pathList.index(apath)	
			famlist = []
			#print "apath=", apath, "total fam", len(self.path2Fam[p])
			for fam in self.path2Fam[p]:
				#print "  >>>check fam", fam, "famMapped-total", len(self.famMapped)
				if fam in self.famMapped:
					if what == "name":
						famlist.append(self.famList[fam])
					else:
						famlist.append(fam)
					#print "   >>>Found in the famMapped"
			return famlist
		else:
			return []

	#return the index of families
	def GetPath2FamUni(self, apath):
		if apath in self.pathList:
			p = self.pathList.index(apath)	
			return self.path2FamUni[p]
		else:
			return []
	#return the index of families
	def GetPath2Fam(self, apath):
		if apath in self.pathList:
			p = self.pathList.index(apath)	
			return self.path2Fam[p]
		else:
			return []

	#note: this function only works when the givenspe is defined (see KEGG2html.py)
	def OrthMapBasedOnKO(self, spe=""):
		if spe == "":
			return
		if spe not in self.orgList:
			return
		self.famMapped = []
		s = self.orgList.index(spe) 
		print "spe", spe, "total gene", len(self.orgGeneList[s])
		for g in range(len(self.orgGeneList[s])):
			for fam in self.orgGene2Fam[s][g]:
				if fam not in self.famMapped:
					self.famMapped.append(fam)
		print "total family=", len(self.famList), "total mapped=", len(self.famMapped)
		#raw_input("type enter to continue")

	def OrthMap(self, famidxlist=[], famnamelist=[], famcount = []):
		if len(famidxlist) == 0 and len(famnamelist) == 0:
			return []
		if famcount:
			self.famCount = [0] * self.famTot
		#get only the orthologs that are mappped to a pathway
		if len(famidxlist) > 0:
			famlist = famidxlist
			famref = self.famList
		else:
			famlist = famnamelist
			famref = self.famName
		orthtot0 = len(famlist)
		orthtotfind = 0
		orthtotmap = 0 
		self.famMapped = []
		pathmap = [0] * self.pathTot
		for i in range(len(famlist)):
			orth = famlist[i]
			if orth in famref:
				idx = famref.index(orth)
				if famcount:
					self.famCount[idx] = famcount[i] 
				orthtotfind += 1
				if len(self.fam2Path[idx]) > 0:
					self.famMapped.append(idx) 
					orthtotmap += 1
					for path in self.fam2Path[idx]:
						pathmap[path] = 1
		print "original ortholog=%d found-in-the-fam-list=%d found-in-the-fam-mapped-to-pathway=%d" % (orthtot0, orthtotfind, orthtotmap)
		self.pathMapped = []
		for p in range(self.pathTot):
			if pathmap[p]:
				self.pathMapped.append(p)

	#assign orthologs to pathways
	#strategy 1: first assign unique ones -- then other orthologs 
	def Orth2PathUni(self, famidxlist=[], famnamelist=[], famcount=[]):
		self.OrthMap(famidxlist = famidxlist, famnamelist = famnamelist, famcount=famcount)
		orthtotmap = len(self.famMapped)
		#sort the orthologs based on their "uniqueness"
		fam2path = [] 
		for i in self.famMapped:
			if len(self.fam2Path[i]) == 0:
				continue
			fam2path.append((len(self.fam2Path[i]), i)) 
                #print "fam2path = ", len(fam2path)
                #raw_input("type enter to continue..")

		fam2pathsorted = sorted(fam2path, key=operator.itemgetter(0))
		pathsort = map(operator.itemgetter(0), fam2pathsorted)
		famsort = map(operator.itemgetter(1), fam2pathsorted)

		#for i in range(len(famsort)):
		#	print "fam0 %d fam %d [%s; %s] topath %d" % (i, famsort[i], self.famList[famsort[i]], self.famName[famsort[i]], pathsort[i])

		#pathfam & pathfam0: the number of fam families assigned to each pathway
		pathfam = [0] * self.pathTot
		pathfam0 = [0] * self.pathTot

		#pathfam0: the number of fam assigned to each pathway (using all multiple assignments)
		for fam in famsort:
			for p in self.fam2Path[fam]:
				pathfam0[p] = pathfam0[p] + 1

		#pathfam: the number of fam assigned to each pathway considering the "uniqueness" of fam to each pathway
		maxhit = pathsort[-1] 
		print "the maximum number of pathways a family is assigned to=%d" % maxhit
		hit = 1
		unassigned = len(famsort)
		beg = 0
		annpath = []
		famassign = [-1] * orthtotmap
		while (hit <= maxhit) and (unassigned > 0):
			#print "check hit=%d beg=%d" % (hit, beg)
			for k in range(beg, orthtotmap):
				if famassign[k] != -1:
					continue
				if pathsort[k] > hit:
					break 
				fam = famsort[k]
				maxsaturate = -1 
				maxsaturate_p = 0
				#print "check k=%d fam=%d %s %s" % (k, fam, self.famList[fam], self.famName[fam])
				for p in self.fam2Path[fam]:
					saturate = 1.0 * pathfam[p] / len(self.fam2Path[fam])
					if(saturate > maxsaturate):
						maxsaturate = saturate
						maxsaturate_p = p
				famassign[k] = maxsaturate_p
				pathfam[maxsaturate_p] = pathfam[maxsaturate_p] + 1
			beg = k

			unassigned = 0
			for k in range(orthtotmap):
				if famassign[k] == -1:
					unassigned = unassigned + 1
			#print "try hit=%d unassigned=%d (tot=%d)" % (hit, unassigned, orthtotmap)
			#raw_input()
			hit = hit + 1

		annpath0 = 0
		annpath = 0
		self.pathMappedOpt = []
		for k in range(self.pathTot):
			if pathfam0[k] != 0:
				annpath0 = annpath0 + 1
			if pathfam[k] != 0:
				annpath = annpath + 1
				self.pathMappedOpt.append(k)

		print "total pathway %d (%d) is found, compared to %d (%d)" % (annpath, len(self.pathMappedOpt), annpath0, len(self.pathMapped))
		print "%-50s %s\t%s" % ("#pathway", "fam-assigned(all)", "fam-assigned(unique)[weight]")
		for p in range(self.pathTot):
			if pathfam0[p] == 0 and pathfam[p] == 0:
				continue
			tmp = self.pathList[p] + "[" + self.pathName[p] + "]"
			weight = .0
			for k in range(orthtotmap):
				ks = famsort[k]
				if famassign[k] != p:
					continue
				weight = weight + 1.0 / len(self.fam2Path[ks])
			print "%-50s %d\t%d[%.1f]" % (tmp, pathfam0[p], pathfam[p], weight)
		return self.pathMappedOpt

	#Parsinomy approach to pathway inference
	def Orth2PathMin(self, famidxlist=[], famnamelist=[], famcount=[], mpsfile="test.mps", glpsol=""):
		# write mps file (the input for glpsol, the integer programming package)
		print "now write mps file.."
		self.WriteMps(famidxlist=famidxlist, famnamelist=famnamelist, famcount=famcount, mpsfile=mpsfile)

		# run glpsol
		global glpsol0
		if glpsol == "":
			glpsol = glpsol0
		lpout = mpsfile + ".LPout"
		command = glpsol + " " + mpsfile + " -o " + lpout
		print "now run command = %s" % command
		os.system(command)

		# check the result
		self.GetLPOut(lpout)

		return self.pathMappedOpt

	#output mps file for integer programming (most parsinomy pathway inference)
	def WriteMps(self, famidxlist=[], famnamelist=[], famcount=[], mpsfile="test.mps"):
		try:
			file = open(mpsfile, "w")
		except IOError:
			print "open file %s error" % mpsfile
			sys.exit()
	        str = "%-14s%s\n" % ("NAME", "PATH")
       	 	file.write(str)
	
		self.OrthMap(famidxlist=famidxlist, famnamelist=famnamelist, famcount=famcount)
		orthtotmap = len(self.famMapped)

        	#write ROWS
		file.write("ROWS\n")
		file.write(" N  NUM\n")
		for orth in self.famMapped:
                	str = " G  F%s\n" % self.famList[orth];
			#note 1: use a different idx of family for mps file
                	#note 2: when use E, there is no feasible solution
                	#G>=1, mean each family has to be assigned to at least one pathway 
			file.write(str)

        	#write COLUMNS
		#the same column (pathway) needs to be organized in the same block 
        	file.write("COLUMNS\n")
		pathvalid = [0] * self.pathTot
		for p in self.pathMapped:
			#note: use a different pathway idx in mps
			pathname = "P%s" % self.pathList[p]
		        str = "    %-10s%-10s%10d\n" %(pathname, "NUM", 1)
			file.write(str)
			for orth in self.famMapped:
       	                 	famname = "F%s" % self.famList[orth];
				for path in self.fam2Path[orth]:
					if path == p:
		                        	str = "    %-10s%-10s%10d\n" %(pathname, famname, 1)
						file.write(str)	
						pathvalid[p] = 1

        	#write RHS
		file.write("RHS\n")
		for orth in self.famMapped:
                	famname = "F%s" % self.famList[orth]
                	str = "    %-10s%-10s%10.1f\n" % ("RHS1", famname, 1.0);
			file.write(str)

        	#write bounds
		file.write("BOUNDS\n")
		for p in range(self.pathTot):
			if pathvalid[p]:
				path = self.pathList[p]
	                	pathname = "P%s" % path;
       	         		str = " BV %-10s%-10s\n" %("BND1", pathname);
                		#all variants are binary (1 keep the pathway; 0 pathway not necessary)
				file.write(str)

        	file.write("ENDATA\n")
		file.close()		

        	print "End of PrintMPS"
		
	def GetLPOut(self, lpoutfile="test.mps.LPout"):
		try:
			file = open(lpoutfile, "r")
		except IOError:
			print "open file %s error" % lpoutfile 
			sys.exit(1)

		keeppath = []
		for aline in file:
			aline = aline.strip()
			cols = aline.split()
			if len(cols) < 2:
				continue
			if cols[0] == "Columns:":
				Columns = int(cols[1])
			elif cols[0] == "Objective:":
				MINimum = int(cols[3])
			elif cols[0] == "No." and cols[1] == "Column":	
				for aline2 in file:
					if aline2[0] == '"'-'"': 
						continue
					aline2 = aline2.strip()
					cols2 = aline2.split()
					if len(cols2) < 1:
						break
					if cols2[3] == '"'1'"':
						keeppath.append(cols2[1][1:])
						#check with WriteMps: the pathway idx used in mps is to add "P" before the pathList
		file.close()
		if len(keeppath) != MINimum:
			print "reading %s error: minimum %d read %d" % (lpoutfile, MINimum, len(keeppath)) 
			sys.exit()

		self.pathMappedOpt = []
		for path in keeppath:
			if path not in self.pathList:
				print "Error: unknown pathList %s" % path
				sys.exit()
			pathidx = self.pathList.index(path)
			self.pathMappedOpt.append(pathidx)	

		print "total pathways mappd: before inference %d, after inference %d" % (len(self.pathMapped), len(self.pathMappedOpt))

	#add the pathways with many functions annotated, even they were considered as redundant ones!
	def PopulatePath(self, pathmapped = [], par = 0.7):
		famvalid = [0] * self.famTot
		for fam in self.famMapped:
			famvalid[fam] = 1

		addpath = 0
		for p in range(self.pathTot):
			if len(self.path2Fam[p]) == 0:
				continue
			if not p in pathmapped:
				add = 0
				for f in self.path2Fam[p]:
					if famvalid[f] == 1:
						add += 1
				#pathways with most functions annotated should be added back -- even it is a redundant one
				print "pathway", p, self.pathList[p], self.pathName[p], "path2fam", len(self.path2Fam[p]), " real-family", add
				if add >= len(self.path2Fam[p]) * par:
					pathmapped.append(p)
					addpath += 1
					print "this pathway is added back"
				#else:
				#	print "this pathway does not have enough functions"
				#raw_input("type enter to continue")

		print "added pathway =", addpath

	#remove the pathways with too few functions annotated (e.g., 2, use par), 
	#even when their associated families are annotated(but NOT the unique ones)
        #not unique families assigned to this pathway? ubiquitous families have to be assigned to at least one of the pathways, right
	def RemoveSparsePath(self, pathmapped = [], par = 2):			
		famvalid = [0] * self.famTot
		for fam in self.famMapped:
			famvalid[fam] = 1

		delpath = 0
		for p in range(self.pathTot):
			if len(self.path2Fam[p]) == 0:
				continue
			if p in pathmapped:
				add = 0
				uni = 0
				for f in self.path2Fam[p]:
					if famvalid[f] == 1:
						add += 1
						if len(self.fam2Path[f]) == 1: 
							uni += 1
				#pathways with few functions annotated are removed
				if uni == 0 and add <= par:
					pathmapped.remove(p)
					delpath += 1
					print "pathway", p, self.pathList[p], self.pathName[p], "path2fam", len(self.path2Fam[p]), " real-family", add, " is removed from the list!!"
					#raw_input("type enter to continue")

		print "deleted pathway =", delpath

	def DiffPathMap(self, maps, tags):
		maps.insert(0, self.pathMapped)
		tags.insert(0, "Ori")
		mapnum = len(maps)
		pathvalid = intmatrix(self.pathTot, mapnum)
		for m in range(mapnum):
			for p in maps[m]:
				pathvalid[p][m] = 1
		print "#Summary for the pathway inference"
		#print description line
		str = "%-5s %-10s %-70s %-5s %-5s" % ("ID", "List", "Name", "Fam", "Fam-found")
		for atag in tags:
			str += " %-3s" % atag
		print str + " Same/Diff"
		#print each pathway
		totsame = 0
		totdiff = 0
		famvalid = [0] * self.famTot
		for fam in self.famMapped:
			famvalid[fam] = 1

		for p in range(self.pathTot):
			add = sum(pathvalid[p])
			if add == 0:
				continue
			add = sum(pathvalid[p][1:])
			if add == mapnum - 1 or add == 0:
				label = "Same" 
				if add != 0:
					totsame += 1
			else:
				label = "Diff"
				totdiff += 1
			add = 0
			for f in self.path2Fam[p]:
				if famvalid[f] == 1:
					add += 1
			str = "%-5d %-10s %-70s %-5d %-5d" % (p + 1, self.pathList[p], self.pathName[p], len(self.path2Fam[p]), add)
			for m in range(len(maps)):
				str += " %-3d" % pathvalid[p][m] 
			print str + " " + label
		#print total number line
		str = "%-5s %-10s %-70s %-5s %-5s" % ("#total", "", "", "", "")
		for m in range(mapnum):
			str += " %-3d" % len(maps[m]) 
		print str
		#print functional diveristy line
		#str = "#functional-diversity [max: log(%d)=%.3f]" % (self.pathTot, math.log(self.pathTot * 1.0))
		#str = "%-99s" % str
		#for m in range(mapnum):
		##	str += " %.3f" % math.log(len(maps[m]) * 1.0)
		#print str
		print "#total match=%d diff=%d" % (totsame, totdiff)

	def WriteReport(self, minpath, reportfile, detailfile):
		na = True
		if self.whichDB == "KEGG":
			keggmap = [] 
			keggdir = "/dataomics/kegg/kegg-curr"
			mapfile = keggdir + "/pathway/" + self.speID.lower() + "/map.list"
			if os.path.exists(mapfile):
				file = open(mapfile, "r")
				for aline in file:
					m = re.match('"'^[^\d]+(?P<id>\d+)'"', aline)
					if m:
						id = m.group('"'id'"')
						idx = self.pathList.index(id)
						keggmap.append(idx)
				na = False
			tags = ["kegg", "naive", "minpath"]
			maps = [keggmap, self.pathMapped, minpath]
		else:
			seedmap = []	
			if self.whichDB == '"'SEED'"':
				tags = ["seed", "naive", "minpath"]
			else:
				tags = ["any", "naive", "minpath"]
			maps = [seedmap, self.pathMapped, minpath]
		mapnum = len(maps)
		pathvalid = intmatrix(self.pathTot, mapnum)
		for m in range(mapnum):
			for p in maps[m]:
				pathvalid[p][m] = 1
		famvalid = [0] * self.famTot
		for fam in self.famMapped:
			famvalid[fam] = 1

		file = open(reportfile, "w")
		if detailfile:
			detail = open(detailfile, "w")
		for p in range(self.pathTot):
			add = sum(pathvalid[p])
			if add == 0:
				continue
			add = 0
			for f in self.path2Fam[p]:
				if famvalid[f] == 1:
					add += 1
			if na:
				tmp = "n/a"
			else:
				tmp = pathvalid[p][0]
			print >> file, "path", self.pathList[p], tags[0], tmp, " naive", pathvalid[p][1], " minpath", pathvalid[p][2], " fam0 ", len(self.path2Fam[p]), " fam-found ", add, " name ", self.pathName[p]

			if not (detailfile and pathvalid[p][2]):
				continue
			#print details
			print >> detail, "path", self.pathList[p], "fam0", len(self.path2Fam[p]), "fam-found", add, "#", self.pathName[p]
			for f in self.path2Fam[p]:
				if famvalid[f] == 1 and self.famCount:
					print >> detail, "  ", self.famID[f], "hits", self.famCount[f], "#", self.famName[f]
				elif famvalid[f] == 1:
					print >> detail, "  ", self.famID[f], "#", self.famName[f]
					
		print
		print "Results are saved in file:", reportfile 
		file.close()
		if detailfile:
			print "Details are saved in file:", detailfile
			detail.close()

# this function reads in KO/fig assignment, then map KO/fig families to the pathways
# last update by Yuzhen Ye on July 3, 2009
def Orth2Path(infile = "demo.ko", whichdb = "KEGG", mpsfile = "test.mps", reportfile = "test.minpath", detailfile = "", mapfile=""):
	try:
		file = open(infile, "r")
	except IOError:
		sys.exit( "open file error " + infile)

	orthlist, orthcount = [], []
	add = 0
	for aline in file:
		#aline = aline.strip()
		#tmp = aline.split("\t")
		tmp = aline.strip().split()
		if len(tmp) >= 2:
			add = add + 1
			if tmp[1] not in orthlist:
				orthlist.append(tmp[1])
				orthcount.append(1)
				#print "ko-%d=%s" % (len(orthlist), tmp[1])
			else:
				idx = orthlist.index(tmp[1])
				orthcount[idx] += 1
	file.close()

	print "total input orth=%d  unique=%d" % (add, len(orthlist))

	test = MinPath(whichdb = whichdb, mapfile = mapfile) #default pathwaydb: KEGG

	if whichdb == "KEGG":
		map = test.Orth2PathMin(famidxlist=orthlist, famnamelist=[], famcount=orthcount, mpsfile=mpsfile)
	elif whichdb == "SEED":
		map = test.Orth2PathMin(famidxlist=[], famnamelist=orthlist, famcount=orthcount, mpsfile=mpsfile)
	else:
		map = test.Orth2PathMin(famidxlist=[], famnamelist=orthlist, famcount=orthcount, mpsfile=mpsfile)
	#KEGG by ids, and fig by names
		
	map_add = map[:]
	par = 0.5
	test.PopulatePath(pathmapped = map_add, par = par)

	test.WriteReport(map_add, reportfile, detailfile)

	os.system("rm test.mps*")

if __name__ == '"'__main__'"':
	kofile, figfile, anyfile, mapfile, mpsfile, reportfile, detailfile = "", "", "", "", "test.mps", "test.minpath", ""
	for i in range(len(sys.argv)):
		if sys.argv[i] == "-ko":
			kofile = sys.argv[i + 1]
		elif sys.argv[i] == "-fig":
			figfile = sys.argv[i + 1]
		elif sys.argv[i] == "-any":
			anyfile = sys.argv[i + 1]
		elif sys.argv[i] == "-map":
			mapfile = sys.argv[i + 1]
		elif sys.argv[i] == "-report":
			reportfile = sys.argv[i + 1]
		elif sys.argv[i] == "-details":
			detailfile = sys.argv[i + 1]
		elif sys.argv[i] == "-mps":
			mpsfile = sys.argv[i + 1]
	if kofile:
		Orth2Path(infile = kofile, mpsfile = mpsfile, reportfile = reportfile, detailfile = detailfile)
	elif figfile:
		Orth2Path(infile = figfile, mpsfile = mpsfile, reportfile = reportfile, detailfile = detailfile, whichdb = "SEED")
	elif anyfile and mapfile:
		Orth2Path(infile = anyfile, mpsfile = mpsfile, reportfile = reportfile, detailfile = detailfile, whichdb = "ANY", mapfile=mapfile)
	else:
		print "Usage: python MinPath.py <-ko filename>/<-fig filename>/<-any annfile> [-map mapfile] [-report filename] [-details detailed-output]"
		print "Note: your input file can contain functional annotations in either of the following"
		print "   -ko file: annotation in KEGG KO families"
		print "   -fig file: annotation in SEED fig families"
		print "   -any file: annotation in any families, then you must specify -map, the pathway-function mapping file"
		print "Example 1: python MinPath.py -ko demo.ko -report demo.ko.minpath"
		print "Example 2: python MinPath.py -ko demo.ko -report demo.ko.minpath -details demo.ko.minpath.details"
		print "Example 3: python MinPath.py -fig demo.fig -report demo.fig.minpath"
		print "Example 4: python MinPath.py -fig demo.fig -report demo.fig.minpath -details demo.fig.minpath.details"
		print "Example 5: python MinPath.py -any demo.ec -map ec2path -report demo.ec.minpath -details demo.ec.minpath.details"
		sys.exit(1)' > MinPath1.2.py

} 

function genes2kronaFunction {

	echo '#!/usr/bin/env python
#################################

import sys,csv, numpy as np

def ReadLengths(f):
    hin = open(f)
    d = {}
    for line in hin:
        line = line.rstrip()
        [g,l] = line.rsplit()
        l = float(l)
        d[g] = l
    return d

def ReadLimits(f):
    limit = {}
    hin = open(f)
    for line in hin:
        line = line.rstrip()
        limit[line.rsplit()[-1]] = ""
    hin.close()
    return limit

def ReadMap(f):
    d = {}
    hin = open(f)
    hincsv = csv.reader(hin, delimiter = '"'\t'"')
    for row in hincsv:
        ann = row[1]
        parent = row[0]
        try: d[ann].append(parent)
        except KeyError: d[ann] = [parent]
    hin.close()
    return d

def ReadCoverage(f, lengths):
    d = {}
    try: hin = open(f, '"'r'"')
    except TypeError: return {}
    hincsv = csv.reader(hin, delimiter = '"'\t'"')
    for row in hincsv:
        g = row[0]
        try: c = float(row[1])
        except ValueError: continue
        if len(lengths.keys())>0:
            try: l = lengths[g]
            except KeyError: sys.exit("No length found for gene "+g+"\n")
        else: l = 1
        c = float(c)/l
        d[g] = c
    hin.close()
    return d

def ReadHierarchy(f):
    d = {}
    hin = open(f)
    hincsv = csv.reader(hin, delimiter = '"'\t'"')
    l = []
    for row in hincsv:
        id = row[0]
        name = row[1]
        try: hiers = row[2]
        except IndexError: hiers = "Unknown"
        d[id] = hiers.split("|")+[name]
        l.append(len(d[id]))
    hin.close()
    return (max(l),d)

def Calculate(hier_c, operation):
    hier_sdev = {}
    if operation == "sum": function = np.sum
    elif operation == "mean" or operation == "meanhalf": function = np.mean
    for hier, l in hier_c.iteritems():
        hier_sdev[hier] = 0.0
        if operation == "meanhalf":
            l.sort()
            l = l[len(l)/2:]
        hier_c[hier] = function(l)
        hier_sdev[hier] = np.std(l)
    return (hier_sdev,hier_c)

def CalcHierarchy(mapping, annotations, coverage, operation, limit, verbose):
    ## Iterate annotations, and sum coverage if available
    ann_c = {}
    if len(coverage.keys()) > 0: cov = True
    else: cov = False
    for annotation, l in annotations.iteritems():
        covlist = []
        if cov:
            for gene in l:
                try: gene_cov = coverage[gene]
                #except KeyError: sys.exit("ERROR: Could not find coverage for gene "+str(gene)+". Are you sure you have coverage information?\n")
                except KeyError: gene_cov=0
		covlist.append(gene_cov)
            ann_c[annotation] = np.mean(covlist)
        else: ann_c[annotation] = len(l)
    ## Transfer annotation sums to nearest parent in mapping, if limit is supplied skip parents in limit
    hier_c = {}
    for annotation, count in ann_c.iteritems():

        try: parents = mapping[annotation]
        except KeyError:
            if verbose: sys.stderr.write("WARNING: Could not find hierarchy parent for "+annotation+"\n")
            continue
        for parent in parents:
            if limit and not parent in limit:
                if verbose: sys.stderr.write("Skipping parent "+ parent+"\n")
                continue
            try: hier_c[parent].append(count)
            except KeyError: hier_c[parent] = [count]
    (hier_sdev,hier_c) = Calculate(hier_c, operation)
    return (hier_sdev,hier_c)

def main():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", type=str,
        help="Tab-delimited file with gene_ids in first column and gene_annotations in second column")
    parser.add_argument("-m", "--mapfile", type=str,
        help="Tab-delimited file mapping each gene_annotation to the nearest parent hierarchy level (e.g. pathway to enzyme)")
    parser.add_argument("-H", "--hierarchy", type=str,
        help="Hierarchy file for parent levels")
    parser.add_argument("-l", "--limit", type=str,
        help="Limit calculations to only this list of parent hierarchies. E.g. a list of predicted pathways only")
    parser.add_argument("-O", "--operation", type=str, default="mean",
        help="Specify how to do calculations, either '"'mean'"' (default), '"'sum'"' or '"'meanhalf'"' where averages will be calculated from the most abundant half of annotations")
    parser.add_argument("-o", "--outfile", type=str,
        help="Write results to outfile. Defaults to stdout")
    parser.add_argument("-n", "--name", type=str,
        help="OPTIONAL: Name to assign to outfile. Defaults to name of infile")
    parser.add_argument("-c", "--coverage", type=str,
        help="OPTIONAL: Supply a file of coverage for each gene in the infile")
    parser.add_argument("-v", "--verbose", action="store_true",
        help="Run in verbose mode")
    parser.add_argument("-s", "--singlecol", action="store_true",
        help="Write only one annotation column for the first parent hierarchy (e.g. pathway)")
    parser.add_argument("-d", "--sdev", type=str,
        help="Write the standard deviation of each first parent hierarchy level to this file")
    parser.add_argument("-L", "--lengthnorm", type=str,
        help="Provide file with lengths for genes to normalize coverage by")
    parser.add_argument("-g", "--missing", type=str, default="missing_pathways.txt",
                        help="File name for missing pathways in heirarchy file")
    args = parser.parse_args()

    err_f = open(args.missing, "w")

    if not args.infile or not args.mapfile or not args.hierarchy: sys.exit(parser.print_help())

    if args.limit: limit = ReadLimits(args.limit)
    else: limit = []

    if args.lengthnorm: lengths = ReadLengths(args.lengthnorm)
    else: lengths = {}

    ## Read the mapping of hierarchies
    mapping = ReadMap(args.mapfile)
    annotations = ReadMap(args.infile) ## Read annotations the same way as above, then get length of the list for counts
    coverage = ReadCoverage(args.coverage, lengths)
    (max_hier, hierarchy) = ReadHierarchy(args.hierarchy)
    (hier_sdev, hier_counts) = CalcHierarchy(mapping, annotations, coverage, args.operation, limit, args.verbose)

    if args.outfile: hout = open(args.outfile, '"'w'"')
    else: hout = sys.stdout
    houtcsv = csv.writer(hout, delimiter = '"'\t'"')

    ## Set name for sample, if not specified use the basename of the input file
    if args.name: name = args.name
    else: name = (args.infile).split("/")[-1]

    out = [name]
    if args.singlecol:out.insert(0,"X")
    else:
        for i in range(1,max_hier+1): out.append("Level"+str(i))
    houtcsv.writerow(out)

    if args.sdev:
        sdevout = open(args.sdev, '"'w'"')
        sdevout.write("X.sdev\t"+name+"\n")

    for hier,count in hier_counts.iteritems():
        out = [count]
        try: h = hierarchy[hier]
        except KeyError: h = ["Unknown"]
        if args.singlecol: out.insert(0,hier)
        else:
            try:
                out+=hierarchy[hier]
            except KeyError:
                err_f.write("%s\n" % hier)

        try: sdevout.write(hier+"\t"+str(hier_sdev[hier])+"\n")
        except NameError: pass
        houtcsv.writerow(out)
    hout.close()

    try: sdevout.close()
    except NameError: pass


if __name__ == "__main__": main()' > genes.to.kronaTable.py

}

function ectopwyFunction {
	echo '12DICHLORETHDEG-PWY	1.1.2.7	methanol dehydrogenase
12DICHLORETHDEG-PWY	1.2.1.4	acetaldehyde dehydrogenase
12DICHLORETHDEG-PWY	3.8.1.3 	2-haloacid dehalogenase
12DICHLORETHDEG-PWY	3.8.1.5	DHLAXANAU-MONOMER
14DICHLORBENZDEG-PWY 	1.3.1.32	chloromaleylacetate reductase
14DICHLORBENZDEG-PWY	3.1.1.45	dienelactone hydrolase
2AMINOBENZDEG-PWY	6.2.1.32	anaerobic aminobenzoate-CoA ligase
2AMINOBENZDEG-PWY 	6.2.1.32 	benzoate-CoA ligase
2ASDEG-PWY	1.13.11	3SC23OALCAL-MONOMER
2OXOBUTYRATECAT-PWY	1.2.7.11 	2-oxoacid:ferredoxin oxidoreductase
3-HYDROXYPHENYLACETATE-DEGRADATION-PWY	1.13.11.15	3,4-dihydroxyphenylacetate dioxygenase subunit
3-HYDROXYPHENYLACETATE-DEGRADATION-PWY	1.13.11.15	MONOMER-2843
3-HYDROXYPHENYLACETATE-DEGRADATION-PWY	1.2.1.60	subunit of 5-carboxymethyl-2-hydroxymuconic semialdehyde dehydrogenase
3-HYDROXYPHENYLACETATE-DEGRADATION-PWY	4.1.1.68	MONOMER-604
3-HYDROXYPHENYLACETATE-DEGRADATION-PWY	4.1.2.52	MONOMER-2849
3-HYDROXYPHENYLACETATE-DEGRADATION-PWY	4.2.1.163	2-hydroxyhepta-2,4-dienedioate hydratase subunit
3-HYDROXYPHENYLACETATE-DEGRADATION-PWY	5.3.3.10	MONOMER-2845
4AMINOBUTMETAB-PWY 	1.2.1.24	NAD<sup>+</sup>-dependent succinate semialdehyde dehydrogenase
4AMINOBUTMETAB-PWY 	2.6.1.19 	4-aminobutyrate aminotransferase
4AMINOBUTMETAB-PWY 	2.6.1.19	4-aminobutyrate aminotransferase
4-HYDROXYMANDELATE-DEGRADATION-PWY 	1.1.5 	(S)-mandelate dehydrogenase
4-HYDROXYMANDELATE-DEGRADATION-PWY 	1.2.1.64 	4-hydroxybenzaldehyde dehydrogenase
4-HYDROXYMANDELATE-DEGRADATION-PWY 	1.2.1.96 	NADP<sup>+</sup>-benzaldehyde dehydrogenase
4-HYDROXYMANDELATE-DEGRADATION-PWY 	4.1.1.7	benzoylformate decarboxylase
4-HYDROXYMANDELATE-DEGRADATION-PWY 	5.1.2.2	mandelate racemase
4TOLCARBDEG-PWY 	1.1.1.M5 	4-sulfobenzyl alcohol dehydrogenase subunit
4TOLCARBDEG-PWY 	1.2.1 	TsaD
6-HYDROXYCINEOLE-DEGRADATION-PWY	1.1.1.241	MONOMER-5721
6-HYDROXYCINEOLE-DEGRADATION-PWY	1.14.13.51	MONOMER-5722
ACETATEUTIL-PWY 	2.3.1.8	phosphate acetyltransferase
ACETATEUTIL-PWY 	2.3.1.8 	phosphate acetyltransferase / phosphate propionyltransferase
ACETATEUTIL-PWY 	2.7.2.1 	acetate kinase
ACETOACETATE-DEG-PWY	2.8.3.9	&beta; complex 
AEROBACTINSYN-PWY	1.14.13.59	L-lysine 6-monooxygenase subunit
AEROBACTINSYN-PWY	2.3.1.102	N<sup>6</sup>-hydroxylysine O-acetyltransferase subunit
AEROBACTINSYN-PWY	6.3.2.38	citrate:N<SUP>6</SUP>-acetyl-N<SUP>6</SUP>-hydroxy-L-lysine ligase
AEROBACTINSYN-PWY	6.3.2.39	aerobactin synthase
ALADEG-PWY	1.4.5.1 	D-amino acid dehydrogenase
ALANINE-DEG3-PWY	2.6.1.2	alanine aminotransferase subunit
ALANINE-SYN2-PWY 	2.6.1.2	alanine aminotransferase 1
ALANINE-SYN2-PWY 	2.6.1.2	alanine aminotransferase 2
ALKANEMONOX-PWY	1.14.14.5 	alkanesulfonate monooxygenase, FMNH2-dependent
ALKANEMONOX-PWY	1.5.1.38	NADPH-dependent FMN reductase
ALLANTOINDEG-PWY 	3.5.2.5	allantoinase
ALLANTOINDEG-PWY 	3.5.3.4	allantoicase trimer
ALLANTOINDEG-PWY 	4.3.2.3	YIR032C-MONOMER
ALLANTOINDEG-PWY 	6.3.4.6 	urea amidolyase
ALL-CHORISMATE-PWY 	1.14.13	2-octaprenylphenol hydroxylase
ALL-CHORISMATE-PWY 	1.14.13	OCTAPRENYL-METHOXYPHENOL-OH-MONOMER
ALL-CHORISMATE-PWY 	1.14.99	OCTAPRENYL-METHYL-METHOXY-BENZOQ-OH-MON
ALL-CHORISMATE-PWY 	1.3.1.28	2,3-dihydro-2,3-dihydroxybenzoate dehydrogenase
ALL-CHORISMATE-PWY 	2.1.1.222 	bifunctional 3-demethylubiquinone-8 3-O-methyltransferase and 2-octaprenyl-6-hydroxyphenol methylase
ALL-CHORISMATE-PWY 	2.2.1.9	2-succinyl-5-enolpyruvyl-6-hydroxy-3-cyclohexene-1-carboxylate synthase
ALL-CHORISMATE-PWY 	2.4.2.18 	anthranilate synthase component I 
ALL-CHORISMATE-PWY 	2.4.2.18 	anthranilate synthase component II
ALL-CHORISMATE-PWY 	2.5.1.39	4OHBENZOATE-OCTAPRENYLTRANSFER-MONOMER
ALL-CHORISMATE-PWY 	2.5.1.54	2-dehydro-3-deoxyphosphoheptonate aldolase
ALL-CHORISMATE-PWY 	2.5.1.74	DMK-MONOMER
ALL-CHORISMATE-PWY 	2.5.1.90	octaprenyl diphosphate synthase
ALL-CHORISMATE-PWY 	2.6.1.57 	branched-chain amino-acid aminotransferase
ALL-CHORISMATE-PWY 	2.7.7 	holo-EntB dimer 
ALL-CHORISMATE-PWY 	2.7.8.7 	ENTD-MONOMER
ALL-CHORISMATE-PWY 	2.7.8.7	holo-[acyl-carrier-protein] synthase
ALL-CHORISMATE-PWY 	2.7.8.7	holo-[acyl carrier protein] synthase 2
ALL-CHORISMATE-PWY 	3.1.2.28	1,4-dihydroxy-2-naphthoyl-CoA thioesterase
ALL-CHORISMATE-PWY 	3.3.2.1 	holo-EntB dimer
ALL-CHORISMATE-PWY 	4.1.1.48 	PRAI-IGPS
ALL-CHORISMATE-PWY 	4.1.1.98	3-octaprenyl-4-hydroxybenzoate carboxy-lyase
ALL-CHORISMATE-PWY 	4.1.3.36	1,4-dihydroxy-2-naphthoyl-CoA synthase
ALL-CHORISMATE-PWY 	4.2.1.113	O-SUCCINYLBENZOATE-COA-SYN-MONOMER
ALL-CHORISMATE-PWY 	4.2.1.20 	tryptophan synthase, &alpha; subunit
ALL-CHORISMATE-PWY 	4.2.1.20 	tryptophan synthase, &alpha; subunit 
ALL-CHORISMATE-PWY 	4.2.1.20 	tryptophan synthase, &beta; subunit dimer
ALL-CHORISMATE-PWY 	4.2.3.4	AROB-MONOMER
ALL-CHORISMATE-PWY 	4.2.99.20	EG12438-MONOMER
ALL-CHORISMATE-PWY 	5.4.4.2	isochorismate synthase 1
ALL-CHORISMATE-PWY 	5.4.4.2	isochorismate synthase 2
ALL-CHORISMATE-PWY 	6.2.1.26	o-succinylbenzoate-CoA ligase
ALL-CHORISMATE-PWY 	6.3.2.14 	2,3-dihydroxybenzoate-AMP ligase
ALL-CHORISMATE-PWY 	6.3.2.14 	holo [EntF peptidyl-carrier protein]
AMMASSIM-PWY 	6.3.1.2	glutamine synthetase
AMMOXID-PWY	1.14.99.39	ammonia monooxygenase &alpha; subunit 
AMMOXID-PWY	1.7.2.6	cytochrome P460 monomer
AMMOXID-PWY	1.7.2.6	Cytochrome P460 monomer
ANAEROFRUCAT-PWY 	1.1.1.27	lactate dehydrogenase subunit
ANAEROFRUCAT-PWY 	1.1.1.27	L-lactate dehydrogenase subunit
ANAEROFRUCAT-PWY 	1.2.1.12	glyceraldehyde-3-phosphate dehydrogenase
ANAEROFRUCAT-PWY	1.2.1.12	glyceraldehyde-3-phosphate dehydrogenase subunit
ANAEROFRUCAT-PWY 	1.2.1.12	testis-specific glyceraldehyde-3-phosphate dehydrogenase
ANAEROFRUCAT-PWY 	2.7.1.11	6-phosphofructokinase, liver type
ANAEROFRUCAT-PWY 	2.7.1.11	6-phosphofructokinase, muscle type
ANAEROFRUCAT-PWY	2.7.1.11	6-phosphofructokinase subunit
ANAEROFRUCAT-PWY 	2.7.1.11	6-phosphofructokinase type C
ANAEROFRUCAT-PWY 	2.7.1.1 	glucokinase
ANAEROFRUCAT-PWY 	2.7.1.1 	hexokinase-1
ANAEROFRUCAT-PWY 	2.7.1.1 	hexokinase-2
ANAEROFRUCAT-PWY 	2.7.1.1 	hexokinase-3
ANAEROFRUCAT-PWY 	2.7.1.40	pyruvate kinase M1/M2
ANAEROFRUCAT-PWY 	2.7.1.40	pyruvate kinase R/L
ANAEROFRUCAT-PWY	2.7.1.40	pyruvate kinase subunit
ANAEROFRUCAT-PWY 	2.7.2.3	phosphoglycerate kinase 1
ANAEROFRUCAT-PWY 	2.7.2.3	phosphoglycerate kinase 2
ANAEROFRUCAT-PWY	2.7.2.3	phosphoglycerate kinase subunit
ANAEROFRUCAT-PWY 	4.1.2.13	fructose-bisphosphate aldolase A
ANAEROFRUCAT-PWY 	4.1.2.13	fructose-bisphosphate aldolase B
ANAEROFRUCAT-PWY 	4.1.2.13	fructose-bisphosphate aldolase C
ANAEROFRUCAT-PWY	4.1.2.13	fructose bisphosphate aldolase subunit
ANAEROFRUCAT-PWY 	4.2.1.11	alpha enolase
ANAEROFRUCAT-PWY 	4.2.1.11	beta enolase
ANAEROFRUCAT-PWY	4.2.1.11	enolase subunit
ANAEROFRUCAT-PWY 	4.2.1.11	gamma enolase
ANAEROFRUCAT-PWY 	5.3.1.1	triosephosphate isomerase
ANAEROFRUCAT-PWY	5.3.1.1	triosephosphate isomerase subunit
ANAEROFRUCAT-PWY 	5.3.1.9	glucose-6-phosphate isomerase
ANAEROFRUCAT-PWY	5.3.1.9	MONOMER-13046
ANAEROFRUCAT-PWY 	5.4.2.11 	bisphosphoglycerate mutase
ANAGLYCOLYSIS-PWY	1.2.1.12	Gap
ANAGLYCOLYSIS-PWY	2.7.1.11	phosphofructokinase
ANAGLYCOLYSIS-PWY	2.7.1.1 	glucokinase subunit
ANAGLYCOLYSIS-PWY	2.7.1.40	pyruvate kinase
ANAGLYCOLYSIS-PWY	2.7.2.3 	Pgk
ANAGLYCOLYSIS-PWY	4.1.2.13	MONOMER-364
ANAGLYCOLYSIS-PWY	4.2.1.11	Eno dimer
ANAGLYCOLYSIS-PWY	5.3.1.9	MONOMER-381
ANAGLYCOLYSIS-PWY	5.4.2.11	2,3-bisphosphoglycerate-dependent phosphoglycerate mutase
ARGASEDEG-PWY 	1.2.1.88	AT5G62530-MONOMER
ARGASEDEG-PWY 	1.2.1.88 	&Delta;<sup>1</sup>-pyrroline-5-carboxylate dehydrogenase, mitochondrial
ARGASEDEG-PWY	1.2.1.88	ROCABACSU-MONOMER
ARGASEDEG-PWY	2.6.1.13	ornithine &delta;-aminotransferase
ARGASEDEG-PWY	2.6.1.13	ROCDBACSU-MONOMER
ARGASEDEG-PWY	3.5.3.1	ROCFBACSU-MONOMER
ARGDEG-IV-PWY	1.2.1.54	MONOMER-11560
ARGDEG-IV-PWY	1.4.3.2	MONOMER-11561
ARGDEG-IV-PWY	3.5.3.7	guanidinobutyrase subunit
ARGDEG-IV-PWY	4.1.1.75	MONOMER-11559
ARGDEGRAD-PWY 	2.1.3.3	ArcB
ARGDEGRAD-PWY 	2.1.3.3	AT1G75330-MONOMER
ARGDEGRAD-PWY 	2.1.3.3	MONOMER-17454
ARGDEGRAD-PWY 	2.1.3.3	MONOMER-533
ARGDEGRAD-PWY 	2.7.2.2	carbamate kinase subunit
ARGDEGRAD-PWY 	2.7.2.2	MONOMER-17457
ARGDEGRAD-PWY 	2.7.2.2	MONOMER-534
ARGDEGRAD-PWY 	2.7.2.2	MONOMER-663
ARGDEGRAD-PWY	3.5.3.6	ArcA
ARGDEGRAD-PWY	3.5.3.6	arginine deiminase subunit
ARGDEGRAD-PWY	3.5.3.6	MONOMER-507
ARGDEG-V-PWY	1.13.12.1	MONOMER-11564
ARGDEG-V-PWY	1.13.12.1	MONOMER-11570
ARGDEG-V-PWY	3.5.1.4	MONOMER-11565
ARGDEG-V-PWY	3.5.3.7	guanidinobutyrase subunit
ARGDEG-V-PWY	3.5.3.7	MONOMER-11566
ARG-GLU-PWY 	3.5.3.1	arginase subunit
ARGORNPROST-PWY 	1.21.4.1 	D-proline reductase complex
ARGORNPROST-PWY 	1.5.1.2	ProC
ARGORNPROST-PWY 	2.1.3.3	ORNCARBAMTRANSFERST-MONOMER
ARGORNPROST-PWY 	2.6.1.13	ORNAMTRANST-MONOMER
ARGORNPROST-PWY 	3.5.3.6	ARGDEIMST-MONOMER
ARGORNPROST-PWY 	4.3.1.12	ornithine cyclodeaminase subunit
ARGORNPROST-PWY 	5.1.1.12	ornithine racemase monomer
ARGORNPROST-PWY 	5.4.3.5	D-ornithine aminomutase component E subunit 
ARG+POLYAMINE-SYN 	1.2.1.38	N-ACETYLGLUTPREDUCT-MONOMER
ARG+POLYAMINE-SYN 	2.1.3.3	ornithine carbamoyltransferase chain F
ARG+POLYAMINE-SYN 	2.1.3.3	ornithine carbamoyltransferase chain I
ARG+POLYAMINE-SYN 	2.3.1.1	N-acetylglutamate synthase
ARG+POLYAMINE-SYN 	2.6.1.11 	N-acetylornithine aminotransferase / N-succinyldiaminopimelate aminotransferase
ARG+POLYAMINE-SYN 	2.7.2.8	acetylglutamate kinase
ARG+POLYAMINE-SYN 	3.5.1.16	acetylornithine deacetylase
ARG+POLYAMINE-SYN 	4.1.1.17	ornithine decarboxylase, biosynthetic
ARG+POLYAMINE-SYN 	4.1.1.18	lysine decarboxylase 1
ARG+POLYAMINE-SYN 	4.1.1.18	lysine decarboxylase 2
ARG+POLYAMINE-SYN 	4.1.1.19	arginine decarboxylase, biosynthetic
ARG+POLYAMINE-SYN 	4.1.1.50	&alpha; cleavage product of SpeD 
ARG+POLYAMINE-SYN 	4.3.2.1	ARGSUCCINLYA-MONOMER
ARG-PRO-PWY	1.5.1.2	pyrroline-5-carboxylate reductase subunit
ARG-PRO-PWY	2.6.1.13	YLR438W-MONOMER
ARG-PRO-PWY	3.5.3.1	arginase subunit
ARGSPECAT-PWY	2.5.1.22	spermine synthase
ARGSPECAT-PWY 	4.1.1.50	S-adenosylmethionine decarboxylase
ARGSYNBSUB-PWY	1.2.1.38	ARGCBACSU-MONOMER
ARGSYNBSUB-PWY	2.1.3.3	ArgF
ARGSYNBSUB-PWY	2.3.1.1	acetylglutamate synthase
ARGSYNBSUB-PWY	2.3.1.1	L-glutamate N-acetyltransferase
ARGSYNBSUB-PWY	2.3.1.35	acetylornithine acetyltransferase
ARGSYNBSUB-PWY	2.3.1.35 	ArgJ
ARGSYNBSUB-PWY	2.6.1.11	acetylornithine aminotransferase
ARGSYNBSUB-PWY	2.6.1.11	ARGDBACSU-MONOMER
ARGSYNBSUB-PWY	2.7.2.8	ArgB
ARGSYNBSUB-PWY	2.7.2.8 	YER069W-MONOMER
ARGSYNBSUB-PWY	4.3.2.1	ARGHBACSU-MONOMER
ARGSYNBSUB-PWY	4.3.2.1	argininosuccinate lyase
ARGSYNBSUB-PWY	6.3.4.5	ArgG
ARGSYNBSUB-PWY	6.3.4.5	arginosuccinate synthetase
ARGSYNBSUB-PWY	6.3.5.5 	carbamoyl-phosphate synthetase A
ARO-PWY 	1.1.1.25	ARODBACSU-MONOMER
ARO-PWY 	1.1.1.25	AROE-MONOMER
ARO-PWY 	1.1.1.25 	shikimate dehydrogenase
ARO-PWY 	1.1.1.282 	shikimate dehydrogenase / quinate dehydrogenase
ARO-PWY 	2.5.1.19	AROA-MONOMER
ARO-PWY 	2.5.1.19	AROEBACSU-MONOMER
ARO-PWY 	2.7.1.71	AROIBACSU-MONOMER
ARO-PWY	2.7.1.71	shikimate kinase
ARO-PWY 	2.7.1.71	shikimate kinase II
ARO-PWY 	4.2.1.10	3-dehydroquinate dehydratase
ARO-PWY 	4.2.1.10	AroC
ARO-PWY 	4.2.3.4	AROBBACSU-MONOMER
ARO-PWY 	4.2.3.5	AROFBACSU-MONOMER
ARO-PWY 	4.2.3.5	chorismate synthase
ASPARAGINE-BIOSYNTHESIS 	6.3.5.4 	asparagine synthetase B
ASPARAGINE-DEG1-PWY-1	3.1.1.47 	60 kDa lysophospholipase
ASPARAGINE-DEG1-PWY-1	3.4.19.5 	isoaspartyl peptidase/L-asparaginase
ASPARAGINE-DEG1-PWY-1	3.5.1.26	N(4)-(beta-N-acetylglucosaminyl)-L-asparaginase
ASPARAGINE-DEG1-PWY 	3.4.19 	isoaspartyl dipeptidase / asparaginase III
ASPARAGINE-DEG1-PWY 	3.5.1.38 	asparaginase I
ASPARAGINE-DEG1-PWY 	3.5.1.38 	asparaginase II
ASPARAGINE-DEG1-PWY	3.5.1.38 	L-asparaginase subunit
ASPARAGINE-DEG1-PWY	3.5.1.38 	MONOMER-9561
ASPARTATE-DEG1-PWY	2.6.1.1	aspartate aminotransferase subunit
ASPARTATE-DEG1-PWY 	2.6.1.1 	cysteine aminotransferase subunit
ASPARTATE-DEG1-PWY	2.6.1.1	MONOMER-13012
ASPASN-PWY 	2.6.1.1 	aspartate aminotransferase, PLP-dependent
ASPASN-PWY 	4.3.1.1	aspartate ammonia-lyase
ASPSYNII-PWY	3.5.5.4 	&beta;-cyano-L-alanine hydrolase
ASPSYNII-PWY	3.5.5.4 	MONOMER-17621
ASPSYNII-PWY	3.5.5.4 	MONOMER-17622
ASPSYNII-PWY	4.4.1.9	AT3G61440-MONOMER
ASPSYNII-PWY	4.4.1.9	&beta;-cyano-L-alanine synthase
ASPSYNII-PWY	4.4.1.9	MONOMER-16250
ASPSYNII-PWY	4.4.1.9	MONOMER-16251
AST-PWY	1.2.1.71	aldehyde dehydrogenase
AST-PWY	2.3.1.109	arginine succinyltransferase &alpha; subunit 
AST-PWY	2.3.1.109	ARGSUCCTRAN-MONOMER
AST-PWY	2.6.1.13 	succinylornithine transaminase subunit
AST-PWY	2.6.1.81 	succinylornithine transaminase
AST-PWY	3.5.1.96	SUCCGLUDESUCC-MONOMER
AST-PWY	3.5.3.23	succinylarginine dihydrolase
BENZCOA-PWY 	1.1.1.368	6-hydroxycyclohex-1-ene-1-carbonyl-CoA dehydrogenase
BENZCOA-PWY 	1.17.5.1	phenylacetyl-CoA dehydrogenase
BENZCOA-PWY 	1.2.1.39	phenylacetaldehyde dehydrogenase
BENZCOA-PWY 	1.2	phenylglyoxylate:acceptor oxidoreductase
BENZCOA-PWY 	1.3.7.8	benzoyl-CoA reductase &alpha; subunit 
BENZCOA-PWY 	1.3.7.9	4-hydroxybenzoyl-CoA reductase &alpha; subunit 
BENZCOA-PWY 	2.6.1.57	phenylalanine transaminase
BENZCOA-PWY 	3.7.1.21	&beta;-oxoacyl-CoA hydrolase subunit
BENZCOA-PWY 	4.1.1.43	phenylpyruvate decarboxylase
BENZCOA-PWY 	4.2.1.100	cyclohexa-1,5-diene-1-carboxyl-CoA hydratase subunit
BENZCOA-PWY 	6.2.1.27	4-hydroxybenzoate-CoA ligase
BENZCOA-PWY 	6.2.1.30	phenylacetate-CoA ligase (anaerobic)
BETA-ALA-DEGRADATION-I-PWY	1.2.1	malonate-semialdehyde dehydrogenase
BETA-ALA-DEGRADATION-I-PWY	2.6.1.19	&beta;-alanine aminotransferase subunit
BETSYN-PWY	1.1.99.1	choline dehydrogenase
BETSYN-PWY	1.1.99.1	MONOMER-8641
BETSYN-PWY	1.2.1.8	betaine aldehyde dehydrogenase
BGALACT-PWY	3.2.1.23 	&beta;-galactosidase
BGALACT-PWY	3.2.1 	&beta;-galactosidase
BIOTIN-BIOSYNTHESIS-PWY 	2.1.1.197	malonyl-[acp] methyltransferase
BIOTIN-BIOSYNTHESIS-PWY 	2.3.1.41 	&beta;-ketoacyl-ACP synthase I
BIOTIN-BIOSYNTHESIS-PWY 	2.3.1.47	8-amino-7-oxononanoate synthase
BIOTIN-BIOSYNTHESIS-PWY 	2.6.1.62	adenosylmethionine-8-amino-7-oxononanoate aminotransferase
BIOTIN-BIOSYNTHESIS-PWY 	2.8.1.6	biotin synthase
BIOTIN-BIOSYNTHESIS-PWY 	3.1.1.85 	pimeloyl-[acp] methyl ester esterase
BIOTIN-BIOSYNTHESIS-PWY 	3.1.1.85	pimeloyl-[acyl-carrier protein] methyl ester esterase
BIOTIN-BIOSYNTHESIS-PWY 	4.2.1.59 	3-hydroxy-acyl-[acyl-carrier-protein] dehydratase
BIOTIN-BIOSYNTHESIS-PWY 	6.3.3.3	dethiobiotin synthetase
BRANCHED-CHAIN-AA-SYN-PWY 	2.2.1.6	acetohydroxy acid synthase I, large subunit 
BRANCHED-CHAIN-AA-SYN-PWY 	2.3.3.13	2-isopropylmalate synthase
BRANCHED-CHAIN-AA-SYN-PWY 	2.6.1.57 	tyrosine aminotransferase
BSUBPOLYAMSYN-PWY 	2.5.1.16 	spermidine synthase
BSUBPOLYAMSYN-PWY	2.5.1.16	spermidine synthase
CALVIN-PWY 	1.2.1.13 	chloroplastic glyceraldehyde 3-phosphate dehydrogenase
CALVIN-PWY	3.1.3.11	fructose 1,6-bisphosphatase
CAMALEXIN-SYN 	1.14.13.125 	CYP79B3
CAMALEXIN-SYN 	1.14.13.125 	L-tryptophan N-monooxygenase
CAMALEXIN-SYN 	4.99.1.6	AT2G30770-MONOMER
CARNMET-PWY	1.3.99	CROBETREDUCT-MONOMER
CARNMET-PWY	2.8.3.21	&gamma;-butyrobetainyl-CoA:carnitine CoA transferase
CARNMET-PWY	4.2.1.149	CARNRACE-MONOMER
CARNMET-PWY	6.2.1	CAIC-MONOMER
CAROTENOID-PWY 	1.13.99 	MONOMER-12386
CAROTENOID-PWY 	1.14.99.45	carotene &epsilon;-monooxygenase
CAROTENOID-PWY 	1.14.99	AT1G31800-MONOMER
CAROTENOID-PWY 	1.14.99	&beta;-ring hydroxylase subunit
CAROTENOID-PWY 	1.3.5.5	15-cis-phytoene desaturase
CAROTENOID-PWY 	1.3.5.6	&zeta;-carotene desaturase
CAROTENOID-PWY 	5.2.1.12	AT1G10830-MONOMER
CAROTENOID-PWY 	5.2.1.13 	carotenoid isomerase
CAROTENOID-PWY 	5.2.1 	carotenoid isomerase
CAROTENOID-PWY 	5.5.1.18	lycopene &epsilon;-cyclase
CAROTENOID-PWY 	5.5.1 	lycopene &beta; cyclase
CAROTENOID-PWY 	5.5.1 	lycopene &epsilon; cyclase
CATECHOL-ORTHO-CLEAVAGE-PWY	5.3.3.4 	muconolactone isomerase
CHLOROPHYLL-SYN	1.3.1.33	protochlorophyllide oxidoreductase A
CHLOROPHYLL-SYN	1.3.1.33	protochlorophyllide oxidoreductase B
CHLOROPHYLL-SYN	1.3.1.33	protochlorophyllide oxidoreductase C
CHLOROPHYLL-SYN	2.1.1.11	AT4G25080-MONOMER
CHLOROPHYLL-SYN	2.1.1.11	Mg-protoporphyrin IX methyltransferase
CHLOROPHYLL-SYN	6.6.1.1	magnesium chelatase subunit D 
CHLOROPHYLL-SYN	6.6.1.1	magnesium chelatase subunit H 
CHOLINE-BETAINE-ANA-PWY	1.1.1 	peroxisomal choline monooxygenase
CHOLINE-BETAINE-ANA-PWY	1.1.99.1	choline dehydrogenase
CHOLINE-BETAINE-ANA-PWY 	1.2.1.8 	&alpha;-aminoadipic semialdehyde dehydrogenase
CHOLINE-BETAINE-ANA-PWY	1.2.1.8	MONOMER-16787
COA-PWY-1	2.7.1.33	pantothenate kinase 1
COA-PWY-1	2.7.1.33	pantothenate kinase 2
COA-PWY-1	2.7.1.33	pantothenate kinase 3
COA-PWY-1	2.7.7.3 	bifunctional coenzyme A synthase
COA-PWY-1	4.1.1.36	phosphopantothenoylcysteine decarboxylase
COA-PWY-1	6.3.2.5 	phosphopantothenate--cysteine ligase (ATP)
COA-PWY	2.7.1.24	dephospho-CoA kinase
COA-PWY	2.7.7.3	pantetheine-phosphate adenylyltransferase
COA-PWY	4.1.1.36	phosphopantothenoylcysteine decarboxylase
COA-PWY	6.3.2.5	phosphopantothenate-cysteine ligase
COBALSYN-PWY	2.7.8.26	COBS-MONOMER
COBALSYN-PWY	3.1.3.73	predicted adenosylcobalamin phosphatase/&alpha;-ribazole phosphatase
CODH-PWY 	1.2.7.4	carbon monoxide dehydrogenase/acetyl-CoA synthase &alpha; subunit dimer
CODH-PWY	2.1.1.258	methyltetrahydrofolate:corrinoid/iron-sulfur protein methyltransferase
CODH-PWY	2.3.1.169	carbon monoxide dehydrogenase/acetyl-CoA synthase &beta; subunit
CODH-PWY 	3.5.4.9 	methenylTHF cyclohydrolase/methyleneTHF dehydrogenase subunit
CODH-PWY 	6.3.4.3	formyltetrahydrofolate synthetase subunit
COLANSYN-PWY 	1.1.1.22	UDP-glucose 6-dehydrogenase
CRNFORCAT-PWY	1.5.3.1	MONOMER-11002
CRNFORCAT-PWY	1.5.3.1	sarcosine oxidase large subunit 
CRNFORCAT-PWY	1.5.8.3	sarcosine dehydrogenase &alpha; subunit 
CRNFORCAT-PWY	3.5.2.10	creatininase subunit
CRNFORCAT-PWY	3.5.3.3	creatinase subunit
CRNFORCAT-PWY	3.5.3.3	MONOMER-11001
CYANCAT-PWY	4.2.1.104	cyanase
CYANCAT-PWY	4.2.1.1	carbonic anhydrase 1
CYCLOHEXANOL-OXIDATION-PWY	1.1.1.245	MONOMER-3121
CYCLOHEXANOL-OXIDATION-PWY	1.1.1.258	MONOMER-3124
CYCLOHEXANOL-OXIDATION-PWY	1.14.13.22	MONOMER-3122
CYCLOHEXANOL-OXIDATION-PWY 	1.2.1.63	MONOMER-3125
CYCLOHEXANOL-OXIDATION-PWY	3.1.1	MONOMER-3123
CYSTSYN-PWY	2.5.1.47	O-acetylserine (thiol) lyase
DAPLYSINESYN-PWY 	1.17.1.8	4-hydroxy-tetrahydrodipicolinate reductase
DAPLYSINESYN-PWY	1.17.1.8	4-hydroxy-tetrahydrodipicolinate reductase
DAPLYSINESYN-PWY	2.6.1.17	MONOMER-6501
DAPLYSINESYN-PWY	3.5.1.18	N-succinyl-L,L-diaminopimelate desuccinylase subunit
DAPLYSINESYN-PWY 	4.3.3.7	4-hydroxy-tetrahydrodipicolinate synthase
DAPLYSINESYN-PWY	5.1.1.7	MONOMER-6381
DENITRIFICATION-PWY	1.7.2.1	NirS
DENITRIFICATION-PWY	1.7.2.4	nitrous oxide reductase
DENITRIFICATION-PWY	1.7.2.5	nitric oxide reductase
DENITRIFICATION-PWY	1.7.5.1	respiratory nitrate reductase
DESULFONATION-PWY 	1.14.12 	2-aminobenzenesulfonate dioxygenase system / orthanilate dioxygenase system
DESULFONATION-PWY 	1.14.12	benzenesulfonate dioxygenase system, oxygenase component 
DETOX1-PWY-1 	1.11.1.6 	catalase
DETOX1-PWY	1.11.1.6 	CAT1 subunit
DETOX1-PWY	1.11.1.6 	CAT2 subunit
DETOX1-PWY	1.11.1.6 	CAT3 subunit
DETOX1-PWY	1.11.1.6 	catalase II
DETOX1-PWY	1.11.1.6 	hydroperoxidase I
DETOX1-PWY-1	1.15.1.1	extracellular superoxide dismutase [Cu-Zn]
DETOX1-PWY-1	1.15.1.1	superoxide dismutase [Cu-Zn]
DETOX1-PWY-1	1.15.1.1	superoxide dismutase [Mn]
DETOX1-PWY	1.15.1.1	superoxide dismutase (Cu-Zn)
DETOX1-PWY	1.15.1.1	superoxide dismutase (Fe)
DETOX1-PWY	1.15.1.1	superoxide dismutase (Mn)
DETOX1-PWY	3.5.4.2 	cryptic adenine deaminase
DHGLUCONATE-PYR-CAT-PWY	1.1.1.43	MONOMER-12749
DHGLUCONATE-PYR-CAT-PWY	1.1.99.3	gluconate 2-dehydrogenase complex
DHGLUCONATE-PYR-CAT-PWY	2.7.1.13	MONOMER-12748
DISSULFRED-PWY	1.8.5.M1	[DsrC]-trisulfide reductase
DISSULFRED-PWY 	1.8.99.2	adenylylsulfate reductase, &alpha; complex 
DISSULFRED-PWY	1.8.99.2	adenylylsulfate reductase &alpha; subunit 
DISSULFRED-PWY 	1.8.99.5	dissimilatory sulfite reductase
DISSULFRED-PWY	1.8.99.5	dissimilatory sulfite reductase &alpha; subunit 
DISSULFRED-PWY 	1.8.99.5	MONOMER-12408
DISSULFRED-PWY	1.8.99.5	sulfite reductase, dissimilatory &alpha; subunit 
DISSULFRED-PWY	2.7.7.4	ATP sulfurylase asubunit
DTDPRHAMSYN-PWY	1.1.1.133	MONOMER-18754
DTDPRHAMSYN-PWY	2.7.7.24	MONOMER-18751
DTDPRHAMSYN-PWY	4.2.1.46	MONOMER-18752
DTDPRHAMSYN-PWY	5.1.3.13	MONOMER-18753
ECASYN-PWY 	1.1.1.336	UDP-N-acetyl-D-mannosamine dehydrogenase
ECASYN-PWY	2.4.1.180	UDPMANACATRANS-MONOMER
ECASYN-PWY	2.4.1.325	G7800-MONOMER
ECASYN-PWY	2.7.8.33	undecaprenyl-phosphate &alpha;-N-acetylglucosaminyl transferase
ECASYN-PWY 	5.1.3.14	UDP-N-acetylglucosamine 2-epimerase
ERGOSTEROL-SYN-PWY 	1.1.1.170	C-3 sterol dehydrogenase
ERGOSTEROL-SYN-PWY 	1.1.1.270	3-keto sterol reductase
ERGOSTEROL-SYN-PWY 	1.1.1.34	HMG-CoA reductase 1
ERGOSTEROL-SYN-PWY 	1.1.1.34	HMG-CoA reductase 2
ERGOSTEROL-SYN-PWY 	1.14.13.70	cytochrome P450 51
ERGOSTEROL-SYN-PWY 	1.14.14.17	squalene monooxygenase
ERGOSTEROL-SYN-PWY 	1.14.19.20	YLR056W-MONOMER
ERGOSTEROL-SYN-PWY 	1.14.19.41	sterol C-22 desaturase
ERGOSTEROL-SYN-PWY 	1.3.1.70	YNL280C-MONOMER
ERGOSTEROL-SYN-PWY 	1.3.1.71	YGL012W-MONOMER
ERGOSTEROL-SYN-PWY 	2.1.1.41	sterol 24-C-methyltransferase subunit
ERGOSTEROL-SYN-PWY 	2.3.3.10	HmgS
ERGOSTEROL-SYN-PWY 	2.5.1.10 	farnesyl diphosphate synthase
ERGOSTEROL-SYN-PWY 	2.5.1.21 	squalene synthetase
ERGOSTEROL-SYN-PWY 	2.7.1.36	Erg12
ERGOSTEROL-SYN-PWY 	2.7.4.2	YMR220W-MONOMER
ERGOSTEROL-SYN-PWY 	4.1.1.33	Erg19
ERGOSTEROL-SYN-PWY 	5.3.3.2	YPL117C-MONOMER
ERGOSTEROL-SYN-PWY 	5.4.99.7	2,3-oxidosqualene-lanosterol cyclase
ETHYL-PWY	1.14.17.4	MONOMER-15536
ETHYL-PWY	1.14.17.4	MONOMER-15544
ETHYL-PWY	1.14.17.4	MONOMER-16249
ETHYL-PWY 	4.4.1.14	1-aminocyclopropane-1-carboxylate synthase
ETHYL-PWY	4.4.1.14	MONOMER-15063
ETHYL-PWY	4.4.1.14	MONOMER-15066
ETOH-ACETYLCOA-ANA-PWY 	1.1.1.1 	alcohol/aldehyde dehydrogenase
ETOH-ACETYLCOA-ANA-PWY 	1.1.1.1	MONOMER-13298
ETOH-ACETYLCOA-ANA-PWY 	1.2.1.10	MONOMER-13299
FAO-PWY	1.3.8	ACYLCOADEHYDROG-MONOMER
FAO-PWY	2.3.1.16 	mitochondrial 3-ketoacyl-CoA thiolase
FAO-PWY 	2.3.1.16 	Non-specific lipid-transfer protein
FAO-PWY	2.3.1.16 	trifunctional enzyme beta subunit, mitochondrial precursor
FAO-PWY 	2.3.1.223 	peroxisomal 3-ketoacyl-CoA thiolase
FASYN-ELONG-PWY	1.3.1.39 	enoyl-[acyl-carrier-protein] reductase (NADPH)
FASYN-ELONG-PWY 	2.3.1.86 	3-oxoacyl-[acyl-carrier-protein] synthase I, chloroplastic
FASYN-ELONG-PWY 	2.3.1.86 	&beta;-hydroxyacyl-ACP dehydratase/isomerase
FERMENTATION-PWY	1.1.1.28	D-lactate dehydrogenase - fermentative
FERMENTATION-PWY	1.1.7 	formate dehydrogenase H 
FERMENTATION-PWY	4.1.1.31	phosphoenolpyruvate carboxylase
FLUORENE-DEG-9-ONE-PWY	1.1.1.256	MONOMER-3092
FOLSYN-PWY 	1.5.1.3	dihydrofolate reductase
FOLSYN-PWY 	2.5.1.15	dihydropteroate synthase
FOLSYN-PWY 	2.6.1.85	aminodeoxychorismate synthase component 1 
FOLSYN-PWY 	2.7.6.3	H2PTERIDINEPYROPHOSPHOKIN-MONOMER
FOLSYN-PWY 	3.5.4.16	GTP cyclohydrolase I
FOLSYN-PWY 	3.5.4.9 	bifunctional 5,10-methylene-tetrahydrofolate dehydrogenase/ 5,10-methylene-tetrahydrofolate cyclohydrolase
FOLSYN-PWY 	3.6.1.67 	dihydroneopterin triphosphate pyrophosphohydrolase
FOLSYN-PWY 	4.1.2.25 	dihydroneopterin aldolase
FOLSYN-PWY 	4.1.3.38	aminodeoxychorismate lyase
FORMASS-PWY	1.2.1.46	formaldehyde dehydrogenase subunit
FUC-RHAMCAT-PWY 	2.7.1.51	L-fuculokinase
FUC-RHAMCAT-PWY 	2.7.1.5	RHAMNULOKIN-MONOMER
FUC-RHAMCAT-PWY 	4.1.2.17	L-fuculose-phosphate aldolase
FUC-RHAMCAT-PWY 	4.1.2.19	rhamnulose-1-phosphate aldolase
FUC-RHAMCAT-PWY 	5.1.3.29	L-fucose mutarotase
FUC-RHAMCAT-PWY 	5.1.3.32	L-rhamnose mutarotase
FUC-RHAMCAT-PWY 	5.3.1.14 	L-rhamnose isomerase
FUC-RHAMCAT-PWY 	5.3.1.25 	L-fucose isomerase
GALACTCAT-PWY	2.7.1.178 	2-dehydro-3-deoxyglucono/2-dehydro-3-deoxygalactonokinase
GALACTCAT-PWY	2.7.1.178 	DEHYDDEOXGALACTKIN-MONOMER
GALACTCAT-PWY	4.1.2.55 	DEHYDDEOXPHOSGALACT-ALDOL-MONOMER
GALACTCAT-PWY	4.2.1.140 	GALACTONATE-DEHYDRATASE-MONOMER
GALACTCAT-PWY	4.2.1.140 	MONOMER-12925
GALACTCAT-PWY 	4.2.1.140 	MONOMER-4862
GALACT-GLUCUROCAT-PWY 	1.1.1.57	MANNONOXIDOREDUCT-MONOMER
GALACT-GLUCUROCAT-PWY 	1.1.1.58	ALTRO-OXIDOREDUCT-MONOMER
GALACT-GLUCUROCAT-PWY 	2.7.1.178 	2-keto-3-deoxygluconokinase
GALACT-GLUCUROCAT-PWY 	3.2.1.31	&beta;-D-glucuronidase
GALACT-GLUCUROCAT-PWY 	4.2.1.7	ALTRODEHYDRAT-MONOMER
GALACT-GLUCUROCAT-PWY 	4.2.1.8	MANNONDEHYDRAT-MONOMER
GALDEG-PWY	1.1.1.48	D-galactose dehydrogenase subunit
GALDEG-PWY	1.1.1.48	MONOMER-12877
GALDEG-PWY	1.1.1.48	MONOMER-12878
GALDEG-PWY	1.1.1.48	MONOMER-12879
GALDEG-PWY	3.1.1.25	lactonase
GALLATE-DEGRADATION-I-PWY 	1.13.11 	protocatechuate 3,4-dioxygenase
GALLATE-DEGRADATION-I-PWY	3.1.1.57	MONOMER-10723
GALLATE-DEGRADATION-I-PWY	4.1.3.17	HMG aldolase
GAMMAHEXCHLORDEG-PWY	1.1.1	LINCPSEPA-MONOMER
GAMMAHEXCHLORDEG-PWY	1.13.11.66	LINEPSEPA-MONOMER
GAMMAHEXCHLORDEG-PWY	1.3.1.32	MONOMER-14639
GAMMAHEXCHLORDEG-PWY	3.8.1.5	haloalkane dehalogenase
GAMMAHEXCHLORDEG-PWY	4.5.1	chlorocyclohexene dehydrochlorinase
GDPRHAMSYN-PWY	1.1.1.281 	GDP-4-dehydro-<small>D</small>-rhamnose reductase
GDPRHAMSYN-PWY	4.2.1.47	MONOMER-12853
GLCMANNANAUT-PWY 	2.7.1.60	N-acetylmannosamine kinase
GLCMANNANAUT-PWY 	4.1.3.3	N-acetylneuraminate lyase
GLCMANNANAUT-PWY 	4.1.3.3	MONOMER-18990
GLCMANNANAUT-PWY 	5.1.3.9	N-acetylmannosamine-6-phosphate 2-epimerase subunit
GLNSYN-PWY 	6.3.1.2	glutamine synthetase, type I
GLNSYN-PWY 	6.3.1.2	glutamine synthetase, type III
GLNSYN-PWY 	6.3.1.2	glutamine synthetase, type IV
GLUCARGALACTSUPER-PWY 	1.1.1.60 	TSA-REDUCT-MONOMER
GLUCARGALACTSUPER-PWY 	2.7.1.165	glycerate kinase I
GLUCARGALACTSUPER-PWY 	2.7.1.165	glycerate kinase II
GLUCARGALACTSUPER-PWY 	4.1.2.20	&alpha;-dehydro-&beta;-deoxy-D-glucarate aldolase
GLUCARGALACTSUPER-PWY 	4.2.1.42	GALACTARDEHYDRA-MONOMER
GLUCARGALACTSUPER-PWY 	5.1.2 	D-glucarate dehydratase
GLUCONEO-PWY	1.1.1.38 	malate dehydrogenase, NAD-requiring
GLUCONEO-PWY	1.1.1.40	malate dehydrogenase
GLUCONEO-PWY	4.1.1.49	PEPCARBOXYKIN-MONOMER
GLUCONEO-PWY 	4.1.1.49	phosphoenolpyruvate carboxykinase (ATP)
GLUCONSUPER-PWY	2.7.1.12	D-gluconate kinase, thermostable
GLUCOSE1PMETAB-PWY	1.1.5.2	quinoprotein glucose dehydrogenase
GLUCOSE1PMETAB-PWY	3.1.1 	GLUCONOLACT-MONOMER
GLUCOSE1PMETAB-PWY	3.1.3.10	&alpha;-D-glucose-1-phosphatase
GLUCOSE1PMETAB-PWY	3.1.3.23 	phosphosugar phosphatase
GLUCOSE1PMETAB-PWY	3.1.3.23 	sugar phosphatase
GLUCUROCAT-PWY 	1.1.1.57	D-mannonate dehydrogenase
GLUCUROCAT-PWY 	2.7.1.178 	2-keto-3-deoxygluconate kinase
GLUCUROCAT-PWY 	4.1.3.16 	MONOMER-4906
GLUCUROCAT-PWY 	4.2.1.8	D-mannonate dehydratase
GLUCUROCAT-PWY 	5.3.1.12	UXAC-MONOMER
GLUDEG-II-PWY 	1.3.8.1	MONOMER-13470
GLUDEG-II-PWY 	2.3.1.9	acetyl-CoA acetyltransferase
GLUDEG-II-PWY 	4.1.3.25 	citramalate lyase, active form
GLUDEG-II-PWY	4.2.1.150	MONOMER-13469
GLUDEG-II-PWY 	4.2.1.34	component I 
GLUDEG-II-PWY 	4.3.1.2	methylaspartase subunit
GLUDEG-II-PWY 	5.4.99.1	component S of glutamate mutase 
GLUDEG-I-PWY 	1.2.1.24	succinate semialdehyde dehydrogenase
GLUDEG-I-PWY 	2.6.1.19 	4-aminobutyrate aminotransferase, mitochondrial
GLUDEG-I-PWY	4.1.1.15	glutamate decarboxylase 1
GLUDEG-I-PWY	4.1.1.15	glutamate decarboxylase 2
GLUTAMATE-DEG1-PWY	1.4.1.2	glutamate dehydrogenase &alpha; subunit
GLUTAMATE-DEG1-PWY	1.4.1.2	glutamate dehydrogenase &beta; subunit
GLUTAMATE-DEG1-PWY 	1.4.1.2	NAD-dependent glutamate dehydrogenase
GLUTAMATE-DEG1-PWY 	1.4.1.2	NAD-glutamate dehydrogenase subunit
GLUTAMATE-DEG1-PWY	1.4.1.2	NAD-glutamate dehydrogenase subunit
GLUTAMINDEG-PWY 	1.4.1.13 	glutamate synthase
GLUTAMINDEG-PWY	3.5.1.38 	glutaminase
GLUTAMINDEG-PWY	3.5.1.38 	glutaminase A
GLUTAMINDEG-PWY	3.5.1.38 	glutaminase B monomer
GLUTAMINDEG-PWY 	3.5.1.38 	glutaminase, liver isoform
GLUTAMINDEG-PWY	3.5.1.38 	glutaminase subunit
GLUTAMINDEG-PWY 	3.5.2.3 	CAD protein
GLUTAMINDEG-PWY 	4.3.2.M2 	imidazole glycerol phosphate synthase, HisH subunit 
GLUTAMINDEG-PWY 	6.3.4.2 	CTP synthase
GLUTAMINDEG-PWY 	6.3.5.2 	GMP synthetase
GLUTAMINDEG-PWY 	6.3.5.4 	asparagine synthetase subunit
GLUTAMINDEG-PWY 	6.3.5.5 	carbamoyl phosphate synthetase, &beta; chain 
GLUTAMINEFUM-PWY	1.4.1.13 	glutamate synthase large subunit 
GLUTAMINEFUM-PWY	1.4.1.13 	glutamate synthase protomer
GLUTAMINEFUM-PWY	1.4.1.13 	glutamate synthase subunit
GLUTATHIONESYN-PWY	6.3.2.2	AT4G23100-MONOMER
GLUTATHIONESYN-PWY	6.3.2.2	&gamma;-glutamate-cysteine ligase
GLUTATHIONESYN-PWY	6.3.2.3	glutathione synthetase
GLUTATHIONESYN-PWY	6.3.2.3	GSH2 subunit
GLUT-REDOX-PWY	1.8.1.7	glutathione reductase (NADPH)
GLYCEROLMETAB-PWY	2.7.1.121	dihydroxyacetone kinase
GLYCGREAT-PWY	2.1.1.2	guanidinoacetate N-methyltransferase
GLYCGREAT-PWY	2.1.1.2	MONOMER-7641
GLYCGREAT-PWY	2.1.4.1	Glycine amidinotransferase
GLYCGREAT-PWY	2.1.4.1	L-arginine:glycine amidinotransferase subunit
GLYCINE-SYN2-PWY 	1.4.4.2 	aminomethyltransferase
GLYCINE-SYN2-PWY 	1.4.4.2 	glycine cleavage complex
GLYCINE-SYN2-PWY 	1.4.4.2	glycine decarboxylase
GLYCINE-SYN2-PWY 	1.4.4.2 	lipoamide dehydrogenase
GLYCINE-SYN2-PWY 	2.1.2.10	aminomethyltransferase
GLYCINE-SYN2-PWY 	2.1.2.10 	glycine decarboxylase
GLYCLEAV-PWY	1.4.4.2 	GCVT-MONOMER
GLYCLEAV-PWY 	1.4.4.2 	lipoamide dehydrogenase
GLYCLEAV-PWY 	1.8.1.4 	glycine cleavage system
GLYCLEAV-PWY	2.1.2.10 	glycine decarboxylase
GLYCOCAT-PWY	2.4.1.1	glycogen phosphorylase
GLYCOCAT-PWY	2.4.1.1	maltodextrin phosphorylase
GLYCOCAT-PWY	2.4.1.25	AMYLOMALT-MONOMER
GLYCOCAT-PWY 	2.7.1.8 	GLUCOKIN-MONOMER
GLYCOCAT-PWY	3.2.1.196	EG10381-MONOMER
GLYCOCAT-PWY	3.2.1.20	maltodextrin glucosidase
GLYCOGENSYNTH-PWY	2.4.1.18	GLYCOGEN-BRANCH-MONOMER
GLYCOGENSYNTH-PWY	2.4.1.21	glycogen synthase
GLYCOGENSYNTH-PWY	2.7.7.27	glucose-1-phosphate adenylyltransferase
GLYCOL-GLYOXDEG-PWY 	1.1.1.77	L-1,2-propanediol oxidoreductase
GLYCOL-GLYOXDEG-PWY 	1.1.1	tartronate semialdehyde reductase 2
GLYCOL-GLYOXDEG-PWY 	1.1.99.14	glycolate oxidase
GLYCOL-GLYOXDEG-PWY 	1.2.1.21 	aldehyde dehydrogenase A, NAD-linked
GLYCOL-GLYOXDEG-PWY 	2.3.3.9	MALSYNG-MONOMER
GLYCOL-GLYOXDEG-PWY 	4.1.1.47	glyoxylate carboligase
GLYCOLYSIS	1.2.1.12	MONOMER-548
GLYCOLYSIS	2.7.1.11	MONOMER-545
GLYCOLYSIS	2.7.1.40	pyruvate kinase subunit
GLYCOLYSIS	2.7.2.3	MONOMER-549
GLYCOLYSIS	4.1.2.13	MONOMER-546
GLYCOLYSIS	4.2.1.11	MONOMER-551
GLYCOLYSIS	5.3.1.1	MONOMER-547
GLYCOLYSIS	5.3.1.9	MONOMER-544
GLYCOLYSIS	5.4.2.11	2,3-bisphosphoglycerate-dependent phosphoglycerate mutase
GLYCOLYSIS	5.4.2.12	2,3-bisphosphoglycerate-independent phosphoglycerate mutase
GLYCOLYSIS-E-D 	1.1.1.38 	multifunctional 2-keto-3-deoxygluconate 6-phosphate aldolase and 2-keto-4-hydroxyglutarate aldolase and oxaloacetate decarboxylase
GLYCOLYSIS-E-D 	4.2.1.12	PGLUCONDEHYDRAT-MONOMER
GLYCOLYSIS-E-D 	5.3.1.9	phosphoglucose isomerase
GLYCOLYSIS-TCA-GLYOX-BYPASS 	1.1.1.37	malate dehydrogenase
GLYCOLYSIS-TCA-GLYOX-BYPASS 	1.1.1.42	isocitrate dehydrogenase
GLYCOLYSIS-TCA-GLYOX-BYPASS 	1.1.5.4	EG12069-MONOMER
GLYCOLYSIS-TCA-GLYOX-BYPASS 	2.3.3.16 	citrate synthase
GLYCOLYSIS-TCA-GLYOX-BYPASS 	2.3.3.9	MALATE-SYNTHASE
GLYCOLYSIS-TCA-GLYOX-BYPASS 	4.1.3.1	isocitrate lyase
GLYCOLYSIS-TCA-GLYOX-BYPASS 	4.2.1.2 	fumarase A
GLYCOLYSIS-TCA-GLYOX-BYPASS 	4.2.1.2 	fumarase B
GLYCOLYSIS-TCA-GLYOX-BYPASS 	4.2.1.2	fumarase C
GLYCOLYSIS-TCA-GLYOX-BYPASS 	4.2.1.99	bifunctional aconitate hydratase 2 and 2-methylisocitrate dehydratase
GLYCOLYSIS-TCA-GLYOX-BYPASS 	6.2.1.5	succinyl-CoA synthetase, &alpha; subunit 
GLYSYN-ALA-PWY	2.6.1.44	alanine--glyoxylate aminotransferase 1
GLYSYN-PWY 	2.1.2.1	serine hydroxymethyltransferase
GLYSYN-PWY	2.1.2.1	serine hydroxymethyltransferase, cytosolic
GLYSYN-PWY	2.1.2.1	serine hydroxymethyltransferase, mitochondrial
GOLPDLCAT-PWY 	1.1.1.202	1,3-propanediol dehydrogenase
GOLPDLCAT-PWY 	1.1.1.202	1,3-propanediol dehydrogenasesubunit
GOLPDLCAT-PWY 	1.1.1.202	1,3-propanediol dehydrogenase subunit
GOLPDLCAT-PWY 	1.1.1.6	glycerol dehydrogenase subunit
GOLPDLCAT-PWY 	2.7.1.29	dihydroxyacetone kinase subunit
GOLPDLCAT-PWY 	4.2.1.30	glycerol dehydratase large subunit 
GOLPDLCAT-PWY 	4.2.1.30	glycerol dehydratase subunit
HEME-BIOSYNTHESIS-II 	1.3.3.4	protoporphyrinogen oxidase
HEME-BIOSYNTHESIS-II	4.1.1.37	uroporphyrinogen decarboxylase
HEME-BIOSYNTHESIS-II 	4.99.1.1	ferrochelatase
HEXITOLDEGSUPER-PWY 	1.1.1.140	sorbitol-6-phosphate dehydrogenase
HEXITOLDEGSUPER-PWY 	1.1.1.17	mannitol-1-phosphate 5-dehydrogenase
HEXITOLDEGSUPER-PWY 	1.1.1.251	L-galactitol-1-phosphate 5-dehydrogenase
HEXITOLDEGSUPER-PWY 	4.1.2.40	tagatose-1,6-bisphosphate aldolase 1
HEXITOLDEGSUPER-PWY 	4.1.2.40	tagatose-1,6-bisphosphate aldolase 2
HEXPPSYN-PWY	2.5.1.83	hexaprenyl diphosphate synthase subunit
HISDEG-PWY	3.5.2.7	HUTIBACSU-MONOMER
HISDEG-PWY	3.5.3.8	formiminoglutamate formiminohydrolase subunit
HISDEG-PWY	4.2.1.49	HutU
HISHP-PWY	1.17.3	MONOMER-11668
HISHP-PWY 	3.5.2.7	MONOMER-11650
HISHP-PWY 	4.3.1.3	histidase subunit
HISTDEG-PWY	1.1.1.111	IMILACTDEHYDROG-MONOMER
HISTDEG-PWY	2.6.1.38	HISTTRANSAM-MONOMER
HISTSYN-PWY	2.4.2.17	AT1G09795-MONOMER
HISTSYN-PWY	2.4.2.17	AT1G58080-MONOMER
HISTSYN-PWY	2.4.2.17	MONOMER-463
HISTSYN-PWY 	2.6.1.9 	histidinol-phosphate aminotransferase
HISTSYN-PWY	2.6.1.9	MONOMER-470
HISTSYN-PWY 	3.1.3.15 	AT4G39120-MONOMER
HISTSYN-PWY	3.1.3.15	histidinol phosphate phosphatase subunit
HISTSYN-PWY	3.5.4.19 	AT1G31860-MONOMER
HISTSYN-PWY	3.5.4.19 	MONOMER-465
HISTSYN-PWY	4.2.1.19	AT3G22425-MONOMER
HISTSYN-PWY	4.2.1.19	MONOMER-469
HISTSYN-PWY	4.3.2.M2 	imidazole-glycerol-phosphate synthase
HISTSYN-PWY	5.3.1.16	AT2G36230-MONOMER
HISTSYN-PWY	5.3.1.16	MONOMER-466
HOMOSER-THRESYN-PWY	2.7.1.39	homoserine kinase monomer
HOMOSER-THRESYN-PWY	4.2.3.1	MONOMER-14634
HSERMETANA-PWY 	2.3.1.31	MONOMER-9364
HSERMETANA-PWY 	2.3.1.31	MONOMER-9423
HSERMETANA-PWY 	2.3.1.31	MONOMER-9447
HSERMETANA-PWY 	2.5.1.49	G18NG-10215-MONOMER
HSERMETANA-PWY 	2.5.1.49	MONOMER-9366
HSERMETANA-PWY 	2.5.1.49	MONOMER-9421
HYDROXYPRODEG-PWY	1.2.1 	&Delta;<sup>1</sup>-pyrroline-5-carboxylate dehydrogenase
HYDROXYPRODEG-PWY	1.5.5	proline dehydrogenase 2
HYDROXYPRODEG-PWY 	2.6.1.23 	aspartate aminotransferase 2
HYDROXYPRODEG-PWY	2.6.1.23	MONOMER-12108
HYDROXYPRODEG-PWY	4.1.3.16 	4-hydroxy-2-oxoglutarate aldolase
HYDROXYPRODEG-PWY	4.1.3.16 	MONOMER-12109
ILEUDEG-PWY	1.1.1.178	3-hydroxy-2-methylbutyryl-CoA dehydrogenase subunit
ILEUDEG-PWY	1.1.1.35 	3-hydroxy-2-methylbutyryl-CoA dehydrogenase subunit
ILEUDEG-PWY 	1.3.8.5 	2-methylacyl-CoA dehydrogenase
ILEUDEG-PWY 	2.3.1.9	mitochondrial acetyl-CoA acetyltransferase
ILEUDEG-PWY 	4.2.1.17	short chain enoyl-CoA hydratase
ILEUSYN-PWY 	1.1.1.86	acetohydroxyacid isomeroreductase
ILEUSYN-PWY 	2.2.1.6	acetohydroxyacid synthase
ILEUSYN-PWY 	3.5.99 	2-iminobutanoate/2-iminopropanoate deaminase
ILEUSYN-PWY 	4.2.1.9	dihydroxyacid dehydratase
KDO-NAGLIPASYN-PWY	2.3.1.242	PALMITOTRANS-MONOMER
KDO-NAGLIPASYN-PWY 	2.4.99.13 	KDO transferase
KDOSYN-PWY 	2.4.99.13 	(KDO)-lipid IVA 3-deoxy-D-manno-octulosonic acid transferase
KETOGLUCONMET-PWY 	1.1.1.264 	IDONDEHYD-MONOMER
KETOGLUCONMET-PWY 	1.1.1.69 	5-keto-D-gluconate 5-reductase
KETOGLUCONMET-PWY	1.1.1.79 	GhrB
KETOGLUCONMET-PWY 	2.7.1.12	D-gluconate kinase, thermosensitive
LACTOSECAT-PWY	2.7.1.144	D-tagatose-6-phosphate kinase subunit
LACTOSECAT-PWY	3.2.1.85	MONOMER-2746
LACTOSECAT-PWY	4.1.2.40	MONOMER-2744
LACTOSECAT-PWY	5.3.1.26	galactose-6-phosphate isomerase lacA subunit 
LACTOSEUTIL-PWY	1.1.99.13	lactose 3-dehydrogenase
LACTOSEUTIL-PWY 	3.2.1.23 	&alpha;-D-3-ketoglucoside 3-ketoglucohydrolase
LARABITOLUTIL-PWY	1.1.1.9	MONOMER-13244
LARABITOLUTIL-PWY	2.7.1.17	MONOMER-13243
LEU-DEG2-PWY 	1.3.8.4 	2-methylacyl-CoA dehydrogenase
LEU-DEG2-PWY	1.3.8.4	isovaleryl-CoA dehydrogenase subunit
LEU-DEG2-PWY 	1.8.1.4 	branched-chain &alpha;-keto acid dehydrogenase complex
LEU-DEG2-PWY 	2.6.1.6 	branched-chain-amino-acid aminotransferase
LEU-DEG2-PWY 	4.1.3.4	hydroxymethylglutaryl-CoA lyase
LEU-DEG2-PWY	4.2.1.18	methylglutaconyl-CoA hydratase subunit
LEU-DEG2-PWY	4.2.1.18	MONOMER-16071
LEU-DEG2-PWY	6.4.1.4	methylcrotonyl-CoA carboxylase
LEU-DEG2-PWY	6.4.1.4	methylcrotonyl-CoA carboxylase &alpha;-subunit 
LEUSYN-PWY 	1.1.1	&beta;-isopropylmalate dehydrogenase
LEUSYN-PWY 	4.2.1.35	isopropylmalate isomerase
LIPAS-PWY	3.1.1.3	AT5G04040-MONOMER
LIPASYN-PWY 	3.1.1.32 	acylhydrolase
LIPASYN-PWY	3.1.1.32	AT3G03301-MONOMER
LIPASYN-PWY	3.1.1.32	phospholipase A1
LIPASYN-PWY	3.1.1.4 	AT2G06925-MONOMER
LIPASYN-PWY	3.1.1.4 	AT4G29460-MONOMER
LIPASYN-PWY	3.1.1.4 	AtPLA IVA
LIPASYN-PWY	3.1.4.3	AT3G03530-MONOMER
LIPASYN-PWY	3.1.4.3	AT3G03540-MONOMER
LPSSYN-PWY 	2.3.1.129	UDP-N-acetylglucosamine acyltransferase
LPSSYN-PWY 	2.3.1.191	UDP-3-O-(R-3-hydroxymyristoyl)-glucosamine N-acyltransferase
LPSSYN-PWY 	2.3.1.241	LAUROYLACYLTRAN-MONOMER
LPSSYN-PWY 	2.3.1 	MYRISTOYLACYLTRAN-MONOMER
LPSSYN-PWY 	2.4.1.182	lipid A disaccharide synthase
LPSSYN-PWY 	2.4.1.44	UDP-D-galactose:(glucosyl)lipopolysaccharide-1,6-D-galactosyltransferase
LPSSYN-PWY 	2.4.1.56 	lipopolysaccharide core biosynthesis; heptosyl transferase IV; probably hexose transferase
LPSSYN-PWY 	2.4.1.58	UDP-glucose:(glucosyl)LPS &alpha;-1,2-glucosyltransferase
LPSSYN-PWY 	2.4.1.73	UDP-D-glucose:(glucosyl)LPS &alpha;-1,3-glucosyltransferase
LPSSYN-PWY 	2.4.1	lipopolysaccharide glucosyltransferase I
LPSSYN-PWY 	2.4	ADP-heptose:LPS heptosyltransferase I
LPSSYN-PWY 	2.4	lipopolysaccharide core heptosyltransferase III
LPSSYN-PWY 	2.5.1.55	3-deoxy-D-manno-octulosonate 8-phosphate synthase
LPSSYN-PWY 	2.7.1.130	TETRAACYLDISACC4KIN-MONOMER
LPSSYN-PWY 	2.7.1	lipopolysaccharide core heptose (II) kinase
LPSSYN-PWY 	2.7.1	lipopolysaccharide core heptose (I) kinase
LPSSYN-PWY 	2.7.7.38	3-deoxy-D-manno-octulosonate cytidylyltransferase
LPSSYN-PWY 	2	ADP-heptose:LPS heptosyltransferase II
LPSSYN-PWY 	3.5.1.108	UDPACYLGLCNACDEACETYL-MONOMER
LPSSYN-PWY 	3.6.1.54	EG12666-MONOMER
LPSSYN-PWY 	5.3.1.13	D-arabinose 5-phosphate isomerase
LYSINE-AMINOAD-PWY	1.1.1.87 	homo-isocitrate dehydrogenase
LYSINE-AMINOAD-PWY	1.2.1.95	L-2-aminoadipate reductase
LYSINE-AMINOAD-PWY	1.5.1.10	saccharopine dehydrogenase (NADP+, L-glutamate-forming)
LYSINE-AMINOAD-PWY	1.5.1.7	YIR034C-MONOMER
LYSINE-AMINOAD-PWY	2.3.3.14	homocitrate synthase
LYSINE-AMINOAD-PWY	4.2.1.36	homoaconitase
MALATE-ASPARTATE-SHUTTLE-PWY 	1.1.1.37	cytosolic malate dehydrogenase
MALATE-ASPARTATE-SHUTTLE-PWY	1.1.1.37	malate dehydrogenase subunit
MALATE-ASPARTATE-SHUTTLE-PWY 	1.1.1.37	mitochondrial malate dehydrogenase
MALATE-ASPARTATE-SHUTTLE-PWY 	2.6.1.1 	aspartate aminotransferase 1
MALATE-ASPARTATE-SHUTTLE-PWY	2.6.1.1	aspartate aminotransferase subunit
MALTOSECAT-PWY	2.4.1.8	maltose phosphorylase
MALTOSECAT-PWY	5.4.2.6	&beta;-phosphoglucomutase
MALTOSECAT-PWY 	5.4.2.6	MONOMER-5821
MANNOSYL-CHITO-DOLICHOL-BIOSYNTHESIS	2.4.1.117	UDP-glucose:dolichyl-phosphate glucosyltransferase
MANNOSYL-CHITO-DOLICHOL-BIOSYNTHESIS	2.4.1.131	DP-Man:Man<small>3</small>GlcNAc<small>2</small>-PP-dolichol &alpha;-1,2-mannosyltransferase
MANNOSYL-CHITO-DOLICHOL-BIOSYNTHESIS	2.4.1.141	N-acetylglucosaminyldiphosphodolichol N-acetylglucosaminyltransferase
MANNOSYL-CHITO-DOLICHOL-BIOSYNTHESIS	2.4.1.142	Dol-PP-GlcNAc2:Man transferase
MANNOSYL-CHITO-DOLICHOL-BIOSYNTHESIS	2.4.1.256	dolichyl-P-Glc:Glc<small>2</small>Man<small>9</small>GlcNAc<small>2</small>-PP-dolichol &alpha;-1,2-glucosyltransferase
MANNOSYL-CHITO-DOLICHOL-BIOSYNTHESIS	2.4.1.257 	ALG2 mannosyltransferase
MANNOSYL-CHITO-DOLICHOL-BIOSYNTHESIS	2.4.1.258	dolichyl-P-Man:Man<small>5</small>GlcNAc<small>2</small>-PP-dolichol &alpha;-1,3-mannosyltransferase
MANNOSYL-CHITO-DOLICHOL-BIOSYNTHESIS	2.4.1.259 	Dol-PP-GlcNAc2:Man6:Man transferase
MANNOSYL-CHITO-DOLICHOL-BIOSYNTHESIS	2.4.1.260	dolichyl-P-Man:Man<small>7</small>GlcNAc<small>2</small>-PP-dolichol &alpha;-1,6-mannosyltransferase
MANNOSYL-CHITO-DOLICHOL-BIOSYNTHESIS	2.4.1.265	(glucosyl)-(mannosyl)9-(N-acetylglucosaminyl)2-diphosphodolichol glucosyltransferase
MANNOSYL-CHITO-DOLICHOL-BIOSYNTHESIS	2.4.1.267	(mannosyl)9-(N-acetylglucosaminyl)2-diphosphodolichol glucosyltransferase
MANNOSYL-CHITO-DOLICHOL-BIOSYNTHESIS	2.4.1.83	dolichyl-phosphate mannose synthase
MANNOSYL-CHITO-DOLICHOL-BIOSYNTHESIS	2.4.99.18	oligosaccharyl transferase complex
MANNOSYL-CHITO-DOLICHOL-BIOSYNTHESIS	2.7.8.15	UDP-N-acetylglucosamine&mdash;dolichyl-phosphate N-acetylglucosaminephosphotransferase
M-CRESOL-DEGRADATION-PWY	1.14.13.23	MONOMER-2949
METHGLYUT-PWY 	1.1.1.283	METHYLGLYREDUCT-MONOMER
METHGLYUT-PWY 	1.1.1.75 	L-1,2-propanediol dehydrogenase / glycerol dehydrogenase
METHGLYUT-PWY 	1.1.1	L-glyceraldehyde 3-phosphate reductase
METHGLYUT-PWY 	1.1.1	methylglyoxal reductase
METHGLYUT-PWY 	1.1.1 	MONOMER0-148
METHGLYUT-PWY 	1.1.1 	MONOMER0-149
METHGLYUT-PWY 	1.1.1 	NADPH-dependent aldehyde reductase
METHGLYUT-PWY 	1.1.5	DLACTDEHYDROGFAD-MONOMER
METHGLYUT-PWY 	3.1.2.6	glyoxalase II
METHGLYUT-PWY 	4.4.1.5	glyoxalase I
METHIONINE-DEG1-PWY 	2.5.1.6	methionine adenosyltransferase 1
METHIONINE-DEG1-PWY 	2.5.1.6	methionine adenosyltransferase 2
METHYLGALLATE-DEGRADATION-PWY 	1.13.11	3-O-methylgallate 3,4-dioxygenase
METHYLGALLATE-DEGRADATION-PWY 	1.13.11 	protocatechuate 4,5-dioxygenase
MET-SAM-PWY 	2.5.1.6	methionine adenosyltransferase
MGLDLCTANA-PWY	1.1.1.78	methylglyoxal reductase &alpha; subunit 
MGLDLCTANA-PWY	1.1.1.78	MONOMER-12889
MGLDLCTANA-PWY	1.1.2.4 	mitochondrial D-lactate dehydrogenase
MGLDLCTANA-PWY 	1.2.1.3 	fatty aldehyde dehydrogenase
N2FIX-PWY	1.18.6.1	[FeMo]-nitrogenase complex
N2FIX-PWY	1.18.6.1	[FeV]-nitrogenase complex
N2FIX-PWY	1.18.6.1	[MoFe]-nitrogenase complex
NAD-BIOSYNTHESIS-II	2.7.7.1 	MONOMER-8322
NAD-BIOSYNTHESIS-II	2.7.7.1 	NadR DNA-binding transcriptional repressor and NMN adenylyltransferase
NAD-BIOSYNTHESIS-II	3.1.3.5 	acid phosphomonoesterase
NAD-BIOSYNTHESIS-II	3.1.3.5 	NAD pyrophosphatase
NAD-BIOSYNTHESIS-III	2.4.2.12	MONOMER-12230
NAD-BIOSYNTHESIS-III	2.4.2.12	MONOMER-8262
NAD-BIOSYNTHESIS-III	2.7.7.1	nicotinamide mononucleotide adenylyltransferase
NADPHOS-DEPHOS-PWY-1	1.6.1.2 	NAD(P) transhydrogenase
NADPHOS-DEPHOS-PWY-1 	2.7.1.23	NAD kinase
NADPHOS-DEPHOS-PWY-1	2.7.1.23	putative inorganic polyphosphate/ATP-NAD kinase
NADSYN-PWY 	2.4.2.19	nicotinate-nucleotide pyrophosphorylase [carboxylating]
NADSYN-PWY 	2.7.7.18 	nicotinamide mononucleotide adenylyltransferase 2
NADSYN-PWY 	2.7.7.18 	nicotinamide mononucleotide adenylyltransferase 3
NADSYN-PWY 	6.3.5.1	glutamine-dependent NAD(+) synthetase
NONMEVIPP-PWY	1.1.1.267	MONOMER-15827
NOPALINEDEG-PWY	1.5	nopaline oxidase subunit A 
NPGLUCAT-PWY	1.1.1.119 	glucose dehydrogenase subunit
NPGLUCAT-PWY	1.2.99.8	D-glyceraldehyde dehydrogenase
NPGLUCAT-PWY	2.7.1.165	glycerate 2-kinase monomer
NPGLUCAT-PWY 	2.7.1.40	pyruvate kinase subunit
NPGLUCAT-PWY 	4.2.1.140 	MONOMER-4882
OANTIGEN-PWY 	5.4.99.9	GALPMUT-MONOMER
OCTOPINEDEG-PWY	1.5	octopine oxidase subunit A 
ORN-AMINOPENTANOATE-CAT-PWY 	4.3.1.12	ornithine cyclodeaminase
ORNARGDEG-PWY 	1.2.1.19	&gamma;-aminobutyraldehyde dehydrogenase
ORNARGDEG-PWY 	1.2.1.79	succinate-semialdehyde dehydrogenase (NADP<sup>+</sup>)
ORNARGDEG-PWY 	1.2.1	&gamma;-glutamyl-&gamma;-aminobutyraldehyde dehydrogenase
ORNARGDEG-PWY 	2.6.1.29	putrescine aminotransferase
ORNARGDEG-PWY 	3.5.1.94	&gamma;-glutamyl-&gamma;-aminobutyrate hydrolase
ORNARGDEG-PWY 	3.5.3.11	agmatinase
ORNARGDEG-PWY 	4.1.1.17	ornithine decarboxylase, degradative
ORNARGDEG-PWY 	4.1.1.19	arginine decarboxylase, degradative
ORNARGDEG-PWY 	6.3.1.11	glutamate-putrescine ligase
OXIDATIVEPENT-PWY	1.1.1.44	6-phosphogluconate dehydrogenase, decarboxylating
OXIDATIVEPENT-PWY	1.1.1.44	MONOMER-6842
OXIDATIVEPENT-PWY	1.1.1.47 	GDH/6PGL endoplasmic bifunctional protein
OXIDATIVEPENT-PWY	1.1.1.49	glucose-6-phosphate 1-dehydrogenase
OXIDATIVEPENT-PWY	3.1.1.31	6-phosphogluconolactonase
P101-PWY	1.2.1.11	MONOMER-1062
P101-PWY	2.3.1.178	MONOMER-802
P101-PWY	2.6.1.46 	EctB
P101-PWY	2.7.2.4	MONOMER-1061
P101-PWY	4.2.1.108	MONOMER-803
P105-PWY	1.2.1.79	succinate-semialdehyde dehydrogenase
P105-PWY	2.3.3.9	malate synthase
P105-PWY	4.1.1.71	2-oxoglutarate decarboxylase
P105-PWY	4.2.1.2	fumarase
P108-PWY	1.1.1.37	malate dehydrogenase subunit
P108-PWY	2.1.3.1	methylmalonyl-CoA carboxyltransferase
P108-PWY	2.8.3	MONOMER-13083
P108-PWY 	5.1.99.1	methylmalonyl CoA epimerase
P108-PWY 	5.4.99.2	methylmalonyl-CoA mutase small subunit 
P122-PWY	1.1.1.28	D-lactate dehydrogenase subunit
P122-PWY	1.1.1.49	glucose-6-phosphate 1-dehydrogenase subunit
P122-PWY	1.2.1.87 	acetaldehyde dehydrogenase
P122-PWY	2.3.1.8	MONOMER-13062
P122-PWY	2.7.1.4	MONOMER-13042
P122-PWY	4.1.2.22 	phosphoketolase
P122-PWY	5.3.1.9	MONOMER-13041
P122-PWY	5.4.2.11	2,3-bisphosphoglycerate-dependent phosphoglycerate mutase
P124-PWY	4.1.2.9 	D-xylulose 5-phosphate/D-fructose 6-phosphate phosphoketolase
P124-PWY	4.1.2.9 	fructose-6-phosphate phosphoketolase &alpha; subunit 
P124-PWY 	5.3.1.9	MONOMER-15668
P125-PWY 	1.1.1 	2,3-butanediol dehydrogenase subunit
P125-PWY 	2.2.1.6	&alpha;-acetolactate synthase subunit
P125-PWY 	2.2.1.6	catabolic acetolactate synthase
P125-PWY 	4.1.1.5	&alpha;-acetolactate decarboxylase
P125-PWY 	4.1.1.5	&alpha;-acetolactate decarboxylase subunit
P142-PWY 	2.3.1.8	PTACLOS-MONOMER
P142-PWY 	2.3.1.8	PTAPROP-MONOMER
P142-PWY 	2.7.2.1	ACEKINPROP-MONOMER
P142-PWY 	2.7.2.1	acetate kinase
P161-PWY 	2.3.1.8	MONOMER-13296
P161-PWY 	2.7.2.1	MONOMER-13297
P161-PWY	4.2.1.112	MONOMER-941
P162-PWY	1.1.1.399	2-oxoglutarate reductase
P162-PWY	2.8.3.12	glutaconate CoA-transferase
P162-PWY	4.1.1.70	GcdA 
P162-PWY	4.2.1.167	HgdA 
P164-PWY	1.17.1.4	purine hydroxylase
P164-PWY 	1.17.1.4	xanthine dehydrogenase
P164-PWY	1.17.1.4	xanthine dehydrogenase &alpha; subunit 
P164-PWY 	1.2.1.2	NAD-dependent formate dehydrogenase subunit
P164-PWY	1.21.4.2	glycine reductase complex component A 
P164-PWY 	2.1.2.4	glycine formimidoyltransferase subunit
P164-PWY 	2.7.2.1	MONOMER-13133
P164-PWY 	3.5.1	MONOMER-13123
P164-PWY 	3.5.2	MONOMER-13129
P164-PWY 	3.5.4.3	MONOMER-13119
P164-PWY 	3.5.4.8	MONOMER-13121
P164-PWY 	4.3.1.4	MONOMER-13125
P165-PWY 	1.1.1.205	IMP dehydrogenase
P165-PWY 	1.17.1.4	xanthine dehydrogenase
P165-PWY 	3.5.2.17 	S-allantoin synthase
P165-PWY 	3.5.2.5	allantoinase
P165-PWY 	3.5.2.5	AT4G04955-MONOMER
P165-PWY 	3.5.3.4	allantoicase dimer
P165-PWY 	3.5.4.15	MONOMER-11783
P165-PWY 	3.5.4.17 	MONOMER-11790
P165-PWY 	3.5.4.3	MONOMER-11785
P165-PWY 	4.3.2.3	ureidoglycolate urea-lyase
P181-PWY	1.1.1.328	NAD(P)H-nicotine blue oxidoreductase
P181-PWY	1.14.13.10	2,6-dihydroxypyridine 3-monooxygenase
P181-PWY	1.2.1.79	succinate-semialdehyde dehydrogenase (NADP<small><sup>+</sup></small>)
P181-PWY	1.5.3.19	4-methylaminobutyrate oxidase (demethylating)
P181-PWY	1.5.3.21	4-methylaminobutanoate oxidase
P181-PWY	1.5.99.14	6-hydroxypseudooxynicotine dehydrogenase
P181-PWY	1.5.99.4	NdhA 
P181-PWY	3.5.1.111	2-ketoglutaramate amidase
P181-PWY	3.7.1.19	2,6-dihydroxypseudooxynicotine hydrolase
P184-PWY	1.1.1.312	&alpha;-hydroxy-&gamma;-carboxymuconate &epsilon;-semialdehyde dehydrogenase subunit
P184-PWY	1.1.1.38 	HMG aldolase
P184-PWY 	3.1.1.57	2-pyrone-4,6-dicarboxylate hydrolase
P184-PWY	3.1.1.57	MONOMER-3182
P184-PWY 	4.1.3.17	HMG aldolase
P184-PWY 	4.2.1.83	MONOMER-3183
P184-PWY	4.2.1.83	subunit of 4-oxalomesaconate hydratase
P185-PWY	1.2.1.12	glyceraldehyde-3-phosphate dehydrogenase subunit
P185-PWY	2.2.1.3	dihydroxyacetone synthase subunit
P185-PWY	2.7.1.29	dihydroxyacetone kinase subunit
P185-PWY	2.7.2.3	MONOMER-13169
P201-PWY	4.99.1	MONOMER-1002
P201-PWY	4.99.1	nitroglycerin reductase
P21-PWY	2.2.1.1	MONOMER-583
P21-PWY	5.1.3.1	MONOMER-582
P221-PWY	1.18.1.1 	alkane hydroxylase system
P221-PWY	1.18.1.1 	rubredoxin-NAD(+) reductase
P221-PWY	1.1.99	MONOMER-1081
P221-PWY 	1.2.1.3 	MONOMER-1082
P221-PWY	6.2.1.3	medium-chain acyl-CoA synthetase
P222-PWY	1.8.5.4	sulfide:quinone oxidoreductase
P23-PWY	1.1.1.37	MONOMER-11859
P23-PWY	1.1.1.42	MONOMER-11847
P23-PWY	1.3.5.4	MONOMER-11860
P23-PWY	2.3.3.8	ATP citrate synthase &alpha; subunit 
P23-PWY	2.7.9.2	pyruvate, water dikinase
P23-PWY	4.1.1.31	MONOMER-11858
P23-PWY	4.2.1.2	MONOMER-11855
P23-PWY	6.2.1.5	MONOMER-11861
P241-PWY	1.1.1 	NAD-dependent threo-isocitrate dehydrogenase
P241-PWY	4.1.3 	(R)-homocitrate synthase
P241-PWY	4.2.1.114 	methanogen homoaconitase
P261-PWY	1.1.1.337 	MONOMER-2264
P261-PWY	1.1.1.375 	L-2-hydroxycarboxylate dehydrogenase [NAD(P)+]
P261-PWY	3.1.3.71	MONOMER-2263
P261-PWY	4.1.1.79	sulfopyruvate decarboxylase
P261-PWY	4.4.1.19	ComA
P281-PWY 	1.13.11.4	gentisate 1,2-dioxygenase subunit
P281-PWY 	5.2.1.4	MONOMER-10866
P282-PWY	1.7.2	nitrite oxidoreductase
P282-PWY	1.7.2	nitrite oxidoreductase &alpha; subunit 
P283-PWY	1.12.99.6	[NiFe]-hydrogenase (membrane bound)
P283-PWY	1.12.99.6	uptake hydrogenase
P2-PWY	2.3.1.49	MONOMER-6581
P2-PWY	2.4.2.52	G6339-MONOMER
P2-PWY 	2.4.2.52	MONOMER-14203
P2-PWY	2.7.7.61	G6340-MONOMER
P2-PWY	3.1.2.16	citrate lyase deacetylase subunit
P2-PWY	6.2.1.22	[citrate [pro-3S]-lyase] ligase
P302-PWY	1.1.1.140	MONOMER-44
P302-PWY	1.1.1	MONOMER-42
P303-PWY	1.7.2.1	nitrite reductase
P303-PWY	1.7.2.1	nitrite reductase monomer
P303-PWY	1.7.2.7	hydrazine synthase
P303-PWY	1.7.2.8	hydrazine dehydrogenase
P303-PWY	1.7.2.8	hydrazine oxidoreductase A
P303-PWY	1.7.2.8	hydrazine oxidoreductase B
P321-PWY	1.3.7.8	BadD 
P341-PWY	1.2.7.6	MONOMER-11816
P341-PWY	2.7.1.146	ADP-dependent phosphofructokinase subunit
P341-PWY	2.7.1.147	ADP-dependent glucokinase subunit
P341-PWY	2.7.1.40	MONOMER-11819
P341-PWY	4.1.2.13	fructose-1,6-bisphosphate aldolase subunit
P341-PWY	4.2.1.11	enolase subunit
P341-PWY	5.3.1.1	triose phosphate isomerase subunit
P341-PWY	5.3.1.9	phosphoglucose isomerase subunit
P341-PWY	5.4.2.12	2,3-bisphosphoglycerate-independent phosphoglycerate mutase
P342-PWY	1.13.11.M3	2,3,5-trihydroxytoluene-1,2-oxygenase
P342-PWY	1.14.13.6	orcinol hydroxylase
P342-PWY	3.7.1.6	acetylpyruvate hydrolase
P343-PWY	1.13.11.37	hydroxyquinol 1,2-dioxygenase
P343-PWY	1.13.11.37	hydroxyquinol 1,2-dioxygenase II
P343-PWY	1.14.13.219	resorcinol 4-hydroxylase (NADPH)
P343-PWY	1.3.1.32	maleylacetate reductase I
P343-PWY	1.3.1.32	maleylacetate reductase II
P344-PWY	3.5.1.4	&alpha; subunit of amidase
P344-PWY 	3.5.1.4	amidase
P344-PWY 	3.5.5.7	aliphatic nitrilase subunit
P344-PWY 	4.2.1.84	nitrile hydratase
P345-PWY	3.5.1.19	nicotinamidase
P345-PWY	4.2.1.84	nitrile hydratase
P345-PWY	4.2.1	aldoxime dehydratase
P345-PWY	4.2.1	MONOMER-2282
P381-PWY 	1.13.11.79	5,6-dimethylbenzimidazole synthase
P381-PWY 	1.14.13.83	MONOMER-14881
P381-PWY 	1.14.13.83	MONOMER-84
P381-PWY	1.16.8.1	cob(II)yrinate a,c-diamide reductase
P381-PWY 	1.16.8.1	cob(II)yrinate a,c-diamide reductase subunit
P381-PWY 	1.3.1.54	MONOMER-113
P381-PWY 	2.1.1.130	CobI
P381-PWY 	2.1.1.131	MONOMER-85
P381-PWY 	2.1.1.132	precorrin-6B C<small><sup>5,15</sup></small>-methyltransferase
P381-PWY 	2.1.1.133	MONOMER-86
P381-PWY 	2.1.1.152	MONOMER-87
P381-PWY 	2.4.2.21	nicotinate-nucleotide--dimethylbenzimidazole phosphoribosyltransferase subunit
P381-PWY 	2.5.1.17	CobO
P381-PWY 	2.7.1.156 	CobP
P381-PWY 	2.7.8.26	adenosyl-cobalamin (5-phosphate) synthase
P381-PWY 	4.1.1.81	threonine-phosphate decarboxylase
P381-PWY 	4.2.1.24	porphobilinogen synthase subunit
P381-PWY 	5.4.99.61	precorrin-8x methylmutase subunit
P381-PWY 	6.3.1.10	MONOMER-142
P381-PWY 	6.3.5.10	CobQ
P381-PWY 	6.3.5.11 	cobyrinate a,c-diamide synthase
P381-PWY 	6.6.1.2	CobN 
P3-PWY	1.3.1.57	phloroglucinol reductase subunit
P3-PWY	1.3.8.1	MONOMER-21
P3-PWY	1.97.1.2	pyrogallol transhydroxylase subunit
P3-PWY	2.3.1.8	MONOMER-8
P3-PWY	2.3.1	MONOMER-7
P3-PWY	2.7.2.1	MONOMER-9
P3-PWY	2.8.3.1 	MONOMER-6
P3-PWY	3.7.1	MONOMER-5
P3-PWY	4.1.1.59	gallate decarboxylase
P3-PWY	4.2.1.150	enoyl-CoA hydratase
P401-PWY	4.4.1.9	subunit of &beta;-cyano-L-alanine synthase
P41-PWY 	1.1.1.27	MONOMER-601
P41-PWY 	1.2.1 	pyruvate dehydrogenase complex
P41-PWY 	2.3.1.8	MONOMER-602
P41-PWY 	2.7.2.1	MONOMER-603
P421-PWY	1.1.3	MONOMER-2121
P42-PWY	1.1.1.37	MONOMER-11922
P42-PWY	1.2.7.1	pyruvate oxidoreductase &alpha; subunit 
P42-PWY	1.2.7.3	2-oxoglutarate ferredoxin oxidoreductase &delta; subunit 
P42-PWY	1.3.5.4	fumarate reductase
P42-PWY	4.2.1.2	MONOMER-11923
P42-PWY	6.2.1.5	MONOMER-11928
P42-PWY	6.4.1.1	pyruvate carboxylase subunit B 
P441-PWY 	1.2.1.12	glyceraldehyde 3-phosphate dehydrogenase
P441-PWY 	2.3.1.54 	PYRUVFORMLY-CPLX
P441-PWY 	2.7.1.11 	6-phosphofructokinase II
P441-PWY 	2.7.1.40	pyruvate kinase I
P441-PWY 	2.7.1.40	pyruvate kinase II
P441-PWY 	2.7.9.2	phosphoenolpyruvate synthetase
P441-PWY	3.5.1.25	N-acetylglucosamine-6-phosphate deacetylase subunit
P441-PWY	3.5.99.6	glucosamine-6-phosphate deaminase subunit
P441-PWY 	4.1.2.13 	fructose bisphosphate aldolase class II
P441-PWY	4.1.3.3	N-acetylneuraminate lyase subunit
P441-PWY 	4.2.1.11	enolase
P441-PWY 	5.4.2.11	2,3-bisphosphoglycerate-dependent phosphoglycerate mutase
P441-PWY 	5.4.2.12	2,3-bisphosphoglycerate-independent phosphoglycerate mutase
P461-PWY 	1.1.1.140	sorbitol-6-phosphate 2-dehydrogenase
P461-PWY	1.1.1.17	MONOMER-13098
P461-PWY	1.1.1.1	MONOMER-13097
P461-PWY	1.1.1.27	L-lactate dehydrogenase subunit
P461-PWY 	1.2.1.12	glyceraldehyde-3-phosphate dehydrogenase subunit
P461-PWY	2.3.1.54	pyruvate formate-lyase subunit
P461-PWY 	2.7.1.40	pyruvate kinase subunit
P461-PWY 	4.2.1.11	enolase subunit
P481-PWY 	1.5.1.42 	2,5-diketocamphane 1,2-monooxygenase complex
P483-PWY	3.11.1.2	MONOMER-16712
P483-PWY	3.11.1.2	PhnA
P541-PWY 	2.1.1.161 	glycine/sarcosine/N,N-dimethylglycine methyltransferase
P541-PWY	2.1.1.161	MONOMER-8543
P541-PWY	2.1.1.161	MONOMER-8544
P541-PWY	2.1.1.161	MONOMER-8565
P541-PWY	2.1.1.20	MONOMER-8542
P541-PWY	2.1.1.20	MONOMER-8561
P541-PWY	2.1.1.20	MONOMER-8562
P542-PWY 	1.1.99.1	MONOMER-244
P542-PWY 	1.1.99.1	MONOMER-8621
P542-PWY 	1.2.1.8	MONOMER-246
P542-PWY	3.1.6.6	MONOMER-223
P561-PWY	1.14.13	stachydrine N-demethylase
P561-PWY	1.5.3	MONOMER-18942
P561-PWY	1.5.3	MONOMER-2221
P561-PWY	5.1.1	hydroxyproline betaine 2-epimerase / proline betaine racemase
P562-PWY	1.1.1.18	myo-inositol 2-dehydrogenase
P562-PWY 	1.1.1.18 	inositol dehydrogenase
P562-PWY	1.2.1	malonate-semialdehyde dehydrogenase
P562-PWY	2.7.1.92	5-dehydro-2-deoxygluconokinase
P562-PWY	4.1.2.29	5-dehydro-2-deoxyphosphogluconate aldolase
P562-PWY	4.2.1.44	myo-inosose-2 dehydratase
P601-PWY	1.1.1.327	5-exo-hydroxycamphor dehydrogenase
P601-PWY 	1.14.13.160	2-oxo-&Delta;<sup>3</sup>-4,5,5-trimethylcyclopentenyl acetyl coA 1,2-monooxygenase
P601-PWY 	1.14.15.1 	cytochrome P-450cam
P601-PWY 	6.2.1.38	2-oxo-&Delta;<sup>3</sup>-4,5,5-trimethylcyclopentenyl acetyl coA synthase
P621-PWY	3.5.1.117	6-aminohexanoate oligomer hydrolase polypeptide I 
P621-PWY	3.5.1.46	6-aminohexanoate linear dimer hydrolase subunit
P621-PWY	3.5.2.12	6-aminohexanoate cyclic dimer subunit
P641-PWY	1.16.1.1	mercuric ion reductase subunit
P641-PWY	4.99.1.2	MONOMER-3304
P661-PWY	1.13.11	MONOMER-1863
P661-PWY	1.14.12	hydrolase subunit of dioxin dioxygenase 
P662-PWY	1.13.11	MONOMER-1905
P662-PWY	1.14.12	dioxin dioxygenase &alpha; subunit 
P662-PWY	3.7.1.8	H1 subunit of 2-hydroxy-6-oxo-6-phenyl-2,4-hexadienoic acid hydrolase 
PANTO-PWY	1.1.1.169	MONOMER-18022
PANTO-PWY	2.1.2.11	3-methyl-2-oxobutanoate hydroxymethyltransferase
PANTO-PWY	2.7.1.33	pantothenate kinase
PANTO-PWY	6.3.2.1	pantoate &beta;-alanine ligase
PANTOSYN-PWY 	1.1.1.169	2-dehydropantoate 2-reductase
PANTOSYN-PWY 	2.1.2.11	3-methyl-2-oxobutanoate hydroxymethyltransferase
PANTOSYN-PWY 	2.7.1.24	EG12312-MONOMER
PANTOSYN-PWY 	2.7.1.33	pantothenate kinase
PANTOSYN-PWY 	2.7.7.3	phosphopantetheine adenylyltransferase
PANTOSYN-PWY 	4.1.1.11	aspartate 1-decarboxylase
PANTOSYN-PWY 	6.3.2.1	pantothenate synthetase
PANTOSYN-PWY 	6.3.2.5 	fused 4-phosphopantothenoylcysteine decarboxylase and phosphopantothenoylcysteine synthetase
PARATHION-DEGRADATION-PWY	3.1.8.1	MONOMER-3322
PCEDEG-PWY	1.97.1.8	MONOMER-1921
PCEDEG-PWY	1.97.1	MONOMER-1922
PCEDEG-PWY	1.97.1 	tetrachloroethene dehalogenase
PCEDEG-PWY	1.97.1 	tetrachloroethene-reductive dehalogenase
PCEDEG-PWY	1.97.1 	tetrachloroethene reductive dehalogenase subunit
PCEDEG-PWY	1.97.1	vinyl chloride reductase catalytic subunit 
PCPDEG-PWY	1.14.13.50	pentachlorophenol hydroxylase
PENTOSE-P-PWY 	1.1.1.44	6-phosphogluconate dehydrogenase, decarboxylating
PENTOSE-P-PWY 	1.1.1.49	GLU6PDEHYDROG-MONOMER
PENTOSE-P-PWY 	2.2.1.1	transketolase I
PENTOSE-P-PWY 	2.2.1.1	transketolase II
PENTOSE-P-PWY 	2.2.1.2	transaldolase B
PENTOSE-P-PWY 	3.1.1.31	6PGLUCONOLACT-MONOMER
PENTOSE-P-PWY	3.1.1.31	6-phosphogluconolactonase
PENTOSE-P-PWY 	5.1.3.1	ribulose-5-phosphate 3-epimerase
PENTOSE-P-PWY 	5.3.1.6	ribose-5-phosphate isomerase A
PEPTIDOGLYCANSYN-PWY 	1.3.1.98	UDP-N-acetylenolpyruvoylglucosamine reductase
PEPTIDOGLYCANSYN-PWY	2.4.1.227	NACGLCTRANS-MONOMER
PEPTIDOGLYCANSYN-PWY 	2.5.1.7	UDP-N-acetylglucosamine enolpyruvoyl transferase
PEPTIDOGLYCANSYN-PWY	2.7.8.13	PHOSNACMURPENTATRANS-MONOMER
PEPTIDOGLYCANSYN-PWY 	5.1.1.3	glutamate racemase
PEPTIDOGLYCANSYN-PWY 	6.3.2.10	UDP-NACMURALGLDAPAALIG-MONOMER
PEPTIDOGLYCANSYN-PWY 	6.3.2.13	UDP-NACMURALGLDAPLIG-MONOMER
PEPTIDOGLYCANSYN-PWY 	6.3.2.4	DALADALALIGA-MONOMER
PEPTIDOGLYCANSYN-PWY 	6.3.2.4	D-alanine-D-alanine ligase B
PEPTIDOGLYCANSYN-PWY 	6.3.2.8	UDP-N-acetylmuramate-alanine ligase
PEPTIDOGLYCANSYN-PWY 	6.3.2.9	UDP-N-acetylmuramoyl-L-alanine:D-glutamate ligase
PHENYLALANINE-DEG1-PWY	1.14.16.1	MONOMER-12067
PHENYLALANINE-DEG1-PWY	1.5.1.34 	Dihydropteridine reductase
PHENYLALANINE-DEG1-PWY	4.2.1.96	Pterin-4-alpha-carbinolamine dehydratase
PHENYLALANINE-DEG1-PWY	4.2.1.96	pterin-4-&alpha;-carbinolamine dehydratase 2
PHESYN 	2.6.1.57 	AROJBACSU-MONOMER
PHESYN	4.2.1.51	PheA
PHESYN 	5.4.99.5 	AroA
PHESYN 	5.4.99.5	AroH
PHOSLIPSYN2-PWY 	2.7.7.14	CDP-ethanolamine synthase
PHOSLIPSYN2-PWY 	2.7.7.14	ethanolamine-phosphate cytidylyltransferase
PHOSLIPSYN2-PWY 	2.7.8.5	AT3G55030-MONOMER
PHOSLIPSYN2-PWY 	2.7.8.5	MONOMER-12277
PHOSLIPSYN2-PWY 	2.7.8.5	phosphatidylglycerolphosphate synthase
PHOSLIPSYN2-PWY	2.7.8.8	phosphatidylserine synthase
PHOSLIPSYN2-PWY 	3.1.3.27	MONOMER-12278
PHOSLIPSYN2-PWY 	4.1.1.65	AT5G57190-MONOMER
PHOSLIPSYN2-PWY 	4.1.1.65	phosphatidylserine decarboxylase
PHOSLIPSYN-PWY 	2.7.7.41	CDP-diglyceride synthetase
PHOSLIPSYN-PWY 	2.7.8.5	PHOSPHAGLYPSYN-MONOMER
PHOSLIPSYN-PWY 	2.7.8.8	phosphatidylserine synthase
PHOSLIPSYN-PWY 	2.7.8	cardiolipin synthase 1
PHOSLIPSYN-PWY 	2.7.8	cardiolipin synthase 2
PHOSLIPSYN-PWY 	3.1.3.27 	PGPPHOSPHAB-MONOMER
PHOSLIPSYN-PWY 	3.1.3.27	phosphatidylglycerophosphatase A
PHOSLIPSYN-PWY 	3.1.3.27	phosphatidylglycerophosphatase C
PHOSLIPSYN-PWY 	4.1.1.65	phosphatidylserine decarboxylase, heterodimer
PHOSPHONOTASE-PWY	2.6.1.37	2-aminoethylphosphonate aminotransferase
PHOSPHONOTASE-PWY	3.11.1.1	phosphonoacetaldehyde hydrolase subunit
PLPSAL-PWY	1.4.3.5	pyridoxine-5-phosphate oxidase
PLPSAL-PWY	2.7.1.35	pyridoxal kinase
POLYAMINSYN3-PWY 	2.5.1.16	spermidine synthase subunit
POLYAMINSYN3-PWY 	2.5.1.22	spermine synthase subunit
POLYAMINSYN3-PWY 	3.5.1.53	At2g27450-monomer
POLYAMINSYN3-PWY 	3.5.1.53	N-carbamoylputrascine amidohydrolase subunit
POLYAMINSYN3-PWY 	3.5.3.12	agmatine deiminase subunit
POLYAMINSYN3-PWY 	3.5.3.12	agmatine iminohydrolase
POLYAMINSYN3-PWY 	4.1.1.17	ornithine decarboxylase subunit
POLYAMINSYN3-PWY 	4.1.1.19	arginine decarboxylase
POLYAMINSYN3-PWY 	4.1.1.19	arginine decarboxylase, biosynthetic
POLYAMINSYN3-PWY 	4.1.1.50	S-adenosylmethionine decarboxylase &alpha; subunit 
POLYISOPRENSYN-PWY 	2.5.1.31	undecaprenyl diphosphate synthase
POLYISOPRENSYN-PWY 	3.6.1.27 	undecaprenyl pyrophosphate phosphatase
POLYISOPRENSYN-PWY 	3.6.1.27	undecaprenyl pyrophosphate phosphatase
PPGPPMET-PWY	2.7.6.5	RELA-MONOMER
PPGPPMET-PWY	3.1.7.2 	SPOT-MONOMER
PPGPPMET-PWY	3.6.1.40	guanosine-5-triphosphate, 3-diphosphate pyrophosphatase
PROPFERM-PWY 	1.1.1.28	D-lactate dehydrogenase subunit
PROPFERM-PWY 	1.2.7.1	PFORPROP-MONOMER
PROPFERM-PWY 	1.3.1.95	acryloyl-CoA reductase monomer
PROPFERM-PWY 	1.4.1.2	glutamate dehydrogenase subunit
PROPFERM-PWY 	2.6.1.2	MONOMER-12754
PROPFERM-PWY 	2.8.3.1	propanoyl-CoA transferase
PROPFERM-PWY 	4.2.1.54	lactyl-CoA dehydratase
PROPIONMET-PWY 	5.1.99.1	ethylmalonyl-CoA/methylmalonyl-CoA epimerase
PROPIONMET-PWY	5.4.99.2	MONOMER-13591
PROPIONMET-PWY	6.4.1.3	propionyl-CoA carboxylase &alpha; subunit 
PROSYN-PWY	1.2.1.41	glutamate-5-semialdehyde dehydrogenase
PROSYN-PWY	1.5.1.2	pyrroline-5-carboxylate reductase
PROSYN-PWY	2.7.2.11	glutamate 5-kinase
PROUT-PWY	1.2.1.88	1-pyrroline-5-carboxylate dehydrogenase
PROUT-PWY	1.2.1.88	YHR037W-MONOMER
PROUT-PWY	1.5.5.2	YLR142W-MONOMER
PRPP-PWY 	2.1.2.3 	AICARTRANSIMPCYCLO-CPLX
PRPP-PWY 	2.4.2.17	ATP phosphoribosyltransferase
PRPP-PWY 	2.6.1.9	histidinol-phosphate aminotransferase
PRPP-PWY 	2.7.4.23	EG10723-MONOMER
PRPP-PWY 	2.7.4.3 	adenylate kinase
PRPP-PWY 	2.7.4.8 	guanylate kinase
PRPP-PWY 	2.7.6.1	PRPPSYN-MONOMER
PRPP-PWY 	3.5.4.19 	HISTCYCLOPRATPPHOS
PRPP-PWY 	4.2.1.19 	imidazoleglycerol-phosphate dehydratase / histidinol-phosphatase
PRPP-PWY 	5.3.1.16	PRIBFAICARPISOM-MONOMER
PRPP-PWY 	5.4.99.18	N<sup>5</sup>-carboxyaminoimidazole ribonucleotide mutase
PRPP-PWY 	6.3.2.6	phosphoribosylaminoimidazole-succinocarboxamide synthase
PRPP-PWY 	6.3.4.18	N<sup>5</sup>-carboxyaminoimidazole ribonucleotide synthetase
PRPP-PWY 	6.3.4.4	adenylosuccinate synthetase
PWY0-1061 	2.6.1.2	glutamate-pyruvate aminotransferase
PWY0-1061 	2.6.1.66 	valine-pyruvate aminotransferase
PWY0-1061 	2.8.1.7 	L-cysteine desulfurase
PWY0-1061 	5.1.1.1	alanine racemase 1, PLP-binding, biosynthetic
PWY0-1061 	5.1.1.1	alanine racemase 2
PWY0-1182	3.2.1.28	cytoplasmic trehalase
PWY0-1241	2.7.7.70 	fused heptose 7-phosphate kinase/heptose 1-phosphate adenyltransferase
PWY0-1241	3.1.3.82	D,D-heptose 1,7-bisphosphate phosphatase
PWY0-1241	5.1.3.20	ADP-L-glycero-D-mannoheptose-6-epimerase
PWY0-1241	5.3.1.28	D-sedoheptulose 7-phosphate isomerase
PWY0-1261	2.7.1.170	anhydro-N-acetylmuramic acid kinase
PWY0-1261	3.2.1.52	&beta;-N-acetylhexosaminidase
PWY0-1261 	3.4.17.13	L,D-carboxypeptidase A
PWY0-1261 	3.4 	G7147-MONOMER
PWY0-1261	3.5.1.28	anhydro-N-acetylmuramoyl-L-alanine amidase
PWY0-1261	3.5.1.28	N-acetyl-anhydromuramyl-L-alanine amidase
PWY0-1261	4.2.1.126	N-acetylmuramic acid 6-phosphate etherase
PWY0-1261	6.3.2.45	EG12440-MONOMER
PWY0-1264	6.3.4.15	bifunctional biotin-[acetyl-CoA-carboxylase] ligase and transcriptional repressor
PWY0-1264	6.3.4.15 	biotin holocarboxylase synthetase
PWY0-1264	6.4.1.2	acetyl-CoA carboxyltransferase, &alpha; subunit 
PWY0-1264	6.4.1.2 	biotin carboxylase
PWY0-1275 	2.8.1.8	lipoyl synthase
PWY0-1277 	3.7.1.14	2-hydroxy-6-oxononatrienedioate hydrolase
PWY0-1277 	4.1.3.39	MHPELY-MONOMER
PWY0-1277 	4.2.1.80	2-hydroxypentadienoate hydratase
PWY0-1295	2.4.2.2 	uridine phosphorylase
PWY0-1296 	5.4.2.7 	phosphoglucomutase-2
PWY0-1297 	3.5.4.4	adenosine deaminase
PWY0-1297 	4.1.2.4	DEOXYRIBOSE-P-ALD-MONOMER
PWY0-1297 	5.4.2.7 	phosphopentomutase
PWY0-1298 	2.4.2.2 	thymidine phosphorylase / uracil phosphorylase
PWY0-1300	3.2.1.170	&alpha;-mannosidase
PWY0-1301	3.2.1.22	&alpha;-galactosidase
PWY0-1305	4.1.1.15	glutamate decarboxylase A subunit
PWY0-1305	4.1.1.15	glutamate decarboxylase B subunit
PWY0-1305	4.1.1.15	YMR250W-MONOMER
PWY0-1306	1.1.1.M8	L-galactonate 5-dehydrogenase
PWY0-1306	1.1.1.M8	L-galactonate oxidoreductase
PWY0-1309	3.2.1.86	monoacetylchitobiose-6-phosphate hydrolase
PWY0-1309	3.5.1.105	chito-oligosaccharide mono-deacetylase
PWY0-1313 	6.2.1.1 	acetyl-CoA synthetase
PWY0-1313 	6.2.1.1	acetyl CoA synthetase 1
PWY0-1313 	6.2.1.1	acetyl CoA synthetase 2
PWY0-1314	2.7.1.56	1-phosphofructokinase
PWY0-1315	1.1.1.77	L-lactaldehyde reductase
PWY0-1317 	1.2.1.22 	bifunctional L-rahmnulose-phosphate aldolase/L-lactaldehyde dehydrogenase
PWY0-1319 	2.3.1.15	glycerol-3-phosphate acyltransferase
PWY0-1321 	1.1.5.6	formate dehydrogenase N, subcomplex
PWY0-1321 	1.1.5.6	formate dehydrogenase-O, &alpha; subunit 
PWY0-1321 	1.7.5.1 	nitrate reductase A
PWY0-1321 	1.7.5.1 	nitrate reductase Z
PWY0-1325 	6.3.1.1	asparagine synthetase A
PWY0-1329 	1.10.3.10	cytochrome bo terminal oxidase
PWY0-1329 	1.3.5.1	succinate:quinone oxidoreductase subcomplex
PWY0-1334 	1.10.3.14	cytochrome bd-I terminal oxidase subunit I 
PWY0-1334 	1.6.5.3 	NADH:quinone oxidoreductase I
PWY0-1336 	1.3.5.4	fumarate reductase flavoprotein 
PWY0-1337	3.1.2	thioesterase III
PWY0-1338	2.1.2.13 	fused UDP-L-Ara4N formyltransferase and UDP-GlcA dehydrogenase
PWY0-1338	2.4.2.43	4-amino-4-deoxy-L-arabinose (L-Ara4N) transferase
PWY0-1338	2.4.2.53	undecaprenyl phosphate-L-Ara4FN transferase
PWY0-1338	2.6.1.87	G7166-MONOMER
PWY0-1355 	1.7.2.3	trimethylamine N-oxide reductase TorCA
PWY0-1356 	1.8.5.3	dimethyl sulfoxide reductase, chain A 
PWY0-1415 	1.3.98.3	coproporphyrinogen III dehydrogenase
PWY0-1433 	1.5.1.50 	G6862-MONOMER
PWY0-1433	5.1.99.7	dihydroneopterin triphosphate 2-epimerase
PWY0-1465 	1.1.1.83 	D-malate / 3-isopropylmalate dehydrogenase (decarboxylating)
PWY0-1466	3.2.1.28	periplasmic trehalase
PWY0-1471	1.1.1.298	3-hydroxy acid dehydrogenase
PWY0-1471	1.1.1.298	predicted malonic semialdehyde reductase
PWY0-1471	1.14.99.46	pyrimidine oxygenase
PWY0-1477	4.3.1.7	ethanolamine ammonia-lyase, &alpha; subunit 
PWY0-1479	2.7.7.56	RNase PH
PWY0-1479	3.1.13.1 	polynucleotide phosphorylase
PWY0-1479	3.1.13.1	ribonuclease II
PWY0-1479	3.1.13.1	RNase BN
PWY0-1479	3.1.13.5 	RNase D
PWY0-1479	3.1.13.5 	RNase T
PWY0-1479	3.1.26.12	ribonuclease E
PWY0-1479	3.1.26.5	RNase P
PWY0-1507	2.6.1.62	7,8-diaminopelargonic acid aminotransferase
PWY0-1507	2.8.1.6	biotin synthase
PWY0-1517 	2.7.1 	6-phosphofructokinase I
PWY0-1533 	3.1.4.55	5-phospho-&alpha;-D-ribosyl 1,2-cyclic phosphate phosphodiesterase
PWY0-1544 	1.5.5.2 	fused PutA DNA-binding transcriptional repressor / proline dehydrogenase / 1-pyrroline-5-carboxylate dehydrogenase
PWY0-1546	3.4.13.18 	peptidase D
PWY0-1546	5.1.1.20	L-Ala-D/L-Glu epimerase
PWY0-1568 	1.6.5.9 	NADH:quinone oxidoreductase II
PWY0-1569	2.3.1.245	3-hydroxy-2,4-pentadione 5-phosphate thiolase
PWY0-1569	2.7.1.189	autoinducer-2 kinase
PWY0-1569	5.3.1.32	phospho-AI-2 isomerase
PWY0-1576 	1.12.99.6	hydrogenase 2
PWY0-1577 	1.12.99.6	hydrogenase 1, oxygen tolerant hydrogenase
PWY0-1581 	1.1.5.3	glycerol-3-phosphate dehydrogenase, anaerobic
PWY0-1584 	1.1.5.3	glycerol-3-phosphate dehydrogenase, aerobic
PWY0-1584	1.7.5.1	periplasmic nitrate reductase
PWY0-1586	2.4.1.129	biosynthetic peptidoglycan transglycosylase
PWY0-1586	2.4.1.129	G7322-MONOMER
PWY0-1586	2.4.1.129 	peptidoglycan glycosyltransferase - peptidoglycan D,D transpeptidase
PWY0-1586	2.4.1.129 	peptidoglycan glycosyltransferase / peptidoglycan D,D-transpeptidase
PWY0-1586	3.4.17.13	L,D-transpeptidase ErfK
PWY0-1586	3.4.17.13	L,D-transpeptidase LdtD
PWY0-1586	3.4.17.13	L,D-transpeptidase LdtE
PWY0-1586	3.4.17.13	L,D-transpeptidase YbiS
PWY0-1586	3.4.17.13	L,D-transpeptidase YcfS
PWY0-1586	3.4 	EG12867-MONOMER
PWY0-1586	3.5.2.6 	D-alanyl-D-alanine carboxypeptidase, fraction A; penicillin-binding protein 5
PWY0-1587	2.3.1.234	N<sup>6</sup>-L-threonylcarbamoyladenine synthase
PWY0-1587	2.7.7.87	threonylcarbamoyl-AMP synthase
PWY0-162 	2.7.4.6	nucleoside diphosphate kinase
PWY0-166 	1.17.4.1	ribonucleoside diphosphate reductase 1
PWY0-166 	2.1.1.45	thymidylate synthase
PWY0-166 	2.7.4.12 	DTMPKI-MONOMER
PWY0-166 	2.7.4.12 	thymidylate kinase
PWY0-166 	3.5.4.13	dCTP deaminase
PWY0-166 	3.6.1.23 	deoxyuridine triphosphatase
PWY0-166 	3.6.1.23 	pyrimidine deoxynucleoside triphosphate pyrophosphohydrolase
PWY0-166 	3.6.1.9 	nucleoside triphosphate pyrophosphohydrolase
PWY-0	2.3.1.57	MONOMER-27
PWY0-301	3.1.1	L-ascorbate 6-phosphate lactonase
PWY0-301	4.1.1.85	3-keto-L-gulonate 6-phosphate decarboxylase
PWY0-301	5.1.3.22	L-xylulose 5-phosphate 3-epimerase
PWY0-301	5.1.3.4	L-ribulose 5-phosphate 4-epimerase
PWY0-321	1.14.13.149	ring 1,2-phenylacetyl-CoA epoxidase
PWY0-321	2.3.1.174 	&beta;-ketoadipyl-CoA thiolase
PWY0-321	4.2.1.17	2,3-dehydroadipyl-CoA hydratase
PWY0-321	5.3.3.18	1,2-epoxyphenylacetyl-CoA isomerase
PWY0-321	6.2.1.30	phenylacetate-CoA ligase (aerobic)
PWY0-381 	2.7.1.30	glycerol kinase
PWY0-381 	3.1.4.46 	glycerophosphodiester phosphodiesterase, cytosolic
PWY0-381 	3.1.4.46	glycerophosphoryl diester phosphodiesterase, periplasmic
PWY0-41	1.1.1.154 	MONOMER-13523
PWY0-41	1.1.1.350	ureidoglycolate dehydrogenase
PWY0-41	2.1.3.5	MONOMER-13524
PWY0-41 	3.5.2.5	allantoinase
PWY0-41 	3.5.3.26	S-ureidoglycine aminohydrolase
PWY0-41 	3.5.3.9	allantoate amidohydrolase
PWY0-42	2.3.3.5 	2-methylcitrate synthase
PWY0-42	2.3.3.5	MONOMER-63
PWY0-42	4.1.3.30	2-methylisocitrate lyase
PWY0-42	4.1.3.30	2-methylisocitrate lyase subunit
PWY0-42	4.2.1.79	G6199-MONOMER
PWY0-42	4.2.1.79	MONOMER-64
PWY0-42	4.2.1.99	aconitate hydratase 1
PWY0-42	4.2.1.99	aconitate hydratase 2
PWY0-42 	6.2.1.17 	acetyl-CoA synthetase (AMP-forming)
PWY0-42	6.2.1.17	MONOMER-62
PWY0-42	6.2.1.17	propionyl-CoA synthetase
PWY0-43	2.8.3	G7517-MONOMER
PWY0-43	4.1.1.M5	methylmalonyl-CoA decarboxylase
PWY0-43 	5.4.99.2	methylmalonyl-CoA mutase
PWY0-44	2.7.1.55	EG11956-MONOMER
PWY0-44	5.1.3	EG11957-MONOMER
PWY0-44 	5.3.1 	allose-6-phosphate isomerase / ribose-5-phosphate isomerase B
PWY0-461 	4.1.1.18	MONOMER-15076
PWY0-501	2.3.1.181	lipoyl(octanoyl) transferase
PWY0-522	6.3.1.20	lipoate-protein ligase
PWY0-522 	6.3.1.20	lipoate-protein ligase A
PWY0-541	2.1.1.79	cyclopropane fatty acyl phospholipid synthase
PWY0-662	2.7.6.1	ribose-phosphate pyrophosphokinase 1
PWY0-662	2.7.6.1	ribose-phosphate pyrophosphokinase 2
PWY0-662	2.7.6.1	ribose-phosphate pyrophosphokinase 3
PWY0-781 	1.17.1.8	4-hydroxy-tetrahydrodipicolinate reductase
PWY0-781 	1.2.1.11	aspartate semialdehyde dehydrogenase
PWY0-781 	1.5.99 	L-aspartate oxidase
PWY0-781 	2.1.1.13 	cobalamin-dependent methionine synthase
PWY0-781 	2.1.1.14	HOMOCYSMET-MONOMER
PWY0-781 	2.3.1.117	tetrahydrodipicolinate succinylase
PWY0-781 	2.3.1.46	homoserine O-succinyltransferase
PWY0-781 	2.4.2.19	quinolinate phosphoribosyltransferase (decarboxylating)
PWY0-781 	2.5.1.48 	O-succinylhomoserine(thiol)-lyase / O-succinylhomoserine lyase
PWY0-781 	2.5.1.72	quinolinate synthase
PWY0-781 	2.7.2.4 	aspartate kinase / homoserine dehydrogenase
PWY0-781 	2.7.2.4	aspartate kinase III
PWY0-781 	2.7.7.18	NICONUCADENYLYLTRAN-MONOMER
PWY0-781 	3.5.1.18	N-succinyl-L-diaminopimelate desuccinylase subunit
PWY0-781 	4.1.1.20	diaminopimelate decarboxylase
PWY0-781 	4.3.3.7	4-hydroxy-tetrahydrodipicolinate synthase
PWY0-781 	5.1.1.7	diaminopimelate epimerase
PWY0-845 	1.1.1.290	erythronate-4-phosphate dehydrogenase
PWY0-845 	1.2.1.72	erythrose 4-phosphate dehydrogenase
PWY0-845 	2.2.1.7	1-deoxyxylulose-5-phosphate synthase
PWY0-845 	2.6.99.2	pyridoxine 5-phosphate synthase
PWY0-845 	2.7.1.35	pyridoxal kinase 2
PWY0-845 	2.7.1.35 	pyridoxal kinase I / hydroxymethylpyrimidine kinase
PWY0-862	2.3.1.41	MONOMER-19312
PWY0-862	5.3.3.14 	3-hydroxyacyl-[acyl-carrier-protein] dehydratase/isomerase
PWY0-862	5.3.3.14	trans-2-decenoyl-[acyl-carrier-protein] isomerase
PWY0-881 	2.3.1.86 	MALONYL-COA-ACP-TRANSACYL-MONOMER
PWY0-881 	6.3.4.14 	 
PWY0-901	2.7.9.3	selenide, water dikinase
PWY0-901	2.9.1.1	selenocysteine synthase
PWY-1001	2.4.1.12	cellulose synthase
PWY-101	1.10.3.9	photosystem II
PWY-101	1.97.1.12	photosystem I
PWY-102	1.14.11.13	MONOMER-11640
PWY-102 	1.14.11.15 	gibberellin 2&beta;,3&beta;-hydroxylase
PWY-102	1.14.11 	gibberellin 2-oxidase
PWY-102	1.14.11	gibberellin 2-oxidase
PWY-1061	2.1.1	homogalacturonan methyltransferase
PWY-1061	2.1.1	MONOMER-2485
PWY-1061	2.4.1.43	&alpha; 1,4-galacturonosyltransferase
PWY-1081	3.1.1.11	MONOMER-14860
PWY-1081	3.1.1.11	MONOMER-16134
PWY-1081	3.1.1.11	pectin methylesterase
PWY-1081	3.2.1.15	MONOMER-14861
PWY-1081	3.2.1.15	MONOMER-14863
PWY-1081	3.2.1.15	polygalacturonase
PWY-1121	1.1.1	&omega;-hydroxyacid dehydrogenase
PWY-1121 	1.1.1	&omega;-hydroxy fatty acid &omega;-alcohol dehydrogenase
PWY-1121	1.14.13	very-long-chain fatty acid &omega;-hydroxylase
PWY-1121 	6.2.1.34 	4-coumarate coenzyme A ligase 1
PWY-1121 	6.2.1.3 	long-chain acyl-CoA synthetase 1
PWY-1121 	6.2.1.3 	long-chain acyl-CoA synthetase 3
PWY-1121 	6.2.1.3 	long-chain acyl-CoA synthetase 4
PWY-1121 	6.2.1.3 	long chain acyl-CoA synthetase 5
PWY-1121 	6.2.1.3 	long-chain acyl-CoA synthetase 8
PWY-1121 	6.2.1.3 	long-chain acyl-CoA synthetase 9
PWY-116	2.4.1.111	AT5G26310-MONOMER
PWY-116	2.4.1.111	AT5G66690-MONOMER
PWY-116	2.4.1.111	MONOMER-11622
PWY-116	2.4.1.111	MONOMER-594
PWY-116 	2.4.1 	UDPG:coniferyl alcohol glucosyltransferase
PWY-116	3.2.1.126	glucosidase subunit I 
PWY-116	3.2.1.126	MONOMER-11624
PWY-116	3.2.1.126	MONOMER-11625
PWY-116	3.2.1.126	MONOMER-11629
PWY-1186 	2.3.3	methylthioalkylmalate synthase
PWY1-2	1.4.1.1	alanine dehydrogenase subunit
PWY1-2	1.4.1.1	L-alanine dehydrogenase subunit
PWY-12	2.3.1.216	PCS subunit
PWY-1263	2.6.1.77	MONOMER-2761
PWY-1264	1.4.99.2	taurine dehydrogenase subunit
PWY-1269	2.5.1.55	3-deoxy-8-phosphooctulonate synthase subunit
PWY-1269	2.7.7.38	MONOMER-11935
PWY-1269	3.1.3.45	3-deoxy-D-manno-octulosonate 8-phosphate phosphatase
PWY-1269	5.3.1.13	arabinose-5-phosphate isomerase
PWY-1281	2.3.3.15	sulfoacetaldehyde acetyltransferase
PWY1-3	1.1.1.36	acetoacetyl-CoA reductase subunit
PWY1-3 	1.1.1.36	(R)-3-hydroxybutyryl-CoA dehydrogenase
PWY-13 	2.3.1.151	MONOMER-11400
PWY-13 	2.3.1.220	BPS subunit
PWY1-3 	2.3.1.9	MONOMER-13585
PWY1-3	2.3.1	poly-&beta;-hydroxybutyrate polymerase subunit
PWY-1361	1.14.13.208	benzoyl-CoA 2,3-epoxidase
PWY-1361	1.2.1.77	3,4-dehydroadipyl-CoA semialdehyde dehydrogenase monomer
PWY-1361	4.1.2.44	2,3-epoxybenzoyl-CoA dihydrolase
PWY-13 	6.2.1.25	MONOMER-11385
PWY-13 	6.2.1.37	MONOMER-11399
PWY-13 	6.2.1.37 	MONOMER-11580
PWY-1422 	2.1.1.295	MONOMER-13899
PWY-1422	2.1.1.95	tocopherol methyltransferase
PWY-1422	2.5.1.115	AT2G18950-MONOMER
PWY-1422	2.5.1.115	MONOMER-13902
PWY-1422	5.5.1.24	tocopherol cyclase
PWY-142 	4.1.99.11	benzylsuccinate synthase
PWY-1541 	1.14.11.17	taurine dioxygenase
PWY-1541 	1.4.99.2	taurine dehydrogenase
PWY-1541 	2.3.1.8	phosphate acetyltransferase
PWY-1541 	2.3.3.15	sulfoacetaldehyde acetyltransferase
PWY-1541 	2.6.1.55	MONOMER-2966
PWY-1541 	2.6.1.77	taurine-pyruvate aminotransferase
PWY-1541 	2.6.1.77	taurine:pyruvate aminotransferase
PWY-1622	1.1.1.37	malate dehydrogenase subunit
PWY-1622	1.1.1.81 	hydroxypyruvate reductase subunit
PWY-1622 	2.1.2.1	serine hydroxymethyltransferase
PWY-1622 	2.1.2.1	serine hydroxymethyltransferase subunit
PWY-1622	2.6.1.45	serine-glyoxylate aminotransferase
PWY-1622	2.7.1.165	glycerate 2-kinase
PWY-1622	4.1.1.31	phosphoenolpyruvate carboxylase subunit
PWY-1622 	4.1.3.24	malyl coenzyme A lyase
PWY-1622	4.2.1.11	enolase
PWY-1622	5.4.2.11	2,3-bisphosphoglycerate-dependent phosphoglycerate mutase
PWY-1622 	6.2.1.9	malate thiokinase &alpha; subunit 
PWY-1641	1.14.13.25	soluble methane monooxygenase
PWY-1701 	1.4.9.1	methylamine dehydrogenase large subunit 
PWY-1723	1.5.1	MtdB
PWY-1781	1.2.1	malonate semialdehyde dehydrogenase
PWY-1781	2.6.1.18	&beta;-alanine transaminase
PWY-1782 	2.3.1.72	IAA-myo-inositol synthase
PWY-1782 	2.4.1.121	IAA-glucose synthase
PWY-1782 	2.4.1.156	indol-3-ylacetyl-myo-inositol galactoside synthase
PWY-1801	1.1.1.284 	AdhI
PWY-1801	1.1.1.284 	formaldehyde dehydrogenase
PWY-1801	3.1.2.12	S-formylglutathione hydrolase
PWY-1801 	3.1.2.12 	S-formylglutathione hydrolase / S-lactoylglutathione hydrolase
PWY-1801	3.1.2.12	S-formylglutathione hydrolase subunit
PWY-181	1.1.1.29	hydroxypyruvate reductase
PWY-181	1.1.3.15	peroxisomal (S)-2-hydroxy-acid oxidase
PWY-181	2.6.1.44 	alanine transaminase
PWY-181	2.6.1.45	AT2G13360-MONOMER
PWY-181	2.6.1.4	glyoxylate aminotransferase
PWY-181	2.7.1.31	AT1G80380-MONOMER
PWY-181	3.1.3.18	phosphoglycolate phosphatase
PWY-1822	3.5.1	indole-3-acetyl-ala hydrolase
PWY-1822	3.5.1	indole-3-acetyl-leu hydrolase
PWY-1861 	4.1.2.43	hexulose-6-phosphate synthase
PWY-1861 	4.1.2.43	hexulose-6-phosphate synthase subunit
PWY-1861 	5.3.1.27	6-phospho-3-hexuloisomerase subunit
PWY-1861 	5.3.1.27	6-phspho-3-hexuloisomerase subunit
PWY-1881	1.2.1.2	formate dehydrogenase
PWY-1881	1.2.1.2	NAD-dependent formate dehydrogenase &alpha; subunit 
PWY-1881	1.2.1.2	NAD-dependent formate dehydrogenase subunit
PWY-1881	1.2.1.2	NAD-linked formate dehydrogenase protomer
PWY-1882 	1.1.1.284 	FlhA
PWY-1882 	1.1.2.7	cytochrome c-dependent methanol dehydrogenase
PWY-1882 	1.1.2.7	methanol dehydrogenase (cytochrome c)
PWY-1882 	1.2.1.2	NAD-dependent formate dehydrogenase &alpha; subunit 
PWY-1882 	1.4.9.1	methylamine dehydrogenase large subunit 
PWY-1882 	1.5.1	MtdB
PWY-1882 	1.5.1 	NADP-dependent methylene tetrahydromethanoptrin dehydrogenase subunit
PWY-1882 	3.1.2.12	FghA
PWY-1882 	3.5.4.27	methenyltetrahydropterin cyclohydrolase subunit
PWY-1882 	4.2.1.147	formaldehyde-activating enzyme
PWY-1882 	4.4.1.22	Gfa
PWY-1901	1.21.3.6	aureusidin synthase
PWY-1901	2.4.1.286	chalcone 4-O-glucosyltransferase
PWY-1901 	2.4.1	6-deoxychalcone 4-glucosyltransferase
PWY-1921	3.1.1	1-O-IAGlu hydrolase
PWY-1921	3.1.1	6-O-IAGlu hydrolase
PWY1A0-6325	1.1.1	ActVI-1
PWY1A0-6325	1.14.14	ActVA-5
PWY1A0-6325	1.3.1	ActIII
PWY1A0-6325	2.3.1	type II minimal act PKS
PWY1A0-6325	3.1.2 	ActIV
PWY1A0-6325	4.2.1	ActVII
PWY1A0-6325	5.5.1	ActVI-3
PWY1F-353	1.14.15.7	chloroplastic choline monooxygenase
PWY1F-353	1.2.1.8	betaine aldehyde dehydrogenase
PWY1F-467 	1.14.13.11	cinnamate 4-hydroxylase
PWY1F-823 	1.1.1.219	dihydroflavonol 4-reductase
PWY1F-823	1.1.1.219	dihydroflavonol 4-reductase
PWY1F-823	1.1.1.219	dihydroflavonol-4-reductase
PWY1F-823 	1.14.11.9	flavanone 3&beta;-hydroxylase
PWY1F-823	1.14.11.9	MONOMER-11942
PWY1F-FLAVSYN	2.3.1.170	MONOMER-4908
PWY1F-FLAVSYN	2.3.1.170	MONOMER-4921
PWY1F-FLAVSYN	2.3.1.170	MONOMER-4923
PWY1F-FLAVSYN	2.3.1.74	naringenin chalcone synthase
PWY1F-FLAVSYN	5.5.1.6	chalcone isomerase
PWY1G-0	2.3.1.189	G185E-4968-MONOMER
PWY1G-0	2.4.1.250	D-inositol-3-phosphate glycosyltransferase
PWY1G-0	3.5.1.103	N-acetyl-1-D-myo-inositol-2-amino-2-deoxy-&alpha;-D-glucopyranoside deacetylase
PWY1G-0	6.3.1.13	L-cysteine : 1 D-myo-inositol 2-amino-2-deoxy-&alpha;-D-glucopyranoside ligase
PWY1G-126	1.8.1.15	mycothione reductase subunit
PWY1G-1	3.5.1.115	G185E-5244-MONOMER
PWY1G-1	3.5.1.115	MONOMER-9708
PWY1G-170	1.1.1.306	mycothiol-dependent formaldehyde dehydrogenase subunit
PWY-2002	2.1.1.150	MONOMER-4964
PWY-2002	2.1.1.150	MONOMER-7441
PWY-2002	2.1.1.154	MONOMER-4981
PWY-2002	5.5.1.6	chalcone isomerase
PWY-2002	5.5.1.6	MONOMER-4386
PWY-2002 	5.5.1.6	MONOMER-4742
PWY-2002 	5.5.1.6	MONOMER-4743
PWY-2002 	5.5.1.6	MONOMER-4744
PWY-2002	5.5.1.6	MONOMER-4745
PWY-2002	5.5.1.6	MONOMER-4761
PWY-2002	5.5.1.6	MONOMER-4762
PWY-2055 	1.1.1.348	MONOMER-7261
PWY-2055 	1.14.13.28	MONOMER-6705
PWY-2055 	1.14.13.28	MONOMER-7281
PWY-2055 	1.14.13.85	MONOMER-6708
PWY-2055 	1.14.13.85	MONOMER-6709
PWY-2055 	1.14.13.85	MONOMER-6710
PWY-2055 	2.1.1.212	MONOMER-5148
PWY-2055 	2.1.1.212	MONOMER-7461
PWY-2055 	2.1.1.46	MONOMER-5152
PWY-2055 	2.5.1.36	MONOMER-6706
PWY-2055 	2.5.1.36	MONOMER-6707
PWY-2083 	1.14.11.22	flavone synthase I
PWY-2083 	1.14.13.52	isoflavone 3-hydroxylase
PWY-2083 	1.14.13	flavone synthase II
PWY-2083 	2.1.1.46 	flavonoid 4-O-methyltransferase
PWY-2 	1.2.1 	betaine aldehyde / aminoaldehyde dehydrogenase
PWY-2	1.4.3.10	copper-containing amine oxidase
PWY-2	1.4.3.10	putrescine oxidase subunit
PWY-2	1.4.3.21 	copper-containing amine oxidase
PWY-2161 	6.3.2.17 	bifunctional folylpolyglutamate synthetase / dihydrofolate synthetase
PWY-2161	6.3.2.17	folylpoly-&alpha;-glutamate synthetase
PWY-2161	6.3.2.17	folylpolyglutamate synthetase
PWY-2161	6.3.2.17	tetrahydrofolate synthase
PWY-2161B	3.4.19.9	&gamma;-glutamyl hydrolase
PWY-2201 	1.5.1.20	5,10-methylenetetrahydrofolate reductase
PWY-2201 	1.5.1.20	methylene tetrahydrofolate reductase
PWY-2201 	1.5.7.1	ferredoxin-dependent methylenetetrahydrofolate reductase
PWY-2201	2.1.1.19	MONOMER-5201
PWY-2201	3.5.1.10	formyltetrahydrofolate deformylase
PWY-2201 	3.5.4.9	folate cyclohydrolase subunit
PWY-2201 	3.5.4.9 	FolD
PWY-2201 	6.3.3.2	5-formyltetrahydrofolate cyclo-ligase
PWY-2201	6.3.3.2	MONOMER-5261
PWY-2201 	6.3.4.3	formate-tetrahydrofolate ligase subunit
PWY-2221 	1.1.1.119 	glucose dehydrogenase subunit
PWY-2221	1.1.1.119 	glucose dehydrogenase subunit
PWY-2221	1.2.1.3	NAD<sup>+</sup>-dependent glyceraldehyde-3-phosphate dehydrogenase subunit
PWY-2221	2.7.1.178 	MONOMER-4885
PWY-2221 	4.1.2.55 	2-keto-3-deoxy-6-phosphogluconate aldolase
PWY-2229 	1.1.1.348	(-)-sophorol reductase
PWY-2229 	1.3.1.45	MONOMER-6024
PWY-2229 	4.2.1.139	cis-(-)-7,2-dihydroxy-4,5-methylenedioxyisoflavanol dehydratase
PWY-2261	1.6.5	MONOMER-5129
PWY-2261	1.6.5	MONOMER-5130
PWY-2261	1.8.5.1	MONOMER-5131
PWY-2261	1.8.5.1	MONOMER-5132
PWY-2301	3.1.3.25	AT1G31190-MONOMER
PWY-2301 	3.1.3.25	inositol monophosphatase 1
PWY-2301	3.1.3.25	inositol monophosphatase subunit
PWY-2301 	3.1.3.25 	L-galactose-1-phosphate phosphatase
PWY-2301 	3.1.3.25	myo-inositol monophosphatase 2
PWY-2301 	5.5.1.4	inositol-1-phosphate synthase subunit
PWY-2343	2.3.1.115	MONOMER-5302
PWY-2343	2.4.1.170	MONOMER-5242
PWY-2343	3.2.1	MONOMER-5202
PWY-2343	3.2.1	MONOMER-5223
PWY-2345	2.3.1.115	MONOMER-5301
PWY-2345	2.4.1.170	MONOMER-5241
PWY-2381	1.7.1	MONOMER-2123
PWY-2381	4.3.3	MONOMER-2124
PWY-241 	1.1.1.40	NADP malic enzyme
PWY-241 	1.1.1.82	malate dehydrogenase (NADP -linked)
PWY-241 	2.7.9.1	pyruvate orthophosphate dikinase
PWY-241 	4.1.1.31	phosphoenolpyruvate carboxylase
PWY-241 	4.2.1.1	&beta;-carbonic anhydrase
PWY-2421	1.14.13	MONOMER-18108
PWY-2421	1.14.13	MONOMER-18109
PWY-2421	1.1	MONOMER-8442
PWY-2421	3.5.2.20	MONOMER-8443
PWY-2463 	1.14.13.53	isoflavone 2-hydroxylase
PWY-2463 	1.14.13.53	MONOMER-5945
PWY-2464	1.14.13.52	isoflavone 3-hydroxylase
PWY-2467 	1.14.13	MONOMER-6103
PWY-2467 	2.1.1.270	MONOMER-6105
PWY-2467 	4.2.1	isoflavene synthase
PWY-2501	1.13.11	fatty acid &alpha;-dioxygenase
PWY-2501	1.2.1.3 	70-kDa subunit 
PWY-2503 	1.14.12.10 	benzoate 1,2-dioxygenase
PWY-2504 	1.13.11.11 	tryptophan dioxygenase subunit
PWY-2504 	1.2.1.28 	MONOMER-2661
PWY-2504 	3.7.1.3	kynureninase subunit
PWY-2504 	4.1.1.44	subunit of 4-carboxymuconolactone decarboxylase
PWY-2504 	5.5.1.2	subunit of 3-carboxymuconate cycloisomerase
PWY-2541	1.1.1 	3&beta;-hydroxysteroid dehydrogenase/C4-decarboxylase
PWY-2541	1.1.1	sterone reducatase
PWY-2541	1.14.13	4,4-dimethyl-9&beta;,19-cyclopropylsterol-4&alpha;-methyl oxidase
PWY-2541	1.14.13	4-&alpha;-methyl-&Delta;7-sterol-4&alpha;-methyl oxidase
PWY-2541	1.14.13.70	MONOMER-7161
PWY-2541	1.14.13.70	obtusifoliol 14&alpha;-demethylase
PWY-2541 	1.14.14.17	squalene monooxygenase
PWY-2541	1.14.19.20	&Delta;<sup>7</sup>-sterol-C5-desaturase
PWY-2541	1.14.19.41	AT2G34500-MONOMER
PWY-2541	1.14.19.41	C-22 sterol desaturase
PWY-2541	1.3.1.21	sterol &Delta;<sup>7</sup> reductase
PWY-2541	1.3.1.72 	&Delta;<sup>5</sup>-sterol-&Delta;<sup>24</sup>- reductase
PWY-2541	1.3.1 	&Delta;<sup>8,14</sup>sterol &Delta;14-sterol reductase
PWY-2541	2.1.1.143	AT1G76090-MONOMER
PWY-2541	2.1.1.143	MONOMER-15391
PWY-2541	2.1.1.143 	sterol methyltransferase
PWY-2541	2.1.1.41	sterol methyltransferase
PWY-2541 	2.5.1.21 	squalene synthase
PWY-2541	5.3.3.5	&Delta;8-&Delta;7-sterol isomerase
PWY-2541	5.4.99.8	cycloartenol synthase
PWY-2541	5.5.1 	cycloeucalenol cycloisomerase
PWY-2561	2.3.1.115	MONOMER-6327
PWY-2561	2.4.1.170	MONOMER-6326
PWY-2582 	1.14.13 	AT4G36380-MONOMER
PWY-2582 	1.14.13	steroid 22&alpha;-hydroxylase
PWY-2601	1.1.1	MONOMER-6042
PWY-2622	5.4.99.16	MONOMER-5647
PWY-2622	5.4.99.16	MONOMER-5661
PWY-2622	5.4.99.16	MONOMER-6161
PWY-2622	5.4.99.16	trehalose synthase subunit
PWY-2661	3.2.1.141	MONOMER-5702
PWY-2661	3.2.1.141	MONOMER-6022
PWY-2661	3.2.1.68	isoamylase
PWY-2661	5.4.99.15	MONOMER-5701
PWY-2661	5.4.99.15	MONOMER-6021
PWY-2681	2.5.1.112 	adenylate isopentenyltransferase
PWY-2681	2.5.1.112	adenylate isopentenyltransferase
PWY-2681	2.5.1.27	adenylate isopentenyltransferase
PWY-2681	3.2.2	MONOMER-15646
PWY-2701	2.3.1.115	MONOMER-6322
PWY-2701	2.4.1.170	MONOMER-6321
PWY-2721	2.4.1.216	MONOMER-5861
PWY-2722	2.4.1.64	trehalose phosphorylase subunit
PWY-2722	5.4.2.6	MONOMER-6101
PWY-2723	2.4.1.231	trehalose phosphorylase subunit
PWY-2724	1.1.3.20 	AT3G23410-MONOMER
PWY-2724	1.1.3 	long-chain fatty alcohol oxidase
PWY-2724	1.14.14	cytochrome P450 52A3-A
PWY-2781	2.5.1.75	AT2G27760-MONOMER
PWY-2821	1.14.13.124	L-phenylalanine N-monooxygenase
PWY-282	1.2.1	alcohol-forming fatty acyl-CoA reductase
PWY-2821	2.4.1.195	MONOMER-6102
PWY-282	2.3.1.75	AT5G37300-MONOMER
PWY-282	4.1.99.5	fatty aldehyde decarbonylase
PWY-283	6.2.1.25	benzoate-CoA ligase (aerobic)
PWY-283 	6.2.1.25 	benzoate-CoA ligase (anaerobic)
PWY-283	6.2.1.25	MONOMER-926
PWY-2841	1.5.99	cytokinin oxidase
PWY-2861	2.3.1.115	MONOMER-6182
PWY-2861	2.4.1.170	MONOMER-6181
PWY-2901 	2.4.1	cytokinin UDP glycosyltransferase
PWY-2901 	2.4.1	trans-zeatin-O-glucoside UDP glycosyltransferase [multifunctional]
PWY-2902	2.4.1.215	CISZOG1-MONOMER
PWY-2902	2.4.1.215	CISZOG2-MONOMER
PWY-2902	2.4.1 	UDP glucose:cytokinin glycosyltransferase
PWY-2902	2.4.1 	zeatin-O-glucosyltransferase
PWY-2921	1.14.13.119	5-epi-aristolochene-1,3-dihydroxylase
PWY-2921	4.2.3 	5-epi-aristolochene synthase
PWY-2921	4.2.3.61	EAS12-MONOMER
PWY-2921	4.2.3.61	EAS34-MONOMER
PWY-2921	4.2.3.61	EAS37-MONOMER
PWY-2941	1.17.1.8	4-hydroxy-tetrahydrodipicolinate reductase
PWY-2941	2.3.1.89	MONOMER-6602
PWY-2941	3.5.1.47	MONOMER-6604
PWY-2941	4.1.1.20	diaminopimelate decarboxylase subunit
PWY-2941	4.3.3.7	4-hydroxy-tetrahydrodipicolinate synthase
PWY-2941	5.1.1.7	MONOMER-6606
PWY-2941	5.1.1.7	MONOMER-6841
PWY-2942	1.4.1.16	meso-diaminopimelate D-dehydrogenase subunit
PWY-2942	1.4.1.16	meso-diaminopimelate dehydrogenase subunit
PWY-2961	4.2.3.21	VS1-MONOMER
PWY-2981 	1.14.13.144	9&beta;-pimara-7,15-diene oxidase
PWY-2981	4.2.3.28	DTC1-MONOMER
PWY-2981	4.2.3.30 	ent-sandaracopimaradiene synthase
PWY-2981	4.2.3.30	MONOMER-13866
PWY-2981	4.2.3.33	DTC2-MONOMER
PWY-2981	4.2.3.34	MONOMER-13868
PWY-2981	4.2.3.35	MONOMER-6803
PWY-2981	5.5.1.14	CYC1-MONOMER
PWY2OL-4 	2.5.1.1	geranyl diphosphate synthase
PWY2OL-4 	4.2.3.25	(S)-linalool synthase
PWY2OL-4 	4.2.3.26	(R)-linalool synthase
PWY-3001 	2.2.1.6	acetolactate synthase II, large subunit, N-ter fragment (pseudogene) 
PWY-3001 	2.7.1.39 	homoserine kinase
PWY-3001 	4.2.3.1	THRESYN-MONOMER
PWY-301	1.3.8.11	AliB
PWY-301	6.2.1	AliA
PWY-3022 	1.14.13.68	CYP71E7
PWY-3041	4.2.3.105 	myrcene synthase
PWY-3041	4.2.3.106	1,8-cineole synthase
PWY-3041	4.2.3.106	monoterpene synthese [multifunctional]
PWY-3041	4.2.3.108 	1,8-cineole synthase
PWY-3041	4.2.3.113	MONOMER-14948
PWY-3041 	4.2.3.115 	3-carene synthase
PWY-3041	4.2.3 	monoterpene synthase 2
PWY-3061	1.1.1.207 	menthol dehydrogenase
PWY-3061	1.1.1.208 	neomenthol dehydrogenase
PWY-3061	1.14.13.104	MONOMER-6783
PWY-3061	1.14.13.47	(-)-limonene-3-hydroxylase
PWY-3061	1.14.13.47	MONOMER-6761
PWY-3061	1.3.1.81	(+)-pulegone reductase
PWY-3061	1.3.1.82	MONOMER-6684
PWY-3081	1.1.1.87 	homoisocitrate dehydrogenase subunit
PWY-3081	1.2.1	[LysW]-L-2-aminoadipate 6-phosphate reductase
PWY-3081 	1.2.1 	[LysW]-L-2-aminoadipate/[LysW]-L-glutamate phosphate reductase
PWY-3081	2.3.3.14	MONOMER-6722
PWY-3081	2.6.1.39	L-2-aminoadipate aminotransferase
PWY-3081 	2.6.1 	[LysW]-aminoadipate semialdehyde/glutamate semialdehyde transaminase
PWY-3081	2.6.1	[LysW]-L-2-aminoadipate semialdehyde transaminase
PWY-3081 	2.7.2.8	[LysW}-glutamate/[LysW]-aminoadipate kinase
PWY-3081	3.5.1	[L-2-aminoadipate carrier protein]-L-lysine lysine-hydrolase
PWY-3081 	3.5.1	L-2-aminoadipate carrier protein]-lysine/ornithine hydrolase
PWY-3081	4.2.1.36	homoaconitase
PWY-3081	6.3.2.43	L-2-aminoadipate--LysW ligase
PWY-3097 	1.14.13.52	MONOMER-5323
PWY-3097 	1.14.13.89	isoflavone-2-hydroxylase
PWY-3097 	1.14.13.89 	isoflavone-2-hydroxylase / isoflavone 2-hydroxylase
PWY-3097 	1.3.1.51	MONOMER-6703
PWY-3097 	2.1.1.150	MONOMER-5151
PWY-3097 	2.1.1.46	MONOMER-5153
PWY-3097 	4.2.1.105	2-hydroxyisoflavanone dehydratase
PWY-3101	1.14.11.19 	anthocyanidin synthase/flavonol synthase
PWY-3101 	1.14.11.23 	flavonol synthase/flavanone 3-dioxygenase
PWY-3101	1.14.13.21	flavonoid 3-monooxygenase
PWY-31	1.6.6	canaline reductase
PWY-3121	3.2.1.21	linamarase
PWY-3121	3.2.1.21	MONOMER-6861
PWY-3121 	4.1.2.46	(S)-hydroxynitrile lyase
PWY-3121	4.1.2.47 	(S)-hydroxynitrile lyase
PWY-3	1.2.1.54 	aminobutyraldehyde dehydrogenase subunit
PWY-31	3.5.3.1	MONOMER-11543
PWY-3161	1.13.12.3	MONOMER-7661
PWY-3161	1.13.12.3	MONOMER-7863
PWY-3161	3.5.1.4	MONOMER-7681
PWY-3161	3.5.1.4	MONOMER-7862
PWY-3162	1.1.99 	tryptophan side chain oxidase / tryptophan 2-dioxygenase
PWY-3181	1.2.3.1 	aldehyde oxidase / indole-3-acetaldehyde oxidase
PWY-3181 	1.2.3.9 	aldehyde oxidase / indole-3-acetaldehyde oxidase / abscisic aldehyde oxidase
PWY-3181	4.1.1.28	MONOMER-7901
PWY-3181	4.1.1.28	MONOMER-8127
PWY-3181 	4.1.1.28	tryptophan decarboxylase subunit
PWY-321	1.11.2.3	fatty acid peroxygenase (fatty acid hydroperoxide-dependent)
PWY-321	1.14.13.205	long-chain fatty acid &omega;-hydroxylase
PWY-321	1.14.13	long-chain fatty acid &omega;-hydroxylase
PWY-321	3.3.2	fatty acid epoxide hydrolase
PWY-321 	6.2.1.3	AT1G49430-MONOMER
PWY-3261	4.2.1.76	MONOMER-18755
PWY-3	2.6.1	MONOMER-17
PWY-3	2.6.1	polyamine aminotransferase subunit
PWY-3282 	1.4.7.1 	glutamate synthase (ferredoxin-dependent)
PWY-3282 	6.3.1.2	glutamine synthetase
PWY-3282 	6.3.1.2	glutamine synthetase, type II
PWY-3301	2.3.1.103	sinapoylglucose:sinapoylglucose O-sinapoyltransferase
PWY-3301	2.3.1.91	AT5G09640-MONOMER
PWY-3301	2.3.1.91	MONOMER-7342
PWY-3301	2.3.1.91	MONOMER-7921
PWY-3301	2.3.1.91	MONOMER-8001
PWY-3301	2.3.1.92	AT2G22990-MONOMER
PWY-3301	2.3.1.92	MONOMER-7343
PWY-3301	2.4.1.120	AT3G21560-MONOMER
PWY-3301	2.4.1.120	MONOMER-7341
PWY-3301	2.4.1.120	MONOMER-7865
PWY-3301	2.4.1.120	MONOMER-7981
PWY-3301	3.1.1.49	MONOMER-7365
PWY-3341	1.2.1.41 	&delta;1-pyrroline-5-carboxylate synthetase
PWY-3341	1.5.1.2	pyrroline-5-carboxylate reductase
PWY-3341	2.6.1.13	MONOMER-7802
PWY-3461	1.3.1.78	arogenate dehydrogenase monomer
PWY-3481 	1.3.1.78	arogenate dehydrogenase
PWY-3481 	1.3.1.78	AT5G34930-MONOMER
PWY-3481 	2.6.1.79	MONOMER-8101
PWY-3481 	4.2.1.91	MONOMER-8126
PWY-3481 	4.2.1.91	MONOMER-8128
PWY-3481 	4.2.1.91	MONOMER-8129
PWY-3481 	5.4.99.5	AT1G69370-MONOMER
PWY-3481 	5.4.99.5	AT3G29200-MONOMER
PWY-3481 	5.4.99.5	AT5G10870-MONOMER
PWY-3502 	2.4.2.19	YFR047C-MONOMER
PWY-3502 	2.7.7.18 	nicotinate nucleotide adenylyltransferase
PWY-3502 	3.5.1.19	nicotinamidase
PWY-3502 	3.6.1.22	NADH pyrophosphatase 1
PWY-3502 	6.3.4.21	nicotinate phosphoribosyl transferase
PWY-3502 	6.3.5.1	MONOMER3O-845
PWY-3581	1.14.13.71	N-methylcoclaurine 3&prime;-monooxygenase
PWY-3581	2.1.1.116	3-hydroxy-N-methylcoclaurine 4-O-methyltransferase
PWY-3581	2.1.1.128	norcoclaurine 6-O-methyltransferase
PWY-3581	4.2.1.78	MONOMER-13859
PWY-3581	4.2.1.78	norcoclaurine synthase
PWY-3602	1.1.1.108	L-carnitine dehydrogenase
PWY-361	1.1.1.195 	cinnamyl alcohol dehydrogenase 4
PWY-361	1.1.1.195 	cinnamyl alcohol dehydrogenase 5
PWY-361 	1.14.13.36	p-coumaroyl shikimate/quinate 3-hydroxylase
PWY-361	1.2.1.44	cinnamoyl coenzyme A reductase
PWY-361	2.3.1 	hydroxycinnamoyl-CoA: shikimate/quinate hydroxycinnamoyltransferase
PWY-361 	6.2.1.12	4-coumarate coenzyme A ligase 2
PWY-361 	6.2.1.12 	4-coumarate coenzyme A ligase 3
PWY-361 	6.2.1.12	AT3G21230-MONOMER
PWY-3621 	1.1.1.108	L-carnitine dehydrogenase
PWY-3621	1.14.11.1	&gamma;-butyrobetaine dioxygenase
PWY-3641	1.14.13.M4	carnitine monooxygenase
PWY-3641	1.14.13.M4	carnitine monooxygenase complex
PWY-3661-1	1.5.8 	dimethylglycine dehydrogenase
PWY-3661-1	1.5.8 	sarcosine dehydrogenase
PWY-3661	1.5.3.1	sarcosine oxidase &alpha; subunit 
PWY-3661	1.5.8.4	dimethylglycine dehydrogenase, mitochondrial
PWY-3661	2.1.1.5	glycine betaine-homocysteine methyltransferase subunit
PWY-3661	2.1.1.5	MONOMER-8525
PWY-3661	2.1.2.1	serine hydroxymethyltransferase 1
PWY-3721	1.1.3.17	choline oxidase
PWY-3722	1.1.1	MONOMER-8601
PWY-3722	1.2.1.8	betaine aldehyde dehydrogenase subunit
PWY-3781	1.10.2.2	ubiquinol-cytochrome c oxidoreductase
PWY-3781 	1.3.5.1	succinate dehydrogenase complex
PWY-3781	1.9.3.1	cytochrome-c oxidase
PWY-381	1.7.1.1	assimilatory nitrate reductase (NADH)
PWY-3821 	2.7.1.6	AT3G06580-MONOMER
PWY-3821 	2.7.1.6	galactokinase
PWY-3821 	2.7.7.64 	MONOMER-15711
PWY-3841 	1.4.4.2 	glycine decarboxylase
PWY-3841	1.5.1.20	MTHFR2 subunit
PWY-3841	1.5.1.20	MTHFR subunit
PWY-3841	1.5.1.5 	5,10-methylenetetrahydrofolate dehydrogenase [multifunctional]
PWY-3841	6.3.3.2	AT5G13050-MONOMER
PWY-3841	6.3.3.2	MONOMER-9382
PWY-3841 	6.3.4.3	10-formyltetrahydrofolate synthetase
PWY-3841	6.3.4.3	10-formyltetrahydrofolate synthetase
PWY-3861	1.1.1.255	mannitol dehydrogenase
PWY-3861	2.7.1.7	mannokinase
PWY-3861 	5.3.1.8	phosphomannose isomerase
PWY-3861	5.3.1.8	phosphomannose isomerase
PWY-3881	1.1.1.224	mannose 6-phosphate reductase
PWY-3881	3.1.3.22	MONOMER-9341
PWY-3901 	1.14.21.5 	benzylisoquinoline synthase
PWY-3901	1.14.21.5	MONOMER-9451
PWY-3901	1.21.3.3	berberine bridge enzyme
PWY-3901	1.3.3.8	MONOMER-16686
PWY-3901	1.3.3.8	MONOMER-9454
PWY-3901	2.1.1.117	MONOMER-9444
PWY-3901	2.1.1.117	scoulerine 9-O-methyltransferase
PWY-3941	2.6.1.18	MONOMER-9450
PWY-3981	1.2.1	3-aminopropionaldehyde dehydrogenase
PWY-3981	1.4.3.22	diamine oxidase
PWY-3982	1.3.1.2	dihydropyrimidine dehydrogenase (NADP<small><sup>+</sup></small>)
PWY-3982	3.5.1.6	&beta;-ureidopropionase
PWY-3982	3.5.2.2	dihydropyrimidinase
PWY-3982	3.5.2.2	MONOMER-9448
PWY3DJ-11470	1.3.1.27	prostaglandin reductase 1
PWY3DJ-11470	2.7.1.91	sphingosine kinase 1
PWY3DJ-11470	2.7.1.91	sphingosine kinase 2
PWY3DJ-11470	3.5.1.23	N-acylsphingosine amidohydrolase 1
PWY3DJ-11470	3.5.1.23	N-acylsphingosine amidohydrolase 2
PWY3DJ-11470	3.5.1.23	N-acylsphingosine amidohydrolase (alkaline ceramidase) 3
PWY3DJ-11470	4.3.2	sphingosine phosphate lyase 1
PWY3DJ-35471	3.1.1.17	MONOMER-13233
PWY3DJ-35471	3.1.1.17	regucalcin
PWY3DJ-35471	3.2.1	UDP-glucuronidase
PWY3O-1109 	2.6.1.57 	aromatic amino acid aminotransferase I
PWY3O-355 	2.3.1 	fatty acid synthase, &beta; subunit
PWY3O-4106	2.7.1.173 	nicotinamide riboside kinase
PWY3O-4106 	2.7.7.1 	nicotinamide mononucleotide adenylyltransferase 1
PWY3O-4106 	2.7.7.1 	YGR010W-MONOMER
PWY3O-4108 	1.1.1.1 	alcohol dehydrogenase IV
PWY3O-4108 	2.6.1.58 	aromatic amino acid aminotransferase II
PWY3O-450 	2.7.1.32	choline kinase &alpha;1
PWY3O-450	2.7.1.32	YLR133W-MONOMER
PWY3O-450 	2.7.7.15	choline-phosphate cytidylyltransferase A
PWY3O-450	2.7.7.15	YGR202C-MONOMER
PWY3O-450 	2.7.8.2	cholinephosphotransferase 1
PWY3O-6	1.1.1.116	MONOMER3O-32
PWY3O-6	1.1.1.117	D-arabinose dehydrogenase, heavy chain 
PWY3O-6	1.1.3.37	D-arabinono-1,4-lactone oxidase
PWY-4002	2.6.1.14 	MONOMER-9722
PWY-4002	3.5.1	MONOMER-9724
PWY-401	2.4.1.184	galactolipid:galactolipid galactosyltransferase
PWY-401	2.4.1.241	AT3G11670-MONOMER
PWY-401	2.4.1.241	UDP-galactose:MGDG galactosyltransferase
PWY-401	2.4.1.46	UDP-galactose:DAG galactosyltransferase
PWY-4021	2.1.1 	&beta;-alanine N-methyltransferase
PWY-40	3.5.3.11	agmatinase
PWY-40	4.1.1.19	SPEABACSU-MONOMER
PWY-4041	2.3.2.4	MONOMER-10043
PWY-4041	2.3.2.4	MONOMER-10044
PWY-4041	3.4.19.13 	&gamma;-glutamyltransferase 1
PWY-4041 	3.4.19.13 	&gamma;-glutamyltransferase 5
PWY-4041	3.5.2.9	5-oxo-L-prolinase subunit
PWY-4041	6.3.2.2	glutamate--cysteine ligase
PWY-4041	6.3.2.2	glutamate--cysteine ligase catalytic subunit 
PWY-4041	6.3.2.3	Glutathione synthetase
PWY-4041	6.3.2.3	glutathione synthetase large subunit 
PWY-4061	2.5.1.18	glutathione S-transferase A2
PWY-4061	2.5.1.18	glutathione S-transferase A heterodimer
PWY-4061	2.5.1.18 	Maleylacetoacetate isomerase
PWY-4081	1.11.1.12	glutathione peroxidase 4
PWY-4081 	1.11.1.9	glutathione peroxidase 1
PWY-4081 	1.11.1.9	glutathione peroxidase 2
PWY-4081 	1.11.1.9	glutathione peroxidase 3
PWY-4081 	1.8.1.7	glutathione reductase subunit
PWY-4081 	1.8.1.7	MONOMER-5144
PWY-4081 	1.8.1.7	MONOMER-5145
PWY-4101	1.1.1.14 	D-sorbitol dehydrogenase
PWY-4101	1.1.1.14	MONOMER-11719
PWY-4101 	1.1.1.14 	MONOMER-13193
PWY-4101	1.1.1.14	MONOMER-17164
PWY-4101	1.1.1.56 	sorbitol dehydrogenase
PWY-4101	1.1.1.9 	glucitol dehydrogenase
PWY-4101	1.1.1.9 	L-iditol 2-dehydrogenase
PWY-4121	3.5.1.78 	CFGSS-MONOMER
PWY-4121 	3.5.1.78 	CFTRS-MONOMER
PWY-4121	6.3.1.9 	fused glutathionylspermidine amidase / glutathionylspermidine synthetase
PWY-4161 	1.14.13.137	indole monooxygenase
PWY-4161 	1.14.13.138	indolin-2-one monooxygenase
PWY-4161 	1.14.13.139	3-hydroxyindolin-2-one monooxygenase
PWY-4161 	1.14.13.140	HBOA monooxygenase
PWY-4161 	1.14.20.2	DIMBOA-glucoside dioxygenase
PWY-4161 	2.1.1.241	TRIBOA-glucoside methyltransferase
PWY-4161 	2.4.1.202	DIBOA 2-<small>D</small>-glucosyltransferase
PWY-4161 	4.1.2.8	indole-3-glycerol phosphate lyase
PWY-4181	1.11.1.17	MONOMER-10183
PWY-4181	1.8.1.16	glutathione amide reductase subunit
PWY-4201	2.1.1.146	CVOMT1 subunit
PWY-4201	2.1.1.146	EOMT1 subunit
PWY-4201	2.1.1.146	MONOMER-13884
PWY-4201	2.1.1.146 	MONOMER-17506
PWY-4201 	2.1.1.279 	MONOMER-15433
PWY-4201	2.1.1.279 	S-adenosyl-L-methionine:(iso)eugenol-O-methyltransferase
PWY-4202	2.1.1.137	arsenic methyltransferase
PWY-4203	2.1.1.273 	S-adenosyl-L-methionine:benzoic acid-salicylic acid carboxymethyltransferase
PWY-4203 	2.3.1.196	MONOMER-10246
PWY-4203	2.3.1.224	MONOMER-10302
PWY-4221 	2.1.2.11	ketopantoate hydroxymethyltransferase
PWY-4221 	6.3.2.1	pantoate:&beta;-alanine ligase
PWY-4242 	1.1.1.358 	2-dehydropantolactone reductase
PWY-4242 	2.7.1.24	AT2G27490-MONOMER
PWY-4242 	2.7.1.33	AT4G32180-MONOMER
PWY-4242 	2.7.1.33	pantothenate kinase 1
PWY-4242 	2.7.7.3	AT2G18250-MONOMER
PWY-4242 	4.1.1.36	AT3G18030-MONOMER
PWY-4242 	6.3.2.5	AT1G12350-MONOMER
PWY-4261	2.7.1.30	AT1G80460-MONOMER
PWY-4	2.7.1.44	MONOMER-18
PWY-4302	1.10.3.11	ubiquinone oxidase
PWY-4302	1.10.3.11	ubiquinone oxidase 1
PWY-4321	1.1.1 	&gamma;-hydroxybutyrate dehydrogenase
PWY-4321 	1.2.1.24	NAD-dependent succinate semialdehyde dehydrogenase
PWY-4321	1.2.1.24	NAD-dependent succinate semialdehyde dehydrogenase
PWY-4321 	2.6.1.96	&gamma;-aminobutyrate transaminase (pyruvate dependent)
PWY-4321 	2.6.1.96	MONOMER-15561
PWY-4321 	2.6.1.96	MONOMER-15562
PWY-4321	4.1.1.15	AT1G65960-MONOMER
PWY-4321	4.1.1.15	AT5G17330-MONOMER
PWY-4321	4.1.1.15	MONOMER-15560
PWY-43	3.5.1.53	MONOMER-17350
PWY-43	3.5.3.12	MONOMER-17349
PWY-43	4.1.1.19	arginine decarboxylase (biosynthetic)
PWY-4381	2.3.1.180	&beta;-ketoacyl-[acyl-carrier-protein] synthase
PWY-4421	2.4.1	UDPG:glucosyltransferase
PWY-4441	3.2.1.182	&beta;-glucosidase
PWY-4521	1.20.2.1	arsenite oxidoreductase
PWY-4601	1.20.99.1	respiratory arsenate reductase large subunit 
PWY-4621	1.20.4.1	cytoplasmic arsenate reductase subunit
PWY-46	4.1.1.17	ornithine decarboxylase
PWY-4661	2.7.1.159 	myo-inositol polyphosphate kinase
PWY-4661	2.7.1.64	MONOMER-10821
PWY-4661	5.5.1.4	myo-inositol-3-monophosphate synthase
PWY-4681 	2.5.1	DMAPP:isoflavone dimethylallyl transferase
PWY-4702	3.1.3.25 	phytase / 4-phytase
PWY-4702	3.1.3.62 	phytase / 3-phytase
PWY-4702	3.1.3.62 	phytase / 4-phytase
PWY-4722	1.5.8.3	sarcosine dehydrogenase subunit
PWY-4722	3.5.1.59	N-carbamoylsarcosine amidohydrolase subunit
PWY-4722	3.5.2.14	N-methylhydantoin amidohydrolase &alpha; subunit 
PWY-4741	3.5.3.16	MONOMER-11026
PWY-4762 	2.1.1.103	phosphoethanolamine N-methyltransferase
PWY-4762 	2.1.1	ethanolamine N-methyltransferase
PWY-4762 	2.1.1	MONOMER-8301
PWY-4762 	2.7.1.82	MONOMER-7807
PWY-4762 	2.7.1.82	MONOMER-8201
PWY-4762 	2.7.7.15	cholinephosphate cytidylyltransferase subunit
PWY-4762 	2.7.7.15	MONOMER-8323
PWY-4762 	2.7.8.2	MONOMER-8324
PWY-4762 	3.1.3.75	MONOMER-7824
PWY-4762 	3.1.4.4	MONOMER-8325
PWY-4762 	3.1.4.4	PLD &alpha;
PWY-4762 	3.1.4.4	PLD &beta;
PWY-4762 	3.1.4.4	PLD &delta;
PWY-4762 	3.1.4.4	PLD &gamma;
PWY-4762 	3.1.4.4	PLD &zeta;
PWY-4762 	4.1.1	MONOMER-7804
PWY-4762 	4.1.1	MONOMER-8241
PWY-4762 	4.1.1	serine decarboxylase subunit
PWY-4781	3.1.3.62	multiple inositol-polyphosphate phosphatase
PWY-4801	2.3.1	ALS subunit
PWY-481	1.1.1.311	1-phenylethanol dehydrogenase subunit
PWY-481	1.17.99.2	&alpha; subunit of ethylbenzene dehydrogenase 
PWY-481	1.17.99.2	EBDA subunit 
PWY-481	6.4.1.8	acetophenone carboxylase subunit 
PWY-4821 	1.1.1.22	UDP-glucose dehydrogenase
PWY-4821	4.1.1.35	UDP-xylose synthase 1
PWY-4841	1.13.99.1	AT1G14520-MONOMER
PWY-4841	1.13.99.1	AT2G19800-MONOMER
PWY-4841	1.13.99.1	AT4G26260-MONOMER
PWY-4841	2.7.1.43	glucuronokinase
PWY-4841	2.7.1.43	MONOMER-11167
PWY-4841	2.7.1.43	MONOMER-11206
PWY-4841	2.7.7.44	MONOMER-11216
PWY490-3	1.7.7.1	SYNPCC7942_1240-MONOMER
PWY490-3	1.7.7.2	MONOMER-13441
PWY490-4	6.1.1 	aspartyl-tRNA synthetase, non-discriminating
PWY-4921	3.5.3.15	MONOMER-11208
PWY-4921	3.5.3.15	peptidylarginine deiminase I
PWY-4921	3.5.3.15	peptidylarginine deiminase II
PWY-4921	3.5.3.15	peptidylarginine deiminase III
PWY-4921	3.5.3.15	peptidylarginine deiminase IV
PWY-4921	3.5.3.15	peptidylarginine deiminase type VI
PWY-4922	2.1.1.108	MONOMER-11169
PWY-4922	2.3.1	6HM subunit
PWY-4942	2.1.1.79	MONOMER-11215
PWY-4961	2.5.1.51	&beta;-pyrazole-1-ylalanine synthase
PWY-4981	1.5.1.2	MONOMER-11347
PWY-4981	2.6.1.13	MONOMER-11346
PWY-4981	3.5.1.20	MONOMER-11345
PWY-4981	3.5.3.6	MONOMER-11344
PWY-4985	2.5.1.52	mimosine synthase
PWY4FS-11	1.3.2	MONOMER4FS-17
PWY4FS-11 	5.1.3 	GDP-L-gulose-3-epimerase [multifunctional]
PWY4FS-12 	2.7.7.69	GDP-L-galactose phosphorylase
PWY4FS-5 	2.1.1.103	phosphoethanolamine N-methyltransferase
PWY4FS-5 	2.1.1.71	MONOMER4FS-5
PWY4FS-5 	2.1.1.71	MONOMER4FS-6
PWY4FS-5 	2.1.1.71	MONOMER4FS-7
PWY4FS-5 	2.1.1.71	MONOMER4FS-8
PWY4FS-5 	2.7.1.32	MONOMER-12865
PWY4FS-5 	2.7.1.32	MONOMER-12866
PWY4FS-5 	2.7.7.57	MONOMER4FS-3
PWY4FS-5 	2.7.7.57	MONOMER4FS-4
PWY4FS-5 	2.7.7	MONOMER4FS-9
PWY4FS-6 	2.7.8.1 	aminoalcoholphosphotransferase
PWY4FS-7 	2.3.1.15	chloroplastidic glycerol-3-phosphate O-acyltransferase
PWY4FS-7 	2.3.1.51	AT4G30580-MONOMER
PWY4FS-8 	2.3.1.15	AT4G00400-MONOMER
PWY4FS-8 	2.3.1.15	cytoplasmic glycerol-3-phosphate acyltransferase
PWY4FS-8 	2.3.1.51	MONOMER-12148
PWY4FS-8 	2.3.1.51	MONOMER-12149
PWY4FS-8 	2.3.1.51	MONOMER-12150
PWY4FS-8 	2.7.7.41	AT1G62430-MONOMER
PWY4FS-8 	2.7.7.41 	AT5G04490-MONOMER
PWY4FS-8 	2.7.7.41	MONOMER-12282
PWY4LZ-257	1.1.1.28	D-lactate dehydrogenase
PWY4LZ-257 	1.2.1.10 	alcohol dehydrogenase / acetaldehyde dehydrogenase
PWY4LZ-257 	1.2.7.1	pyruvate:ferredoxin oxidoreductase
PWY4LZ-257 	2.3.1.54	pyruvate formate lyase
PWY4LZ-257 	2.3.1.8	phosphate acetyltransferase
PWY4LZ-257 	2.7.2.1	acetate kinase
PWY4LZ-257 	4.1.1.1	pyruvate decarboxylase
PWY-5004 	1.14.13.39	eNOS
PWY-5004 	1.14.13.39	iNOS
PWY-5004 	1.14.13.39	nNOS
PWY-5004 	1.2.1.41 	&delta;-1-pyrroline-5-carboxylate synthase
PWY-5004 	1.2.1.41 	&delta;1-pyrroline-5-carboxylate synthetase
PWY-5004 	1.4.7.1 	glutamate synthase (ferredoxin-dependent)
PWY-5004 	1.5.5.2	mitochondrial proline dehydrogenase 1
PWY-5004 	1.5.5.2	proline dehydrogenase 1
PWY-5004 	1.5.5.2	proline dehydrogenase 2
PWY-5004 	2.1.3.3	Ornithine carbamoyltransferase, mitochondrial
PWY-5004 	2.6.1.13	Ornithine aminotransferase
PWY-5004 	2.6.1.13	ornithine-&delta;-aminotransferase
PWY-5004 	3.5.1.38 	glutaminase, kidney isoform
PWY-5004 	3.5.3.1	arginase 1
PWY-5004 	3.5.3.1	arginase 2
PWY-5004 	3.5.3.1	AT4G08900-MONOMER
PWY-5004 	4.3.2.1	argininosuccinate lyase
PWY-5004 	4.3.2.M2 	imidazole-glycerol-phosphate synthase
PWY-5004 	4.3.3.6 	glutaminase 
PWY-5004 	6.3.4.16 	ammonia-dependent carbamoyl-phosphate synthase
PWY-5004 	6.3.4.2 	CTP synthase 1
PWY-5004 	6.3.4.2 	CTP synthase 2
PWY-5004 	6.3.4.5	Argininosuccinate synthase
PWY-5005 	2.3.1.47	8-amino-7-oxononanoate synthase
PWY-5005 	2.6.1.105	lysine--8-amino-7-oxononanoate transaminase
PWY-5005 	2.8.1.6	BIOBBACSU-MONOMER
PWY-5005 	6.2.1.14	6-carboxyhexanoate-CoA ligase
PWY-5005 	6.2.1.14	6-carboxyhexanoate--CoA ligase monomer
PWY-5005 	6.3.3.3	dethiobiotin synthetase subunit
PWY-5021	2.5.1.53	MONOMER-11507
PWY-5022	1.1.1.61	MONOMER-13486
PWY-5022	2.8.3	4-hydroxybutyrate CoA-transferase subunit
PWY-5022	5.3.3.3 	4-hydroxybutyryl-CoA dehydratase subunit
PWY-5024	2.6.1	arginine-&alpha;-ketoglutarate transaminase subunit
PWY-5024 	3.5.3.7	guanidinobutyrase subunit
PWY-5025	4.2.1.84	nitrile hydratase subunit
PWY-5026	3.5.5.1	nitrilase subunit
PWY-5028	3.5.1.68	MONOMER-11615
PWY-5028	3.5.2.7	MONOMER-11617
PWY-5028	3.5.3.13	formiminoglutamate deiminase subunit
PWY-5028	4.2.1.49	urocanase subunit
PWY-5028	4.3.1.3	histidase subunit
PWY-5029	1.1.1.111	MONOMER-11648
PWY-5029	2.6.1.38	MONOMER-11647
PWY-5030 	1.5.1.5 	C1-THF synthase
PWY-5030	4.3.1.4 	formimidoyltransferase-cyclodeaminase dimer
PWY-5031	2.6.1.58	histidine-pyruvate aminotransferase subunit
PWY-5035 	1.14.11.12 	gibberellin 20-oxidase
PWY-5035 	1.14.11 	gibberellin 20-oxidase
PWY-5035 	1.14.11	gibberellin 20-oxidase
PWY-5035	1.14.11	MONOMER-17926
PWY-5035	1.14.11	MONOMER-17927
PWY-5036 	1.14.11	gibberellin 20-oxidase
PWY-5037	2.1.1.159	theobromine synthase
PWY-5038 	2.1.1.158	7-methylxanthosine synthase
PWY-5038 	2.1.1.160 	caffeine synthase
PWY-5038 	3.2.2.25	N-methyl nucleosidase
PWY-5041 	2.1.1.13 	MONOMER-9225
PWY-5041	2.1.1.14	MONOMER-11750
PWY-5041	2.5.1.6	MONOMER-11743
PWY-5041	3.3.1.1	S-adenosylhomocysteine hydrolase
PWY-5045	2.3.1.146	pinosylvin synthase
PWY-5046	1.2.1 	branched-chain &alpha;-keto acid dehydrogenase complex E2 component
PWY-5046	1.2.1 	branched-chain &alpha;-keto acid dehydrogenase E1 component
PWY-5047	1.14.11	gibberellin 13-hydroxylase
PWY-5047	1.14.11 	gibberellin 2&alpha;-hydroxylase / gibberellin 1,2-desaturase
PWY-5047	1.14.13	gibberellin 20-oxidase
PWY-5049 	1.1.1.237 	hydroxyphenylpyruvate reductase
PWY-5049	2.3.1.140	MONOMER-16537
PWY-5049 	2.3.1.140 	rosmarinic acid synthase
PWY-5052 	1.14.11.15	AT1G80330-MONOMER
PWY-5052 	1.14.11.15 	gibberellin 3-oxidase
PWY-5052 	1.14.11 	gibberellin 20-oxidase
PWY-5052 	1.14.11	gibberellin 20-oxidase
PWY-5052 	1.14.11	MONOMER-11682
PWY-5053 	1.14.11	gibberellin 3&beta;-hydroxylase / ent-kaurenoate oxidase
PWY-5053 	4.2.3.19	ent-kaurene synthase
PWY-5053 	4.2.3.19	ent-kaurene synthase B
PWY-5053 	5.5.1.13	CYC2-MONOMER
PWY-5053 	5.5.1.13	ent-copalyl diphosphate synthase
PWY-5053 	5.5.1.13	MONOMER-7941
PWY-5054	1.1.1.200	sorbitol-6-phosphate 1-dehydrogenase
PWY-5054	3.1.3.50	MONOMER-11698
PWY-5059	2.3.1	PDCHSX subunit
PWY-5059	2.3.1	PSCHS subunit
PWY-5059	2.3.1	SbCHS2 subunit
PWY-5059	5.5.1.6	chalcone isomerase I
PWY-5059	5.5.1.6	chalcone isomerase II
PWY-5060	1.14.11.22	flavone synthase I
PWY-5060 	1.14.13.21	flavonoid 3-hydroxylase
PWY-5060	1.14.13.21	flavonoid 3-hydroxylase
PWY-5061 	1.14.13.86 	isoflavone synthase
PWY-5061 	1.14.13 	isoflavone synthase
PWY-5061 	4.2.1 	2-hydroxyisoflavanone dehydratase
PWY-5062 	1.1.1.291	2-(hydroxymethyl)glutarate dehydrogenase subunit
PWY-5062 	1.17.1.5	nicotinate hydroxylase 50 kD subunit 
PWY-5062 	1.17.3.3	MONOMER-11655
PWY-5062 	1.17.99	MONOMER-11653
PWY-5062 	1.3.7.1	6-hydroxynicotinate reductase subunit
PWY-5062 	3.5.2.18	enamidase subunit
PWY-5062 	4.1.3.32	2,3-dimethylmalate lyase subunit
PWY-5062 	4.2.1.85	dimethylmaleate hydratase subunit A 
PWY-5062 	5.3.3.6	3-methylitaconate &Delta;-isomerase
PWY-5062 	5.4.99.4	&alpha;-methyleneglutarate mutase subunit
PWY-5063	1.3.1.83	geranylgeranyl diphosphate reductase
PWY-5064 	2.5.1.62	chlorophyll synthetase
PWY-5067	2.4.1.11	glycogen synthase
PWY-5067	2.4.1.18	1,4-&alpha;-glucan branching enzyme
PWY-5067	2.4.1.186	glycogenin 1
PWY-5067	2.4.1.186	glycogenin 2
PWY-5068	1.1.1.294	chlorophyll b reductase
PWY-5068	1.1.1.294	chlorphyll b reductase
PWY-5068	1.1.1.294	MONOMER-13499
PWY-5068	1.17.7.2	7<sup>1</sup>-hydroxychlorophyll a reductase
PWY-5068	1.17.7.2	MONOMER-13500
PWY-5070 	1.14.11.15 	gibberellin 3-oxidase (multifunctional)
PWY-5070 	1.14.11	gibberellin 2&beta;,3&beta;-hydroxylase
PWY-5070 	1.14.11 	gibberellin 3-oxidase
PWY-5071 	1.14.13	4-coumaroyl-4-hydroxyphenyllactate 3-hydroxylase
PWY-5071 	1.14.13	hydroxycinnamoyl-4-hydroxyphenyllactate 3-hydroxylase
PWY-5071 	1.14.13	MONOMER-11772
PWY-5071 	1.14.16.2	tyrosine 3-monooxygenase
PWY-5071 	2.6.1.57 	subunits of TAT1
PWY-5071 	2.6.1.57 	subunit TAT 2
PWY-5074	1.1.1.88	hydroxymethylglutaryl-CoA reductase
PWY-5075	5.4.3.7	MONOMER-11830
PWY-5076 	1.1.1 	methylglyoxal reductase
PWY-5078 	2.6.1.42 	branched-chain amino acid aminotransferase
PWY-5078 	4.1.1.72 	keto-isocaproate decarboxylase subunit
PWY-5080 	1.1.1.330	very-long-chain 3-oxoacyl-CoA reductase
PWY-5080 	1.3.1.38 	2-enoyl thioester reductase
PWY-5080 	1.3.1.38 	peroxisomal trans-2-enoyl-CoA reductase
PWY-5080 	4.2.1.119 	peroxisomal multifunctional enzyme type 2
PWY-5082 	1.1.1.1 	alcohol dehydrogenase I
PWY-5082 	1.1.1.1 	alcohol dehydrogenase II
PWY-5082 	1.1.1.1 	alcohol dehydrogenase III
PWY-5082 	1.1.1.1 	alcohol dehydrogenase V
PWY-5082 	4.1.1.74 	2-oxo acid decarboxylase
PWY-5083	1.6.1.2 	MONOMER-11876
PWY-5083	1.6.1.2 	MONOMER-11877
PWY-5083 	1.6.5.3	NADH-ubiquinone oxidoreductase
PWY-5083	2.7.1.23	AT1G21640-MONOMER
PWY-5083	2.7.1.86 	NADK1 subunit
PWY-5083	2.7.1.86	NADK3 subunit
PWY-5083	3.1.3.2 	MONOMER-11869
PWY-5083	3.1.3.2 	MONOMER-11878
PWY-5083	3.1.3	MONOMER-11870
PWY-5084	1.2.1 	dihydrolipoyl dehydrogenase component
PWY-5084 	1.2.4.2 	2-oxoglutarate dehydrogenase complex, (E2) inner core
PWY-5084	1.2.4.2 	dihydrolipoyltranssuccinylase
PWY-5084 	1.8.1.4 	2-ketoglutarate dehydrogenase complex
PWY-5084 	1.8.1.4 	2-oxoglutarate dehydrogenase complex
PWY-5084	2.3.1.61 	2-oxoglutarate decarboxylase, thiamine-requiring
PWY-5084 	2.3.1.61 	2-oxoglutarate dehydrogenase complex E1 dimer
PWY-5086 	1.3.1.75	3,8-divinyl protochlorophyllide a 8-vinyl-reductase (NADPH)
PWY-5086	2.5.1.62	MONOMER-11751
PWY-5096	2.6.1.2	alanine aminotransferase subunit
PWY-5096 	4.1.1.1 	pyruvate ferredoxin oxidoreductase &alpha; subunit 
PWY-5096 	6.2.1.13	acetate-CoA ligase (ADP-forming)
PWY-5097	1.17.1.8	4-hydroxy-tetrahydrodipicolinate reductase
PWY-5097	2.6.1.83	MONOMER-15639
PWY-5097	4.3.3.7	4-hydroxy-tetrahydrodipicolinate synthase
PWY-5097	5.1.1.7	diaminopimelate epimerase
PWY-5098 	1.3.7.12	AT4G37000-MONOMER
PWY-5098	3.1.1.14	AT1G19670-MONOMER
PWY-5098	3.1.1.14	AT5G43860-MONOMER
PWY-5100 	1.1.1.27	Ldh
PWY-5100	1.2.7.1	pyruvate:ferredoxin oxidoreductase &alpha; subunit 
PWY-5100 	2.3.1.8	Pta
PWY-5100 	2.7.2.1	AckA
PWY-5101 	1.1.1.383 	ketol-acid reductoisomerase
PWY-5101 	1.1.1	3-isopropylmalate dehydrogenase subunit
PWY-5101 	2.2.1.6	acetohydroxy-acid synthase
PWY-5101	2.3.1.182	(R)-citramalate synthase subunit
PWY-5101	2.3.1.182	(R)-citratemalate synthase subunit
PWY-5103	5.4.99.1	glutamate mutase &sigma; subunit 
PWY-5104	1.2.7	&alpha;-ketobutyrate synthase &alpha; subunit 
PWY-5104 	2.6.1.42	branched-chain-amino-acid aminotransferase subunit
PWY-5104 	4.2.1.9	dihydroxy-acid dehydratase subunit
PWY-5105 	2.4.1.185	flavanone 7-O-glucosyltransferase
PWY-5105 	2.4.1.236	1,2 rhamnosyltransferase
PWY-5105 	2.4.1	flavanone 7-O-glucoside-6\-O-rhamnosyltransferase"
PWY-5107	2.7.4	MONOMER-11913
PWY-5109	1.1.1.178	2-methylacetoacetyl-coenzyme A reductase subunit
PWY-5109	1.3.8.5	2-methyl branched-chain acyl-CoA dehydrogenase
PWY-5110	2.1.1.7	MONOMER-11926
PWY-5110	2.1.1.7	MONOMER-11927
PWY-5113 	4.1.1 	UDP-D-apiose/UDP-D-xylose synthase subunit
PWY-5114 	1.1.1.22	MONOMER-9001
PWY-5114 	4.1.1.35	AT3G53520-MONOMER
PWY-5114 	4.1.1.35	AT5G59290-MONOMER
PWY-5114 	4.1.1.35	MONOMER-11101
PWY-5114 	4.1.1.35	subunit of UDP-glucuronate decarboxylase
PWY-5114 	4.1.1.35 	UDP-D-glucuronate decarboxylase subunit
PWY-5114 	4.2.1.76	AT1G78570-MONOMER
PWY-5114 	5.1.3.2	UDP-galactose 4-epimerase
PWY-5114 	5.1.3.2	UDP-glucose 4-epimerase
PWY-5114 	5.1.3.5	AT1G30620-MONOMER
PWY-5114 	5.1.3.5	MONOMER-74
PWY-5114 	5.1.3.6	subunit of UDP-GlcA 4-epimerase
PWY-5116	2.1.1.232	MONOMER-11973
PWY-5118	2.1.1.231	MONOMER-11996
PWY-5119	2.1.1.75	MONOMER-11997
PWY-5120	2.5.1.29	farnesyltranstransferase
PWY-5120 	2.5.1.29 	short chain isoprenyl diphosphate synthase
PWY-5121 	1.1.1.267	1-deoxy-D-xylulose 5-phosphate reductoisomerase
PWY-5121 	1.1.1.267	AT5G62790-MONOMER
PWY-5121 	1.17.7.3	EG10370-MONOMER
PWY-5121 	1.17.7.4	1-hydroxy-2-methyl-2-(E)-butenyl 4-diphosphate reductase
PWY-5121 	2.2.1.7	AT4G15560-MONOMER
PWY-5121 	2.5.1.10	AT4G17190-MONOMER
PWY-5121 	2.5.1.10 	AT5G47770-MONOMER
PWY-5121 	2.5.1.1	AT2G34630-MONOMER
PWY-5121 	2.5.1.1 	farnesyl diphosphate synthase
PWY-5121 	2.7.1.148	4-diphosphocytidyl-2-C-methylerythritol kinase
PWY-5121 	2.7.7.60	4-diphosphocytidyl-2C-methyl-D-erythritol synthase
PWY-5121 	2.7.7.60	AT2G02500-MONOMER
PWY-5121 	4.6.1.12	2-C-methyl-D-erythritol 2,4-cyclodiphosphate synthase
PWY-5121 	4.6.1.12	AT1G63970-MONOMER
PWY-5121 	5.3.3.2	AT3G02780-MONOMER
PWY-5121 	5.3.3.2	AT5G16440-MONOMER
PWY-5121 	5.3.3.2	isopentenyl diphosphate isomerase
PWY-5121	5.3.3.2	MONOMER-12170
PWY-5122 	2.5.1.10 	farnesyl diphosphate synthase
PWY-5122 	2.5.1.1 	FPPSYN-MONOMER
PWY-5122	2.5.1.1	geranyl pyrophosphate synthase, large subunit 
PWY-5123	2.5.1.10	MONOMER-11954
PWY-5125 	1.14.11.19	anthocyanidin synthase
PWY-5125 	2.4.1.115	anthocyanidin 3-O-glucosyltransferase
PWY-5129	1.1.1.102 	3-ketodihydrosphinganine reductase
PWY-5129	1.1.1.102	MONOMER-11981
PWY-5129	1.14.18.5	AT1G14290-MONOMER
PWY-5129	1.14.18.5	AT1G69640-MONOMER
PWY-5129	1.14.18.5	MONOMER-11984
PWY-5129	1.14.18.5	sphingolipid C4-hydroxylase
PWY-5129	1.14.18.7	AT2G34770-MONOMER
PWY-5129	1.14.19.17	sphingolipid 4-desaturase
PWY-5129	1.14.19.29	sphingolipid 8-desaturase
PWY-5129	2.3.1.50	MONOMER-11980
PWY-5129	2.3.1.50	serine palmitoyltransferase subunit LCB1 
PWY-5129	2.3.1	ceramide synthase
PWY-5129	2.3.1	MONOMER-15535
PWY-5129	2.4.1.80	MONOMER-11987
PWY-5129	2.4.1	MONOMER-11988
PWY-5129	2.7.1.138	AT5G51290-MONOMER
PWY-5129	2.7.1	AT2G37940-MONOMER
PWY-5129	2.7.1	MONOMER-11989
PWY-5133 	2.3.1 	phlorisobutyrophenone synthase / phlorisovalerophenone synthase
PWY-5134 	2.3.1 	naringenin chalcone synthase
PWY-5134 	2.5.1	MONOMER-12001
PWY-5134 	2.5.1	MONOMER-12002
PWY-5134 	2.5.1	MONOMER-12004
PWY-5134 	2.5.1	phlorisobutyrophenone dimethylallyltransferase
PWY-5136 	5.1.2.3 	multifunctional protein MFP I
PWY-5137	5.3.3.8	3,2-trans-enoyl-CoA isomerase 2
PWY-5137 	5.3.3.8	delta-3, delta2-enoyl-CoA isomerase
PWY-5137 	5.3.3.8 	MFP a
PWY-5137 	5.3.3.8 	multifunctional protein MFP II
PWY-5137 	5.3.3.8 	multifunctional protein MFP III
PWY-5138	1.3.1.34	MONOMER-12061
PWY-5138	4.2.1.119 	very-long-chain 3-hydroxyacyl-CoA dehydratase
PWY-5140	1.21.3.7	&Delta;<sup>9</sup>-tetrahydrocannabinolate synthase
PWY-5140	1.3.3	cannabichromenate synthase
PWY-5140	1.3 	cannabidiolate synthase
PWY-5140	2.5.1 	MONOMER-12031
PWY-5140	4.4.1.26	olivetolic acid cyclase
PWY-5140	6.2.1	hexanoyl-CoA synthetase
PWY-5142	3.1.2.14	AT3G25110-MONOMER
PWY-5151	1.2.3.13	MONOMER-12071
PWY-5151	2.6.1.57 	MONOMER-12070
PWY-5152 	1.14.13.88 	flavonoid 3,5-hydroxylase
PWY-5153	1.14.11	MONOMER-12083
PWY-5153	2.4.1.115	MONOMER-12085
PWY-5154	2.1.3.9	acetylornithine transcarbamylase subunit
PWY-5154	3.5.1.16	acetylcitrulline deacetylase subunit
PWY-5155 	4.1.1.11 	L-tyrosine decarboxylase subunit
PWY-5156 	1.14.19.2	stearoyl-[acp] 9-desaturase
PWY-5156 	2.3.1.180	&beta;-ketoacyl-[acp] synthase III
PWY-5156 	2.3.1.41	3-ketoacyl-[acp] synthase 2, chloroplastic
PWY-5156 	2.3.1.86 	MONOMER-10721
PWY-5156 	2.3.1.86 	MONOMER-16561
PWY-5156 	6.2.1.3 	long-chain acyl-CoA synthetase 6
PWY-5156 	6.2.1.3 	long-chain acyl-CoA synthetase 7
PWY-5156 	6.4.1.2	biotin carboxyl carrier protein 
PWY-5159	1.2.1.26	&alpha;-ketoglutaric semialdehyde dehydrogenase (4-hydroxy-L-proline-induced isozyme)
PWY-5159 	1.2.1 	&alpha;-ketoglutaric semialdehyde dehydrogenase subunit
PWY-5159 	1.2.1 	MONOMER-18709
PWY-5159	1.4.99	D-hydroxyproline dehydrogenase
PWY-5159	1.4.99	D-hydroxyproline dehydrogenase &gamma; subunit 
PWY-5159	3.5.4.22	&Delta;<sup>1</sup>-pyrroline-4-hydroxy-2-carboxylate deaminase
PWY-5159	5.1.1.8	4-hydroxyproline epimerase
PWY-5159	5.1.1.8	hydroxyproline 2-epimerase subunit
PWY-5161	1.14.13	MONOMER-12118
PWY-5161	2.4.1	6-deoxychalcone 4-glucosyltransferase
PWY-5162	4.2.1.80	MONOMER-17215
PWY-5168 	1.14.13	ferulate 5-hydroxylase
PWY-5168	1.2.1 	hydroxycinnamaldehyde dehydrogenase
PWY-5169	3.5.2.15	cyanurate amidohydrolase
PWY-5172	2.3.3.8	ATP-citrate synthase
PWY-5173 	1.2.1 	pyruvate dehydrogenase complex
PWY-5173 	2.3.3.8	ATP-citrate lyase 
PWY-5174	5.3.99.8	capsanthin-capsorubin synthase
PWY-5175 	5.5.1.18	lycopene &epsilon;-cyclase
PWY-5176	1.14.13.14	MONOMER-12159
PWY-5176	2.4.1.114	MONOMER-12160
PWY-5176	3.2.1.21	coumarinic acid glucoside &beta;-glucosidase
PWY-5176	5.2.1	MONOMER-12161
PWY-5177	1.3.8.6 	glutaryl-CoA dehydrogenase
PWY-5177 	4.2.1.150 	enoyl-CoA hydratase
PWY-5178 	1.14.13.dw 	toluene methyl-monooxygenase
PWY-5181 	1.14.13.2	4-hydroxybenzoate hydroxylase
PWY-5181 	1.17.99.1	&alpha; subunit of 4-cresol dehydrogenase
PWY-5183 	1.13.11.2	MONOMER-11102
PWY-5183	1.13.11.2	MONOMER-11210
PWY-5183 	1.13.11.2	MONOMER-11354
PWY-5183 	1.14.12.10	benzoate 1,2-dioxygenase hydroxylase component 
PWY-5183 	1.14.12.11	toluene dioxygenase iron-sulfur protein 
PWY-5183 	1.14.13.M1	toluene-4-monooxygenase
PWY-5183 	1.14.13	toluene 2-monooxygenase
PWY-5183 	1.14.13 	toluene 3-monooxygenase
PWY-5183 	1.2.1.10	MONOMER-11388
PWY-5183 	1.2.1.10	MONOMER-354
PWY-5183 	1.3.1.19	MONOMER-11352
PWY-5183 	1.3.1.25	MONOMER-2964
PWY-5183 	3.7.1	2-hydroxy-6-oxohepta-2,4-dienoate hydrolase
PWY-5183 	3.7.1	MONOMER-11103
PWY-5183 	3.7.1	MONOMER-11224
PWY-5183 	4.1.3.39	MONOMER-11226
PWY-5183 	4.1.3.39	MONOMER-11384
PWY-5183 	4.1.3.39	MONOMER-353
PWY-5183 	4.2.1.80	2-oxopent-4-enoate hydratase
PWY-5183 	4.2.1.80	MONOMER-11225
PWY-5183 	4.2.1.80	MONOMER-352
PWY-5184 	1.1.1.35	3-hydroxyacyl-CoA dehydrogenase
PWY-5184 	1.3.8.3	benzylsuccinyl-CoA dehydrogenase subunit
PWY-5184 	2.3.1.M17	benzoylsuccinyl-CoA thiolase
PWY-5184 	2.8.3.15	succinyl-CoA:benzylsuccinate CoA-transferase BbsE subunit 
PWY-5184 	4.1.99.11	benzylsuccinate synthase &alpha; subunit 
PWY-5184 	4.2.1.17	putative phenylitaconyl-CoA hydratase subunit
PWY-5194	2.1.1.107 	siroheme synthase
PWY-5194	2.1.1.107	uroporphyrinogen-III C-methyltransferase
PWY-5194	4.99.1.4	sirohydrochlorin ferrochelatase
PWY-5195	1.3.1.92	artemisinic aldehyde &Delta;<sup>11(13)</sup>-reductase
PWY-5195	4.2.3.24	MONOMER-12182
PWY-5195	4.2.3.24	MONOMER-12183
PWY-5197	1.1.98.5	F420-dependent methylglyoxal reductase
PWY-5197	1.2.1.22	lactaldehyde dehydrogenase subunit
PWY-5197	4.1.2.17	fuculose-1-phosphate aldolase subunit
PWY-5198	1.5.98.1	Mtd
PWY-5198	2.5.1.77	7,8-didemethyl-8-hydroxy-5-deazariboflavin synthase subunit 1 
PWY-5198	2.5.1.77	MONOMER-12180
PWY-5198	2.7.7.68	2-phospho-L-lactate guanylyltransferase
PWY-5198	2.7.8.28	LPPG:Fo 2-phospho-L-lactate transferase subunit
PWY-5199	6.3.2.34 	F420:glutamyl ligase
PWY-5203	1.14.99.43	&beta;-amyrin 24-hydroxylase
PWY-5203	2.4.1.262	MONOMER-12198
PWY-5203	2.4.1.272	soyasapogenol B glucuronide galactosyltransferase
PWY-5203	2.4.1.273	soyasaponin III rhamnosyltransferase
PWY-5203 	5.4.99.39	MONOMER-12194
PWY-5203 	5.4.99.39	MONOMER-12195
PWY-5203 	5.4.99.39	MONOMER-12196
PWY-5207	1.12.98.3	F420 non-reducing hydrogenase I
PWY-5207	1.5.98	F420H2 dehydrogenase subunit FpoD 
PWY-5207	1.8.98.1	CoB-CoM heterodisulfide reductase
PWY-5207	1.8.98.1	CPLX-502
PWY-5207	1.8.98.1	heterodisulfide oxidoreductase
PWY-5209	1.5.98	Mer
PWY-5209	1.5.99 	Mtd
PWY-5	2.1.3	ornithine carbamoyltransferase
PWY-5247 	2.1.1.247	methylated [methylamine-specific corrinoid protein]:coenzyme M methyltransferase
PWY-5254	2.5.1.131	(4-{4-[2-(&gamma;-L-glutamylamino)ethyl]phenoxymethyl}furan-2-yl)methanamine synthase
PWY-5254	2.6.1.108	MONOMER-18938
PWY-5254	2.7.4.31	[5-(aminomethyl)furan-3-yl]methyl phosphate kinase
PWY-5254	4.2.3.153	MONOMER-18937
PWY-5254	6.3.4.24	tyramine--L-glutamate ligase
PWY-5257 	1.1.1.10	L-xylulose reductase subunit
PWY-5257 	1.1.1.10	MONOMER-18763
PWY-5257 	1.1.1.11	MONOMER-12213
PWY-5257 	1.1.1.11	MONOMER-12215
PWY-5257 	1.1.1.12	MONOMER-13195
PWY-5257 	1.1.1.179	D-xylose dehydrogenase subunit
PWY-5257 	1.1.1.21	MONOMER-13190
PWY-5257 	1.1.1.21	MONOMER-13197
PWY-5257 	1.1.1.376	MONOMER-13199
PWY-5257 	1.1.1.56	MONOMER-12235
PWY-5257 	1.1.1.56	ribitol dehydrogenase subunit
PWY-5257 	1.1.1.67 	MONOMER-12217
PWY-5257 	1.1.1.9	MONOMER-12199
PWY-5257 	1.1.1.9	MONOMER-12201
PWY-5257 	1.1.1.9	MONOMER-12203
PWY-5257 	1.1.1.9	MONOMER-13192
PWY-5257 	1.1.1	D-arabinose 1-dehydrogenase subunit
PWY-5257 	1.2.1.26	&alpha;-ketoglutarate semialdehyde dehydrogenase
PWY-5257 	1.2.1.26	&alpha;-ketoglutarate semialdehyde dehydrogenase subunit
PWY-5257 	2.3.3.9	apparent malate synthase
PWY-5257 	2.7.1.16	L-ribulokinase
PWY-5257 	2.7.1.17	D-xylulose kinase subunit
PWY-5257 	2.7.1.17	D-xyulokinase subunit
PWY-5257 	2.7.1.17	MONOMER-12200
PWY-5257 	2.7.1.17	MONOMER-12202
PWY-5257 	2.7.1.17	MONOMER-12204
PWY-5257 	2.7.1.17	MONOMER-12216
PWY-5257 	2.7.1.17	MONOMER-12218
PWY-5257 	2.7.1.17 	xylulokinase
PWY-5257 	2.7.1.47	DarK
PWY-5257 	2.7.1.47	D-ribulokinase subunit
PWY-5257 	2.7.1.47	MONOMER-12236
PWY-5257 	2.7.1.53	L-xylulose kinase
PWY-5257 	3.1.1.15	MONOMER-13200
PWY-5257 	4.1.2.28 	CP4-6 prophage; probable 2-keto-3-deoxygluconate aldolase
PWY-5257 	4.1.2.28	KpLE2 phage-like element; 2-dehydro-3-deoxy-D-pentonate aldolase
PWY-5257 	4.2.1.141	2-keto-3-deoxy-D-arabinonate dehydratase subunit
PWY-5257 	4.2.1.141	MONOMER-16376
PWY-5257 	4.2.1.25	L-arabonate dehydratase subunit
PWY-5257 	4.2.1.43	L-2-keto-3-deoxyarabonate dehydratase subunit
PWY-5257 	4.2.1.5	D-arabinonate dehydratase subunit
PWY-5257 	4.2.1.82	CP4-6 prophage; D-xylonate dehydratase
PWY-5257 	4.2.1.82	D-xylonate dehydratase subunit
PWY-5257 	4.2.1.82	KpLE2 phage-like element; D-xylonate dehydratase
PWY-5257 	5.1.3.4	EG12287-MONOMER
PWY-5257 	5.1.3.4	L-ribulose 5-phosphate 4-epimerase
PWY-5257 	5.3.1.25 	L-fucose isomerase
PWY-5257 	5.3.1.4	L-arabinose isomerase
PWY-5257 	5.3.1.5	xylose isomerase
PWY-5258 	2.1.1.251	3-(methylthio)propanoate:coenzyme M methyltransferase
PWY-5259 	2.1.1.251	methylthiol:coenzyme M methyltransferase
PWY-5261	2.1.1.252	MONOMER-12245
PWY-5261	2.1.1.253	MONOMER-12244
PWY-5265 	1.3.1.98	MONOMER-12254
PWY-5265 	2.3.2.16	MONOMER-15452
PWY-5265 	2.3.2.17	FemA
PWY-5265 	2.3.2.18	FemB monomer
PWY-5265	2.4.1.129	monofunctional peptidoglycan glycosyltransferase
PWY-5265	2.4.1.129 	PBP2
PWY-5265	2.4.1.227	UDP-N-acetylglucosamine--N-acetylmuramyl-(pentapeptide) pyrophosphoryl-undecaprenol N-acetylglucosamine transferase
PWY-5265 	2.5.1.7	UDP-N-acetylglucosamine 1-carboxyvinyltransferase 1
PWY-5265 	2.5.1.7	UDP-N-acetylglucosamine 1-carboxyvinyltransferase 2
PWY-5265	2.7.8.13	phospho-N-acetylmuramoyl-pentapeptide-transferase
PWY-5265	3.4.16.4	PBP1
PWY-5265	3.4.16.4	PBP2a
PWY-5265	3.4.16.4	PBP3
PWY-5265	3.5.2.6 	PBP4
PWY-5265 	5.1.1.3	MONOMER-15462
PWY-5265 	6.3.2.10	MONOMER-12255
PWY-5265 	6.3.2.4	D-alanine--D-alanine ligase
PWY-5265 	6.3.2.7	MONOMER-12253
PWY-5265 	6.3.2.8	MONOMER-15459
PWY-5265 	6.3.2.9	UDP-N-acetylmuramoyl-L-alanine--D-glutamate ligase
PWY-5266 	1.1.1.M33	<I>p</I>-cumic alcohol dehydrogenase
PWY-5266 	1.13.11 	<I>p</I>-cumate 2,3-dioxygenase
PWY-5266 	1.13.11.M5 	2,3-dihydroxy-<I>p</I>-cumate 3,4-dioxygenase
PWY-5266 	1.14.13.dv	<I>p</I>-cymene methyl-monooxygenase
PWY-5266 	1.2.1.29	MONOMER-324
PWY-5266 	1.3.1.58	MONOMER-348
PWY-5266 	3.7.1	HOMODA hydrolase
PWY-5266 	4.1.1.M8	HCOMODA decarboxylase
PWY-5267	3.2.1.147	AT5G25980-MONOMER
PWY-5267	3.2.1.147	AT5G26000-MONOMER
PWY-5267	3.2.1.147	myrosinase
PWY-5267	3.2.1.147	myrosinase 
PWY-5269	2.7.8.41	AT4G04870-MONOMER
PWY-5269	2.7.8.5	mitochondrial CDP-diacylglycerol--glycerol-3-phosphate 3-phosphatidyltransferase
PWY-5270	1.1.1.218 	codeinone reductase
PWY-5270	1.1.1.248	salutaridine oxidoreductase
PWY-5270	1.14.11.31	thebaine 6-O-demethylase
PWY-5270	1.14.11.32	codeine O-demethylase
PWY-5270	1.14.21.4	MONOMER-12299
PWY-5270	2.3.1.150	salutaridinol-7-O-acetyltransferase
PWY-5270	5.1.99	reticuline epimerase
PWY-5271	1.14.13.93	MONOMER-12323
PWY-5272 	2.4.1.177 	indole-3-acetic acid glucosyltransferase / abscisic acid glucosyltransferase
PWY-5272	2.4.1.177 	MONOMER-12304
PWY-5272	2.4.1.263	AT1G05530-MONOMER
PWY-5272	2.4.1.263	AT1G05560-MONOMER
PWY-5272	2.4.1.263	AT2G23250-MONOMER
PWY-5272	2.4.1.263	AT3G21780-MONOMER
PWY-5272	3.2.1.175	abscisic acid glucose ester &beta;-glucosidase
PWY-5274	1.8.2.3	flavocytochrome c cytochrome subunit 
PWY-5274	1.8.2.3	flavocytochrome c cytochrome subunit 
PWY-5276	1.8.2.1	sulfite dehydrogenase
PWY-5277	2.8.1.3	thiosulfate&mdash;thiol sulfurtransferase
PWY-5278	1.8.99.2	adenylylsulfate reductase &alpha; subunit 
PWY-5278	2.7.7.4	dissimilatory sulfate adenylyltransferase
PWY-5279 	1.8.99.2	APS reductase &alpha; subunit 
PWY-5279	2.7.7.5	adenylylsulfate:phosphate adenylyltransferase subunit
PWY-5283	1.4.3.3	MONOMER-17731
PWY-5284 	2.3.1.172	anthocyanin 5-O-glucoside-6-O-malonyltransferase
PWY-5284 	2.3.1.215 	anthocyanidin 3-O-glucoside 6-O-acyltransferase
PWY-5284 	2.4.1 	anthocyanin 5-O-glucosyltransferase
PWY-5287	1.14.13.56	MONOMER-12347
PWY-5287	1.14.13.57	MONOMER-12349
PWY-5287 	1.14.13 	(S)-cis-N-methylstylopine 14-hydroxylase/(S)-N-methylcanadine 14-hydroxylase
PWY-5287	1.14.21.1	MONOMER-12341
PWY-5287	1.14.21.1	stylopine synthase
PWY-5287	1.14.21.2	MONOMER-12340
PWY-5287	1.14.21.2	(S)- cheilanthifoline synthase
PWY-5287	1.14.21.e 	MONOMER-13846
PWY-5287	1.5.3.12	dihydrobenzophenanthridine oxidase
PWY-5287	2.1.1.119	MONOMER-12348
PWY-5287	2.1.1.120	MONOMER-12350
PWY-5287	2.1.1	MONOMER-12342
PWY-5288	1.13.99 	&beta;-carotene ketolase
PWY-5288	1.13.99 	&beta;-ring carotenoid 4-ketolase
PWY-5288	1.13.99 	MONOMER-12401
PWY-5288	1.14.13.129	&beta;-carotene hydroxylase
PWY-5290	1.1.1.324 	acyclic monoterpene primary alcohol oxidoreductase
PWY-5290	1.14.13.152	MONOMER-12352
PWY-5290	1.14.13.74	MONOMER-12353
PWY-5290	1.3.3.9	MONOMER-12354
PWY-5290	2.1.1.50	MONOMER-13906
PWY-5290	2.1.1.50	MONOMER-13907
PWY-5290	2.4.1.323	7-deoxyloganetic acid glucosyltransferase
PWY-5290	4.3.3.2	MONOMER-11582
PWY-5292	1.14.11.20	MONOMER-12361
PWY-5292	1.14.13.73	MONOMER-12358
PWY-5292	1.14.13.73	MONOMER-7709
PWY-5292	2.1.1.94	MONOMER-12359
PWY-5292	2.1.1.99	MONOMER-12360
PWY-5292	2.3.1.107	MONOMER-12364
PWY-5294 	1.8.2.1	sulfite dehydrogenase &alpha; subunit 
PWY-5294 	1.8.5.4	sulfide:quinone oxidoreductase
PWY-5297	6.3.5.M1	siroheme amidase
PWY-5301	1.1.1.273	MONOMER-12406
PWY-5301	1.14.13.75	MONOMER-17775
PWY-5301	1.14.13.91	MONOMER-12413
PWY-5301	1.3.1.73	MONOMER-12397
PWY-5301	1.5.1.32	MONOMER-12395
PWY-5301	2.3.1.160	MONOMER-12394
PWY-5301	3.1.1.78	MONOMER-7725
PWY-5301	3.1.1.80	acetylajmalan esterase
PWY-5301	3.2.1.105	strictosidine glucosidase
PWY-5301	3.2.1.125	MONOMER-12400
PWY-5303	2.8.1	DsrE3A thiosulfate-carrier protein
PWY-5303	2.8.1	TusA sulfur/thiosulfate-carrier protein
PWY-5304 	1.13.11.55	sulfur oxygenase reductase hexamer
PWY-5304 	1.8.2.1	sulfite dehydrogenase
PWY-5304 	1.8.5.2	thiosulfate:quinone oxidoreductase
PWY-5304	1.8.99.2	MONOMER-12420
PWY-5304	2.7.7.5	MONOMER-12421
PWY-5305 	1.13.12	MONOMER-12416
PWY-5305	1.2.1	MONOMER-12424
PWY-5305	2.1.1	MONOMER-12425
PWY-5306 	1.8.99.2	MONOMER-12409
PWY-5306 	2.7.7.4	MONOMER-12410
PWY-5308	1.8.2.1	sulfite dehydrogenase
PWY-5308	1.8.99.2	MONOMER-12414
PWY-5308	2.7.7.4	MONOMER-12415
PWY-5310 	2.3.1.171	MONOMER-12396
PWY-5310 	2.3.1 	anthocyanin 5-(6-hydroxycinnamoyltransferase)
PWY-5310 	2.4.1.298	anthocyanidin 5-O-glucosyltransferase
PWY-5310 	2.4.1.298 	anthocyanin 3-glucosyltransferase
PWY-5310 	2.4.1.298 	MONOMER-12013
PWY-5312 	2.3.1.214	MONOMER-12281
PWY-5312 	2.4.1.116	MONOMER-12017
PWY-5313 	2.4.1.115	anthocyanidin 3-O-glucosyltransferase
PWY-5313 	2.4.1.115	flavonoid 3-O-glucosyltransferase
PWY-5313 	2.4.1.298	flavonoid 5-O-glucosyltransferase
PWY-5313 	2.4.1	anthocyanidin 3-O-glucoside 2-O-glucosyltransferase
PWY-5313 	2.4.1	anthocyanidin 5,3-O-glycosyltransferase
PWY-5315	2.1.1.53	MONOMER-12438
PWY-5318	1.1.1.236	MONOMER-13851
PWY-5318	1.1.1.236	MONOMER-13853
PWY-5318	1.1.1.236	tropinone reductase (TR-II)
PWY-5319	3.1.1.35	MONOMER-12433
PWY-5320 	2.4.1.91	flavonol 3-O-glucosyltransferase
PWY-5321 	2.4.1.237 	flavonol 7-O-glucosyltransferase / abscisic acid glucosyltransferase
PWY-5321 	2.4.1.91 	flavonol-3-O-rhamnosyltransferase
PWY-5321 	2.4.1	flavonol 7-O-rhamnosyltransferase
PWY-5326	1.8.3.1	SO subunit
PWY-5327 	1.13.12.2	L-lysine monooxygenase subunit
PWY-5327 	1.2.1	MONOMER-12336
PWY-5327	1.3.8.1	short-chain acyl-CoA dehydrogenase
PWY-5327 	1.4.1.11	L-erythro-3,5-diaminohexanoate dehydrogenase subunit
PWY-5327 	1.4.1.11	MONOMER-12294
PWY-5327 	1.4.3.14	L-lysine &alpha;-oxidase subunit
PWY-5327	1.4.3.3 	D-amino acid oxidase subunit
PWY-5327 	1.5.1.1 	ketimine reductase mu-crystallin
PWY-5327 	1.5.1.1 	MONOMER-12351
PWY-5327 	1.5.1.9 	&alpha;-aminoadipic semialdehyde synthase, mitochondrial
PWY-5327 	1.5.1.9 	lysine-ketoglutarate reductase subunit
PWY-5327 	1.5.3.7	L-pipecolate oxidase
PWY-5327 	1.5.3.7	MONOMER-12355
PWY-5327 	2.3.1.247	3-keto-5-aminohexanoate cleavage enzyme subunit
PWY-5327 	2.3.1.247	3-keto,5-aminohexanoate cleavage enzyme subunit
PWY-5327 	2.3.1.8	MONOMER-12291
PWY-5327 	2.3.1.9	acetyl-CoA acetyltransferase
PWY-5327 	2.3.1.9	acetyl-CoA C-acetyltransferase
PWY-5327 	2.3.1.9	cytosolic acetyl-CoA acetyltransferase
PWY-5327 	2.3.1	MONOMER-8142
PWY-5327 	2.3.1	MONOMER-8145
PWY-5327 	2.3.1	MONOMER-8167
PWY-5327 	2.3.1	MONOMER-8221
PWY-5327 	2.6.1.36	L-lysine 6-aminotransferase subunit
PWY-5327 	2.6.1.36	MONOMER-12388
PWY-5327 	2.6.1.39	&alpha;-aminoadipate aminotransferase subunit
PWY-5327 	2.6.1.48	MONOMER-12326
PWY-5327 	2.6.1.48	MONOMER-8144
PWY-5327 	2.6.1.71	MONOMER-12465
PWY-5327 	2.6.1	N<sup>6</sup>-acetyllysine aminotransferase subunit
PWY-5327 	2.6.1	MONOMER-8143
PWY-5327	2.8.3.13	succinyl-CoA-glutarate CoA-transferase
PWY-5327 	2.8.3.9	butyryl-CoA:acetoacetate CoA-transferase
PWY-5327 	2.8.3.9	MONOMER-12296
PWY-5327 	3.5.1.30	MONOMER-12325
PWY-5327 	4.1.1.18	lysine decarboxylase subunit
PWY-5327 	4.3.1.14	MONOMER-12286
PWY-5327 	4.3.1.14	MONOMER-12295
PWY-5327 	5.1.1.5	MONOMER-12337
PWY-5327 	5.4.3.2	L-lysine 2,3-aminomutase subunit
PWY-5327 	5.4.3.2	lysine 2,3-aminomutase
PWY-5327 	5.4.3.3	L-&beta;-lysine-5,6-aminomutase &beta; subunit 
PWY-5327 	5.4.3.3	lysine 5,6-aminomutase &alpha; subunit 
PWY-5328 	1.13.11.20	cysteine dioxygenase
PWY-5328 	1.13.11.20	MONOMER-8844
PWY-5328 	1.8.2.1 	mitochondrial sulfite oxidase
PWY-5328 	2.1.1.10	S-methylmethionine--homocysteine S-methyltransferase
PWY-5328 	2.1.1.13 	methionine synthase
PWY-5328 	2.1.1.5	betaine-homocysteine S-methyltransferase subunit
PWY-5328 	2.1.1.5	betaine--homocysteine S-methyltransferase 1
PWY-5328 	5.1.99.1	methylmalonyl CoA epimerase
PWY-5328 	5.1.99.1	methylmalonyl-CoA racemase subunit
PWY-5328 	5.4.99.2	methylmalonyl-CoA mutase
PWY-5328 	5.4.99.2	MONOMER-8610
PWY-5328 	6.4.1.3	propionyl-CoA carboxylase
PWY-5328 	6.4.1.3	propionyl-CoA carboxylase &alpha; subunit 
PWY-5329	2.8.1.2	3-mercaptopyruvate sulfurtransferase
PWY-5329	2.8.1.2	MONOMER-12473
PWY-5331	1.13.11.19	cysteamine dioxygenase
PWY-5331	1.8.1.3	hypotaurine dehydrogenase
PWY-5331	4.1.1.29	cysteine sulfinic acid decarboxylase
PWY-5331	4.1.1.29	MONOMER-18720
PWY-5331	4.1.1.29	sulfinoalanine decarboxylase
PWY-5332	1.12.98.4	82 kDa subunit 
PWY-5332	1.12.98.4	sulfur reductase molybdopterin subunit 
PWY-5335 	1.8.2.1	sulfite dehydrogenase
PWY-5337	2.4.1.123	MONOMER-12620
PWY-5337	2.4.1.82	MONOMER-12495
PWY-5337	2.4.1 	MONOMER-12496
PWY-5337	2.4.2 	stachyose synthase
PWY-5338 	2.4.1 	galactosylononitol synthase / stachyose synthase
PWY-5339	2.4.1	MONOMER-12483
PWY-5339	2.4.1	MONOMER-12484
PWY-5340	2.7.1.25	AT2G14750-MONOMER
PWY-5340	2.7.1.25	AT4G39940-MONOMER
PWY-5340	2.7.7.4 	3-phosphoadenosine 5-phosphosulfate synthase 1
PWY-5340	2.7.7.4 	3-phosphoadenosine 5-phosphosulfate synthase 2
PWY-5340 	2.7.7.4	APS3
PWY-5340 	2.7.7.4	ATP sulfurylase
PWY-5340 	2.7.7.4	sulfate adenylyltransferase
PWY-5342 	2.4.1 	ajugose synthase [multifunctional]
PWY-5343	2.4.1	MONOMER-12486
PWY-5344	2.3.1.31	homoserine O-acetyltransferase
PWY-5345 	1.2.1.11	aspartate semialdehyde dehydrogenase subunit
PWY-5345 	2.1.1.13 	G18NG-11090-MONOMER
PWY-5345 	2.1.1.14	G18NG-10711-MONOMER
PWY-5345 	2.1.1.14	MONOMER-14559
PWY-5345 	2.1.1.14	YER091C-MONOMER
PWY-5345 	2.6.1.1	aspartate aminotransferase
PWY-5345 	2.7.2.4	aspartate kinase II
PWY-5345 	2.7.2.4	aspartate kinase III
PWY-5345 	2.7.2.4	aspartokinase &alpha; subunit 
PWY-5349 	1.14.11	4-coumaroyl 2-hydroxylase
PWY-5350	2.8.1.1	inner membrane thiosulfate sulfurtransferase
PWY-5350	2.8.1.1	thiosulfate sulfurtransferase
PWY-5350	2.8.1.1	thiosulfite-cleaving enzyme subunit
PWY-5350 	2.8.1.2 	EG11600-MONOMER
PWY-5350	2.8.1 	thiosulfate sulfurtransferase
PWY-5352	1.8.2.M1	thiosulfate reductase (cytochrome)
PWY-5353	1.14.19.30	acyl-lipid 5-desaturase
PWY-5353	1.14.19.47	acyl-lipid 6-desaturase
PWY-5353 	2.3.1	C18 &Delta;6-polyunsaturated fatty acyl-CoA elongase
PWY-5353	2.3.1	C18-&Delta;6-polyunsaturated fatty acyl-CoA elongase
PWY-5355	1.7.3.1	nitroalkane oxidase subunit
PWY-5358	1.7.2.2 	tetrathionate reductase
PWY-5360 	1.8.5.c	thiosulfate reductase (quinone)
PWY-5360	1.8.99.5	dissimilatory sulfite reductase &alpha; subunit 
PWY-5360 	1.97.1	tetrathionate reductase &alpha; subunit 
PWY-5361	1.14.19.10	acyl-CoA 5-desaturase
PWY-5362	1.14.19.26	acyl-[acyl-carrier protein] 6-desaturase
PWY-5363	1.14.11.22	flavone synthase I
PWY-5364	1.12.98.4	hydrogenase gamma subunit 
PWY-5364	1.12.98.4	polysulfide reductase &alpha; subunit 
PWY-5364	1.12.98.4	sulfhydrogenase II
PWY-5365	1.14.13.102	MONOMER-12571
PWY-5365	2.1.1.69	MONOMER-12575
PWY-5365	2.5.1	MONOMER-12568
PWY-5366	1.14.19.2	acyl-[acyl-carrier-protein] 9-desaturase
PWY-5366 	3.1.2.14	acyl-[acp] thioesterase
PWY-5367	1.14.19.11	acyl-[acp] 4-desaturase
PWY-5367	1.14.19.2 	acyl-[acyl-carrier-protein] 4/9-desaturase
PWY-5367	2.3.1.179	3-oxoacyl-[acyl-carrier-protein] synthase
PWY-5368	1.14.19.34	oleoyl-lipid 12-(E) desaturase
PWY-5368	1.14.99	(9Z,12E)-octadecadienoyl-lipid 9-hydroxylase
PWY-5370 	2.8.4.1	methyl-coenzyme M reductase
PWY-5372 	1.2.1.43	NADP-dependent formate dehydrogenase &alpha; subunit 
PWY-5373	1.14.19.14	&Delta;9 acyl-lipid conjugase
PWY-5374	1.14.19.16	linoleoyl-lipid &Delta;12 conjugase (11E,13Z-forming)
PWY-5375	1.14.19.33	acyl-lipid &Delta;12 conjugase
PWY-5375	1.14.19.33	acyl-lipid &Delta;12-conjugase
PWY-5375	1.14.19.34 	acyl-lipid &Delta;12-conjugase
PWY-5377 	5.4.99.40 	AT1G78960-MONOMER
PWY-5377 	5.4.99.54 	multifunctional triterpene synthase
PWY-5379	2.4.1.123 	MONOMER-12680
PWY-5379 	2.4.1 	MONOMER-12627
PWY-5380 	2.4.1.123 	galactinol synthase / fagopyritol synthase
PWY-5381	3.2.2.14	MONOMER-12666
PWY-5381	3.5.1.19	MONOMER-12671
PWY-5381	3.6.1.22	MONOMER-12667
PWY-5381	6.3.4.21	MONOMER-12669
PWY-5382	1.12.1.2	bidirectional hydrogenase
PWY-5382	1.12.1.2	[NiFe]-hydrogenase (soluble)
PWY-5384	2.4.1.7	MONOMER-12645
PWY-5384	2.4.1.7	sucrose phosphorylase subunit
PWY-5384	2.7.1.4	MONOMER-12642
PWY-5384	5.4.2.2	MONOMER-15669
PWY-5386	3.1.2.6	YDR272W-MONOMER
PWY-5386	4.4.1.5	glyoxalase I
PWY-5386	4.4.1.5	glyoxalase I dimer
PWY-5386	4.4.1.5	glyoxalase I monomer
PWY-5388	2.4.1.196	MONOMER-12672
PWY-5389 	1.13.11.53 	acireductone dioxygenase [iron(II)-requiring]
PWY-5390 	1.14.11.23	flavonol synthase
PWY-5390	2.4.1.159	flavonol-3-O-glucoside L-rhamnosyltransferase
PWY-5390	2.4.1.91	flavonoid 3-O-glucosyltransferase
PWY-5392	1.1.1	oxalosuccinate reductase subunit
PWY-5392	1.2.7.1	pyruvate:ferredoxin oxidoreductase &alpha; subunit 
PWY-5392	1.2.7.3	2-oxoglutarate synthase
PWY-5392	4.1.3.34	citryl-coA lyase subunit
PWY-5392	6.2.1.18	citryl-CoA synthetase small subunit 
PWY-5392	6.4.1.7	2-oxoglutarate carboxylase small subunit 
PWY-5393	1.3.1	MONOMER-12714
PWY-5393	2.3.1.212	4-hydroxybenzalacetone synthase
PWY-5393	2.3.1.212	BAS subunit
PWY-5393	2.3.1.212	benzalacetone synthase
PWY-5393 	2.3.1.212 	polyketide type III synthase
PWY-5394	1.13.11.29	4,5-DOPA dioxygenase extradiol
PWY-5398	2.4.1.271	crocetin glucosyltransferase
PWY-5398	2.4.1.271	UDP-glucose:crocetin 8,8-glucosyltransferase
PWY-5398	2.4.1.330	crocetin glucosyl ester glucosyltransferase
PWY-5399	1.10.3	tyrosinase
PWY-5405 	1.10.3	tyrosinase
PWY-5405 	1.13.11.29	MONOMER-12753
PWY-5405 	2.3.1	1-O-hydroxycinnamoyl-&beta;-glucose:amaranthin O-hydroxycinnamoyl transferase
PWY-5405 	2.4.1	MONOMER-12776
PWY-5405 	2.4.1	MONOMER-12777
PWY-5405 	2.4.1	MONOMER-12783
PWY-5405 	2.4.1	MONOMER-12784
PWY-5405 	2.4.1	MONOMER-12825
PWY-5412 	5.5.1.12 	chloroplastic isopimara-7,15-diene synthase
PWY-5415 	1.13.11.2	C23OALCAL-MONOMER
PWY-5416 	1.14.13.108	abietadiene hydroxylase
PWY-5416 	1.2.1.74	abietadienal dehydrogenase
PWY-5416 	4.2.3.18 	abietadiene synthase
PWY-5416 	4.2.3.32 	levopimaradiene/abietadiene synthase
PWY-5416 	4.2.3 	levopimaradiene/abietadiene synthase
PWY-5417	2.3.1.174	&beta;-ketoadipyl CoA thiolase subunit
PWY-5418 	1.14.13.7	phenol 2-monooxygenase
PWY-5420 	1.2.1.10	acetaldehyde dehydrogenase  (acylating)
PWY-5420 	4.1.3.39	4-hydroxy-2-oxovalerate aldolase
PWY-5422 	1.14.13	abietadienol hydroxylase
PWY-5424 	4.2.3.107 	R-linalool synthase
PWY-5424 	4.2.3.117 	camphene synthase
PWY-5424 	4.2.3.119 	terpinolene synthase
PWY-5424 	4.2.3.120 	(4S)-limonene synthase
PWY-5424 	4.2.3.120 	pinene synthase
PWY-5424 	4.2.3.15	myrcene synthase
PWY-5424 	4.2.3.20 	(E)-&alpha;-bisabolene synthase
PWY-5424 	4.2.3.20 	pinene synthase
PWY-5424 	4.2.3.38	(E)-&alpha;-bisabolene synthase
PWY-5424 	4.2.3.52 	(4S)-limonene synthase
PWY-5424 	4.2.3.52 	&beta;-phellandrene synthase
PWY-5424 	4.2.3 	longifolene synthase
PWY-5425	4.2.3.76	MONOMER-14954
PWY-5430 	1.13.11.2	catechol 2,3-dioxygenase
PWY-5430 	3.7.1.9	MONOMER-12763
PWY-5	4.3.2.1	argininosuccinate lyase
PWY-5433 	1.1.1.21 	aldo-keto reductase family 4 member C9
PWY-5433 	1.1.1	aldehyde reductase
PWY-5433 	1.1.1	aldehyde reductase (chloroplastic)
PWY-5433 	1.1.1 	aldehyde reductase (cytosolic)
PWY-5433 	1.13.11.12	linoleate 13S-lipoxygenase
PWY-5433 	1.13.11.12	lipoxygenase
PWY-5433	1.13.11.12	lipoxygenase
PWY-5433 	1.13.11 	13S/9S lipoxygenase
PWY-5433 	1.13.11 	9S-lipoxygenase
PWY-5433 	1.13.11 	9S-lipoxygenase 1
PWY-5433 	1.3.1.42	12-oxophytodienoate reductase
PWY-5433 	1.3.3.6	MONOMER-14969
PWY-5433 	2.3.1.16 	AT2G33150-MONOMER
PWY-5433 	2.3.1.195	acetyl CoA:(Z)-3-hexen-1-ol acetyltransferase
PWY-5433 	4.2.1.121	9-divinyl ether synthase
PWY-5433 	4.2.1.92	allene oxide synthase
PWY-5433 	5.3.99.6	allene oxide cyclase
PWY-5433 	6.2.1	AT1G20510-MONOMER
PWY-5434	4.2.3.48	MONOMER-12826
PWY-5436 	4.1.2.48 	low-specificity L-threonine aldolase
PWY-5439 	1.10.3 	tyrosinase
PWY-5441 	2.1.1.10 	homocysteine S-methyltransferase
PWY-5441 	2.1.1.12	AT5G49810-MONOMER
PWY-5450	1.14.12.3	benzene 1,2 dioxygenase terminal oxygenase component 
PWY-5450	1.3.1.19	cis-benzene glycol dehydrogenase subunit
PWY-5453 	1.1.1	aldose reductase
PWY-5456	1.2.1.23	MONOMER-12919
PWY-5456	1.2.1.23	MONOMER-12930
PWY-5458	1.1.2.3	L-lactate dehydrogenase subunit
PWY-5458	1.2.1.22	MONOMER-12910
PWY-5461	1.11.1.7	MONOMER-12917
PWY-5462	1.2.1.49	MONOMER-12918
PWY-5464	1.1.1.37	MONOMER-17398
PWY-5464	1.1.1.37	MONOMER-17400
PWY-5464	1.1.1.37	MONOMER-17411
PWY-5464 	1.1.1.41	NAD+-dependent isocitrate dehydrogenase
PWY-5464 	1.2.1.12	MONOMER-12868
PWY-5464 	1.2.1.12	NAD-dependent glyceraldehyde-3-phosphate dehydrogenase
PWY-5464 	1.2.1.9	cytosolic glyceraldehyde 3-phosphate dehydrogenase
PWY-5464 	1.2.1.9	cytosolic NADP-dependent glyceraldehyde-3-phosphate dehydrogenase
PWY-5464 	1.3.5.1	succinate dehydrogenase complex
PWY-5464 	2.7.1.11	6-phosphofructokinase subunit
PWY-5464 	2.7.1.11	cytosolic 6-phosphofructokinase
PWY-5464 	2.7.1.40	cytosolic pyruvate kinase
PWY-5464 	2.7.1.90	diphosphate:fructose-6-phosphate 1-phosphotransferase
PWY-5464 	2.7.1.90	diphosphate&mdash;fructose-6-phosphate 1-phosphotransferase
PWY-5464 	2.7.2.3	cytosolic 3-phosphoglycerate kinase
PWY-5464 	4.1.2.13	cytosolic fructose-1,6-biphosphate aldolase
PWY-5464 	4.1.2.13	cytosolic fructose-bisphosphate aldolase
PWY-5464 	4.1.2.13	fructose-1,6-biphosphate aldolase
PWY-5464 	4.1.2.13	fructose-1,6-phosphate aldolase
PWY-5464 	4.1.2 	cytosolic fructose-1,6-bisphosphate aldolase
PWY-5464 	4.2.1.11	cytoplasmic enolase
PWY-5464 	4.2.1.11	cytosolic enolase
PWY-5464 	4.2.1.11	MONOMER-12858
PWY-5464 	4.2.1.2	fumarase
PWY-5464 	5.3.1.1	cytosolic triose phosphate isomerase
PWY-5464 	5.3.1.1	cytosolic triose-phosphate isomerase
PWY-5464 	5.3.1.1	cytosolic triosephosphate isomerase
PWY-5464 	5.4.2.12	cytoplasmic 2,3-bisphosphoglycerate-independent phosphoglycerate mutase
PWY-5464 	5.4.2.12	cytosolic cofactor-independent phosphoglycerate mutase
PWY-5464 	6.2.1.5	MONOMER-12976
PWY-5466	1.1.1	(-)-secoisolariciresinol dehydrogenase
PWY-5466	1.23.1.2 	(+)-pinoresinol/(+)-lariciresinol reductase
PWY-5466	1.23.1.4 	(+)-pinoresinol/(+)-lariciresinol reductase
PWY-5467	2.1.1	MONOMER-12933
PWY-5468	2.3.1.93	MONOMER-12941
PWY-5468	2.3.1.93	MONOMER-13855
PWY-5468	4.1.1.18	MONOMER-12938
PWY-5469	1.14.13	(+)-piperitol/(+)-sesamin synthase
PWY-5470 	1.3.3.8	(S)-tetrahydroprotoberberine oxidase
PWY-5470	2.1.1.118 	columbamine O-methyltransferase
PWY-5470	2.1.1.118	MONOMER-12942
PWY-5472	1.14.21 	berbamunine synthase
PWY-5473	1.14.13	MONOMER-12983
PWY-5474 	2.3.1.110	tyramine hydroxycinnamoyltransferase
PWY-5474 	2.3.1.110	tyramine N-hydroxycinnamoyltransferase
PWY-5474	2.3.1.110	tyramine N-hydroxycinnamoyltransferase
PWY-5474 	4.1.1.25 	MONOMER-12764
PWY-5474 	4.1.1.28 	MONOMER-8382
PWY-5478 	1.10.3	MONOMER-13001
PWY-5478 	1.10.3	tellimagrandin II O2 oxidoreductase subunit
PWY-5478 	2.3.1.143	MONOMER-13000
PWY-5478 	2.3.1.90	MONOMER-12995
PWY-5478 	2.3.1	&beta;-glucogallin: 1,2,3,4,6-pentagalloyl-&beta;-D-glucose (2,4-O-digalloyl)-galloyltransferase
PWY-5478 	2.3.1	&beta;-glucogallin: 3-O-trigalloyl-1,2,4,6-tetra-O-galloylglucose O-galloyltransferase
PWY-5478 	2.3.1	galloyltransferase C subunit
PWY-5478 	2.3.1	MONOMER-12996
PWY-5478 	2.3.1	MONOMER-12997
PWY-5478 	2.3.1	MONOMER-12999
PWY-5478 	2.3.1	MONOMER-13005
PWY-5478 	2.3.1	MONOMER-13006
PWY-5478 	2.4.1.136	gallate 1-&beta;-glucosyltransferase
PWY-5478 	2.4.1.136	MONOMER-12993
PWY-5478 	2.4.1.136	MONOMER-12994
PWY-5479	1.14.13	MONOMER-13010
PWY-5479	2.1.1	MONOMER-13011
PWY-5480 	2.3.1.54	pyruvate formate-lyase subunit
PWY-5481 	1.1.1.27	L-lactate dehydrogenase subunit
PWY-5486 	4.1.1.1 	pyruvate decarboxylase
PWY-5486 	4.1.1.1 	pyruvate decarboxylase 1
PWY-5486 	4.1.1.1 	pyruvate decarboxylase 2
PWY-5486 	4.1.1.1 	pyruvate decarboxylase 3
PWY-5487	1.13.11.66	hydroquinone 1,2-dioxygenase
PWY-5487	1.13.11.66	hydroquinone dioxygenase
PWY-5487	1.14.13.166 	4-nitrophenol 4-monooxygenase
PWY-5487	1.2.1.61	4-hydroxymuconic-semialdehyde dehydrogenase
PWY-5487	1.2.1.61	MONOMER-13027
PWY-5487	1.3.1.32	MONOMER-13029
PWY-5487	1.6.5.7 	4-benzoquinone reductase
PWY-5488	1.13.11.37	MONOMER-13015
PWY-5488	1.14.13.29 	4-nitrophenol 2-monooxygenase
PWY-5488	1.14.13.29	4-nitrophenol 2-monooxygenase
PWY-5488 	1.3.1.32	maleylacetate reductase
PWY-5488	1.3.1.32	MONOMER-13024
PWY-5489	3.1.8.1	methyl parathion hydrolase
PWY-5490 	3.1.8.1	organophosphate acid anhydrase
PWY-5491	3.1.3	phosphodiesterase subunit
PWY-5497	1.2.7.1	pyruvate-ferredoxin oxidoreductase subunit
PWY-5497	1.5.1.5	methylenetetrahydrofolate dehydrogenase subunit
PWY-5497	2.1.2.1	MONOMER-13128
PWY-5497	2.3.1.8	MONOMER-13134
PWY-5497 	3.5.4.9	MONOMER-13130
PWY-5499	1.1.1.107	pyridoxal 4-dehydrogenase subunit
PWY-5499	1.1.3.12	MONOMER-13147
PWY-5499	1.14.12.4	MONOMER-13153
PWY-5499	2.6.1.30	pyridoxamine-pyruvate aminotransferase subunit
PWY-5499	3.1.1.27	4-pyridoxolactonase subunit
PWY-5499	3.5.1.29	MONOMER-13154
PWY-5499	4.1.1.51	MONOMER-13152
PWY-5505 	1.4.1.14	glutamate synthase (NADH)
PWY-5505 	1.4.1.14	MONOMER-7421
PWY-5505 	1.4.1.3 	glutamate dehydrogenase 1
PWY-5505 	1.4.1.3 	glutamate dehydrogenase 2
PWY-5505 	1.4.1.4	glutamate dehydrogenase
PWY-5505 	1.4.1.4	NADP-dependent glutamate dehydrogenase 2 subunit
PWY-5505 	1.4.1.4	NADP-specific glutamate dehydrogenase 1 subunit
PWY-5505 	6.3.1.2	glutamine synthetase, type II
PWY-5505 	6.3.1.2	glutamine synthetase, type III
PWY-5505 	6.3.1.2	glutamine synthetase, type IV
PWY-5506	1.11.1.6 	catalase (peroxisomal)
PWY-5506	1.1.3.13	alcohol oxidase
PWY-5507	1.16.1.3 	MONOMER-13237
PWY-5507 	1.16.8.1	flavodoxin A
PWY-5507 	1.3.1.106	cobalt-precorrin-6A reductase
PWY-5507 	1.3.1.106	MONOMER-13224
PWY-5507 	1.3.1.76	precorrin-2 dehydrogenase
PWY-5507 	2.1.1.151	cobalt-precorrin-2 C(20)-methyltransferase subunit
PWY-5507 	2.1.1.195	cobalt-precorrin-6A synthase
PWY-5507 	2.1.1.196	cobalt-precorrin-6B (C<sup>15</sup>)-methyltransferase [decarboxylating]
PWY-5507 	2.1.1.271	MONOMER-13218
PWY-5507 	2.1.1.272	MONOMER-17924
PWY-5507 	2.1.1.289	cobalt-precorrin-7 (C<sup>5</sup>)-methyltransferase
PWY-5507 	2.4.2.21	nicotinate-nucleotide--dimethylbenzimidazole phosphoribosyltransferase subunit
PWY-5507 	2.7.1.156 	CobU
PWY-5507 	2.7.1.177	MONOMER-13230
PWY-5507 	2.7.8.26	adenosylcobalamin 5-phosphate synthase
PWY-5507 	3.1.3.73	&alpha;-ribazole-5-phosphate phosphatase
PWY-5507 	3.7.1.12	MONOMER-13219
PWY-5507 	4.1.1.81	L-threonine-O-3-phosphate decarboxylase subunit
PWY-5507 	4.99.1.3	cobaltochelatase
PWY-5507 	4.99.1.4 	sirohydrochlorin cobaltochelatase
PWY-5507 	4.99.1.4 	uroporphyrinogen-III C-methyltransferase subunit
PWY-5507 	5.4.3.8	glutamate-1-semialdehyde 2,1-aminomutase subunit
PWY-5507 	5.4.99.60	cobalt-precorrin-8 methylmutase
PWY-5507 	6.3.1.10	MONOMER-13228
PWY-5507 	6.3.5.10	adenosyl-cobyrate synthase subunit
PWY-5507 	6.3.5.11	cobyrinate a,c-diamide synthase
PWY-5508	3.1.3.73	MONOMER-13245
PWY-5512 	5.1.3.7 	UDP-glucose 4-epimerase
PWY-5512	5.1.3.7	UDP-N-acetylglucosamine C4-epimerase subunit
PWY-5514	2.3.1.4	MONOMER-13186
PWY-5514	2.7.1.1 	glucokinase
PWY-5514	2.7.7.23	MONOMER-13188
PWY-5514	3.5.99.6	glucosamine 6-phosphate isomerase 1
PWY-5514	5.1.3.7	MONOMER-13189
PWY-5514	5.3.1.9	glucose-6-phosphate isomerase
PWY-5514	5.4.2.3	phosphoacetylglucosamine mutase
PWY-5521 	1.1.99.32 	L-sorbose/L-sorbosone dehydrogenase 64.5 kD subunit 
PWY-5523	1.5.1.39 	FMN reductase [NAD(P)H]
PWY-5525 	1.1.1.19	aldo-keto reductase family 1, member A4 (aldehyde reductase)
PWY-5525	1.1.1.45	L-gulonate 3-dehydrogenase subunit
PWY-5525	1.1.1 	L-xylulose reductase subunit
PWY-5525	4.1.1.34	MONOMER-13239
PWY-5529 	1.1.1.396	bacteriochlorophyllide a dehydrogenase
PWY-5529 	1.3.1.75	3,8-divinyl protochlorophyllide a 8-vinyl-reductase (NADPH)
PWY-5529 	1.3.7.13 	chlorophyllide a reductase
PWY-5529 	1.3.7.15	chlorophyllide a reductase
PWY-5529 	1.3.7.7	light-independent protochlorophyllide reductase
PWY-5529 	1.3.98.3	coproporphyrinogen dehydrogenase
PWY-5529 	2.1.1.11	MONOMER-13249
PWY-5529 	2.5.1.bk	bacteriochlorophyll a synthase
PWY-5529 	4.2.1.165	chlorophyllide a hydratase
PWY-5529 	6.6.1.1	BchI-BchD complex 
PWY-5530	1.1.99.28	glucose-fructose oxidoreductase
PWY-5530	2.7.1.12	MONOMER-13275
PWY-5530	3.1.1.17	gluconolactonase subunit
PWY-5532	2.4.2.15 	MONOMER-19658
PWY-5532	2.4.2.1	adenosine phosphorylase
PWY-5532	2.4.2.2 	uridine phosphorylase
PWY-5532	2.4.2.57	AMP phosphorylase
PWY-5532	2.7.1.av	&alpha;-D-ribose-1-phosphate 5-kinase (ADP)
PWY-5532	2.7.1.aw	cytidine kinase
PWY-5532	4.1.1.39	RbcL dimer
PWY-5532	5.3.1.29	MONOMER-13273
PWY-5533	6.4.1.6	acetone carboxylase &beta; subunit 
PWY-5534	1.1.1.268	2-(R)-hydroxypropyl-CoM dehydrogenase subunit
PWY-5534	1.1.1.269	2-(S)-hydroxypropyl-CoM dehydrogenase subunit
PWY-5534	1.14.13.69	alkene monooxygenase
PWY-5534	1.8.1.5	2-oxopropyl-CoM reductase subunit
PWY-5534	4.4.1.23	2-hydroxypropyl-CoM lyase
PWY-5537	1.2.1 	pyruvate dehydrogenase complex
PWY-5537 	2.8.3.18	MONOMER-13300
PWY-5537	6.2.1.5	succinyl-CoA synthetase &alpha; subunit 
PWY-5538	1.2.7.1	pyruvate:ferredoxin oxidoreductase subunit
PWY-5538 	2.8.3.18	acetyl:succinate CoA-transferase
PWY-5538	6.2.1.5	succinyl-CoA synthetase &alpha;-2 subunit 
PWY-5600	1.2.7.1	MONOMER-13306
PWY-5600	2.3.1.8	MONOMER-13305
PWY-5600	2.7.2.1	MONOMER-13304
PWY-561 	1.1.1.35 	peroxisomal fatty acid &beta;-oxidation multifunctional protein MFP2
PWY-561 	1.1.1.37	AT3G47520-MONOMER
PWY-561 	1.1.1.37	malate dehydrogenase
PWY-561 	2.3.3.16 	citrate synthase (mitochondrial)
PWY-561 	2.3.3.9	malate synthase
PWY-561 	4.1.3.1	isocitrate lyase
PWY-5634 	1.14.11.26 	deacetoxycephalosporin C synthase / deacetoxycephalosporin C hydroxylase
PWY-5634 	1.14.11.26	MONOMER-13408
PWY-5634 	1.14.20.1	MONOMER-13407
PWY-5634 	1.21.3.1	MONOMER-13364
PWY-5634 	1.21.3.1	MONOMER-13365
PWY-5634 	2.1.3.7	MONOMER-13420
PWY-5634 	2.3.1.175	MONOMER-13419
PWY-5634 	5.1.1.17	MONOMER-13378
PWY-5634 	5.1.1.17	MONOMER-13381
PWY-5	6.3.4.5	argininosuccinate synthetase
PWY-5634 	6.3.2.26	MONOMER-13331
PWY-5636	1.14.13.31	MONOMER-13318
PWY-5637	1.7.1	MONOMER-13319
PWY-5637 	1.7.1	nitrobenzene nitroreductase
PWY-5637 	5.4.4.1	hydroxylaminobenzene mutase
PWY-5637	5.4.4.1	MONOMER-13322
PWY-5640	1.14.12.23 	nitroarene dioxygenase complex
PWY-5641	1.13.11.2	MONOMER-13375
PWY-5641	1.14.12.23	2-nitrotoluene dioxygenase complex
PWY-5641	3.7.1	2-hydroxy-6-oxohepta-2,4-dienoate hydrolase subunit
PWY-5642	1.13.11	2,4,5-trihydroxytoluene dioxygenase
PWY-5642	1.14.12.24	2,4-dinitrotoluene dioxygenase complex
PWY-5642	1.14.13.210	MONOMER-13337
PWY-5642	5.3.2 	4-hydroxy-2-keto-5-methyl-6-oxo-3-hexenoate isomerase/hydrolase
PWY-5643	1.14.12.23	MONOMER-13342
PWY-5644	1.13.11.74	MONOMER-13352
PWY-5644	1.2.1.32	MONOMER-13353
PWY-5644	1.7.1	MONOMER-13350
PWY-5644	3.5.99.5	2-aminomuconate deaminase
PWY-5644	4.1.1.77	MONOMER-13355
PWY-5644	4.2.1	MONOMER-13356
PWY-5644	5.4.4.1	MONOMER-13351
PWY-5645	1.13.11.74 	2-aminophenol 1,6-dioxygenase
PWY-5645	1.2.1.32	2-aminomuconate semialdehyde dehydrogenase
PWY-5645	3.5.1.120	2-aminomuconate deaminase
PWY-5645	4.1.1.77 	MONOMER-14737
PWY-5645	4.2.1.80	2-oxopent-4-enoate hydratase
PWY-5645	5.3.2.6 	2-hydroxymuconate tautomerase
PWY-5647	1.13.11.6	3-hydroxyanthranilate 3,4-dioxygenase subunit
PWY-5647 	1.2.1.32	MONOMER-13361
PWY-5647 	3.5.99.5	MONOMER-13362
PWY-5647 	4.1.1.45	MONOMER-13360
PWY-5647 	4.1.1.77	MONOMER-13363
PWY-5648	1.14.12.1	MONOMER-13368
PWY-5655 	1.13.11.6	3-hydroxyanthranilate 3,4-dioxygenase subunit
PWY-5655 	1.2.1.32	MONOMER-13372
PWY-5655 	3.5.99.5	2-aminomuconate deaminase subunit
PWY-5655 	4.1.1.45	MONOMER-13371
PWY-5656	2.4.1.217	MONOMER-13380
PWY-5656	3.1.3.70	MONOMER-13379
PWY-5658	2.4.1.269	mannosylglycerate synthase subunit
PWY-5659	2.7.7.13 	MONOMER-13383
PWY-5659	5.4.2.8	MONOMER-13382
PWY-5660	1.14.13.146	MONOMER-13409
PWY-5660	1.14.13.76	MONOMER-13401
PWY-5660	1.14.13.77	MONOMER-13398
PWY-5660	1.14.13	MONOMER-17465
PWY-5660	1.14.99.37	taxadiene 5&alpha;-hydroxylase
PWY-5660	2.3.1.162	MONOMER-13400
PWY-5660	2.3.1.162	MONOMER-14815
PWY-5660	2.3.1.166	MONOMER-13413
PWY-5660	2.3.1.167	MONOMER-13414
PWY-5660	2.3.1.167	MONOMER-17464
PWY-5660	2.3.1	MONOMER-13417
PWY-5660	2.3.1	MONOMER-13418
PWY-5660	4.2.3.17	MONOMER-14814
PWY-5660 	4.2.3.17	taxadiene synthase
PWY-5661 	2.7.1.1 	AT2G19860-MONOMER
PWY-5661 	2.7.7.34 	GDP-mannose pyrophosphorylase &alpha; subunit 
PWY-5661	2.7.7.34	MONOMER-13386
PWY-5661 	5.4.2.2	phosphoglucomutase-1
PWY-5662	2.4.1.266	glucosyl-3-phosphoglycerate synthase
PWY-5662	3.1.3.85	MONOMER-13395
PWY-5663 	3.5.4.16	GTP cyclohydrolase I
PWY-5663 	4.2.3.12	6-pyruvoyl tetrahydrobiopterin synthase trimer
PWY-5663 	4.2.3.12	6-pyruvoyl-tetrahydropterin synthase subunit
PWY-5664 	1.1.1.153	sepiapterin reductase
PWY-5664	1.1.1.220	MONOMER-13405
PWY-5665	1.14.14.1	MONOMER-13425
PWY-5665	2.1.1	MONOMER-13426
PWY-5665	3.7.1	MONOMER-13423
PWY-5666	2.4.1	glycosterol rhamnosyltransferase
PWY-5666	2.4.1	MONOMER-13432
PWY-5666	2.4.1	MONOMER-13434
PWY-5667 	2.3.1.51	1-acylglycerol-3-phosphate O-acyltransferase
PWY-5667	2.3.1.51	1-acyl-sn-glycerol-3-phosphate acyltransferase beta
PWY-5667	2.3.1.51	1-acyl-sn-glycerol-3-phosphate acyltransferase delta
PWY-5667	2.3.1.51	1-acyl-sn-glycerol-3-phosphate acyltransferase epsilon
PWY-5667	2.3.1.51	1-acyl-sn-glycerol-3-phosphate acyltransferase gamma
PWY-5667	2.7.7.41	phosphatidate cytidylyltransferase 1
PWY-5667	2.7.7.41	phosphatidate cytidylyltransferase 2
PWY-5668 	2.7.8.5	CDP-diacylglycerol--glycerol-3-phosphate 3-phosphatidyltransferase
PWY-5668 	3.1.3.27	phosphatidylglycerophosphatase
PWY-5670	1.14.14.17	MONOMER-13443
PWY-5670	1.14.14.17	MONOMER-14485
PWY-5670	2.5.1.21 	squalene synthase
PWY-5672	1.14.13.183	dammarenediol 12-hydroxylase
PWY-5672	1.14.13.184	protopanaxadiol 6-hydroxylase
PWY-5672	2.4.1.314	MONOMER-13452
PWY-5672	4.2.1.125	MONOMER-13444
PWY-5674	1.7.2.2	nitrite reductase (dissimilatory)
PWY-5674	1.9.6.1	periplasmic nitrate reductase (cytochrome)
PWY-5675	1.4.1.4	glutamate dehydrogenase (NADP+)
PWY-5675	1.7.1.3	assimilatory nitrate reductase (NADPH)
PWY-5675 	1.7.1.4 	assimilatory nitrite reductase
PWY-5675	6.3.1.2	glutamine synthetase, type II (octamer)
PWY-5675	6.3.1.2	glutamine synthetase, type II (tetramer)
PWY-5676	1.1.1.36	3-hydroxybutyryl-CoA dehydrogenase subunit
PWY-5676 	2.3.1.8	MONOMER-13471
PWY-5676 	2.7.2.1	MONOMER-13473
PWY-5676 	2.8.3.1 	MONOMER-13474
PWY-5677	1.1.1.61	4-hydroxybutyrate dehydrogenase (NAD)
PWY-5677	1.2.1.76	succinate semialdehyde dehydrogenase subunit
PWY-5677	2.8.3.18	succinyl-CoA:CoA transferase
PWY-5677	2.8.3	MONOMER-13467
PWY-5677	5.3.3.3 	4-hydroxybutyryl-CoA dehydratase
PWY-5679	1.14.11.21	clavaminate synthase
PWY-5679	2.5.1.66	N<sup>2</sup>-(carboxyethyl)arginine synthase subunit
PWY-5679	3.5.3.22	proclavaminate amidino hydrolase subunit
PWY-5679	6.3.3.4	MONOMER-13483
PWY-5690	1.1.1.41	isocitrate dehydrogenase [NAD]
PWY-5690 	1.3.5.1	succinate dehydrogenase (ubiquinone)
PWY-5691	1.7.3.3	factor-independent urate hydroxylase
PWY-5694	3.5.3.4	MONOMER-3152
PWY-5695 	3.1.3.5 	cytosolic 5-nucleotidase II
PWY-5701	1.14.13.116	geranylhydroquinone 3-hydroxylase
PWY-5701	2.5.1.93	MONOMER-13508
PWY-5701	2.5.1.93	MONOMER-13509
PWY-5704	3.5.1.5	urease &alpha; subunit 
PWY-5705 	3.5.3.9	allantoate amidohydrolase subunit
PWY-5705	4.3.2.3	ureidoglycolate lyase
PWY-5706	4.4.1.4	MONOMER-13494
PWY-5710	2.1.1.104	MONOMER-13528
PWY-5710	2.6.1	MONOMER-13529
PWY-5710	4.1.2.41 	hydroxycinnamoyl-CoA hydratase lyase
PWY-5717 	1.1.1.348	MONOMER-5044
PWY-5717 	1.1.1.348	MONOMER-5841
PWY-5717 	1.1.1.348	MONOMER-5884
PWY-5717 	1.14.13.52	MONOMER-5361
PWY-5717 	1.14.21.8	pseudobaptigenin synthase
PWY-5717 	1.3.1.45	MONOMER-5043
PWY-5717 	1.3.1.45	MONOMER-5681
PWY-5717 	1.3.1.45	MONOMER-5883
PWY-5717 	2.3.1.115	MONOMER-6244
PWY-5717 	2.4.1.170	MONOMER-6243
PWY-5717 	4.2.1.139	7,2-dihydroxy-4-methoxy-isoflavanol dehydratase
PWY-5717 	4.2.1.139	cis-(-)-7,2-dihydroxy-4,5-methylenedioxyisoflavanol dehydratase
PWY-5723 	2.7.1.19	chloroplastic phosphoribulokinase
PWY-5723 	4.1.1.39 	ribulose bisphosphate carboxylase/oxygenase
PWY-5724 	3.5.1.84 	allophanate hydrolase subunit
PWY-5724 	3.5.1.84	MONOMER-12131
PWY-5724 	3.5.1.84	MONOMER-13538
PWY-5724 	3.5.2.15	cyanurate amidohydrolase subunit
PWY-5724 	3.5.4.42	N-isopropylammelide isopropylaminohydrolase subunit
PWY-5724 	3.5.4.q	hydroxyatrazine ethylaminohydrolase
PWY-5724 	3.8.1.8	AtzA
PWY-5724 	3.8.1.8	triazine hydrolase
PWY-5725	4.2.3.46 	&alpha;-farnesene synthase
PWY-5725 	4.2.3.46 	AT4G16740-MONOMER
PWY-5725	4.2.3.46	(E,E)-&alpha;-farnesene synthase
PWY-5725	4.2.3.46	MONOMER-13546
PWY-5725	4.2.3.47	(E)-&beta;-farnesene synthase
PWY-5725	4.2.3.47	MONOMER-13549
PWY-5725	4.2.3.47	MONOMER-14489
PWY-5725	4.2.3 	sesquiterpene synthase
PWY-5729	1.3.1	MONOMER-13553
PWY-5729	1.3.1	MONOMER-13554
PWY-5731	3.8.1.8	triazine hydrolase subunit
PWY-5733 	4.2.3.23	germacrene A synthase
PWY-5733	4.2.3.23	germacrene A synthase
PWY-5733	4.2.3.23	MONOMER-13576
PWY-5733 	4.2.3.23	MONOMER-14837
PWY-5733	4.2.3.60	germacrene C synthase
PWY-5733	4.2.3.71	germacrene B synthase
PWY-5733	4.2.3.71	MONOMER-14864
PWY-5733	4.2.3.75	germacrene D synthase
PWY-5733	4.2.3.75	MONOMER-13577
PWY-5733	4.2.3.75	MONOMER-17487
PWY-5733	4.2.3.77	MONOMER-13556
PWY-5733	4.2.3.77	MONOMER-13578
PWY-5737	1.14.20.3	carbapenem synthase
PWY-5737	2.3.1.226	carboxymethylproline synthase subunit
PWY-5737	6.3.3.6	carbapenam synthetase subunit
PWY-5741	1.3.1.86 	crotonyl-CoA carboxylase/reductase
PWY-5741	1.3.8.12	(2S)-methylsuccinyl-CoA dehydrogenase
PWY-5741	4.1.3.24	malyl-CoA lyase subunit
PWY-5741	4.2.1.148	MONOMER-13587
PWY-5741	5.4.99.63	(2R)-ethylmalonyl-CoA mutase
PWY-5742 	1.2.1.54 	MONOMER-16
PWY-5742	2.6.1.84	L-arginine:pyruvate transaminase subunit
PWY-5742	3.5.3.7	guanidinobutyrase subunit
PWY-5742	4.1.1.75	MONOMER-13584
PWY-5744 	1.2.1.75	malonyl-CoA reductase subunit
PWY-5744 	1.3.1.84	acrylyl-CoA reductase (NADPH)
PWY-5744 	4.1.3.25 	L-malyl CoA lyase
PWY-5747	2.3.3.5	MONOMER-13635
PWY-5747	4.1.3.30	MONOMER-13636
PWY-5747	4.2.1.117	MONOMER-13619
PWY-5748	1.5.1	MONOMER-13624
PWY-5748	2.1.1	MONOMER-13626
PWY-5748	2.6.1	MONOMER-17488
PWY-5748	2.6.1	MONOMER-17489
PWY-5749	2.8.3.22 	MONOMER-13628
PWY-5749	4.1.3.25	MONOMER-13630
PWY-5749	4.2.1.56	MONOMER-13629
PWY-5749	6.2.1.4	(GDP-forming) succinyl-CoA synthetase
PWY-5750	2.3.3.16 	MONOMER-13633
PWY-5750	4.1.1.6	cis-aconitate decarboxylase
PWY-5750	6.4.1.1	pyruvate carboxylase monomer
PWY-5751	1.1.1	MONOMER-13638
PWY-5751	1.1.1	MONOMER-13639
PWY-5751	1.1.1 	MONOMER-13644
PWY-5751	1.4.3.21	MONOMER-13652
PWY-5751	4.1.1.53	MONOMER-13640
PWY-5751	4.1.1.53	MONOMER-13641
PWY-5751	4.1.1.53	MONOMER-13642
PWY-5752	2.3.1.145	MONOMER-13650
PWY-5755	4.1.3.40	chorismate lyase
PWY-5755 	4.1.3.40	MONOMER-18799
PWY-5756	2.4.1.17	MONOMER-13651
PWY-5757	1.11.1.23	(S)-2-hydroxypropylphosphonic acid epoxidase subunit
PWY-5757	1.1.1.309	MONOMER-13659
PWY-5757	2.1.1.308	MONOMER-13660
PWY-5757	4.1.1.82	MONOMER-13658
PWY-5757	5.4.2.9	phosphoenolpyruvate phosphomutase
PWY-5759	2.4.1	triterpene sapogenin carboxylic acid glucosyltransferase
PWY-5760 	1.2.1 	aldehyde dehydrogenase
PWY-5760	1.5.3.13 	non-specific polyamine oxidase
PWY-5765	2.1.1	MONOMER-13679
PWY-5765 	2.1.1 	MONOMER-17507
PWY-5765 	2.1.1 	MONOMER-17508
PWY-5765	2.1.1	phloroglucinol O-methyltransferase
PWY-5766	1.4.1.3 	glutamate dehydrogenase subunit
PWY-5766	1.4.1.3 	trimer complex
PWY-5768 	1.2.1.4 	aldehyde dehydrogenase, mitochondrial
PWY-5768 	1.2.1.4 	potassium-activated aldehyde dehydrogenase, mitochondrial
PWY-5770	1.10.3.i	dihydrophenazinedicarboxylate synthase
PWY-5770	2.5.1	PhzA/B
PWY-5770	2.6.1.86	MONOMER-13690
PWY-5770	3.3.2.15	2-amino-4-deoxychorismate hydrolase
PWY-5770	5.3.3.17	trans-2,3-dihydro-3-hydroxy-anthranilate isomerase subunit
PWY-5773	1.14.13	MONOMER-13698
PWY-5773	2.1.1	desoxyhemigossypol-6-O-methyltransferase
PWY-5773 	4.2.3.13	MONOMER-13696
PWY-5773 	4.2.3.13	MONOMER-13697
PWY-5774	2.4.1	triterpene sapogenin UDP-glucosyltransferase
PWY-5780	4.2.1 	phenolic oxidative coupling protein
PWY-5782	1.1.1	L-sorbosone dehydrogenase subunit
PWY-5782	1.1.99.21	D-sorbitol dehydrogenase large subunit 
PWY-5782	1.1.99.32	L-sorbose dehydrogenase, FAD dependent
PWY-5787	1.10.3	lacquer laccase
PWY-5787	1.10.3	MONOMER-13742
PWY-5787	1.10.3	MONOMER-13743
PWY-5787	1.10.3	MONOMER-17509
PWY-5788 	6.3	indole-3-acetyl amido synthetase
PWY-5788 	6.3	indole-3-acetyl-amido synthetase
PWY-5789	1.1.1.298	NADP-dependent malonate semialdehyde reductase
PWY-5789	1.1.1.35 	crotonyl-CoA hydratase/(S)-3-hydroxybutyryl-CoA dehydrogenase
PWY-5789	1.1.1.co	NADPH-dependent succinate semialdehyde reductase
PWY-5789	1.2.1.75 	succinyl-CoA reductase (NADPH) subunit
PWY-5789	1.2.1.76	succinyl-CoA reductase
PWY-5789	2.3.1.9	acetyl-CoA C-acetyltransferase
PWY-5789	4.2.1.116	3-hydroxypropionyl-CoA dehydratase
PWY-5789	4.2.1.120	4-hydroxybutyryl-CoA dehydratase
PWY-5789	5.1.99.1	methylmalonyl-CoA epimerase
PWY-5789	5.4.99.2	methylmalonyl-CoA mutase
PWY-5789	6.2.1.36	3-hydroxypropionyl-CoA synthetase (AMP-forming)
PWY-5789	6.2.1.40	4-hydroxybutyryl-CoA synthetase (AMP-forming)
PWY-5789	6.4.1.2 	acyl CoA carboxylase holoenzyme
PWY-5793	2.4.1	MONOMER-13720
PWY-5793	2.4.1	MONOMER-13722
PWY-5794	2.3.1.187 	malonate decarboxylase acetyl-[acp]:malonate-[acp] transferase component
PWY-5794	2.3.1 	[acyl-carrier protein] S-malonyltransferase
PWY-5794	4.1.1.89 	malonate decarboxylase
PWY-5796	2.7.7.66	MONOMER-14204
PWY-5800	2.4.1	MONOMER-13758
PWY-5800	2.4.2.24	1,4-&beta;-<small>D</small>-xylan synthase
PWY-5800	2.4.2	MONOMER-13757
PWY-5808	2.3.1.156	type III polyketide synthase
PWY-5808	2.5.1	phlorisobutyrophenone dimethylallyltransferase
PWY-581	3.5.1.4	indole-3-acetamide amidohydrolase
PWY-5813 	5.5.1.8	(+)-bornyl diphosphate synthase
PWY-5815	2.5.1 	rubber transferase
PWY-5818	1.14.11.52	validamycin A dioxygenase
PWY-5818	1.3.1	valienone 7-phosphate dehydrogenase
PWY-5818	2.4.1.338	validoxylamine A glucosyltransferase
PWY-5818	2.5.1.bp	validamine 7-phosphate valienyltransferase
PWY-5818	2.6.1.M1	validone 7-phosphate aminotransferase
PWY-5818	2.7.1.ax	C7-cyclitol kinase
PWY-5818	2.7.7.ab	valienol-1-phosphate guanylyltransferase
PWY-5818	3.1.3.101	validoxylamine A 7-phosphate phosphatase
PWY-5818	4.2.3.152	MONOMER-13772
PWY-5818	5.1.3.33	2-epi-5-epi-valiolone epimerase
PWY-5823 	1.1.1.341	MONOMER-13793
PWY-5823 	1.1.1.342	CDP-paratose synthase monomer
PWY-5823 	1.17.1.1	CDP-4-dehydro-6-deoxyglucose reductase
PWY-5823 	2.7.7.33	&alpha;-D-glucose-1-phosphate cytidylyltransferase monomer
PWY-5823 	4.2.1.45	CDP-D-glucose-4,6-dehydratase monomer
PWY-5823 	5.1.3.10	MONOMER-13795
PWY-5825	3.2.1.161	7-O--D-apiofuranosyl-(1-6)-&beta;-D-glucopyranosidase
PWY-5826	2.3.2.2	MONOMER-13775
PWY-5826	2.6.1	MONOMER-13789
PWY-5829	1.1.1.183	MONOMER-13798
PWY-5829	1.1.1.183	MONOMER-13799
PWY-5829	3.1.7.11	C-geraniol synthase
PWY-5829 	3.1.7.11 	MONOMER-12827
PWY-5829	3.1.7.11	MONOMER-12834
PWY-5829	3.1.7.11	MONOMER-13792
PWY-5835	2.3.1.84	MONOMER-13797
PWY-5835	2.3.1.84	MONOMER-13800
PWY-5835	3.1.1	MONOMER-13805
PWY-5835	3.1.1	MONOMER-17521
PWY-5835	3.1.1	MONOMER-17522
PWY-5835	3.1.1	MONOMER-17523
PWY-5835	3.1.1	MONOMER-17524
PWY-5835	3.1.1	MONOMER-17525
PWY-5835	3.1.1	MONOMER-17526
PWY-5835	3.1.1	MONOMER-17527
PWY-5837	5.4.4.2	MONOMER-13723
PWY-5838 	2.1.1.163 	bifunctional 2-octaprenyl-6-methoxy-1,4-benzoquinone methylase and S-adenosylmethionine:2-DMK methyltransferase
PWY-5838 	2.1.1.163	S-adenosylmethionine:2-demethylmenaquinol methyltransferase
PWY-5838 	2.5.1.90	solanesyl diphosphate synthase
PWY-5840 	2.1.1.163	S-adenosylmethionine:2-demethylmenaquinol methyltransferase
PWY-5840 	2.5.1.74	MONOMER-13813
PWY-5840 	2.5.1 	heptaprenyl diphosphate synthase
PWY-5840 	4.1.1.71 	MONOMER-13808
PWY-5840 	4.1.3.36	MONOMER-13812
PWY-5840 	4.2.1.113	MONOMER-13796
PWY-5840 	5.4.4.2	isochorismate synthase monomer
PWY-5840 	6.2.1.26	MONOMER-13811
PWY-5845 	2.5.1.74	MONOMER-13819
PWY-5846	1.14.21	MONOMER-13814
PWY-5846	2.1.1	MONOMER-13815
PWY-5846	2.3.1	MONOMER-13816
PWY-5846	3.5.1 	MONOMER-17529
PWY-5848	1.1.1	cinchoninone:NADPH oxidoreductase I
PWY-5848	1.1.1	cinchoninone:NADPH oxidoreductase II
PWY-5850 	2.1.1.163	MONOMER-13830
PWY-5850 	2.5.1.83	hexaprenyl diphosphate synthase
PWY-5856	1.14.99	MONOMER-13881
PWY-5856	2.5.1.39	MONOMER-13869
PWY-5859 	1.1.1.318 	anol synthase (multifunctional)
PWY-5859	1.1.1.318	MONOMER-13833
PWY-5859	1.1.1.318	MONOMER-13834
PWY-5859	1.1.1.318	MONOMER-13842
PWY-5859	1.1.1.318	MONOMER-13843
PWY-5859	1.1.1.319	MONOMER-13832
PWY-5859	1.1.1.319	MONOMER-13837
PWY-5859	2.3.1.224	MONOMER-13831
PWY-5859 	2.3.1.224	MONOMER-17540
PWY-5861 	2.2.1.9	MONOMER-13857
PWY-5861 	6.2.1.26	MONOMER-13858
PWY-5862 	2.5.1.84	solanesyl diphosphate synthase
PWY-5862 	4.1.3.36	naphthoate synthase trimer
PWY-5863 	1.6.5.12	demethylphylloquinone dehydrogenase
PWY-5863 	2.1.1.fq	2-phytyl-1,4-naphthoquinone methyltransferase
PWY-5863 	2.2.1.9	MONOMER-13721
PWY-5863 	2.2.1.9	MONOMER-13897
PWY-5863 	2.5.1.130	DHNA phytyl transferase
PWY-5863 	2.5.1.130	MONOMER-13822
PWY-5863 	3.1.2.28	DHNA-CoA thioesterase
PWY-5863 	3.1.2.28	MONOMER-15406
PWY-5863 	4.1.3.36	MONOMER-13896
PWY-5863 	4.2.1.113	MONOMER-11678
PWY-5863 	5.4.4.2	AT1G74710-MONOMER
PWY-5863 	5.4.4.2	MONOMER-2342
PWY-5863 	6.2.1.26	MONOMER-11669
PWY-5863 	6.2.1.26	MONOMER-13898
PWY-5864 	1.13.11.27	4-hydroxyphenylpyruvate dioxygenase
PWY-5864 	1.13.11.27	MONOMER-13901
PWY-5864 	2.5.1.117	MONOMER-3784
PWY-5864 	2.6.1.57 	tryptophan aminotransferase
PWY-5867 	2.1.1.279	MONOMER-13876
PWY-5867 	2.1.1.279	MONOMER-13887
PWY-5871	1.14.99	MONOMER-13880
PWY-5871	2.1.1.64 	MONOMER-13874
PWY-5871	2.5.1.39	AT4G23660-MONOMER
PWY-5872	1.14.13	ubiquinone biosynthesis monooxygenase COQ6
PWY-5872	1.14.99	ubiquinone biosynthesis protein COQ7 homolog
PWY-5872	2.1.1.201	ubiquinone biosynthesis methyltransferase COQ5
PWY-5872	2.1.1.222 	hexaprenyldihydroxybenzoate methyltransferase
PWY-5872	2.5.1.39	4-hydroxybenzoate polyprenyltransferase
PWY-5872 	2.5.1.91	decaprenyl diphosphate synthase
PWY-5874	1.14.14.18	heme oxygenase 1
PWY-5874	1.3.1.24	biliverdin reductase A
PWY-5875	1.1.1	MONOMER-17444
PWY-5875	1.14.99	diaponeurosporene oxidase
PWY-5875	1.3.8	dehydrosqualene desaturase
PWY-5875	2.3.1	MONOMER-13883
PWY-5875	2.4.1	MONOMER-13882
PWY-5875	2.5.1.96 	dehydrosqualene synthase
PWY-5876	1.14.21	MONOMER-13877
PWY-5876 	2.1.1 	coclaurine N-methyltransferase
PWY-5884	1.2.1.84	alcohol-forming fatty acyl-CoA reductase
PWY-5884	2.3.1.75	MONOMER-13890
PWY-5885	2.3.1.75 	wax synthase/diacyglycerol acyltransferase [multifunctional]
PWY-5887	4.2.3.37	MONOMER-13904
PWY-5903 	1.3.1.28	2,3-dihydro-2,3-dihydroxybenzoate dehydrogenase
PWY-5903	2.7.7.58 	bacillibactin synthase multienzyme complex
PWY-5903 	3.3.2.1	apo-DhbB pentamer
PWY-5903 	5.4.4.2	MONOMER-13807
PWY-5905	1.14.99.29	deoxyhypusine hydroxylase
PWY-5905	2.5.1.46	deoxyhypusine synthase
PWY-5905	2.5.1.46	deoxyhypusine synthase subunit
PWY-5905	2.5.1.46	MONOMER-13916
PWY-5907 	2.5.1.45 	deoxyhypusine synthase
PWY-5907	2.5.1.45	MONOMER-13925
PWY-5907	2.5.1.45	MONOMER-13926
PWY-5910 	1.1.1.34	AT1G76490-MONOMER
PWY-5910 	2.3.1.9	acetoacetyl-CoA thiolase
PWY-5910 	2.3.3.10	AT4G11820-MONOMER
PWY-5910 	2.5.1.29	AT1G49530-MONOMER
PWY-5910 	2.5.1.29	AT2G18640-MONOMER
PWY-5910 	2.5.1.29	AT2G23800-MONOMER
PWY-5910 	2.5.1.29	AT3G14550-MONOMER
PWY-5910 	2.5.1.29	AT4G36810-MONOMER
PWY-5910 	2.5.1.29	YPL069C-MONOMER
PWY-5910 	2.7.1.36	AT5G27450-MONOMER
PWY-5910 	4.1.1.33	mevalonate diphosphate decarboxylase subunit
PWY-5912	1.1.1.285	3&prime;&prime;-deamino-3&prime;&prime;-oxonicotianamine reductase
PWY-5912	2.5.1.43	nicotianamine synthase 1
PWY-5912	2.5.1.43	nicotianamine synthase 2
PWY-5912	2.5.1.43	nicotianamine synthase 3
PWY-5912	2.6.1.80	nicotianamine aminotransferase A
PWY-5912	2.6.1.80	nicotianamine aminotransferase B
PWY-5912	2.6.1.80	nicotianamine transferase 1
PWY-5915	1.3.7.2	MONOMER-13950
PWY-5915	1.3.7.3	MONOMER-13952
PWY-5917	1.14.15.20	heme oxygenase 1
PWY-5917 	1.14.15.20	MONOMER-18992
PWY-5917	1.3.7.5	phycocyanobilin:ferredoxin oxidoreductase
PWY-5918 	1.2.1.70	glutamyl-tRNA reductase
PWY-5918 	1.3.3.3	coproporphyrinogen III oxidase
PWY-5918 	1.3.3.4	protoporphyrinogen IX oxidase isozyme II
PWY-5918 	1.3.3.4 	protoporphyrinogen oxidase
PWY-5918 	1.3.3.4	protoporphyrinogen oxidase
PWY-5918 	2.5.1.61	OHMETHYLBILANESYN-MONOMER
PWY-5918 	2.5.1.61	porphobilinogen deaminase
PWY-5918 	4.1.1.37	MONOMER-11784
PWY-5918 	4.1.1.37	uroporphyrinogen decarboxylase
PWY-5918 	4.2.1.24	porphobilinogen synthase
PWY-5918 	4.2.1.75	MONOMER-11788
PWY-5918 	4.2.1.75	UROGENIIISYN-MONOMER
PWY-5918 	4.2.1.75	uroporphyrinogen III synthase
PWY-5918 	4.99.1.1	PROTOHEME-FERROCHELAT-MONOMER
PWY-5918 	5.4.3.8	glutamate-1-semialdehyde aminotransferase
PWY-5918 	6.1.1.17	glutamate-tRNA ligase
PWY-5920 	1.3.3.3	coproporphyrinogen-III oxidase
PWY-5920 	1.3.3.3	MONOMER-11601
PWY-5920 	1.3.3.4	Protoporphyrinogen oxidase
PWY-5920 	2.3.1.37	5-aminolevulinate syntase 2
PWY-5920 	2.3.1.37	5-aminolevulinate synthase 1
PWY-5920 	2.3.1.37	5-aminolevulinate synthase 2
PWY-5920 	2.3.1.37	5-aminolevulinic acid synthase 1
PWY-5920 	2.5.1.61	MONOMER-11598
PWY-5920 	2.5.1.61	porphobilinogen deaminase
PWY-5920 	4.1.1.37	MONOMER-11600
PWY-5920 	4.1.1.37	uroporphyrinogen decarboxylase
PWY-5920 	4.2.1.24	&delta;-aminolevulinate dehydratase subunit
PWY-5920 	4.2.1.24	delta-aminolevulinic acid dehydratase
PWY-5920 	4.2.1.75	MONOMER-11599
PWY-5920 	4.2.1.75	uroporphyrinogen-III synthase
PWY-5920 	4.99.1.1	ferrochelatase
PWY-5920 	4.99.1.1	MONOMER-11602
PWY-5921	6.1.1.17 	glutamyl-tRNA<sup>Glx</sup> synthetase
PWY-5921	6.1.1.24 	Glutamyl-tRNA<sup>Glx</sup> synthetase
PWY-5921	6.3.5.7	glutamyl-tRNA(Gln) amidotransferase subunit E 
PWY-5921 	6.3.5.7 	glutamyl-tRNA<sup>Gln</sup> amidotransferase subunit A 
PWY-5921	6.3.5.7	glutamyl-tRNA<sup>Gln</sup> amidotransferase subunit A 
PWY-5923 	3.3.2.8	limonene-1,2-epoxide hydrolase
PWY-5924 	1.14.13.107	limonene 1,2-monooxygenase
PWY-5925	1.14.11.24	MONOMER-13953
PWY-5925	1.14.11.25	MONOMER-13972
PWY-5925	1.14.11.25	mugineic-acid 3-dioxygenase
PWY-5927 	1.1.1.296	NAD+-dependent dihydrocarveol dehydrogenase
PWY-5927 	1.14.13.105	monocyclic monoterpene ketone monooxygenase
PWY-5927 	1.3.99.25	carvon reductase
PWY-5927 	3.1.1.83	&epsilon;-lactone hydrolase
PWY-5928	1.14.13.48	MONOMER-15424
PWY-5928	1.14.13.48	MONOMER-6743
PWY-5928 	4.2.3.16	(4S)-limonene synthase
PWY-5929	2.1.1.38	MONOMER-13985
PWY-5930	4.2.3.36	terpentetriene synthase monomer
PWY-5930	5.5.1.15	MONOMER-13987
PWY-5934	1.16.1.7	AT5G49740-MONOMER
PWY-5934	1.16.1.7	MONOMER-13991
PWY-5934	1.16.1.7	MONOMER-13992
PWY-5935	3.1.7.9 	tuberculosinol synthase
PWY-5935	5.5.1.16	G185E-7653-MONOMER
PWY-5936	2.4.1.168	xyloglucan glycosyltransferase 4
PWY-5936	2.4.1	&alpha;-1,2-fucosyltransferase
PWY-5936	2.4.1	xyloglucan &alpha;-1,2-fucosyltransferase
PWY-5936	2.4.2.39	xyloglucan 6-xylosyltransferase 1
PWY-5936	2.4.2.39	xyloglucan 6-xylosyltransferase 2
PWY-5936	2.4.2.39	xyloglucan 6-xylosyltransferase 5
PWY-5940	2.1.4.2	L-arginine:inosamine-phosphate amidinotransferase
PWY-5940	2.4.2.27	dTDP-L-dihydrostreptose:streatidine-6-phosphate dihydrostreptosyltransferase subunit
PWY-5940	2.6.1.50	L:glutamine:2-keto-myo-inositol aminotransferase subunit
PWY-5940	2.6.1.56	MONOMER-15285
PWY-5940	2.7.1.65	MONOMER-14012
PWY-5940	2.7.1.72	MONOMER-14016
PWY-5940	3.1.3.39	MONOMER-14015
PWY-5940	3.1.3.40	MONOMER-15286
PWY-5941	2.4.1.1	glycogen phosphorylase, brain isoform
PWY-5941	2.4.1.1	glycogen phosphorylase, liver isoform
PWY-5941	2.4.1.1	glycogen phosphorylase, muscle isoform
PWY-5941	2.4.1.25 	glycogen debranching enzyme
PWY-5942	1.3.99.26	all-trans-&zeta;-carotene desaturase
PWY-5942	1.3.99.29	phytoene desaturase
PWY-5942 	1.3.99.30	phytoene desaturase (3,4-didehydrolycopene-forming)
PWY-5943	5.5.1.19	chloroplastic/chromoplastic lycopene &beta; cyclase
PWY-5943	5.5.1.19	lycopene &beta;-cyclase
PWY-5947 	5.5.1 	lycopene &beta;-cyclase
PWY-5950	4.1.99.16 	MONOMER-14022
PWY-5951 	1.1.1.4 	2,3-butanediol dehydrogenase
PWY-5951	1.1.1.4	-(R,R)-2,3-butanediol dehydrogenase subunit
PWY-5957	2.5.1.43	MONOMER-14045
PWY-5957 	2.5.1.43	nicotianamine synthase
PWY-5958	2.1.1.111	MONOMER-14024
PWY-5958	2.3.1.159	MONOMER-14039
PWY-5958	4.1.3.27	anthranilate synthase
PWY-5959 	1.14.13.175	MONOMER-14044
PWY-5960 	2.1.1.109	O-methyltransferase I
PWY-5960 	2.1.1.110	O-methyltransferase II
PWY-5961 	1.1.1.349	MONOMER-14027
PWY-5961 	1.1.1.352	hydroxyaverantin dehydrogenase
PWY-5961 	1.1.1.353	versiconal hemiacetal acetate reductase
PWY-5961 	1.13.12.20	anthrone oxidase
PWY-5961 	1.14.13.174	(1S)-averantin monooxygenase
PWY-5961 	1.3.1	MONOMER-14041
PWY-5961 	2.3.1.221	MONOMER-17895
PWY-5961 	3.1.1.94	MONOMER-14036
PWY-5961 	4.2.1.142	5-oxoaverantin cyclase subunit
PWY-5961 	4.2.1.143	versiconal cyclase subunit
PWY-5963	2.8.1.9	molybdenum cofactor sulfurtransferase
PWY-5966 	2.3.1.41 	beta-keto-acyl synthase homolog
PWY-5966 	2.3.1.86 	fatty acid synthase
PWY-5966 	2.3.1.86 	fatty acid synthase, &alpha; subunit
PWY-5967	2.5.1.50	MONOMER-14061
PWY-5968	2.4.1.177	MONOMER-14070
PWY-5968	2.4.1.177	MONOMER-14073
PWY-5968	2.4.1.177	naringenin 2-hydroxylase
PWY-5970	2.3.1 	fatty acid synthase
PWY-5971 	1.3.1.9	enoyl-[acyl-carrier-protein] reductase [NADH], chloroplastic
PWY-5971 	2.3.1.86 	3-hydroxyacyl-[acp] dehydrase
PWY-5972	1.1.1	MONOMER-14115
PWY-5972	1.3.1	MONOMER-14116
PWY-5972	2.3.1.199	long-chain 3-oxoacyl-CoA synthase
PWY-5972	2.3.1.199	long-chain fatty acid elongase 1
PWY-5972	2.3.1.199 	long chain fatty acid elongase 6
PWY-5972	4.2.1.119	MONOMER-14114
PWY-5973 	1.3.1.9	enoyl-[acyl-carrier-protein] reductase
PWY-5973 	2.3.1.86 	3-oxoacyl-[acyl-carrier-protein] reductase
PWY-5973 	3.1.2.14 	multifunctional acyl-CoA thioesterase I and protease I and lysophospholipase L1
PWY-5975	1.3.1 	enone oxidoreductase
PWY-5975	2.1.1	MONOMER-14074
PWY-5975	2.4.1	MONOMER-14075
PWY-5976	3.2.1	MONOMER-14082
PWY-5976	3.2.1	MONOMER-14083
PWY-5976	4.1.2.11	MONOMER-14084
PWY-5979 	4.2.1.144 	3-amino-5-hydroxybenzoate synthase
PWY-5980	2.4.2.41	xylogalacturonan &beta;-1,3-xylosyltransferase
PWY-5981	2.3.1.51	MONOMER-14088
PWY-5981	2.3.1	MONOMER-14087
PWY-5981 	2.7.7.41	phosphatidate cytidylyltransferase
PWY-5982	1.2.1.73	sulfoacetaldehyde dehydrogenase subunit
PWY-5983	2.4.1.245	trehalose synthase subunit
PWY-5987	1.14.19.13	acyl-CoA 15-desaturase
PWY-5987	1.14.19.1	acyl-CoA 9-desaturase
PWY-5987	1.14.19.6	palmitoleoyl-CoA 12-desaturase
PWY-5987	2.1.1	MONOMER-14104
PWY-5987	2.3.1	MONOMER-15594
PWY-5987	2.3.1	MONOMER-15595
PWY-5989 	2.3.1.86 	3-oxoacyl-[acyl-carrier-protein] reductase, chloroplastic
PWY-5990 	1.14.13.117 	valine N-monooxygenase (oxime forming)
PWY-5990 	1.14.13.118 	isoleucine N-monooxygenase (oxime forming)
PWY-5990 	1.14.13.118 	valine N-monooxygenase (oxime forming)
PWY-5990 	2.4.1.63	UDP-glucose:acetone cyanohydrin &beta;-glucosyltransferase
PWY-5992	1.14.13	AT5G48000-MONOMER
PWY-5992	1.14.21	AT5G47990-MONOMER
PWY-5992	5.4.99.31	AT5G48010-MONOMER
PWY-5993 	1.1.1	UDP-3-dehydro-&alpha;-D-glucose dehydrogenase
PWY-5993 	2.1.1.315	27-O-demethyl-rifamycin SV methyltransferase subunit
PWY-5993 	2.3.1	rifamycin polyketide synthase
PWY-5993 	2.3.1	rifamycin polyketide synthase A subunit
PWY-5993 	2.5.1	MONOMER-14079
PWY-5993 	2.7.1.179	MONOMER-14080
PWY-5994 	2.3.1.86 	hydroxyacyl-thioester dehydratase
PWY-5994 	3.1.2.2	acyl-CoA thioesterase 2 (mitochondrial)
PWY-5994 	3.1.2.2 	acyl-CoA thioesterase 4, proxisomal
PWY-5995 	1.14.19.22	acyl-lipid &omega;-6 desaturase
PWY-5995 	1.14.19.23	acyl-lipid &omega;-6 desaturase (ferredoxin)
PWY-5995 	2.3.1.23	lysophospholipid acyltransferase 1
PWY-5996	1.14.19.1	acyl-CoA 9-desaturase 1
PWY-5996 	3.1.2.2	acyl-CoA thioesterase 1, cytosolic
PWY-5997 	1.14.19.25	acyl-lipid  &omega;-3 desaturase (endoplasmic reticulum)
PWY5F9-12	1.13.11.39	2,3-dihydroxybiphenyl 1,2-dioxygenase
PWY5F9-12	1.14.12.18	biphenyl 2,3-dioxygenase, ferredoxin component 
PWY5F9-12	1.14.12.18	ethylbenzene dioxygenase
PWY5F9-12	1.3.1.56	cis-2,3-dihydrobiphenyl-2,3-diol dehydrogenase
PWY5F9-12	1.3.1.56	cis-3-phenylcyclohexa-3,5-diene-1,2-diol dehydrogenase
PWY5F9-12	3.7.1.8	2-hydroxy-6-oxo-6-phenylhexa-2,4-dienoate hydrolase
PWY5F9-3233	1.14.12	phthalate 3,4-dioxygenase ferredoxin subunit 
PWY5F9-3233	1.3.1	phthalate 3,4-dihydrodiol dehydrogenase
PWY5F9-3233	4.1.1.69	3,4-dihydroxyphthalate decarboxylase
PWY-6001	1.14.19.6	oleoyl-CoA 12-desaturase
PWY-6002 	3.2.1.21	&beta;-glucosidase
PWY-6002	3.2.1.21	lotaustralin beta-glucosidase
PWY-6002	3.2.1.21	MONOMER-14132
PWY-6002	4.1.2.46	aliphatic (R)-hydroxynitrile lyase
PWY-6005	5.4.99.53	AT5G42600-MONOMER
PWY-6007	4.2.1.124	AT4G15340-MONOMER
PWY-6008	5.4.99.57	AT4G15370-MONOMER
PWY-6010	2.3.1.115	isoflavone-7-O-&beta;-glucoside 6-O-malonyltransferase
PWY-6010	2.4.1.236	MONOMER-14145
PWY-6010	2.4.1.81	flavone 7-O-&beta;-glucosyltransferase
PWY-6010	2.4.2.25	flavone apiosyltransferase
PWY-6011	3.2.1.117	MONOMER-14142
PWY-6011	3.2.1.118	MONOMER-14143
PWY-6011	4.1.2.10	MONOMER-14144
PWY-601 	1	CYP83B1 monooxygenase
PWY-6012-1	2.7.8.7	L-aminoadipate-semialdehyde dehydrogenase-phosphopantetheinyl transferase
PWY-6012	3.1.4.14	EG11095-MONOMER
PWY-6013	1.14.19.39	linoleate 12-acetylenase
PWY-6014	1.14.18	MONOMER-14156
PWY-6015	2.1.1.153	MONOMER-14161
PWY-6015	2.4.1.105	MONOMER-14160
PWY-6019	2.7.1.83	MONOMER-14165
PWY-6019	4.2.1.70	MONOMER-14164
PWY-6021	1.5.1.42 	nitrilotriacetate monooxygenase
PWY-6024	2.4.1.106	MONOMER-14169
PWY-6024	2.4.1	isovitexin 2\-O-arabinosyltransferase"
PWY-6024	2.4.1	isovitexin 7-O-galactosyltransferase
PWY-6024	2.4.1	MONOMER-14168
PWY-6024	2.4.1	MONOMER-14170
PWY-6024	2.4.1	MONOMER-14172
PWY-6027	1.1.1.195 	MONOMER-14183
PWY-6027	1.2.1.44	MONOMER-14182
PWY-6028	2.3.1.190	acetoin dehydrogenase complex
PWY-6029 	1.17.1.3	leucoanthocyanidin reductase
PWY-6030	1.14.16.4	tryptophan 5-hydroxylase 1
PWY-6030	1.14.16.4	tryptophan hydroxylase 2
PWY-6030	2.1.1.4	acetylserotonin O-methyltransferase
PWY-6030	2.1.1.4	acetylserotonin O-methyltransferase
PWY-6030	2.3.1.5 	serotonin N-acetyltransferase
PWY-6030	2.3.1.87	aralkylamine N-acetyltransferase
PWY-6030 	4.1.1.28	aromatic L-amino acid decarboxylase
PWY-6030	4.1.1.28	L-DOPA decarboxylase
PWY-6032	1.1.1.278	MONOMER-14295
PWY-6032	1.3.1.22	MONOMER-14269
PWY-6032	1.3.1	MONOMER-14261
PWY-6035	1.3.1.77	anthocyanidin reductase
PWY-6035	1.3.1.77	MONOMER-14199
PWY-6038	4.1.3.6 	citrate lyase
PWY-6038	4.1.3.6 	citrate lyase &alpha; subunit hexamer
PWY-6038	4.1.3.6 	citrate lyase &beta; subunit hexamer
PWY-6038	4.1.3.6 	citrate lyase complex
PWY-6039	1.14.13.36	MONOMER-14223
PWY-6039	1.14.13.36	MONOMER-14259
PWY-6039	1.14.13.36	p-coumaroyl ester 3-hydroxylase
PWY-6039	2.1.1.104	MONOMER-14260
PWY-6039	2.3.1.133	MONOMER-14258
PWY-6039	2.3.1 	hydroxycinnamoyl-coenzyme A shikimate/quinate hydroxycinnamoyl transferase
PWY-6040	1.14.13.36	MONOMER-14201
PWY-6040	1.14.13.36	MONOMER-14219
PWY-6040 	2.3.1.133	MONOMER-15052
PWY-6040 	2.3.1.133	MONOMER-15053
PWY-6041	1.13.11.3	protocatechuate 3,4-dioxygenase type II
PWY-6041	1.3.1.32	MONOMER-14207
PWY-6041	3.1.1.92	MONOMER-14206
PWY-6041	5.5.1.2	3-carboxymuconate lactonizing enzyme type 2
PWY-6044	1.14.13.111	methanesulfonate monooxygenase hydroxylase component 
PWY-6046 	4.4.1.3	dimethylsulfoniopropionate lyase
PWY-6046	4.4.1.3	dimethylsulfoniopropionate lyase
PWY-6046 	4.4.1.3	DMSP lyase monomer
PWY-6046	4.4.1.3	MONOMER-14221
PWY-6046	4.4.1.3	MONOMER-14241
PWY-6047	1.8.3.4	MONOMER-14227
PWY-6049 	1.3.8	MONOMER-16784
PWY-6049 	2.1.1.269	dimethylsulfoniopropionate demthylase
PWY-6049 	4.2.1.155	methylthioacryloyl-CoA hydrolase
PWY-6049 	6.2.1.44	3-methylmercaptopropionyl-CoA ligase
PWY-6051 	2.4.1.91 	flavonol 3-O-glucosyltransferase
PWY-6051 	2.4.1 	UDP glucose:cytokinin glycosyltransferase
PWY-6051	2.4.1	UDP glucose:dinitrotoluene glycosyltransferase
PWY-6051 	2.4.1	UDP-glucose:flavonol-3-O-glycoside-7-O-glucosyltransferase
PWY-6051	2.4.1	UDP-indole-3-butyric acid glucosyltransferase
PWY-6054	1.2.1.3	MONOMER-14236
PWY-6054	2.1.1.12	methionine S-methyltransferase subunit
PWY-6056	2.3.1	dimethylsulfoniopropanoate:acyl-CoA transferase
PWY-6057	1.8.2.4	dimethylsulfide dehydrogenase
PWY-6059	1.14.13.M46	dimethylsulfoxide monooxygenase
PWY-6059	1.14.14.ah	methanesulfonate monooxygenase
PWY-6059	1.14.14.ai	dimethylsulfone monooxygenase
PWY-6060	2.3.1.187 	MONOMER-14254
PWY-6060	4.1.1.89 	malonyl-S-ACP:biotin-protein carboxyltransferase &alpha; subunit 
PWY-6060	6.2.1.35 	MONOMER-14253
PWY-6061	1.1.1.181	3&beta;-hydroxy-&Delta;<sup>5</sup>-C27 steroid oxidoreductase
PWY-6061	1.14.14.23	cytochrome P450 7A1
PWY-6061	1.14.18.8	cytochrome P450 8B1
PWY-6061	1.17.99	peroxisomal acyl-coenzyme A oxidase 2
PWY-6061	1.3.1.3	&Delta;<sup>4</sup>-3-oxosteroid 5&beta;-reductase
PWY-6061	2.3.1.176	sterol carrier protein 2
PWY-6061	2.3.1.65	bile acid-CoA:amino acid N-acyltransferase
PWY-6061	4.2.1.107 	peroxisomal multifunctional enzyme type 2
PWY-6061	5.1.99.4	&alpha;-methylacyl-CoA racemase
PWY-6061	6.2.1.28 	bile acid CoA ligase
PWY-6061 	6.2.1.7 	very long-chain acyl-CoA synthetase
PWY-6064 	2.1.1.42	flavonol 3-O-methyltransferase
PWY-6068	2.4.1.220	MONOMER-14292
PWY-6071 	1.1.1.35	3-hydroxyadipyl-CoA dehydrogenase (NAD<sup>+</sup>)
PWY-6071 	1.14.13.149	ring 1,2-phenylacetyl-CoA epoxidase
PWY-6071 	1.2.1.39	phenylacetaldehyde dehydrogenase
PWY-6071 	1.4.3.21	copper-containing amine oxidase
PWY-6071 	2.3.1.223 	&beta;-ketoadipyl-CoA thiolase
PWY-6071 	3.3.2.12 	oxepin-CoA hydrolase/3-oxo-5,6-dehydrosuberyl-CoA semialdehyde dehydrogenase
PWY-6072 	1.14.13.131	NADH-dependent dimethylsulfide monooxygenase
PWY-6072 	1.8.1.17	MONOMER-14331
PWY-6073	1.1.1.132	GDP-mannose 6-dehydrogenase
PWY-6073	1.1.1.132	MONOMER-14358
PWY-6073	2.4.1.33	MONOMER-14359
PWY-6073	5.1.3.37	mannuronan C-5 epimerase
PWY-6074 	1.1.1.170	NAD(P)-dependent 3&beta;-hydroxy-4&alpha;-carboxy-sterol 3-dehydrogenase (decarboxylating)
PWY-6074 	1.3.1.70	lamin B receptor
PWY-6076	1.14.14.24	cytochrome P450 2R1
PWY-6076	1.14.15.18	cytochrome p450 27B1
PWY-6077	1.14.13.40	2-amninobenzoyl-CoA monooxygenase/reductase subunit
PWY-6077	6.2.1.32	aerobic 2-aminobenzoate-CoA ligase type 2
PWY-6077	6.2.1.32	aerobic 2-aminobenzoate-CoA ligase type I
PWY-6079	1.14.12.1	anthranilate dioxygenase reductase component 
PWY-6080	1.1.1	MONOMER-14370
PWY-6080	1.1.1	MONOMER-14379
PWY-6080	1.17.99	p-ethylphenol methylhydroxylase subunit 1 
PWY-6080	1.3.7.9	4-hydroxybenzoyl-CoA reductase, &alpha; subunit 
PWY-6081	1.14.12	MONOMER-14388
PWY-6082	1.1.1.132	GDP-mannose 6-dehydrogenase
PWY-6082	2.4.1.33	mannuronosyl transferase
PWY-6082	5.1.3.37	mannuronan C5 epimerase
PWY-6083 	1.14.12	chlorobenzene dioxygenase
PWY-6084	5.2.1.10	MONOMER-14432
PWY-6086 	1.14.11	2,4-dichlorophenoxyacetate dioxygenase
PWY-6086 	1.14.13.7 	2,4-dichlorophenol 6-monooxygenase I
PWY-6086 	1.14.13.7 	2,4-dichlorophenol 6-monooxygenase II
PWY-6087 	3.1.1.45	dienelactone hydrolase
PWY-6087 	3.1.1.45	dienelactone hydrolase I
PWY-6087	5.5.1.7	chloromuconate cycloisomerase II
PWY-6088 	1.3.1 	1,6-dihydroxycyclohexa-2,4-diene-1-carboxylate dehydrogenase
PWY-6089 	1.3.1.32	chloromaleylacetate reductase
PWY-6089 	1.3.1.32	maleylacetate reductase
PWY-6089 	1.3.1.32	maleylacetate reductase 1
PWY-6089 	3.1.1.45	dienelactone hydrolase
PWY-6089 	5.5.1.7 	chloromuconate cycloisomerase
PWY-6091 	1.14.12	chlorobenzene dioxygenase
PWY-6094 	5.5.1 	chloromuconate cycloisomerase
PWY-6095	5.4.99.37	MONOMER-14418
PWY-6098	4.2.1.129 	hydroxyhopane synthase
PWY-6098	5.4.99.8	MONOMER-14421
PWY-6100	1.14.11.1	&gamma;-butyrobetaine dioxygenase
PWY-6100	1.14.11.1	&gamma;-butyrobetaine hydroxylase subunit
PWY-6100	1.14.11.8	&epsilon;-N-trimethyllysine hydroxylase subunit
PWY-6100	1.14.11.8	trimethyllysine dioxygenase
PWY-6100	1.2.1.47	MONOMER-14430
PWY-6102 	5.5.1.7 	chloromuconate cycloisomerase I
PWY-6104 	1.1.1.90	benzyl alcohol dehydrogenase
PWY-6104 	1.14.13.dw	xylene methyl-monooxygenase
PWY-6104 	1.2.1.28	benzaldehyde dehydrogenase
PWY-6105	1.3.1.96	squalene synthase-like 2
PWY-6105	1.3.1.97	squalene synthase-like 3
PWY-6105	2.5.1.103 	MONOMER-14640
PWY-6105	2.5.1.103	Squalene synthase-like 1
PWY-6107	1.3.1.32	MONOMER-14438
PWY-6107	3.1.1.45	trans-dienelactone hydrolase
PWY-6109 	5.4.99.39	MONOMER-14451
PWY-6109 	5.4.99.39 	triterpene synthase
PWY-6109 	5.4.99.41	MONOMER-14452
PWY-6109 	5.4.99.41 	triterpene synthase
PWY-6109	5.4.99.8	MONOMER-14445
PWY-6109	5.4.99.8	MONOMER-14449
PWY-6111	2.3.1.21	carnitine O-palmitoyltransferase 1, brain isoform
PWY-6111	2.3.1.21	carnitine O-palmitoyltransferase 1, liver isoform
PWY-6111	2.3.1.21	carnitine O-palmitoyltransferase 1, muscle isoform
PWY-6111	2.3.1.21	carnitine O-palmitoyltransferase 2
PWY-6111	2.3.1.21	carnitine palmitoyltransferase 1a
PWY-6111	2.3.1.21	carnitine palmitoyltransferase 2
PWY-6111	2.3.1.7	carnitine O-acetyltransferase
PWY-6111	2.3.1.7	MONOMER-14440
PWY-6113 	1.1.1.M9 	meromycolic acid 3-oxoacyl-[acyl-carrier-protein] reductase
PWY-6113 	1.1 	G185E-4279-MONOMER
PWY-6113 	1.1	Rv2509
PWY-6113 	1.3.1.9 	fatty acid synthase
PWY-6113 	1.3.1.M4	meromycolic acid enoyl-[acyl-carrier-protein] reductase
PWY-6113 	2.1.1.79 	methoxy mycolic acid synthase 2
PWY-6113 	2.1.1	cyclopropane mycolic acid synthase 2
PWY-6113 	2.1.1	cyclopropane synthase
PWY-6113 	2.1.1	methoxy mycolic acid synthase 1
PWY-6113 	2.1.1	methoxy mycolic acid synthase 3
PWY-6113 	2.1.1	methoxymycolic acid synthase 4
PWY-6113	2.3.1.86 	acyl-carrier-protein S-malonyltransferase
PWY-6113 	2.3.1.M1	meromycolic acid 3-oxoacyl-(acyl carrier protein) synthase I
PWY-6113 	2.3.1.M2	meromycolic acid 3-oxoacyl-(acyl carrier protein) synthase II
PWY-6113 	2.3.1	Pks13
PWY-6113 	3.1.3.12	mycolyl-trehalose-6-phosphate phosphatase
PWY-6113 	3.1.3.12	trehalose-6-phosphate phosphatase OtsB1
PWY-6113 	4.2.1.M1 	meromycolic acid 3-hydroxyacyl-[acyl-carrier-protein] dehydratase I
PWY-6113 	4.2.1.M2 	meromycolic acid 3-hydroxyacyl-[acyl-carrier-protein] dehydratase II
PWY-6113 	5.3.3.14	2-trans-enoyl-ACP isomerase
PWY-6113 	6.2.1	acyl-CoA synthetase
PWY-6113	6.4.1.2	accD6
PWY-6113	6.4.1.2	acetyl/propionyl-CoA carboxylase alpha chain
PWY-6113 	6.4.1.3	propionyl-CoA carboxylase beta chain 4
PWY-6113 	6.4.1.3	propionyl-coa carboxylase beta chain 5
PWY-6115	5.4.99.8	MONOMER-14473
PWY-6116	2.4.1.246	MONOMER-14460
PWY-6116	3.1.3.79	MONOMER-14461
PWY-6117	1.5.3.13 	N<sup>1</sup>-acetylpolyamine oxidase (3-acetamamidopropanal -forming)
PWY-6117	1.5.3.17 	polyamine oxidase 1 (spermidine-forming)
PWY-6117	2.3.1.57	spermidine/spermine N<sup>1</sup>-acetyltransferase 1
PWY-6118	1.1.1.8	AT2G41540-MONOMER
PWY-6118	1.1.1.8	cytoplasmic glycerol-3-phosphate dehydrogenase (NAD+)
PWY-6118 	1.1.5.3	AT3G10370-MONOMER
PWY-6118	1.1.5.3	mitochondrial glycerol-3-phosphate dehydrogenase
PWY-6120	1.3.1.43	MONOMER-8134
PWY-6122	2.1.2	phosphoribosylglycinamide formyltransferase 2
PWY-6124 	4.3.2.2	adenylosuccinate lyase
PWY-6125 	1.17.4.1	ribonucleoside-diphosphate reductase 2
PWY-6126 	1.1.98.f	anaerobic ribonucleoside-triphosphate reductase
PWY-6126 	4.3.2.2	adenylosuccinate lyase
PWY-6128	4.2.3.67	cis-muuroladiene synthase
PWY-6129	1.3.1.94	polyprenol reductase
PWY-6129	2.5.1.87	dehydrodolichyl diphosphate syntase complex
PWY-6129	2.5.1.87	dehydrodolichyl diphosphate synthase complex
PWY-6129	2.7.1.108	dolichol kinase
PWY-6131	2.7.1.29	dihydroxyacetone kinase 2
PWY-6132	5.4.99.7	MONOMER-14484
PWY-6134 	1.14.16.1	phenylalanine hydroxylase subunit
PWY-6139	2.5.1.56	N-acetylneuraminate synthase subunit
PWY-6139	3.2.1.183	MONOMER-14542
PWY-6140	2.5.1.132	MONOMER-14547
PWY-6140	2.7.7.ac	MONOMER-14549
PWY-6140	3.1.3.aa	MONOMER-14548
PWY-6143	2.3.1.202	HP0327-MONOMER
PWY-6143	2.5.1.97	pseudaminic acid synthase
PWY-6143	2.6.1.92	PseC monomer
PWY-6143	2.7.7.81	HP0326-MONOMER
PWY-6143	3.6.1.57	UDP-2,4-diacetamido-2,4,6-trideoxy-&beta;-<small>L</small>-altropyranose hydrolase
PWY-6145 	1.6.2.2 	cytidine monophosphate-N-acetylneuraminate hydroxylase system
PWY-6145 	2.5.1.132	MONOMER-14550
PWY-6145 	2.5.1.56	N-acetylneuraminate synthase subunit
PWY-6145 	2.5.1.57	N-acetylneuraminate 9-phosphate synthase subunit
PWY-6145 	2.5.1.57	MONOMER-14514
PWY-6145 	2.5.1.57	sialic acid synthase
PWY-6145 	2.5.1	7-O-acetyl-N-acetylneuraminate synthase subunit
PWY-6145 	2.7.1.60 	bifunctional UDP-N-acetylglucosamine 2-epimerase/N-acetylmannosamine kinase
PWY-6145 	2.7.7.43	cytidine 5-monophosphate N-acetylneuraminate synthetase subunit
PWY-6145 	2.7.7.43	MONOMER-14544
PWY-6145 	2.7.7.ac	MONOMER-14551
PWY-6145 	3.1.3.29	N-acetylneuraminate-9-phosphatase
PWY-6145 	3.1.3.29	MONOMER-14517
PWY-6145 	3.2.1.183	MONOMER-14541
PWY-6146 	1.2.7.1	pyruvate synthase subunit &alpha; 
PWY-6146 	2.7.8.38	archaetidylserine synthase
PWY-6146 	4.1.1.31	phosphoenolpyruvate carboxylase monomer
PWY-6146 	4.2.1.11	MONOMER-14533
PWY-6146 	4.2.1.1	carbonic anhydrase subunit
PWY-6146 	5.3.1.1	triosephosphate isomerase monomer
PWY-6146 	5.4.2.12	2,3-bisphosphoglycerate-independent phosphoglycerate mutase 1
PWY-6146 	5.4.2.12	2,3-bisphosphoglycerate-independent phosphoglycerate mutase 2
PWY-6146 	6.4.1.1	pyruvate carboxylase subunit A 
PWY-6147 	2.7.6.3 	dihydropterin pyrophosphokinase/dihydropteroate  synthase
PWY-6148	1.5.99.15	dihydromethanopterin reductase (acceptor)
PWY-6148	2.4.2.54	&beta;-ribofuranosylhydroxybenzene 5-phosphate synthase
PWY-6148	2.5.1.105	MONOMER-14570
PWY-6148 	2.7.6.3	6-hydroxymethyl-7,8-dihydropterin pyrophosphokinase
PWY-6148	6.3.1	4-(&beta;-D-ribofuranosyl)hydroxybenzene 5-phosphate--L-aspartate ligase
PWY-6151	4.4.1.21	S-ribosylhomocysteine lyase monomer
PWY-6153 	3.2.2.9 	5-methylthioadenosine/S-adenosylhomocysteine nucleosidase
PWY-6153	4.4.1.21	S-ribosylhomocysteine lyase
PWY-6153	4.4.1.21	S-ribosylhomocysteine lyase monomer
PWY-6154 	3.2.2.9 	5-methylthioadenosine/S-adenosylhomocysteine nucleosidase
PWY-6154	4.4.1.21	S-ribosylhomocysteine lyase monomer
PWY-6157	2.3.1.184	acyl-homoserine-lactone synthase
PWY-6158	2.7.3.2	creatine kinase BB isoform
PWY-6158	2.7.3.2	creatine kinase MB isoform
PWY-6158	2.7.3.2	creatine kinase MM isoform
PWY-6158	2.7.3.2	sarcomeric mitochondrial creatine kinase
PWY-6158	2.7.3.2	ubiquitous mitochondrial creatine kinase
PWY-6164 	4.2.3.4 	pentafunctional AROM polypeptide
PWY-6165 	1.1.1.25	shikimate dehydrogenase
PWY-6165 	1.4.1.24	dehydroquinate synthase II
PWY-6165 	2.2.1.11 	6-deoxy-5-ketofructose 1-phosphate synthase
PWY-6165 	2.7.1.71	MONOMER-18801
PWY-6165 	4.1.2.13 	D-fructose 1,6-bisphosphate aldolase
PWY-6165	4.2.1.10	MONOMER-14591
PWY-6165 	5.3.1.1 	triosephosphate isomerase
PWY-6167	1.1.1.302	MONOMER-14596
PWY-6167	2.5.1.78	6,7-dimethyl-8-ribityllumazine synthase monomer
PWY-6167	2.5.1.9	riboflavin synthase monomer
PWY-6167	2.7.1.161	MONOMER-13864
PWY-6167	2.7.7.2	FAD synthase
PWY-6167	3.5.1.102	2-amino-5-formylamino-6-ribosylaminopyrimidin-4(3H)-one 5-monophosphate deformylase
PWY-6167	3.5.4.29	GTP cyclohydrolase III
PWY-6167	4.1.99.12	3,4-dihydroxy-2-butanone 4-phosphate synthase monomer
PWY-6168	1.1.1.302	MONOMER3O-27
PWY-6168	2.5.1.78	6,7-dimethyl-8-ribityllumazine synthase monomer
PWY-6168	2.5.1.9	riboflavine synthetase monomer
PWY-6168	2.7.1.26	riboflavin kinase
PWY-6168	2.7.7.2	FAD synthetase
PWY-6168	3.5.4.25	GTP cyclohydrolase II
PWY-6168	3.5.4	DRAP deaminase
PWY-6168	4.1.99.12	MONOMER-14611
PWY-6173	4.1.1.22	histidine decarboxylase
PWY-6173	4.1.1.22	L-histidine decarboxylase C-terminally truncated subunit
PWY-6173	4.1.1.22	pyridoxal 5-phosphate-dependent histidine decarboxylase
PWY-6173	4.1.1.22	pyruvoyl-dependent histidine decarboxylase
PWY-6174	1.1.1.34	MONOMER-14624
PWY-6174	2.3.1.9	acetyl-CoA C-acetyltransferase monomer
PWY-6174 	2.7.1.36	mevalonate kinase
PWY-6174	2.7.1.36	mevalonate kinase monomer
PWY-6174	2.7.4.26	MONOMER-14619
PWY-6174	4.1.1.99	MONOMER-18695
PWY-6174	5.3.3.2	isopentenyl-diphosphate &delta;-isomerase
PWY-6176 	1.8.3	MONOMER-13526
PWY-6176 	4.4.1.4	MONOMER-13501
PWY-6176	4.4.1.4	MONOMER-13504
PWY-6176 	5.3.99.M1	sulfenic acid isomerase
PWY-6178	1.14.14	REUT_A1585-MONOMER
PWY-6181	1.4.3.22	amiloride-sensitive amine oxidase
PWY-6181	1.4.3.22	diamine oxidase
PWY-6181	2.1.1.8	histamine N-methyltransferase
PWY-6181	2.1.1.8	MONOMER-14646
PWY-6182 	1.13.11.1	MONOMER-14642
PWY-6182 	5.3.3.4	MONOMER-14645
PWY-6182 	5.5.1.1	MONOMER-14643
PWY-6184 	1.14.13.1	salicylate 1-hydroxylase
PWY-6185 	1.13.11.1	MONOMER-14435
PWY-6185	5.3.3.4	methylmuconolactone isomerase monomer
PWY-6185	5.4.99.14	4-methylmuconolactone methylisomerase monomer
PWY-6185 	5.5.1.7 	muconate cycloisomerase I
PWY-6192 	1.13.11	chlorocatechol 1,2-dioxygenase
PWY-6192 	1.14.12 	chlorobenzene dioxygenase
PWY-6192 	1.3.1	chlorobenzene dihydrodiol dehydrogenase
PWY-6192 	3.1.1	dienelactone hydrolase
PWY-6192 	5.5.1.7 	chloromuconate cycloisomerase
PWY-6193	1.3.1.32	maleylacetate reductase 1
PWY-6193	3.1.1.45	MONOMER-14667
PWY-6193	5.5.1	5-chloromuconolactone dehalogenase
PWY-6193	5.5.1.7	chloromuconate cycloisomerase 2
PWY-6196	5.1.1.18	serine racemase
PWY-6200	1.13.11.37	1,2,4-trihydroxybenzene 1,2-dioxygenase
PWY-6200	1.14.14	chlorophenol-4-monooxygenase
PWY-6200	1.3.1.32	MONOMER-14680
PWY-6200	1.6.5.7	hydroxybenzoquinone reductase
PWY-6200	4.5.1	MONOMER-14675
PWY-6210	1.13.11.74	2-aminophenol 1,6-dioxygenase &alpha; subunit 
PWY-6210	1.13.11.74	2-aminophenol-1,6-dioxygenase &alpha; subunit 
PWY-6210	1.2.1.32	2-aminomucoate semialdehyde dehydrogenase subunit
PWY-6210	1.2.1.32	2-aminomuconic 6-semialdehyde dehydrogenase monomer
PWY-6210	3.5.99.5	2-aminomuconate deaminase monomer
PWY-6210	3.5.99.5	2-aminomuconate deaminase subunit
PWY-6210	4.1.1.77	4-oxalocrotonate decarboxylase monomer
PWY-621	2.7.1.1 	hexokinase
PWY-621	3.2.1.26 	invertase
PWY-621	3.2.1.26 	MONOMER-19181
PWY-621	3.2.1.26 	neutral/alkaline invertase
PWY-6215	1.14.13.2	MONOMER-14755
PWY-6215	3.1.2.23	4-hydroxybenzoyl-CoA thioesterase monomer
PWY-6215	3.8.1.7	4-chlorobenzoate-coenzyme A dehalogenase monomer
PWY-6215	6.2.1.33	4-chlorobenzoate:coenzyme A ligase monomer
PWY-6216 	1.13.11.8 	protocatechuate 4,5-dioxygenase
PWY-6217 	1.14.12	3-chlorobenzoate-3,4-dioxygenase
PWY-6221	1.14.12.13	2-halobenzoate dioxygenase multicomponent enzyme system
PWY-622	2.4.1.18	1,4-alpha-glucan branching enzyme
PWY-622	2.4.1.18	starch branching enzyme I
PWY-622	2.4.1.18	starch branching enzyme II
PWY-622	2.4.1.21	soluble starch synthase I
PWY-622	2.4.1.21	soluble starch synthase III
PWY-622	2.4.1.21	starch synthase I
PWY-622	2.4.1.21	starch synthase II
PWY-622	2.4.1.21	starch synthase III
PWY-622	2.4.1.242	granule-bound starch synthase I
PWY-622	2.7.7.27	ADP-glucose pyrophosphorylase
PWY-622	2.7.7.27	ADP-glucose pyrophosphorylase large subunit 
PWY-622	2.7.7.27	ADP-glucose pyrophosphorylase small subunit 
PWY-6223	1.13.11.4	gentisate 1,2-dioxygenase
PWY-622	3.2.1.68	isoamylase
PWY-622	3.2.1.68	isoamylase catalytic subunit
PWY-6223	3.7.1.20	fumarylpyruvate hydrolase
PWY-6223	5.2.1.4	glutathion-dependent maleylpyruvate isomerase
PWY-6223	5.2.1.4	mycothiol-dependent maleylpyruvate isomerase
PWY-6224	1.14.13.172	salicylate-5-hydroxylase large subunit 
PWY-622	5.3.1.9	MONOMER-12904
PWY-6232 	2.1.1.42	phenylpropanoid/flavonoid 3-O-methyltransferase
PWY-6234 	6.3	jasmonoyl-amino acid synthetase
PWY-6235	2.8.2	hydroxyjasmonate sulfotransferase
PWY-6239	2.4.1.236	MONOMER-14785
PWY-6239	2.4.1.81	MONOMER-14784
PWY-6241	1.11.1.8	thyroid peroxidase
PWY-6241	3.4.22.15	cathepsin L
PWY-6241	3.4.22.1	cathepsin B
PWY-6241	3.4.22.38	cathepsin K
PWY-6241	3.4.23.5	cathepsin D
PWY-6243	4.2.3.81	trans-&alpha;-bergamotene synthase
PWY-6243	4.2.3 	sesquiterpene synthase
PWY-6244	2.5.1.92	Z,Z-farnesyl pyrophosphate synthase
PWY-6254 	4.2.3.50 	santalene and bergamontene synthase
PWY-6257	4.2.3.94	MONOMER-14858
PWY-6258	4.2.3 	MONOMER-14859
PWY-6258	4.2.3 	patchoulol synthase subunit
PWY-6260 	1.21.99.4 	type I iododthyronine deiodinase
PWY-6261	2.4.1.17	MONOMER-14907
PWY-6261	2.4.1.17	MONOMER-14909
PWY-6261	2.4.1.17	MONOMER-14910
PWY-6261	2.4.1.17	MONOMER-14911
PWY-6261	2.4.1.17	UDP-glucuronosyltransferase 1-8
PWY-6261	2.4.1.17	UDP-glucuronosyltransferase 1A3
PWY-6261	2.8.2.1	sulfotransferase 1A1
PWY-6261 	2.8.2.1	sulfotransferase 1A2
PWY-6265	1.1.1.326	zerumbone synthase
PWY-6265	1.14.13.150	&alpha;-humulene 10-hydroxylase
PWY-6265	4.2.3.57 	alpha-humulene synthase
PWY-6266 	1.14.11.23	flavonol synthase
PWY-6266 	1.14.13.21	MONOMER-12072
PWY-6266 	2.1.1.267	myricetin O-methyltransferase
PWY-6266 	2.1.1.42 	AT5G54160-MONOMER
PWY-6266 	2.1.1.76	MONOMER-14285
PWY-6266 	2.1.1 	caffeoyl-CoA O-methyltransferase
PWY-6266 	2.3.1.173	MONOMER-12519
PWY-6266 	2.4.1.159	MONOMER-12683
PWY-6266 	2.8.2.25	AT3G45070-MONOMER
PWY-6266 	2.8.2.25	MONOMER-14676
PWY-6266 	2.8.2.25	MONOMER-14678
PWY-6266 	2.8.2.26	MONOMER-14702
PWY-6266 	2.8.2.27	MONOMER-14703
PWY-6268 	2.5.1.17	CobA
PWY-6268 	2.5.1.17	MONOMER-14887
PWY-6269	2.7.7.62	MONOMER-14884
PWY-6269	2.7.8.26	MONOMER-14886
PWY-6269	3.1.3.73	MONOMER-14882
PWY-6269	3.5.1.90	adenosylcobinamide amidohydrolase
PWY-6269	3.5.1.90 	MONOMER-14883
PWY-6269	6.3.1.10	MONOMER-14885
PWY-6270 	4.2.3.27	isoprene synthase
PWY-6271	4.2.3.84 	&beta;-eudesmol synthase
PWY-6	2.7.1.52	MONOMER-11586
PWY-6273 	2.7.8.29	phosphatidylserine synthase 2
PWY-6275	4.2.3.57	MONOMER-14906
PWY-6275	4.2.3.57	MONOMER-14912
PWY-6275	4.2.3.57	MONOMER-14914
PWY-6275	4.2.3.57	MONOMER-14917
PWY-6275	4.2.3.57	MONOMER-14918
PWY-6275	4.2.3.57 	MONOMER-16135
PWY-6275	4.2.3 	(-)-&beta;-caryophyllene synthase
PWY-6275	4.2.3 	(E)-beta-caryophyllene/beta-elemene synthase
PWY-6277 	2.1.2.2	phosphoribosylglycinamide formyltransferase 1
PWY-6277 	2.1.2 	phosphoribosylglycinamide formyltransferase 2
PWY-6277 	2.4.2.14	amidophosphoribosyl transferase
PWY-6	2.7.7.30	MONOMER-11585
PWY-6277 	6.3.3.1	phosphoribosylformylglycinamide cyclo-ligase
PWY-6277 	6.3.4.13	GLYCRIBONUCSYN-MONOMER
PWY-6277 	6.3.5.3	FGAMSYN-MONOMER
PWY-6278	4.2.3.74	MONOMER-14922
PWY-6279	1.14.13	beta-carotene hydroxylase
PWY-6279	4.2.1.131	MONOMER-14923
PWY-6279	5.5.1	1-hydroxylycopene cyclase
PWY-6281	2.7.1.164	MONOMER-14956
PWY-6281	2.7.1.164	MONOMER-14957
PWY-6281	2.7.9.3	MONOMER-14953
PWY-6281	2.7.9.3	selenophosphate synthetase 2
PWY-6281	2.9.1.2	O-phosphoseryl-tRNA:selenocysteinyl-tRNA synthase dimer
PWY-6281	2.9.1.2	selenocysteine synthase
PWY-6281	6.1.1.11	seryl-tRNA synthetase
PWY-6284 	2.3.1.179 	&beta;-ketoacyl-ACP synthase II
PWY-6285 	2.3.1.180 	&beta;-ketoacyl-ACP synthase III
PWY-6285 	6.2.1.3 	fatty acyl-CoA synthetase
PWY-6285 	6.2.1.3 	short chain acyl-CoA synthetase
PWY-6286	1.14.15.9	spheroidene monooxygenase
PWY-6286	1.14.99	hydroxyneurosporene desaturase
PWY-6286 	2.1.1.210 	1-hydroxy-carotenoid methylase
PWY-6286	2.1.1.210 	hydroxyneurosporene methyltransferase
PWY-6286 	4.2.1.131	acyclic carotenoid 1,2-hydratase
PWY-6286	4.2.1.131	MONOMER-14930
PWY-6287	1.3.99.28	phytoene desaturase
PWY-6288	2.4.1.276	zeaxanthin glucosyltransferase
PWY-6289	2.7.7.M9	3,4-dihydroxybenzoate-AMP ligase
PWY-6289	4.2.1.118	3-dehydroshikimate dehydratase
PWY-6290	4.2.3.128	MONOMER-14947
PWY-6291	4.2.3.73	MONOMER-14951
PWY-6291	4.2.3.86 	(+)-valencene synthase
PWY-6291	4.2.3 	selina-4,11-diene / intermedeol synthase
PWY-6292 	2.5.1.6	methionine adenosyltransferase I, &alpha; subunit
PWY-6292 	3.3.1.1	S-adenosylhomocysteine hydrolase subunit
PWY-6292 	4.2.1.22	cystathionine beta-synthase
PWY-6292 	4.2.1.22	cystathionine &beta;-synthase subunit
PWY-6292 	4.4.1.1	cystathionine &gamma;-lyase subunit
PWY-6294	4.2.3.66	MONOMER-14960
PWY-6294	4.2.3.66 	selinene synthase
PWY-6297	2.4.1	MONOMER-14966
PWY-6297	2.4.1	MONOMER-14967
PWY-6303	2.1.1.278	753205-MONOMER
PWY-6303	2.1.1.278	AT5G55250-MONOMER
PWY-6303	3.1.1.1	methyl indole-3-acetate methylesterase
PWY-6304	4.2.3.8	MONOMER-14979
PWY-6304	4.2.3.8	MONOMER-14980
PWY-6305	3.5.3.1	MONOMER-14987
PWY-6305	3.5.3.1	MONOMER-14988
PWY-6305	4.1.1.17	MONOMER-14989
PWY-6305	4.1.1.17	MONOMER-14990
PWY-6305	4.1.1.17	MONOMER-14991
PWY-6305	4.1.1.19	arginine decarboxylase
PWY-6305	4.1.1.19	MONOMER-14981
PWY-6305	4.1.1.19	MONOMER-14983
PWY-6305	4.1.1.19	MONOMER-14984
PWY-6305	4.1.1.19	MONOMER-14985
PWY-6307	1.1.1.2	aldehyde reductase
PWY-6307	4.1.1.28	aromatic L-amino acid decarboxylase
PWY-6308	2.5.1.73	Sep-tRNA:Cys-tRNA synthase
PWY-6308	6.1.1.27	O-phosphoseryl-tRNA ligase monomer
PWY-6309 	1.2.1.32	MONOMER-15010
PWY-6309 	2.6.1.7 	kynurenine/alpha-aminoadipate aminotransferase
PWY-6309 	2.6.1.7 	kynurenine--oxoglutarate transaminase 3
PWY-6310	2.3.1	aloesone synthase
PWY-6312	2.3.1	MONOMER-15019
PWY-6312	2.3.1	octaketide synthase
PWY-6312	2.3.1	type III polyketide synthase
PWY-6313 	2.8.2.1	sulfotransferase 1A3/1A4
PWY-6314	2.3.1	&alpha; pyrone polyketide synthase
PWY-6316	2.3.1.74	MONOMER-15028
PWY-6316	2.3.1 	aromatic polyketide synthase
PWY-6317	2.7.1.6	GALACTOKIN-MONOMER
PWY-6317	2.7.7.12	galactose 1-phosphate uridyl transferase
PWY-6317	2.7.7.12	galactose-1-phosphate uridylyltransferase
PWY-6317	5.1.3.3	galactose-1-epimerase
PWY-6317	5.1.3.3	MONOMER-6141
PWY-6318 	1.13.11 	4-hydroxyphenylpyruvate dioxygenase
PWY-6318 	1.4.3.4	monoamine oxidase B
PWY-6318	2.3.1.14	MONOMER-15067
PWY-6318	2.3.1.192	MONOMER-15068
PWY-6318	2.6.1.1	mitochondrial aspartate aminotransferase
PWY-6318	4.1.1.53	aromatic L-amino acid decarboxylase
PWY-6321 	1.2.1.16 	NAD(P)+-dependent succinate semialdehyde dehydrogenase
PWY-6321 	2.6.1 	4-aminobutyrate aminotransferase
PWY-6322	1.1.1.309	phophonoacetaldehyde reductase
PWY-6322	1.1.1	hydroxymethylphosphonate dehydrogenase
PWY-6322	1.13.11.72	hydroxyethylphosphonate dioxygenase
PWY-6322	1.2.1	phosphonoformaldehyde dehydrogenase
PWY-6322	2.1.1.326	N-acetyl-demethylphophinothricin P-methyltransferase
PWY-6322	2.3.1.183	demethyl-phosphinothricin N-acetyltransferase
PWY-6322	2.3.1	N-acetyl demethylphosphinothricinyl-transferase
PWY-6322	2.3.1 	nonribosomal peptide synthetase PhsA
PWY-6322	2.3.1 	nonribosomal peptide synthetase PhsB
PWY-6322	2.3.1	nonribosomal peptide synthetase PhsC
PWY-6322	2.3.3.b	2-phosphinomethylmalate synthase
PWY-6322	2.7.1	CMP-5-phosphonoformate--3-phosphoglycerate phosphonoformyl transferase
PWY-6322	2.7.7.ad	phosphonoformate--CTP cytidylyltransferase
PWY-6322	2.7.8.23	carboxyphosphonoenolpyruvate phosphonomutase
PWY-6322	3.1.1	N-acetyl-phophinothricin tripeptide acetyl hydrolase
PWY-6322	3.1.2	phosphinothricin tripeptide thioesterase
PWY-6322	4.1.1.82	phosphonopyruvate decarboxylase
PWY-6322	4.2.1.166	phosphinomethylmalate isomerase
PWY-6322	4.2.1	2-phophono-formylglycerate enolase
PWY-6322	5.4.2.9	phosphoenolpyruvate phosphomutase
PWY-6323	2.3.1.144	MONOMER-15060
PWY-6323	2.3.1.144	MONOMER-15061
PWY-6323	2.3.1.144	MONOMER-15062
PWY-6324	1.13.12.17	flavin-dependent monooxygenase RebC
PWY-6324	1.13.12.17	P450 oxygenase RebP
PWY-6324	1.14.19.9	tryptophan 7-halogenase
PWY-6324	1.21.98.2 	dichlorochromopyrrolate synthase
PWY-6324	1.4.3.23	7-chloro-L-tryptophan oxidase
PWY-6324	2.1.1.164	demethylrebeccamycin--D-glucose O-methyltransferase
PWY-6324	4.1.99 	RebG N-glycosyl transferase
PWY-6325	1.14.13.87	MONOMER-15064
PWY-6325	2.1.1.65	MONOMER-15065
PWY-6326	4.1.1.28	MONOMER-15092
PWY-6326	4.1.1.28	MONOMER-15093
PWY-6326	4.1.1.28	MONOMER-15094
PWY-6326	4.3.3.2	MONOMER-15079
PWY-6326	4.3.3.2	MONOMER-15080
PWY-6328	1.2.1	glutarate semialdehyde dehydrogenase
PWY-6328	1.5.1	MONOMER-15078
PWY-6328	2.6.1.48	5-aminovalerate aminotransferase
PWY-6328	2.6.1	L-lysine aminotransferase
PWY-6330	4.1.1.1	MONOMER-15096
PWY-6333 	1.1.1.1 	MONOMER-11088
PWY-6337	1.21.3.3	MONOMER-15118
PWY-6337 	1.21.3 	tetrahydroprotoberberine synthase
PWY-6338 	1.1.1.312	LigC
PWY-6338 	1.1.1.38 	LigK
PWY-6338 	1.2.1.67	vanillin dehydrogenase
PWY-6338 	3.1.1.57	2-pyrone-4,6-dicarboxylate hydrolase
PWY-6338 	4.2.1.83	4-oxalomesaconate hydratase monomer
PWY-6339 	1.13.11.57	gallate dioxygenase
PWY-6339	2.1.1	syringate O-demethylase
PWY-6339 	2.1.1	vanillate/3-O-methylgallate O-demethylase
PWY-6339 	4.2.1.83	4-oxalmesaconate hydratase monomer
PWY-6339 	5.3.2.8	MONOMER-16014
PWY-6340	1.14.13	5,5-dehydrodivanillate O-demethylase
PWY-6340	4.1.1	5-carboxyvanillate decarboxylase
PWY-6340	4.1.1	MONOMER-15125
PWY-6342 	1.1.1.1 	alcohol dehydrogenase 2
PWY-6342	1.1.1.1	alcohol dehydrogenase 4
PWY-6342 	1.2.1.3 	mitochondrial aldehyde dehydrogenase
PWY-6342 	1.4.3.4	monoamine oxidase A
PWY-6342 	2.1.1.28	phenylethanolamine N-methyltransferase
PWY-6342 	2.1.1.6	catechol O-methyltransferase
PWY-6343	4.2.1.101 	feruloyl-CoA hydratase/lyase
PWY-6343	6.2.1.12 	feruloyl-CoA synthetase
PWY-6345	1.4.3	L-amino acid oxidase
PWY-6346	1.4.3	L-amino acid oxidase
PWY-6346	2.1.1.139	MONOMER-15145
PWY-6346	2.1.1	MONOMER-15146
PWY-6346	2.4.1	MONOMER-15148
PWY-6348	1.11.1.7 	MONOMER-15160
PWY-6348	1.11.1.7 	MONOMER-15162
PWY-6348	3.1.3.2 	MONOMER-15153
PWY-6348	3.1.3.2 	MONOMER-15154
PWY-6348	3.1.3.2 	MONOMER-15155
PWY-6348	3.1.3.2 	purple acid phosphatase
PWY-6348	3.1.3.2 	purple acid phosphatase [multifunctional]
PWY-6349	1.3.7.11	2,3-bis-O-geranylgeranyl-sn-glycero-phospholipid reductase
PWY-6350 	1.1.1.261	(NAD(P)-dependent glycerol-1-phosphate dehydrogenase
PWY-6350 	1.3.7.11	2,3-bis-O-geranylgeranyl-sn-glycero-phospholipid reductase
PWY-6350 	2.5.1.41	phosphoglycerol geranylgeranyltransferase
PWY-6350 	2.5.1.42	digeranylgeranylglyceryl phosphate synthase
PWY-6350 	2.7.7.67	CDP-archaeol synthase
PWY-6350	2.7.8.39 	archaetidylinositol phosphate synthase
PWY-6351 	2.7.1.149 	phosphatidylinositol-5-phosphate 4-kinase
PWY-6351 	2.7.1.67	phosphatidylinositol 4 kinase alpha
PWY-6351 	2.7.1.67	phosphatidylinositol 4-kinase beta
PWY-6351 	2.7.1.67	phosphatidylinositol 4-kinase type 2-alpha
PWY-6351 	2.7.1.67	phosphatidylinositol 4-kinase type-2 beta
PWY-6351 	2.7.1.68	phosphatidylinositol-4-phosphate 5-kinase type 1 alpha
PWY-6351 	2.7.1.68	phosphatidylinositol-4-phosphate 5-kinase, type 1 beta
PWY-6351 	2.7.1.68	phosphatidylinositol-4-phosphate 5-kinase, type 1 gamma
PWY-6351 	2.7.8.11	CDP-diacylglycerol--inositol 3-phosphatidyltransferase
PWY-6351 	2.7.8.11	phosphatidylinositol synthase
PWY-6351 	3.1.4.11	phosphoinositide phospholipase C 1
PWY-6351 	3.1.4.11	phosphoinositide phospholipase C 2
PWY-6351 	3.1.4.11	phospholipase C&beta;2
PWY-6351 	3.1.4.11	phospholipase C&delta;1
PWY-6351 	3.1.4.11	phospholipase C &epsilon;1
PWY-6351 	3.1.4.11	phospholipase C&gamma;2
PWY-63	5.1.3.5	UDP-xylose 4-epimerase subunit
PWY-6352	2.7.1.137	phosphoinositide-3-kinase, class 3
PWY-6352	2.7.1.150	FYVE finger-containing phosphoinositide kinase
PWY-6352	2.7.1.153	phosphatidylinositol 3-kinase class IA p110&alpha;/p85&alpha;
PWY-6352	2.7.1.153	phosphatidylinositol 3-kinase class IA p110&beta;/p85&gamma;
PWY-6352	2.7.1.153	phosphatidylinositol 3-kinase class IA p110&delta;/p85&beta;
PWY-6352	2.7.1.153	phosphatidylinositol 3-kinase, class IB, p110&gamma;/p101
PWY-6352	2.7.1.153	phosphatidylinositol 3-kinase, class IB, p110&gamma;/p87
PWY-6352	2.7.1.154	phosphatidylinositol 3-kinase C2 domain-containing beta polypeptide
PWY-6352	2.7.1.154	phosphatidylinositol 3-kinase C2 domain-containing gamma polypeptide
PWY-6352	2.7.1.154	phosphoinositide-3-kinase, class 2, alpha polypeptide
PWY-6352 	3.1.3 	phosphatidylinositide phosphatase SAC1
PWY-6352	3.1.3	phosphatidylinositol 3,5-bisphosphate 5-phosphatase
PWY-6353 	1.1.1.205	IMP dehydrogenase
PWY-6353 	1.1.1.205	inosine-5-monophosphate dehydrogenase 1
PWY-6353 	1.1.1.205	inosine-5-monophosphate dehydrogenase 2
PWY-6353 	1.17.1.4 	xanthine oxidase
PWY-6353 	1.17.1.4 	xanthine oxidoreductase
PWY-6353	2.4.2.15 	MONOMER-14469
PWY-6353 	3.1.3.1 	alkaline phosphatase
PWY-6353 	3.1.3.5	cytosolic 5-nucleotidase 1A
PWY-6353 	3.1.3.5 	pyrimidine nucleotidase
PWY-6353 	3.1.3.99 	purine nucleotidase
PWY-6353 	3.5.4.3	G7502-MONOMER
PWY-6353 	3.5.4.3	guanine deaminase
PWY-6358 	2.7.1.127 	inositol polyphosphate multikinase
PWY-6358 	2.7.1.127	inositol-trisphosphate 3-kinase A
PWY-6358 	2.7.1.127	inositol-trisphosphate 3-kinase B
PWY-6358 	2.7.1.127	inositol-trisphosphate 3-kinase C
PWY-6358 	3.1.3.25	inositol monophosphatase 3
PWY-6358 	3.1.3.57	inositol polyphosphate 1-phosphatase
PWY-6358 	3.1.3.62 	multiple inositol polyphosphate phosphatase 1
PWY-6362 	2.7.1.158	inositol-1,3,4,5,6-pentakisphosphate 2-kinase
PWY-6363 	3.1.3.36 	inositol polyphosphate 5-phosphatase OCRL1
PWY-6363	3.1.3.36 	phosphatidylinositol/ inositol polyphosphate 5-phosphatase
PWY-6363	3.1.3.56	inositol-polyphosphate 5-phosphatase
PWY-6363 	3.1.3.56 	inositol polyphosphate 5-phosphatase K
PWY-6363 	3.1.3.56 	phosphatidylinositol-4,5-bisphosphate 5-phosphatase, A
PWY-6363 	3.1.3.56 	synaptojanin 1
PWY-6363 	3.1.3.56 	synaptojanin 2
PWY-6363 	3.1.3.56 	type II inositol-1,4,5-trisphosphate 5-phosphatase
PWY-6363 	3.1.3.56	type I inositol-1,4,5-trisphosphate 5-phosphatase
PWY-6364 	3.1.3.56	type I inositol-1,4,5-trisphosphate 5-phosphatase 1
PWY-6365 	2.7.1.140 	inositol 1,3,4,6-tetrakisphosphate 5-kinase
PWY-6366 	2.7.1.159 	inositol-polyphosphate kinase/phosphatase
PWY-6367 	3.1.3.78	type II phosphatidylinositol 4,5-bisphosphate 4-phosphatase
PWY-6367 	3.1.3.78	type I phosphatidylinositol 4,5-bisphosphate 4-phosphatase
PWY-6367 	3.1.3.95 	myotubularin-related protein 14
PWY-6367 	3.1.3.95 	myotubularin-related protein 3
PWY-6368	3.1.3.36 	phosphatidylinositide phosphatase SAC2
PWY-6368	3.1.3.36 	type IV inositol polyphosphate 5-phosphatase
PWY-6368	3.1.3.66	type II inositol-polyphosphate 4-phosphatase
PWY-6368	3.1.3.66	type I inositol-polyphosphate 4-phosphatase
PWY-6368 	3.1.3.67 	phosphatidylinositol-3,4,5-trisphosphate 3-phosphatase PTEN
PWY-6368 	3.1.3.86 	phosphatidylinositol-3,4,5-trisphosphate 5-phosphatase 1
PWY-6368 	3.1.3.86 	phosphatidylinositol-3,4,5-trisphosphate 5-phosphatase 2
PWY-6368	3.1.3	phosphoinositide lipid phosphatase PLIP
PWY-6369	2.7.4.21 	inositol hexakisphosphate and diphosphoinositol-pentakisphosphate kinase 1
PWY-6369	2.7.4.21 	inositol hexakisphosphate/diphosphoinositol-pentakisphosphate kinase 2
PWY-6369	2.7.4.21 	inositol hexakisphosphate kinase 1
PWY-6369	2.7.4.21 	inositol hexakisphosphate kinase 2
PWY-6369	2.7.4.21 	inositol hexakisphosphate kinase 3
PWY-6370 	1.1.1.213	3&alpha;-hydroxysteroid dehydrogenase (A-specific)
PWY-6370	1.6.5.4	MONOMER-15193
PWY-6370	1.8.1.9 	NADPH thioredoxin reductase
PWY-6370	1.8.5.1	glutaredoxin-1
PWY-6370	1.8.5.1 	protein disulfide isomerase
PWY-6370 	2.5.1.18 	glutathione transferase, &Omega; class
PWY-6371 	2.7.1.150	1-phosphatidylinositol-3-phosphate 5-kinase
PWY-6371 	2.7.1.158	inositol (1,3,4,5,6)-pentakisphosphate 2-kinase
PWY-6371 	2.7.8.11	YPR113W-MONOMER
PWY-6371 	3.1.3	phosphatidylinositol 3,5-bisphosphate 5-phosphatase
PWY-6371 	3.1.3 	phosphoinositide phosphatase
PWY-6372	2.7.1.140	inositol-1,3,4,6-tetrakis phosphate 5-kinase
PWY-6372	2.7.1.64	myo-inositol 3-kinase
PWY-6372	5.5.1.4	MONOMER-15198
PWY-6373	1.1.1.59	MONOMER-15202
PWY-6373	1.2.1.18 	malonate-semialdehyde dehydrogenase (acetylating)
PWY-6373	2.3.1	AcuN
PWY-6374 	1.3.1.28	vibriobactin-specific 2,3-dihydro-2,3-dihydroxybenzoate dehydrogenase
PWY-6374	2.7.7.58	vibriobactin-specific 2,3-dihydroxybenzoate-AMP ligase
PWY-6374 	3.3.2.1	holo-[VibB protein dimer]
PWY-6374 	5.4.4.2	vibriobactin-specific isochorismate synthase
PWY-6374	6.3.2	Vibrio cholerae apo-VibF protein
PWY-6375	4.1.1.18	SCO2782-MONOMER
PWY-6375	6.3.3	desferrioxamine E synthetase
PWY-6376 	4.1.1.18	MONOMER-15208
PWY-6377	1.14.13.30	cytochrome P450 4F2
PWY-6381	4.1.1.18	MONOMER-15217
PWY-6385	1.3.1.98	G185E-4610-MONOMER
PWY-6385	2.4.1.129 	MONOMER-15251
PWY-6385	2.4.1.227	G185E-6361-MONOMER
PWY-6385	2.5.1.7	G185E-5489-MONOMER
PWY-6385	2.7.8	G185E-6364-MONOMER
PWY-6385	5.1.1.3	G185E-5517-MONOMER
PWY-6385	6.3.2.10	G185E-6365-MONOMER
PWY-6385	6.3.2.13	G185E-6366-MONOMER
PWY-6385	6.3.2.4	G185E-7236-MONOMER
PWY-6385	6.3.2.8	G185E-6360-MONOMER
PWY-6385	6.3.2.9	G185E-6363-MONOMER
PWY-6389 	1.1.1.304 	acetoin reductase subunit
PWY-6389 	1.1.1.304 	L-glycol dehydrogenase
PWY-6390 	1.1.1.76 	diacetyl (acetoin) reductase
PWY-6391	1.1.1.304 	meso-butanediol dehydrogenase [(R)-acetoin-forming]
PWY-6395 	2.1.1.280 	MONOMER-15249
PWY-6395 	2.1.1.280	MONOMER-15257
PWY-6395 	2.3.1.30	serine acetyltransferase
PWY-6395 	2.5.1.47	MONOMER-16911
PWY-6395 	4.4.1.8	MONOMER-16912
PWY-6401	1.11.2	MONOMER-15266
PWY-6401	1.11.2	MONOMER-15267
PWY-6401	1.11.2	quinol vinyl ether peroxidase
PWY-6402 	1.11.2.2	myeloperoxidase
PWY-6402 	1.14.14.1	cytochrome P450 1A2
PWY-6402 	1.14.14.1	cytochrome P450 1B1
PWY-6402 	1.14.14.1	cytochrome P450 2C19
PWY-6402 	1.4.3.4	monoamine oxidase A
PWY-6403	2.5.1	D-galactose-2,6-sulfurylase I
PWY-6403	2.5.1	D-galactose-2,6-sulfurylase II
PWY-6404 	1.1.1.333	G185E-8087-MONOMER
PWY-6404 	1.1.98.3	decaprenylphospho-&beta;-D-ribofuranose 2-dehydrogenase
PWY-6404 	2.3.1 	antigen 85A
PWY-6404 	2.3.1 	antigen 85B
PWY-6404 	2.3.1 	antigen 85C
PWY-6404 	2.4.1.288	arabinogalactan UDP-galactofuranosyltransferase 2
PWY-6404 	2.4.1.289	G185E-7539-MONOMER
PWY-6404 	2.4.1	arabinogalactan galactofuranosyl transferase 1
PWY-6404 	2.4.2.45	DPPR synthase
PWY-6404 	2.4.2.46	G185E-8088-MONOMER
PWY-6404 	2.4.2.47	G185E-6921-MONOMER
PWY-6404 	2.5.1.68	Z-farnesyl diphosphate synthase
PWY-6404 	2.5.1.86	Z-decaprenyl diphosphate synthase
PWY-6404 	2.7.8.35	G185E-5476-MONOMER
PWY-6404 	3.6.1	G185E-6342-MONOMER
PWY-6404 	5.1.3.2	G185E-7913-MONOMER
PWY-6404 	5.4.99.9	G185E-8105-MONOMER
PWY-6405 	5.4.2.4 	phosphoglycerate mutase 1
PWY-6405 	5.4.2.4 	phosphoglycerate mutase 2
PWY-6406	4.2.99.21	MONOMER-15305
PWY-6406	5.4.4.2	MONOMER-15306
PWY-6406	5.4.4.2 	salicylate synthase
PWY-6407	2.7.7 	non-ribosomal peptide synthetase HMWP2 
PWY-6407	2.7.7 	salicyl-AMP ligase
PWY-6408	2.3.1 	PchE
PWY-6408	2.7.7	non-ribosomal peptide synthetase PchD
PWY-6409	1.14.13.195	MONOMER-15307
PWY-6409	2.6.1.76	L-2,4-diaminobutyrate:2-ketoglutarate 4-aminotransferase monomer
PWY-6409	6.3.2	non-ribosomal peptide synthetase PvdD 
PWY-641	1.14.11.19	anthocyanidin synthase
PWY-641 	1.17.1.3	leucoanthocyanidin reductase
PWY-641	1.3.1.77	anthocyanidin reductase 1
PWY-641	1.3.1.77	anthocyanidin reductase 2
PWY-6411	3.2.1.193	3-O-&beta;-D-(1&rarr;2)-glucose-hydrolyzing ginsenosidase
PWY-6411	3.2.1.195	20-O-&beta;-D-(1&rarr;6)-glucose-hydrolyzing ginsenosidase
PWY-6412	3.2.1.195	ginsenoside Rb1 &beta;- glucosidase
PWY-641	2.4.1.M12	epicatechin-specific glucosyltransferase
PWY-6413	3.2.1.193	MONOMER-15322
PWY-6415	1.1.1.22	MONOMER-15323
PWY-6415	1.1.1.365	MONOMER-15325
PWY-6415	3.1.1	aldonolactonase
PWY-6416	1.1.1.24	MONOMER-15333
PWY-6416	4.2.1.10	3-dehydroquinate dehydratase monomer
PWY-6416	4.2.1.10	MONOMER-15335
PWY-6416	4.2.1.118	MONOMER-15329
PWY-6416 	4.2.1.118	MONOMER-15334
PWY-6418	2.3.1.208	4-hydroxycoumarin synthase
PWY-6419 	1.1.1 	quinate/shikimate dehydrogenase
PWY-6420	1.3.3.11	pyrroloquinoline-quinone synthase monomer
PWY-6421	1.20.4.3	mycoredoxin 1
PWY-6421	1.8.1.15	mycothione reductase monomer
PWY-6421	2.8.4.2	mycothiol-dependent arsenate reductase
PWY-6422	1.4.1	anabolic L-arginine dehydrogenase
PWY-6423	3.4.23.38	plasmepsin-1
PWY-6423	3.4.23.39	plasmepsin-2
PWY-6423	4.99.1.8	heme ligase monomer
PWY-6426	1.17.99.4	uracil/thymine oxidase
PWY-6426	3.5.1.95	MONOMER-15397
PWY-6426	3.5.2.1	barbiturase
PWY-6430 	1.3.1.2	(NADP)-dependent dihydropyrimidine dehydrogenase
PWY-6430 	3.5.1.6	MONOMER-15400
PWY-6430 	3.5.2.2	dihydropyrimidinase
PWY-6431	1.2.1.64	MONOMER-15419
PWY-6431	3.7.1	MONOMER-15418
PWY-6432	2.3.1.217 	curcumin synthase
PWY-6432	2.3.1.218	diketide-CoA synthase
PWY-6433	1.14.19.25	acyl-lipid &omega;-3 desaturase
PWY-6433	2.3.1.199	3-ketoacyl-CoA synthase
PWY-6435 	6.2.1.12	MONOMER-13687
PWY-6435	6.2.1.12	MONOMER-15421
PWY-6436	1.14.13.48 	(-)-limonene-7-hydroxylase
PWY-6437 	4.2.3.10	MONOMER-15442
PWY-6437 	4.2.3.10	MONOMER-15444
PWY-6438	2.3.1.218	p-coumaroyl-diketide-CoA synthase
PWY-6440	1.5.3.14	polyamine oxidase (propane-1,3-diamine forming)
PWY-6441	1.5.3.16 	spermine oxidase
PWY-6441	1.5.3.17 	AT1G65840-MONOMER
PWY-6441	1.5.3.17 	non-specific polyamine oxidase
PWY-6442	1.14.99	CYP98A8
PWY-6443	2.3.1.196	benzoyl-CoA:benzyl alcohol/phenylethanol benzoyltransferase
PWY-6443 	6.2.1 	cinnamate:CoA ligase
PWY-6444	1.2.1.28 	AT5G04580-MONOMER
PWY-6444	1.2.1.3 	benzaldehyde dehydrogenase
PWY-6444 	4.3.1.25 	phenylalanine ammonia lyase
PWY-6445	4.2.3.119 	monoterpene synthase
PWY-6446 	4.2.1	MONOMER-15509
PWY-6446 	6.2.1	cinnamate:CoA ligase
PWY-6447	2.5.1.28	MONOMER-15449
PWY-6447	4.2.3.51	&beta;-phellandrene synthase 1
PWY-6448	2.3.1.64	MONOMER-15451
PWY-6451	4.2.3.107	MONOMER-15470
PWY-6453	1.14.13.204	acyl-CoA &omega;-hydroxylase
PWY-6454	1.1.1.28	MONOMER-15467
PWY-6454	3.4.13.22	VanX D-alanyl-D-alanine dipeptidase
PWY-6454	6.1.2.1	VanA D-alanine--D-lactate ligase
PWY-6454	6.1.2.1	VanB D-alanine--D-lactate ligase
PWY-6454	6.1.2.1	VanD D-alanine--D-lactate ligase
PWY-6455	3.4.16.4 	VanXY D,D-dipeptidase/D,D-carboxypeptidase
PWY-6455	5.1.1.1 	serine/alanine racemase
PWY-6455	6.3.2.35	VanC1 D-alanine--D-serine ligase
PWY-6455	6.3.2.35	VanC2/3 D-alanine--D-serine ligase
PWY-6455	6.3.2.35	VanC4 D-alanine--D-serine ligase
PWY-6455	6.3.2.35	VanE D-alanine--D-serine ligase
PWY-6455	6.3.2.35	VanG D-Alanine--D-Serine ligase
PWY-6456	2.6.1	MONOMER-15485
PWY-6456	3.1.3	MONOMER-15486
PWY-6457	4.3.1.25 	phenylalanine ammonia lyase
PWY-6457	6.2.1	cinnamate:CoA ligase
PWY-6458	4.2.1	cinnamoyl-CoA hydratase
PWY-6462	2.3.2.16 	N-acetylmuramoyl-N-glucosaminyl-diphosphoundecaprenyl-pentapeptide:alanine alanyltransferase
PWY-6462	2.3.2	N-acetylmuramoyl-N-glucosaminyl-diphosphoundecaprenyl-hexapeptide:alanine alanyltransferase
PWY-6463	2.3.2.10	UDP-N-acetylmuramoylpentapeptide-lysine N<sup>6</sup>-alanyltransferase
PWY-6464	1.1.2.6	polyvinyl alcohol dehydrogenase (cytochrome)
PWY-6466	4.3.3.6 	pyridoxal 5-phosphate synthase complex
PWY-6467	2.4.99.14 	3-deoxy-D-manno-2-octulosonate transferase
PWY-6470	3.4.16	MONOMER-15534
PWY-6470	3.4.17.14	DdcY D-alanyl-D-alanine carboxypeptidase
PWY-6471	3.4.16.4	PBP5
PWY-6471 	6.3.1.12	MONOMER-15487
PWY-6471 	6.3.1	peptidoglycan asparagine synthase
PWY-6475	1.3.5.5	15-cis phytoene desaturase
PWY-6475	1.3.5.6	&zeta;-carotene desaturase
PWY-6475	5.2.1.12	15-cis-&zeta;-carotene isomerase
PWY-6476	2.7.7.76	molybdenum cofactor cytidylyltransferase
PWY-6477	2.1.1 	gibberellin carboxyl methyltransferase
PWY-6478	2.7.1.168	MONOMER-15575
PWY-6478	2.7.7.71	MONOMER-15576
PWY-6478	3.1.3.83	MONOMER-15574
PWY-6478	5.3.1.28	MONOMER-15573
PWY-6482	2.1.1.98	diphthine synthase
PWY-6482	2.5.1.108	2-(3-amino-3-carboxypropyl)histidine synthase
PWY-6483 	2.7.1.91	sphingoid base kinase
PWY-6483	3.5.1.23	MONOMER-15591
PWY-6493	1.1.1.332	chanoclavine-I dehydrogenase
PWY-6493	2.1.1.261	4-dimethylallyltryptophan N-methyltransferase
PWY-6493	2.5.1.34	4-dimethylallyltryptophan synthase
PWY-6494	1.14.13	gibberellin 16&alpha;, 17-epoxidase
PWY-6495	1.14.13	elymoclavine monooxygenase
PWY-6498 	1.14.18.1 	tyrosinase
PWY-6498	5.3.3.12	L-dopachrome tautomerase
PWY-6501 	1.1.1.203	D-uronate dehydrogenase
PWY-6501 	1.1.1.203	uronate dehydrogenase
PWY-6501 	5.5.1.27	D-galactarolactone cycloisomerase
PWY-6502 	3.6.1.55 	8-oxo-dGTP pyrophosphohydrolase
PWY-6502	3.6.1.55 	HS02879-MONOMER
PWY-6503 	1.1.1.332	chanoclavine-I dehydrogenase
PWY-6503 	1.14.13	elmoclavine monooxygenase
PWY-6503 	1.5.1.46	agroclavine dehydrogenase
PWY-6503 	2.1.1.261	MONOMER-17450
PWY-6503 	2.5.1.34	tryptophan dimethylallyltransferase
PWY-6505	1.13.11.6	MONOMER-15634
PWY-6505 	1.14.14.8	anthranilate hydroxylase
PWY-6505	4.1.1.45	MONOMER-15635
PWY-6507	1.1.1.127	2-dehydro-3-deoxy-<small>D</small>-gluconate 5-dehydrogenase
PWY-6507	5.3.1.17	5-dehydro-4-deoxy-<small>D</small>-glucuronate isomerase
PWY-6509	1.1.99.36 	NDMA-dependent methanol dehydrogenase
PWY-6510	1.1.1.244	NAD-dependent methanol dehydrogenase
PWY-6511	2.1.1.243	2-ketoarginine:SAM methyltransferase
PWY-6511	2.1.1.243	MONOMER-18913
PWY-6511	2.6.1	MONOMER-15642
PWY-6511	2.6.1	MONOMER-18912
PWY-6512	1.12.1.3	NADP-dependent hydrogenase
PWY-6515	1.3.1	MONOMER-15673
PWY-6515	2.3.1	MONOMER-15677
PWY-6516 	1.1.1.127 	predicted 2-keto-3-deoxy-D-gluconate dehydrogenase
PWY-6516 	1.1.1.2 	aldehyde reductase, NADPH-dependent
PWY-6516 	1.1.1.2 	NADP<sup>+</sup>-dependent aldehyde reductase
PWY-6516 	1.1.1.365	MONOMER-15601
PWY-6516 	1.1.1.365	MONOMER-15602
PWY-6516 	1.1.1.372 	MONOMER-15607
PWY-6516 	1.2.1.26	&alpha;-ketoglutarate semialdehyde dehydrogenase (D-glucarate-induced isozyme)
PWY-6516 	1.2.1.26	&alpha;-ketoglutarate semialdehyde dehydrogenase subunit
PWY-6516 	2.7.1.178 	MONOMER-15644
PWY-6516 	4.1.2.54	MONOMER-15605
PWY-6516 	4.1.2.55 	MONOMER-15645
PWY-6516 	4.2.1.146	MONOMER-15603
PWY-6516 	4.2.1.146	MONOMER-15604
PWY-6516 	4.2.1.40	D-glucarate dehydratase subunit
PWY-6516 	4.2.1.41	5-dehydro-4-deoxyglucarate dehydratase subunit
PWY-6516 	4.2.1.41	MONOMER-15617
PWY-6516 	4.2.1.42	MONOMER-15616
PWY-6516 	5.3.1.17	4-deoxy-L-threo-5-hexosulose-uronate ketol-isomerase
PWY-6516 	5.3.1.17	5-dehydro-4-deoxy-<small>D</small>-glucuronate isomerase
PWY-6516 	5.3.1.17	MONOMER-15648
PWY-6516 	5.4.1.4	D-galactarolactone isomerase
PWY-6517	2.7.1.59	N-acetyl-D-glucosamine kinase
PWY-6517	2.7.1.59	MONOMER-19002
PWY-6517 	3.5.1.25	N-acetylglucosamine-6-phosphate deacetylase
PWY-6517 	3.5.99.6	glucosamine-6-phosphate deaminase
PWY-6518	1.1.1.159	7&alpha;-hydroxysteroid dehydrogenase
PWY-6518	1.1.1.176	12&alpha;-hydroxysteroid dehydrogenase
PWY-6518	1.1.1.201	7&beta;-hydroxysteroid dehydrogenase
PWY-6518	1.1.1.238	12&beta;-hydroxysteroid dehydrogenase
PWY-6518	1.1.1.391	NAD-dependent bile acid 3&beta;-dehydrogenase
PWY-6518	1.1.1.52	3&alpha;-hydroxysteroid dehydrogenase
PWY-6518 	3.5.1.24	choloylglycine hydrolase
PWY-6518	3.5.1.24	choloylglycine hydrolase
PWY-6520	2.5.1.85	AT1G17050-MONOMER
PWY-6520	2.5.1.85	AT1G78510-MONOMER
PWY-6524 	2.4.1.166	MONOMER-15687
PWY-6525	2.4.1	MONOMER-15688
PWY-6526	1.1.1.144	perillyl alcohol dehydrogenase
PWY-6526	1.14.13	limonene hydroxylase
PWY-6526	1.14.13	perillyl aldehyde monooxygenase
PWY-6527 	2.7.1.6	MONOMER-15709
PWY-6527	2.7.7.12	MONOMER-15710
PWY-6527 	2.7.7.64 	MONOMER-9002
PWY-6527 	2.7.7.64 	UDP-galactose pyrophosphorylase
PWY-6527	3.2.1.22	MONOMER-15706
PWY-6527	3.2.1.22	MONOMER-15708
PWY-6527	3.2.1.22	MONOMER-16141
PWY-6527	3.2.1	MONOMER-16136
PWY-6527	3.2.1	MONOMER-16140
PWY-6529	1.13.11.49	chlorite O2-lyase
PWY-6529	1.97.1.1	chlorate reductase
PWY-6530	1.13.11.49	chlorite O2-lyase
PWY-6530	1.97.1	perchlorate reductase
PWY-6531	1.1.1.17	mannitol-1-phosphate 5-dehydrogenase
PWY-6531	1.1.1.67	MONOMER-15725
PWY-6531	2.7.1.4	hexokinase
PWY-6531	3.1.3.22	24323504
PWY-6531	3.1.3.22	MONOMER-15724
PWY-6533	1.14.12	&gamma;-glutamylanilide dioxygenase complex
PWY-6533	3.5.1.ac	N-glutamylanilide hydrolase
PWY-6533	6.3.1.d	N-glutaminylanilide synthase
PWY-6534	1.2.1.39	NAD-dependent phenylacetaldehyde dehydrogenase
PWY-6534	1.4.99	quinohaemoprotein amine dehydrogenase
PWY-6534	1.4.99 	quinohemoprotein amine dehydrogenase
PWY-6535 	1.2.1.16 	NAD(P)-dependent succinate-semialdehyde dehydrogenase
PWY-6535	1.2.1.24	NAD-dependent succinate semialdehyde dehydrogenase
PWY-6535 	2.6.1.19	4-aminobutyrate transaminase subunit
PWY-6536	1.2.1.16 	NAD(P)-dependent succinate semialdehyde dehydrogenase
PWY-6536	2.6.1.19	4-aminobutyrate aminotransferase
PWY-6537	1.2.1.79	NADP-dependent succinate semialdehyde dehydrogenase
PWY-6538	1.14.13.128	MONOMER-15753
PWY-6538	1.14.13.178	MONOMER-15754
PWY-6538	1.14.13.178	MONOMER-15755
PWY-6538	1.14.13.179	theobromine demethylase
PWY-6538 	1.17.1.4	MONOMER-15740
PWY-6538 	1.17.3.2	MONOMER-15741
PWY-6538	1.5.1.39 	glutathione-S-transferase NdmE 
PWY-6539	1.1.1.M1	sulfenic acid dehydrogenase
PWY-6539	4.4.1.4	cysteine sulfoxide lyase
PWY-6540	1.1.1.314 	germacrene A alcohol dehydrogenase/hydroxylase
PWY-6540	1.14.13.120	costunolide synthase
PWY-6540	1.14.13.123	germacrene A oxidase
PWY-6540	1.14.13.123	MONOMER-13607
PWY-6540	1.14.13	germacrene A acid 8&beta;-hydroxylase
PWY-6540	4.2.3.23	germacrene A synthase
PWY-6544 	1.1.1.357	MONOMER-7322
PWY-6544 	1.1.1 	brassinolide synthase
PWY-6544 	1.1.1 	brassinosteroid-6-oxidase
PWY-6544 	1.14.99	23&alpha; hydroxylase
PWY-6544 	1.1.99.M1	3&beta;-hydroxysteroid dehydrogenase/&Delta;<sup>5</sup>-&Delta;<sup>4</sup>-isomerase
PWY-6544 	1.3.1	3-oxo-5&alpha;-steroid 4-dehydrogenase
PWY-6544 	5.3.3.1	MONOMER-7321
PWY-6545 	1.17.4.1	B12-dependent ribonucleotide reductase
PWY-6545 	1.17.4.1	ribonucleotide reductase
PWY-6545	2.1.1.148	flavin-dependent thymidylate synthase
PWY-6545 	2.7.4.6	nucleoside-diphosphate kinase
PWY-6545	3.6.1.15 	nucleoside-triphosphatase
PWY-6545	3.6.1.23 	dUTPase subunit
PWY-6546 	2.4.1 	UDP glucose:cytokinin glycosyltransferase
PWY-6546	2.8.2	brassinosteroid sulfotransferase
PWY-6549	1.1.1.42	NADP+-dependent isocitrate dehydrogenase
PWY-6549 	4.1.1.31	phosphoenolpyruvate carboxylase
PWY-6550	1.13.11.39	2-aminobiphenyl-2,3-diol dioxygenase complex
PWY-6550	1.14.12.22	carbazole 1,9a-dioxygenase complex
PWY-6550	3.7.1.13	2-hydroxy-6-oxo-6-(2-aminophenyl)-hexa-2,4dienoate hydrolase
PWY-6554 	2.7.1.151 	inositol polyphosphate multiple-kinase 2&alpha;
PWY-6555 	2.7.1.151 	inositol polyphosphate multiple-kinase 2&beta;
PWY-6555 	2.7.1.159	myo-inositol-1,3,4-trisphosphate 5/6-kinase
PWY-6555 	2.7.1 	polyphosphate 2-kinase
PWY-6559 	1.5.1.43	carboxynorspermidine dehydrogenase
PWY-6559	1.5.1.43	MONOMER-17345
PWY-6559 	4.1.1 	carboxynorspermidine decarboxylase
PWY-6559	4.1.1	carboxyspermidine decarboxylase subunit
PWY-6564 	2.4.1.133	xylosylprotein 4-&beta;-galactosyltransferase
PWY-6564 	2.4.1.134	&beta;-1,3-galactosyltransferase 6
PWY-6564 	2.4.1.135	galactosylgalactosylxylosylprotein 3-&beta;-glucuronosyltransferase 1
PWY-6564 	2.4.1.135	galactosylgalactosylxylosylprotein 3-&beta;-glucuronosyltransferase 2
PWY-6564 	2.4.1.135	galactosylgalactosylxylosylprotein 3-&beta;-glucuronosyltransferase 3
PWY-6564 	2.4.1.223	glucuronyl-galactosyl-proteoglycan 4-&alpha;-N-acetylglucosaminyltransferase 1
PWY-6564 	2.4.1.224	N-acetylglucosaminyl-proteoglycan 4-&beta;-glucuronosyltransferase EXTL1
PWY-6564 	2.4.1.225 	exostosin complex
PWY-6564 	2.4.2.26	xylosyltransferase 1
PWY-6564 	2.4.2.26	xylosyltransferase 2
PWY-6564 	2.8.2.29	[heparan sulfate]-glucosamine 3-O-sulfotransferase 2
PWY-6564 	2.8.2.30	[heparan sulfate]-glucosamine 3-O-sulfotransferase 3A1
PWY-6564 	2.8.2.30	[heparan sulfate]-glucosamine 3-O-sulfotransferase 3B1
PWY-6564 	2.8.2	[heparan-sulfate]-6-O-sulfotransferase 1
PWY-6564 	2.8.2	[heparan-sulfate]-6-O-sulfotransferase 2
PWY-6564 	2.8.2	[heparan-sulfate]-6-O-sulfotransferase 3
PWY-6564 	2.8.2 	[heparan sulfate]-glucosamine 3-O-sulfotransferase 1
PWY-6564 	2.8.2	[heparan sulfate]-uronosyl-2-O-sulfotransferase 1
PWY-6564 	3.1.1 	bifunctional heparan sulfate N-deacetylase/N-sulfotransferase 1
PWY-6564 	3.1.1 	bifunctional heparan sulfate N-deacetylase/N-sulfotransferase 2
PWY-6564 	3.1.1 	bifunctional heparan sulfate N-deacetylase/N-sulfotransferase 3
PWY-6564 	3.1.1 	bifunctional heparan sulfate N-deacetylase/N-sulfotransferase 4
PWY-6564 	5.1.3.17	heparosan-N-sulfate-glucuronate 5-epimerase
PWY-6565 	1.2.1.11	aspartate-&beta;-semialdehyde dehydrogenase subunit
PWY-6565 	4.1.1.86 	L-2,4-diaminobutyrate aminotransferase/decarboxylase
PWY-6565 	4.1.1.96	carboxynorspermidine decarboxylase
PWY-6569 	2.8.2.17	chondroitin 6-O-sulfotransferase 1
PWY-6569 	2.8.2.17	chondroitin 6-sulfotransferase 2
PWY-6569 	2.8.2.5	chondroitin 4-sulfotransferase 1
PWY-6569 	2.8.2.5	chondroitin 4-sulfotransferase 2
PWY-6569 	2.8.2.5	chondroitin 4-sulfotransferase 3
PWY-6569 	2.8.2	uronyl 2-sulfotransferase
PWY-6571 	2.4.1.175 	chondroitin sulfate N-acetylgalactosaminyltransferase 1
PWY-6571 	2.4.1.175 	chondroitin sulfate N-acetylgalactosaminyltransferase 2
PWY-6571 	2.4.1.175 	chondroitin sulfate synthase 1
PWY-6571 	2.4.1.175 	chondroitin sulfate synthase 2
PWY-6571 	2.4.1.175 	chondroitin sulfate synthase 3
PWY-6571 	2.4.1.226	chondroitin sulfate glucuronyltransferase
PWY-6571 	2.8.2.33	N-acetylgalactosamine 4-sulfate 6-O-sulfotransferase
PWY-6571 	2.8.2.35	dermatan 4-sulfotransferase 1
PWY-6571 	5.1.3.19	dermatan-sulfate epimerase
PWY-6572	3.1.6.10	chondro-6-sulfatase
PWY-6572	4.2.2.20 	chondroitin AC lyase
PWY-6573	3.1.6.12	arylsulfatase B
PWY-6573	3.1.6.4	N-acetylgalactosamine-6-sulfatase
PWY-6573	3.2.1.35	hyaluronidase-4
PWY-6573	3.2.1.52	&beta;-hexosaminidase A
PWY-6575	1.1.1.216	NADP+-dependent farnesol dehydrogenase 1
PWY-6575	1.14.13.202	farnesyl methyl ester epoxidase
PWY-6575	1.2.1.94	fatty aldehyde dehydrogenase 3
PWY-6575	2.1.1.325	juvenile hormone acid methyltransferase
PWY-6575	3.1.3 	farnesyl phosphatase
PWY-6576	3.1.6.13	iduronate 2-sulfatase
PWY-6576	3.2.1.35	hyaluronidase-1
PWY-6576	3.2.1.35	hyaluronidase PH-20
PWY-6576	3.2.1.76	&alpha;-L-iduronidase
PWY-6577	1.1.1.216	farnesol dehydrogenase
PWY-6577	1.8.3.6	farnesylcysteine lyase
PWY-6577	2.7.1	farnesol kinase
PWY-6577	2.7.1	prenyl alcohol kinase
PWY-6577	2.7.4	farnesyl monophosphate kinase
PWY-6579 	2.4.2.1 	purine nucleoside phosphorylase
PWY-6579 	2.4.2.8	hypoxanthine phosphoribosyltransferase
PWY-6579 	2.7.1.73	inosine-guanosine kinase
PWY-6580	3.1.3.64	1-phosphatidyl-1D-myo-inositol 3-phosphate phosphatase
PWY-6581	1.14.15.9	acyclic carotenoid 2-ketolase
PWY-6581	1.3.99.27	carotenoid 3,4-desaturase
PWY-6585 	3.1.2 	3-oxoacyl-CoA hydrolase
PWY-6585	3.1.2	methylketone synthase 2
PWY-6591	1.11.1.M1 	manganese-oxidizing peroxidase
PWY-6591	1.11.1.M1 	manganese peroxidase 1
PWY-6591 	1.16.3.c	manganese-oxidizing multicopper oxidase
PWY-6591 	1.16.3.c	MnxG manganese-oxidizing multicopper oxidase
PWY-6592 	1.11.1.M1 	manganese-oxidizing peroxidase
PWY-6592	1.16.3.c	CotA manganese-oxidizing multicopper oxidase
PWY-6592	1.16.3.c	MnxG manganese oxidizing multicopper oxidase
PWY-6592	1.16.3.c	MnxG manganese-oxidizing multicopper oxidase
PWY-6592	1.16.3.c	MofA manganese-oxidizing multicopper oxidase
PWY-6593	1.2.1.81	sulfoacetaldehyde dehydrogenase (acylating)
PWY-6593	6.2.1	sulfoacetate-CoA ligase (AMP-forming)
PWY-6595 	3.2.2.1 	ribonucleoside hydrolase
PWY-6596 	3.1.3.99 	5-nucleotidase
PWY-6599	2.4.2.8	hypoxanthine-guanine phosphoribosyltransferase
PWY-6599	3.2.2.1	purine nucleosidase
PWY-6603	1.14.19.38 	acyl-lipid 6-desaturase/acetylenase
PWY-6603	1.14.19.47	acyl-lipid 6-desaturase
PWY-6604 	1.1.1.35	3-hydroxybutyryl-CoA dehydrogenase
PWY-6604 	1.1.1.35	(S)-3-hydroxybutyryl-CoA dehydrogenase
PWY-6604 	1.1.1	NADH-dependent butanol dehydrogenase A
PWY-6604 	1.1.1	NADH-dependent butanol dehydrogenase B
PWY-6604 	1.1.1	NADPH-dependent butanol dehydrogenas
PWY-6604 	1.2.7.1	pyruvate:ferredoxin oxidoreductase
PWY-6604 	1.3.8.1	BCDCLOS-MONOMER
PWY-6604 	2.3.1.19	phosphotransbutyrylase subunit
PWY-6604 	2.3.1.9	acetyl-CoA C-acetyltransferase
PWY-6604 	2.7.2.7	butyrate kinase monomer
PWY-6604 	4.1.1.4	acetoacetate decarboxylase
PWY-6604 	4.2.1.17 	3-hydroxybutyryl-CoA dehydratase
PWY-6605 	2.4.2.7	adenine phosphoribosyltransferase
PWY-6605	2.4.2.7	AT1G27450-MONOMER
PWY-6605 	2.4.2.8 	adenine-guanine phosphoribosyltransferase
PWY-6605	3.2.2.7	adenosine nucleosidase
PWY-6606	3.2.2.1	MONOMER-15864
PWY-6607 	3.2.2.1 	purine nucleosidase
PWY-6608 	1.17.1.4	xanthine dehydrogenase
PWY-6608 	3.1.3.2 	acid phosphatase / phosphotransferase
PWY-6608 	3.1.3.5 	5-nucleotidase
PWY-6608 	3.1.3.5 	broad specificity 5(3)-nucleotidase and polyphosphatase
PWY-6608 	3.1.3.5 	UMP phosphatase
PWY-6608 	3.5.4.3	guanine deaminase
PWY-6609	2.4.2.8	MONOMER-862
PWY-6609	3.5.4.4	MONOMER-841
PWY-6610	2.4.2.7	adenine phosphoribosyltransferase
PWY-6610	2.4.2.8	hypoxanthine guanine phosphoribosyltransferase
PWY-6610	3.5.4.2	adenine deaminase
PWY-66	1.1.1.271	AT1G73250-MONOMER
PWY-6612 	1.5.1.3 	dihydrofolate reductase/thymidylate synthase
PWY-6612 	1.5.1.3 	MONOMER-9361
PWY-6612 	2.6.1.85	AT2G28880-MONOMER
PWY-6612 	2.6.1.85	MONOMER-8741
PWY-6612 	3.5.4.16	GTP cyclohydrolase I
PWY-6612 	3.6.1.67	MONOMER-13396
PWY-6612 	4.1.2.25 	7,8-dihydroneopterin aldolase [multifunctional]
PWY-6612 	4.1.2.25	MONOMER-8764
PWY-6612 	4.1.3.38	ADCL subunits
PWY-6612 	4.1.3.38	AT5G57850-MONOMER
PWY-6612 	6.3.2.12	AT5G41480-MONOMER
PWY-6616	1.1.1.310	S-sulfolactate dehydrogenase subunit
PWY-6616	1.1.1.338	MONOMER-15868
PWY66-161 	1.14.14 	cytochrome P450 2E1
PWY-6616	4.4.1.24	R-sulfolactate sulfo-lyase B subunit 
PWY-6617	3.2.2.4	AMP nucleosidase
PWY-6618	2.7.1.73	MONOMER-11773
PWY-6619	2.7.1.20	adenosine kinase 1
PWY-6619	2.7.1.20	adenosine kinase 2
PWY66-201	1.14.13.148 	dimethylaniline monooxygenase (n-oxide forming) 3
PWY66-201 	1.14.14.1 	cytochrome P450 2A6
PWY66-201	2.1.1.49	indolethylamine N-methyltransferase
PWY66-201 	2.4.1.17	UDP glycosyltransferase 1 family, polypeptide A4
PWY66-21 	6.2.1.1	acetyl-coA synthetase
PWY66-221 	1.2.3.1 	aldehyde oxidase
PWY-6622 	4.1.99.5 	aldehyde decarbonylase
PWY-6622	4.1.99.5	aldehyde decarbonylase
PWY-6622 	6.2.1.l 	long-chain-fatty-acid&mdash;[acyl-carrier-protein] ligase
PWY-6623 	2.4.1 	hydroxybenzoate UDP-glucosyltransferase / quercetin UDP-glucosyltransferase
PWY66-241	1.14.14.1	cytochrome P450 2B6
PWY-6626	1.1.1	1,3-dihydroxyacetone/glyceraldehyde reductase
PWY-6626	2.7.1	MONOMER-15888
PWY-6626	2.7.7	MONOMER-15889
PWY-6627	1.1.1	chloroethylmalonyl-CoA dehydrogenase (decarboxylating)
PWY-6627	2.3.1	salinosporamide A polyketide synthase
PWY-6627	2.4.2.1	MONOMER-15880
PWY-6627	2.5.1.94	adenosyl-chloride synthase
PWY-6627	2.6.1	MONOMER-15891
PWY-6627	3.1.3	MONOMER-15881
PWY-6627	4.1.1	5-chloro-4-hydroxy-2-oxopentanoate decarboxylase (acylating)
PWY-6627	4.2.1	dihydroxy-acid dehydratase
PWY-6627	4.2.1	MONOMER-15884
PWY-6627	5.4.99.5	chorismate mutase
PWY-6628 	5.4.99.5 	chorismate mutase / prephenate dehydrogenase
PWY66-301	1.14.16.2	tyrosine 3-monooxygenase
PWY66-301	1.14.17.1	dopamine &beta;-monooxygenase
PWY-6630 	5.4.99.5 	chorismate mutase / prephenate dehydratase
PWY66-3 	1.14.13.70	cytochrome P450 51A1
PWY-6631	2.1.1	MONOMER-15895
PWY-6632	1.17.3.2	xanthine oxidase
PWY-6633	1.14.13.212	MONOMER-18616
PWY-6633	1.17.5.2	caffeine dehydrogenase
PWY-6634	1.1.1.308	sulfopropanediol 3-dehydrogenase
PWY66-341 	1.14.19.20	lathosterol oxidase
PWY-6634	1.1.1	(R)-sulfopropanediol 2-dehydrogenase
PWY-6634	1.1.1	(S)-sulfopropanediol 2-dehydrogenase
PWY66-341 	1.3.1.21	&Delta; 7-sterol reductase
PWY66-341 	1.3.1.72	24-dehydrocholesterol reductase precursor
PWY66-3 	5.3.3.5	3-&beta;-hydroxysteroid-&Delta;<sup>8</sup>-&Delta;<sup>7</sup>-isomerase
PWY-6636	4.1.1.91	MONOMER-15904
PWY66-366	2.7.1.26	riboflavin kinase
PWY66-366	2.7.7.2	FAD synthetase
PWY66-367 	1.1.1.30	mitochondrial D-&beta;-hydroxybutyrate dehydrogenase
PWY66-367 	4.1.1.4	MONOMER-12921
PWY66-367 	4.1.3.4	mitochondrial hydroxymethylglutaryl-CoA lyase
PWY66-368	1.1.1.30	3-hydroxybutyrate dehydrogenase type 2
PWY66-368	2.8.3.5	succinyl-CoA:3-ketoacid-coenzyme A transferase 1, mitochondrial
PWY66-373	2.7.1.3	ketohexokinase
PWY66-373	3.2.1.26 	sucrase-isomaltase
PWY66-374	1.14.99.1	prostaglandin G/H synthase 1
PWY66-374	5.3.99.2	hematopoietic prostaglandin D synthase
PWY66-374	5.3.99.2	prostaglandin-D2 synthase
PWY66-374	5.3.99.3	prostaglandin E synthase
PWY66-374	5.3.99.3	prostaglandin E synthase 2 truncated form
PWY66-374	5.3.99.3	prostaglandin E synthase 3
PWY66-374	5.3.99.4	prostacyclin synthase
PWY66-374	5.3.99.5	thromboxane-A synthase
PWY66-375	3.3.2 	leukotriene A4 hydrolase
PWY66-375	3.4.13.19 	Dipeptidase 1
PWY66-375	3.4.13.19 	Dipeptidase 2
PWY66-375	4.4.1.20	Leukotriene C4 synthase
PWY66-378 	5.3.3.1 	3 &beta;-hydroxysteroid dehydrogenase/ &Delta; 5-->4-isomerase type 1
PWY66-381 	1.14.14.ag 	steroid 17-alpha-hydroxylase/17,20 lyase
PWY66-382 	1.14.13 	steroid 21-monooxygenase
PWY66-382 	1.14.15.4	steroid 11&beta;-hydroxylase
PWY-6638	4.4.1.25	L-cysteate sulfo-lyase subunit
PWY66-387	1.14.11 	phytanoyl-CoA dioxygenase, peroxisomal
PWY66-387	4.1.2	2-hydroxyacyl-CoA lyase 1
PWY66-387 	6.2.1.3 	long-chain-fatty-acid-CoA ligase 1
PWY66-388	1.14.18.6	sphingoliid fatty acid 2-hydroxylase
PWY66-388 	1.2.1.19 	4-trimethylaminobutyraldehyde dehydrogenase
PWY66-388 	4.1.2	2-hydroxyacyl-CoA lyase 1
PWY66-388 	6.2.1.3 	long-chain acyl-CoA synthetase 3
PWY66-388 	6.2.1.3 	long-chain acyl-CoA synthetase 4
PWY66-388 	6.2.1.3 	long-chain acyl-CoA synthetase 5
PWY66-388 	6.2.1.3 	long-chain acyl-CoA synthetase 6
PWY66-388 	6.2.1.3 	long-chain acyl-CoA synthetase 8 ACSBG2
PWY66-388 	6.2.1.3 	long-chain acyl-CoA synthetase ACSBG1
PWY66-391 	1.3.3.6	peroxisomal acyl-coenzyme A oxidase 1
PWY66-391 	1.3.3.6 	peroxisomal acyl-coenzyme A oxidase 2
PWY66-392	1.13.11 	arachidonate 12-lipoxygenase
PWY66-392 	1.13.11	arachidonate 15-lipoxygenase
PWY66-392 	1.13.11 	arachidonate 5-lipoxygenase
PWY66-394 	1.13.11 	prostaglandin G/H synthase 2
PWY66-398	1.1.1.41	isocitrate dehydrogenase (NAD)
PWY66-398	2.3.3.16 	citrate synthase
PWY66-398	4.2.1.2	fumarate hydratase monomer
PWY66-398	6.2.1.4	(GDP-forming) succinate-CoA ligase
PWY66-398	6.2.1.5	(ADP-forming) succinyl-CoA ligase
PWY66-399	3.1.3.11	fructose-1,6-bisphosphatase 1
PWY66-399	3.1.3.11	fructose-1,6-bisphosphatase 2
PWY66-399	3.1.3.58 	glucose-6-phosphatase
PWY66-399	3.1.3.58 	glucose-6-phosphatase 2
PWY66-399	3.1.3.58 	glucose-6-phosphatase 3
PWY66-399	4.1.1.32	phosphoenolpyruvate carboxykinase, cytosolic [GTP]
PWY66-399	6.4.1.1	pyruvate carboxylase, mitochondrial
PWY-6640 	1.13.11.4	MONOMER-15921
PWY-6640	1.14.13.209	MONOMER-15920
PWY-6640	2.3.1 	MONOMER-15914
PWY-6640	2.3.1	MONOMER-15917
PWY66-409 	2.4.2.8	hypoxanthine-guanine phosphoribosyltransferase
PWY66-409 	2.7.1.20	adenosine kinase
PWY-6641 	1.1.1	MONOMER-15905
PWY-6641 	2.3.1.8	phosphate acetyltransferase
PWY-6641 	2.3.3.15	sulfoacetaldehyde acetyltransferase
PWY-6641 	2.6.1.1	MONOMER-15911
PWY-6641 	4.1.1.79	sulfopyruvate decarboxylase
PWY-6641 	4.4.1.25	L-cysteate sulfo-lyase
PWY-6642	1.1.1.337	MONOMER-15913
PWY-66	4.2.1.47	GDP-mannose 4,6-dehydratase subunit
PWY66-421 	6.3.2.11	carnosine synthase 1
PWY66-422	5.1.3.3	aldose 1-epimerase
PWY-6642	2.6.1.1	MONOMER-15912
PWY66-423	2.7.1.105 	6-phosphofructo-2-kinase/fructose-2,6-bisphosphatase 3
PWY66-423	3.1.3.46 	6-phosphofructo-2-kinase/fructose-2,6-bisphosphatase 1
PWY66-423	3.1.3.46 	6-phosphofructo-2-kinase/fructose-2,6-bisphosphatase 2
PWY66-423	3.1.3.46 	6-phosphofructo-2-kinase/fructose-2,6-bisphosphatase 4
PWY66-423	3.1.3.46	fructose-2,6-bisphosphatase TIGAR
PWY-6642	4.4.1.24	3-sulfolactate sulfo-lyase small subunit 
PWY66-425 	1.5.1 	pyrroline-5-carboxylate reductase 1
PWY66-426 	4.4.1.1	cystathionine gamma-lyase
PWY66-428 	4.3.1.19	L-serine dehydratase/L-threonine deaminase
PWY66-428 	4.3.1.19	serine dehydratase-like
PWY-6643	2.5.1.76	cysteate synthase
PWY-6643	2.6.1.52 	phosphoserine aminotransferase
PWY-6643	4.1.1.79	sulfopyruvate decarboxylase
PWY-6644	1.2.1.69	fluoroacetaldehyde dehydrogenase
PWY-6644	2.2.1.8	4-fluorothreonine transaldolase
PWY-6644	2.4.2.1	5-fluoro-5-deoxy-adenosine phosphorylase
PWY-6644	2.5.1.63	5-fluoro-5-deoxyadenosine synthase
PWY-6644	5.3.1	5-fluoro-5-deoxy-ribose 1-phosphate isomerase
PWY-6645	4.2.1.133	MONOMER-15960
PWY-6646	3.8.1.3	fluoroacetate dehalogenase
PWY-6649	1.1.99.14	glycolate oxidase
PWY-6649	4.1.3.14	L-erythro-3-hydroxyaspartate aldolase
PWY-6649	4.3.1.20	L-erythro-3-hydroxyaspartate ammonia-lyase
PWY-6650	1.2.1.94	farnesal dehydrogenase
PWY-6650	2.1.1.325	juvenile hormone acid methyltransferase
PWY66-5 	1.1.1.34	3-hydroxy-3-methylglutaryl-coenzyme A reductase
PWY66-5 	1.14.14.17	squalene monooxygenase
PWY66-5 	2.3.3.10	hydroxymethylglutaryl-CoA synthase
PWY66-5 	2.3.3.10	mitochondrial hydroxymethylglutaryl-CoA synthase
PWY66-5 	2.5.1.1 	geranylgeranyl diphosphate synthase
PWY66-5 	2.5.1.1 	E-farnesyl diphosphate synthase
PWY66-5 	2.5.1.21 	farnesyl-diphosphate farnesyltransferase
PWY66-5 	2.7.1.36	mevalonate kinase
PWY66-5 	2.7.4.2	phosphomevalonate kinase
PWY-6653	5.5.1.13	MONOMER-15965
PWY66-5 	4.1.1.33	diphosphomevalonate decarboxylase
PWY-6654	2.7.1.169	MONOMER-15970
PWY-6654	6.3.2.36	MONOMER-15971
PWY-6655	2.4.1.251	GlcA-&beta;-(1,2)-D-Man-&alpha;-(1,3)-D-Glc-&beta;-(1,4)-D-Glc-&alpha;-1-diphosphoundecaprenol 4-&beta;-mannosyltransferase
PWY-6655	2.4.1.252	GDP-mannose:cellobiosyl-diphosphopolyprenol &alpha;-mannosyltransferase
PWY-6655	2.4.1.264	MONOMER-15982
PWY-6655	2.4.1	&alpha;-D-glucopyranosyl-diphosphoundecaprenol glucosyltransferase
PWY-6655	2.5.1 	xanthan ketal pyruvate transferase
PWY-6655	2.7.8.31	undecaprenyl-phosphate glucose phosphotransferase
PWY66-5 	5.3.3.2	isopentenyl-diphosphate delta-isomerase 1
PWY66-5 	5.3.3.2	isopentenyl-diphosphate delta-isomerase 2
PWY66-5 	5.4.99.7	lanosterol synthase
PWY-6657	4.2.1.119 	(R)-specific enoyl-CoA hydratase
PWY-6658	2.4.1.252	MONOMER-15989
PWY-6658	2.4.1	MONOMER-15988
PWY-6658	2.7.8.31	MONOMER-15987
PWY-6659	1.1.1	short-chain dehydrogenase/reductase
PWY-6659	1.14.11	fusicocca-2,10(14)-diene-8&beta;,16-diol dioxygenase
PWY-6659	2.5.1.29 	fusicocca-2,10(14)-diene synthase
PWY-6659	2.5.1	diterpene glucoside O-glucose prenyltransferase
PWY-6659 	4.2.3 	phomopsene synthase
PWY-6662 	1.14.13.182	MONOMER-16010
PWY-6662 	2.3.1.230	2-aminobenzoyloctanoate synthase
PWY-6662 	2.3.1.bq	anthraniloyl-CoA malonyltransferase
PWY-6662 	3.1.2.f	2-aminobenzoylacetyl-CoA thioesterase
PWY-6662 	4.1.3.27	anthranilate synthase &beta; subunit 
PWY-6662 	4.1.3.27	PhnB anthranilate synthase subunit 
PWY-6662 	6.2.1.32	MONOMER-16009
PWY-6664	2.7.7.74 	bifubnctional CTP:inositol-1-phosphate cytidylyltransferase/CDP-inositol:inositol-1-phosphate transferase
PWY-6665	2.1.1.240	resveratrol 3,5-O-dimethyltransferase
PWY-6666	1.14.13.218	5-methylphenazine-1-carboxylate 1-monooxygenase
PWY6666-1	3.5.1.99	fatty acid amide hydrolase
PWY-6666	2.1.1.fr	5-methyl-phenazine-1-carboxylate N-methyltransferase
PWY6666-2 	1.4.3.4	amine oxidase A
PWY6666-2 	2.1.1.6	catechol O-methyltransferase, membrane-bound form
PWY-6667	1.13.11	MONOMER-16023
PWY-6668 	1.14.13	AT3G25180-MONOMER
PWY-6668	4.2.3.144	geranyllinalool synthase
PWY-6669	4.2.3.23	MONOMER-16036
PWY-6669	4.2.3.23	MONOMER-16040
PWY-6670	1.3.99	MONOMER-16061
PWY-6672 	4.1.3.26 	hydroxymethylglutaryl-CoA lyase subunit
PWY-6672	4.2.1.57	MONOMER-16076
PWY-6672	6.4.1.5	cis-geranyl-CoA carboxylase &alpha;-subunit 
PWY-6673	2.3.1.98	MONOMER-16057
PWY-6673	2.3.1.98	MONOMER-16156
PWY-6673 	2.3.1 	hydroxycinnamoyl-CoA: shikimate/quinate hydroxycinnamoyltransferase
PWY-6673	2.3.1	MONOMER-16127
PWY-6673	6.2.1.12	MONOMER-16125
PWY-6676 	1.8.5.4	sulfide:quinone reductase
PWY-6676 	1.8.5	sulfite dehydrogenase
PWY-6676 	1.8.99.2	dissimilatory adenylyl-sulfate reductase
PWY-6676 	1.8.99.5	siroheme sulfite reductase, &alpha; subunit 
PWY-6676 	2.7.7.4	dissimilatory sulfate adenylyltransferase
PWY-6676 	2.8.1	DsrC dimer (A. vinosum)
PWY-6676 	2.8.1	DsrEFH hexamer
PWY-6676 	2.8.1	rhodanese 2599
PWY-6676 	2.8.1	TusA sulfur-carrier protein
PWY-6678	1.1.1.347	citrol dehydrogenase
PWY-6678	1.2.1.86	citral dehydrogenase
PWY-6679	2.3.1	2,6-dideoxy-&alpha;-L-ribohexopyranosyl-O-glycosyltransferase
PWY-6679	2.3.1	jadomycin polyketide synthase complex
PWY-6679	2.3.1	UWM6 synthase
PWY-6679	4.2.1	UWM6 dehydratase
PWY-6679	6.4.1.2	acetyl CoA carboxylase
PWY-6681	1.13.11.59	MONOMER-16133
PWY-6681	1.2.1.82 	MONOMER-16130
PWY-6681	2.5.1	MONOMER-16129
PWY-6682	1.1.1	1-hydroxy-2-phosphorylethylphosphonate dehydrogenase
PWY-6682	1.1.1.309	phophonoacetaldehyde reductase
PWY-6682	1.14.11	2-hydroxyethylphosphonate 1-hydroxylase
PWY-6682	2.1.1.M5	phosphonate O-methyltransferase
PWY-6682	2.3.2	[1-(2-amino-4-methylpentanamido)ethenyl]phosphonate glycyltransferase
PWY-6682	2.3.2	1-amino-2-phosphorylethylphosphonate leucinyltransferase
PWY-6682	2.6.1	1-oxo-2-phosphorylethylphosphonate transaminase
PWY-6682	2.7.1	1,2-dihydroxyethylphosphonate kinase
PWY-6682	4.1.1.82	phosphonopyruvate decarboxylase
PWY-6682	5.4.2.9	phosphoenolpyruvate phosphomutase
PWY-6683	1.8.4.10	phosphoadenylyl-sulfate reductase (thioredoxin)
PWY-6683	2.7.7.4	sulfate adenylyltransferase
PWY-6685	2.4.1.268	glucosylglycerate synthase
PWY-6686	2.4.1.266	glucosyl-3-phosphoglycerate synthase
PWY-6686	2.4.1.270	MONOMER-16132
PWY-6686	3.1.3	glucosyl-3-phosphoglycerate phosphatase
PWY-6687	2.4.1.268	glucosylglycerate synthase
PWY-6688 	1.21.99.3	type III iododthyronine deiodinase
PWY-6688 	1.21.99.3	type I iodothyronine deiodinase
PWY-6688 	1.21.99.4	type II iododthyronine deiodinase
PWY-6689	2.7.1.160	tRNA 2-phosphotransferase
PWY-6689	3.1.3.84	ADP-ribose 1-phosphate phosphatase
PWY-6689	3.1.4.37 	cyclic phosphodiesterase
PWY-6689	4.6.1.16	tRNA-splicing endonuclease
PWY-6689	6.5.1.M3 	tRNA ligase
PWY-6690 	1.13.11.16	2,3-dihydroxyphenylpropionate 1,2-dioxygenase
PWY-6690 	1.14.12.19	3-phenylpropionate dioxygenase system
PWY-6690 	1.14.13.127	3-(3-hydroxyphenyl)propionate 2-hydroxylase
PWY-6690 	1.3.1.87	2,3-dihydroxy-2,3-dihydrophenylpropionate dehydrogenase
PWY-6691	1.14.13.110	geranylgeraniol 18-hydroxylase
PWY-6691	2.5.1.29	MONOMER-16157
PWY-6691	3.1.7.5	geranylgeranyl diphosphate diphosphatase I
PWY-6691	3.1.7.5	geranylgeranyl diphosphate diphosphatase II
PWY-6692	1.10.2.2	cytochrome bc1
PWY-6692	1.16.9.1	iron-rusticyanin reductase
PWY-6692	1.9.3.1	aa3-type cytochrome oxidase
PWY-6693 	1.1.1.21	aldose reductase
PWY-6693 	1.1.1 	MONOMER-13196
PWY-6693	1.1.1	MONOMER-17163
PWY-6693	1.1.1	MONOMER-18677
PWY-6693	1.1.1	MONOMER-18678
PWY-6694	1.2.7.10	oxalate oxidoreductase trimer
PWY-6695	2.8.3.16	formyl-coenzyme A transferase
PWY-6695	4.1.1.8	oxalyl-CoA decarboxylase
PWY-6696	1.2.1.17	glyoxylate dehydrogenase (acylating)
PWY-6696	1.2.1.2	NAD-dependent formate dehydrogenase
PWY-6696	2.8.3.2	oxalate CoA-transferase
PWY-6696	4.1.1.8	oxalyl-CoA decarboxylase
PWY-6697	1.2.3.4	oxalate oxidase
PWY-6698	4.1.1.2	oxalate decarboxylase
PWY-6699	3.7.1.1	oxaloacetate acetylhydrolase
PWY-6700	1.17.99.6	epoxyqueuosine reductase
PWY-6700	1.7.1.13	7-cyano-7-deazaguanine reductase
PWY-6700	2.4.2.29	tRNA-guanine transglycosylase
PWY-6700	2.4.99.17	EG10812-MONOMER
PWY-6703	3.5.4.16	GTP cyclohydrolase I
PWY-6703	4.1.2.50	6-carboxy-5,6,7,8-tetrahydropterin synthase
PWY-6703	4.3.99.3	7-carboxy-7-deazaguanine synthase
PWY-6703	6.3.4.20	7-cyano-7-deazaguanine synthase
PWY-6704	1.1.1.366	MONOMER-16192
PWY-6710	1.14.13	CYP77A4
PWY-6711	2.4.2.48	7-cyano-7-deazaguanine tRNA-ribosyltransferase
PWY-6711	2.6.1.97	archaeosine synthase
PWY-6713	1.1.1.173	NAD<sup>+</sup>-dependent L-rhamnose 1-dehydrogenase
PWY-6713	1.1.1.378 	NAD(P)<sup>+</sup>-dependent L-rhamnose 1-dehydrogenase
PWY-6713	1.2.1.22 	NAD(P)<sup>+</sup> L-lactaldehyde dehydrogenase
PWY-6713	1.2.1.22	NAD<sup>+</sup>-dependent L-lactaldehyde dehydrogenase
PWY-6713	3.1.1.65	MONOMER-16238
PWY-6713	4.1.2.53	MONOMER-16242
PWY-6713	4.1.2.53	MONOMER-16243
PWY-6713	4.2.1.90	MONOMER-16239
PWY-6714	1.1.1.378 	NAD(P)<sup>+</sup>-dependent L-rhamnose 1-dehydrogenase
PWY-6714	1.1.1.401	2-dehydro-3-deoxy-L-rhamnonate dehydrogenase
PWY-6714	4.2.1.90	MONOMER-16232
PWY-6717	3.2.1.37	exo-1,4-&beta;-xylosidase
PWY-6717 	3.2.1.8	endo-1,4-&beta;-xylanase XynD
PWY-6718	1.1.1.313	MONOMER-16212
PWY-6718	1.1.1.313	sulfoacetaldehyde reductase
PWY-6720	1.7.1	ToyE
PWY-6720	2.4.2	preQ0 phosphoribosyltransferase
PWY-6720	3.1.3	MONOMER-16220
PWY-6720	4.3.2	MONOMER-16219
PWY-6720	6.3.4	MONOMER-16218
PWY-6721	4.2.1	toyocamycin nitrile hydratase
PWY-6724	2.4.1.25	disproportionating enzyme
PWY-6724	2.7.9.4 	glucan, water dikinase
PWY-6724	2.7.9.5	phosphoglucan, water dikinase
PWY-6724	3.1.3	glucan phosphatase
PWY-6724	3.2.1.2	&beta;-amylase
PWY-6724	3.2.1.68	isoamylase
PWY-6728	1.1.1.42	isocitrate dehydrogenase, NADP
PWY-6728	2.3.3.16 	citrate synthase
PWY-6728	4.1.3.24	&beta;-methylmalyl-CoA lyase
PWY-6728	4.2.1.148	mesaconyl-CoA hydratase
PWY-6728	4.3.1.2	methylaspartate ammonia-lyase
PWY-6728	5.4.99.1	glutamate mutase
PWY-6728	5.4.99.2	methylmalonyl-CoA mutase
PWY-6728	6.4.1.3	propionyl-CoA carboxylase
PWY-6730	2.1.1.165	methyl halide transferase
PWY-6730	2.1.1.165	MONOMER-16271
PWY-6730	2.1.1.165	MONOMER-16279
PWY-6730	2.1.1.165	MONOMER-16280
PWY-6731	2.4.1.19	cyclodextrin glucanotransferase subunit
PWY-6731	2.4.1.1	maltodextrin phosphorylase subunit
PWY-6731	3.2.1.54	MONOMER-16275
PWY-6731 	5.4.2.2	phosphoglucomutase subunit
PWY-6733	1.1.1	AT1G68540-MONOMER
PWY-6733	1.1.1	AT4G35420-MONOMER
PWY-6733 	1.14.13.205	long-chain fatty acid &omega;-hydroxylase
PWY-6733	1.14.13.205 	long-chain fatty acid &omega;-hydroxylase
PWY-6733	1.14.13.206	laurate 7-monooxygenase
PWY-6733 	1.14.13 	fatty acid &omega;-hydroxylase
PWY-6733	1.2.1.84	alcohol-forming fatty acyl-CoA reductase
PWY-6733 	3.1.2.2 	acyl-CoA thioesterase
PWY-6733 	4.2.99.20 	PHYLLO [multifunctional]
PWY-6733	6.2.1.3 	fatty acyl-CoA synthetase
PWY-6735	2.4.1.19	MONOMER-16282
PWY-6735	3.2.1 	cyclodextrinase
PWY-6736	2.1.1.9	MONOMER-16293
PWY-6736	2.1.1.9	MONOMER-16294
PWY-6736	2.1.1.9	MONOMER-16295
PWY-6736 	2.1.1.9	MONOMER-16766
PWY-6736	2.1.1	MONOMER-16290
PWY-6736	2.1.1	MONOMER-16292
PWY-6737	2.4.1.1	&alpha;-glucan/maltodextrin phosphorylase
PWY-6737	2.4.1.25	MONOMER-16304
PWY-6737	3.2.1.41 	amylopullulanase
PWY-6738	2.1.1.129	MONOMER-16307
PWY-6739	1.1.1.142	MONOMER-16309
PWY-6739	1.1.1.143	MONOMER-16308
PWY-6744	1.12.1.4	bifurcating [FeFe]-hydrogenase (soluble)
PWY-6745	2.3.2.15	MONOMER-16319
PWY-6745	2.3.2.15	MONOMER-16321
PWY-6748	1.7.5.2	MONOMER-13457
PWY-6748	1.7.5.2	nitric oxide reductase (menaquinol)
PWY-6749	2.3.1	MONOMER-16346
PWY-6749	2.3.1	MONOMER-16351
PWY-6749	2.5.1.101	N,N-diacetyllegionaminate synthase
PWY-6749	2.6.1.16	L-glutamine-D-fructose-6-phosphate isomerase subunit 
PWY-6749	2.6.1	MONOMER-16350
PWY-6749	2.7.7.82	CMP-legionaminic acid synthase
PWY-6749	2.7.7	glucosamine-1-phosphate guanylyltransferase
PWY-6749	3.2.1	MONOMER-16352
PWY-6749	4.2.1	MONOMER-16348
PWY-6749	5.4.2.10	MONOMER-16332
PWY-6751 	1.12.1.2	NADH-dependent hydrogenase
PWY-6751 	1.12.7.2	ech hydrogenase
PWY-6752	1.10.3.1	MONOMER-16342
PWY-6752	1.10.3.1	MONOMER-16349
PWY-6752	1.10.3.1	polyphenol oxidase
PWY-6753	2.4.2.44	S-methyl-5-thioinosine phosphorylase
PWY-6753	2.4.2.44 	purine nucleoside phosphorylase
PWY-6753	3.5.4.2 	adenosine deaminase
PWY-6759	1.12.7.2	ech hydrogenase
PWY-6759	1.12.7.2	Mbh hydrogenase
PWY-6765	1.12.1.3	hydrogenase delta subunit 
PWY-6765	1.12.1.5 	sulfhydrogenase II &delta; subunit 
PWY-6767	1.14.99.44	diapolycopene oxygenase
PWY-6767	1.2.99	4,4-diapolycopenedial dehydrogenase
PWY-6767	1.3.8.2 	diapophytoene desaturase
PWY-6769	3.1.1.86	rhamnogalacturonan acetylesterase
PWY-6769	3.2.1.171	rhamnogalacturonan hydrolase
PWY-6769	3.2.1.173	rhamnogalacturonan galacturonohydrolase
PWY-6769	3.2.1.174	rhamnogalacturonan rhamnohydrolase
PWY-6769	4.2.2.23	rhamnogalacturonan endolyase
PWY-6771	3.2.1.172	unsaturated &alpha;-galacturonyl hydrolase
PWY-6771	4.2.2.23	rhamnogalacturonan endolyase
PWY-6771	4.2.2.24	rhamnogalacturonan exolyase
PWY-6772	1.1.7 	formate dehydrogenase H
PWY-6773	2.4.1.34	callose synthase
PWY-6773	2.4.1.34	MONOMER-16406
PWY-6778	2.4.1.31	MONOMER-16426
PWY-6780	1.12.7.2	carbon monoxide induced hydrogenase
PWY-6780	1.2.7.4	carbon monoxide dehydrogenase
PWY-6781	3.1.1.42	MONOMER-16433
PWY-6784	3.1.1.73 	endo-1,4-&beta;-xylanase XynY
PWY-6784	3.1.1.73 	endo-1,4-&beta;-xylanase XynZ
PWY-6784	3.2.1.4 	endo-&beta;-1,4-glucanase XghA
PWY-6784	3.2.1.8 	cellulase J
PWY-6784	3.2.1.8	endo-1,4-&beta;-xylanase XynA
PWY-6784	3.2.1.8	endo-1,4-&beta;-xylanase XynC
PWY-6786	1.3.1	alkenal/one oxidoreductase
PWY-6787	1.14.11.22	MONOMER-16536
PWY-6787	1.14.11.23	MONOMER-16535
PWY-6787	1.14.11.9	MONOMER-16523
PWY-6787	2.3.1.74	MONOMER-16518
PWY-6788	3.2.1.21	&beta;-glucosidase Bgl1
PWY-6788	3.2.1.21	&beta;-glucosidase Bgl2
PWY-6788	3.2.1.4	cellulase EG I
PWY-6788	3.2.1.4	cellulase EG II
PWY-6788	3.2.1.4	cellulase EG III
PWY-6788	3.2.1.4	cellulase EG IV
PWY-6788	3.2.1.4	cellulase EG V
PWY-6788	3.2.1.4	cellulase VII
PWY-6788	3.2.1.91	cellulose 1,4-&beta;-cellobiosidase Cbh2
PWY-6788	3.2.1.91	cellulose 1,4-&beta;-cellobiosidase CbhI
PWY-6789	3.2.1.32	endo-1,3-&beta;-xylanase
PWY-6789	3.2.1.72	xylan 1,3-&beta;-xylosidase
PWY-6790	3.2.1.55	&alpha;-N-arabinofuranosidase A
PWY-6791 	3.2.1.120	oligoxyloglucan &beta;-glycosidase
PWY-6791	3.2.1.151	xyloglucan-specific endo-beta-1,4-glucanase A
PWY-6791	3.2.1.23 	&beta;-galactosidase
PWY-6791	3.2.1.51	&alpha;-fucosidase A
PWY-6792 	1.14.11	4-coumaroyl 2-hydroxylase
PWY-6792 	1.14.11	AT3G13610-MONOMER
PWY-6792 	2.1.1.104 	AT1G67990-MONOMER
PWY-6792 	2.1.1.104	caffeoyl-CoA O-methyltransferase
PWY-6794	2.7.7.51	MONOMER-16551
PWY-6797	3.1.4.56	7,8-dihydro-D-neopterin 2,3-cyclic phosphate phosphodiesterase
PWY-6797	3.5.4.39	GTP cyclohydrolase IV
PWY-6797 	4.1.2.25	dihydroneopterin aldolase
PWY-6799	2.3.1	AT2G04540-MONOMER
PWY-6799	6.2.1	MONOMER-16560
PWY-6799	6.2.1	MONOMER-16562
PWY-6801	2.3.1.196	MONOMER-16566
PWY-6801	2.3.1.84	acetyl-CoA:alcohol O-acyltransferase
PWY-6801	2.3.1.84	MONOMER-16565
PWY-6802	4.1.1.25	MONOMER-16586
PWY-6802	4.1.1.25	MONOMER-16587
PWY-6803	2.3.1.23	MONOMER-16595
PWY-6804	2.7.8.2	MONOMER-16594
PWY-6805	2.4.1.20	cellobiose phosphorylase
PWY-6805	3.2.1.91	cellulose 1,4-&beta;-cellobiosidase CbhA
PWY-6805	3.2.1 	cellulase R
PWY-6806	1.14.99 	carotenoid cleavage dioxygenase
PWY-6806	1.14.99	MONOMER-16606
PWY-6807	3.2.1.155	xyloglucan-specific exo-&beta;-1,4-glucanase
PWY-6807	3.2.1.23 	&beta;-galactosidase
PWY-6809 	5.3.99.9	MONOMER-12165
PWY-6809 	5.3.99.9	MONOMER-16625
PWY-681	1.14.14.22	dibenzothiophene-5,5-dioxide monooxygenase
PWY-6812	3.2.1.150	oligoxyloglucan reducing-end-specific cellobiohydrolase
PWY-681	3.13.1.3	2&prime;-hydroxybiphenyl-2-sulfinate desulfinase
PWY-6813	3.2.1.136	glucuronoarabinoxylan endo-1,4-&beta;-xylanase
PWY-6815	3.2.1.178	&beta;-porphyranase A
PWY-6815	3.2.1.178	&beta;-porphyranase B
PWY-6815 	3.2.1.81	&beta;-agarase A
PWY-6815 	3.2.1.81	&beta;-agarase B
PWY-6816	3.2.1.159	&alpha;-neoagaro-oligosaccharide hydrolase
PWY-6816 	3.2.1.159	neoagaro-oligosaccharide 1,3-&alpha;-3,6-anhydro-L-galactosidase
PWY-6816	3.2.1.81	&beta;-agarase A
PWY-6817	3.2.1.162	&lambda;-carrageenase
PWY-6818	2.3.1	L-ornithine N 3-hydroxy acyltransferase
PWY-6818	2.3.1	MONOMER-16664
PWY-6818	2.3.1	ornithine lipid synthase
PWY-6820	1.14.13	MONOMER-16661
PWY-6821	3.1.6	neocarrabiose sulfate/neotetraose disulfate 4-sulfatase
PWY-6821	3.2.1.83	&kappa;-carrageenase
PWY-6821	3.2.1	neocarratetraose 4-O-monosulfate &beta;-hydrolase
PWY-6822	3.2.1.157	&iota;-carrageenase
PWY-6823	2.10.1.1 	AT5G20990-MONOMER
PWY-6823	2.10.1.1 	gephyrin
PWY-6823	2.10.1.1	molybdopterin molybdenumtransferase
PWY-6823	2.7.7.75	molybdopterin adenylyltransferase
PWY-6823	2.7.7.80	molybdopterin-synthase adenylyltransferase
PWY-6823	2.8.1.11 	MOCS3 adenylyltransferase/sulfurtransferase
PWY-6823	2.8.1.11	molybdopterin synthase sulfurtransferase
PWY-6823	2.8.1.12	molybdopterin synthase
PWY-6823	2.8.1.12	MPT synthase large subunit 
PWY-6823	4.6.1.b	cyclic pyranopterin monophosphate synthase
PWY-6823	5.5.1.n	GTP 3,8-cyclase
PWY-6824	1.23.1.4 	MONOMER-16668
PWY-6825	2.1.1.17	MONOMER-16693
PWY-6825	2.1.1.17	phosphatidylethanolamine methyltransferase
PWY-6825	2.1.1.71	MONOMER-16694
PWY-6825	2.1.1.71 	phosphatidylethanolamine N-methyltransferase
PWY-6825	2.1.1.71	phospholipid methyltransferase
PWY-6826	2.7.8.24	MONOMER-16699
PWY-6827	3.2.1.180 	unsaturated &beta;-glucuronyl hydrolase
PWY-6827	4.2.2.25	gellan lyase
PWY-6828	2.1.1.224	23S rRNA (adenine<sup>2503</sup>-C<sup>8</sup>)-methyltransferase
PWY-6829	2.1.1.203 	multisite-specific tRNA:(cytosine-C<sup>5</sup>)-methyltransferase
PWY-6829	2.1.1.205	tRNA (cytidine<sup>32</sup>/guanosine<sup>34</sup>-2-O)-methyltransferase
PWY-6829	2.1.1.211	tRNA<small><sup>Ser</sup></small> (uridine<small><sup>44</sup></small>-2&prime;-O)-methyltransferase
PWY-6829	2.1.1.214	tRNA (guanine<sup>10</sup>-N<sup>2</sup>)-methyltransferase
PWY-6829	2.1.1.216	tRNA (guanine<small><sup>26</sup></small>-N<small><sup>2</sup></small>)-dimethyltransferase
PWY-6829	2.1.1.220	tRNA (adenine<sup>58</sup>-N<sup>1</sup>)-methyltransferase
PWY-6829	2.1.1.225	tRNA:m<small><sup>4</sup></small>X modification enzyme
PWY-6829	2.1.1.229	tRNA (carboxymethyluridine<small><sup>34</sup></small>-5-O)-methyltransferase
PWY-6830 	1.12.98.2	Hmd
PWY-6830 	1.2.7.4	acetyl-CoA decarbonylase/synthase complex &alpha;2&epsilon;2 component
PWY-6830 	1.2.99.5	FmdA 
PWY-6830 	1.2.99.5	formylmethanofuran dehydrogenase, molybdenum enzyme
PWY-6830 	1.2.99.5	formylmethanofuran dehydrogenase, tungsten enzyme
PWY-6830 	1.5.98.2	Mer
PWY-6830 	2.1.1.245	acetyl-CoA decarbonylase/synthase complex  &gamma;&delta; component
PWY-6830 	2.1.1.248	methylamine--corrinoid protein Co-methyltransferase
PWY-6830 	2.1.1.249	dimethylamine--corrinoid protein Co-methyltransferase
PWY-6830 	2.1.1.250	trimethylamine--corrinoid protein Co-methyltransferase
PWY-6830 	2.1.1.86	methyl-H<SUB>4</SUB>MPT:coenzyme M methyltransferase complex
PWY-6830 	2.1.1.90 	methanol-5-hydroxybenzimidazolylcobamide Co-methyltransferase
PWY-6830 	2.1.1.90 	methanol-CoM methyltransferase complex
PWY-6830 	2.1.1.90 	MTAAMBARK-MONOMER
PWY-6830 	2.1.1	CPLX-322
PWY-6830 	2.1.1	CPLX-481
PWY-6830 	2.3.1.101	Ftr
PWY-6830 	2.3.1.169	acetyl-CoA decarbonylase/synthase complex &beta; subunit
PWY-6830 	2.3.1.8	PTAMSARC-MONOMER
PWY-6830 	2.3.1 	FTRMBARK-MONOMER
PWY-6830 	2.7.2.1	Ack
PWY-6830 	2.8.4.1	&alpha; subunit 
PWY-6830 	2.8.4.1	methyl-coenzyme M reductase I
PWY-6830 	2.8.4.1	methyl-coenzyme M reductase II
PWY-6830 	3.5.4.27	Mch
PWY-6830 	3.5.4 	Mch
PWY-6831	1.14.13	aminopyrrolnitrin oxygenase subunit
PWY-6831	1.14.19.9	tryptophan 7-halogenase subunit
PWY-6832	1.2.1	MONOMER-16711
PWY-6834	2.5.1.104	agmatine aminopropyltransferase subunit
PWY-6834	2.5.1.104	polyamine aminopropyltransferase
PWY-6834	3.5.3.24	aminopropylagmatine ureohydrolase
PWY-6834	4.1.1.19	arginine decarboxylase &alpha;&beta; subunit
PWY-6835	2.1.1	MONOMER-16745
PWY-6836	4.2.3.81 	MONOMER-16738
PWY-6836	4.2.3.81 	MONOMER-16739
PWY-6836	4.2.3.82	MONOMER-16740
PWY-6837	1.3.3.6	acyl-CoA oxidase
PWY-6837	1.3.3.6	AT1G06290-MONOMER
PWY-6837	1.3.3.6	AT3G51840-MONOMER
PWY-6837	1.3.3.6	AT4G16760-MONOMER
PWY-6837 	5.3.3.8	delta3, delta2-enoyl-CoA isomerase
PWY-6837	5.3.3	AT5G43280-MONOMER
PWY-6839	4.1.1.82	MONOMER-16743
PWY-6839	5.4.2.9	phosphoenolpyruvate mutase subunit
PWY-6840	6.3.2.23	MONOMER-16747
PWY-6840	6.3.2.2	MONOMER-16746
PWY-6841 	2.3.2 	MONOMER-16320
PWY-6841	2.3.2	MONOMER-16749
PWY-6842	2.3.1	MONOMER-16763
PWY-6842	2.5.1.18	AT1G17170-MONOMER
PWY-6842	2.5.1.18	AT1G78380-MONOMER
PWY-6842	2.5.1.18	AT2G47730-MONOMER
PWY-6842	3.4.17	AT5G44070-MONOMER
PWY-6842	3.4.19.9	AT4G29210-MONOMER
PWY-6842	4.4.1.13 	MONOMER-16765
PWY-6845	1.14.13.39	MONOMER-16748
PWY-6848	1.13.11.24	MONOMER-16757
PWY-6848	1.13.11.24	quercetin 2,3-dioxygenase
PWY-6848	3.1.1	phenol carboxylic acid acyl esterase
PWY-6848	3.2.1	MONOMER-16764
PWY-6852	1.14.13.101	MONOMER-16767
PWY-6852	1.14.13.101	MONOMER-16771
PWY-6852	1.14.13.101	MONOMER-16772
PWY-6854	1.16.1.7	NADH:Fe(III)EDTA oxidoreductase
PWY-6855	2.7.1.147	ADP-dependent glucokinase
PWY-6855	3.2.1.ak 	chitinase, containing dual catalytic domains
PWY-6855	3.2.1	exo-&beta;-D-glucosaminidase
PWY-6855	3.5.1.33 	diacetylchitobiose deacetylase
PWY-6855	3.5.99.6	archaeal glucosamine-6-phosphate deaminase
PWY-6857	1.13.11.63	&beta;,&beta;-carotene 15,15-dioxygenase
PWY-6857	3.1.1.3	pancreatic triacylglycerol lipase
PWY-6857	3.1.1 	liver carboxylesterase 1
PWY-6859	3.1.3.1 	farnesyl diphosphatase
PWY-6861	1.1.1 	11-cis retinol dehydrogenase
PWY-6861 	1.1.1.300	retinol dehydrogenase 11
PWY-6861 	1.1.1.300	retinol dehydrogenase 3
PWY-6861 	1.1.1.300	retinol dehydrogenase 8
PWY-6861 	1.1.1.315 	retinol dehydrogenase 10
PWY-6861 	1.1.1 	retinol dehydrogenase 12
PWY-6861 	1.1 	retinol dehydrogenase 9
PWY-6861 	2.3.1.135	lecithin retinol acyltransferase
PWY-6861	3.1.1.64	retinoid isomerohydrolase
PWY-6863 	1.1.1.1 	aldehyde-alcohol dehydrogenase
PWY-6863 	1.1.1.1 	ethanol dehydrogenase / alcohol dehydrogenase
PWY-6863 	1.1.1 	aldehyde/alcohol dehydrogenase
PWY-6863	1.2.7.1	PYRUFLAVREDUCT-MONOMER
PWY-6863 	2.3.1.16 	FadI component of anaerobic fatty acid oxidation complex
PWY-6863 	2.3.1.16 	fatty acid oxidation complex, &beta; component
PWY-6871 	1.1.1.85	3-isopropylmalate dehydrogenase
PWY-6871 	4.1.1.1 	&alpha;-ketoisovalerate decarboxylase
PWY-6872	1.1.1.105	epidermal retinol dehydrogenase 2
PWY-6872	1.1.1.105 	retinol dehydrogenase 16
PWY-6872	1.2.1.36	retinal dehydrogenase 1
PWY-6872	1.2.1.36	retinal dehydrogenase 2
PWY-6872	1.2.1.36	retinal dehydrogenase 3
PWY-6873	1.1.5.5	alcohol dehydrogenase II
PWY-6873 	2.3.1.75 	wax ester synthase/acyl-CoA diacylglycerol acyltransferase bifunctional enzyme
PWY-6875 	1.17.3.2 	xanthine dehydrogenase
PWY-6876	1.1.1.80	secondary alcohol dehydrogenase
PWY-6876 	3.1.2.11 	butyrate--acetoacetate CoA-transferase
PWY-6886 	1.18.1.2 	ferredoxin-NADP+ oxidoreductase
PWY-6886 	1.3.1.44	mitochondrial trans-2-enoyl-CoA reductase
PWY-6886	1.8.1.4 	pyruvate dehydrogenase complex
PWY-6886	2.7.1.40	pyruvate kinase
PWY-6886 	4.2.1.11	enolase
PWY-6886 	5.4.2.12	phosphoglycerate mutase
PWY-6887	5.5.1.13	ent-copalyl diphosphate synthase
PWY-6888	5.5.1.17 	MONOMER-15969
PWY-6895 	2.2.1.7	1-deoxyxylulose-5-phosphate synthase
PWY-6895 	2.7.1.49 	4-amino-2-methyl-5-phosphomethylpyrimidine kinase
PWY-6895 	2.7.4.16	thiamine monophosphate kinase
PWY-6895 	2.7.7.73	[sulfur-carrier protein ThiS] adenylyltransferase
PWY-6895 	2.8.1.10	thiazole synthase
PWY-6895 	4.1.99.17	phosphomethylpyrimidine synthase
PWY-6895 	5.3.99.10	thiazole tautomerase
PWY-6896	2.7.1.89	thiamin kinase
PWY-6897 	2.5.1.3	thiamine phosphate synthase
PWY-6897	2.7.1.50	hydroxyethylthiazole kinase
PWY-6898 	2.7.6.2	thiamine pyrophosphokinase
PWY-6898	2.7.6.2	thiamine pyrophosphokinase
PWY-6899	3.5.1	N-formyl-4-amino-5-aminomethyl-2-methylpyrimidine deformylase
PWY-6899 	3.5.99.2 	HMP/HMP-P kinase and thiaminase II
PWY-6899	3.5.99.2	thiaminase II
PWY-6900	4.4.1.4	MONOMER-16811
PWY-6900	5.3.99.M1	MONOMER-16812
PWY-6901	1.2.1.12	glyceraldehyde-3-phosphate dehydrogenase
PWY-6901	5.3.1.5	xylose isomerase
PWY-6902 	3.2.1.14	chitinase C1
PWY-6902 	3.2.1.14	chitinase C2
PWY-6902	3.2.1.14	MONOMER-16643
PWY-6902	3.2.1.14	MONOMER-16645
PWY-6902	3.2.1.14	MONOMER-16820
PWY-6902 	3.2.1.52	chitobiase
PWY-6902	3.2.1.52	periplasmic &beta;-N-acetylglucosaminidase
PWY-6902	3.2.1.am	endo-chitodextrinase (N,N-diacetylchitobiose-forming)
PWY-6906	2.4.1.280	MONOMER-16878
PWY-6906	2.7.1.59	MONOMER-16883
PWY-6906	2.7.1.8	MONOMER-16888
PWY-6906	3.5.1.25	N-acetylglucosamine-6-phosphate deacetylase subunit
PWY-6906	3.5.99.6	MONOMER-16884
PWY-6907	2.5.1.3	thiamine phosphate synthase
PWY-6907	3.1.3.100	thiamine monophosphate phosphatase
PWY-6908 	2.7.1.49 	4-amino-2-methyl-5-phosphomethylpyrimidine kinase
PWY-6914	1.14.13.103	MONOMER-16839
PWY-6914	2.5.1.70	MONOMER-16837
PWY-6914	2.5.1.70	MONOMER-16838
PWY-6914	2.5.1.71	MONOMER-16840
PWY-6915	1.1.1.340	11&beta;-hydroxy-1-deoxypentalante dehydrogenase
PWY-6915	1.14.11.35	1-deoxypentalenate dioxygenase
PWY-6915	1.14.11	pentalenolactone D dioxygenase
PWY-6915	1.14.13.170	1-deoxy-11-oxopentalenate oxygenase
PWY-6915	1.14.13	pentalenene oxygenase
PWY-6915	1.14.19.8	MONOMER-16836
PWY-6915	4.2.3.7	pentalenene synthase
PWY-6917	1.11.2.3	oat seed peroxygenase
PWY-6917	1.13.11 	9S-lipoxygenase
PWY-6919	1.1.1.340	11&beta;-hydroxy-1-deoxypentalante dehydrogenase
PWY-6919	1.14.11.35	1-deoxypentalenate 11-hydroxylase
PWY-6919	1.14.11	neopentalenolactone D dioxygenase
PWY-6919	1.14.13.171	1-deoxy-11-oxopentalenate monooxygenase (neopentalenolactone D producing)
PWY-6919	1.14.13	pentalenene oxygenase
PWY-6919	1.14.15.11	1-deoxypentalenate monooxygenase
PWY-6919	4.2.3.7	pentalenene synthase
PWY-6920	2.3.1 	curcuminoid synthase
PWY-6920 	6.2.1.3 	medium-chain acyl-CoA synthetase
PWY-6920	6.2.1.3	MONOMER-16866
PWY-6922	2.3.1	ornithine N<sup>&delta;</sup>-acetyltransferase
PWY-6923	3.5.5.2	MONOMER-16942
PWY-6926	2.5.1.67	chrysanthemyl diphosphate synthase
PWY-6927	3.1.1.14	pheophytinase
PWY-6928 	1.13.11.25	3,4-dihydroxy-9,10-secoandrosta-1,3,5(10)-triene-9,17-dione 4,5-dioxygenase
PWY-6928 	1.1.3.6 	cholesterol oxidase
PWY-6928 	1.14.13.141	cholest-4-en-3-one 26-monooxygenase [(25S)-3-oxocholest-4-en-26-oate forming]
PWY-6928 	1.14.13 	3-ketosteroid 9&alpha;-hydroxylase complex
PWY-6928 	1.14.14.12	3-hydroxy-9,10-seconandrost-1,3,5(10)-triene-9,17-dione 4-hydroxylase monomer
PWY-6928 	1.2.1.87	propanal dehydrogenase
PWY-6928 	1.3.99	3-oxo-23,24-bisnorchol-4-en-22-oyl-CoA dehydrogenase
PWY-6928 	1.3.99.4	3-ketosteroid-&Delta;<sup>1</sup>-dehydrogenase
PWY-6928 	2.3.1	3-oxo-23,24-bisnorchol-4-en-17-ol-22-oyl-CoA lyase
PWY-6928 	3.7.1.17	4,5-9,10-diseco-3-hydroxy-5,9,17-trioxoandrosta-1(10),2-diene-4-oate hydrolase
PWY-6928 	4.1.3.43	4-hydroxy-2-oxohexanoate aldolase
PWY-6928 	4.2.1.132	2-hydroxyhexa-2,4-dienoate hydratase
PWY-6928 	4.2.1	3-oxo-23,24-bisnorchol-4,17(20)-dien-22-oyl-CoA hydratase
PWY-6928 	6.2.1.41	3-[(3aS,4S,7aS)-7a-methyl-1,5-dioxo-octahydro-1H-inden-4-yl]propanoyl:CoA ligase
PWY-6928 	6.2.1.42	3-oxocholest-4-en-26-oate&ndash;CoA ligase
PWY-6930	2.3.1	phenolic glucoside malonyltransferase
PWY-6937 	1.1.1.50 	3&alpha;-hydroxysteroid dehydrogenase/carbonyl reductase
PWY-6937 	1.1.1.51 	3&beta;-hydroxysteroid dehydrogenase
PWY-6937 	1.13.11.25	3,4-dihydroxy-9,10-secoandrosta-1,3,5(10)-triene-9,17-dione 4,5-dioxygenase
PWY-6937 	1.14.13.142	3-ketosteroid 9&alpha;-hydroxylase reductase component
PWY-6937 	1.2.1.87	propanal dehydrogenase
PWY-6937 	1.3.99.4	3-oxosteroid 1-dehydrogenase
PWY-6937 	1.3.99.5	3-keto-5&alpha;-steroid &Delta;4-dehydrogenase
PWY-6937 	1.3.99	5-hydroxy-3-[(3aS,4S,5R,7aS)-7a-methyl-1,5-dioxo-octahydro-1H-inden-4-yl]propanoyl-CoA dehydrogenase
PWY-6937 	3.7.1.17	4,5-9,10-diseco-3-hydroxy-5,9,17-trioxoandrosta-1(10),2-diene-4-oate hydrolase
PWY-6937 	4.1.3.43	4-hydroxy-2-oxohexanoate aldolase
PWY-6937 	4.2.1.132	2-hydroxyhexa-2,4-dienoate hydratase
PWY-6938	4.2.1.93	ATP-dependent NAD(P)H-hydrate dehydratase
PWY-6938	5.1.99.6	NADHX epimerase
PWY-6940	1.14.19.25	acyl-lipid &omega;-3 desaturase
PWY-6941	1.14.14.11	styrene monooxygenase monooxygenase component
PWY-6941	1.2.1.39	phenylacetaldehyde dehydrogenase
PWY-6941	5.3.99	epoxystyrene isomerase
PWY-6942	2.1.1.234	dTDP-3-amino-3,4,6-trideoxy-&alpha;-D-glucopyranose N,N-dimethyltransferase
PWY-6942 	2.6.1.106 	dTDP-2,6-dideoxy-D-glycero-hex-2-enos-4-ulose transaminase
PWY-6942	2.6.1.106	dTDP-3-oxo-3,4,6-trideoxy-&alpha;-D-glucopyranose transaminase
PWY-6942	2.6.1.33	dTDP-4-amino-4,6-dideoxy-<small>D</small>-glucose transaminase
PWY-6942	2.7.7.24	glucose-1-phosphate thymidylyltransferase
PWY-6942	4.2.1.46	dTDP-glucose 4,6-dehydratase
PWY-6942	4.3.1.30	dTDP-4-amino-4,6-dideoxy-D-glucose deaminase
PWY-6945 	5.3.3.1 	3 &beta;-hydroxysteroid dehydrogenase
PWY-6947 	1.13.11.25	3,4-dihydroxy-9,10-secoandrosta-1,3,5(10)-triene-9,17-dione 4,5-dioxygenase
PWY-6947 	1.1.3.6 	cholesterol oxidase
PWY-6947 	1.14.13.141	cholest-4-en-3-one 26-monooxygenase [(25S)-3-oxocholest-4-en-26-oate forming]
PWY-6947 	1.14.13.142	3-ketosteroid 9&alpha;-hydroxylase
PWY-6947 	1.14.13.221	cholest-4-en-3-one C26 monooxygenase [(25R)-3-oxocholest-4-en-26-oate forming]
PWY-6947 	1.14.14.12	3-hydroxy-9,10-seconandrost-1,3,5(10)-triene-9,17-dione 4-monooxygenase
PWY-6947 	1.2.1.87	propanal dehydrogenase
PWY-6947 	1.3.99	3-oxo-23,24-bisnorchol-4-en-22-oyl-CoA dehydrogenase
PWY-6947 	1.3.99	3-oxochol-4-en-24-oyl-CoA dehydrogenase
PWY-6947 	1.3.99	3-oxocholest-4-en-26-oyl-CoA dehydrogenase
PWY-6947 	1.3.99.4	3-oxosteroid 1-dehydrogenase
PWY-6947 	1.3.99	5-hydroxy-3-[(3aS,4S,5R,7aS)-7a-methyl-1,5-dioxo-octahydro-1H-inden-4-yl]propanoyl-CoA dehydrogenase
PWY-6947 	2.3.1.16 	steroid 3-ketoacyl-CoA thiolase
PWY-6947 	2.3.1	3-oxo-23,24-bisnorchol-4-en-17-ol-22-oyl-CoA lyase
PWY-6947 	3.7.1.17	4,5-9,10-diseco-3-hydroxy-5,9,17-trioxoandrosta-1(10),2-diene-4-oate hydrolase
PWY-6947 	4.1.3.43	4-hydroxy-2-oxohexanoate aldolase
PWY-6947 	4.2.1.132	G185E-7813-MONOMER
PWY-6947 	4.2.1	3-oxo-23,24-bisnorchol-4,17(20)-dien-22-oyl-CoA hydratase
PWY-6947 	4.2.1	3-oxocholest-4,24-dien-26-oyl-CoA hydratase
PWY-6947 	5.1.99.4	&alpha;-methylacyl-CoA racemase
PWY-6947 	6.2.1.41	3-[(3aS,4S,7aS)-7a-methyl-1,5-dioxo-octahydro-1H-inden-4-yl]propanoyl:CoA ligase
PWY-6947 	6.2.1.42	3-oxocholest-4-en-26-oate--CoA ligase
PWY-6948	1.1.1.145	sitosterol dehydrogenase
PWY-6948	6.4.1	MONOMER-17007
PWY-695	1.1.1.288	AT1G52340-MONOMER
PWY-695	1.13.11.51	AT1G78390-MONOMER
PWY-695	1.13.11.51	AT3G14440-MONOMER
PWY-695	1.13.11.51	AT3G24220-MONOMER
PWY-695	1.13.11.51	AT4G18350-MONOMER
PWY-6954 	1.13.11.56	1,2-dihydroxynaphthalene dioxygenase subunit
PWY-6954 	1.13.11	protocatechuate 2,3-dioxygenase
PWY-6954 	1.14.12.12 	naphthalene 1,2-dioxygenase complex
PWY-6954 	1.14.12.12	naphthalene 1,2-dioxygenase complex
PWY-6954 	1.14.12.1	anthranilate dioxygenase oxygenase component 
PWY-6954 	1.14.13.2	MONOMER-11506
PWY-6954 	1.14.13.2	MONOMER-11534
PWY-6954 	1.1.5	MONOMER-11528
PWY-6954 	1.2.1.64	MONOMER-11532
PWY-6954 	1.2.1.65	MONOMER-12810
PWY-6954 	1.2.1.85	2-hydroxymuconate-6-semialdehyde dehydrogenase
PWY-6954 	1.2.1.85	2-hydroxymuconate semialdehyde dehydrogenase
PWY-6954 	1.2.1.96	MONOMER-11531
PWY-6954 	1.3.1.29	cis-1,2-dihydro-1,2-dihydroxynaphthalene-1, 2-dehydrogenase subunit
PWY-6954 	4.1.1.77	4-oxalocrotonate decarboxylase
PWY-6954 	4.1.1.7	MONOMER-11530
PWY-6954 	4.1.2.45	MONOMER-12809
PWY-6954 	4.2.1.10	catabolic 3-dehydroquinate dehydratase
PWY-6954 	4.2.1.118	MONOMER-37
PWY-6954 	5.1.2.2	MONOMER-11527
PWY-6954 	5.3.2.6	2-hydroxymuconate tautomerase
PWY-6954 	5.3.2.6	4-oxalocrotonate tautomerase
PWY-6954 	5.99.1.4	2-hydroxychromene-2-carboxylate isomerase
PWY-6955	1.13.11.30	L-dopa 2,3-dioxygenase
PWY-6955	2.1.1	N-demethyllincomycin methyltransferase
PWY-6956 	1.14.13.1	MONOMER-12811
PWY-6956 	1.14.13.1	MONOMER-14356
PWY-6957 	1.13.11.1	subunit of catechol 1,2-dioxygenase
PWY-6957 	2.3.1.174 	&beta;-ketoadipyl CoA thiolase
PWY-6957 	2.8.3.6	&alpha; subunit of &beta;-ketoadipate succinyl-CoA transferase 
PWY-6957 	3.1.1.24	subunit of 3-oxoadipate enol-lactone hydrolase
PWY-6957 	4.1.3.39	4-hydroxy-2-oxovalerate aldolase
PWY-6957 	5.3.3.4	muconolactone isomerase subunit
PWY-6957 	5.5.1.1	muconate lactonizing enzyme subunit
PWY-6958 	1.14.19.30	acyl-lipid 5-desaturase
PWY-6958 	1.14.19.47	acyl-lipid 6-desaturase
PWY-6960	1.11.1.11	<small>L</small>-ascorbate peroxidase
PWY-6961	1.1.1.130	2,3-diketo-L-gulonate reductase
PWY-6961	4.1.1.85	3-keto-L-gulonate 6-phosphate decarboxylase
PWY-6962 	1.14.13	dimethylamine monooxygenase
PWY-6962 	1.14.13.M46 	trimethylamine monooxygenase
PWY-6962 	1.5.99.5	N-methyl glutamate dehydrogenase
PWY-6962 	1.5.99.5	methylglutamate dehydrogenase
PWY-6962 	2.1.1.21	N-methyl-L-glutamate synthase
PWY-6962 	4.1.2.32	trimethylamine-oxide aldolase
PWY-6962 	6.3.4.12 	glutamate--methylamine ligase
PWY-6962 	6.3.4.12	glutamate--methylamine ligase
PWY-6963 	6.3.1.2	glutamine synthetase, type II
PWY-6964	1.4.7.1 	glutamate synthase (ferredoxin)
PWY-6968	4.1.2.32	trimethylamine-oxide aldolase
PWY-6969	1.1.1.37	G185E-5411-MONOMER
PWY-6969	1.1.1.42	isocitrate dehydrogenase subunit
PWY-6969	1.2.7.3	2-oxoglutarate:ferredoxin oxidoreductase
PWY-6969	1.3.5.1	G185E-7593-MONOMER
PWY-6969	2.3.3.16 	G185E-5050-MONOMER
PWY-6969	2.3.3.9	MONOMER-11953
PWY-6969	4.1.3.1	isocitrate lyase subunit
PWY-6969	4.2.1.2	G185E-5263-MONOMER
PWY-6969	6.2.1.5	Succinyl-CoA ligase [ADP-forming]
PWY-6970	1.2.1.51	pyruvate dehydrogenase (NADP<sup>+</sup>)
PWY-6971	1.14.13	8,8a-deoxyoleandolide monooxygenase
PWY-6971	2.1.1.239	L-olivosyl-oleandolide 3-O-methyltransferase
PWY-6971	2.3.1	8,8a-deoxyoleandolide synthase
PWY-6971	2.4.1	L-oleandrosyl-oleandolide desosaminyltransferase
PWY-6971	2.4.1	oleandolide olivosyltransferase
PWY-6972	2.4.1	oleandomycin glycosyltransferase
PWY-6972	3.2.1	glucosyl-oleandomycin glucohydrolase
PWY-6977 	1.14.13.154	erythromycin C-12 monooxygenase
PWY-6977 	1.14.13.188	6-deoxyerythronolide B 6S-monooxygenase
PWY-6977 	2.1.1.234	dTDP-3-amino-3,4,6-trideoxy-&alpha;-<small>D</small>-glucopyranose N,N-dimethyltransferase
PWY-6977 	2.1.1.254	erythromycin 3-o-methyltransferase
PWY-6977 	2.3.1.94	6-deoxyerythronolide B synthase monomer
PWY-6977 	2.4.1.278	3-L-mycarosyl erythronolide B desosaminyltransferase complex
PWY-6977 	2.4.1.328	dTDP-L-mycarosyl: erythronolide B mycarosyltransferase
PWY-6977 	2.6.1.106	dTDP-3-oxo-3,4,6-trideoxy-&alpha;-D-glucopyranose transaminase
PWY-6977 	2.6.1.33	dTDP-4-amino-4,6-dideoxy-<small>D</small>-glucose transaminase
PWY-6977 	4.3.1.30	dTDP-4-amino-4,6-dideoxy-D-glucose deaminase
PWY-6978	2.5.1.39	MONOMER-17095
PWY-6978	4.1.1.98	4-hydroxy-3-polyprenylbenzoate decarboxylase
PWY-6981 	2.3.1.4	glucosamine-6-phosphate N-acetyltransferase 1
PWY-6981 	2.3.1.4	MONOMER-17117
PWY-6981	2.4.1.16	chitin synthase 2
PWY-6981 	2.4.1.1	glycogen phosphorylase
PWY-6981 	2.4.1.1	glycogen phosphorylase b
PWY-6981 	2.4.1.1	glycogen phosphorylase <I>b
PWY-6981 	2.6.1.16	MONOMER-17113
PWY-6981 	2.7.1.1 	YCL040W-MONOMER
PWY-6981 	2.7.7.23	MONOMER-17119
PWY-6981 	2.7.7.23	YDL103C-MONOMER
PWY-6981 	3.2.1.28	acid trehalase
PWY-6981 	3.2.1.28	neutral trehalase
PWY-6981 	3.2.1.28	trehalase (bound form)
PWY-6981 	3.2.1.33 	glycogen debranching enzyme
PWY-6981 	3.2.1.3	glucan 1,4-&alpha;-glucosidase
PWY-6981 	5.3.1.9	phosphoglucoisomerase subunit
PWY-6981 	5.4.2.2	phosphoglucomutase-1
PWY-6981 	5.4.2.2	phosphoglucomutase-2
PWY-6981 	5.4.2.3	MONOMER-17118
PWY-6981 	5.4.2.3	YEL058W-MONOMER
PWY-6983	1.1.1.325	sepiapterin reductase (L-threo-7,8-dihydrobiopterin forming)
PWY-6984	2.3.1.200	lipoyl amidotransferase
PWY-6984	3.5.1	MONOMER-17108
PWY-6984	6.3.1.20	lipoate&mdash;protein ligase
PWY-6986	1.1.1.126	(4S,5S)-4,5-dihydroxy-2,6-dioxohexanoate reductase
PWY-6986	4.2.2.11	alginate lyase II
PWY-6986	4.2.2.26	oligoalginate lyase
PWY-6986	4.2.2.3	alginate lyase III
PWY-6987	2.3.1.181	octanoyl transferase
PWY-6987	2.3.1.204	GcvH:[lipoyl domain] amidotransferase
PWY-6987	2.8.1.8	lipoyl synthase (lipoic acid synthetase)
PWY-6989	1.1.1	3-exo-hydroxycamphor dehydrogenase
PWY-6989	1.5.1.42 	3,6-diketocamphane 1,6-monooxygenase complex
PWY-6990	1.1.1 	(+)-borneol dehydrogenase
PWY-6990	3.1.3	(+)-bornyl-diphosphate diphosphatase
PWY-6991	1.1.1.227	(-)-borneol dehydrogenase
PWY-699	1.1.1 	brassinolide synthase
PWY-699 	1.14.99 	AT3G13730-MONOMER
PWY-6991	3.1.3	(-)-bornyl diphosphate hydrolase
PWY-6991	5.5.1.22	(-)-bornyl diphosphate synthase
PWY-6992	1.1.1.292	1,5-anhydro-D-fructose reductase
PWY-6993	1.13.11.9	2,5-dihydroxypyridine dioxygenase
PWY-6993	1.14.13.163	6-hydroxy-3-succinoyl-pyridine hydroxylase
PWY-6993	1.2.1.83	3-succinoylsemialdehyde-pyridine dehydrogenase
PWY-6993	1.3.99 	MONOMER-17156
PWY-6993	1.3.99	nicotine dehydrogenase
PWY-6993	1.4.3.24	pesudooxynicotine oxidase
PWY-6993	1.5.99	3-succinoylpyridine monooxygenase complex
PWY-6993	3.5.1.106	N-formylmaleamate deformylase
PWY-6993	3.5.1.107	maleamate amidohydrolase
PWY-6993	5.2.1.1	maleate isomerase
PWY-6994	1.4.1	pyrrolysine synthase
PWY-6994	5.4.99.58	methylornithine synthase
PWY-6995	1.1.3	5-hydroxymethylfurfural oxidase
PWY-6996	1.3.1	daidzein reductase
PWY-6996	1.3.1	dihydrodaidzein reductase
PWY-6997 	1.3.99.8	furoyl-CoA dehydrogenase
PWY-6997 	6.2.1.31 	MONOMER-17170
PWY-6998	2.7.7	MONOMER-17174
PWY-6999 	1.5.1.39 	methylxanthine N1-demethylase
PWY-6999 	1.5.1.39 	methylxanthine N3-demethylase
PWY-7000	1.1.1.355	2-oxokanamycin reductase
PWY-7000	1.1.3 	paromamine 6-oxidase
PWY-7000	1.14.11.37	kanamycin B dioxygenase
PWY-7000	2.4.1.301	paromamine glycosyltransferase
PWY-7000	2.6.1.94 	glutamate--6-dehydroparomanine aminotransferase
PWY-7002	1.13.11.66	hydroquinone 1,2-dioxygenase
PWY-7002	1.14.13.84	4-hydroxyacetophenone monooxygenase subunit
PWY-7002	1.14.13.84	MONOMER-17196
PWY-7002	3.1.1	MONOMER-17199
PWY-7006	1.13.11	4-amino-3-hydroxy-benzoate 2,3-dioxygenase subunit
PWY-7006	1.2.1.85	MONOMER-17212
PWY-7006	3.5.99	2-amino-5-carboxymuconic 6-semialdehyde deaminase subunit
PWY-7006	4.1.1.77	MONOMER-17214
PWY-7006	5.3.2.6	4-oxalocrotonate tautomerase
PWY-7007	1.3.3.6	acyl coA dehydrogenase
PWY-7007 	3.1.2.20 	acyl-CoA thioesterase
PWY-7007 	3.1.2.20 	thioesterase II
PWY-7007 	4.1.1.56 	methylketone synthase 1
PWY-7009	1.13.11.35 	MONOMER-17253
PWY-7009 	1.14.13 	salicylate monooxygenase
PWY-7010 	1.14.13 	2-hydroxybiphenyl-3-monooxygenase [multifunctional]
PWY-7011 	1.13.11.M5 	2,3-dihydroxybiphenyl-1,2-dioxygenase [multifunctional]
PWY-7011 	3.7.1 	2-hydroxy-6-oxo-6-phenyl-2,4-hexadienoate hydrolase [multifunctional]
PWY-7013	1.1.1.1	propanol dehydrogenase
PWY-7013 	2.3.1.8 	phosphate propanoyltransferase
PWY-7013 	2.7.2.15	MONOMER-17272
PWY-7013	4.2.1.28	propanediol dehydratase small subunit 
PWY-7014 	1.1.1.329	2-deoxy-scyllo-inosamine dehydrogenase
PWY-7014 	2.4.1.283 	2-deoxystreptamine glucosyltransferase
PWY-7014 	2.6.1.100 	L-glutamine:2-deoxy-scyllo-inosose aminotransferase
PWY-7014 	3.5.1.112	2-N-acetylparomamine deacetylase
PWY-7014 	4.2.3.124	2-deoxy-scyllo-inosose synthase
PWY-7016 	1.1.3.44 	paromamine 6-oxidase
PWY-7016 	2.6.1.95 	glutamate--6-dehydroparomanine aminotransferase
PWY-7016 	3.5.1.113 	MONOMER-17259
PWY-7018	1.1.3	6-hydroxyparomomycin oxidase
PWY-7018	2.4.1	UDP-GlcNAc:ribosylparomamine N-acetylglucosaminyltransferase
PWY-7018	2.4.2	paromamine phosphoribosyltransferase
PWY-7018	2.6.1	glutamate--6-oxoparomomycin aminotransferase
PWY-7018	3.1.3	5-phosphoribosylparomamine phosphatase
PWY-7018	3.5.1	2-acetyl-6-hydroxyparomomycin deacetylase
PWY-7020 	1.1.3.43	paromamine 6-oxidase
PWY-7020 	1.14.14.13	4-(L-&gamma;-glutamylamino)butanoyl-[BtrI acyl-carrier protein] monooxygenase
PWY-7020 	1.1.99.38	2-deoxy-scyllo-inosamine dehydrogenase (SAM-dependent)
PWY-7020 	2.3.2.19	ribostamycin:(S)-2-hydroxy-4-glutamyl-4-aminobutanoyl-[BtrI acyl-carrier protein] (S)-2-hydroxy-4-glutamyl-4-aminobutanoate transferase
PWY-7020 	2.4.1.283	2-deoxystreptamine N-acetyl-D-glucosaminyltransferase
PWY-7020 	2.4.2.49	neamine phosphoribosyltransferase
PWY-7020 	2.6.1.101 	L-glutamine:2-deoxy-scyllo-inosose aminotransferase
PWY-7020 	2.6.1.93	glutamate--6-dehydroparomanine aminotransferase
PWY-7020 	3.1.3.88	5-phosphoribostamycin phosphatase
PWY-7020 	3.5.1.112	2-N-acetylparomamine deacetylase
PWY-7020 	4.1.1.95	L-glutamyl-[BtrI acyl-carrier protein] decarboxylase
PWY-7020 	4.2.3.124	2-deoxy-scyllo-inosose synthase
PWY-7020 	4.3.2.6	&gamma;-L-glutamyl-butirosin B &gamma;-glutamyl cyclotransferase
PWY-7020 	6.2.1.39	[BtrI acyl-carrier protein]--L-glutamate ligase
PWY-7021 	1.1.1.329	2-deoxy-scyllo-inosamine dehydrogenase
PWY-7021 	2.4.1.283	2-deoxystreptamine N-acetyl-D-glucosaminyltransferase
PWY-7021 	2.4.1.285	UDP-GlcNAc:ribostamycin N-acetylglucosaminyltransferase
PWY-7021 	2.4.2.49	neamine phosphoribosyltransferase
PWY-7021 	2.6.1.101 	L-glutamine:2-deoxy-scyllo-inosose aminotransferase
PWY-7021 	3.1.3.88	5-phosphoribostamycin phosphatase
PWY-7021 	4.2.3.124	MONOMER-17230
PWY-7024 	1.1.1.298 	malonyl CoA reductase
PWY-7024 	2.8.3.22	succinyl CoA:L-malate CoA transferase
PWY-7024 	4.2.1.148	&beta;-methylmalyl-CoA dehydratase
PWY-7024 	4.2.1.153	mesaconyl-C4-CoA hydratase
PWY-7024 	5.1.99.1	methylmalonyl-CoA epimerase
PWY-7024 	5.4.1.3	mesaconyl-C1-CoA-C4-CoA transferase
PWY-7024 	5.4.99.2	methylmalonyl-CoA mutase
PWY-7024 	6.2.1.36 	propionyl-CoA synthase
PWY-7024 	6.4.1.2	acetyl-CoA carboxylase complex
PWY-7024 	6.4.1.3	propionyl-CoA carboxylase
PWY-7025	2.1.1	MONOMER-17296
PWY-7025	2.4.2	paromamine D-xylosyltransferase
PWY-7027	1.1.1	MONOMER-17300
PWY-7027	3.1.2 	OleA
PWY-7028	2.3.1.203	UDP-4-amino-4,6-dideoxy-&alpha;-D-N-acetyl-D-glucosamine N-acetyltransferase subunit
PWY-7028	2.6.1.34	MONOMER-17319
PWY-7028	2.6.1.34	MONOMER-17337
PWY-7028	4.2.1.135	MONOMER-17318
PWY-7028	4.2.1.135	MONOMER-17336
PWY-7029	2.3.1	fatty-acyl-[acp] elongase/decarboxylase
PWY-7031	2.4.1.290	N,N-diacetylbacillosaminyl-diphospho-undecaprenol &alpha;-1,3-N-acetylgalactosaminyltransferase
PWY-7031	2.4.1.291	N-acetylgalactosamine-N,N-diacetylbacillosaminyl-diphospho-undecaprenol 4-&alpha;-N-acetylgalactosaminyltransferase
PWY-7031	2.4.1.292	GalNAc-&alpha;-(1&rarr;4)-GalNAc-&alpha;-(1&rarr;3)-diNAcBac-PP-undecaprenol &alpha;-1,4-N-acetyl-D-galactosaminyltransferase
PWY-7031	2.4.1.293	GalNAc5-diNAcBac-PP-undecaprenol &beta;-1,3-glucosyltransferase
PWY-7031	2.4.99.19	undecaprenyl-diphosphooligosaccharide--protein glycotransferase
PWY-7031	2.7.8.36	MONOMER-17321
PWY-7032 	1.2.1.80	acyl-[acp] reductase
PWY-7032 	4.1.99.5	aldehyde decarbonylase
PWY-7033	1.2.1	acyl-CoA reductase (NADH)
PWY-7033	4.1.99	aldehyde decarbonylase
PWY-7033	4.1.99	MONOMER-17303
PWY-7033	6.2.1.3 	MONOMER-17305
PWY-7035	1.14.14 	aldehyde decarbonylase
PWY-7036 	1.1.1.330	very-long-chain 3-oxoacyl-CoA reductase
PWY-7036 	1.1.1.35 	peroxisomal bifunctional enzyme
PWY-7036 	1.3.1.93	very-long-chain enoyl-CoA reductase
PWY-7036 	2.3.1.199	very-long-chain 3-ketoacyl-CoA synthase
PWY-7036	2.3.1.199	very-long-chain 3-oxoacyl-CoA synthase
PWY-7036 	4.2.1.134	very-long-chain 3-hydroxyacyl-CoA dehydratase
PWY-7036 	4.2.1.17 	mitochondrial enoyl-CoA hydratase
PWY-7036 	4.2.1.17 	peroxisomal fatty acid &beta;-oxidation multifunctional protein AIM1
PWY-7037	2.4.1	bacterial oligosaccharyltransferase PglO
PWY-7037	2.4.1	MONOMER-17339
PWY-7037	2.4.1	MONOMER-17340
PWY-7037	2.4.1	MONOMER-17341
PWY-7037 	2.7.8.36 	bifunctional acetyltransferase / phosphoglycosyltransferas PglB
PWY-7039	2.7.1.107	AT4G30340-MONOMER
PWY-7039	2.7.1.107	AT5G63770-MONOMER
PWY-7039	2.7.1.107	MONOMER-17355
PWY-7039	2.7.4	MONOMER-17354
PWY-7039	3.1.4.11	MONOMER-17353
PWY-7040	1.4.3	MONOMER-17361
PWY-7042 	1.1.1	11-cis-3-hydroxyretinol dehydrogenase
PWY-7042 	1.1.1	all-trans-3-hydroxyretinal dehydrogenase
PWY-7043	1.13.11.65	carotenoid isomerooxygenase
PWY-7043	5.1.2	MONOMER-17371
PWY-7044	1.13.11.64	5-nitrosalicylate 1,2-dioxygenase
PWY-7044	3.1.1.91	2-oxo-3-(5-oxofuran-2-ylidene)propanoate lactonase
PWY-7044	3.5.99.8	5-nitroanthranilate aminohydrolase
PWY-7044	3.7.1.20	3-fumarylpyruvate hydrolase
PWY-7045	1.1.1	mithramycin dehydrogenase
PWY-7045	1.1.1	MtmTII
PWY-7045	1.13.99	MtMOI
PWY-7045	1.14.13	premithramycin B monooxygenase
PWY-7045	2.1.1	demethylpremithracinone o-methyltransferase
PWY-7045	2.1.1	premithramycin A3 methyltransferase
PWY-7045	2.3.1	mithramycin polyketide synthase
PWY-7045	2.4.1	3A-deolivosylpremithramycin B:dTDP-D-olivose D-olivosyltransferase
PWY-7045	2.4.1	MtmGIII
PWY-7045	2.4.1	MtmGIV
PWY-7045	2.4.1	premithramycin A3:dTDP-D-olivose D-olivosyltransferase
PWY-7045	4.2.1	3,5,7,9,11,13,15,17,19-nonaoxoicosanoate cyclase/dehydratase
PWY-7046	1.3.7.9	4-hydroxybenzoyl-CoA reductase HbaD subunit 
PWY-7046	4.1.2 	4-coumaroyl-CoA hydratase/aldolase
PWY-7046	6.2.1.12	MONOMER-17392
PWY-7046	6.2.1.25 	4-hydroxybenzoate-CoA ligase / benzoate-CoA ligase
PWY-7049 	6.2.1.3 	long-chain-fatty-acid--CoA ligase 1
PWY-7050	2.3.1	icosapentaenoate synthase
PWY-7052	3.4.15.6	cyanophycinase
PWY-7052	6.3.2.30 	cyanophycin synthase
PWY-7053 	1.1.1.107 	3-oxoacyl-[acyl-carrier-protein] reductase
PWY-7053	1.14.19.31	acyl-lipid 4-desaturase
PWY-7053	2.3.1	C20 &Delta;5-polyunsaturated fatty acyl-CoA elongase
PWY-7055	2.1.1	MONOMER-17439
PWY-7055	2.4.1	MONOMER-17432
PWY-7056	3.2.1	MONOMER-17430
PWY-7057	2.4.1.104	MONOMER-17442
PWY-7057 	2.4.1.104	UDPG:o-dihydroxycoumarin 7-O-glucosyltransferase
PWY-7057	3.2.1	MONOMER-17441
PWY-7058 	2.1.1	esculetin-O-methyltransferase
PWY-7058	2.1.1	MONOMER-17459
PWY-7058 	2.4.1.126	hydroxycinnamate-4-O-&beta;-glucosyltransferase
PWY-7058 	2.4.1 	flavonol 3-O-glucosyltransferase
PWY-7058	2.4.1	MONOMER-17453
PWY-7059	1.3.1.100	chanoclavine-I aldehyde dehydrogenase
PWY-7059	1.5.1.44	festuclavine synthase
PWY-7059	2.3.1.205	fumigaclavine C acetyltransferase
PWY-7059	2.5.1.100	fumigaclavine A dimethylallyltransferase
PWY-7064	5.4.3.10	MONOMER-17466
PWY-7064	5.4.3.10	Phenylalanine aminomutase
PWY-7065	1.14.13 	MONOMER-17468
PWY-7065	1.14.13	MONOMER-17469
PWY-7066	1.14.13.134 	11&alpha;-hydroxy-&beta;-amyrin dehydrogenase
PWY-7066	1.14.13.173 	11-oxo-&beta;-amyrin 30-oxidase
PWY-7067 	1.14.13 	&alpha;/&beta;-amyrin 28-monooxygenase
PWY-7068 	1.14.13 	&beta;-amyrin 28-oxidase
PWY-7068 	1.14.13 	MONOMER-17476
PWY-7069	1.14.13.201	MONOMER-17477
PWY-7069 	5.4.99.39 	&alpha;/&beta;-amyrin synthase
PWY-7070	1.14.13	ent-kaurenoic acid NADPH:oxygen oxidoreductase
PWY-7071	2.4.1	MONOMER-17484
PWY-7071	2.4.1	UDPG glucosyltransferase 76G1
PWY-7071	2.4.1	UDPG glucosyltransferase 85C2
PWY-7072	1.17.98.b	hydroxysqualene dehydroxylase
PWY-7072	2.1.1	MONOMER-17493
PWY-7072	2.5.1.103	presqualene diphosphate synthase
PWY-7072	4.2.1.129 	squalene hopene cyclase
PWY-7072	4.2.3.ea	hydroxysqualene synthase
PWY-7074	2.4.1	MONOMER-17497
PWY-7074	3.2.1	2-phenylethyl &beta;-D-glucopyranoside glucosidase; subunit 
PWY-7075	2.3.1.224	MONOMER-17504
PWY-7075	2.3.1.224	MONOMER-17505
PWY-7076 	2.1.1 	orcinol O-methyltransferase
PWY-7077	2.7.1	MONOMER-17510
PWY-7077	3.5.1.25	MONOMER-17511
PWY-7077	3.5.99	D-galactosamine-6-phosphate deaminase/isomerase
PWY-7079	1.14.13.43	questin monooxygenase
PWY-7079	1	dihydrogeodin oxidase subunit
PWY-7079 	1 	silochrin oxidase subunit
PWY-7079	2.1.1.283	EOMT
PWY-7080	1.21.3.4	sulochrin oxidase, subunit 
PWY-7080	1.21.3.5	sulochrin oxidase subunit
PWY-7081	1.13.11.37	MONOMER-17532
PWY-7082	1.10.2.2	cytochrome cm552
PWY-7082	1.7.2.5	cytochrome c554
PWY-7084 	1.14.99.39 	particulate methane monooxygenase hydroxylase component
PWY-7084 	1.14.99.39	soluble ammonia monooxygenase
PWY-7084	1.7.2.1	nitrite reductase
PWY-7084	1.7.2.1	nitrite reductase (NO-forming)
PWY-7084	1.7.2.5	cytochrome c-&beta;
PWY-7084	1.7.2.5	nitric oxide reductase
PWY-7084 	1.7.2.6	hydroxylamine oxidoreductase
PWY-7085	1.14.13	triethylamine monooxygenase [multifunctional]
PWY-7087	2.1.1.255	MONOMER-17559
PWY-7087	2.1.1.255	MONOMER-17561
PWY-7087	4.2.3.118	2-methylisoborneol synthase subunit
PWY-7087	4.2.3.118	MONOMER-17562
PWY-7088	1.14.13.41	tyrosine N-monooxygenase
PWY-7089	3.2.1.21	taxiphyllin &beta;-glucosidase
PWY-7091	3.2.1	MONOMER-17612
PWY-7092 	3.2.1.21	&beta;-diglucosidase
PWY-7092 	3.2.1.21	&beta;-glucosidase
PWY-7092 	4.1.2.46	(S)-hydroxynitrile lyase
PWY-7093	3.2.1.119	MONOMER-17571
PWY-7093	3.2.1.119	MONOMER-17579
PWY-7093	4.1.2.10	(R)-mandelonitrile lyase subunit
PWY-7094	1.1.1.35 	fatty acid oxidation complex &alpha; subunit dimer
PWY-7094	1.3.8.7	medium-chain acyl-CoA dehydrogenase
PWY-7094	2.3.1.16	fatty acid oxidation complex &beta; subunit dimer
PWY-7094	2.3.1.207	&beta;-ketodecanoyl-[acyl-carrier-protein] synthase
PWY-7094	6.2.1.3	long-chain-fatty-acid--CoA ligase
PWY-7094	6.2.1.3	medium-chain-fatty-acid--CoA ligase
PWY-7095 	2.4.1.178	MONOMER-17558
PWY-7096	1.3.1.9 	enoyl-[acyl-carrier-protein] reductase (NADH)
PWY-7096	1.3.1.M3 	enoyl-[acyl-carrier-protein] reductase (NADH)
PWY-7098	1.14.13.82	vanillate O-demethylase oxygenase component 
PWY-7098	1.2.1.67	vanillin dehydrogenase
PWY-7100	2.4.1	MONOMER-17606
PWY-7101	1.13.11.68	9-cis-&beta;-carotene 9,10-cleavage dioxygenase
PWY-7101 	1.13.11.69 	10-apo-&beta;-carotenal 13,14-cleaving dioxygenase
PWY-7101	5.2.1.14	&beta;-carotene isomerase
PWY-7102 	4.2.3.133 	&delta;-selinene synthase
PWY-7102 	4.2.3 	&gamma;-humulene synthase
PWY-7105 	4.4.1 	tetraketide synthase/olivetol synthase
PWY-7110 	1.1.1	dTDP-3-N,N-dimethylamino-4-oxo-2,3,6-trideoxy-&alpha;-L-glucose 4-ketoreductase
PWY-7110 	1.14.13.154	erythromycin D C-12 hydroxylase
PWY-7110 	2.4.1	dTDP-L-megosamine:erythromycin C L-megosaminyltransferase
PWY-7110 	5.3.1	dTDP-3-N,N-dimethylamino-4-oxo-2,3,6-trideoxy-&alpha;-D-glucose 5-epimerase
PWY-7111 	1.1.1.86 	acetohydroxy acid isomeroreductase
PWY-7111 	2.2.1.6	acetolactate synthase / acetohydroxybutanoate synthase, catalytic subunit 
PWY-7111 	4.2.1.9	dihydroxy acid dehydratase
PWY-7112 	2.3.1.80	cysteine-S-conjugate N-acetyltransferase
PWY-7112 	2.3.1.80	MONOMER-10112
PWY-7112 	2.5.1.18	glutathione S-transferase A1
PWY-7112 	3.4.11.2 	aminopeptidase M
PWY-7112 	3.4.11.2 	cysteinylglycinase / cytosolic leucyl aminopeptidase
PWY-7112 	3.4.19.9	&gamma;-glutamyl hydrolase
PWY-7112 	3.4.19.9	pteroyl-&gamma;-glutamyl hydrolase
PWY-7113	3.2.1.161	MONOMER-17662
PWY-7114	3.2.1.149	&beta;-primerverosidase (Shuixian)
PWY-7114	3.2.1.149	&beta;-primeverosidase (Assam)
PWY-7114	3.2.1.149	&beta;-primeverosidase (Yabukita)
PWY-7115	1.1.1.37	MONOMER-17671
PWY-7115	1.1.1.38 	MONOMER-17672
PWY-7115	1.1.1.38 	MONOMER-17673
PWY-7115	2.6.1.1	MONOMER-17669
PWY-7115	2.6.1.1	MONOMER-17670
PWY-7115	2.6.1.2	MONOMER-17674
PWY-7115	4.1.1.31	MONOMER-17666
PWY-7115	4.1.1.31	MONOMER-17667
PWY-7117	1.1.1.82	MONOMER-17678
PWY-7117	2.6.1.2	MONOMER-17680
PWY-7117	2.7.9.1	MONOMER-17681
PWY-7117	4.1.1.31	MONOMER-17675
PWY-7117	4.1.1.49	phosphoenolpyruvate carboxykinase (ATP)
PWY-7117	4.2.1.1	&beta;-carbonic anhydrase
PWY-7118 	1.1.1.38 	NAD-dependent malic enzyme, mitochondrial
PWY-7118	3.5.1.41	chitin deacetylase 1
PWY-7118	3.5.1.41	chitin deacetylase 2
PWY-7119	2.7.1.91	sphingoid long chain base kinase 4
PWY-7119	2.7.1.91	sphingoid long chain base kinase 5
PWY-7119	3.1.3	sphingoid base phosphate phosphatase
PWY-7119	3.1.4	inositol phosphosphingolipids phospholipase C
PWY-7119	3.5.1.23	dihydroceramidase
PWY-7119	3.5.1.23	phytoceramidase
PWY-7119	4.1.2.27	sphingoid base phosphate lyase
PWY-7120	2.3.1.188	AT5G41040-MONOMER
PWY-7124 	1.10.3.9	photosystem II
PWY-7124 	1.10.3.9	photosystem II monomer
PWY-7124 	1.10.9.1	cytochrome b6f complex
PWY-7124 	1.14.11.41 	2-oxoglutarate-dependent ethylene/succinate-forming enzyme
PWY-7124 	1.18.1.2 	ferredoxin-NADP oxidoreductase
PWY-7124 	1.2.1.59 	glyceraldehyde 3-phosphate dehydrogenase 2
PWY-7124 	1.97.1.12	photosystem I
PWY-7124 	1.97.1.12	photosystem I complex
PWY-7124 	2.2.1.1	transketolase
PWY-7124 	2.3.3.16 	citrate synthase
PWY-7124 	2.7.1.19	PRKSYN-MONOMER
PWY-7124 	2.7.2.3	PGKSYN-MONOMER
PWY-7124 	3.1.3.11 	D-fructose 1,6-bisphosphatase class 2/sedoheptulose 1,7-bisphosphatase
PWY-7124 	3.1.3.11	fructose-1,6-bisphosphatase class 1
PWY-7124 	4.1.1.39 	ribulose bisphosphate carboxylase/oxygenase
PWY-7124 	4.1.2.13 	FBAASYN-MONOMER
PWY-7124 	4.1.2.13 	FBABSYN-MONOMER
PWY-7124 	4.2.1.11	enolase
PWY-7124 	5.1.3.1	RPESYN-MONOMER
PWY-7124 	5.3.1.1	TpiA
PWY-7124 	5.3.1.6	RPIASYN-MONOMER
PWY-7127	2.7.7	MONOMER-17706
PWY-7128	1.13.11.9	2,5-dihydroxypyridine dioxygenase
PWY-7128	1.14.13.163	6-hydroxy-3-succinoyl-pyridine hydroxylase
PWY-7128	1.4.3.M1	6-hydroxypseudooxynicotine oxidase
PWY-7128	1.5.99.4	nicotine dehydrogenase
PWY-7128	3.5.1.106	N-formylmaleamate deformylase
PWY-7128	3.5.1.107	maleamate amidohydrolase
PWY-7128	5.2.1.1	maleate isomerase
PWY-7129	2.4.1.91	MONOMER-17718
PWY-7129	2.4.1	MONOMER-17717
PWY-7129	2.4.1	MONOMER-17720
PWY-7129	2.4.1	quercetin-7-O-glucosyltransferase
PWY-7129	2.4.1 	quercetin-O-glucosyltransferase
PWY-7130	1.1.1.370 	scyllo-inositol 2-dehydrogenase
PWY-7130	1.1.1	D-idonate 5-dehydrogenase
PWY-7130	1.1.1	L-gluconate 5-dehydrogenase
PWY-7130	2.7.1.178 	2-dehydro-3-deoxygalactonokinase
PWY-7130	4.1.2.55 	2-dehydro-3-deoxy-6-phosphogalactonate aldolase
PWY-7130	4.2.1	D-idonate dehydratase
PWY-7131	2.5.1.101	N,N-diacetyllegionaminate synthase
PWY-7131	2.7.7.82	CMP-legionaminic acid synthase
PWY-7131	3.2.1.184	UDP-N,N-diacetylbacillosamine 2-epimerase
PWY-7133	3.2.1	quercetin-3-O-&beta;-glucosidase
PWY-7133	3.2.1	quercetin-4-O-&beta;-glucosidase
PWY-7134	3.2.1	flavonol 3-O-&beta;-heterodisaccharide glucosidase
PWY-7134	3.2.1	flavonol-3-O-glycosidase
PWY-7134	3.2.1	rutin hydrolase
PWY-7135	2.1.1	MONOMER-17755
PWY-7135	2.1.1	MONOMER-17757
PWY-7135	2.1.1	MONOMER-17758
PWY-7135	2.1.1	MONOMER-17759
PWY-7135	3.2.1	MONOMER-17756
PWY-7135	4.3.3.3	MONOMER-17752
PWY-7136	1.1.1.347	geraniol dehydrogenase
PWY-7136	1.2.1.86	geranial dehydrogenase
PWY-7136	5.4.4.4 	linalool dehydratase/isomerase
PWY-7138	1.1.1	noscapine synthase
PWY-7138	1.14.21.5	MONOMER-17768
PWY-7138	2.1.1.117	MONOMER-17764
PWY-7138	2.1.1.117	MONOMER-17766
PWY-7138 	2.1.1.122 	MONOMER-13847
PWY-7138	2.1.1.89 	MONOMER-17765
PWY-7138	2.1.1	MONOMER-17770
PWY-7139	2.4.1	lignan 1,6-glucosyltransferase
PWY-7139	2.4.1	lignan glucosyltransferase
PWY-7142	4.2.1.66	cyanide hydratase
PWY-7143 	2.4.1.91 	flavonol 3-O-glucosyltransferase
PWY-7143 	2.4.1	curcumin glucoside 1,6-glucosyltransferase
PWY-7145 	2.4.1	flavonoid glucoside 1,6-glucosyltransferase
PWY-7147	1.14.15.12	pimeloyl-[acp] synthase
PWY-7150 	2.1.1.82	MONOMER-14287
PWY-7150 	2.4.1	flavonol 2-O-glucosyltransferase
PWY-7151 	1.14.11	flavonol 6-hydroxylase
PWY-7151 	2.1.1 	3,7dimethylquercetin/quercetagetin 4-O-methyltransferase
PWY-7151 	2.1.1 	methylquercetagetin 6-O-methyltransferase
PWY-7151 	2.1.1	methylquercetagetin-O-methyltransferase
PWY-7151 	2.1.1	methyl-quercetin/quercetagetin 3/5-O-methyltransferase
PWY-7151 	2.1.1	methyl-quercetin/quercetagetin glucoside 2/3-Omethyltransferase
PWY-7151 	2.4.1	flavonol 3-O-glucosyltransferase
PWY-7152	1.14.19.12	acyl-lipid &omega;-13 desaturase
PWY-7153	1.10.3.15 	grixazone synthase
PWY-7153	4.1.2.56	2-amino-4,5-dihydroxy-6-one-heptanoic acid-7-phosphate synthase subunit
PWY-7153	4.1.99.20	MONOMER-17798
PWY-7156 	1.14.19.20	sterol C-5 desaturase
PWY-7157	1.14.13	flavonol 6-hydroxylase
PWY-7157	2.1.1	flavonol 6-O-methyltransferase
PWY-7157	2.1.1	flavonol 7-O-methyltransferase
PWY-7157 	2.4.1.237	flavonol 7-O-glucosyltransferase
PWY-7157	2.4.1	flavonol 3-O-glucosyltransferase
PWY-7157	2.4.1 	flavonol 7O-glucosyltransferase
PWY-7158	1.14.16.1	MONOMER-17816
PWY-7158	1.14.16.1	MONOMER-17817
PWY-7158	4.2.1.96	MONOMER-17818
PWY-7158	4.2.1.96	MONOMER-17819
PWY-7159 	1.14.13.81	magnesium-protoporphyrin IX monomethyl ester [oxidative] cyclase 1
PWY-7159 	1.14.13.81	magnesium-protoporphyrin IX monomethyl ester [oxidative] cyclase 2
PWY-7159	1.3.7.7	light-independent protochlorophyllide reductase
PWY-7161 	2.1.1 	flavonol 3/5-methyltransferase
PWY-7163 	2.1.1 	flavonol 3-O-methyltransferase
PWY-7163 	2.1.1 	flavonol 7/4-methyltransferase
PWY-7165 	1.1.1.274	MONOMER-17832
PWY-7165 	1.1.1.346	2,5-diketo-D-gluconate reductase A
PWY-7165 	1.1.5.2	quinoprotein glucose dehydrogenase
PWY-7165 	1.1.99.3	D-gluconate dehydrogenase dehydrogenase subunit 
PWY-7166	2.4.1	flavonoid 3-O-glucosyltransferase
PWY-7167	4.3.99.4	choline trimethylamine-lyase
PWY-7169 	6.2.1 	cinnamate:CoA ligase
PWY-7170	1.14.15.20	AT1G58300-MONOMER
PWY-7170	1.14.15.20	AT1G69720-MONOMER
PWY-7170	1.14.15.20	AT2G26670-MONOMER
PWY-7170	1.3.7.4	AT3G09150-MONOMER
PWY-7172 	2.3.1	flavonol 3-O-glucoside 3\-O-coumaroyltransferase"
PWY-7172 	2.3.1	flavonol 3-O-glucoside 4\-O-coumaroyltransferase"
PWY-7172 	2.3.1	flavonol 3-O-glucoside 6\-O-acyltransferase"
PWY-7173 	2.4.1.239	flavonol 3-O-glucoside glucosyltransferase
PWY-7173 	2.4.1.240	flavonol 3-O-diglucoside glucosyltransferase
PWY-7173 	2.4.1.91 	flavonol 3-O-glucosyltransferase
PWY-7174	2.5.1.49	O-acetylhomoserine aminocarboxypropyltransferase
PWY-7174	4.4.1	1-methylthio-xylulose 5-phosphate sulfo-lyase
PWY-7174	5.3.3	5-methylthioribulose-1-phosphate isomerase
PWY-7175	1.14.13	carotenoid 2,2-&beta;-ionone ring hydroxylase
PWY-7175 	1.14.13	carotenoid 2-hydroxylase
PWY-7175	1.14.13	zeaxanthin 2,2-&beta;-hydroxylase
PWY-7177 	3.6.1.5 	ectonucleoside triphosphate diphosphohydrolase 1
PWY-7177 	3.6.1.5 	ectonucleoside triphosphate diphosphohydrolase 2
PWY-7177	3.6.1.5 	nucleoside triphosphate phosphohydrolase
PWY-7177	6.3.4.2 	CTP synthetase
PWY-7178	1.1.1.175	MONOMER-13208
PWY-7179-1 	2.4.2.1 	purine nucleoside phosphorylase
PWY-7179-1 	3.5.4.4	adenosine deaminase
PWY-7181 	2.4.2.3 	pyrimidine-nucleoside phosphorylase
PWY-7181	2.4.2.3 	thymidine phosphorylase
PWY-7181 	2.4.2.3 	uridine phosphorylase 1
PWY-7182 	2.5.1.1	farnesyl pyrophosphate synthase
PWY-7182 	4.2.3.26	linalool synthase
PWY-7184 	2.7.4.6	nucleoside diphosphate kinase
PWY-7186 	2.1.1	MONOMER-17458
PWY-7186 	2.4.1.128 	UDP-glucose:scopoletin O-&beta;-D-glucosyltransferase
PWY-7189 	1.14.13	flavanone 2-hydroxylase
PWY-7189 	2.4.1	2-hydroxyflavanone dibenzoylmethane tautomer glucosyltransferase
PWY-7189 	4.2.1	2-hydroxyflavanone-C-glucoside dehydratase
PWY-7192 	2.4.1.234	flavonol 3-O-&beta;-D-galactosyltransferase
PWY-7192 	2.4.1	flavonol 3-O-&beta;-D-galactoside-2\-O-glucosyltransferase"
PWY-7193	2.7.1.aw 	uridine-cytidine kinase 1
PWY-7193	2.7.1.aw 	uridine-cytidine kinase 2
PWY-7194 	3.5.4.1	cytosine deaminase
PWY-7195 	3.2.2.3 	ribonucleoside hydrolase 2 (pyrimidine-specific)
PWY-7195 	3.2.2 	ribonucleoside hydrolase (pyrimidine-specific)
PWY-7195 	3.5.4.1 	creatinine deiminase subunit
PWY-7196 	2.7.1.48	AT5G40870-MONOMER
PWY-7196 	2.7.1.48 	UDK-MONOMER
PWY-7196	2.7.1.aw 	uridine kinase
PWY-7196 	2.7.4.14 	cytidylate kinase
PWY-7196 	3.2.2.3	ribonucleoside hydrolase 1 (pyrimidine-specific)
PWY-7196 	3.2.2.3 	ribonucleoside hydrolase 3
PWY-7196 	3.5.4.5	cytidine deaminase
PWY-7198	3.5.4.30	dCTP deaminase (dUMP-forming)
PWY-7200 	2.1.1.45	MONOMER-14571
PWY-7200 	2.7.1.145 	thymidine kinase
PWY-7200 	2.7.1.145 	thymidine kinase 1
PWY-7200 	2.7.1.145 	thymidine kinase 1a
PWY-7200 	2.7.1.145 	thymidine kinase 1b
PWY-7200 	2.7.1.74 	deoxyadenosine/deoxycytidine kinase
PWY-7200 	2.7.1.74 	thymidine kinase 2
PWY-7200 	2.7.4.25 	UMP/CMP kinase
PWY-7200 	2.7.4.6	nucleoside-diphosphate kinase 1
PWY-7200 	3.5.4.5	cytidine deaminase
PWY-7200 	3.5.4.5	cytidine/deoxycytidine deaminase
PWY-7204	1.1.1.65	pyridoxal reductase
PWY-7204	1.4.3.5	pyridoxine 5-phosphate oxidase
PWY-7204	1.4.3.5	pyridoxine/pyridoxamine 5-phosphate oxidase
PWY-7204	2.7.1.35	pyridoxal kinase
PWY-7205 	2.7.4.14 	uridylate kinase
PWY-7206	3.6.1 	5-hydroxy-CTP diphosphatase
PWY-7206	3.6.1.9 	dCTP diphosphatase 1
PWY-7208 	2.4.2.9	uracil phosphoribosyltransferase
PWY-7209 	1.3.1.2	dihydropyrimidine dehydrogenase
PWY-7209 	1.3.1.2	dihydropyrimidine dehydrogenase [NADP(+)]
PWY-7209 	3.5.1.6	&beta;-ureidopropionase
PWY-7209 	3.5.2.2	dihydropyrimidinase
PWY-7210	3.5.4.12	dCMP deaminase
PWY-7210	3.5.4.12	deoxycytidylate deaminase
PWY-7210	3.6.1.6 	ectonucleoside triphosphate diphosphohydrolase 4
PWY-7211 	1.3.5.2	dihydroorotate dehydrogenase
PWY-7211 	1.3.5.2	dihydroorotate dehydrogenase, type 2
PWY-721	1.14.12.16	MONOMER-2404
PWY-7211 	2.1.1.45	thymidylate synthase
PWY-7211 	2.1.3.2	aspartate carbamoyltransferase, catalytic subunit 
PWY-7211 	2.1.3.2	AT3G20330-MONOMER
PWY-7211 	2.4.2.10	orotate phosphoribosyltransferase
PWY-7211 	2.4.2.10 	UMP synthase
PWY-7211 	2.7.4.12 	thymidylate kinase
PWY-7211 	2.7.4.22 	UMP/CMP kinase
PWY-7211 	2.7.4.22 	UMP kinase
PWY-7211 	3.5.2.3	dihydroorotase
PWY-7211 	3.6.1.15 	cancer-related nucleoside-triphosphatase
PWY-7211 	3.6.1.15 	nucleoside-triphosphatase
PWY-7211 	3.6.1.15	nucleoside-triphosphatase
PWY-7211 	3.6.1.23 	deoxyuridine 5-triphosphate nucleotidohydrolase
PWY-7211 	3.6.1.23 	mitochondrial deoxyuridine 5-triphosphate nucleotidohydrolase
PWY-7211 	3.6.1.5 	apyrase
PWY-721	1.3.99.17	quinoline 2-oxidoreductase &alpha; subunit 
PWY-7211 	4.1.1.23	orotidine-5-phosphate decarboxylase
PWY-7212	3.2.1.167	MONOMER-17918
PWY-7213 	2.4.1.253	baicalein 7-O-glucuronosyltransferase
PWY-7213 	2.4.1	flavonoid 7-O-glucosyltransferase
PWY-7213 	3.2.1.167	baicalin-&beta;-D-glucuronidase
PWY-7214	1.11.1.7	flavonoid peroxidase 1
PWY-7214	1.11.1.7	flavonoid peroxidase 2
PWY-7214	1.11.1.7	flavonoid peroxidase 3
PWY-7218 	1.1.1.35 	FadJ component of anaerobic fatty acid oxidation complex
PWY-7218 	1.1.1.35 	fatty acid oxidation complex, &alpha; component
PWY-7218 	1.1.1.36	MONOMER-16780
PWY-7218 	2.3.1.16 	&beta;-ketothiolase
PWY-7218 	2.3.1.9	acetyl-CoA acetyltransferase subunit
PWY-7218 	2.7.2.3	PGK
PWY-7218 	3.1.3.11	fructose-1,6-bisphosphatase
PWY-7218 	3.1.3.11	fructose-1,6-bisphosphatase I
PWY-7218 	3.1.3.11	fructose 1,6-bisphosphatase II
PWY-7218 	3.1.3.74 	EG11239-MONOMER
PWY-7218 	4.1.2.13	fructose bisphosphate aldolase class I
PWY-7218 	5.3.1.1	triosephosphate isomerase
PWY-722	1.13.11.9	2,5-dihydroxypyridine dioxygenase
PWY-722	1.13.11.9	2,5-dihydroxypyridine oxygenase subunit
PWY-722	1.14.13.114	6-hydroxynicotinate 3-monooxygenase
PWY-722	1.17.2.1	nicotinate dehydrogenase (cytochrome)
PWY-722	3.5.1.106	MONOMER-15549
PWY-722	3.5.1.107	maleamate amidohydrolase
PWY-7224 	2.7.1.145 	deoxynucleoside kinase
PWY-7224 	2.7.1.74 	deoxycytidine kinase
PWY-7224	2.7.1.74 	mitochondrial deoxyguanosine kinase
PWY-7224 	2.7.4.13 	adenylate kinase
PWY-7224 	2.7.4.13 	adenylate kinase isoenzyme 5
PWY-7224 	2.7.4.13 	guanylate kinase
PWY-722	5.2.1.1	maleate cis-trans isomerase
PWY-7228 	1.1.1.205	MONOMER-510
PWY-7228 	1.17.4.1	ribonucleoside-diphosphate reductase
PWY-7228 	2.7.4.13 	guanylate kinase
PWY-7228 	6.3.5.2 	GMP synthase
PWY-7229 	4.3.2.2	Ade13
PWY-7230 	1.14.13	3-hexaprenyl-benzoate-5-monoxygenase
PWY-7230 	2.1.1.64 	hexaprenyldihydroxybenzoate methyltransferase
PWY-723	1.13.12.16	nitronate monooxygenase
PWY-7233 	2.5.1 	hexaprenyl-diphosphate hexaprenyltransferase
PWY-7234	3.5.4.10	MONOMER-14617
PWY-7234 	4.3.2.2 	adenylosuccinate lyase
PWY-7234	6.3.4.23	FAICAR synthetase monomer
PWY-7235 	1.14.99	MONOMER3O-164
PWY-7235 	1.18.1.2 	ferredoxin oxidoreductase
PWY-7235 	2.1.1.201	YML110C-MONOMER
PWY-7236	1.14.21.9	mycocyclosin synthase
PWY-7236	2.3.2.21	cyclo(L-tyrosyl-L-tyrosyl) synthase
PWY-7237	1.1.1.370	scyllo-inositol 2-dehydrogenase (NAD+)
PWY-7237 	1.2.1.27 	(methyl)malonate-semialdehyde dehydrogenase
PWY-7237 	2.7.1.92	5-dehydro-2-deoxy-D-gluconate kinase
PWY-7237 	3.7.1.22	3D-(3,5/4)-trihydroxycyclohexane-1,2-dione hydrolase
PWY-7237 	4.1.2.29	2-deoxy-5-keto-D-gluconate 6-phosphate aldolase
PWY-7237 	4.2.1.44	scyllo-inosose dehydratase
PWY-7237 	5.3.1.30	5-deoxy-D-glucuronate isomerase
PWY-7237	5.3.99.11	scyllo-inosose epimerase/isomerase
PWY-7238 	2.4.1.14	sucrose phosphate synthase
PWY-7238	2.4.1.1	&alpha;-glucan phosphorylase
PWY-7238	2.4.1.25	disproportionating enzyme
PWY-7238	2.7.1.1 	hexokinase
PWY-7241	1.1.1.18	MONOMER-17948
PWY-7241	1.1.1	MONOMER-17949
PWY-724 	1.2.1.11	aspartate semialdehyde dehydrogenase
PWY-7241 	5.3.3	D-tagaturonate epimerase
PWY-7241	5.3.3	MONOMER-17951
PWY-724 	2.1.1.14	AT3G03780-MONOMER
PWY-724 	2.1.1.14	methionine synthase
PWY-724 	2.5.1	cystathionine &gamma;-synthase
PWY-724 	2.6.1.83	AT4G33680-MONOMER
PWY-724 	2.7.1.39	MONOMER-1961
PWY-724 	2.7.2.4	AT3G02020-MONOMER
PWY-7243 	3.1.1.11	pectin methylesterase A
PWY-7243	3.2.1.40	MONOMER-18216
PWY-7243 	3.2.1.67	galacturan 1,4-&alpha;-galacturonidase
PWY-7243	4.2.2.22	pectate trisaccharide-lyase
PWY-7243	4.2.2.9	exopolygalacturonate lyase W
PWY-724 	4.1.1.20	AT3G14390-MONOMER
PWY-724 	4.1.1.20	AT5G11880-MONOMER
PWY-7245 	1.1.1.42	NADP-dependent isocitrate dehydrogenase
PWY-7245 	1.1.1.49	glucose-6-phosphate dehydrogenase
PWY-7245 	1.2.1.4 	magnesium-activated aldehyde dehydrogenase, cytosolic
PWY-7245 	1.6.5.9	NADH:ubiquinone oxidoreductase 1 (external)
PWY-7245 	1.6.5.9	NADH:ubiquinone oxidoreductase 2 (external)
PWY-7245 	1.6.5.9	NADH:ubiquinone oxidoreductase (internal)
PWY-7245 	2.7.1.23 	NAD(+)/NADH kinase
PWY-7245 	2.7.1.86 	NAD(+)/NADH kinase
PWY-7246	3.2.1.82	exopolygalacturonan hydrolase X
PWY-7246 	4.2.2 	oligogalacturonate lyase
PWY-7248 	5.3.1.12	hexuronate isomerase
PWY-7250 	2.8.1.7	cysteine desulfurase
PWY-7250 	2.8.1.7	cysteine desulfurase, mitochondrial
PWY-7250	2.8.1.7	L-cysteine desulfurase
PWY-7250	3.6.4	chaperone for [Fe-S] cluster biosynthesis
PWY-7251 	4.2.1.128 	lupeol/lupan-3&beta;,20-diol synthase
PWY-7251	5.4.99.35 	taraxerol synthase
PWY-7251	5.4.99.36	isomultiflorenol synthase
PWY-7251	5.4.99.39	&beta;-amyrin synthase
PWY-7251 	5.4.99.39	MONOMER-14474
PWY-7251 	5.4.99.39	MONOMER-17474
PWY-7251 	5.4.99.41	MONOMER-14455
PWY-7251	5.4.99.49 	glutinol synthase
PWY-7251	5.4.99.50 	friedelin synthase
PWY-7251	5.4.99.55 	&delta;-amyrin synthase
PWY-7252 	1.1.1.234 	dihydroflavonol-4-reductase/flavanone 4-reductase
PWY-7252 	1.1.1.234	flavanone 4-reductase
PWY-7252 	2.4.1	3-deoxyanthocyanidin 5-O-glucosyltransferase
PWY-7253 	1.1.1.234	dihydroflavonol 4-reductase/flavanone 4-reductase
PWY-7254	2.3.3.16 	citrate synthase
PWY-7254	2.8.3.18	succinyl-CoA:acetate CoA-transferase
PWY-7255	1.14.99.50	hercynine oxygenase
PWY-7255	2.1.1.44	<small>L</small>-histidine N<small><sup>&alpha;</sup></small>-methyltransferase
PWY-7255	3.5.1.118	MONOMER-17986
PWY-7255	6.3.2.2	MONOMER-17988
PWY-7260 	2.4.1.299	anthocyanidin 3-O-glucoside 5-O-glucosyltransferase
PWY-7260 	2.4.1.299 	anthocyanidin 3-O-glucoside 5-O-glucosyltransferase (acyl-glucose dependent)
PWY-7261	2.3.1	anthocyanidin 3-O-glucoside-6\-O-malyltransferase"
PWY-7265	2.3.1	1-O-hydroxycinnamoyl-&beta;-glucose:betanin O-hydroxycinnamoyltransferase
PWY-7270 	2.5.1.6	methionine adenosyltransferase
PWY-7270 	2.6.1	2-oxo-4-methylthiobutanoate-glutamine aminotransferase
PWY-7270 	2.6.1.88 	AT3G19710-MONOMER
PWY-7270 	3.1.3.87 	dehydratase-enolase-phosphatase
PWY-7270 	3.5.1.111	2-oxoglutaramate amidase
PWY-7270 	5.3.1.23	5-methylthioribose-1-phosphate isomerase
PWY-7274	1.14.11	MONOMER-18020
PWY-7274	2.3.1.30	MONOMER-18015
PWY-7274	2.5.1.47 	O-ureido-L-serine synthase
PWY-7274	3.5.3.25	N<sup>&omega;</sup>-hydroxy-L-arginine amidinohydrolase
PWY-7274	5.1.1.19	O-ureido-serine racemase
PWY-7274	6.3.3.5	O-ureido-D-serine cyclase
PWY-7275	2.3.3	2-benzylmalate synthase
PWY-7275	4.2.1	benzylmalate isomerase
PWY-7277 	1.1.1.102	3-ketodihydrosphingosine reductase
PWY-7277 	1.14.19.17	sphingolipid 4-desaturase 1
PWY-7277 	2.3.1.50	serine palmitoyltransferase IA
PWY-7277 	2.3.1.50	serine palmitoyltransferase IIA
PWY-7277 	2.7.8.27	sphingomyelin synthase 1
PWY-7277 	2.7.8.27	sphingomyelin synthase 2
PWY-7277 	3.1.4.12	ectonucleotide pyrophosphatase/phosphodiesterase 7
PWY-7277 	3.1.4.12	sphingomyelin phosphodiesterase 1, acid lysosomal
PWY-7277 	3.1.4.12	sphingomyelin phosphodiesterase 2, neutral
PWY-7277 	3.1.4.12	sphingomyelin phosphodiesterase 3, neutral
PWY-7277 	3.1.4.12	sphingomyelin phosphodiesterase 4
PWY-7277 	3.5.1.109	sphingomyelin deacylase
PWY-7279	1.10.2.2	ubiquinol-cytochrome C oxidoreductase
PWY-7279	1.9.3.1	cytochrome c oxidase
PWY-7280	2.3.1	SPLCAT1 subunit 1 
PWY-7282	1.1.1.65	pyridoxine 4-dehydrogenase
PWY-7282	1.4.3.5	pyridoxamine 5-phosphate oxidase/pyridoxine 5-phosphate oxidase
PWY-7282	2.7.1.35	pyridoxine kinase
PWY-7283 	2.1.1.228	tRNA (guanine<small><sup>37</sup></small>-N<small><sup>1</sup></small>)-methyltransferase
PWY-7283	2.1.1.282	tRNA wybutosine-synthesizing protein 3
PWY-7283	2.3.1.231 	tRNA wybutosine-synthesizing protein 4
PWY-7283	2.5.1.114	tRNA wybutosine-synthesizing protein 2
PWY-7283	4.1.3.44	tRNA-4-demethylwyosine synthase
PWY-7285 	2.1.1.228	tRNA (guanine<small><sup>37</sup></small>-N<small><sup>1</sup></small>)-methyltransferase
PWY-7285 	2.1.1 	MONOMER-18045
PWY-7285 	4.1.3.44	tRNA-4-demethylwyosine synthase
PWY-7286 	2.1.1.282 	tRNA isowyosine<sup>37</sup> 7-methyltransferase 1
PWY-7286 	2.1.1.282 	tRNA isowyosine<sup>37</sup> 7-methyltransferase 2
PWY-7286	2.5.1.114	tRNA 4-demethylwyosine &alpha;-amino-&alpha;-carboxypropyltransferase
PWY-7287 	1.1.1	&beta;-hydroxy-L-tyrosine-[NovH] dehydrogenase
PWY-7287 	1.1.1	dTDP-4-dehydro-5-methyl-L-rhamnose 4-ketoreductase
PWY-7287 	1.13.12 	3-dimethylallyl-4-hydroxyphenylpyruvate oxygenase
PWY-7287 	1.14.13	MONOMER-18085
PWY-7287 	1.3.1.12	prephenate dehydrogenase
PWY-7287	2.1.1.284	8-demethylnovobiocic acid C8-methyltransferase
PWY-7287	2.1.1.285	demethyldecarbamoyl novobiocin O-methyltransferase
PWY-7287 	2.1.1	dTDP-4-dehydro-&beta;-L-rhamnose C5-methyltransferase
PWY-7287	2.1.3.12	decarbamoylnovobiocin carbamoyltranferase
PWY-7287	2.4.1.302	4-O-demethyl-L-noviosyl transferase
PWY-7287 	2.5.1.111	4-hydroxyphenylpyruvate dimethylallyl transferase
PWY-7287 	2.7.7.24	glucose-1-phosphate thymidylyltransferase
PWY-7287 	4.2.1.46	dTDP-glucose 4,6-dehydratase
PWY-7287 	5.1.3.13	dTDP-4-dehydro-6-deoxy-&alpha;-D-glucose 3,5-epimerase
PWY-7287	6.3.1.15	MONOMER-18056
PWY-7289	2.5.1.113	CysO-thiocarboxylate-dependent cysteine synthase
PWY-7289	2.7.7 	[sulfur carrier protein CysO] adenylyltransferase/sulfurtransferase
PWY-7289	3.4.17	G185E-5513-MONOMER
PWY-7290	2.4.1.306	MONOMER-18063
PWY-7290	2.4.1.307	MONOMER-18062
PWY-7290	2.4.1.308	MONOMER-18061
PWY-7290	2.4.1.309	MONOMER-18060
PWY-7290	2.7.8.33	UDP-N-acetylglucosamine&mdash;undecaprenyl-phosphate N-acetylglucosaminephosphotransferase
PWY-7290	5.1.3.26	MONOMER-18098
PWY-7292	3.1.2 	acyl-CoA thioesterase
PWY-7294 	1.1.1.26	MONOMER-18075
PWY-7294 	1.2.1	glycolaldehyde oxidoreductase
PWY-7294 	2.3.3.9	malate synthase subunit
PWY-7294 	4.1.2.28 	KDG-aldolase subunit
PWY-7295 	1.1.1.376 	aldose 1-dehydrogenase
PWY-7295 	4.2.1.25 	D-xylonate dehydratase
PWY-7297	1.14.17	MONOMER-18080
PWY-7297	1.14.17	MONOMER-18081
PWY-7297	4.1.1.25	tyrosine decarboxylase 2
PWY-7298	1.14.11	MONOMER-18078
PWY-7298	2.1.1	flavonoid 8-O-methyltransferase
PWY-7298	2.1.1 	phenylpropanoid and flavonoid O-methyltransferase
PWY-7299 	5.3.3.1 	3 &beta;-hydroxysteroid dehydrogenase/ &Delta; 5-->4-isomerase type 2
PWY-7300	1.1.1	MONOMER-18102
PWY-7300	1.14.14	cytochrome P450 306a1
PWY-7300	1.14.14	cytochrome P450 306A1
PWY-7300	1.14.15	cytochrome P450 302a1, mitochondrial
PWY-7300	1.14.15	cytochrome P450 315a1, mitochondrial
PWY-7300	1.14.15	cytochrome P450 CYP302A1
PWY-7300	1.14.15	cytochrome P450 CYP315A1
PWY-7300	1.14.19.21	MONOMER-18094
PWY-7300	1.14.19.21	MONOMER-18095
PWY-7300	1.14.99.22	cytochrome P450 314a1
PWY-7300	1.14.99.22	cytochrome P450 314a1, mitochondrial
PWY-7300	1.1.99	3-oxoecdysteroid 3&beta;-reductase
PWY-7300	1.3.1	MONOMER-18100
PWY-7305 	1.1.1.62 	3-keto-steroid reductase
PWY-7305 	1.1.1.64	testosterone 17-beta-dehydrogenase 3
PWY-7305 	1.14.15.5 	cyp11B2
PWY-7305 	1.14.15.6	cholesterol side chain cleavage enzyme, mitochondrial
PWY-7305 	5.3.3.1 	3 &beta;-hydroxysteroid dehydrogenase/ &Delta; 5-->4-isomerase type 1
PWY-7307 	5.3.3	&Delta;3,5-&Delta;2,4-dienoyl-CoA isomerase
PWY-7308	4.2.1.84	nitrile hydratase
PWY-7310	4.1.2.55 	2-dehydro-3-deoxy-phosphogluconate aldolase
PWY-7310	4.3.1.29	D-glucosaminate 6-phosphate ammonia-lyase
PWY-7317 	1.1.1.133	DTDPDEHYRHAMREDUCT-MONOMER
PWY-7317 	1.1.1.133	MONOMER-17551
PWY-7317 	1.1.1.266	dTDP-4-dehydro-6-deoxyglucose reductase [NAD(P)H]
PWY-7317	1.1.1.266	MONOMER-18134
PWY-7317	1.1.1.339	dTDP-6-deoxy-<small>L</small>-talose 4-dehydrogenase (NAD<small><sup>+</sup></small>)
PWY-7317 	2.3.1.197	dTDP-3-amino-3,6-dideoxy-&alpha;-D-galactopyranose 3-N-acetyltransferase
PWY-7317 	2.3.1.209	dTDP-4-amino-4,6-dideoxy-<small>D</small>-glucose acetyltransferase
PWY-7317 	2.3.1.210	TDPFUCACTRANS-MONOMER
PWY-7317 	2.3.1	dTDP-3-amino-3,6-dideoxy-&alpha;-D-glucopyranose:acetyl-CoA acetyltransferase
PWY-7317 	2.6.1.33	dTDP-4-amino-4,6-dideoxy-<small>D</small>-glucose transaminase
PWY-7317 	2.6.1.59	dTDP-4-dehydro-6-deoxy-D-glucose transaminase
PWY-7317 	2.6.1.89	dTDP-3-amino-3,6-dideoxy-&alpha;-<small>D</small>-glucopyranose transaminase
PWY-7317 	2.6.1.90	dTDP-3-amino-3,6-dideoxy-&alpha;-D-galactopyranose transaminase
PWY-7317 	2.7.7.24	dTDP-glucose pyrophosphorylase
PWY-7317 	2.7.7.24	dTDP-glucose pyrophosphorylase 2
PWY-7317 	2.7.7.24	glucose-1-phosphate thymidylyltransferase
PWY-7317 	2.7.7.24	MONOMER-17548
PWY-7317 	4.2.1.46	dTDP-glucose 4,6-dehydratase
PWY-7317 	4.2.1.46	dTDP-glucose 4,6-dehydratase 2
PWY-7317 	4.2.1.46	MONOMER-17549
PWY-7317 	5.1.3.13	DTDPDEHYDRHAMEPIM-MONOMER
PWY-7317 	5.1.3.13	MONOMER-17550
PWY-7317 	5.3.2.3 	bifunctional ketoisomerase / N-acetyltransferase FdtD
PWY-7317 	5.3.2.3	dTDP-6-deoxy-3,4-keto-hexulose isomerase
PWY-7317 	5.3.2.4	TDP-4-oxo-6-deoxy-&alpha;-<small>D</small>-glucose-3,4-oxoisomerase
PWY-7317 	5.4.99.59	dTDP-fucopyranose mutase
PWY-7321	1.1.1	MONOMER-18158
PWY-7321	1.1.1	MONOMER-18159
PWY-7321	1.1.3.16	MONOMER-16688
PWY-7321	1.1.3.16	MONOMER-18162
PWY-7321	1.1.3.16	MONOMER-18163
PWY-7321	1.14.14	cytochrome P450 18a1
PWY-7321	2.7.1	ecdysteroid 22-kinase
PWY-7321	3.1.3	ecdysteroid-phosphate phosphatase
PWY-7323 	1.1.1.187 	MONOMER-13570
PWY-7323 	1.1.1.271	GDP-fucose synthase
PWY-7323 	1.1.1.281 	MONOMER-12852
PWY-7323 	1.1.1.356	GDP-L-collitose synthase subunit
PWY-7323 	2.6.1.102	GDP-4-dehydro-6-deoxy-D-mannose-4-aminotransferase subunit
PWY-7323 	2.6.1.54	GDP-4-keto-6-deoxy-D-mannose-3-dehydrase subunit
PWY-7323 	2.7.7.13	MANNPGUANYLTRANGDP-MONOMER
PWY-7323 	4.2.1.47	GDP-D-mannose 4,6-dehydratase subunit
PWY-7323 	4.2.1.47	GDP-mannose 4,6-dehydratase
PWY-7323 	4.2.1.47	MONOMER-13569
PWY-7323 	4.2.1.47	MONOMER-13574
PWY-7323 	5.3.1.8	MANNPISOM-MONOMER
PWY-7323 	5.3.1.8 	MONOMER-13384
PWY-7323 	5.4.2.2 	phosphomannomutase
PWY-7323 	5.4.2.8	PHOSMANMUT-MONOMER
PWY-7325	1.14.13	apigenin synthase
PWY-7325	1.14.13	genkwanin 6-hydroxylase
PWY-7325	2.1.1	flavonoid 4-O-methyltransferase
PWY-7325	2.1.1	flavonoid 6-O-methyltransferase
PWY-7325	2.1.1	flavonoid 7-O-methyltransferase
PWY-7328 	1.1.1.22	UDP-glucose 6-dehydrogenase
PWY-7328 	2.7.7.64 	UTP:glucose-1-phosphate uridylyltransferase
PWY-7328 	5.1.3.2	MONOMER-6142
PWY-7328 	5.1.3.2	UDP-glucose 4-epimerase
PWY-7328 	5.1.3.6	UDP-glucuronate 4-epimerase
PWY-7328 	5.4.2.2	PHOSPHOGLUCMUT-MONOMER
PWY-7332 	1.1.1.136	UDP-N-acetyl-&alpha;-D-glucosamine 6-dehydrogenase
PWY-7332 	1.1.1.335	UDP-N-acetyl-2-amino-2-deoxyglucuronate dehydrogenase
PWY-7332 	1.1.1	UDP-N-acetyl-&alpha;-D-fucosamine dehydrogenase
PWY-7332 	1.1.1	UDP-N-acetyl-&alpha;-D-quinovosamine dehydrogenase
PWY-7332 	1.1.1	UDP-N-acetyl-&beta;-L-rhamnosamine dehydrogenase
PWY-7332 	1.1.1 	UDP-N-acetylgalactosamine 6-dehydrogenase
PWY-7332 	2.3.1.201	UDP-2-acetamido-3-amino-2,3-dideoxy-glucuronate N-acetyltransferase
PWY-7332 	2.6.1.16	L-glutamine:D-fructose-6-phosphate aminotransferase
PWY-7332 	2.6.1.98	UDP-2-acetamido-2-deoxy-ribo-hexuluronate aminotransferase
PWY-7332 	2.7.7.23 	fused N-acetylglucosamine-1-phosphate uridyltransferase and glucosamine-1-phosphate acetyltransferase
PWY-7332	3.5.3	UDP-2,3-diacetamido-2,3-dideoxy-&alpha;-D-mannuronate amidotranferase
PWY-7332	3.5.3	UDP-N-acetyl-&beta;-L-fucosamine amidotransferase
PWY-7332 	4.2.1.115	UDP-N-acetylglucosamine 4,6-dehydratase (configuration-inverting)
PWY-7332 	4.2.1.115	UDP-N-acetylglucosamine 4,6-dehydratase (inverting)
PWY-7332 	4.2.1.135	UDP-N-acetylglucosamine 4,6-dehydratase (configuration-retaining)
PWY-7332 	5.1.3.23	UDP-2,3-diacetamido-2,3-dideoxyglucuronic acid 2-epimerase
PWY-7332 	5.1.3.28	UDP-2-acetamido-2,6-dideoxy-&beta;-L-talose 2-epimerase
PWY-7332 	5.1.3	UDP-N-acetyl-&beta;-L-rhamnosamine 2-epimerase
PWY-7332 	5.4.2.10	phosphoglucosamine mutase
PWY-7334	1.1.1	UDP-N-acetyl-&alpha;-D-quinovosamine dehydrogenase
PWY-7335 	5.1.3.14	UDP-N-acetyl glucosamine 2-epimerase
PWY-7336 	5.1.3 	UDP-N-acetyl-&alpha;-D-glucosaminouronate C4-epimerase
PWY-7338 	1.3.1.34	2,4-dienoyl-CoA reductase
PWY-7340 	1.1.1.M19 	peroxisomal multifunctional enzyme type 2
PWY-7340 	1.3.3.6	fatty-acyl coenzyme A oxidase
PWY-7340 	2.3.1.16 	3-oxoacyl CoA thiolase
PWY-7340 	5.3.3 	&Delta;3-&Delta;2 -enoyl-CoA isomerase
PWY-7341 	1.1.1.206	MONOMER-13848
PWY-7341 	1.1.1.206	tropinone reductase (TR-I)
PWY-7341 	1.14.11.11	MONOMER-12447
PWY-7341 	1.14.11.14	hyoscyamine 6&beta;-hydroxylase
PWY-7341 	1.14.11.14	MONOMER-16825
PWY-7341 	1.4.3.22	diamine oxidase
PWY-7341 	2.1.1.53	MONOMER-12434
PWY-7341 	2.3.1.185	MONOMER-16824
PWY-7342 	1.4.3.22	N-methylputrescine oxidase
PWY-7342 	2.1.1.53	MONOMER-12436
PWY-7342 	2.1.1.53	MONOMER-12437
PWY-7342 	2.4.2.19	MONOMER-12622
PWY-7342 	3.2.2	MONOMER-12670
PWY-7343	2.7.7.64 	UDP-glucose pyrophosphorylase subunit
PWY-7343	2.7.7.64 	UTP--glucose-1-phosphate uridylyltransferase
PWY-7343	2.7.7.64 	UTP-glucose-1-phosphate uridylyltransferase
PWY-7344 	5.1.3.2 	UDP-galactose 4-epimerase
PWY-7345 	1.1.1.1	alcohol dehydrogenase
PWY-7345 	1.1.1.1	MONOMER-15098
PWY-7345 	1.1.1.1	MONOMER-15099
PWY-7345 	1.1.1.1	MONOMER-15105
PWY-7345 	1.1.1.27	MONOMER-9224
PWY-7345 	2.4.1.13	sucrose synthase
PWY-7345 	2.7.1.11	chloroplastic 6-phosphofructokinase
PWY-7345 	2.7.1.40	plastidic pyruvate kinase
PWY-7345 	2.7.1.4 	AT4G29130-MONOMER
PWY-7345	2.7.1.4	fructokinase I
PWY-7345	2.7.1.4	fructokinase-Ia
PWY-7345	2.7.1.4	fructokinase-Ib
PWY-7345	2.7.1.4	fructokinase II
PWY-7345	2.7.1.4	fructokinase-II
PWY-7345	2.7.1.4	fructokinase (recombinant)
PWY-7345 	2.7.2.3	plastidic 3-phosphoglycerate kinase
PWY-7345 	2.7.7.64 	UDP-glucose pyrophosphorylase
PWY-7345 	3.1.3.11	fructose-1,6-bisphosphatase
PWY-7345	4.1.1.1	MONOMER-15104
PWY-7345 	4.1.1.1	pyruvate decarboxylase
PWY-7345 	4.1.2.13	chloroplastic fructose-bisphosphate aldolase
PWY-7345 	4.2.1.11	plastidic enolase
PWY-7345 	5.3.1.9	cytosolic glucose-6-phosphate isomerase
PWY-7345 	5.3.1.9	plastidic phosphoglucose isomerase
PWY-7345 	5.4.2.12	plastidic cofactor-independent phosphoglycerate mutase
PWY-7345 	5.4.2.2	StpPGM
PWY-7346	1.1.1.22	UDP-glucose 6-dehydrogenase
PWY-7347	3.1.3.24 	sucrose-phosphate synthase/phosphatase
PWY-7351	1.5.1.11	octopine dehydrogenase
PWY-7351	1.5.1.17	alanopine dehydrogenase
PWY-7351	1.5.1.17 	strombine/alanopine dehydrogenase
PWY-7351	1.5.1.23	tauropine dehydrogenase
PWY-7351	1.5.1.26	&beta;-alanopine dehydrogenase
PWY-7352	1.1.1.362	aklaviketone reductase
PWY-7352	1.1.1	daunorubicin polyketide synthase reductase
PWY-7352	1.13.12	aklanoate anthrone monooxygenase
PWY-7352	1.14.13.180	aklavinone 12-hydroxylase
PWY-7352	1.14.13.181	MONOMER-18172
PWY-7352	1.14.13 	daunorubicin polyketide synthase aromatase subunit
PWY-7352	2.1.1.288	aklanonic acid methyltransferase
PWY-7352	2.1.1.292	carminomycin 4-O-methyltransferase
PWY-7352	2.3.1	daunorubicin polyketide synthase complex
PWY-7352	2.4.1	dTDP-daunosamine transferase
PWY-7352	3.1.1	rhodomycinone D methylesterase
PWY-7352	3.1.1	rhodomycinone D methyl-esterase
PWY-7352	4.2.1	12-deoxyaklanoate synthase
PWY-7352	5.5.1.23	aklanonic acid methyl ester cyclase
PWY-7354	1.1.1.362	aklaviketone reductase
PWY-7354	1.1.1	aclacinomycin polyketide synthase reductase
PWY-7354	1.13.12	aklanoate anthrone monooxygenase
PWY-7354	1.3.3.14 	aclacinomycin oxidase
PWY-7354	2.1.1.288	aklanonic acid methyltransferase
PWY-7354	2.3.1	aclacinomycin polyketide synthase complex
PWY-7354	2.4.1.326	L-rhodosaminyltransferase
PWY-7354	2.4.1.327	MONOMER-18193
PWY-7354	4.2.1	aclacinomycin polyketide synthase aromatase
PWY-7354	5.5.1.23	aklanonic acid methyl ester cyclase
PWY-7355	1.14.13.181	MONOMER-18171
PWY-7363	2.1.1.121	MONOMER-18209
PWY-7363	2.1.1.291	(R,S)-reticuline 7-O-methyltransferase
PWY-7363	2.1.1	MONOMER-18210
PWY-7366	3.1.4 	phosphatidylglycerol phospholipase C
PWY-7367	3.1.1 	Neuropathy Target Esterase
PWY-7369	2.7.4.15	thiamine diphosphate kinase
PWY-7371	1.21.98.1	dehypoxanthinylfutalosine cyclase
PWY-7371	3.2.2.30	aminodeoxyfutalosine nucleosidase
PWY-7373 	1.21.98.1	dehypoxanthinylfutalosine cyclase
PWY-7373 	2.5.1.83	hexaprenyl-diphosphate synthase
PWY-7373 	3.2.2.16 	methylthioadenosine nucleosidase
PWY-7373 	4.2.1.151	chorismate dehydratase
PWY-7374	1.21.98.1	dehypoxanthinylfutalosine cyclase
PWY-7374	2.5.1.120	aminodeoxyfutalosine synthase
PWY-7374	3.2.2.26	futalosine hydrolase
PWY-7374	3.5.4.40	6-aminodeoxyfutalosine deaminase
PWY-7374	4.2.1.151	chorismate dehydratase
PWY-7375	2.1.1.56	mRNA (guanine-N<small><sup>7</sup></small>)-methyltransferase
PWY-7375	2.7.7.50	mRNA capping enzyme &alpha; subunit
PWY-7375	3.1.3.33	mRNA capping enzyme &beta; subunit
PWY-7379 	2.1.1.56	mRNA cap guanine-N7 methyltransferase
PWY-7379	2.1.1.57	cap1 methyltransferase
PWY-7379 	3.1.3.33 	mRNA capping enzyme
PWY-7382	2.3.1.181	octanoyl transferase
PWY-7382	2.3.1	octanoyl-CoA-protein transferase
PWY-7382	2.8.1.8	lipoyl synthase
PWY-7385 	1.1.1.94 	glycerol-3-phosphate dehydrogenase, biosynthetic
PWY-7385	3.1.3.21	DL-glycerol-3-phosphatase
PWY-7385	3.1.3.21 	sugar phosphatase
PWY-7387	1.2.1.3	aldehyde dehydrogenase (NAD<small><sup>+</sup></small>)
PWY-7387	2.6.1.77	hypotaurine&mdash;pyruvate aminotransferase
PWY-7388	2.3.1.86 	malonyl-CoA:ACP transferase
PWY-7388	6.4.1.2	Acetyl-CoA carboxylase, mitochondrial
PWY-7389 	1.1.1.37	MONOMER-18271
PWY-7389 	1.1.1.38 	NAD<sup>+</sup>-dependent malic enzyme, mitochondrial
PWY-7389 	1.3.5.4	fumarate reductase complex
PWY-7389 	1.8.1.4 	pyruvate dehydrogenase complex
PWY-7389 	2.6.1.1	MONOMER-18269
PWY-7389 	2.6.1.2	MONOMER-18253
PWY-7389 	2.7.1.40	MONOMER-18252
PWY-7389 	2.8.3	MONOMER-18292
PWY-7389 	2.8.3 	succinyl-CoA:acetate CoA-transferase
PWY-7389 	4.1.1.32	phosphoenolpyruvate carboxykinase
PWY-7389 	4.2.1.2	fumarase subunit
PWY-7389 	5.1.1.1	alanine racemase subunit
PWY-7389 	5.1.99.1	MONOMER-18294
PWY-7389 	5.4.99.2	methylmalonyl-CoA mutase subunit
PWY-7389 	6.2.1.5	MONOMER-18307
PWY-7389 	6.4.1.3	MONOMER-18296
PWY-7391	1.1.1.34 	acetyl-CoA acetyltransferase/HMG-CoA reductase
PWY-7391	2.3.3.10	hydroxymethylglutaryl-CoA synthase
PWY-7392 	2.5.1.29	geranylgeranyl diphosphate synthase
PWY-7393 	1.3.99.31	phytoene desaturase
PWY-7393	5.5.1.19	lycopene &beta;-cyclase
PWY-7394	1.14.13.113	FAD-dependent urate hydroxylase
PWY-7394	1.14.13.113	MONOMER-15359
PWY-7394	3.5.2.17	hydroxyisourate hydrolase
PWY-7394	4.1.1.97	OHCU decarboxylase
PWY-7395	2.7.1.11 	6-phosphofructokinase
PWY-7395	3.5.1.25	N-acetylglucosamine-6-phosphate deacetylase
PWY-7395	3.5.99	D-galactosamine-6-phosphate deaminase/isomerase
PWY-7395	4.1.2.40	D-tagatose-1,6-bisphosphate aldolase
PWY-7396 	2.3.3.7 	malate synthase 1
PWY-7396	2.3.3.7 	malate synthase 2
PWY-7397	2.3.1.74	naringenin-chalcone synthase
PWY-7397 	4.3.1.23 	phenylalanine/tyrosine ammonia-lyase
PWY-7397 	5.5.1.6	chalcone isomerase 1
PWY-7397 	6.2.1.12	4-coumarate--CoA ligase 1
PWY-7398 	1.14.14 	4-coumarate 3-monooxygenase
PWY-7399	3.1.4.57	cyclic phosphate dihydrolase
PWY-7400	2.1.3.3	L-ornithine carbamoyltransferase
PWY-7400	6.3.2	glutamate--LysW ligase
PWY-7402 	1.3.7.8 	cyclohex-1-ene-1-carbonyl-CoA dehydrogenase
PWY-7402 	1.3.8.11	cyclohexane-1-carbonyl-CoA dehydrogenase
PWY-7402 	3.7.1.21	6-oxocyclohex-1-ene-1-carbonyl-CoA hydratase
PWY-7402 	4.2.1.100	cyclohexa-1,5-dienecarbonyl-CoA hydratase
PWY-7403	1.1.1	aminoalcohol dehydrogenase
PWY-7403	1.14.12	tetramethylpyrazine oxygenase system
PWY-7403	3.5.1	(Z)-N,N-(but-2-ene-2,3-diyl)diacetamide hydrolase
PWY-7404	2.7.8	MONOMER-18328
PWY-7405	1.14.13	MONOMER-18334
PWY-7405	2.3.1	polyketide synthase complex type II &beta;-ketoacyl-[acp] synthase subunit 
PWY-7405	2.5.1	2-methyl-4-hydroxy-quinoline farnesyltransferase
PWY-7405	6.2.1.32	MONOMER-18329
PWY-7407	1.14.13.da	aurachin C monooxygenase/isomerase
PWY-7407	1.14.13	MONOMER-18337
PWY-7407	1.14.13	MONOMER-18347
PWY-7407	2.3.1	4-hydroxy-2-methylquinoline synthase complex
PWY-7407	2.5.1	2-methyl-4-hydroxyquinoline farnesyltransferase
PWY-7407	6.2.1.32	MONOMER-18342
PWY-7409 	3.1.1.34 	triacylglycerol lipase
PWY-7409 	3.1.1.4 	triacylglycerol lipase
PWY-7409	3.1.1.5	phospholipase B/lyso-phospholipase
PWY-7409 	3.1.4.46 	glycerophosphocholine phosphodiesterase
PWY-7410	1.14.14.ae	myrcene hydroxylase
PWY-7410	4.2.3.15 	geranyl diphosphate synthase / myrcene synthase
PWY-7411	1.1.1.101	1-acyl dihydroxyacetone phosphate reductase
PWY-7411 	1.1.1.8	glycerol 3-phosphate dehydrogenase
PWY-7411 	2.3.1.42 	glycerol-3-phosphate 1-O-acyltransferase
PWY-7411 	2.3.1.42 	glycerol-3-phosphate acyltransferase
PWY-7412	1.14.13	mycinamicin IV hydroxylase/epoxidase
PWY-7412	1.14.13	mycinamicin VIII monooxygenase complex
PWY-7412	2.1.1.237	mycinamicin III 3-O-methyltransferase
PWY-7412	2.1.1.238	mycinamicin VI 2-O-methyltransferase
PWY-7412	2.3.1	protomycinolide IV synthase complex
PWY-7412	2.4.1	mycinamicin VII 6-deoxyallosyltransferase
PWY-7412	2.4.1	protomycinolide IV desosaminyltransferase
PWY-7413	1.1.1.364	dTDP-4-dehydro-6-deoxy-&alpha;-D-gulose 4-ketoreductase
PWY-7413	2.7.7.24	glucose-1-phosphate thymidylyltransferase
PWY-7413	4.2.1.46	dTDP-glucose 4,6-dehydratase
PWY-7413	5.1.3.27	dTDP-4-dehydro-6-deoxy-D-glucose 3-epimerase
PWY-7414	2.1.1.235	dTDP-3-amino-3,6-dideoxy-&alpha;-<small>D</small>-glucopyranose N,N-dimethyltransferase
PWY-7414	2.6.1.89	dTDP-3-amino-3,6-dideoxy-&alpha;-<small>D</small>-glucopyranose transaminase
PWY-7414 	2.7.7.24	glucose-1-phosphate thymidylyltransferase
PWY-7414 	4.2.1.46	dTDP-glucose 4,6-dehydratase
PWY-7414	5.3.2.4	TDP-4-oxo-6-deoxy-&alpha;-<small>D</small>-glucose-3,4-oxoisomerase
PWY-7415	1.1.1 	5-O-mycaminosyltylactone C20-hydroxylase
PWY-7415	1.14.13.186	20-oxo-5-O-mycaminosyltylactone 23-monooxygenase
PWY-7415	2.1.1 	demethylmacrocin O-methyltransferase
PWY-7415	2.1.1 	macrocin O-methyltransferase
PWY-7415	2.3.1	MONOMER-18394
PWY-7415	2.3.1	tylactone synthase complex
PWY-7415	2.4.1.316	tylactone mycaminosyltransferase complex
PWY-7415	2.4.1.317	O-mycaminosyltylonolide 6-deoxyallosyltransferase
PWY-7415	2.4.1.318	demethyllactenocin mycarosyltransferase
PWY-7417 	2.3.1.51	1-acyl-sn-gylcerol-3-phosphate acyl transferase
PWY-7417 	2.3.1.51 	acyl-CoA:lyso-phospholipid acyltransferase
PWY-7417 	2.3.1.51	oleoyl-CoA: lysophosphatidate acyltransferase
PWY-7417 	2.3.1.51 	triacylglycerol lipase
PWY-7419	1.14.11	MONOMER-18404
PWY-7419	1.14.13	MONOMER-18401
PWY-7419	2.3.1	MONOMER-18402
PWY-7419	2.3.3	MONOMER-18399
PWY-7419	3.1.4	MONOMER-18403
PWY-7419	5.4.2.9	MONOMER-18398
PWY-7420 	2.3.1.22 	acyl-CoA:diacylglycerol acyltransferase
PWY-7420 	3.1.1.23	monoacylglycerol lipase
PWY-7421 	1.14.13.185	MONOMER-18405
PWY-7421 	2.3.1.240 	narbonolide/10-deoxymethynolide polyketide synthase complex
PWY-7422 	2.3.1.239 	editing thioesterase
PWY-7422 	2.4.1.277	narbonolide/10-deoxymethynolide desosaminyltransferase complex
PWY-7423	1.14.19.15 	acyl-CoA 11-desaturase/conjugase
PWY-7423	1.2.1.84 	pheromone gland-specific fatty-acyl reductase
PWY-7424	2.3.1.26	acyl-CoA cholesterol acyltransferase
PWY-7424	2.3.1.26	acyl-CoA sterol acyltransferase
PWY-7424 	3.1.1.13 	steryl ester hydrolase
PWY-7424	3.1.1.13	steryl ester hydrolase 1
PWY-7424	3.1.1.13	steryl ester hydrolase 2
PWY-7425	1.3.1.103	2-haloacrylate reductase
PWY-7425	3.8.1.10 	S-2-haloacid dehalogenase
PWY-7426	2.4.1.101	Alpha-1,3-mannosyl-glycoprotein 2-beta-N-acetylglucosaminyltransferase
PWY-7426	2.4.1.143 	Alpha-1,3-mannosyl-glycoprotein 4-beta-N-acetylglucosaminyltransferase A soluble form
PWY-7426	2.4.1.143	&alpha;-1,6-mannosyl-glycoprotein 2-beta-N-acetylglucosaminyltransferase
PWY-7426	2.4.1.144	Beta-1,4-mannosyl-glycoprotein 4-beta-N-acetylglucosaminyltransferase
PWY-7426	2.4.1.145	Alpha-1,3-mannosyl-glycoprotein 4-beta-N-acetylglucosaminyltransferase B
PWY-7426	2.4.1.145	Alpha-1,3-mannosyl-glycoprotein 4-beta-N-acetylglucosaminyltransferase C
PWY-7426	2.4.1.155	Alpha-1,6-mannosylglycoprotein 6-beta-N-acetylglucosaminyltransferase A
PWY-7426	2.4.1.155	Alpha-1,6-mannosylglycoprotein 6-beta-N-acetylglucosaminyltransferase B
PWY-7426	2.4.1.155	alpha-1,6-mannosylglycoprotein 6-beta-N-acetylglucosaminyltransferase V
PWY-7426	3.2.1.114	Alpha-mannosidase 2x
PWY-7426	3.2.1.170 	Alpha-mannosidase 2
PWY-7428	4.2.1	2-haloacrylate hydratase
PWY-7429	1.20.9.1	arsenite oxidoreductase
PWY-7430	1.13.11.17	indole 2,3-dioxygenase
PWY-7431	1.2.1.53	MONOMER-18433
PWY-7431	1.4.3.4	tyramine oxidase
PWY-7431	1.4.99	MONOMER-18430
PWY-7431	4.2.1.87	MONOMER-18429
PWY-7431	4.2.1.88	synephrine dehydratase subunit
PWY-7432	2.6.1.57 	aromatic amino acid transaminase
PWY-7433	2.4.1.102	&beta;-1,3-galactosyl-O-glycosyl-glycoprotein &beta;-1,6-N-acetylglucosaminyltransferase
PWY-7433	2.4.1.122	glycoprotein-N-acetylgalactosamine 3-&beta;-galactosyltransferase 1
PWY-7433	2.4.1.146	&beta;-1,3-galactosyl-O-glycosyl-[glycoprotein] &beta;-1,3-N-acetylglucosaminyltransferase 3
PWY-7433 	2.4.1.41	polypeptide N-acetylgalactosaminyltransferase 1
PWY-7433	2.4.99.4	CMP-N-acetylneuraminate-&beta;-galactosamide-&alpha;-2,3-sialyltransferase 1
PWY-743	3.5.5.8	thiocyanate hydrolase &alpha; subunit 
PWY-7434	2.4.1.149	N-acetyllactosaminide &beta;-1,3-N-acetylglucosaminyltransferase
PWY-7434	2.4.1.150	N-acetyllactosaminide &beta;-1,6-N-acetylglucosaminyl-transferase
PWY-7434	2.4.1.152	&alpha;-(1,3)-fucosyltransferase 9
PWY-7434	2.4.1.38 	&beta;-1,4-galactosyltransferase 1
PWY-7434	2.4.99.1	&beta;-galactoside &alpha;-2,6-sialyltransferase 1
PWY-7434	2.4.99 	CMP-N-acetylneuraminate-&beta;-1,4-galactoside &alpha;-2,3-sialyltransferase
PWY-7435	2.4.1.147	mucin core 3 beta 3-GlcNAc-transferase
PWY-7435	2.4.1.148	mucin core 4 beta 6-GlcNAc-transferase
PWY-7436 	2.1.1.295	2-methyl-6-phytyl-1,4-hydroquinone methyltransferase
PWY-7436 	2.1.1.95	MONOMER-3282
PWY-7436	2.5.1.116	homogentisate geranylgeranyl transferase
PWY-7436 	5.5.1.24	tocopherol cyclase
PWY-7437	2.4.1.255	UDP-N-acetylglucosamine--peptide N-acetylglucosaminyltransferase
PWY-7437	3.2.1.169 	bifunctional protein O-GlcNAcase/histone acetyltransferase
PWY-7438 	1.1.1.384	dTDP-3,4-diketo-2,6-dideoxy-D-glucose 3-ketoreductase
PWY-7438 	1.1.1	D-oliose 4-ketoreductase
PWY-7438 	1.1.1	dTDP-4,6-dihydroxy-2-methyloxan-3-one 4-ketoreductase
PWY-7438 	1.1.1	dTDP-4-dehydro-3-methyl-2,6-dideoxy-&alpha;-L-galactose 4-ketoreductase
PWY-7438 	1.1.1	dTDP-D-mycarose 4-ketoreductase
PWY-7438 	1.17.1	dTDP-(2R,6S)-6-hydroxy-2-methyl-3-oxo-3,6-dihydro-2H-pyran-4-olate 3-ketoreductase (dTDP-4-dehydro-2,6-dideoxy-&alpha;-D-allose-forming)
PWY-7438 	1.17.1	MegBII aldo-keo reductase
PWY-7438 	2.1.1.324	dTDP-4-amino-2,3,4,6-tetradeoxy-D-glucose N,N-dimethyltransferase
PWY-7438 	2.1.1 	D-olivose 4-ketoreductase
PWY-7438 	2.1.1 	dTDP-3-amino-4-oxo-2,3,6-trideoxy-&alpha;-D-glucose N,N-dimethyltransferase
PWY-7438 	2.1.1	dTDP-4-dehydro-2,6-dideoxy-&alpha;-D-allose 3-C-methyltransferase
PWY-7438 	2.6.1.110	dTDP-4-dehydro-2,3,6-trideoxy-D-glucose 4-aminotransferase
PWY-7438 	2.7.7.24	glucose-1-phosphate thymidylyltransferase
PWY-7438 	4.2.1.164	dTDP-4-dehydro-2,6-dideoxy-D-glucose 3-dehydratase
PWY-7438 	4.2.1.46	dTDP-glucose 4,6-dehydratase
PWY-7438 	5.1.3	dTDP-4-dehydro-3-C-methyl-2,6-dideoxy-&alpha;-L-allose 5-epimerase
PWY-7438 	5.1.3	dTDP-4-oxo-2,6-dideoxy-D-glucose 3,5-epimerase
PWY-7440	1.1.1	dTDP-&beta;-L-4-epi-vancosamine 4-reductase
PWY-7440	2.1.1	dTDP-3-amino-4-dehydro-2,3,6-trideoxy-&alpha;-D-glucose C3-methyltransferase
PWY-7440	2.6.1	dTDP-3,4-didehydro-2,6-dideoxy-&alpha;-D-glucose aminotransferase
PWY-7440	5.1.3	dTDP-2,3,6-trideoxy-3-C-methyl-3-amino-4-dehydro-&alpha;-D-glucose 5-epimerase
PWY-7442	1.5.4.1	pyrimidodiazepine synthase
PWY-7443	1.1.1.275	(+)-trans-carveol dehydrogenase
PWY-7443	1.14.13.80	(R)-limonene 6-monooxygenase
PWY-7443	4.2.3.20	(4R)-limonene synthase
PWY-7444	2.4.1.189	MONOMER-14786
PWY-7444	2.4.1.190	MONOMER-14787
PWY-7444	2.4.1.191	MONOMER-14788
PWY-7445	1.11.1.7	luteolin diglucuronide peroxidase
PWY-7445	3.2.1.31	luteolin glucuronides &beta;-glucuronidase
PWY-7446	1.1.1.373 	3-sulfolactaldehyde reductase
PWY-7446	2.7.1.184	6-deoxy-6-sulfofructose kinase
PWY-7446	4.1.2.57	6-deoxy-6-sulfofructose-1-phosphate aldolase
PWY-7446	5.3.1.31 	sulfoquinovose isomerase
PWY-7447	1.13.11.78	2-amino-1-hydroxyethylphosphonate dioxygenase
PWY-7447	1.14.11.46	MONOMER-18484
PWY-7448	2.3.1	epicatechin:1-O-galloyl-&beta;-D-glucose O-galloyltransferase
PWY-7449	2.3.1.213	cyanidin 3-O-(6-O-glucosyl-2-O-xylosylgalactoside) 6-O-hydroxycinnamoyltransferase
PWY-7449	2.4.1.294	cyanidin 3-O-galactosyltransferase
PWY-7449	2.4.2.50	cyanidin 3-O-galactoside 2-O-xylosyltransferase
PWY-7450	2.3.1	anthocyanin 5-O-glucoside 6\-O-malonyltransferase"
PWY-7450 	2.3.1	anthocyanin acyltransferase
PWY-7450	2.3.1 	sinapoyl-&beta;-D-glucose:anthocyanin sinapoyltransferase
PWY-7450	2.4.1	acyl-glucose dependent anthocyanin glucosyltransferase
PWY-7450	2.4.1	anthocyanin 5-O-glucosyltransferase
PWY-7450 	2.4.2 	cyanidin 3-O-glucoside 2-O\-xylosyltransferase"
PWY-7452 	2.3.1.171	anthocyanidin 3-O-glucoside 6\-O-malonyltransferase"
PWY-7452	2.3.1.171	anthocyanidin 3-O-glucoside 6\-O-malonyltransferase"
PWY-7452	2.3.1 	anthocyanidin 3-O-glucoside-3\,6\"-O-dimalonyltransferase"
PWY-7455	1.1.1.213 	Aldo-keto reductase family 1 member C1
PWY-7455 	1.1.1.213 	Aldo-keto reductase family 1 member C3
PWY-7455	1.1.1.213 	Aldo-keto reductase family 1 member C4
PWY-7455	1.1.1.52 	3&alpha; hydroxysteroid dehydrogenase III
PWY-7455 	1.3.1.22	3-oxo-5-alpha-steroid 4-dehydrogenase 1
PWY-7455 	1.3.1.22	3-oxo-5-alpha-steroid 4-dehydrogenase 2
PWY-7456	2.4.1.281	4-O-&beta;-D-mannosyl-D-glucose phosphorylase subunit
PWY-7456	2.4.1.319	&beta;-1,4-mannooligosaccharide phosphorylase subunit
PWY-7456	3.2.1.100	MONOMER-18520
PWY-7456	5.1.3.11	MONOMER-18524
PWY-7458 	2.4.1	anthocyanidin 3-O- glucoside 6\-O-rhamnosyltransferase"
PWY-7458 	2.4.1 	anthocyanidin 3-O-glucoside 7-O-glucosyltransferase (acyl-glucose dependent)
PWY-7459	2.4.1.230	kojibiose phosphorylase
PWY-7459	5.4.2.6	&beta;-phosphoglucomutase
PWY-7460	2.3.1	meso-tartrate p-coumaroyltransferase
PWY-7460	2.3.1	trans-coutarate acetyltransferase
PWY-7461	2.3.1.131	sinapoyl-CoA: glucarate sinapoyltransferase
PWY-7461	2.3.1.132	sinapoyl-CoA: glucarolactone sinapoyltransferase
PWY-7461	2.3.1	glucarate acyltransferase (acyl-glucose dependent)
PWY-7462	1.13.11	3-mercaptopropionate dioxygenase
PWY-7462	1.8.1.4 	dihydrolipoamide dehydrogenase
PWY-7462	3.13.1.4	3-sulfinopropionyl-CoA desulfinase
PWY-7462	6.2.1.5 	succinyl-CoA synthetase
PWY-7464 	2.3.1	anthocyanin 7-O-glucoside acyltransferase (acyl-glucose dependent)
PWY-7464 	2.4.1	anthocyanin 7-O-(6-O-(p-hydroxybenzoyl)-glucoside) glucosyltransferase (acyl-glucose dependent)
PWY-7465	1.13.11	3-mercaptopropionate dioxygenase
PWY-7465	2.3.1	acyl-CoA:3-sulfinopropionate CoA-transferase
PWY-7466 	1.14.14.1	cytochrome P450 2E1
PWY-7467	2.6.1	UDP-acetamido-4-amino-6-deoxygalactopyranose transaminase subunit
PWY-7467	2.7.8	MONOMER-18558
PWY-7467	2.7.8	undecaprenyl-phosphate 2-acetamido-4-amino-2,4,6-trideoxy-&alpha;-D-galactose phosphotransferase
PWY-7468 	2.4.1	salicylic acid / benzoic acid glucosyltransferase
PWY-7470	2.3.1.23	MONOMER-18572
PWY-7470	2.3.1.23	MONOMER-18574
PWY-7471	1.1.1.254	D-carnitine dehydrogenase subunit
PWY-7472	5.1.2.M1	carnitine racemase
PWY-7473 	1.14.13	cytochrome P450 dependent &beta;-amyrin monooxygenase
PWY-7473 	2.4.1	avenacin A glucosyltransferase
PWY-7475 	2.3.1	SCPL desacyl avenacin acyltransferase
PWY-7475 	2.4.1	UDP-glucose: acyl glucosyltransferase
PWY-7476 	2.1.1.111	anthranilate N-methyltransferase
PWY-7477 	1.14.13 	cytochrome P450 dependent ent-sandaracopimardiene 7&beta;-hydroxylase
PWY-7478	1.14.13 	cytochrome P450 dependent ent-sandaracopimardiene 9&beta;-hydroxylase
PWY-7480	1.13.11.14	MONOMER-18584
PWY-7480	4.1.1	MONOMER-18585
PWY-7481	1.14.13.143	ent-isokaurene-C2-hydroxylase
PWY-7481	4.2.3.103	ent-isokaurene synthase
PWY-7482	1.14.13	cyclooctat-9-en-5,7-diol C18-monooxygenase
PWY-7482	1.14.13	octat-9-en-7-ol 5-monooxygenase
PWY-7482	4.2.3.146	cyclooctat-9-en-7-ol synthase
PWY-7483	1.13.12.21	tetracenomycin F1 monooxygenase
PWY-7483	1.14.13.200	tetracenomycin B2  oxygenase
PWY-7483	2.1.1	12-demethyl-elloramycin C12a O-methyltransferase
PWY-7483	2.1.1.305	8-demethyl-8-&alpha;-L-rhamnosyl tetracenomycin-C 2-O-methyltransferase
PWY-7483	2.1.1.306	8-demethyl-8-(2-methoxy-&alpha;-L-rhamnosyl)-tetracenomycin-C 3-O-methyltransferase
PWY-7483	2.1.1.307	8-demethyl-8-(2,3-dimethoxy-&alpha;-L-rhamnosyl)-tetracenomycin-C 4-O-methyltransferase
PWY-7483	2.1.1	tetracenomycin B3 9-O-methyltransferase
PWY-7483	2.1.1	tetracenomycin D3 3-O-methyltransferase
PWY-7483	2.3.1	elloramycin polyketide synthase complex
PWY-7483	2.4.1.331	8-demethyltetracenomycin C rhamnosyltransferase
PWY-7483	4.2.1.154	tetracenomycin F2 cyclase
PWY-7483	4.2.1	elloramycin aromatase
PWY-7483	4.2.1	elloramycin third ring cyclase
PWY-7484	1.14.13.145	ent-cassadiene 11&alpha;-hydroxylase
PWY-7484 	1.14.13 	cytochrome P450 dependent diterpene olefin hydroxylase
PWY-7484	1.14.13	ent-cassadiene-C2-hydroxylase
PWY-7485	1.13.12.21	tetracenomycin F1 monooxygenase
PWY-7485	1.14.13.200	MONOMER-18607
PWY-7485	2.1.1	tetracenomycin B3 8-O-methyl transferase
PWY-7485	2.1.1	tetracenomycin E C-9 O-methyltransferase
PWY-7485	2.3.1	tetracenomycin C polyketide synthase complex
PWY-7485	4.2.1.154	tetracenomycin F2 cyclase
PWY-7485	4.2.1	MONOMER-18609
PWY-7485	4.2.1 	multifunctional aromatase/3-O-methyl transferase TcmN
PWY-7487	2.4.1	(+)-secoisolariciresinol glucosyltransferase
PWY-7489 	1.1.1 	momilactone A synthase
PWY-7489	1.1.1	oryzalexin A synthase
PWY-7489	1.1.1	oryzalexin B synthase
PWY-7490	1.1.1.97	MONOMER-18630
PWY-7490	1.1.1	MONOMER-18631
PWY-7490	1.14.13	cytochrome P450 CYP619C2
PWY-7490	1.14.13	cytochrome P450 CYP619C3
PWY-7490	2.3.1.165	6-methylsalicylate synthase subunit
PWY-7490	4.1.1.52	MONOMER-18627
PWY-7491	2.4.1	(6-methoxy)podophyllotoxin 7-O-glucosyltransferase
PWY-7491	3.2.1	6-methoxypodophyllotoxin 7-O-glucoside glucosidase
PWY-7491	3.2.1	podophyllotoxin 7-O-glucoside glucosidase
PWY-7492	1.14.13	FAD-dependent monooxygenase, PaxM
PWY-7492	2.5.1.29	MONOMER-18634
PWY-7493	1.14.13	cytochrome P450 monooxygenase, PaxQ
PWY-7493	1.14.13	multifunctional cytochrome P450 monooxygenase, PaxP
PWY-7493	2.5.1	paxilline prenyltransferase
PWY-7494	1.1.3.17	choline oxidase
PWY-7494	1.2.1.8	MONOMER-18647
PWY-7495	1.14.99	flavonoid 8-O-hydroxylase
PWY-7495	2.1.1.88	8-hydroxyflavonol 8-O-methyltransferase
PWY-7496	3.5.1	linuron hydrolase subunit
PWY-7496	3.5.1	MONOMER-18671
PWY-7497	1.14.13	parthenolide synthase
PWY-7497	1.14.13	sesquiterpenoid 3&beta;-hydroxylase
PWY-7498	2.1.1.104 	phenylpropanoid/flavonoid O-methyltransferase
PWY-7501	2.7.8	phosphatidylserine synthase 1
PWY-7507 	1.14.13.37 	(S)-cis-N-methylstylopine 14-hydroxylase/(S)-N-methylcanadine 14-monooxygenase
PWY-7507 	1.14.13 	MONOMER-12344
PWY-7507 	1.3.1.107	sanguinarine reductase
PWY-7507 	1.5.3.12	MONOMER-12346
PWY-7509	2.7.8.41	bifunctional cardiolipin / phosphatidylethanolamine synthase
PWY-7510	2.6.1	L-&alpha;-amino acid transaminase
PWY-7510	4.1.1.82	beta subunit 
PWY-7510	4.1.3	phosphonoacetaldehyde aldolase
PWY-7510	5.4.2.9	phosphoenolpyruvate mutase
PWY-7510	6.3.2	L-&alpha;-amino acid ligase
PWY-7511	2.3.2.23	ubiquitin-conjugating enzyme E2
PWY-7511	2.3.2.23	ubiquitin-conjugating enzyme E2D2
PWY-7511	2.3.2.24	(E3-independent) ubiquitin-conjugating enzyme E2 O
PWY-7511	2.3.2.26	HECT-type E3 ubiquitin transferase NEDD4L
PWY-7511	2.3.2.27	RING-type E3 ubiquitin transferase DOA10
PWY-7511	2.3.2.27	RING-type E3 ubiquitin transferase MARCH6
PWY-7511	2.3.2.27	RING-type E3 ubiquitin transferase UBR1
PWY-7511	6.2.1.45	E1 ubiquitin-activating enzyme
PWY-7511	6.2.1.45	MONOMER-16752
PWY-7511	6.2.1.45	MONOMER-16753
PWY-7511	6.2.1.45	MONOMER-16754
PWY-7512	1.14.14	MONOMER-18698
PWY-7512	3.5.99	MONOMER-18699
PWY-7513	1.14.13 	1,3,6,8-tetrahydroxynaphthalene monooxygenase
PWY-7513	1.14.21.7	biflaviolin synthase
PWY-7513	2.3.1.233	1,3,6,8-tetrahydroxynaphthalene synthase
PWY-7513	2.3.1.233	MONOMER-14457
PWY-7514	4.1.1.83	4-hydroxyphenylacetate decarboxylase small subunit 
PWY-7515	1.5.1.49 	&Delta;<sup>1</sup>-pyrroline-2-carboxylate reductase subunit
PWY-7515	4.2.1.77	trans-3-hydroxy-L-proline dehydratase subunit
PWY-7516	4.2.1.43	MONOMER-18708
PWY-7516	4.2.1	MONOMER-18707
PWY-7517	1.14.11	fusicocca-2,10(14)-diene-8&beta;,16-diol dioxygenase
PWY-7517	1.14.13	fusicoccadiene 8-ol C-16 hydroxylase
PWY-7517	1.14.13	fusicoccadiene C-8 hydroxylase
PWY-7517	2.1.1	fusicocca-1,10(14)-diene-3,8&beta;,16-triol O-methyltransferase
PWY-7517	4.2.3.43 	fusicoccadiene synthase
PWY-7518	2.6.1.57 	L-tyrosine:2-oxoglutarate aminotransferase
PWY-7520	1.1.1	quinone reductase
PWY-7520	2.5.1	bisindolylquinone prenyltransferase
PWY-7520	2.6.1.28	L-tryptophan:phenylpyruvate aminotransferase
PWY-7521	1.14.13	MONOMER-18731
PWY-7521	3.1.1	MONOMER-18732
PWY-7523	3.5.3.7	MONOMER-18734
PWY-7524	2.7.1.185	MONOMER-18735
PWY-7524	2.7.1.186	MONOMER-18736
PWY-7524	2.7.4.26	isopentenyl phosphate kinase subunit
PWY-7527 	2.4.2.28	S-methyl-5-thioadenosine phosphorylase
PWY-7527 	2.6.1 	kynurenine aminotransferase I
PWY-7527 	3.5.1.111	&omega;-amidase NIT2
PWY-7527 	5.3.1.23	methylthioribose-1-phosphate isomerase
PWY-7528 	1.13.11.54 	acireductone dioxygenase
PWY-7528 	2.5.1.16	SPEEBACSU-MONOMER
PWY-7528 	2.5.1.6	S-adenosylmethionine synthetase monomer
PWY-7528 	2.6.1.88 	2-oxo-4-methylthiobutanoate-glutamine aminotransferase
PWY-7528 	2.6.1.88 	methionine-oxo-acid transaminase
PWY-7528 	2.7.1.100	MONOMER-10424
PWY-7528 	2.7.1.100	MONOMER-10461
PWY-7528 	2.7.1.100	MONOMER-1295
PWY-7528 	3.1.3.77 	MONOMER-1331
PWY-7528 	3.1.3.87	2-hydroxy-3-keto-5-methylthiopentenyl-1-phosphate phosphatase
PWY-7528 	3.2.2.16	MTA nucleosidase subunit
PWY-7528 	3.2.2.16 	MTA/SAH nucleosidase
PWY-7528 	3.5.1.111	2-oxoglutaramate amidase
PWY-7528 	4.1.1.50	adenosylmethionine decarboxylase proenzyme
PWY-7528 	5.3.1.23	5-methylthioribose-1-phosphate isomerase
PWY-7528 	5.3.2.5	2,3-diketo-5-methylthiopentyl-1-phosphate enolase subunit
PWY-7530	2.4.1.M15	undecaprenyldiphospho-N-acetylgalactosamine 3-&beta;-galactosyltransferase
PWY-7531	1.1.1	MONOMER-18759
PWY-7531	1.1.1	MONOMER-18760
PWY-7531	2.6.1	MONOMER-18757
PWY-7531	3.1.3	MONOMER-18758
PWY-7532	2.3.1	aszonalenin acetyltransferase
PWY-7532	2.5.1	benzodiazepinedione prenyltransferase
PWY-7532	6.3.2	benzodiazepinedione synthetase
PWY-7533	1.14.13	cyclo(L-Phe-L-Ser) hydroxylase
PWY-7533	1.8.3.2	dithiol oxidase
PWY-7533	2.1.1	N-desmethyl-gliotoxin N-methyltransferase
PWY-7533	2.3.2	&gamma;-glutamylcyclotransferase
PWY-7533	2.5.1.18	glutathione S-transferase
PWY-7533	3.4.13.19	dipeptidase
PWY-7533	4.4.1	carbon-sulfate lyase
PWY-7533	6.3.2	cyclo (L-Phe-L-Ser) synthetase
PWY-7534	2.1.1	dithiolgliotoxin S-methyltransferase
PWY-7535	1.14.13.198	dihydromonacolin L hydroxylase
PWY-7535	2.3.1.161	lovastatin nonaketide synthase complex
PWY-7535	2.3.1.244	2-methylbutanoate polyketide synthase
PWY-7535	2.3.1 	monacolin J acid methylbutanoate transferase
PWY-7535	3.1.2.31	dihydromonacolin L-[lovastatin nonaketide synthase] thioesterase
PWY-7536	2.3.1.37	bifunctional 5-aminolevulinate synthase / 5-aminolevulinyl-CoA cyclase
PWY-7536	6.2.1	acyl-CoA ligase
PWY-7537 	1.14.11.38	verruculogen synthase
PWY-7537 	1.14.13.176	tryprostatin B 6-hydroxylase
PWY-7537 	1.14.13.177	fumitremorgin C monooxygenase
PWY-7537 	1.14.21.10	fumitremorgin C synthase
PWY-7537 	2.1.1.293	6-hydroxytryprostatin B O-methyltransferase
PWY-7537 	2.5.1.106	tryprostatin B synthase
PWY-7537 	2.5.1.107	verruculogen O-prenyltransferase
PWY-7537 	2.5.1.110	12&alpha;,13&alpha;-dihydroxyfumitremorgin C prenyltransferase
PWY-7537 	6.3.3	brevianamide F synthase
PWY-7539	2.5.1.15 	folate synthesis bifunctional protein
PWY-7539	3.5.4.25	GTP cyclohydrolase II
PWY-7539	4.1.2.25	dihydroneopterin aldolase
PWY-7539	5.3.1.24 	MONOMER-18791
PWY-7540	1.14.13	cytochrome P450 monooxygenase, AtmQ
PWY-7540	2.5.1	dimethylallyl-diphosphate:paspalinine dimethylallyltransferase, AtmD
PWY-7541 	1.2.1.87	CoA-dependent propanal dehydrogenase
PWY-7542	1.14.13	fumiquinazoline F monooxygenase
PWY-7542	1.5.3	fumiquinazoline C oxidase
PWY-7542	6.3.2	fumiquinazoline A synthetase
PWY-7542	6.3.2	fumiquinazoline F synthetase
PWY-7543	2.5.1	ardeemin FQ prenyltransferase
PWY-7543	6.3.2	ardeemin FQ synthetase
PWY-7544	1.2.5.1	pyruvate oxidase
PWY-7546	2.1.1.314	diphthine synthase
PWY-7546	2.5.1.108	2-(3-amino-3-carboxypropyl)histidine synthase
PWY-7546	3.1.1.97	methylated diphthine methylhydrolase
PWY-7546	6.3.1.14	diphthine&mdash;ammonia ligase
PWY-7547	1.1.1	MONOMER-18834
PWY-7547	1.3.8	L-prolyl-S-peptidyl-carrier protein dehydrogenase
PWY-7547	1.4.3	2-methyl-3-n-amyl-dihydropyrrole oxidase
PWY-7547	2.1.1	4-hydroxy-2,2-bipyrrole-5-carbaldehyde O-methyltransferase
PWY-7547	2.2.1.12 	3-acetyloctanal synthase
PWY-7547	2.6.1	(S)-3-acetyloctanal aminotransferase
PWY-7547	4.2.1 	PigH polyketide synthase
PWY-7547	6.2.1	L-proline--[PigG acyl-carrier protein] ligase
PWY-7547	6.4.1	prodigiosin synthase
PWY-7548	2.2.1	transaldolase, LmbR
PWY-7548	2.7.1	MONOMER-18843
PWY-7548	2.7.7	MONOMER-18842
PWY-7548	3.1.3	MONOMER-18841
PWY-7548	5.3.1	MONOMER-18840
PWY-7549	1.14.13.216	asperlicin C monooxygenase
PWY-7549	6.3.2	asperlicin C-synthetase
PWY-7550	1.14.99.51 	<small>L</small>-histidine N<small><sup>&alpha;</sup></small>-methyltransferase
PWY-7552	1.3.98.M2	MONOMER-18867
PWY-7552	4.1.1	siroheme decarboxylase
PWY-7552	4.1.3	MONOMER-18866
PWY-7554	1.1.1	MONOMER-18859
PWY-7554	1.2.98	MONOMER-18860
PWY-7554	4.1.1	siroheme decarboxylase
PWY-7555	1.21.99.1	&beta;-cyclopiazonate oxidocyclase
PWY-7555	1.21.99.1	&beta;-cyclopiazonate oxidocyclase I
PWY-7555	1.21.99.1	&beta;-cyclopiazonate oxidocyclase II
PWY-7555	1.21.99.1	&beta;-cyclopiazonate oxidocyclase III
PWY-7555	1.21.99.1	&beta;-cyclopiazonate oxidocyclase IV
PWY-7555	1.21.99.1	&beta;-cyclopiazonate oxidocyclase V
PWY-7555	2.5.1	cyclo-acetoacetyl-L-tryptophan dimethylallyltransferase
PWY-7555	6.3.2 	cyclo-acetoacetyl-L-tryptophan synthetase
PWY-7556	6.3.2	microperfuranone synthase
PWY-7557	1.1.1	dehydrodiconiferyl alcohol dehydrogenase
PWY-7557	1.13.11.43 	lignostilbene &alpha;,&beta;-dioxygenase III
PWY-7557	1.13.11 	lignostilbene &alpha;,&beta;-dioxygenase I
PWY-7557	1.2.1.67 	aromatic aldehyde dehydrogenase
PWY-7558	1.14.13	&alpha;-cyclopiazonate hydroxylase
PWY-7559	3.4.19	glutathione hydrolase
PWY-7560	1.1.1.267	1-deoxy-D-xylulose 5-phosphate reductoisomerase, apicoplast
PWY-7560	1.1.1.267	MONOMER-16159
PWY-7560	1.17.7.1	AT5G60600-MONOMER
PWY-7560	2.2.1.7	MONOMER-16158
PWY-7560	2.7.1.148	MONOMER-11959
PWY-7561	1.14.13	pretenellin B N-monooxygenase
PWY-7561	1.14.19	pretenellin A oxidase
PWY-7561	2.3.1	trans-acting enoyl reductase
PWY-7561	6.3.2 	tenellin synthetase
PWY-7562	1.1.1.127	2,5-didehydro-3-deoxy-L-galactonate 5-reductase
PWY-7562	1.1.1.389	2-dehydro-3-deoxy-L-galactonate 5-dehydrogenase
PWY-7562	1.2.1.92	3,6-anhydro-L-galactose dehydrogenase
PWY-7562	2.7.1.178 	2-dehydro-3-deoxygluconate kinase
PWY-7562	4.1.2.55 	2-dehydro-3-deoxy-phosphogluconate aldolase
PWY-7562	5.5.1.25	3,6-anhydro-L-galactonate cycloisomerase
PWY-7563	2.3.1	trans enoyl reductase
PWY-7563	6.3.2 	bassianin synthetase
PWY-7563	6.3.2 	desmethylbassianin synthetase
PWY-7564	2.1.1	MONOMER-18909
PWY-7564	2.1.2	MONOMER-18907
PWY-7564	2.7.6	MONOMER-18910
PWY-7564	3.2.2	5-hydroxymethylcytidine 5-monophosphate nucleosidase subunit
PWY-7565	1.14.13 	preaspyridone A oxidase
PWY-7565	2.3.1	trans enoyl reductase
PWY-7565	6.3.2 	aspyridone synthetase
PWY-7566	1.1.1.380	L-gulonate 5-dehydrogenase
PWY-7568	1.1.1	MONOMER-18914
PWY-7569	2.1.1	MONOMER-18923
PWY-7569	2.4.1	MONOMER-18919
PWY-7569	2.6.1	cytosyl-4-dehydro-&beta;-D-glucuronate transaminase
PWY-7569	3.2.2.10	nucleoside hydrolase
PWY-7569	6.3.2	MONOMER-18920
PWY-7569	6.3.2	MONOMER-18922
PWY-7570	1.1	cytosyl-&beta;-D-glucuronate 4-oxidoreductase
PWY-7570	2.1.1	L-leucyldemethyl-blasticidin S guanidino methyltransferase
PWY-7570	2.4.1	cytosyl-glucuronate synthase
PWY-7570	2.6.1	cytosyl-4-dehydro-&beta;-D-glucuronate transaminase
PWY-7570	3.2.2.10	nucleoside hydrolase
PWY-7570	3.4.13	MONOMER-18933
PWY-7570	4.1.1	cytosinine synthase
PWY-7570	5.4.3	L-arginine 2,3 aminomutase
PWY-7570	6.3.2	&beta;-L-arginine--cytosinine ligase
PWY-7570	6.3.2	demethylblasticidin S--L-leucine ligase
PWY-7571	1.14.13.195	L-ornithine 5-monooxygenase
PWY-7571	2.3.3.10	hydroxymethyl glutaryl-CoA synthase
PWY-7571	3.5.1	N-hydroxyornithine acylase
PWY-7571	4.2.1.18	methylglutaconyl-CoA hydratase
PWY-7571	6.3.2	ferrichrome A synthetase
PWY-7572	1.14.13 	multifunctional cytochrome P450 monooxygenase
PWY-7572	1.14.13	multifunctional cytochrome P450 monooxygenase
PWY-7572	1.14.21	multifunctional cytochrome P450 monooxygenase
PWY-7572	2.5.1	multifunctional aromatic O-prenyltransferase
PWY-7572	2.5.1	multifunctional diprenyltransferase
PWY-7574	1.1.1.59	MONOMER-18957
PWY-7574	1.2.1	malonate semialdehyde dehydrogenase (CoA-acylating)
PWY-7574	3.1.2.4	MONOMER-18956
PWY-7574	4.2.1.119 	D-specific enoyl-CoA hydratase
PWY-7575 	1.14.13	MONOMER-18951
PWY-7575 	2.4.1	MONOMER-18944
PWY-7575 	2.6.1.85	4-amino-4-deoxychorismate synthase
PWY-7575 	2.6.1	MONOMER-18953
PWY-7575	4.1.3.38	MONOMER-18662
PWY-7575	4.1.3.38	MONOMER-18663
PWY-7575 	4.2.1.47	MONOMER-18954
PWY-7575 	6.2.1 	candicidin polyketide synthase complex
PWY-7575 	6.2.1 	MONOMER-18943
PWY-7576	1.19.6.1	[FeFe]-nitrogenase complex
PWY-7576	1.19.6.1	[FeMo]-nitrogenase complex
PWY-7577	2.3.1	N-hydroxyornithine: acetyl CoA N-transacetylase
PWY-7577	6.3.2	ferrichrome synthetase
PWY-7578 	1.14.15.20	heme oxygenase 1
PWY-7578 	1.3.7.5	MONOMER-18997
PWY-7578	4.4.1.31	phycoerythrocyanin &alpha; subunit phycoviolobilin lyase/isomerase
PWY-7579	1.14.15.20	MONOMER-18993
PWY-7579	1.3.7.2	MONOMER-18994
PWY-7579	1.3.7.3	MONOMER-18995
PWY-7579	4.4.1.33	R-phycocyanin V &alpha;-subunit:phycourobilin lyase/isomerase
PWY-7580	1.14.15.20	heme oxygenase
PWY-7580	1.3.7.6	MONOMER-13958
PWY-7581	4.1.3.3	MONOMER-18999
PWY-7581	5.1.3.8	MONOMER-19001
PWY-7582	1.13.11	mercaptosuccinate dioxygenase subunit
PWY-7583	2.3.1	polyunsaturated fatty acid synthase
PWY-7585	2.3.1	docosahexaenoate synthase
PWY-7586	2.4.1 	&beta;-1,4-D-mannosyl-N-acetyl-D-glucosamine phosphorylase
PWY-7588	1.1.1.201	7&beta;-hydroxysteroid dehydrogenase
PWY-7589	1.14.19.27	sn-2 acyl-lipid 9-desaturase
PWY-7590 	1.14.19.23	plastidial acyl-lipid &omega;-6 desaturase
PWY-7590 	1.14.19.35 	chloroplastic acyl-lipid &omega;-3 desaturase
PWY-7590 	1.14.19.42	palmitoyl-monogalactosyldiacylglycerol 7-desaturase
PWY-7591	2.1.1	MONOMER-19028
PWY-7591	4.2.1	MONOMER-19026
PWY-7591	5.5.1.19	MONOMER-19025
PWY-7596 	1.14.19.28	sn-1 stearoyl-lipid 9-desaturase
PWY-7596 	1.14.19.45	sn-1 oleoyl-lipid 12-desaturase
PWY-7596 	1.14.19.46	sn-1 linoleoyl-lipid 6-desaturase
PWY-7596 	2.3.1.51	1-acylglycerol-3-phosphate O-acyltransferase
PWY-7596 	2.3.1.51	MONOMER-19031
PWY-7598 	1.14.19.36	sn-1 &gamma;-linolenoyl-lipid &omega;-3 desaturase
PWY-7599	1.14.11 	MONOMER-19047
PWY-7599	1.14.11	MONOMER-19056
PWY-7599	1.14.13	MONOMER-19044
PWY-7599	1.14.13	MONOMER-19049
PWY-7599	1.1.99	MONOMER-19046
PWY-7599	1.1.99	MONOMER-19053
PWY-7599	1.1.99	MONOMER-19055
PWY-7599	2.3.1	MONOMER-19041
PWY-7599	2.3.1	MONOMER-19054
PWY-7599	2.5.1	MONOMER-19043
PWY-7599	3.1.1 	bifunctional 3,5-dimethylorsellinate hydroxylase / 5,7-dihydroxy-4,6-dimethylphthalide synthase
PWY-7599	5.4.99	MONOMER-19045
PWY-7600	6.3.2	peramine synthetase
PWY-7601 	1.14.19.30	acyl-lipid 5-desaturase
PWY-7601 	1.14.19.4	acyl-lipid 8-desaturase
PWY-7601 	2.3.1.199	C18 &Delta;9-polyunsaturated fatty acyl-CoA elongase
PWY-7602 	2.3.1.199	C18 &Delta;9-polyunsaturated fatty acyl-CoA elongase
PWY-7603	2.5.1.109	brevianamide F reverse prenyltransferase
PWY-7603	2.5.1.109	deoxybrevianamide F synthase
PWY-7603	2.5.1	6-hydroxy-deoxybrevianamide E prenyltransferase
PWY-7604	1.14.13	notoamide E monooxygenase
PWY-7605 	1.14.21	cytochrome P450 roquefortine oxidoreductase
PWY-7605 	2.5.1	roquefortine synthase
PWY-7605 	6.3.2	histidyltryptophanyldiketopiperazine synthetase
PWY-7606 	2.3.1	polyunsaturated long-chain fatty acid elongase 2
PWY-7606 	2.3.1 	polyunsaturated long-chain fatty acid elongase 5
PWY-7608 	1.14.13	cytochrome P450 glandicoline/roquefortine monooxygenase
PWY-7608 	2.1.1	roquefortine/glandicoline O- methyltransferase
PWY-7609 	1.14.13	roquefortine monooxygenase
PWY-7610	1.1.1	MONOMER-19060
PWY-7610	1.1.1	MONOMER-19061
PWY-7610	4.2.1	GDP-D-glycero-&alpha;-D-manno-heptose 4,6-dehydratase subunit
PWY-7610	5.1.3	MONOMER-19059
PWY-7612	1.14.13	cytochrome P450 dependent chaetoglobosin monooxygenase
PWY-7612	1.14.13	cytochrome P450 dependent prochaetoglobosin monooxygenase
PWY-7612	1.14.21	FAD-dependent chaetoglobosin oxidoreductase
PWY-7612	2.3.1	enoyl reductase
PWY-761	2.6.1.76	MONOMER-15538
PWY-7612	6.3.2 	chaetoglobosin synthetase
PWY-7612	6.3.2 	prochaetoglobosin I synthetase
PWY-7613	1.1.1	GDP-6-deoxy-4-keto-D-lyxo-heptose 4-reductase subunit
PWY-7613	4.2.1	MONOMER-19080
PWY-761	4.1.1.86	MONOMER-15539
PWY-7614 	4.4.1.4	alliinase
PWY-7615	1.1.1	pterocarpan:NADPH oxidoreductase
PWY-7615	1.14.13	phaseollin 1&alpha;-hydroxylase
PWY-7615	1.14.13	pisatin demethylase
PWY-7615	1.14.13	pterocarpan 1&alpha;-hydroxylase
PWY-7615	4.2.1.97	phaseollidin hydratase
PWY-7616 	1.1.1.1 	mycothiol-dependent formaldehyde dehydrogenase
PWY-7616	1.1.1.1 	Zn<sup>2+</sup>-dependent alcohol dehydrogenase
PWY-7616	1.1.99.33	G18NG-10091-MONOMER
PWY-7616	1.2.1.3 	acetaldehyde dehydrogenase
PWY-7618 	1.14.18.4	oleate 12-hydroxylase
PWY-7618 	1.14.19.22 	oleoyl-lipid 12-hydroxylase/desaturase
PWY-7619 	1.14.19.37	acyl-CoA 5-desaturase
PWY-7620	4.1.1	naphthalene carboxylase
PWY-7620	6.2.1	2-naphthoate-CoA ligase
PWY-7621	1.1.1	MONOMER-19106
PWY-762	1.14.19.43	palmitoyl-[glycerolipid] 3-(E) desaturase
PWY-7621	2.3.1	3-aminotridec-2-en-4-one synthase
PWY-7622	5.1.3.2	MONOMER-19119
PWY-7622	5.1.3.2	MONOMER-19120
PWY-7622	5.4.99.9	UDP-galactopyranose mutase subunit
PWY-7623	2.4.1	MONOMER-19108
PWY-7624	1.14.13	MONOMER-19116
PWY-7624	1.14.13	MONOMER-19117
PWY-7624	2.3.1	Nys PKS complex
PWY-7624	2.4.1	MONOMER-19118
PWY-7626	1.1.1.385	dihydroanticapsin dehydrogenase
PWY-7626	1.3.1.aa	3-[(5R)-5-hydroxy-7-oxabicyclo[4.1.0]heptan-2-ylidene]-2-oxopropanoate reductase
PWY-7626	2.6.1	3-[(2S,5R)-5-hydroxy-7-oxabicyclo[4.1.0]heptan-2-yl]-2-oxopropanoate aminotransferase
PWY-7626	4.1.1.100	MONOMER-19126
PWY-7626	5.3.3.19	3-[(4R)-4-hydroxycyclohexa-1,5-dien-1-yl]-2-oxopropanoate isomerase
PWY-7626	6.3.2.49	L-alanine-anticapsin ligase
PWY-7627	1.3.98 	2,4,6-trinitrophenol reductase
PWY-7627	1.3.98	MONOMER-19137
PWY-7627	5.3.99	MONOMER-19139
PWY-7628	3.3.2.14	2,4-dinitroanisole o-demethylase
PWY-7629	2.1.1	5-methylthujaplicatin 4-O-methyltransferase
PWY-7629	2.1.1.M19	thujaplicatin 5-O-methyltransferase
PWY-7630 	1.14.21.11	pluviatolide synthase
PWY-7630	1.14.21.11	pluviatolide synthase
PWY-7631	2.1.1	matairesinol 4-O-methyltransferase
PWY-7631	2.1.1	matairesinol 4-O-methyltransferase
PWY-7633	2.4.1	calycosin 7-O-glucosyltransferase
PWY-7634	1.14.13	linoleate 12-epoxygenase
PWY-7635	4.2.1.95	kievitone hydratase
PWY-7636 	1.14.13.129 	carotenoid &beta;-ring 4-dehydrogenase
PWY-7636	1.14.13.129 	carotenoid &beta;-ring 4-dehydrogenase
PWY-7636	1.1.99	carotenoid 4-hydroxy-&beta;-ring 4-dehydrogenase
PWY-7637	1.14.13	carotenoid 2,2-&beta;-ionone ring hydroxylase
PWY-7638	1.14.13	&beta;-carotene hydroxylase
PWY-7638	1.5.3	&beta;-carotene 4-monoketolase
PWY-7639 	2.7.7.77	EG11829-MONOMER
PWY-7640 	1.14.13 	(+)-abscisic acid 8&prime;-hydroxylase
PWY-7640 	1.14.13 	AT2G29090-MONOMER
PWY-7640 	1.14.13 	AT3G19270-MONOMER
PWY-7640 	1.14.13 	AT5G45340-MONOMER
PWY-7641	1.14.19.40	5-hexenoyl-[acp] acetylenase
PWY-7641	6.2.1.l	5-hexenoate:acyl-carrier protein ligase
PWY-7642	1.14.13	abscisic acid 7-hydroxylase
PWY-7643	2.1.1	coniferyl alcohol 9-methyltransferase subunit
PWY-7643	2.1.1	coniferyl alcohol 9-O-methyltransferase
PWY-7644 	3.10.1.1	N-sulfoglucosamine sulfohydrolase
PWY-7644	3.1.6.11	disulfoglucosamine-6-sulfatase
PWY-7644	3.1.6.13	iduronate-2-sulfatase
PWY-7644	4.2.2.7	heparinase I
PWY-7645	4.2.2.1	hyaluronate lyase
PWY-7646 	3.1.6.9	chondro-4-sulfatase
PWY-7646 	3.2.1.179 	unsaturated &beta;-glucuronyl hydrolase
PWY-7646 	4.2.2.20 	chondroitin ABC endo-lyase
PWY-7646	4.2.2.20 	chondroitin-B lyase
PWY-7646 	4.2.2.21	chondroitin sulfate ABC exolyase
PWY-7647	3.2.1	oligo-ulvan &beta;-glucuronyl hydrolase
PWY-7647	4.2.2	ulvan lyase
PWY-7648	1.14.11	leucine 5-hydroxylase
PWY-7649	1.1.1	3-(4-hydroxybenzyl)-malate dehydrogenase
PWY-7649	1.14.11	L-homotyrosine 3-hydroxylase
PWY-7649	2.3.3	2-(4-hydroxybenzyl)-malate synthase
PWY-7649	2.6.1	4-(4-hydroxyphenyl)-2-oxobutanoate transaminase
PWY-7649	5.4.4	2-(4-hydroxybenzyl)-malate isomerase
PWY-7650	1.14.13	echinocandine D hydroxylase
PWY-7650	6.1	linoleyl-AMP ligase
PWY-7650	6.3.2	echinocandin synthetase
PWY-7651 	3.2.1.56 	unsaturated &beta;-glucuronyl hydrolase
PWY-7651 	4.2.2.8 	heparinase II
PWY-7651	4.2.2.8	heparinase III
PWY-7652	3.5.1	echinocandin B deacylase
PWY-7653	1.14.19	griseophenone C halogenase
PWY-7653	2.1.1	2-(2,4-dihydroxy-6-methylbenzoyl)benzene-1,3,5-triol 3-O-methyltransferase
PWY-7653	2.1.1	griseophenone D 9-O-methyltransferase
PWY-7653	2.3.1	norlichexanthone synthase
PWY-7655 	1.14.21	griseophenone synthase
PWY-7655 	1.3.1	dehydrogriseofulvin reductase
PWY-7655 	2.1.1	griseofulvin 5-O-methyltransferase
PWY-7656	1.14.19	(11E)-tetradecenoyl-CoA 9Z desaturase
PWY-7656	1.14.19.5 	acyl-CoA 11-desaturase
PWY-7657	1.1.1.384	dTDP-3,4-didehydro-2,6-dideoxy-&alpha;-D-glucose reductase
PWY-7657	1.1.1	NDP-4-keto-6-deoxyhexose 4-ketoreductase
PWY-7657	2.7.7.24	glucose-1-phosphate thymidylyltransferase
PWY-7657	4.2.1.46	dTDP-glucose 4,6-dehydratase
PWY-7657	5.1.3	dTDP-4-keto-6-deoxyhexose 5-epimerase
PWY-7658	2.4.1.153	UDP-N-acetylglucosamine&ndash;dolichyl-phosphate N-acetylglucosaminyltransferase
PWY-7658	2.4.1.335	UDP-2,3-diacetamido-2,3-dideoxy-&alpha;-D-glucuronate glycosyltransferase
PWY-7658	2.4.1	GDP-2-acetamido-2-deoxy-&beta;-mannuronate glycosyltransferase
PWY-7658	2.4.99.21	dolichyl-monophosphooligosaccharide--protein glycotransferase
PWY-7659	1.14.13	12a-deshydroxy-desmethylanthrotainin 12a-hydroxylase
PWY-7659	1.14.13	desmethylanthrotainin 5-hydroxylase
PWY-7659	1.14.21	viridicatumtoxin synthase
PWY-7659	2.1.1	5-hydroxy-desmethylanthrotainin 8-O-methyltransferase
PWY-7659	2.3.1	4-(3-acetyl-2,5,7,10-tetrahydroxy-4-oxo-3-hydroanthracen-2-yl)-3-oxobutanamide-ACP synthase
PWY-7659	2.5.1.1	geranyl diphosphate synthase
PWY-7659	2.5.1	5-hydroxyanthrotainin geranyltransferase
PWY-7659	6.2.1	malonamate-CoA ligase
PWY-7660	1.1.1	tryptoquialanone dehydrogenase
PWY-7660	1.14.13	15-dimethyl-2-epi-fumiquinazoline A oxidase
PWY-7660	1.14.13	deoxynortryptoquialanone monooxygenase
PWY-7660	1.14.13	fumiquinazoline F monooxygenase
PWY-7660	2.3.1	tryptoquialanol O-acetyltransferase
PWY-7660	3.1.1	deoxynortryptoquialanone lactonohydrolase
PWY-7660	6.2.1 	fumiquinazoline F-indoline-2-2-aminobutyryl-3-ol synthetase
PWY-7661	2.1.1	phosphodolichol-tetrasaccharide glucuronate methyltransferase
PWY-7661	2.4.1.83	MONOMER-19290
PWY-7661	2.4.1	AglE glucuronate-transferase
PWY-7661	2.4.1	AglG glucuronate-transferase
PWY-7661	2.4.1	AglI glucuronate-transferase
PWY-7661	2.4.1	dolichol phosphate hexosyltransferase
PWY-7661	2.4.1	dolichyl-monophosphomannose--[protein]-tetrasaccharide mannosyltransferase
PWY-7661	2.4.99.21	MONOMER-19288
PWY-7662	4.2.1.110 	MONOMER-19294
PWY-7662	4.2.1.111	1,5-anhydro-<small>D</small>-fructose dehydratase
PWY-7662	4.2.2.13	&alpha;-1,4-glucan lyase 1
PWY-7662	5.3.2.7	ascopyrone tautomerase
PWY-7663 	1.1.1.100	3-oxoacyl-[acyl-carrier-protein] reductase
PWY-7663 	1.3.1.9	enoyl-[acyl-carrier-protein] reductase
PWY-7663 	2.3.1.41	3-oxoacyl-[acyl-carrier-protein] synthase
PWY-7663 	4.2.1.59 	3-hydroxyacyl-[acyl-carrier-protein] dehydratase/isomerase
PWY-7665	6.3.2	aureobasidin A synthetase
PWY-7666	2.4.1.241	digalactosyldiacylglycerol synthase
PWY-7666	2.4.1.336	monoglucosyldiacylglycerol synthase
PWY-7666	5.1.3.34	&beta;-monoglucosyldiacylglycerol epimerase
PWY-7667	1.1.1	L-2-amino-8-hydroxdecanoate dehydrogenase
PWY-7667	1.14.13	L-2-aminodecanoate 8-hydroxylase
PWY-7667	1.5.1	1-piperideine 6-carboxylate reductase
PWY-7667	2.1.1	N-hydroxytryptophan O-methyltransferase
PWY-7667	6.3.2	apicidin synthetase
PWY-7668	1.5.1	L-piperideine 6-carboxylate reductase
PWY-7668	2.1.1	N-hydroxy-L-tryptophan O-methyltransferase
PWY-7668	6.3.2	apicidin F synthetase
PWY-7669	2.1.1	trichosetin N-methyltransferase
PWY-7669	2.3.1	equisetin enoylreductase
PWY-7669	6.3.2 	equisetin synthetase
PWY-7670	2.3.1	trans-enoyl reductase
PWY-7670	6.3.2 	fusaridione A synthetase
PWY-7671	1.11.2.5	3-methyl-<small>L</small>-tyrosine peroxygenase
PWY-7671	2.1.1.304	L-tyrosine C<small><sup>3</sup></small>-methyltransferase
PWY-7671	2.1.1.M6	3-hydroxy-5-methyl-L-tyrosine 4O-methyltransferase
PWY-7671	2.3.1 	non-ribosomal peptide synthase SfmC
PWY-7671	6.3.2	non-ribosomal peptide synthase SfmA
PWY-7671	6.3.2 	non-ribosomal peptide synthase SfmB
PWY-7672	2.3.1	fusaric acid synthase
PWY-7673	1.14.13	fusarin cytochrome P450 monooxygenase
PWY-7673	1.14.13	fusarin oxygenase
PWY-7673	2.1.1	carboxy-fusarin C O-methyltransferase
PWY-7673	6.3.2 	fusarin C synthetase
PWY-7674	1.1.3.48	3-deoxy-D-manno-octulosonate 8-oxidase
PWY-7674	2.6.1.109	8-amino-3,8-dideoxy-D-manno-octulosonate transaminase
PWY-7674	2.7.7.90	8-amino-3,8-dideoxy-manno-octulosonate cytidylyltransferase
PWY-7675	2.7.1.166	MONOMER-15507
PWY-7676	2.4.99	8-amino-3,8-dideoxy-D-manno-octulosonate transferase
PWY-7676	2.7.1.166	8-amino-3,8-dideoxy-D-manno-octulosonate kinase
PWY-7677	1.14.13	20-deoxo-20-dihydro-12,13-deepoxyrosamicin 12,13 epoxidase
PWY-7677	1.14.13	MONOMER-19354
PWY-7678	2.4.1.295	anthocyanidin 3-O-sambubioside 5-O-glucosyltransferase
PWY-7679	2.3.1 	anthocyanidin glucoside 6-O-acyltransferase
PWY-7679 	2.4.2.51	anthocyanidin 3-O-glucoside 2-O-xylosyltransferase
PWY-7680	1.14.13.190	ferruginol synthase
PWY-7680	4.2.3.131	multiradiene synthase
PWY-7680	5.5.1.12	copalyl diphosphate synthase
PWY-7681	1.14.12.23	nitroarene dioxygenase complex
PWY-7682	1.13.11	extradiol ring-cleavage dioxygenase
PWY-7682	1.14.13	p-coumaraldehyde 3-hydroxylase
PWY-7683	1.7.2	hemoglobin
PWY-7683	1.7.2	non-symbiotic hemoglobin 1
PWY-7684	2.4.1.208	diacylglycerol &alpha;-glucosyltransferase
PWY-7684	2.4.1.337	&alpha;-monoglucosyldiacylglycerol synthase
PWY-7685	4.1.1.101	malolactic enzyme
PWY-7686	1.1.1.38 	NAD-dependent malic enzyme
PWY-7687	1.14.11	stipitaldehyde synthase
PWY-7687	1.14.13	3-methylorcinaldehyde monooxygenase
PWY-7687	1.14.13	stipitaldehyde monooxygenase
PWY-7687	2.3.1	3-methylorcinaldehyde synthase
PWY-7687	4.1.1.60	stipitatonate decarboxylase
PWY-7688	2.1.1.236	dTDP-3-amino-3,6-dideoxy-&alpha;-<small>D</small>-galactopyranose N,N-dimethyltransferase
PWY-7688	2.3.1	ravidosamine 4-O-acetyltransferase
PWY-7688	2.6.1.90	dTDP-3-amino-3,6-dideoxy-&alpha;-<small>D</small>-galactopyranose transaminase
PWY-7688	2.7.7.24	glucose-1-phosphate thymidylyltransferase
PWY-7688	4.2.1.46	dTDP-glucose 4,6-dehydratase
PWY-7688	5.3.2.3	dTDP-6-deoxy-3,4-keto-hexulose isomerase
PWY-7689	1.1.1 	methylfusarubin monooxygenase
PWY-7689	2.1.1	fusarubinaldehyde O-methyltransferase
PWY-7689	2.3.1	fusarubinaldehyde synthase
PWY-7690	1.4.99	(5E)-3-amino-4-sulfanyl-5-(sulfanylmethylidene)pyrrolidin-2-one dehydrogenase
PWY-7690	1.8.3.M1	holomycin dithiol oxidase
PWY-7690	2.3.1	MONOMER-19320
PWY-7690	4.1.1	MONOMER-19405
PWY-7690	6.3.2	2-azanidyl-3-sulfanylprop-2-enamido]-3-[NRPS] synthase
PWY-7690	6.3.2	HlmE non-ribosomal peptide synthase
PWY-7691	2.1.1	MONOMER-19407
PWY-7692	1.14.13	bikaverin monooxygenase
PWY-7692	2.1.1	bikaverin O-methyltransferase
PWY-7692	2.3.1	pre-bikaverin synthase
PWY-7693	1.1.1	L-seryl-[seryl-carrier protein] dehydrogenase
PWY-7693	1.14.11.M1	6-dehydroguadinomine B 6-hydroxylase
PWY-7693	1.2.98	2-aminomalonyl semialdehyde-[seryl-carrier protein] hydroxylase
PWY-7693	2.1.3	guadinomine carbamoyltransferase
PWY-7693	2.1.4.1	glycine amidinotransferase
PWY-7693	2.3.1	guadinomine polyketide synthase
PWY-7693	2.7.2	carbamate phosphorylase
PWY-7693	6.2.1	L-serine--[seryl-carrier protein] ligase
PWY-7694	1.1.1	D-glyceryl/L-seyl-[acp] dehydrogenase
PWY-7694	1.2.98	2-aminomalonyl semialdehyde-[pcp] hydroxylase
PWY-7694	1.2.98	2-hydroxy-3-oxopropanoyl-[acp] monooxygenase
PWY-7694	2.1.3	(S)-2,3-diaminopropanoate carbamoyltransferase
PWY-7694	2.3.1	1,3-bisphosphoglycerate:acyl-carrier protein glyceryltransferase
PWY-7694	2.3.1	proto-zwittermicin A-L-alanine C-hydroxylase
PWY-7694	2.3.1	zwittermicin A polyketide synthase complex
PWY-7694	2.6.1	L-serine aminase
PWY-7694	3.4.11	proto-zwittermicin A hydrolase
PWY-7694	6.2.1	L-serine--[seryl-carrier protein] ligase
PWY-7695	1.1.1	dimeric rubrofusarin dehydrogenase
PWY-7695	1.14.13	rubrofusarin monooxygenase
PWY-7695	1.14.21 	rubrofusarin laccase
PWY-7695	2.1.1	nor-rubrofusarin O-methyltransferase
PWY-7695	2.3.1	YWA1 synthase
PWY-7695	4.2.1	YWA1 dehydratase
PWY-7696	2.3.1	bikisocoumarin synthase
PWY-7696	2.3.1	YWA1 synthase
PWY-7697	2.4.1	geraniol glucoside 6-O-xylosyltransferase
PWY-7697	2.4.1	geraniol glucosyltransferase
PWY-7698 	1.1.1.90	benzyl alcohol dehydrogenase
PWY-7698 	1.13.11.4	gentisate 1,2-dioxygenase I
PWY-7698 	1.13.11.4	gentisate 1,2-dioxygenase II
PWY-7698 	1.14.13.24	3-hydroxybenzoate 6-hydroxylase I
PWY-7698 	1.14.13	xylenol methylhydroxylase
PWY-7698 	1.2.1.7	MONOMER-2948
PWY-7698 	3.7.1.23	maleylpyruvate hydrolase
PWY-7701	4.1.3.17	4-hydroxy-4-methyl-2-oxoglutarate aldolase
PWY-7702	2.3.1	enoyl reductase
PWY-7702	6.3.2 	sch210971/2 synthetase
PWY-7703	1.14.13	4-hydroxyisophthalate monooxygenase
PWY-7703	1.14.13.M11	4-hydroxy-3-methylbenzoate monooxygenase
PWY-7703 	1.2.1.96	4-hydroxybenzaldehyde dehydrogenase
PWY-7704 	1.14.11	4-methoxycyclopenine synthase
PWY-7704	6.3.2.40	cyclopeptine synthetase
PWY-7705 	6.3.2.40	cyclopeptine synthetase
PWY-7706	1.14.11.M2	N-3-fumaramoyl-(S)-2,3-diaminopropanoyl-[DdaD non-ribosomal peptide synthase] oxygenase
PWY-7706	2.3.1 	DdaD non-ribosomal peptide synthase
PWY-7706	3.1.2	MONOMER-19480
PWY-7706	6.3.2.46	fumarate-(S)-2,3-diaminopropanoate ligase
PWY-7706	6.3.2.47	dapdiamide synthase
PWY-7708	1.14.13	lyngbyatoxin A monooxygenase
PWY-7708	1.14.19	(-)-indolactam synthase
PWY-7708	2.3.1	N-methyl-L-valyl-L-tryptophanol synthase
PWY-7708	2.5.1	MONOMER-19482
PWY-7710	1.19.6.1	[FeMo]-nitrogenase complex reductase component
PWY-7710	2.3.3.14	homocitrate synthase
PWY-7712	3.1.1	trichothecene C-3 esterase
PWY-7713 	1.14.13	trichothecene 7,8-dihydroxylase
PWY-7715 	1.14.13	12,13-epoxytrichothec-9-ene 4-hydroxylase
PWY-7715 	1.14.13	isotrichodermin 15-hydroxylase
PWY-7715 	1.14.13	trichodiene oxygenase
PWY-7715 	1.14.13	trichothecene 4-hydroxylase
PWY-7715 	1.14.13	trichothecene 8-hydroxylase
PWY-7715 	2.3.1	isotrichodermol 3-acetyltransferase
PWY-7715 	2.3.1	trichodermol 15-O-acetyltransferase
PWY-7715 	2.3.1	trichodermol 15-O-acyltransferase
PWY-7715 	2.3.1	trichothecene 15-O-acetyltransferase
PWY-7715 	2.3.1	trichothecene C-4-acetyltransferase
PWY-7715 	2.3.1	trichothecene C-8-acyltransferase
PWY-7715 	3.1.1	trichothecene C-15 esterase
PWY-7715 	3.1.1	trichothecene C-3 esterase
PWY-7715 	4.2.3.6	trichodiene synthase
PWY-7716	2.3.1.164	isopenicillin N acyltransferase
PWY-7716 	2.3.1 	acyl coenzyme A: isopenicillin N acyltransferase
PWY-7716	6.2.1.30	phenylacetyl-CoA ligase
PWY-7717	1.13.11.11 	tryptophan 2,3-dioxygenase
PWY-7717	2.1.1.M8	3-hydroxy-D-kynurenine methyltransferase
PWY-7717	3.5.1.9	arylformamidase
PWY-7717	3.7.1	3-hydroxy-4-methyl-D-kynurenine hydrolase
PWY-7718	1.10.3.4	phenoxazinone synthase
PWY-7718	2.7.7.x	3-hydroxy-4-methylanthranilate adenylyltransferase
PWY-7718	6.3.1	3-hydroxy-4-methyl-anthranilate pentapeptide lactone synthase III
PWY-7718	6.3.1	actinomycin synthase II
PWY-7719	2.3.1	UDP-4-amino-4,6-dideoxy-N-acetyl-&beta;-<small>L</small>-idosamine acetyltransferase
PWY-7719	2.5.1	diacetamido-8-epilegionaminic acid synthase
PWY-7719	2.6.1	UDP-4-amino-4,6-dideoxy-N-acetyl-&beta;-<small>L</small>-idosamine transaminase
PWY-7719	2.7.7	diacetamido-8-epilegionaminic acid cytidylyl transferase
PWY-7719	3.6.1	UDP-2,4-diacetamido-2,4,6-trideoxy-&beta;-L-gulopyranose hydrolase
PWY-7719	4.2.1.115	UDP-N-acetylglucosamine 4,6-dehydratase (configuration-inverting)
PWY-7719	5.1.3	UDP-2,4-diacetamido-2,4,6-trideoxy-&beta;-L-idopyranose C2-epimerase
PWY-7720	4.2.3.145 	ophiobolin F synthase
PWY-7722	1.1.1.390	sulfoquinovose 1-dehydrogenase
PWY-7722	1.2.1.97	3-sulfolactaldehyde dehydrogenase
PWY-7722	3.1.1.99	6-deoxy-6-sulfogluconolactonase
PWY-7722	4.1.2.58	2-dehydro-3,6-dideoxy-6-sulfogluconate aldolase
PWY-7722	4.2.1.162	6-deoxy-D-gluconate 6-sulfonate dehydratase
PWY-7723	1.14.14.3	bacterial luciferase
PWY-7723	1.2.1.50	long-chain acyl-protein thioester reductase
PWY-7723	1.5.1.38	NADPH-flavin oxidoreductase
PWY-7723	1.5.1.42	FMN reductase (NADH)
PWY-7723	3.1.2.2 	activated long-chain acyl hydrolase
PWY-7723	6.2.1.19	long-chain-fatty-acid--protein ligase
PWY-7725 	1.14.19.44	acyl-CoA 5-desaturase
PWY-7727 	1.14.19 	acyl-CoA 4/6/8-desaturase
PWY-7729	2.1.1.M10	5-methoxy-benzimidazole C-methyltransferase
PWY-7729	2.1.1.M9	5-hydroxy-benzimidazole O-methyltransferase
PWY-7729	4.1.99.M1	5-hydroxybenzimidazole synthase
PWY-7731 	1.10.3.9	photosystem II
PWY-7731 	1.10.9.1	Cyt b6f complex
PWY-7731 	1.12.7.2	Fe hydrogenase
PWY-7731	1.6.5.2	NAD(P)H: plastoquinone dehydrogenase
PWY-7731 	1.97.1.12	photosystem I
PWY-7731 	1.97.1.12	Photosystem I iron-sulfur center
PWY-7733	1.13.11	&beta;-hydroxy-L-tryptophan 2,3-dioxygenase
PWY-7733	1.14.13	L-tryptophanyl-[tryptophanyl-carrier protein] &beta;-monooxygenase
PWY-7733	1.17.1	3-hydroxykynurenate reductase/dehydratase
PWY-7733	2.6.1	(2S,3R)-4-(2-aminophenyl)-2-amino-3-hydroxy-4-oxobutanoate aminotransferase
PWY-7733	3.1.2	&beta;-hydroxy-L-tryptophanyl-[tryptophanyl-carrier protein] thioesterase
PWY-7733	6.2.1	tryptophanyl-carrier-protein Swb11
PWY-7734	1.1.1	(2S,3R)-2-[(2-aminophenyl)amino]-3-carboxy-3-hydroxypropanoate dehydrogenase
PWY-7734	1.13.11	&beta;-hydroxy-L-tryptophan 2,3-dioxygenase
PWY-7734	1.14.13	L-tryptophanyl-[tryptophanyl-carrier protein] &beta;-monooxygenase
PWY-7734	2.6.99	(2S,3R)-4-(2-aminophenyl)-2-amino-3-hydroxy-4-oxobutanoate oxidoreductase
PWY-7734	3.1.2	&beta;-hydroxy-L-tryptophanyl-[tryptophanyl-carrier protein] thioesterase
PWY-7734	3.5.1.9	N-formyl-&beta;-hydroxy-L-kynurenine formamidase
PWY-7734	6.2.1	tryptophanyl-carrier-protein Ecm13
PWY-7735	1.8.3	triostin A dithiol oxidase
PWY-7735	2.1.1	echinomycin synthase
PWY-7735	2.3.1	MONOMER-19598
PWY-7735	2.7.7	quinoxaline-2-carboxylate activating enzyme
PWY-7735	6.3.2	TrsI non-ribosomal peptide synthase
PWY-7735	6.3.2	TrsJ non-ribosomal peptide synthase
PWY-7736	1.14.13	stellata-2,6,19-triene monooxygenase
PWY-7736	4.2.3 	stellata-2,6,19-triene synthase
PWY-7737	2.3.1	3-hydroxy-quinaldate carrier protein
PWY-7737	2.7.7	3-hydroxyquinaldate adenylase
PWY-7737	6.3.2	TioR non-ribosomal peptide synthase
PWY-7737	6.3.2	TioS non-ribosomal peptide synthase
PWY-7738	2.3.1.252	mycolipanoate synthase
PWY-7738	2.3.1.M11	Chp2 diacyltrehalose acyltransferase
PWY-7738	2.3.1.M12 	2,3-diacyltrehalose acyltransferase
PWY-7739	1.14.13	3-hydroxy-5-methoxybiphenyl 4-hydroxylase
PWY-7739	2.1.1	3,5-dihydroxybiphenyl O-methyltransferase
PWY-7739	2.1.1	noraucuparin O-methyltransferase
PWY-7739	2.3.1.177	3,5-dihydroxybiphenyl synthase
PWY-7740	2.3.1.M3	mycoketide-CoA synthase
PWY-7741	1.1.1.M10	(phenol)phthiodiol-4-en-one enoyl reductase
PWY-7741	1.1.1.M11	phthiodiolone ketoreductase
PWY-7741	2.1.1.M12	(phenol)phthiotriol methyltransferase
PWY-7741	2.7.7.M6	long-chain fatty acid adenylase
PWY-7741	3.1.2.M1	polyketide thiesterase TesA
PWY-7742 	1.1.1.M10	(phenol)phthiodiolenone enoyl reductase
PWY-7742 	1.1.1.M11	(phenol)phthiodiolone ketoreductase
PWY-7742 	2.1.1.M12	(phenol)phthiotriol methyltransferase
PWY-7742	2.3.1 	4-hydroxybenzate adenylase
PWY-7742	2.3.1.bp	4-hydroxyphenylalkanoate synthase
PWY-7742 	2.3.1.M4	(phenol)carboxyphthiodiolenone synthase
PWY-7742	2.7.7.y	4-hydroxyphenylalkanoate adenylase
PWY-7743 	2.3.1.M9	phenolphthiocerol/phthiocerol/phthiodiolone dimycocerosyl transferase
PWY-7744 	2.3.1.111	mycocerosic acid synthase
PWY-7744 	2.3.1.M9	phenolphthiocerol/phthiocerol/phthiodiolone dimycocerosyl transferase
PWY-7744 	2.7.7.aa	mycocerosic acid adenylyltransferase
PWY-7744 	3.1.2.M1	polyketide synthase thioesterase
PWY-7745 	2.1.1.M13	4-O-(&alpha;-L-rhamnopyranosyl)-hydroxybenzoate methyl ester 2-O-methyltransferase
PWY-7745 	2.1.1.M14	[&alpha;-L-fucopyranosyl-(1&rarr;3)-&alpha;-L-rhamnopyranosyl-(1&rarr;3)-2-O-methyl-&alpha;-L-rhamnopyranosyl] dimycocerosyl phenol-phthiocerol 2-O-methyltransferase
PWY-7745 	2.1.1.M15	[2-O-methyl-&alpha;-L-fucopyranosyl-(1&rarr;3)-&alpha;-L-rhamnopyranosyl-(1&rarr;3)-2-O-methyl-&alpha;-L-rhamnopyranosyl] dimycocerosyl phenol-phthiocerol 4-O-methyltransferase
PWY-7745 	2.1.1.M16	[2,4-di-O-methyl-&alpha;-L-fucopyranosyl-(1&rarr;3)-&alpha;-L-rhamnopyranosyl-(1&rarr;3)-2-O-methyl-&alpha;-L-rhamnopyranosyl] dimycocerosyl phenol-phthiocerol 3-O-methyltransferase
PWY-7745 	2.4.1.M3	dimycocerosyl phenolphthiocerol rhamnosyltransferase
PWY-7745 	2.4.1.M4	mycoside B rhamnosyltransferase
PWY-7745 	2.4.1.M5	[&alpha;-L-rhamnopyranosyl-(1&rarr;3)-2-O-methyl-&alpha;-L-rhamnopyranosyl] dimycocerosyl phenol-phthiocerol fucosyltransferase
PWY-7745 	4.1.3.40	chorismate lyase
PWY-7745	4.1.3.40	chorismate lyase
PWY-7746	2.3.1.M13	2-O-sulfo trehalose long-chain-acyltransferase
PWY-7746	2.3.1.M14	2-n-acyl-2-O-sulfo-trehalose (hydroxy)phthioceranyltransferase
PWY-7746	2.3.1.M15	3-(hydroxy)phthioceranyl-2-palmitoyl(stearoyl)-2-O-sulfo-trehalose (hydroxy)phthioceranyltransferase
PWY-7746	2.3.1.M7	phthioceranic/hydroxyphthioceranic acid synthase
PWY-7746	2.7.7.M7	long-chain fatty acid adenylyltransferase
PWY-7746	2.8.2.37	trehalose 2-sulfotransferase
PWY-7747	1.13.11.39	biphenyl-2,3-diol 1,2-dioxygenase
PWY-7747	1.14.12.18	biphenyl 2,3-dioxygenase
PWY-7747	1.3.1.56	cis-2,3-dihydrobiphenyl-2,3-diol dehydrogenase
PWY-7748	1.14.13.213	bursehernin 5-hydroxylase
PWY-7748	2.1.1.323	pluviatolide 4-O-methyltransferase
PWY-7748	2.1.1.fu	5-demethyl-yatein O-methyltransferase
PWY-7749	1.14.11.50	(-)-deoxypodophyllotoxin synthase
PWY-7749	1.14.13.214	(-)-deoxypodophyllotoxin hydroxylase
PWY-7749	1.14.13.M16	(-)-deoxypodophyllotoxin oxidative demethylase
PWY-7750	1.2.5.3	carbon monoxide dehydrogenase
PWY-7751	2.1.1	demethyl-4-deoxygadusol methyltransferase
PWY-7751	4.2.3.154	desmethyl-4-deoxygadusol synthase
PWY-7751	6.3.4	4-deoxygadusol glycyltransferase
PWY-7752	2.1.1.M22 	gadusol synthase
PWY-7752	4.2.3.152	2-epi-5-epi-valiolone synthase
PWY-7754	1.1.1.395	3&alpha;-hydroxysteroyl-CoA 3-dehydrogenase
PWY-7754	1.1.1.395	3&alpha;-hydroxysteroyl-CoA 3-dehydrogenase 2
PWY-7754	1.3.1.M6	7&alpha;-hydroxy-3-oxochol-4-en-24-oyl-CoA dehydrogenase
PWY-7754	1.3.1.M7	bile acid &Delta;<SUP>6</SUP>-reductase
PWY-7754	1.3.1.M8	bile acid &Delta;<SUP>4</SUP>-reductase
PWY-7754	2.8.3.h 	bile acid coenzyme A transferase
PWY-7754	2.8.3.h	bile acid coenzyme A transferase
PWY-7754	2.8.3.h	bile acid coenzyme A transferse
PWY-7754	3.5.1.74	chenodeoxycholoyltaurine hydrolase
PWY-7754	4.2.1.106	bile acid 7&alpha;-dehydratase
PWY-7754	6.2.1.7	bile acid-coenzyme A ligase
PWY-7755	1.1.1.391	NAD-dependent bile acid 3&beta;-dehydrogenase
PWY-7755	1.1.1.52	NAD-dependent bile acid 3&alpha;-dehydrogenase
PWY-7756	1.1.1.392	MONOMER-19696
PWY-7756	1.1.1.393	3 &beta;-hydroxysteroid dehydrogenase
PWY-7756	1.1.1.393	3&beta;-hydroxysteroid dehydrogenase
PWY-7757	1.14.13.M17	bisphenol A hydroxylase
PWY-7758	2.1.1.fw	bacteriochlorophyllide d C12<sup>1</sup>-methyltransferase
PWY-7758	2.1.1.fx	bacteriochlorophyllide d C8<sup>2</sup>-methyltransferase
PWY-7758	2.5.1.bl	bacteriochlorophyll d synthase
PWY-7758	3.1.1.r	chlorophyllide a hydrolase
PWY-7758	4.2.1.bi	3-vinyl-bacteriochlorophyllide d  3<sup>1</sup>-hydratase
PWY-7759	1.3.7.15 	chlorophyllide a reductase
PWY-7759	2.1.1.fw	8-ethyl-12-methyl-3-vinyl-bacteriochlorophyllide d 12<sup>1</sup>-methyltransferase
PWY-7759	2.1.1.fx	bacteriochlorophyllide d 8<sup>2</sup>-methyltransferase
PWY-7759	2.1.1.fy	bacteriochlorophyllide d C20 methyltransferase
PWY-7759	2.5.1.bl	bacteriochlorophyll c synthase
PWY-7759	3.1.1.r	chlorophyllide a hydrolase
PWY-7759	4.2.1.bi	3-vinyl-bacteriochlorophyllide d  3<sup>1</sup>-hydratase
PWY-7760	1.17.98.c	bacteriochlorophyllide c 7<sup>1</sup>-hydroxylase
PWY-7760	2.1.1.fw	8-ethyl-12-methyl-3-vinyl-bacteriochlorophyllide d 12<sup>1</sup>-methyltransferase
PWY-7760	2.1.1.fx	bacteriochlorophyllide d 8<sup>2</sup>-methyltransferase
PWY-7760	2.1.1.fy	bacteriochlorophyllide d C20 methyltransferase
PWY-7760	2.5.1.bl	bacteriochlorophyll e synthase
PWY-7760	3.1.1.r	chlorophyllide a hydrolase
PWY-7760	4.2.1.bi	3-vinyl-bacteriochlorophyllide d  3<sup>1</sup>-hydratase
PWY-7761	3.5.1.42	NMN aminohydrolase
PWY-7762	1.1.1.396	bacteriochlorophyllide a dehydrogenase
PWY-7762	1.3.1.111	geranylgeranyl-bacteriochlorophyll b reductase
PWY-7762	1.3.7.14	3,8-divinyl-chlorophyllide a reductase
PWY-7762	2.5.1.bk	bacteriohlorophyll b synthase
PWY-7762	4.2.1.165	chlorophyllide a 3<sup>1</sup>-hydratase
PWY-7764	1.3.7.13	divinyl chlorophyllide a 8-vinyl-reductase
PWY-7765	1.14.13.9	kynurenine 3-monooxygenase
PWY-7765	2.1.1.M8 	3-hydroxy-L-kynurenine methyltransferase
PWY-7765	3.7.1.M2	3-hydroxy-4-methyl-L-kynurenine hydrolase
PWY-7766	1.3.3.15	coproporphyrinogen III oxidase (coproporphyrin-forming)
PWY-7766	1.3.99	coproheme III oxidative decarboxylase
PWY-7766	4.1.1.37	uroporphyrinogen III decarboxylase
PWY-7766	4.99.1.9	coproporphyrin ferrochelatase
PWY-7767	1.1.1.345	4-methyl-2-oxopentanoate reductase
PWY-7767	2.8.3.24	2-hydroxyisocaproate CoA transferase
PWY-7767	4.2.1.157	(R)-2-hydroxyisocaproyl-CoA dehydratase
PWY-7769	1.1.1.309	phophonoacetaldehyde reductase
PWY-7769	1.1.1	hydroxymethylphosphonate dehydrogenase
PWY-7769	1.13.11.72	hydroxyethylphosphonate dioxygenase
PWY-7769	1.2.1	phosphonoformaldehyde dehydrogenase
PWY-7769	2.1.1.326	N-acetyl-demethylphophinothricin P-methyltransferase
PWY-7769	2.3.1.183	demethyl-phosphinothricin N-acetyltransferase
PWY-7769	2.3.1	N-acetyl demethylphosphinothricinyl-transferase
PWY-7769	2.3.1	nonribosomal peptide synthetase PhsC
PWY-7769	2.3.3.b	2-phosphinomethylmalate synthase
PWY-7769	2.7.1	CMP-5-phosphonoformate--3-phosphoglycerate phosphonoformyl transferase
PWY-7769	2.7.7.ad	phosphonoformate--CTP cytidylyltransferase
PWY-7769	2.7.7 	nonribosomal peptide synthetase PhsA
PWY-7769	2.7.7 	nonribosomal peptide synthetase PhsB
PWY-7769	2.7.8.23	carboxyphosphonoenolpyruvate phosphonomutase
PWY-7769	3.1.1	N-acetyl-phophinothricyl-L-alanyl-L-leucine acetyl hydrolase
PWY-7769	3.1.2	phosphinothricin tripeptide thioesterase
PWY-7769	4.1.1.82	phosphonopyruvate decarboxylase
PWY-7769	4.2.1.166	phosphinomethylmalate isomerase
PWY-7769	4.2.1	2-phophono-formylglycerate enolase
PWY-7769	5.4.2.9	phosphoenolpyruvate phosphomutase
PWY-7770	1.1.1.397	&beta;-methylindole-3-pyruvate reductase
PWY-7770	1.3.1	Ind4
PWY-7770	2.1.1.47	indolepyruvate C-methyltransferase
PWY-7770	2.1.1.fz	N-demethylindolmycin N-methyltransferase
PWY-7770	2.6.1.99	tryptophan&mdash;pyruvate aminotransferase
PWY-7770	4.3.3 	(2R,4E)-2-amino-5-({[(2S,3R)-2-hydroxy-3-(1H-indol-3-yl)butanamido]methanimidoyl}amino)pent-4-enoate lyase
PWY-7770	4.3.3 	Ind5-Ind6 complex
PWY-7770	6.3.2	4,5-dehydro-L-arginine--indolmycenate ligase
PWY-7771	3.1.3	butachlor esterase
PWY-7772	1.13.11.37	hydroxyquinol 1,2-dioxygenase
PWY-7772	1.14.13.220	resorcinol 4-hydroxylase (NADH)
PWY-7772	1.3.1.32	maleylacetate reductase
PWY-7772	4.1.1.n	&gamma;-resorcylate decarboxylase
PWY-7773	1.13.11.37	hydroxyquinol 1,2-dioxygenase
PWY-7773	1.14.14.27	resorcinol 4-hydroxylase (FADH2)
PWY-7773	1.3.1.32	maleylacetate reductase
PWY-7773	4.1.1.n	&gamma;-resorcylate decarboxylase
PWY-7774	1.1.1.M22	NAD+-dependent secondary alcohol dehydrogenase I
PWY-7774	1.1.1.M22	NAD+-dependent secondary alcohol dehydrogenase II
PWY-7774	1.1.1.M22	NAD+-dependent secondary alcohol dehydrogenase III
PWY-7774	1.14.13.dm	acetone monooxygenase (methylacetate-forming)
PWY-7774	1.14.13.dn	propane 2-monooxygenase
PWY-7774	3.1.1.M3	methyl acetate hydrolase
PWY-7775	1.1.1.M22	MONOMER-19816
PWY-7775	1.14.13.M31 	propane 2-monooxygenase
PWY-7775	1.14.13.M34	acetol monooxygenase
PWY-7776	1.1.1.M24	2-hyroxyethyl-CoM 2-dehydrogenase
PWY-7776	1.14.13.69	alkene monooxygenase
PWY-7776	1.8.1.5	2-oxoehtyl-CoM reductase/carboxylase
PWY-7776	4.4.1.23	epoxyalkane:coenzyme M transferase
PWY-7777	1.14.13.69	alkene monooxygenase
PWY-7777	4.4.1.i	1,2-epoxy-2-methyl-3-butene--glutathione S-transferase
PWY-7778	1.14.13.69	alkene monooxygenase
PWY-7778	1.2.1.98 	2-methyl-propane-1,2-diol dehydrogenase
PWY-7778	3.3.2.10	1,2-epoxy-2-methylpropane hydrolase
PWY-7779	1.1.1.35	(S)-3-hydroxybutanoyl-CoA dehydrogeanse
PWY-7779	1.1.1.400	2-methyl-propane-1,2-diol dehydrogenase
PWY-7779	1.1.1.M26	tert-butoxymethanol dehydrogenase
PWY-7779	1.14.13.M36	MTBE monooxygenase
PWY-7779	1.14.19.ar 	tert-butyl alcohol hydroxylase
PWY-7779	1.2.1.98	2-hydroxy-2-methylpropanal dehydrogenase
PWY-7779	2.3.1.9	acetyl-CoA C-acetyltransferase
PWY-7779	3.7.1	tert-butyl formate deformylase
PWY-7779	5.3.3.20	2-hydroxyisobutanoyl-CoA mutase
PWY-7779	6.2.1.M4	2-hydroxyisobutanoate-CoA ligase
PWY-7780	1.1.2.9	1-butanol dehydrogenase (cytochrome c)
PWY-7780	1.14.13	butanal monooxygenase
PWY-7780	1.14.13.dq	soluble butane monooxygenase
PWY-7780	1.1.5.11	1-butanol dehydrogenase
PWY-7781	1.14.15.M3	II-dihydromenaquinone-9 &omega;-hydroxylase
PWY-7781	2.8.2	&omega;-hydroxy-II-dihydromenaquinone-9 sulfotransferase Stf3
PWY-7782	1.1.1.101	acyl/alkyl-dihydroxyacetone phosphate reductase
PWY-7782	1.14.19	2-acyl-1-alkyl-sn-glycero-3-phosphoethanolamine desaturase
PWY-7782	1.2.1.84	fatty acyl-CoA reductase 1
PWY-7782	2.3.1.42	glycerone-phosphate O-acyltransferase
PWY-7782 	2.3.1.51	1-acyl-sn-glycerol-3-phosphate acyltransferase alpha
PWY-7782	2.5.1.26	alkylglycerone-phosphate synthase
PWY-7782	2.7.1.82	ethanolamine kinase 1
PWY-7782	2.7.7.14	ethanolamine-phosphate cytidylyltransferase
PWY-7782	2.7.8.1	ethanolaminephosphotransferase 1
PWY-7782 	2.7.8.2 	choline/ethanolaminephosphotransferase 1
PWY-7782	3.1.3.4	2-acyl-1-alkyl-sn-glycerol 3-phosphate phosphohydrolase
PWY-7783	3.1.1.4	calcium-independent cytosolic phospholipase A2
PWY-7783	3.1.4.39	alkylglycerophosphoethanolamine phosphodiesterase
PWY-7783 	3.1.4.3	phosphatidylcholine-specific phospholipase C
PWY-7783	3.3.2.2	lysoplasmalogenase
PWY-7786	1.1.1.403	D-threitol dehydrogenase
PWY-7786 	2.7.1.210	D-erythrulose kinase
PWY-7787 	1.1.1.56 	erythritol/L-threitol dehydrogenase
PWY-7787	2.7.1.209	L-erythrulose 4-kinase
PWY-7787	5.1.3.39	L-erythrulose-4-phosphate epimerase
PWY-7787 	5.3.1.34	D-erythrulose-4-phosphate isomerase 1
PWY-7789	1.1.1.402	D-erythritol 1-phosphate dehydrogenase
PWY-7789	2.7.1.az	erythritol kinase (D-erythritol 1-phosphate-forming)
PWY-7789	5.1.3.38	D-erythrulose 1-phosphate 3-epimerase
PWY-7789	5.3.1.33	L-erythrulose 1-phosphate isomerase
PWY-7789	5.3.1.34	D-erythrulose 4-phosphate isomerase
PWY-7790	1.3.98.1	dihydroorotate dehydrogenase, type 1A
PWY-7790	2.1.3.2 	carbamoyl phosphate synthase/aspartate carbamoyltransferase
PWY-7790	2.4.2.10	orotate phosphoribosyltransferase
PWY-7790	3.5.2.3	dihydrooratase
PWY-7790	4.1.1.23	orotidine-5-phosphate decarboxylase
PWY-7791	1.3.1.14	dihydroorotate dehydrogenase
PWY-7791	1.3.1.14	dihydroorotate dehydrogenase, type 1B
PWY-7791	2.1.3.2	aspartate carbamoyltransferase
PWY-7791	2.4.2.10	orotate phosphoribosyltransferase
PWY-7791	3.5.2.3	dihydroorotase
PWY-7791	4.1.1.23	orotidine 5-phosphate decarboxylase
PWY-7791	6.3.5.5 	carbamoyl-phosphate synthetase P
PWY-7793	2.1.1.ga	methanethiol S-methyltransferase
PWY-7793	4.4.1.11	methionine &gamma;-lyase
PWY-7794	3.1.1.s	poly(ethylene terephthalate) hydrolase
PWY-7794	3.1.1.t	mono(ethylene terephthalate) hydrolase
PWY-7795	1.14.12.15	terephthalate 1,2-dioxygenase oxygenase component 
PWY-7795	1.3.1.53	1,2-dihydroxy-3,5-cyclohexadiene-1,4-dicarboxylate dehydrogenase
PWY-7796	4.1.2.43 	bifunctional 3-hexulose-6-phosphate formaldehyde lyase/6-phospho-3-hexuloisomerase
PWY-7797	1.14.15.M5 	nocardicin C N-oxygenase
PWY-7797	2.5.1.38	nocardicin 3-amino-3-carboxypropyltransferase
PWY-7797	3.4.21.4	trypsin
PWY-7797	5.1.1.14	nocardicin C-9-epimerase
PWY-7797	6.3.2	non-ribosomal peptide synthetase NocAB
PWY-7798 	1.1.1.284 	alcohol dehydrogenase class-3
PWY-7798 	1.1.1.284 	S-hydroxymethyl-glutathione dehydrogenase/S-nitrosoglutathione dehydrogenase
PWY-7799	2.3.2.8	arginyl-tRNA--protein transferase 1
PWY-7799 	3.4.11.18	methionine aminopeptidase 1
PWY-7799 	3.4.11.18	methionine aminopeptidase 2
PWY-7799	3.4.19.1	acylaminoacyl-peptidase
PWY-7799	3.5.1.aa	protein N-terminal asparagine amidohydrolase
PWY-7799	3.5.1.ab	protein N-terminal glutamine amidohydrolase
PWY-7800 	2.3.1.bt	N-terminal methionine N<sup>&alpha;</sup>-acetyltransferase NatB
PWY-7800	2.3.1.bu	N-terminal amino acyl N<sup>&alpha;</sup>-acetyltransferase NatA
PWY-7800	2.3.1.bx 	N-terminal methionine N<sup>&alpha;</sup>-acetyltransferase NatC
PWY-7800	2.3.1.bx 	N-terminal methionine N<sup>&alpha;</sup>-acetyltransferase NatE
PWY-7801	2.3.2.6	leucyl, phenylalanyl-tRNA-protein transferase
PWY-7802	2.3.2.p	leucyl<sup>D,E</sup>-transferase
PWY-7803	4.6.1.16	tRNA-splicing endonuclease
PWY-7803	6.5.1.d	tRNA-splicing ligase complex
PWY-7804	4.3.99.c	glyphosate oxidoreductase
PWY-7805	2.3.1.M28	aminoalkylphosphonate N-acetyltransferase
PWY-7805 	2.7.8.37	methylphosphonate degradation complex
PWY-7805	3.6.1.1 	inorganic pyrophosphatase
PWY-7805 	3.6.1.63	RPnTP hydrolase
PWY-7805 	4.7.1.1	carbon-phosphorus lyase core complex, PhnJ subunit
PWY-7806 	1.4.3.19	glycine oxidase
PWY-7807	2.4.2.7	adenine phosphoribosyltransferase
PWY-7807	2.7.4.23	ribose 1,5-bisphosphate phosphokinase
PWY-7807	2.7.8.37	&alpha;-D-ribose 1-methylphosphonate 5-triphosphate synthase subunit PhnI
PWY-7807	3.1.4.55	phosphoribosyl 1,2-cyclic phosphate phosphodiesterase
PWY-7807	3.6.1.1	inorganic diphosphatase
PWY-7807	3.6.1.63	&alpha;-<small>D</small>-ribose 1-methylphosphonate 5-triphosphate diphosphatase
PWY-7807	4.7.1.1	&alpha;-<small>D</small>-ribose 1-methylphosphonate 5-phosphate C-P-lyase
PWY-7808	1.14.13.dr	tetracycline 11a-monooxygenase
PWY-7809 	1.1.1.M31	nonaketamide C9-dehydrogenase
PWY-7809 	1.14.13.ds 	6-methylpretetramide 4,12a-monooxygenase
PWY-7809 	1.14.13.ds	6-methylpretetramide 4-monooxygenase
PWY-7809 	1.14.13.du 	anhydrotetracycline 5,6-monooxygenase
PWY-7809 	1.3.98.d	dehydrooxytetracycline dehydrogenase
PWY-7809 	2.1.1.gb	4-amino-anhydrotetracycline N<sup>&alpha;</sup>-methyltransferase
PWY-7809 	2.1.1.M31	pretetramide C-6-methyltransferase
PWY-7809 	2.3.1.bz	tetracycline polyketide synthase
PWY-7809 	2.3.1.M29	OxyC acyl-binding protein
PWY-7809 	2.6.1.M2	4-dedimethylamino-4-oxo-anhydrotetracycline
PWY-7809 	4.2.1.M10	tetracycline nonaketamide B/C rings cyclase/dehydratase
PWY-7809 	4.2.1.M9	tetracycline nonaketamide D-ring cylase/dehydratase
PWY-7809 	5.5.1	pretetramide synthase
PWY-7809 	6.3.5.M2	malonamoyl-CoA synthetase
PWY-7810	1.14.19.as	tetracycline 7-halogenase
PWY-7814	1.1.1	dTDP-3-amino-4-dehydro-2,3,6-trideoxy-&beta;-L-glucose 4-ketoreductase
PWY-7814	2.6.1	dTDP-3,4-didehydro-2,6-dideoxy-&beta;-L-glucose 3-aminotranferase
PWY-7814	2.7.7.24	glucose-1-phosphate thymidylyltransferase
PWY-7814	4.2.1.46	dTDP-D-glucose 4,6-dehydratase
PWY-7814	4.2.1	dTDP-4-dehydro-&beta;-L-rhamnose 2,3-dehydratase
PWY-7814	5.1.3.13	dTDP-4-dehydrorhamnose 3,5-epimerase
PWY-7815	1.1.1.137	ribitol-5-phosphate 1-dehydrogenase
PWY-7815	2.3.1.M31	D-alanine carrier protein:PG D-alanyltransferase
PWY-7815	2.3.1.M32	PG:teichoic acid D-alanyltransferase
PWY-7815	2.4.1.187	N-acetylglucosaminyldiphosphoundecaprenol N-acetyl-&beta;-<small>D</small>-mannosaminyltransferase
PWY-7815	2.4.1.53	(poly)ribitol-phosphate teichoic acid &beta;-D-glucosyltransferase
PWY-7815	2.7.7.39	glycerol-3-phosphate cytidylyltransferase
PWY-7815	2.7.7.40	ribitol-5-phosphate cytidylyltransferase
PWY-7815	2.7.8.14	teichoic acid poly(ribitol phosphate) polymerase
PWY-7815	2.7.8.33	undecaprenyl-phosphate N-acetylglucosaminyl 1-phosphate transferase
PWY-7815	2.7.8.M1	polyisoprenyl-teichoic acid--peptidoglycan teichoic acid transferase TarT
PWY-7815	2.7.8.M1	polyisoprenyl-teichoic acid--peptidoglycan teichoic acid transferase TarU
PWY-7815	2.7.8.M1	polyisoprenyl-teichoic acid--peptidoglycan teichoic acid transferase TarV
PWY-7815	2.7.8.M2	teichoic acid ribitol-phosphate primase
PWY-7815	2.7.8.q	teichoic acid glycerol-phosphate primase
PWY-7815	5.1.3.14	UDP-N-acetylglucosamine 2-epimerase
PWY-7815	6.2.1.M5	D-alanine--[D-alanyl carrier protein] ligase
PWY-7816	1.1.1.137	ribitol-5-phosphate 1-dehydrogenase
PWY-7816 	2.3.1.M31	D-alanine carrier protein:PG D-alanyltransferase
PWY-7816	2.4.1.187	N-acetylglucosaminyldiphosphoundecaprenol N-acetyl-&beta;-<small>D</small>-mannosaminyltransferase
PWY-7816	2.4.1.70	poly(ribitol phosphate) teichoic acid &alpha;-O-GlcNAc transferase
PWY-7816	2.4.1.M14	poly(ribitol phosphate) teichoic acid &beta;-O-GlcNAc transferase
PWY-7816	2.7.7.39	glycerol-3-phosphate cytidylyltransferase
PWY-7816	2.7.7.40	ribitol-5-phosphate cytidylyltransferase
PWY-7816	2.7.8.14 	teichoic acid ribitol-phosphate primase/polymerase
PWY-7816	2.7.8.33	UDP-N-acetylglucosamine&mdash;undecaprenyl-phosphate N-acetylglucosaminephosphotransferase
PWY-7816	2.7.8.M1	polyisoprenyl-teichoic acid--peptidoglycan teichoic acid transferase MsrR
PWY-7816	2.7.8.M1	polyisoprenyl-teichoic acid--peptidoglycan teichoic acid transferase SA2103
PWY-7816	2.7.8.M1	polyisoprenyl-teichoic acid--peptidoglycan teichoic acid transferase SA908
PWY-7816	2.7.8.q	teichoic acid glycerol-phosphate primase
PWY-7816	2.7.8	teichoic acid glycerol-phosphate transferase
PWY-7816 	6.2.1.M5	D-alanine--[D-alanyl carrier protein] ligase
PWY-7817 	2.3.1.M32	PG:teichoic acid D-alanyltransferase
PWY-7817	2.4.1.315 	processive diacylglycerol &beta;-glucosyltransferase
PWY-7817	2.4.2	undecaprenyl phosphate-N-acetyl-&alpha;-D-glucosamine transferase
PWY-7817	2.7.1.107	diacylglycerol kinase
PWY-7817	2.7.7.64 	UTP--glucose-1-phosphate uridylyltransferase
PWY-7817	2.7.8.M4 	poly(glycerol phosphate) lipoteichoic acid synthase
PWY-7818	1.1.1.137	ribitol-5-phosphate 1-dehydrogenase
PWY-7818	2.3.1.M31	D-alanine carrier protein:PG D-alanyltransferase
PWY-7818	2.3.1.M32	PG:teichoic acid D-alanyltransferase
PWY-7818	2.4.1.337	1,2-diacylglycerol 3-&alpha;-glucosyltransferase
PWY-7818	2.4.1	AATGal-PP-lipid glucosyltransferase
PWY-7818	2.4.1	GalNAc-Rbo-P-Glc-AATGal-PP-undecaprenol N-acetyl-D-galactosamine transferase
PWY-7818	2.4.1	Rbo-P-Glc-AATGal-PP-undecaprenol N-acetyl-D-galactosamine transferase
PWY-7818	2.7.1.32	choline kinase
PWY-7818	2.7.7.15	choline-phosphate cytidylyltransferase
PWY-7818	2.7.7.40	ribitol-5-phosphate cytidylyltransferase
PWY-7818	2.7.8	(Chol-P-GalNAc)-GalNAc-Rbo-P-Glc-AATGal-PP-undecaprenol phosphocholine transferase
PWY-7818	2.7.8	GalNAc-GalNAc-Rbo-P-Glc-AATGal-PP-undecaprenol phosphocholine transferase
PWY-7818	2.7.8	Glc-AATGal-PP-undecaprenol phosphoribityltransferase
PWY-7818	2.7.8.M1	polyisoprenyl-teichoic acid--peptidoglycan teichoic acid transferase Cps2A
PWY-7818	2.7.8.M1	polyisoprenyl-teichoic acid--peptidoglycan teichoic acid transferase LytR
PWY-7818	2.7.8.M1	polyisoprenyl-teichoic acid--peptidoglycan teichoic acid transferase Psr
PWY-7818	2.7.8	type IV lipoteichoic acid polymerase
PWY-7818	6.2.1.M5	D-alanine--[D-alanyl carrier protein] ligase
PWY-7819	2.4.1	poly(glucosyl N-acetylgalactosamine 1-phosphate) glucosyltransferase
PWY-7819	2.4.1	poly(glucosyl N-acetylgalactosamine 1-phosphate) N-acetylgalactosamine 1-phosphate primase/transferase
PWY-7820	1.1.1.22	UDP-glucose 6-dehydrogenase
PWY-7820	2.4.1	UDP-GalNAc:&alpha;-<small>D</small>-GalNAc-diphosphoundecaprenol &alpha;-1,6-N-acetylgalactosaminyltransferase
PWY-7820	2.4.1	UDP-GlcA:&alpha;-<small>D</small>-GalNAc-(1,6)-&alpha;-<small>D</small>-GalNAc-diphosphoundecaprenol &beta;-1,3-glucuronosyltransferase
PWY-7820	2.4.1	UDP-GlcA:&beta;-D-GlcA-(1,3)-&alpha;-<small>D</small>-GalNAc-(1,6)-&alpha;-<small>D</small>-GalNAc-diphosphoundecaprenol &beta;-1,4-glucuronosyltransferase
PWY-7820	2.7.8.40	undecaprenyl-phosphate N-acetylgalactosaminyl 1-phosphate transferase
PWY-7820 	2.7.8.M1	polyisoprenyl-teichoic acid--peptidoglycan teichoic acid transferase TagT
PWY-7820 	2.7.8.M1	polyisoprenyl-teichoic acid--peptidoglycan teichoic acid transferase TagU
PWY-7820 	2.7.8.M1	polyisoprenyl-teichoic acid--peptidoglycan teichoic acid transferase TagV
PWY-7820	2.7.8	teichuronic acid polymerase
PWY-7821	2.3.1.M33	N-acetyl-D-glucosaminyl-tunicaminyl-uracil acyltransferase
PWY-7821	2.4.1.M16	N-acetyl-tunicaminyl-uracil N-acetylglucosmainyltransferase
PWY-7821	3.1.3.5	UMP phosphatase
PWY-7821	3.5.1.M4	N-acetyl-D-glucosaminyl-N-acetyl-tunicaminyl-uracil N-deacetylase
PWY-7821	3.6.1.9	nucleoside triphosphate pyrophosphatase
PWY-7821	3.6.1.M1	UDP-N-acetyl-tunicaminyl-uracil hydrolase
PWY-7821	4.1.99.M2	UDP-N-acetyl-tunicamine-uracil synthase
PWY-7821	4.2.1.M11	UDP-N-acetyl-&alpha;-D-glucosamine 5,6-dehydratase
PWY-7821	5.1.3.7 	UDP-6-deoxy-5,6-ene-GlcNAc 4-epimerase
PWY-7822	1.14.99.q	CBP21
PWY-7822	3.2.1.ak	chitinase B
PWY-7822	3.2.1.al	chitinase A
PWY-7823	3.5.2.M1	chlorzoxazone hydrolase
PWY-7824	1.14.13.124	L-phenylalanine N-monooxygenase
PWY-7824	1.14.13.68 	(Z)-phenylacetaldehyde oxime monooxygenase
PWY-7825	1.1.1.M32	3,5-dihydroxy-1,4-naphthoquinone reductase
PWY-7825	1.17.3.a	juglone hydroxylase 1
PWY-7825	1.17.3.a	juglone hydroxylase 2
PWY-7826	1.14.13	MONOMER-20067
PWY-7826	1.14.19.at 	noroxomaritidine synthase
PWY-7826	1.3.1 	noroxomaritidine reductase
PWY-7826	2.1.1.gc	norbelladine O-methyltranferase
PWY-7826	2.1.1	MONOMER-20068
PWY-7828	3.2.1.ao	1,3-&alpha;-isomaltosidase
PWY-7828	3.2.1.ap	isomaltose glucohydrolase
PWY-821 	1.8.1.2	assimilatory sulfite reductase (NADPH)
PWY-821 	1.8.4.8	Met16
PWY-821 	2.3.1.31	Met2
PWY-821 	2.5.1.49	Met17
PWY-821 	2.7.1.25	Met14
PWY-821 	2.7.7.4	Met3
PWY-821 	4.2.1.22	Cys4
PWY-822	2.4.1.10	sucrose:fructan 6-fructosyltransferase
PWY-822	2.4.1.99	sucrose:sucrose 1-fructosyltransferase
PWY-82	2.7.1.46	AT4G16130-MONOMER
PWY-82	2.7.1.46	MONOMER-452
PWY-82 	2.7.7.37 	MONOMER-11146
PWY-83	2.4.1	UDPG:glucosyltransferase
PWY-841 	2.1.2.2	YDR408C-MONOMER
PWY-841 	2.4.2.14	amidophosphoribosyltransferase
PWY-841 	2.4.2.14	phosphoribosylpyrophosphate amidotransferase
PWY-841 	2.7.4.13 	adenylate kinase 1
PWY-841 	3.5.4.10 	bifunctional purine biosynthesis protein ADE16
PWY-841 	3.5.4.10 	bifunctional purine biosynthesis protein ADE17
PWY-841 	3.5.4.10 	bifunctional purine biosynthesis protein PURH
PWY-841 	4.1.1.21 	multifunctional protein ADE2
PWY-841 	4.1.1.21	phosphoribosylaminoimidazole carboxylase
PWY-841 	6.3.2.6	YAR015W-MONOMER
PWY-841 	6.3.3.1 	bifunctional purine biosynthetic protein ADE5,7
PWY-841 	6.3.4.13 	trifunctional purine biosynthetic protein adenosine-3
PWY-841 	6.3.4.4	Ade12
PWY-841 	6.3.4.4	adenylosuccinate synthase
PWY-841 	6.3.5.3	phosphoribosylformylglycinamidine synthase
PWY-841 	6.3.5.3	YGR061C-MONOMER
PWY-84	2.3.1.95	trihydroxystilbene synthase
PWY-84	2.3.1.95	trihydroxystilbene synthase I
PWY-84	2.3.1	MONOMER-11632
PWY-84	2.3.1	MONOMER-11694
PWY-842	3.2.1.1 	&alpha;-amylase
PWY-842	3.2.1.20	&alpha;-glucosidase
PWY-842	3.2.1.2	&beta;-amylase
PWY-842	3.2.1.41	debranching enzyme
PWY-842	3.2.1 	&alpha;-amylase
PWY-861	1.14.13.41	tyrosine N-monooxygenase
PWY-861	1.14.13.68	4-hydroxyphenylacetaldehyde oxime monooxygenase
PWY-861	2.4.1.85	MONOMER-523
PWY-862	3.2.1	MONOMER-16316
PWY-862	3.2.1	MONOMER-2243
PWY-881	2.4.1	MONOMER-553
PWY-881	3.1.3.12	MONOMER-554
PWY-881 	3.1.3.12	MONOMER-6002
PWY-882	1.1.1.316	L-galactose dehydrogenase
PWY-882	1.1.1.316	MONOMER-2306
PWY-882	1.3.2.3	MONOMER-2307
PWY-882 	2.7.7.13	AT2G39770-MONOMER
PWY-882	2.7.7.69	GDP-L-galactose phosphorylase
PWY-882 	5.1.3.18 	MONOMER-11941
PWY-882 	5.3.1.8	MONOMER-2302
PWY-882 	5.4.2.8	MONOMER-2303
PWY-882 	5.4.2.8	phosphomannomutase
PWY8J2-20	1.14.15.13	pulcherriminic acid synthase
PWY8J2-20	2.3.2.22	cyclodipeptide synthase
PWY8J2-22	1.1.1.361	glucose-6-phosphate 3-dehydrogenase
PWY8J2-22	2.6.1.104	3-oxo-glucose-6-phosphate:glutamate aminotransferase
PWY8J2-22	3.1.3.92	kanosamine-6-phosphate phosphatase
PWY-922	2.7.1.36	MONOMER-11969
PWY-922	2.7.4.2	MONOMER-11968
PWYDQC-4 	1.14.13.168	indole-3-pyruvate monooxygenase
PWYDQC-4	1.14.13.168	indole-3-pyruvate monooxygenase
PWYDQC-4	2.6.1.99	L-tryptophan:pyruvate aminotransferase
PWY-I9	2.5.1.48	cystathionine &gamma;-lyase
PWY-I9 	2.5.1.6	S-adenosylmethionine synthase
PWY-I9	3.2.2.16 	5-methylthioadenosine/adenosylhomocysteine/6-amino-6-deoxyfutalosine nucleosidase
PWY-I9	4.2.1.22	cystathionine &beta;-synthase
PWY-I9	4.4.1.21	S-ribosylhomocysteine lyase
PWYQT-4427	2.4.1.M2	sulfoquinovosyldiacylglycerol synthase
PWYQT-4427	2.4.1.M2	UDP-sulfoquinovose:DAG sulfoquinovosyltransferase
PWYQT-4427	3.13.1.1	UDP-sulfoquinovose synthase
PWYQT-4429 	4.2.1.1	AT3G01500-MONOMER
PWYQT-4450	1.1	isopropylmalate dehydrogenase
PWYQT-4450 	2.3.3	methylthioalkylmalate synthase
PWYQT-4450	5.4.4	isopropylmalate isomerase, large subunit
PWYQT-4450	5.4.4	isopropylmalate isomerase, small subunit
PWYQT-4471 	2.6.1	branched-chain aminotransferase
PWYQT-4472 	1.14.13	CYP79F1 monooxygenase
PWYQT-4472 	2.4.1.195	UGT74B1
PWYQT-4472 	2.8.2 	desulfoglucosinolate sulfotransferase
PWYQT-4472 	2.8.2	desulfoglucosinolate sulfotransferase
PWYQT-4475 	1.14.13	CYP79F2 monooxygenase
PWYQT-4475 	1.14.99	glucosinolate S-oxygenase
PWYQT-4475	1.14.99	glucosinolate S-oxygenase
PWYQT-4475 	1	CYP83A1 monooxygenase
PWYQT-4476 	4.2.1.84 	Nit1
PWYQT-4476 	4.2.1.84 	Nit2
PWYQT-4476 	4.2.1.84 	Nit3
PYRIDNUCSAL-PWY	3.5.1.19 	NICOTINAMID-MONOMER
PYRIDNUCSAL-PWY 	3.6.1.22	MONOMER-8361
PYRIDNUCSAL-PWY	6.3.4.21	NICOTINATEPRIBOSYLTRANS-MONOMER
PYRIDNUCSAL-PWY 	6.3.5.1 	NAD synthetase, NH3-dependent
PYRIDNUCSYN-PWY	1.4.3.16	AT5G14760-MONOMER
PYRIDNUCSYN-PWY	1.4.3.16	MONOMER-12626
PYRIDNUCSYN-PWY	2.4.2.19	AT2G01350-MONOMER
PYRIDNUCSYN-PWY	2.4.2.19	MONOMER-12623
PYRIDNUCSYN-PWY	2.4.2.19	QPT subunit
PYRIDNUCSYN-PWY	2.5.1.72	AT5G50210-MONOMER
PYRIDNUCSYN-PWY 	2.7.7.18	MONOMER-12632
PYRIDNUCSYN-PWY 	6.3.5.1	MONOMER-12633
PYRIDOXSYN-PWY 	1.4.3.5	pyridoxine 5-phosphate oxidase / pyridoxamine 5-phosphate oxidase
PYRUVDEHYD-PWY	1.2.1 	pyruvate dehydrogenase
PYRUVDEHYD-PWY	1.2.1 	pyruvate dehydrogenase E1 component (somatic)
PYRUVDEHYD-PWY	1.2.1 	pyruvate dehydrogenase E2 component
PYRUVDEHYD-PWY	1.2.1 	pyruvate dehydrogenase, E2 subunit
PYRUVDEHYD-PWY 	1.8.1.4 	dihydrolipoyl dehydrogenase
PYRUVDEHYD-PWY 	1.8.1.4 	lipoamide dehydrogenase 
PYRUVDEHYD-PWY	1.8.1.4 	pyruvate dehydrogenase complex
QUINATEDEG-PWY	1.1.5.8	MONOMER-67
QUINATEDEG-PWY	4.2.1.10	periplasmic dehydroquinate dehydratase
QUINATEDEG-PWY 	4.2.1.118	MONOMER-1
REDCITCYC	1.1.1.42	isocitrate dehydrogenase subunit
REDCITCYC	1.1.5.4	malate:quinone oxidoreductase
REDCITCYC	1.2.7.3	2-oxoglutarate synthase
REDCITCYC	1.3.5.4	fumarate reductase
REDCITCYC	2.3.3.16 	citrate synthase (gltA)
REDCITCYC	2.8.3.5	succinyl-CoA:acetoacetate CoA-transferase subunit B 
REDCITCYC	4.2.1.2	fumarase (fumC)
RHAMCAT-PWY	2.7.1.5	rhamnulokinase
RHAMCAT-PWY 	4.1.2.19 	bifunctional L-rhamnulose 1-phosphate aldolase/lactaldehyde dehydrogenase
RHAMCAT-PWY	5.1.3.32	L-rhamnose mutarotase
RIBITOLUTIL-PWY	1.1.1.56	ribitol dehydrogenase subunit
RIBOKIN-PWY	2.7.1.15	ribokinase
RIBOKIN-PWY	5.4.99.62	ribose pyranase
RIBOSYN2-PWY	1.1.1.193	MONOMER-17844
RIBOSYN2-PWY	2.5.1.78	6,7-dimethyl-8-ribityllumazine synthase
RIBOSYN2-PWY	2.5.1.78	chloroplastic 6,7-dimethyl-8-ribityllumazine synthase
RIBOSYN2-PWY	2.5.1.9	riboflavin synthase
RIBOSYN2-PWY	2.7.1.26 	bifunctional riboflavin kinase / FMN adenylyltransferase
RIBOSYN2-PWY	2.7.7.2 	bifunctional riboflavin kinase / FMN adenylyltransferase
RIBOSYN2-PWY	2.7.7.2 	MONOMER-14605
RIBOSYN2-PWY	3.1.3.z 	5-amino-6-(5-phospho-D-ribitylamino)uracil phosphatase
RIBOSYN2-PWY	3.5.4.25 	chloroplastic riboflavin biosynthesis protein ribBA
RIBOSYN2-PWY	3.5.4.25	GTP cyclohydrolase II
RIBOSYN2-PWY	3.5.4.26	diaminohydroxyphosphoribosylaminopyrimidine deaminase
RIBOSYN2-PWY	3.5.4.26 	fused diaminohydroxyphosphoribosylaminopyrimidine deaminase / 5-amino-6-(5-phosphoribosylamino)uracil reductase
RIBOSYN2-PWY	4.1.99.12	3,4-dihydroxy-2-butanone 4-phosphate synthase
RUMP-PWY	1.1.1.343	6-phophogluconate dehydrogenase subunit
RUMP-PWY	1.1.1.49	glucose-6-phosphate dehydrogenase subunit
SALVADEHYPOX-PWY 	3.6.1.45 	USHA-MONOMER
SALVPURINE2-PWY 	2.4.2.1 	xanthosine phosphorylase
SALVPURINE2-PWY 	2.4.2.22 	xanthine-guanine phosphoribosyltransferase
SER-GLYSYN-PWY 	1.1.1.95	D-3-phosphoglycerate dehydrogenase
SER-GLYSYN-PWY 	2.1.2.1 	serine hydroxymethyltransferase
SER-GLYSYN-PWY 	2.1.2.1	serine hydroxymethyltransferase 1
SER-GLYSYN-PWY 	2.1.2.1	serine hydroxymethyltransferase 2
SER-GLYSYN-PWY 	2.6.1.52	phosphoserine aminotransferase
SER-GLYSYN-PWY 	3.1.3.3	phosphoserine phosphatase
SERSYN-PWY	1.1.1.95	3-phosphoglycerate dehydrogenase subunit
SERSYN-PWY	1.1.1.95	AT1G17745-MONOMER
SERSYN-PWY	2.6.1.1 	phosphoserine aminotransferase
SERSYN-PWY	2.6.1.52	AT4G35630-MONOMER
SERSYN-PWY	2.6.1.52	phosphoserine aminotransferase 1
SERSYN-PWY	3.1.3.3	AT1G18640-MONOMER
SHIKIMATEDEG-PWY 	1.1.5.8	quinate/shikimate dehydrogenase
SOPHOROSYLOXYDOCOSANOATE-DEG-PWY	3.1.1.6	acetylesterase
SOPHOROSYLOXYDOCOSANOATE-SYN-PWY	1.14.13 	long-chain fatty acid &omega;/(&omega;-1) hydroxylase
SOPHOROSYLOXYDOCOSANOATE-SYN-PWY	2.3.1	acetyl-CoA:sophorosyloxyfatty acid acetyltransferase
SOPHOROSYLOXYDOCOSANOATE-SYN-PWY	2.4.1 	fatty acid-&beta;-glucosyltransferase
SOPHOROSYLOXYDOCOSANOATE-SYN-PWY	2.4.1 	MONOMER-18771
SOPHOROSYLOXYDOCOSANOATE-SYN-PWY	2.4.1	UDP-glucose:hydroxy fatty acid glucosyltransferase
SPHINGOLIPID-SYN-PWY	1.1.1.102	3-dehydrosphinganine reductase
SPHINGOLIPID-SYN-PWY	1.14.18.5	sphingolipid C4-hydroxylase
SPHINGOLIPID-SYN-PWY	1.14.18.6	ceramide very long chain fatty acid hydroxylase
SPHINGOLIPID-SYN-PWY	2.3.1.24	ceramide synthase
SPHINGOLIPID-SYN-PWY	2.3.1.50	serine palmitoyltransferase
SPHINGOLIPID-SYN-PWY	2.4.1	inositol phosphorylceramide mannosyltransferase
SPHINGOLIPID-SYN-PWY	2.7.1	inositol phosphorylceramide synthase
SPHINGOLIPID-SYN-PWY	2.7.1	mannosyl-inositol-phosphoceramide inositolphosphotransferase
SUCROSEUTIL2-PWY	1.1.1	MONOMER-6110
SUCROSEUTIL2-PWY 	1.1.99.13	D-glucoside-3-dehydrogenase
SUCSYN-PWY	2.4.1.14	MONOMER-11709
SUCSYN-PWY	2.4.1.14	sucrose phosphate synthase
SUCSYN-PWY	3.1.3.24	sucrose-6F-phosphate phosphohydrolase
SUCSYN-PWY 	5.4.2.2	phosphoglucomutase
SUCUTIL-PWY	2.7.1.4	fructokinase I
SUCUTIL-PWY	2.7.1.4	MONOMER-12612
SUCUTIL-PWY	2.7.1.4	MONOMER-12630
SUCUTIL-PWY	3.2.1	MONOMER-12599
SUCUTIL-PWY	3.2.1	MONOMER-12602
SUCUTIL-PWY	3.2.1	MONOMER-12611
SUCUTIL-PWY	3.2.1	MONOMER-12621
SULFATE-CYS-PWY 	1.5.1.30 	assimilatory sulfite reductase (NADPH)
SULFATE-CYS-PWY 	1.8.4.8	3-phospho-adenylylsulfate reductase
SULFATE-CYS-PWY 	2.3.1.30	serine acetyltransferase
SULFATE-CYS-PWY 	2.5.1.47	cysteine synthase B
SULFATE-CYS-PWY 	2.5.1.47	O-acetylserine sulfhydrylase A
SULFATE-CYS-PWY 	2.7.1.25	adenylylsulfate kinase
SULFATE-CYS-PWY 	2.7.7.4	sulfate adenylyltransferase
SULFATE-CYS-PWY 	3.1.3.3	PSERPHOSPHA-MONOMER
SULFMETII-PWY	1.8.4.9	1-Apr
SULFMETII-PWY	1.8.4.9	2-Apr
SULFMETII-PWY	1.8.4.9	3-Apr
SULFMETII-PWY	1.8.4.9	adenosine 5-phosphosulphate reductase subunit
SULFMETII-PWY	1.8.4.9	assimilatory adenylyl-sulfate reductase
SULFMETII-PWY	1.8.7.1	assimilatory sulfite reductase (ferredoxin)
SULFMETII-PWY	2.7.7.4	assimilatory sulfate adenylyltransferase
TAURINEDEG-PWY	2.6.1.55	taurine:2-oxoglutarate aminotransferase
TCA	1.1.1.42	isocitrate dehydrogenase (NADP-dependent)
TEICHOICACID-PWY	2.3.1.M31	D-alanine carrier protein:PG D-alanyltransferase
TEICHOICACID-PWY	2.3.1.M32	PG:teichoic acid D-alanyltransferase
TEICHOICACID-PWY 	2.4.1.187	UDP-N-acetylmannosamine transferase
TEICHOICACID-PWY	2.4.1.52	poly(glycerol-phosphate) &alpha;-glucosyltransferase
TEICHOICACID-PWY 	2.7.7.39	glycerol-3-phosphate cytidylyltransferase
TEICHOICACID-PWY	2.7.8.12	teichoic acid poly(glycerol phosphate) polymerase
TEICHOICACID-PWY 	2.7.8.33	UDP-N-acetylglucosamine&mdash;undecaprenyl-phosphate N-acetylglucosaminephosphotransferase
TEICHOICACID-PWY 	2.7.8.q	teichoic acid glycerol-phosphate primase
TEICHOICACID-PWY 	5.1.3.14	UDP-N-acetylglucosamine 2-epimerase
TEICHOICACID-PWY	6.2.1.M5	D-alanine--[D-alanyl carrier protein] ligase
THIOREDOX-PWY	1.8.1.9	MONOMER-543
THIOREDOX-PWY	1.8.1.9	thioredoxin reductase
THIOSULFOX-PWY 	1.8.7.M1 	thiosulfate dehydrogenase
THIOSULFOX-PWY	1.8.7.M1	thiosulfate dehydrogenase
THISYNARA-PWY 	2.7.1.50 	TMP diphosphorylase/hydroxyethylthiazole kinase
THISYNARA-PWY 	2.7.6.2	thiamine pyrophosphokinase
THISYNARA-PWY 	2.7.6.2	thiamine pyrophosphokinase 1
THISYNARA-PWY 	2.7.6.2	thiamine pyrophosphokinase 2
THISYNARA-PWY 	2.8.1	thiazole synthase
THISYNARA-PWY 	3.1.3.100	thiamine phosphate phosphatase
THISYNARA-PWY 	4.1.99.17	phosphomethylpyrimidine synthase
THISYN-PWY 	2.5.1.3	THIE-MONOMER
THISYN-PWY 	2.7.1.49 	hydroxymethylpyrimidine kinase / phosphohydroxymethylpyrimidine kinase
THISYN-PWY 	2.7.4.16	thiamine monophosphate kinase
THISYN-PWY 	2.8.1.10	1-deoxy-D-xylulose 5-phosphate:thiol sulfurtransferase
THISYN-PWY 	2.8.1 	ThiI
THISYN-PWY 	2.8.1 	ThiS adenylyltransferase
THISYN-PWY 	4.1.99.17	THIC-MONOMER
THISYN-PWY 	4.1.99.19	2-iminoacetate synthase
THREOCAT-PWY 	1.1.1.103	L-threonine dehydrogenase
THREOCAT-PWY 	1.1.1.103	threonine dehydrogenase
THREOCAT-PWY 	1.2.1.10	acetaldehyde dehydrogenase (acylating)
THREOCAT-PWY 	2.3.1.29	2-amino-3-ketobutyrate CoA ligase
THREOCAT-PWY 	2.3.1 	2-ketobutyrate formate-lyase / pyruvate formate-lyase 4
THREOCAT-PWY 	2.7.2.15	propionate kinase
THREOCAT-PWY 	4.1.2 	low-specificity L-threonine aldolase
TOLSULFDEG-PWY	1.14.12.8	4-sulfobenzoate 3,4-dioxygenase oxygenase component 
TOLSULFDEG-PWY 	1.14.99	4-toluene-sulfonate monooxygenase oxygenase component 
TREDEGLOW-PWY	2.7.1.1 	MONOMER-5981
TREDEGLOW-PWY	3.2.1.93	MONOMER-5946
TREDEGLOW-PWY	3.2.1.93	TRE6PHYDRO-MONOMER
TREHALOSESYN-PWY	2.4.1.36	GDP-glucose trehalose-6-phosphate synthase
TREHALOSESYN-PWY 	2.4.1.36 	trehalose-6-phosphate synthase
TRESYN-PWY	2.4.1.15	MONOMER-17100
TRESYN-PWY	2.4.1.15	MONOMER-5983
TRESYN-PWY	2.4.1.15	TREHALOSE6PSYN-MONOMER
TRESYN-PWY	3.1.3.12 	trehalose-6-phosphate phosphatase
TRESYN-PWY	3.1.3.12	TREHALOSEPHOSPHASYN-MONOMER
TRESYN-PWY	3.1.3.12	trehalose phosphatase
TRESYN-PWY	3.1.3.12 	trehalose synthase complex
TRESYN-PWY	3.1.3.12 	YBR126C-MONOMER
TRIGLSYN-PWY	2.3.1.158	AT5G13640-MONOMER
TRIGLSYN-PWY 	2.3.1.158 	lyso-phosphatidylcholine acyltransferase
TRIGLSYN-PWY	2.3.1.158	MONOMER-12153
TRIGLSYN-PWY	2.3.1.20	diacylglycerol acyltransferase
TRIGLSYN-PWY	2.3.1	MONOMER-16590
TRIGLSYN-PWY	2.7.1.174	diacylglycerol kinase
TRIGLSYN-PWY	3.1.3.4 	diacylglycerol pyrophosphate phosphatase 1
TRIGLSYN-PWY	3.1.3.4	G3O-32855-MONOMER
TRIGLSYN-PWY 	3.1.3.4 	lipid phosphatidate phosphatase
TRIGLSYN-PWY	3.1.3.81 	YDR503C-MONOMER
TRNA-CHARGING-PWY	6.1.1.10	methionine-tRNA ligase
TRNA-CHARGING-PWY 	6.1.1.11	serine-tRNA ligase
TRNA-CHARGING-PWY	6.1.1.12	aspartate-tRNA ligase
TRNA-CHARGING-PWY	6.1.1.14	glycine-tRNA ligase
TRNA-CHARGING-PWY	6.1.1.15	proline-tRNA ligase
TRNA-CHARGING-PWY	6.1.1.16	cysteine-tRNA ligase
TRNA-CHARGING-PWY	6.1.1.18	glutamine-tRNA ligase
TRNA-CHARGING-PWY	6.1.1.19	arginine-tRNA ligase
TRNA-CHARGING-PWY	6.1.1.1	tyrosine-tRNA ligase
TRNA-CHARGING-PWY	6.1.1.20	phenylalanine-tRNA ligase
TRNA-CHARGING-PWY	6.1.1.21	histidine-tRNA ligase
TRNA-CHARGING-PWY	6.1.1.22	asparagine-tRNA ligase
TRNA-CHARGING-PWY	6.1.1.2	tryptophan-tRNA ligase
TRNA-CHARGING-PWY	6.1.1.3	threonine-tRNA ligase
TRNA-CHARGING-PWY	6.1.1.4	leucine-tRNA ligase
TRNA-CHARGING-PWY	6.1.1.5	isoleucine-tRNA ligase
TRNA-CHARGING-PWY	6.1.1.6 	lysine-tRNA ligase / Ap4A synthetase / Ap3A synthetase
TRNA-CHARGING-PWY	6.1.1.6	lysine-tRNA ligase, constitutive
TRNA-CHARGING-PWY	6.1.1.7	alanine-tRNA ligase and DNA-binding transcriptional repressor
TRNA-CHARGING-PWY	6.1.1.9	valine-tRNA ligase
TRPCAT-PWY	1.13.11.11 	tryptophan dioxygenase subunit
TRPCAT-PWY	3.5.1.9	kynurenine formamidase I
TRPCAT-PWY	3.5.1.9	kynurenine formamidase II
TRPCAT-PWY	3.5.1.9	MONOMER-7524
TRPIAACAT-PWY	1.2.3.7 	indolepyruvate decarboxylase subunit a 
TRPIAACAT-PWY	2.6.1.27	aromatic amino acid aminotransferase subunit
TRPIAACAT-PWY	2.6.1.27	MONOMER-7806
TRPIAACAT-PWY	4.1.1.74	indole-3-pyruvate decarboxylase subunit
TRPIAACAT-PWY	4.1.1.74	indolepyruvate decarboxylase subunit
TRPKYNCAT-PWY	1.1.1.110	MONOMER-7403
TRPKYNCAT-PWY	1.1.1.110	MONOMER-8131
TRPKYNCAT-PWY	2.6.1.27	MONOMER-7401
TRPKYNCAT-PWY	2.6.1.27	MONOMER-8130
TRPSYN-PWY	2.4.2.18	anthranilate synthase component II
TRPSYN-PWY	4.1.1.48	indole-3-glycerol phosphate synthase
TRPSYN-PWY	4.1.1.48	MONOMER-343
TRPSYN-PWY	4.1.3.27	anthranilate synthase &alpha; subunit 
TRPSYN-PWY	4.1.3.27	anthranilate synthase component I 
TRPSYN-PWY	4.2.1.122	tryptophan synthase &beta;2 dimer
TRPSYN-PWY	4.2.1.122	tryptophan synthase, beta subunit
TRPSYN-PWY	4.2.1.122	trytophan synthase, beta subunit
TRPSYN-PWY	4.2.1.20 	&alpha; subunit of tryptophan synthase
TRPSYN-PWY	4.2.1.20 	&alpha; subunit of tryptophan synthase 
TRPSYN-PWY	4.2.1.20 	tryptophan synthase &beta;1 dimer
TRPSYN-PWY	4.2.1.20 	tryptophan synthase &beta; subunit dimer
TRPSYN-PWY	4.2.1.20 	tryptophan synthase complex
TRPSYN-PWY	4.2.1.20 	tryptophan synthase subunit &alpha;
TRPSYN-PWY	5.3.1.24	phosphoribosylanthranilate isomerase subunit
TRYPTOPHAN-DEGRADATION-1 	1.11.2.M1 	indoleamine 2,3-dioxygenase
TRYPTOPHAN-DEGRADATION-1 	1.1.1.35 	3-hydroxy-2-methylbutyryl-CoA dehydrogenase subunit
TRYPTOPHAN-DEGRADATION-1 	1.1.1.35 	Hydroxyacyl-coenzyme A dehydrogenase, mitochondrial
TRYPTOPHAN-DEGRADATION-1 	1.13.11.11 	tryptophan 2,3-dioxygenase
TRYPTOPHAN-DEGRADATION-1 	1.13.11.11 	YJR078W-MONOMER
TRYPTOPHAN-DEGRADATION-1 	1.13.11.6	3-hydroxyanthranilate 3,4-dioxygenase
TRYPTOPHAN-DEGRADATION-1 	1.14.13.9	kynurenine 3-monooxygenase
TRYPTOPHAN-DEGRADATION-1 	1.14.13.9	YBL098W-MONOMER
TRYPTOPHAN-DEGRADATION-1 	3.5.1.9	arylformamidase
TRYPTOPHAN-DEGRADATION-1 	3.5.1.9	YJL060W-MONOMER
TRYPTOPHAN-DEGRADATION-1 	3.7.1.3	kynureninase
TRYPTOPHAN-DEGRADATION-1 	3.7.1.3	YLR231C-MONOMER
TRYPTOPHAN-DEGRADATION-1 	4.1.1.45	aminocarboxymuconate-semialdehyde decarboxylase
TYRFUMCAT-PWY	1.13.11.27	4-hydroxyphenylpyruvate dioxygenase subunit
TYRFUMCAT-PWY	1.13.11.5	Homogentisate 1,2-dioxygenase
TYRFUMCAT-PWY	1.13.11.5	MONOMER-12040
TYRFUMCAT-PWY	1.13.11.5	MONOMER-12041
TYRFUMCAT-PWY	2.5.1.18 	Maleylacetoacetate isomerase
TYRFUMCAT-PWY 	2.6.1.57 	Tyrosine aminotransferase
TYRFUMCAT-PWY 	2.6.1.57 	tyrosine aminotransferase subunit
TYRFUMCAT-PWY	3.7.1.2	Fumarylacetoacetase
TYRFUMCAT-PWY	3.7.1.2	MONOMER-12048
TYRFUMCAT-PWY	5.2.1.2	MONOMER-12045
TYRFUMCAT-PWY	5.2.1.2	MONOMER-12046
UBISYN-PWY 	4.1.3.40	CHORPYRLY-MONOMER
UDPNACETYLGALSYN-PWY	2.3.1.4	glucosamine-6-phosphate N-acetyltransferase 1
UDPNACETYLGALSYN-PWY	2.3.1.4	glucosamine 6-phosphate N-acetyltransferase
UDPNACETYLGALSYN-PWY	2.6.1.16	AT3G24090-MONOMER
UDPNACETYLGALSYN-PWY	2.6.1.16	glucosamine-6-phosphate synthase subunit
UDPNACETYLGALSYN-PWY	2.6.1.16	glutamine:fructose-6-phosphate amidotransferase I
UDPNACETYLGALSYN-PWY	2.6.1.16	glutamine:fructose 6-phosphate amidotransferase subunit
UDPNACETYLGALSYN-PWY 	2.7.7.23 	N-acetylglucosamine-1-phosphate uridylyltransferase
UDPNACETYLGALSYN-PWY	2.7.7.23	UDP-N-acetylglucosamine pyrophosphorylase subunit
UDPNACETYLGALSYN-PWY 	2.7.7.64 	UDP-sugar pyrophosphorylase
UDPNACETYLGALSYN-PWY	5.4.2.3	phosphoacetylglucosamine mutase
URDEGR-PWY 	3.5.1.116	AT5G43600-MONOMER
URDEGR-PWY 	3.5.1.116	MONOMER-13516
URDEGR-PWY	3.5.1.116	ureidoglycolate amidohydrolase
URDEGR-PWY 	3.5.1.5	urease
URDEGR-PWY	3.5.2.5	MONOMER-11658
URDEGR-PWY 	3.5.3.26	AT4G17050-MONOMER
URDEGR-PWY 	3.5.3.9	allantoate amidohydrolase
URDEGR-PWY 	3.5.3.9	allantoate amidohydrolase monomer
URDEGR-PWY 	3.5.3.9	MONOMER-13525
URDEGR-PWY 	4.3.2.3	ureidoglycolate urea-lyase subunit
URSIN-PWY 	1.7.3.3	factor-independent urate hydroxylase
URSIN-PWY 	2.4.2.1	purine nucleoside phosphorylase
URSIN-PWY 	3.5.2.17	MONOMER-13506
VALDEG-PWY	1.1.1.31	3-hydroxyisobutyrate dehydrogenase subunit
VALDEG-PWY	1.2.1.27	methylmalonate-semialdehyde dehydrogenase subunit
VALDEG-PWY	2.6.1.19 	4-aminobutyrate aminotransferase, mitochondrial
VALDEG-PWY	3.1.2.4	MONOMER-11699
VALDEG-PWY	4.2.1.17	enoyl-CoA hydratase subunit
VALSYN-PWY 	1.1.1.383 	acetohydroxy acid isomeroreductase
VALSYN-PWY 	2.6.1.6 	branched-chain amino acid aminotransferase
VALSYN-PWY 	4.2.1.9	dihydroxy acid dehydratase' > ec.to.pwy
}

function pwyhierarchyFunction {
	echo 'PWY-6519	8-amino-7-oxononanoate biosynthesis I	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|BIOTIN-SYN|7-Keto-8-aminopelargonate-Biosynthesis
PWY-7641	5-hexynoate biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|Acetlylenic-Fatty-Acid-Biosynthesis
PWY-6013	crepenynate biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis
PWY-6333	acetaldehyde biosynthesis I	Biosynthesis|Carbohydrates-Biosynthesis|Acetaldehyde-Biosynthesis
PWY-6401	hispidol and hispidol 4-O-&beta;-D-glucoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|AURONE-SYN
PWY-1901	aurone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|AURONE-SYN
PWY-5305	bixin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TETRATERPENOID-SYN|APOCAROTENOID-SYN
PWY-6681	neurosporaxanthin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TETRATERPENOID-SYN|APOCAROTENOID-SYN
PWY-5398	crocetin esters biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TETRATERPENOID-SYN|APOCAROTENOID-SYN
PWY-7539	6-hydroxymethyl-dihydropterin diphosphate biosynthesis III (Chlamydia)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Folate-Biosynthesis|6-HM-Dihydropterin-PP-Biosynthesis
PWY-6797	6-hydroxymethyl-dihydropterin diphosphate biosynthesis II (archaea)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Folate-Biosynthesis|6-HM-Dihydropterin-PP-Biosynthesis
PWY-6147	6-hydroxymethyl-dihydropterin diphosphate biosynthesis I	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Folate-Biosynthesis|6-HM-Dihydropterin-PP-Biosynthesis
PWY-7766	heme biosynthesis IV (Gram-positive bacteria)	Biosynthesis|Cofactor-Biosynthesis|Porphyrin-Compounds-Biosynthesis|HEME-SYN
PWY-7554	heme d1 biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Porphyrin-Compounds-Biosynthesis|HEME-SYN
HEMESYN2-PWY	heme biosynthesis II (anaerobic)	Biosynthesis|Cofactor-Biosynthesis|Porphyrin-Compounds-Biosynthesis|HEME-SYN
PWY-7552	heme biosynthesis III (from siroheme)	Biosynthesis|Cofactor-Biosynthesis|Porphyrin-Compounds-Biosynthesis|HEME-SYN
HEME-BIOSYNTHESIS-II	heme biosynthesis I (aerobic)	Biosynthesis|Cofactor-Biosynthesis|Porphyrin-Compounds-Biosynthesis|HEME-SYN
PWY-5068	chlorophyll cycle	Biosynthesis|Cofactor-Biosynthesis|Porphyrin-Compounds-Biosynthesis|Chlorophyll-Biosynthesis|Chlorophyll-a-Biosynthesis
PWY-7764	chlorophyll a biosynthesis III	Biosynthesis|Cofactor-Biosynthesis|Porphyrin-Compounds-Biosynthesis|Chlorophyll-Biosynthesis|Chlorophyll-a-Biosynthesis
PWY-5064	chlorophyll a biosynthesis II	Biosynthesis|Cofactor-Biosynthesis|Porphyrin-Compounds-Biosynthesis|Chlorophyll-Biosynthesis|Chlorophyll-a-Biosynthesis
PWY-5086	chlorophyll a biosynthesis I	Biosynthesis|Cofactor-Biosynthesis|Porphyrin-Compounds-Biosynthesis|Chlorophyll-Biosynthesis|Chlorophyll-a-Biosynthesis
DESULFONATION-PWY	benzenesulfonate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Sulfoaromatics-Degradation
PWY-6041	4-sulfocatechol degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Sulfoaromatics-Degradation
PWY-7521	3-(4-sulfophenyl)butanoate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Sulfoaromatics-Degradation
PWY-5165	4-toluenesulfonate degradation II	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Sulfoaromatics-Degradation|4-Toluenesulfonate-Degradation
TOLSULFDEG-PWY	4-toluenesulfonate degradation I	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Sulfoaromatics-Degradation|4-Toluenesulfonate-Degradation
PWY-5641	2-nitrotoluene degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nitroaromatic-Degradation|Nitro-Toluene-Degradation
PWY-6051	2,4,6-trinitrotoluene degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nitroaromatic-Degradation|Nitro-Toluene-Degradation
PWY-5644	4-nitrotoluene degradation II	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nitroaromatic-Degradation|Nitro-Toluene-Degradation|4-Nitrotoluene-Degradation
P421-PWY	4-nitrotoluene degradation I	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nitroaromatic-Degradation|Nitro-Toluene-Degradation|4-Nitrotoluene-Degradation
PWY-5754	4-hydroxybenzoate biosynthesis I (eukaryotes)	Biosynthesis|AROMATIC-COMPOUNDS-BIOSYN|4-Hydroxybenzoate-Biosynthesis
PWY-6435	4-hydroxybenzoate biosynthesis V	Biosynthesis|AROMATIC-COMPOUNDS-BIOSYN|4-Hydroxybenzoate-Biosynthesis
PWY-6431	4-hydroxybenzoate biosynthesis IV	Biosynthesis|AROMATIC-COMPOUNDS-BIOSYN|4-Hydroxybenzoate-Biosynthesis
PWY-5755	4-hydroxybenzoate biosynthesis II (microbes)	Biosynthesis|AROMATIC-COMPOUNDS-BIOSYN|4-Hydroxybenzoate-Biosynthesis
PWY0-1586	peptidoglycan maturation (meso-diaminopimelate containing)	Biosynthesis|Cell-Structure-Biosynthesis|Cell-Wall-Biosynthesis|Peptidoglycan-Biosynthesis
PWY-7723	bacterial bioluminescence	Bioluminescence
PWY-5751	phenylethanol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN
PWY-5787	oligomeric urushiol biosynthesis	Biosynthesis|AROMATIC-COMPOUNDS-BIOSYN
PWY-5886	4-hydroxyphenylpyruvate biosynthesis	Biosynthesis|AROMATIC-COMPOUNDS-BIOSYN
PWY-5901	2,3-dihydroxybenzoate biosynthesis	Biosynthesis|AROMATIC-COMPOUNDS-BIOSYN
PWY-6320	phaselate biosynthesis	Biosynthesis|AROMATIC-COMPOUNDS-BIOSYN
PWY-6323	benzoylanthranilate biosynthesis	Biosynthesis|AROMATIC-COMPOUNDS-BIOSYN
PWY-6406	salicylate biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|BENZENOID-SYN|Salicylate-Biosynthesis
PWY-6457	trans-cinnamoyl-CoA biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|CINNAMATE-SYN
PWY-6458	benzoyl-CoA biosynthesis	Biosynthesis|AROMATIC-COMPOUNDS-BIOSYN
PWY-6621	salicylate glucosides biosynthesis I	Super-Pathways
PWY-6623	salicylate glucosides biosynthesis II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|BENZENOID-SYN|Salicylate-Biosynthesis
PWY-6624	salicylate glucosides biosynthesis III	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|BENZENOID-SYN|Salicylate-Biosynthesis
PWY-6707	gallate biosynthesis	Biosynthesis|AROMATIC-COMPOUNDS-BIOSYN
PWY-6752	o-diquinones biosynthesis	Biosynthesis|AROMATIC-COMPOUNDS-BIOSYN
PWY-6762	salicylate glucosides biosynthesis IV	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|BENZENOID-SYN|Salicylate-Biosynthesis
PWY-6930	phenolic malonylglucosides biosynthesis	Metabolic-Clusters
PWY-7303	3-dimethylallyl-4-hydroxybenzoate biosynthesis	Biosynthesis|AROMATIC-COMPOUNDS-BIOSYN
PWY-981	salicylate biosynthesis II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|BENZENOID-SYN|Salicylate-Biosynthesis
PWY-6164	3-dehydroquinate biosynthesis I	Biosynthesis|AROMATIC-COMPOUNDS-BIOSYN|3-Dehydroquinate-Biosynthesis
PWY-6160	3-dehydroquinate biosynthesis II (archaea)	Biosynthesis|AROMATIC-COMPOUNDS-BIOSYN|3-Dehydroquinate-Biosynthesis
PWY-6244	bergamotene biosynthesis II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN|Bergamotene-Biosynthesis
PWY-6243	bergamotene biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN|Bergamotene-Biosynthesis
PWY-6444	benzoate biosynthesis II (CoA-independent, non-&beta;-oxidative)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|BENZENOID-SYN|Benzoate-Biosynthesis
PWY-6766	salicin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|BENZENOID-SYN|Benzoate-Biosynthesis
PWY-6763	salicortin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|BENZENOID-SYN|Benzoate-Biosynthesis
PWY-5400	amaranthin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|BETALAIN-ALKALOIDS
PWY-5399	betacyanin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|BETALAIN-ALKALOIDS
PWY-5404	betaxanthin biosynthesis (via dopaxanthin)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|BETALAIN-ALKALOIDS
PWY-5439	betacyanin biosynthesis (via dopamine)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|BETALAIN-ALKALOIDS
PWY-5394	betalamic acid biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|BETALAIN-ALKALOIDS
PWY-5403	betaxanthin biosynthesis (via dopamine)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|BETALAIN-ALKALOIDS
PWY-5426	betaxanthin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|BETALAIN-ALKALOIDS
PWY-7198	pyrimidine deoxyribonucleotides de novo biosynthesis IV	Biosynthesis|Nucleotide-Biosynthesis|Deoxyribonucleotide-Biosynthesis|Pyrimid-Deoxyribonucleot-De-Novo-Biosyn
PWY-7678	anthocyanidin sambubioside biosynthesis	Metabolic-Clusters
PWY-7224	purine deoxyribonucleosides salvage	Metabolic-Clusters
PWY-7799	Arg/N-end rule pathway (eukaryotic)	Macromolecule-Modification|Protein-Modification
PWY-1121	suberin monomers biosynthesis	Metabolic-Clusters
PWY-7321	ecdysteroid metabolism (arthropods)	Metabolic-Clusters
PWY-7802	N-end rule pathway II (prokaryotic)	Macromolecule-Modification|Protein-Modification
PWY-4201	volatile cinnamoic ester biosynthesis	Metabolic-Clusters
PWY-6801	volatile esters biosynthesis (during fruit ripening)	Metabolic-Clusters
PWY-7491	podophyllotoxin glucosides metabolism	Metabolic-Clusters
TRNA-CHARGING-PWY	tRNA charging	Metabolic-Clusters
PWY-5406	divinyl ether biosynthesis I	Metabolic-Clusters
PWY-6035	2,3-cis-flavanols biosynthesis	Metabolic-Clusters
PWY-7197	pyrimidine deoxyribonucleotide phosphorylation	Biosynthesis|Nucleotide-Biosynthesis|PYR-NUC-SYN|Pyrimidine-Nucleotide-Salvage
PWY-7656	Spodoptera littoralis pheromone biosynthesis	Metabolic-Clusters
PWY-7210	pyrimidine deoxyribonucleotides biosynthesis from CTP	Biosynthesis|Nucleotide-Biosynthesis|PYR-NUC-SYN|Pyrimidine-Nucleotide-Salvage
PWY-7696	citreoisocoumarin and bikisocoumarin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY-7261	anthocyanidin 3-malylglucoside biosynthesis (acyl-glucose dependent)	Metabolic-Clusters
PWY-7801	N-end rule pathway I (prokaryotic)	Macromolecule-Modification|Protein-Modification
PWY-2902	cytokinin-O-glucosides biosynthesis	Metabolic-Clusters
PWY-6733	sporopollenin precursors biosynthesis	Metabolic-Clusters
PWY-7424	sterol:steryl ester interconversion (yeast)	Metabolic-Clusters
PWY-5286	anthocyanidin sophoroside metabolism	Metabolic-Clusters
PWY-6829	tRNA methylation (yeast)	Biosynthesis|Nucleotide-Biosynthesis|Nucleic-Acid-Processing
PWY-6029	2,3-trans-flavanols biosynthesis	Metabolic-Clusters
PWY-7184	pyrimidine deoxyribonucleotides de novo biosynthesis I	Biosynthesis|Nucleotide-Biosynthesis|Deoxyribonucleotide-Biosynthesis|Pyrimid-Deoxyribonucleot-De-Novo-Biosyn
PWY-7654	(8E,10E)-dodeca-8,10-dienol biosynthesis	Metabolic-Clusters
PWY-7679	anthocyanidin acylglucoside and acylsambubioside biosynthesis	Metabolic-Clusters
OANTIGEN-PWY	O-antigen building blocks biosynthesis (E. coli)	Super-Pathways
PWY-6545	pyrimidine deoxyribonucleotides de novo biosynthesis III	Biosynthesis|Nucleotide-Biosynthesis|Deoxyribonucleotide-Biosynthesis|Pyrimid-Deoxyribonucleot-De-Novo-Biosyn
PWY-7251	pentacyclic triterpene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-7800	Ac/N-end rule pathway	Macromolecule-Modification|Protein-Modification
PWY-6603	dicranin biosynthesis	Metabolic-Clusters
PWY-7423	bombykol biosynthesis	Metabolic-Clusters
PWY-83	monolignol glucosides biosynthesis	Metabolic-Clusters
PWY-5045	pinosylvin metabolism	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|STILBENE-SYN
PWY-7498	phenylpropanoids methylation (ice plant)	Metabolic-Clusters
PWY-7120	esterified suberin biosynthesis	Metabolic-Clusters
PWY-7637	2,2-dihydroxyketocarotenoids biosynthesis	Metabolic-Clusters
PWY-6219	indole-3-acetyl-amide conjugate biosynthesis	Metabolic-Clusters
PWY-6030	serotonin and melatonin biosynthesis	Biosynthesis|HORMONE-SYN
PWY-6076	vitamin D3 biosynthesis	Biosynthesis|HORMONE-SYN
PWY-6241	thyroid hormone biosynthesis	Biosynthesis|HORMONE-SYN
PWY-6575	juvenile hormone III biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN|JH-III-Biosynthesis
PWY-6650	juvenile hormone III biosynthesis II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN|JH-III-Biosynthesis
PWY-7299	progesterone biosynthesis	Biosynthesis|HORMONE-SYN
PWY-7300	ecdysone and 20-hydroxyecdysone biosynthesis	Biosynthesis|HORMONE-SYN
PWY-7305	superpathway of steroid hormone biosynthesis	Super-Pathways
PWY-7306	estradiol biosynthesis II	Biosynthesis|HORMONE-SYN
PWY-7410	ipsdienol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|MONOTERPENOID-SYN
PWY66-301	catecholamine biosynthesis	Biosynthesis|HORMONE-SYN
PWY66-374	C20 prostanoid biosynthesis	Biosynthesis|HORMONE-SYN
PWY66-375	leukotriene biosynthesis	Biosynthesis|HORMONE-SYN
PWY66-377	pregnenolone biosynthesis	Biosynthesis|HORMONE-SYN
PWY66-378	androgen biosynthesis	Biosynthesis|HORMONE-SYN
PWY66-380	estradiol biosynthesis I (via estrone)	Biosynthesis|HORMONE-SYN
PWY66-381	glucocorticoid biosynthesis	Biosynthesis|HORMONE-SYN
PWY66-382	mineralocorticoid biosynthesis	Biosynthesis|HORMONE-SYN
PWY66-392	lipoxin biosynthesis	Biosynthesis|HORMONE-SYN
PWY66-393	aspirin-triggered lipoxin biosynthesis	Biosynthesis|HORMONE-SYN
2ASDEG-PWY	orthanilate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY66-394	aspirin triggered resolvin E biosynthesis	Biosynthesis|HORMONE-SYN
3-HYDROXYPHENYLACETATE-DEGRADATION-PWY	4-hydroxyphenylacetate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY66-395	aspirin triggered resolvin D biosynthesis	Biosynthesis|HORMONE-SYN
4-HYDROXYMANDELATE-DEGRADATION-PWY	4-hydroxymandelate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY66-397	resolvin D biosynthesis	Biosynthesis|HORMONE-SYN
4TOLCARBDEG-PWY	4-toluenecarboxylate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
BENZCOA-PWY	anaerobic aromatic compound degradation (Thauera aromatica)	Super-Pathways
M-CRESOL-DEGRADATION-PWY	<I>m</I>-cresol degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
METHYLGALLATE-DEGRADATION-PWY	methylgallate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
P342-PWY	orcinol degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
P343-PWY	resorcinol degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
P345-PWY	aldoxime degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
P661-PWY	dibenzo-<I>p</I>-dioxin degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
P662-PWY	dibenzofuran degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PARATHION-DEGRADATION-PWY	parathion degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-142	m-xylene degradation (anaerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY0-1319	CDP-diacylglycerol biosynthesis II	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis|CDP-diacylglycerol-Biosynthesis
PWY-2421	indole-3-acetate degradation VIII (bacterial)	Degradation|HORMONE-DEG|PLANT-HORMONE-DEG|AUXINS-DEGRADATION
PWY-5981	CDP-diacylglycerol biosynthesis III	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis|CDP-diacylglycerol-Biosynthesis
PWY-2504	superpathway of aromatic compound degradation via 3-oxoadipate	Super-Pathways
PWY-5667	CDP-diacylglycerol biosynthesis I	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis|CDP-diacylglycerol-Biosynthesis
PWY-481	ethylbenzene degradation (anaerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-5062	superpathway of nicotinate degradation	Super-Pathways
PWY-5163	<I>p</I>-cumate degradation to 2-oxopent-4-enoate	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-5169	cyanurate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-5266	<I>p</I>-cymene degradation	Super-Pathways
PWY-5273	p-cumate degradation	Super-Pathways
PWY-5450	benzene degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-5489	methyl parathion degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-5490	paraoxon degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-6182	superpathway of salicylate degradation	Super-Pathways
PWY-6184	methylsalicylate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-6185	4-methylcatechol degradation (ortho cleavage)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-6210	2-aminophenol degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-6339	syringate degradation	Super-Pathways
PWY-6340	5,5-dehydrodivanillate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-6343	ferulate degradation	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|Vanillin-Biosynthesis
PWY-6532	diphenylamine degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-6533	aniline degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-6550	carbazole degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-6781	chlorogenic acid degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-681	dibenzothiophene desulfurization	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-6941	styrene degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-6954	superpathway of aromatic compound degradation via 2-oxopent-4-enoate	Super-Pathways
PWY-6956	naphthalene degradation to acetyl-CoA	Super-Pathways
PWY-7002	4-hydroxyacetophenone degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-7006	4-amino-3-hydroxybenzoate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-7008	2-hydroxybiphenyl degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-7009	2,2-dihydroxybiphenyl degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-7010	2-propylphenol degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-7011	2-isopropylphenol degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-7081	4-aminophenol degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-721	3-methylquinoline degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-741	<I>p</I>-cymene degradation to <I>p</I>-cumate	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-7430	indole degradation to anthranil and anthranilate	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-7480	2,3-dihydroxybenzoate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-7557	(-)-dehydrodiconiferyl alcohol degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-7698	2,5-xylenol and 3,5-xylenol degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-7747	diphenyl ethers degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-7795	terephthalate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY-7825	juglone degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY0-1277	3-phenylpropanoate and 3-(3-hydroxyphenyl)propanoate degradation	Super-Pathways
PWY5F9-12	biphenyl degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
PWY5F9-3233	phthalate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION
12DICHLORETHDEG-PWY	1,2-dichloroethane degradation	Degradation|CHLORINATED-COMPOUNDS-DEG
GAMMAHEXCHLORDEG-PWY	&gamma;-hexachlorocyclohexane degradation	Degradation|CHLORINATED-COMPOUNDS-DEG
PCEDEG-PWY	tetrachloroethene degradation	Degradation|CHLORINATED-COMPOUNDS-DEG
PCPDEG-PWY	pentachlorophenol degradation	Degradation|CHLORINATED-COMPOUNDS-DEG
PWY-5849	menaquinol-6 biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Menaquinone-Biosynthesis
PWY-5822	trichloroethene degradation	Degradation|CHLORINATED-COMPOUNDS-DEG
PWY-5891	menaquinol-11 biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Menaquinone-Biosynthesis
PWY-5839	menaquinol-7 biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Menaquinone-Biosynthesis
PWY-6085	2,4-dichlorophenoxyacetatedegradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation
PWY-5890	menaquinol-10 biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Menaquinone-Biosynthesis
PWY-6086	4-chloro-2-methylphenoxyacetate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation
PWY-6102	5-chloro-3-methyl-catechol degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation
PWY-5895	menaquinol-13 biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Menaquinone-Biosynthesis
PWY-6107	chlorosalicylate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation
PWY-6178	2,4,6-trichlorophenol degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation
PWY-5844	menaquinol-9 biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Menaquinone-Biosynthesis
PWY-6197	chlorinated phenols degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation
PWY-6200	2,4,5-trichlorophenoxyacetate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation
PWY-5892	menaquinol-12 biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Menaquinone-Biosynthesis
PWY-7496	linuron degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation
PWY-7512	3,5,6-trichloro-2-pyridinol degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation
MENAQUINONESYN-PWY	menaquinol-8 biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Menaquinone-Biosynthesis
PWY-7771	butachlor degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation
PWY-7823	chlorzoxazone degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation
PWY-6104	3-chlorotoluene degradation II	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorotoluene-Degradation|3-Chlorotoluene-Degradation
PWY-6103	3-chlorotoluene degradation I	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorotoluene-Degradation|3-Chlorotoluene-Degradation
PWY-6084	3,5-dichlorocatechol degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorocatechol-Degradation
PWY-6087	4-chlorocatechol degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorocatechol-Degradation
PWY-6093	4,5-dichlorocatechol degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorocatechol-Degradation
PWY-6094	3,4,6-trichlorocatechol degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorocatechol-Degradation
PWY-6193	3-chlorocatechol degradation II (ortho)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorocatechol-Degradation|3-Chlorocatechol-Degradation
PWY-6089	3-chlorocatechol degradation I (ortho)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorocatechol-Degradation|3-Chlorocatechol-Degradation
PWY-6214	3-chlorocatechol degradation III (meta pathway)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorocatechol-Degradation|3-Chlorocatechol-Degradation
PWY-6215	4-chlorobenzoate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorobenzoate-Degradation
PWY-6217	3,4-dichlorobenzoate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorobenzoate-Degradation
PWY-6221	2-chlorobenzoate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorobenzoate-Degradation
PWY-6216	3-chlorobenzoate degradation II (via protocatechuate)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorobenzoate-Degradation|3-Chlorobenzoate-Degradation
PWY-6088	3-chlorobenzoate degradation I (via chlorocatechol)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorobenzoate-Degradation|3-Chlorobenzoate-Degradation
PWY-6228	3-chlorobenzoate degradation III (via gentisate)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorobenzoate-Degradation|3-Chlorobenzoate-Degradation
PWY-7131	CMP-legionaminate biosynthesis II	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|CMP-Sugar-Biosynthesis|CMP-Legionaminate-Biosynthesis
PWY-6749	CMP-legionaminate biosynthesis I	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|CMP-Sugar-Biosynthesis|CMP-Legionaminate-Biosynthesis
PWY-7511	protein ubiquitylation	Macromolecule-Modification|Protein-Modification
PWY-7798	protein S-nitrosylation and denitrosylation	Macromolecule-Modification|Protein-Modification
PWY-7031	protein N-glycosylation (bacterial)	Macromolecule-Modification|Protein-Modification|Protein-Glycosylation
PWY-7433	mucin core 1 and core 2 O-glycosylation	Macromolecule-Modification|Protein-Modification|Protein-Glycosylation
PWY-7437	protein O-[N-acetyl]-glucosylation	Macromolecule-Modification|Protein-Modification|Protein-Glycosylation
PWY-7426	mannosyl-glycoprotein N-acetylglucosaminyltransferases	Macromolecule-Modification|Protein-Modification|Protein-Glycosylation
PWY-7435	mucin core 3 and core 4 O-glycosylation	Macromolecule-Modification|Protein-Modification|Protein-Glycosylation
PWY-7661	protein <I>N-glycosylation (Haloferax volcanii)	Macromolecule-Modification|Protein-Modification|Protein-Glycosylation
MANNOSYL-CHITO-DOLICHOL-BIOSYNTHESIS	protein N-glycosylation (eukaryotic, high mannose)	Macromolecule-Modification|Protein-Modification|Protein-Glycosylation
PWY-7037	protein O-glycosylation (bacterial)	Macromolecule-Modification|Protein-Modification|Protein-Glycosylation
PWY-7434	terminal O-glycans residues modification	Macromolecule-Modification|Protein-Modification|Protein-Glycosylation
PWY-7658	protein <I>N-glycosylation (Methanococcus voltae)	Macromolecule-Modification|Protein-Modification|Protein-Glycosylation
PWY-6279	myxol-2 fucoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN
PWY-6287	neurosporene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN
PWY-5174	capsanthin and capsorubin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN
PWY-5291	canthaxanthin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN
PWY-7636	astaxanthin biosynthesis (flowering plants)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN
PWY-5944	zeaxanthin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN
PWY-5947	lutein biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN
PWY18HP-2	decaprenoxanthin and decaprenoxanthin diglucoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN
PWY-6286	spheroidene and spheroidenone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN
PWY-6581	spirilloxanthin and 2,2-diketo-spirilloxanthin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN
PWY-7175	nostoxanthin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN
PWY-5288	astaxanthin biosynthesis (bacteria, fungi, algae)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN
PWY-1501	mandelate degradation I	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|MANDELATE-DEG
PWY-7623	astaxanthin dirhamnoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN
PWY-6571	dermatan sulfate biosynthesis	Super-Pathways
PWY-5943	&beta;-carotene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN
PWY-6975	superpathway of erythromycin biosynthesis (without sugar biosynthesis)	Super-Pathways
PWY-7638	echinenone and zeaxanthin biosynthesis (Synechocystis)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN
PWY-7186	superpathway of scopolin and esculin biosynthesis	Super-Pathways
PWY-5946	&delta;-carotene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN
PWY-7283	wybutosine biosynthesis	Super-Pathways
PWY-6280	synechoxanthin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN
PWY-7622	UDP-galactofuranose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
PWY-6288	zeaxanthin-&beta;-D-diglucoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN
ALL-CHORISMATE-PWY	superpathway of chorismate metabolism	Super-Pathways
PWY-6809	neoxanthin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN
NADSYN-PWY	NAD biosynthesis II (from tryptophan)	Super-Pathways
PWY-5175	lactucaxanthin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN
PWY-2229	superpathway of pterocarpan biosynthesis (via formononetin)	Super-Pathways
PWY-7591	okenone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN
PWY-5005	biotin biosynthesis II	Super-Pathways
PWY-6475	trans-lycopene biosynthesis II (plants)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|CAROTENOID-SYN|Lycopene-Biosynthesis
PWY-5179	toluene degradation V (aerobic) (<I>via</I> toluene-<I>cis</I>-diol)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|TOLUENE-DEG
PWY-7108	erythromycin A biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis|Macrolides-Biosynthesis
PWY-5313	superpathway of anthocyanin biosynthesis (from cyanidin and cyanidin 3-O-glucoside)	Super-Pathways
PWY-7412	mycinamicin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis|Macrolides-Biosynthesis
PWY-5897	superpathway of menaquinol-11 biosynthesis	Super-Pathways
PWY-7422	methymycin, neomethymycin and novamethymycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis|Macrolides-Biosynthesis
PWY-6126	superpathway of adenosine nucleotides de novo biosynthesis II	Super-Pathways
PWY-6971	oleandomycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis|Macrolides-Biosynthesis
PWY-6285	superpathway of fatty acids biosynthesis (E. coli)	Super-Pathways
PWY-7106	erythromycin D biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis|Macrolides-Biosynthesis
PWY-6415	L-ascorbate biosynthesis V	Super-Pathways
PWY-7421	narbomycin, pikromycin and novapikromycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis|Macrolides-Biosynthesis
PWY-6446	benzoate biosynthesis III (CoA-dependent, non-&beta;-oxidative)	Super-Pathways
PWY-7677	rosamicin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis|Macrolides-Biosynthesis
PWY-6887	kauralexin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHYTOALEXIN-SYN|TERPENOID-PHYTOALEXINS
PWY-7109	megalomicin A biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis|Macrolides-Biosynthesis
PWY-7000	kanamycin biosynthesis	Super-Pathways
PWY-7415	tylosin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis|Macrolides-Biosynthesis
PWY-7200	superpathway of pyrimidine deoxyribonucleoside salvage	Super-Pathways
PWY-7624	nystatin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis|Macrolides-Biosynthesis
PWY-7317	superpathway of dTDP-glucose-derived O-antigen building blocks biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis
PWY-7438	superpathway of dTDP-glucose-derived antibiotic building blocks biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis
PWY-7809	superpathway of tetracycline and oxytetracycline biosynthesis	Super-Pathways
ALANINE-DEG3-PWY	L-alanine degradation III	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|ALANINE-DEG
PWY0-845	superpathway of pyridoxal 5-phosphate biosynthesis and salvage	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Vitamin-B6-Biosynthesis
ALADEG-PWY	L-alanine degradation I	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|ALANINE-DEG
SULFATE-CYS-PWY	superpathway of sulfate assimilation and cysteine biosynthesis	Super-Pathways
ALACAT2-PWY	L-alanine degradation II (to D-lactate)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|ALANINE-DEG
ECASYN-PWY	enterobacterial common antigen biosynthesis	Super-Pathways
PWY1-2	L-alanine degradation IV	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|ALANINE-DEG
PHOSLIPSYN2-PWY	superpathway of phospholipid biosynthesis II (plants)	Super-Pathways
ORNARGDEG-PWY	superpathway of L-arginine and L-ornithine degradation	Super-Pathways
PWY-5053	superpathway of gibberellin GA12 biosynthesis	Super-Pathways
PWY-5742	L-arginine degradation IX (arginine:pyruvate transaminase pathway)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|ARGININE-DEG
PWY-5182	toluene degradation II (aerobic) (via 4-methylcatechol)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|TOLUENE-DEG
ARGASEDEG-PWY	L-arginine degradation I (arginase pathway)	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|GLUTAMATE-SYN
PWY-5717	superpathway of formononetin derivative biosynthesis	Super-Pathways
ARGDEG-PWY	superpathway of L-arginine, putrescine, and 4-aminobutanoate degradation	Super-Pathways
PWY-5903	bacillibactin biosynthesis	Super-Pathways
AST-PWY	L-arginine degradation II (AST pathway)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|ARGININE-DEG
PWY-6145	superpathway of CMP-sialic acids biosynthesis	Super-Pathways
PWY-5024	L-arginine degradation XI	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|ARGININE-DEG
PWY-6612	superpathway of tetrahydrofolate biosynthesis	Super-Pathways
PWY0-823	L-arginine degradation III (arginine decarboxylase/agmatinase pathway)	Biosynthesis|Polyamine-Biosynthesis
PWY-7021	superpathway of neomycin biosynthesis	Super-Pathways
ARG-PRO-PWY	L-arginine degradation VI (arginase 2 pathway)	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|PROLINE-SYN
PWY-7211	superpathway of pyrimidine deoxyribonucleotides de novo biosynthesis	Super-Pathways
ARGDEG-IV-PWY	L-arginine degradation VIII (arginine oxidase pathway)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|ARGININE-DEG
PWY-7332	superpathway of UDP-<I>N-acetylglucosamine-derived O-antigen building blocks biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
ARGDEGRAD-PWY	L-arginine degradation V (arginine deiminase pathway)	Super-Pathways
PWY-7509	cardiolipin and phosphatidylethanolamine biosynthesis (Xanthomonas)	Super-Pathways
PWY-4983	L-citrulline-nitric oxide cycle	Biosynthesis|Metabolic-Regulators|Nitric-Oxide-Biosynthesis
PWY3O-1109	superpathway of 4-hydroxybenzoate biosynthesis (yeast)	Super-Pathways
PWY-7523	L-arginine degradation XII	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|ARGININE-DEG
THISYNARA-PWY	superpathway of thiamine diphosphate biosynthesis III (eukaryotes)	Super-Pathways
ARG-GLU-PWY	L-arginine degradation VII (arginase 3 pathway)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|ARGININE-DEG
FASYN-INITIAL-PWY	superpathway of fatty acid biosynthesis initiation (E. coli)	Super-Pathways
ARGDEG-III-PWY	L-arginine degradation IV (arginine decarboxylase/agmatine deiminase pathway)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|ARGININE-DEG
KDO-NAGLIPASYN-PWY	superpathway of (Kdo)<SUB>2</SUB>-lipid A biosynthesis	Super-Pathways
ARGDEG-V-PWY	L-arginine degradation X (arginine monooxygenase pathway)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|ARGININE-DEG
P165-PWY	superpathway of purines degradation in plants	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Urate-Degradation
ASPARAGINE-DEG1-PWY-1	L-asparagine degradation III (mammalian)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|ASPARAGINE-DEG
POLYAMSYN-PWY	superpathway of polyamine biosynthesis I	Super-Pathways
ASPARAGINE-DEG1-PWY	L-asparagine degradation I	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|ASPARAGINE-DEG
PWY-5507	adenosylcobalamin biosynthesis I (anaerobic)	Super-Pathways
PWY-4002	L-asparagine degradation II	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|ASPARAGINE-DEG
PWY-5838	superpathway of menaquinol-8 biosynthesis I	Super-Pathways
ASPARTATE-DEG1-PWY	L-aspartate degradation I	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|ASPARTATE-DEG
PWY-5920	superpathway of heme biosynthesis from glycine	Super-Pathways
MALATE-ASPARTATE-SHUTTLE-PWY	L-aspartate degradation II	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|ASPARTATE-DEG
PWY-6165	chorismate biosynthesis II (archaea)	Super-Pathways
LCYSDEG-PWY	L-cysteine degradation II	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|CYSTEINE-DEG
PWY-6350	archaetidylinositol biosynthesis	Super-Pathways
CYSTEINE-DEG-PWY	L-cysteine degradation I	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|CYSTEINE-DEG
PWY-6940	icosapentaenoate biosynthesis III (fungi)	Super-Pathways
PWY-5329	L-cysteine degradation III	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|CYSTEINE-DEG
PWY-7229	superpathway of adenosine nucleotides de novo biosynthesis I	Super-Pathways
GLUTDEG-PWY	L-glutamate degradation II	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|GLUTAMATE-DEG
PWY4FS-8	phosphatidylglycerol biosynthesis II (non-plastidic)	Super-Pathways
PWY-5087	L-glutamate degradation VI (to pyruvate)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|GLUTAMATE-DEG
ARO-PWY	chorismate biosynthesis I	Super-Pathways
PWY0-1305	L-glutamate degradation IX (via 4-aminobutanoate)	Super-Pathways
GALACT-GLUCUROCAT-PWY	superpathway of hexuronide and hexuronate degradation	Super-Pathways
GLUTAMATE-DEG1-PWY	L-glutamate degradation I	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|GLUTAMATE-DEG
MET-SAM-PWY	superpathway of S-adenosyl-L-methionine biosynthesis	Super-Pathways
PWY-4321	L-glutamate degradation IV	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|GLUTAMATE-DEG
P381-PWY	adenosylcobalamin biosynthesis II (aerobic)	Super-Pathways
PWY-5766	L-glutamate degradation X	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|GLUTAMATE-DEG
PRPP-PWY	superpathway of histidine, purine, and pyrimidine biosynthesis	Super-Pathways
GLUTAMINEFUM-PWY	L-glutamine degradation II	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|GLUTAMINE-DEG
PWY-4041	&gamma;-glutamyl cycle	Super-Pathways
GLUTAMINDEG-PWY	L-glutamine degradation I	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|GLUTAMINE-DEG
PWY-5114	UDP-sugars interconversion	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
GLYCGREAT-PWY	creatine biosynthesis	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|GLYCINE-DEG
PWY-5416	superpathway of diterpene resin acids biosynthesis	Super-Pathways
GLYCLEAV-PWY	glycine cleavage	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|GLYCINE-DEG
PWY-5850	superpathway of menaquinol-6 biosynthesis I	Super-Pathways
HISHP-PWY	L-histidine degradation VI	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|HISTIDINE-DEG
PWY-5991	superpathway of linamarin and lotaustralin biosynthesis	Super-Pathways
PWY-5030	L-histidine degradation III	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|HISTIDINE-DEG
PWY-6371	superpathway of inositol phosphate compounds	Super-Pathways
HISDEG-PWY	L-histidine degradation I	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|HISTIDINE-DEG
PWY-6555	superpathway of 1D-myo-inositol hexakisphosphate biosynthesis (plants)	Super-Pathways
PWY-5028	L-histidine degradation II	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|HISTIDINE-DEG
PWY-7124	ethylene biosynthesis V (engineered)	Super-Pathways
HISTDEG-PWY	L-histidine degradation IV	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|HISTIDINE-DEG
PWY-7592	arachidonate biosynthesis III (6-desaturase, mammals)	Super-Pathways
PWY-5031	L-histidine degradation V	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|HISTIDINE-DEG
PWY66-341	cholesterol biosynthesis I	Super-Pathways
HOMOCYSDEGR-PWY	L-cysteine biosynthesis III (from L-homocysteine)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|HOMOCYSTEINE-DEG
URSIN-PWY	ureide biosynthesis	Super-Pathways
ILEUDEG-PWY	L-isoleucine degradation I	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|ISOLEUCINE-DEG
BIOTIN-BIOSYNTHESIS-PWY	biotin biosynthesis I	Super-Pathways
PWY-5078	L-isoleucine degradation II	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|ISOLEUCINE-DEG
PWY-4242	pantothenate and coenzyme A biosynthesis III	Super-Pathways
PWY-7767	leucine degradation IV	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|LEUCINE-DEG
PWY-5134	bitter acids biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|TERPENOPHENOLICS-SYN
PWY-5076	L-leucine degradation III	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|LEUCINE-DEG
PWY-5424	superpathway of oleoresin turpentine biosynthesis	Super-Pathways
PWY-5075	L-leucine degradation II	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|LEUCINE-DEG
PWY-5634	superpathway of penicillin, cephalosporin and cephamycin biosynthesis	Super-Pathways
LEU-DEG2-PWY	L-leucine degradation I	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|LEUCINE-DEG
PWY-5862	superpathway of demethylmenaquinol-9 biosynthesis	Super-Pathways
PWY66-425	L-lysine degradation II (L-pipecolate pathway)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|LYSINE-DEG
PWY-6266	superpathway of flavones and derivatives biosynthesis	Super-Pathways
PWY-5298	L-lysine degradation VI	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|LYSINE-DEG
PWY-6569	chondroitin sulfate biosynthesis	Super-Pathways
PWY-5324	L-lysine degradation IX	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|LYSINE-DEG
PWY-7169	hyperxanthone E biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|XANTHONE-SYN
PWY0-461	L-lysine degradation I	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|LYSINE-DEG
PWY-7277	sphingolipid biosynthesis (mammals)	Super-Pathways
LYSINE-DEG1-PWY	L-lysine degradation XI (mammalian)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|LYSINE-DEG
PWY-7393	&beta;-carotene biosynthesis (engineered)	Super-Pathways
PWY-5283	L-lysine degradation V	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|LYSINE-DEG
PWY0-166	superpathway of pyrimidine deoxyribonucleotides de novo biosynthesis (E. coli)	Super-Pathways
PWY-5314	L-lysine degradation VIII	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|LYSINE-DEG
PWY66-5	superpathway of cholesterol biosynthesis	Super-Pathways
PWY-6328	L-lysine degradation X	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|LYSINE-DEG
COLANSYN-PWY	colanic acid building blocks biosynthesis	Super-Pathways
LYSDEGII-PWY	L-lysine degradation III	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|LYSINE-DEG
NAD-BIOSYNTHESIS-II	NAD salvage pathway III	Super-Pathways
PWY-5280	L-lysine degradation IV	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|LYSINE-DEG
PANTOSYN-PWY	pantothenate and coenzyme A biosynthesis I	Super-Pathways
PWY-5311	L-lysine degradation VII	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|LYSINE-DEG
PWY-2055	superpathway of pterocarpan biosynthesis (via daidzein)	Super-Pathways
PWY-5327	superpathway of L-lysine degradation	Super-Pathways
PWY-5178	toluene degradation IV (aerobic) (via catechol)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|TOLUENE-DEG
METHIONINE-DEG1-PWY	L-methionine degradation I (to L-homocysteine)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|METHIONINE-DEG
PWY-5312	superpathway of anthocyanin biosynthesis (from pelargonidin 3-O-glucoside)	Super-Pathways
PWY-701	L-methionine degradation II	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|METHIONINE-DEG
PWY-5433	superpathway of lipoxygenase	Super-Pathways
PWY-5328	superpathway of L-methionine salvage and degradation	Super-Pathways
PWY-5896	superpathway of menaquinol-10 biosynthesis	Super-Pathways
PWY-5082	L-methionine degradation III	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|METHIONINE-DEG
PWY-6125	superpathway of guanosine nucleotides de novo biosynthesis II	Super-Pathways
PHENYLALANINE-DEG1-PWY	L-phenylalanine degradation I (aerobic)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|PHENYLALANINE-DEG
PWY-6284	superpathway of unsaturated fatty acids biosynthesis (E. coli)	Super-Pathways
PWY-7158	L-phenylalanine degradation V	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|PHENYLALANINE-DEG
PWY-6404	superpathway of mycolyl-arabinogalactan-peptidoglycan complex biosynthesis	Super-Pathways
ANAPHENOXI-PWY	L-phenylalanine degradation II (anaerobic)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|PHENYLALANINE-DEG
PWY-6443	benzoate biosynthesis I (CoA-dependent, &beta;-oxidative)	Super-Pathways
PWY-6318	L-phenylalanine degradation IV (mammalian, via side chain)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|PHENYLALANINE-DEG
PWY-6886	1-butanol autotrophic biosynthesis (engineered)	Super-Pathways
PWY-5079	L-phenylalanine degradation III	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|PHENYLALANINE-DEG
PWY-6981	chitin biosynthesis	Super-Pathways
PROUT-PWY	L-proline degradation	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|PROLINE-DEG
PWY-7196	superpathway of pyrimidine ribonucleosides salvage	Super-Pathways
SERDEG-PWY	L-serine degradation	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|SERINE-DEG
PWY-7411	superpathway of phosphatidate biosynthesis (yeast)	Super-Pathways
THREOCAT-PWY	superpathway of L-threonine metabolism	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|THREONINE-DEG
PWY0-781	aspartate superpathway	Super-Pathways
PWY-5436	L-threonine degradation IV	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|THREONINE-DEG
SUCSYN-PWY	sucrose biosynthesis I (from photosynthesis)	Super-Pathways
THRDLCTCAT-PWY	L-threonine degradation III (to methylglyoxal)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|THREONINE-DEG
PHOSLIPSYN-PWY	superpathway of phospholipid biosynthesis I (bacteria)	Super-Pathways
PWY66-428	L-threonine degradation V	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|THREONINE-DEG
PWY-2541	plant sterol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
THREONINE-DEG2-PWY	L-threonine degradation II	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|THREONINE-DEG
PWY-5052	superpathway of gibberellin biosynthesis	Super-Pathways
PWY-5437	L-threonine degradation I	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|THREONINE-DEG
PWY-5181	toluene degradation III (aerobic) (via p-cresol)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|TOLUENE-DEG
PWY-3162	L-tryptophan degradation V (side chain pathway)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|TRYPTOPHAN-DEG
PWY-5478	superpathway of hydrolyzable tannin biosynthesis	Super-Pathways
PWY-5651	L-tryptophan degradation to 2-amino-3-carboxymuconate semialdehyde	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|TRYPTOPHAN-DEG
PWY-5899	superpathway of menaquinol-13 biosynthesis	Super-Pathways
PWY-6309	L-tryptophan degradation XI (mammalian, via kynurenine)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|TRYPTOPHAN-DEG
PWY-6142	gluconeogenesis II (Methanobacterium thermoautotrophicum)	Super-Pathways
TRPIAACAT-PWY	L-tryptophan degradation VII (via indole-3-pyruvate)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|TRYPTOPHAN-DEG
PWY-6897	thiamine salvage II	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Thiamine-Biosynthesis|Thiamin-Salvage
TRYPTOPHAN-DEGRADATION-1	L-tryptophan degradation III (eukaryotic)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|TRYPTOPHAN-DEG
PWY-7020	superpathway of butirocin biosynthesis	Super-Pathways
PWY-5081	L-tryptophan degradation VIII (to tryptophol)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|TRYPTOPHAN-DEG
PWY-7328	superpathway of UDP-glucose-derived O-antigen building blocks biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
PWY-6307	L-tryptophan degradation X (mammalian, via tryptamine)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|TRYPTOPHAN-DEG
PWY-841	superpathway of purine nucleotides de novo biosynthesis I	Super-Pathways
TRPCAT-PWY	L-tryptophan degradation I (via anthranilate)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|TRYPTOPHAN-DEG
THISYN-PWY	superpathway of thiamine diphosphate biosynthesis I	Super-Pathways
TRYPDEG-PWY	L-tryptophan degradation II (via pyruvate)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|TRYPTOPHAN-DEG
ERGOSTEROL-SYN-PWY	superpathway of ergosterol biosynthesis I	Super-Pathways
PWY-3181	L-tryptophan degradation VI (via tryptamine)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|TRYPTOPHAN-DEG
POLYAMINSYN3-PWY	superpathway of polyamine biosynthesis II	Super-Pathways
PWY-5655	L-tryptophan degradation IX	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|TRYPTOPHAN-DEG
PWY-5071	superpathway of rosmarinic acid biosynthesis	Super-Pathways
PWY-6505	L-tryptophan degradation XII (Geobacillus)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|TRYPTOPHAN-DEG
PWY-5184	toluene degradation VI (anaerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|TOLUENE-DEG
TRPKYNCAT-PWY	L-tryptophan degradation IV (via indole-3-lactate)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|TRYPTOPHAN-DEG
PWY-5823	superpathway of CDP-glucose-derived O-antigen building blocks biosynthesis	Super-Pathways
PWY-5151	L-tyrosine degradation II	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|TYROSINE-DEG
PWY-5918	superpathay of heme biosynthesis from glutamate	Super-Pathways
TYRFUMCAT-PWY	L-tyrosine degradation I	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|TYROSINE-DEG
PWY-6937	superpathway of testosterone and androsterone degradation	Super-Pathways
PWY3O-4108	L-tyrosine degradation III	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|TYROSINE-DEG
PWY-7050	icosapentaenoate biosynthesis IV (bacteria)	Super-Pathways
PWY-7514	L-tyrosine degradation IV (to 4-methylphenol)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|TYROSINE-DEG
PWY-7228	superpathway of guanosine nucleotides de novo biosynthesis I	Super-Pathways
VALDEG-PWY	L-valine degradation I	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|VALINE-DEG
PWY-7342	superpathway of nicotine biosynthesis	Super-Pathways
PWY-5057	L-valine degradation II	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|VALINE-DEG
PWY4FS-7	phosphatidylglycerol biosynthesis I (plastidic)	Super-Pathways
PWY-1781	&beta;-alanine degradation II	Degradation|Amino-Acid-Degradation|MISCELLANEOUS-DEG|Beta-Alanine-Degradation
FUC-RHAMCAT-PWY	superpathway of fucose and rhamnose degradation	Super-Pathways
BETA-ALA-DEGRADATION-I-PWY	&beta;-alanine degradation I	Degradation|Amino-Acid-Degradation|MISCELLANEOUS-DEG|Beta-Alanine-Degradation
LPSSYN-PWY	superpathway of lipopolysaccharide biosynthesis	Super-Pathways
PWY-6422	D-arginine degradation	Degradation|Amino-Acid-Degradation|MISCELLANEOUS-DEG|D-Amino-Acid-Degradation
PWY-5415	catechol degradation I (<I>meta</I>-cleavage pathway)	Super-Pathways
PWY-6196	D-serine metabolism	Super-Pathways
PWY-5845	superpathway of menaquinol-9 biosynthesis	Super-Pathways
PWY0-1535	D-serine degradation	Degradation|Amino-Acid-Degradation|MISCELLANEOUS-DEG|D-Amino-Acid-Degradation
PWY-6358	superpathway of D-myo-inositol (1,4,5)-trisphosphate metabolism	Super-Pathways
PWY-7515	trans-3-hydroxy-L-proline degradation	Degradation|Amino-Acid-Degradation|MISCELLANEOUS-DEG|Hydroxyprolines-Degradation
PWY-6544	superpathway of C28 brassinosteroid biosynthesis	Super-Pathways
PWY-5159	trans-4-hydroxy-L-proline degradation II	Degradation|Amino-Acid-Degradation|MISCELLANEOUS-DEG|Hydroxyprolines-Degradation
PWY-7110	superpathway of megalomicin A biosynthesis	Super-Pathways
HYDROXYPRODEG-PWY	trans-4-hydroxy-L-proline degradation I	Degradation|Amino-Acid-Degradation|MISCELLANEOUS-DEG|Hydroxyprolines-Degradation
PWY-7379	mRNA capping II	Super-Pathways
PWY-6991	(-)-camphor biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|MONOTERPENOID-SYN|Camphor-Biosynthesis
PWY-7575	superpathway of candicidin biosynthesis	Super-Pathways
PWY-5789	3-hydroxypropanoate/4-hydroxybutanate cycle	Degradation|C1-COMPOUNDS|CO2-Fixation|Autotrophic-CO2-Fixation
PWY66-3	cholesterol biosynthesis II (via 24,25-dihydrolanosterol)	Super-Pathways
CALVIN-PWY	Calvin-Benson-Bassham cycle	Energy-Metabolism|Photosynthesis
GLUCARGALACTSUPER-PWY	superpathway of D-glucarate and D-galactarate degradation	Super-Pathways
PWY-7784	reductive acetyl coenzyme A pathway II (autotrophic methanogens)	Degradation|C1-COMPOUNDS|CO2-Fixation|Autotrophic-CO2-Fixation
PWY-4221	pantothenate and coenzyme A biosynthesis II (plants)	Super-Pathways
PWY-5743	3-hydroxypropanoate cycle	Degradation|C1-COMPOUNDS|CO2-Fixation|Autotrophic-CO2-Fixation
PWY-5420	catechol degradation II (<I>meta</I>-cleavage pathway)	Super-Pathways
PWY-7024	superpathway of the 3-hydroxypropanoate cycle	Super-Pathways
PWY-5861	superpathway of demethylmenaquinol-8 biosynthesis	Super-Pathways
CODH-PWY	reductive acetyl coenzyme A pathway I (homoacetogenic bacteria)	Degradation|C1-COMPOUNDS|CO2-Fixation|Autotrophic-CO2-Fixation
PWY-6263	superpathway of menaquinol-8 biosynthesis II	Super-Pathways
PWY-7808	tetracycline resistance	Detoxification|Antibiotic-Resistance
PWY-6385	peptidoglycan biosynthesis III (mycobacteria)	Super-Pathways
PWY-7652	echinocandin B degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|ANTIBIOTIC-DEGRADATION
PWY-6565	superpathway of polyamine biosynthesis III	Super-Pathways
PWY-6848	rutin degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|FLAVONOID-DEG
PWY-7156	superpathway of ergosterol biosynthesis II	Super-Pathways
PWY-7134	rutin degradation (plants)	Degradation|SECONDARY-METABOLITE-DEGRADATION|FLAVONOID-DEG
PWY-7392	taxadiene biosynthesis (engineered)	Super-Pathways
PWY-7133	quercetin glucoside degradation (Allium)	Degradation|SECONDARY-METABOLITE-DEGRADATION|FLAVONOID-DEG
PWY-7609	superpathway of roquefortine, meleagrin and neoxaline biosynthesis	Super-Pathways
PWY-6996	daidzin and daidzein degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|FLAVONOID-DEG
PWY0-162	superpathway of pyrimidine ribonucleotides de novo biosynthesis	Super-Pathways
PWY-7445	luteolin triglucuronide degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|FLAVONOID-DEG
PWY66-409	superpathway of purine nucleotide salvage	Super-Pathways
PWY-6990	(+)-camphor biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|MONOTERPENOID-SYN|Camphor-Biosynthesis
CAROTENOID-PWY	superpathway of carotenoid biosynthesis	Super-Pathways
PWY-641	proanthocyanidins biosynthesis from flavanols	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|PROANTHOCYANIDIN-SYN
PWY-4821	UDP-D-xylose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
NOPALINEDEG-PWY	nopaline degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|N-CONTAINING-SECONDARY-CMPD-DEG
PWY-5310	superpathway of anthocyanin biosynthesis (from delphinidin 3-O-glucoside)	Super-Pathways
OCTOPINEDEG-PWY	octopine degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|N-CONTAINING-SECONDARY-CMPD-DEG
PWY-5431	aromatic compounds degradation via &beta;-ketoadipate	Super-Pathways
PWY-31	canavanine degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|N-CONTAINING-SECONDARY-CMPD-DEG
PWY-5864	superpathway of plastoquinol biosynthesis	Super-Pathways
PWY-6667	resveratrol degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|PHENYLPROPANOID-DERIVATIVE-DEG
PWY-6113	superpathway of mycolate biosynthesis	Super-Pathways
PWY-5319	coumarin metabolism (to melilotic acid)	Degradation|SECONDARY-METABOLITE-DEGRADATION|PHENYLPROPANOID-DERIVATIVE-DEG
PWY-6277	superpathway of 5-aminoimidazole ribonucleotide biosynthesis	Super-Pathways
PWY-7214	baicalein degradation (hydrogen peroxide detoxification)	Degradation|SECONDARY-METABOLITE-DEGRADATION|PHENYLPROPANOID-DERIVATIVE-DEG
PWY-6579	superpathway of guanine and guanosine salvage	Super-Pathways
PWY-116	coniferin metabolism	Degradation|SECONDARY-METABOLITE-DEGRADATION|PHENYLPROPANOID-DERIVATIVE-DEG
PWY-6977	superpathway of erythromycin biosynthesis	Super-Pathways
PWY-7113	furcatin degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|PHENYLPROPANOID-DERIVATIVE-DEG
PWY-7194	pyrimidine nucleobases salvage II	Super-Pathways
PWY-7615	pterocarpan phytoalexins modification (maackiain, medicarpin, pisatin, phaseollin)	Degradation|SECONDARY-METABOLITE-DEGRADATION|PHYTOALEXIN-DEG
PWY-7287	novobiocin biosynthesis	Super-Pathways
PWY-7635	kievitone detoxification	Degradation|SECONDARY-METABOLITE-DEGRADATION|PHYTOALEXIN-DEG
PWY-7715	superpathway of trichothecene biosynthesis	Super-Pathways
PWY-5316	nicotine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|PYRROLIDINE-ALKALOIDS
DENOVOPURINE2-PWY	superpathway of purine nucleotides de novo biosynthesis II	Super-Pathways
PWY-5752	piperine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|PYRROLIDINE-ALKALOIDS
PEPTIDOGLYCANSYN-PWY	peptidoglycan biosynthesis I (<I>meso</I>-diaminopimelate containing)	Super-Pathways
PWY-7022	paromamine biosynthesis II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis|Paromamine-Biosynthesis
PWY-5180	toluene degradation I (aerobic) (<I>via</I> o</I>-cresol)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|TOLUENE-DEG
PWY-7014	paromamine biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis|Paromamine-Biosynthesis
PWY-5898	superpathway of menaquinol-12 biosynthesis	Super-Pathways
PWY3DJ-12	ceramide de novo biosynthesis	Biosynthesis|Lipid-Biosynthesis|Sphingolipid-Biosynthesis
PWY-6141	archaetidylserine and archaetidylethanolamine biosynthesis	Super-Pathways
PWY-5129	sphingolipid biosynthesis (plants)	Biosynthesis|Lipid-Biosynthesis|Sphingolipid-Biosynthesis
PWY-6432	curcuminoid biosynthesis	Super-Pathways
SPHINGOLIPID-SYN-PWY	sphingolipid biosynthesis (yeast)	Biosynthesis|Lipid-Biosynthesis|Sphingolipid-Biosynthesis
PWY-6467	Kdo transfer to lipid IV<SUB>A</SUB> III (Chlamydia)	Super-Pathways
6-HYDROXYCINEOLE-DEGRADATION-PWY	1,8-cineole degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|TERPENOID-DEG
PWY-6895	superpathway of thiamine diphosphate biosynthesis II	Super-Pathways
PWY-6413	ginsenoside degradation III	Degradation|SECONDARY-METABOLITE-DEGRADATION|TERPENOID-DEG
PWY-7003	glycerol degradation to butanol	Super-Pathways
PWY-6670	citronellol degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|TERPENOID-DEG
PWY-7208	superpathway of pyrimidine nucleobases salvage	Super-Pathways
PWY-6672	cis-genanyl-CoA degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|TERPENOID-DEG
PWY-7323	superpathway of GDP-mannose-derived O-antigen building blocks biosynthesis	Super-Pathways
PWY-6678	geraniol and nerol degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|TERPENOID-DEG
PWY-7476	superpathway avenacin A biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-6806	carotenoid cleavage	Metabolic-Clusters
PWY0-881	superpathway of fatty acid biosynthesis I (E. coli)	Super-Pathways
PWY-7141	(3S)-linalool biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|MONOTERPENOID-SYN|Linalool-Biosynthesis
ARG+POLYAMINE-SYN	superpathway of arginine and polyamine biosynthesis	Super-Pathways
PWY-4841	UDP-&alpha;-D-glucuronate biosynthesis (from myo-inositol)	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
ENTBACSYN-PWY	enterobactin biosynthesis	Super-Pathways
PWY-6361	1D-myo-inositol hexakisphosphate biosynthesis I(from Ins(1,4,5)P3)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|SUGAR-DERIVS|Cyclitols-Biosynthesis|Phytate-Biosynthesis
PWY-3097	superpathway of isoflavonoids (via naringenin)	Super-Pathways
PWY-1741	indole-3-acetyl-ester conjugate biosynthesis	Activation-Inactivation-Interconversion|Inactivation
PWY-5183	superpathway of aerobic toluene degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|TOLUENE-DEG
PWY-5380	A series fagopyritols biosynthesis	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|CYCLITOLS-DEG
PWY-5724	superpathway of atrazine degradation	Super-Pathways
PWY-6554	1D-myo-inositol hexakisphosphate biosynthesis V (from Ins(1,3,4)P3)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|SUGAR-DERIVS|Cyclitols-Biosynthesis|Phytate-Biosynthesis
PWY-5910	superpathway of geranylgeranyldiphosphate biosynthesis I (via mevalonate)	Super-Pathways
PWY-5379	B series fagopyritols biosynthesis	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|CYCLITOLS-DEG
PWY-6338	superpathway of vanillin and vanillate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Vanillin-Degradation
PWY-6362	1D-myo-inositol hexakisphosphate biosynthesis II (mammalian)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|SUGAR-DERIVS|Cyclitols-Biosynthesis|Phytate-Biosynthesis
PWY-6503	superpathway of ergotamine biosynthesis	Super-Pathways
PWY-7650	echinocandin B biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis|LIPOPEPTIDE-BIOSYNTHESIS
PWY-6928	superpathway of cholesterol degradation I (cholesterol oxidase)	Super-Pathways
DARABITOLUTIL-PWY	D-arabitol degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ALCOHOLS-DEG
PWY-7218	photosynthetic 3-hydroxybutanoate biosynthesis (engineered)	Super-Pathways
GALACTITOLCAT-PWY	galactitol degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ALCOHOLS-DEG
PWY-7341	superpathway of hyoscyamine and scopolamine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|TROPANE-ALKALOIDS
HEXITOLDEGSUPER-PWY	superpathway of hexitol degradation (bacteria)	Super-Pathways
PWY4FS-5	superpathway of phosphatidylcholine biosynthesis	Super-Pathways
LARABITOLUTIL-PWY	xylitol degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ALCOHOLS-DEG
ARGORNPROST-PWY	arginine, ornithine and proline interconversion	Super-Pathways
MANNIDEG-PWY	mannitol degradation I	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ALCOHOLS-DEG
FOLSYN-PWY	superpathway of tetrahydrofolate biosynthesis and salvage	Super-Pathways
P562-PWY	<I>myo</I>-inositol degradation I	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ALCOHOLS-DEG
P281-PWY	3-phenylpropanoate degradation	Super-Pathways
PWY-4101	D-sorbitol degradation I	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ALCOHOLS-DEG
POLYISOPRENSYN-PWY	polyisoprenoid biosynthesis (E. coli)	Super-Pathways
PWY-6531	mannitol cycle	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ALCOHOLS-DEG
PWY-3502	superpathway of NAD biosynthesis in eukaryotes	Super-Pathways
PWY-7237	myo-, chiro- and scillo-inositol degradation	Super-Pathways
PWY-5265	peptidoglycan biosynthesis II (staphylococci)	Super-Pathways
PWY-7241	<I>myo</I>-inositol degradation II	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ALCOHOLS-DEG
PWY-5405	superpathway of betalain biosynthesis	Super-Pathways
PWY-7786	D-threitol degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ALCOHOLS-DEG
PWY-5529	superpathway of bacteriochlorophyll a biosynthesis	Super-Pathways
PWY-7787	L-threitol degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ALCOHOLS-DEG
PWY-5840	superpathway of menaquinol-7 biosynthesis	Super-Pathways
RIBITOLUTIL-PWY	ribitol degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ALCOHOLS-DEG
PWY-5942	trans-lycopene biosynthesis I (bacteria)	Super-Pathways
SORBDEG-PWY	D-sorbitol degradation II	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ALCOHOLS-DEG
PWY-6176	superpathway of Allium flavor precursors	Super-Pathways
PWY-7446	sulfoquinovose degradation I	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|Sulfoquinovose-Degradation
PWY-6947	superpathway of cholesterol degradation II (cholesterol dehydrogenase)	Super-Pathways
PWY-7722	sulfoquinovose degradation II	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|Sulfoquinovose-Degradation
PWY-7235	superpathway of ubiquinol-6 biosynthesis (eukaryotic)	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Ubiquinone-Biosynthesis
PWY-6497	D-galactarate degradation II	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ACIDS-DEG|D-Galactarate-Degradation
PWY-7373	superpathway of demethylmenaquinol-6 biosynthesis II	Super-Pathways
GALACTARDEG-PWY	D-galactarate degradation I	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ACIDS-DEG|D-Galactarate-Degradation
PWY-7537	superpathway of fumitremorgin biosynthesis	Super-Pathways
GALACTUROCAT-PWY	D-galacturonate degradation I	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ACIDS-DEG|D-Galacturonate-Degradation
UBISYN-PWY	superpathway of ubiquinol-8 biosynthesis (prokaryotic)	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Ubiquinone-Biosynthesis
PWY-6491	D-galacturonate degradation III	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ACIDS-DEG|D-Galacturonate-Degradation
PWY-13	superpathway of tetrahydroxyxanthone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|XANTHONE-SYN
PWY-6486	D-galacturonate degradation II	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ACIDS-DEG|D-Galacturonate-Degradation
PWY-5121	superpathway of geranylgeranyl diphosphate biosynthesis II (via MEP)	Super-Pathways
GLUCARDEG-PWY	D-glucarate degradation I	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ACIDS-DEG|D-Glucarate-Degradation
PWY-5417	catechol degradation III (ortho-cleavage pathway)	Super-Pathways
PWY-6499	D-glucarate degradation II	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ACIDS-DEG|D-Glucarate-Degradation
PWY-5860	superpathway of demethylmenaquinol-6 biosynthesis I	Super-Pathways
PWY-5525	D-glucuronate degradation I	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ACIDS-DEG|D-Glucuronate-Degradation
PWY-5993	superpathway of rifamycin B biosynthesis	Super-Pathways
GLUCUROCAT-PWY	superpathway of &beta;-D-glucuronide and D-glucuronate degradation	Super-Pathways
PWY-6234	superpathway of jasmonoyl-amino acid conjugates biosynthesis	Super-Pathways
PWY-7247	&beta;-D-glucuronide and D-glucuronate degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ACIDS-DEG|D-Glucuronate-Degradation
PWY-6374	vibriobactin biosynthesis	Super-Pathways
PWY-6501	D-glucuronate degradation II	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ACIDS-DEG|D-Glucuronate-Degradation
PWY-6564	heparan sulfate biosynthesis	Super-Pathways
PWY-5747	2-methylcitrate cycle II	Degradation|CARBOXYLATES-DEG|Propionate-Degradation|Methyl-Citrate-Cycle
PWY-6957	mandelate degradation to acetyl-CoA	Super-Pathways
PWY0-42	2-methylcitrate cycle I	Degradation|CARBOXYLATES-DEG|Propionate-Degradation|Methyl-Citrate-Cycle
PWY-7149	superpathway polymethylated quercetin/quercetagetin glucoside biosynthesis (Chrysosplenium)	Super-Pathways
PWY-5173	superpathway of acetyl-CoA biosynthesis	Super-Pathways
PWY-7245	superpathway NAD/NADP - NADH/NADPH interconversion (yeast)	Super-Pathways
PWY-5172	acetyl-CoA biosynthesis III (from citrate)	Energy-Metabolism|Acetyl-CoA-Biosynthesis
PWY-7596	superpathway of stearidonate biosynthesis (cyanobacteria)	Super-Pathways
PWY-6970	acetyl-CoA biosynthesis II (NADP-dependent pyruvate dehydrogenase)	Energy-Metabolism|Acetyl-CoA-Biosynthesis
PWY0-1415	superpathway of heme biosynthesis from uroporphyrinogen-III	Super-Pathways
PWY-6946	cholesterol degradation to androstenedione II (cholesterol dehydrogenase)	Degradation|Steroids-Degradation|Cholesterol-Degradation
PWY66-4	cholesterol biosynthesis III (via desmosterol)	Super-Pathways
PWY-7457	sulfite oxidation V (SoeABC)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfite-Oxidation
PWY-1782	superpathway of indole-3-acetate conjugate biosynthesis	Super-Pathways
AMMOXID-PWY	ammonia oxidation I (aerobic)	Energy-Metabolism|CHEMOAUTOTROPHIC-ENERGY-METABOLISM
PWY-4762	superpathway of choline biosynthesis	Super-Pathways
P282-PWY	nitrite oxidation	Energy-Metabolism|Electron-Transfer
PWY-5156	superpathway of fatty acid biosynthesis II (plant)	Super-Pathways
PWY-5285	sulfide oxidation III (persulfide dioxygenase)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfide-Oxidation
PWY-5430	meta cleavage pathway of aromatic compounds	Super-Pathways
PWY-7429	arsenite oxidation II (respiratory)	Energy-Metabolism|Electron-Transfer
PWY-5647	2-nitrobenzoate degradation I	Super-Pathways
P222-PWY	sulfide oxidation I (sulfide-quinone reductase)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfide-Oxidation
PWY-5863	superpathway of phylloquinol biosynthesis	Super-Pathways
P303-PWY	ammonia oxidation II (anaerobic)	Energy-Metabolism|CHEMOAUTOTROPHIC-ENERGY-METABOLISM
PWY-6270	isoprene biosynthesis I	Super-Pathways
PWY-5276	sulfite oxidation I (sulfite oxidoreductase)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfite-Oxidation
PWY-6945	cholesterol degradation to androstenedione I (cholesterol oxidase)	Degradation|Steroids-Degradation|Cholesterol-Degradation
PWY-7082	ammonia oxidation IV (autotrophic ammonia oxidizers)	Super-Pathways
PWY-5815	rubber biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|POLYTERPENOID-SYN
SULFUROX-PWY	sulfur oxidation I (aerobic)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfur-Oxidation
PWY-7565	aspyridone A biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
FESULFOX-PWY	sulfur oxidation II (Fe<sup>+3</sup>-dependent)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfur-Oxidation
PWY-12	pentaketide chromone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
P301-PWY	galena oxidation	Energy-Metabolism|CHEMOAUTOTROPHIC-ENERGY-METABOLISM
PWY-7513	flaviolin dimer and mompain biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY-4521	arsenite oxidation I (respiratory)	Energy-Metabolism|Electron-Transfer
PWY-7695	aurofusarin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY-5308	superpathway of sulfur metabolism (Desulfocapsa sulfoexigens)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Thiosulfate-Disproportionation
PWY-6312	barbaloin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY-3781	aerobic respiration I (cytochrome c)	Energy-Metabolism|Electron-Transfer
PWY-7631	arctigenin and isoarctigenin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|LIGNAN-SYN
PWY0-1353	succinate to cytochrome bd oxidase electron transfer	Energy-Metabolism|Electron-Transfer
PWY-7487	(+)-secoisolariciresinol diglucoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|LIGNAN-SYN
PWY-5358	tetrathionate reduction I (to thiosulfate)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Tetrathionate-Reduction
PWY-6820	diphyllin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|LIGNAN-SYN
PWY0-1544	proline to cytochrome bo oxidase electron transfer	Energy-Metabolism|Electron-Transfer
PWY-7672	fusaric acid biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY-7544	pyruvate to cytochrome bo oxidase electron transfer	Energy-Metabolism|Electron-Transfer
PWY-5466	matairesinol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|LIGNAN-SYN
PWY0-1567	NADH to cytochrome bo oxidase electron transfer II	Energy-Metabolism|Electron-Transfer
PWY-7749	(-)-4-desmethyl-epipodophyllotoxin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|LIGNAN-SYN
PWY0-1329	succinate to cytochrome bo oxidase electron transfer	Energy-Metabolism|Electron-Transfer
PWY-6824	justicidin B biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|LIGNAN-SYN
PWY0-1576	hydrogen to fumarate electron transfer	Energy-Metabolism|Electron-Transfer
PWY-7692	bikaverin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY0-1336	NADH to fumarate electron transfer	Energy-Metabolism|Electron-Transfer
PWY-5670	epoxysqualene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY0-1581	nitrate reduction IX (dissimilatory)	Degradation|Noncarbon-Nutrients|NITROGEN-DEG|Nitrate-Reduction
PWY-5469	sesamin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|LIGNAN-SYN
PWY-101	photosynthesis light reactions	Energy-Metabolism|Photosynthesis
PWY-7630	hinokinin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|LIGNAN-SYN
PWY0-1352	nitrate reduction VIII (dissimilatory)	Degradation|Noncarbon-Nutrients|NITROGEN-DEG|Nitrate-Reduction
PWY-6627	salinosporamide A biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY0-1356	formate to dimethyl sulfoxide electron transfer	Energy-Metabolism|Electron-Transfer
PWY-6174	mevalonate pathway II (archaea)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|ISOPRENOIDS|Isopentenyl-Diphosphate-Biosynthesis
PWY0-1565	D-lactate to cytochrome bo oxidase electron transport	Energy-Metabolism|Electron-Transfer
PWY-7139	sesaminol glucoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|LIGNAN-SYN
PWY0-1321	nitrate reduction III (dissimilatory)	Degradation|Noncarbon-Nutrients|NITROGEN-DEG|Nitrate-Reduction
PWY-5479	6-methoxypodophyllotoxin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|LIGNAN-SYN
PWY0-1573	nitrate reduction VIIIb (dissimilatory)	Degradation|Noncarbon-Nutrients|NITROGEN-DEG|Nitrate-Reduction
PWY-7563	bassianin and desmethylbassianin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY0-1335	NADH to cytochrome bo oxidase electron transfer I	Energy-Metabolism|Electron-Transfer
PWY-5152	leucodelphinidin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|LEUCOANTHOCYANIDIN-SYN
PWY0-1578	hydrogen to trimethylamine N-oxide electron transfer	Energy-Metabolism|Electron-Transfer
PWY1F-823	leucopelargonidin and leucocyanidin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|LEUCOANTHOCYANIDIN-SYN
PWY0-1348	NADH to dimethyl sulfoxide electron transfer	Energy-Metabolism|Electron-Transfer
PWY-7524	mevalonate pathway III (archaea)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|ISOPRENOIDS|Isopentenyl-Diphosphate-Biosynthesis
PWY0-1584	nitrate reduction X (periplasmic, dissimilatory)	Degradation|Noncarbon-Nutrients|NITROGEN-DEG|Nitrate-Reduction
PWY-4801	aloesone biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY-4302	aerobic respiration III (alternative oxidase pathway)	Energy-Metabolism|Electron-Transfer
PWY-7684	&alpha;-diglucosyldiacylglycerol biosynthesis	Biosynthesis|Lipid-Biosynthesis
PWY0-1355	formate to trimethylamine N-oxide electron transfer	Energy-Metabolism|Electron-Transfer
PWY-7738	polyacyltrehalose biosynthesis	Biosynthesis|Lipid-Biosynthesis
PWY-7279	aerobic respiration II (cytochrome c) (yeast)	Energy-Metabolism|Electron-Transfer
PWY-7740	&beta;-D-mannosyl phosphomycoketide biosynthesis	Biosynthesis|Lipid-Biosynthesis
PWY0-1561	glycerol-3-phosphate to cytochrome bo oxidase electron transfer	Energy-Metabolism|Electron-Transfer
PWY-6310	aloesone biosynthesis II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY-7545	pyruvate to cytochrome bd terminal oxidase electron transfer	Energy-Metabolism|Electron-Transfer
PWY-7742	phenolphthiocerol biosynthesis	Biosynthesis|Lipid-Biosynthesis
PWY0-1568	NADH to cytochrome bd oxidase electron transport II	Energy-Metabolism|Electron-Transfer
PWY-7743	dimycocerosyl triglycosyl phenolphthiocerol biosynthesis	Biosynthesis|Lipid-Biosynthesis
PWY0-1334	NADH to cytochrome bd oxidase electron transfer I	Energy-Metabolism|Electron-Transfer
PWY-7744	dimycocerosyl phthiocerol biosynthesis	Biosynthesis|Lipid-Biosynthesis
PWY0-1577	hydrogen to dimethyl sulfoxide electron transfer	Energy-Metabolism|Electron-Transfer
PWY-6316	aromatic polyketides biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY0-1347	NADH to trimethylamine N-oxide electron transfer	Energy-Metabolism|Electron-Transfer
PWY-7746	mycobacterial sulfolipid biosynthesis	Biosynthesis|Lipid-Biosynthesis
PWY0-1582	glycerol-3-phosphate to fumarate electron transfer	Energy-Metabolism|Electron-Transfer
PWY-782	glycolipid desaturation	Biosynthesis|Lipid-Biosynthesis
NPGLUCAT-PWY	Entner-Doudoroff pathway II (non-phosphorylative)	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|Entner-Duodoroff-Pathways
PWY0-1264	biotin-carboxyl carrier protein assembly	Biosynthesis|Lipid-Biosynthesis
ENTNER-DOUDOROFF-PWY	Entner-Doudoroff pathway I	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|Entner-Duodoroff-Pathways
PWY-7689	8-O-methylfusarubin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY-2221	Entner-Doudoroff pathway III (semi-phosphorylative)	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|Entner-Duodoroff-Pathways
PWYQT-4427	sulfoquinovosyl diacylglycerol biosynthesis	Biosynthesis|Lipid-Biosynthesis
GLYCOLYSIS	glycolysis I (from glucose 6-phosphate)	Energy-Metabolism|GLYCOLYSIS-VARIANTS
SOPHOROSYLOXYDOCOSANOATE-SYN-PWY	sophorolipid biosynthesis	Biosynthesis|Lipid-Biosynthesis
PWY-5484	glycolysis II (from fructose 6-phosphate)	Energy-Metabolism|GLYCOLYSIS-VARIANTS
TRIGLSYN-PWY	diacylglycerol and triacylglycerol biosynthesis	Biosynthesis|Lipid-Biosynthesis
ANAGLYCOLYSIS-PWY	glycolysis III (from glucose)	Energy-Metabolism|GLYCOLYSIS-VARIANTS
PWY-6438	phenylphenalenone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY-1042	glycolysis IV (plant cytosol)	Energy-Metabolism|GLYCOLYSIS-VARIANTS
PWY3DJ-11281	sphingomyelin metabolism	Biosynthesis|Lipid-Biosynthesis
P341-PWY	glycolysis V (Pyrococcus)	Energy-Metabolism|GLYCOLYSIS-VARIANTS
PWY-7745	p-HBAD biosynthesis	Biosynthesis|Lipid-Biosynthesis
PWY-6744	hydrogen production I	Energy-Metabolism|Hydrogen-Production
PWY-7741	phthiocerol biosynthesis	Biosynthesis|Lipid-Biosynthesis
PWY-6759	hydrogen production III	Energy-Metabolism|Hydrogen-Production
PWY-7561	tenellin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY-6780	hydrogen production VI	Energy-Metabolism|Hydrogen-Production
PWY-7666	galactolipid biosynthesis II	Biosynthesis|Lipid-Biosynthesis
PWY-6758	hydrogen production II	Energy-Metabolism|Hydrogen-Production
PWY-7420	monoacylglycerol metabolism (yeast)	Biosynthesis|Lipid-Biosynthesis
PWY-6772	hydrogen production V	Energy-Metabolism|Hydrogen-Production
PWY-6818	ornithine lipid biosynthesis	Biosynthesis|Lipid-Biosynthesis
PWY-7731	superpathway of photosynthetic hydrogen production	Super-Pathways
PWY-5393	raspberry ketone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY-6751	superpathway of hydrogen production	Super-Pathways
PWY-6804	diacylglycerol biosynthesis (PUFA enrichment in oilseed)	Biosynthesis|Lipid-Biosynthesis
PWY-6765	hydrogen production IV	Energy-Metabolism|Hydrogen-Production
PWY-6803	phosphatidylcholine acyl editing	Biosynthesis|Lipid-Biosynthesis
PWY-6785	hydrogen production VIII	Energy-Metabolism|Hydrogen-Production
PWY-6453	stigma estolide biosynthesis	Biosynthesis|Lipid-Biosynthesis
PWY66-367	ketogenesis	Energy-Metabolism|OTHER-ENERGY
PWY-6314	plumbagin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY-6118	glycerol-3-phosphate shuttle	Energy-Metabolism|OTHER-ENERGY
PWY-6129	dolichol and dolichyl phosphate biosynthesis	Biosynthesis|Lipid-Biosynthesis
PWY-5302	sulfur disproportionation II (aerobic)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfur-Disproportionation
PWY-5473	hydroxycinnamic acid serotonin amides biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN
PWY66-368	ketolysis	Energy-Metabolism|OTHER-ENERGY
PWY-5474	hydroxycinnamic acid tyramine amides biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN
P203-PWY	sulfur disproportionation I (anaerobic)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfur-Disproportionation
PWY-6442	spermidine hydroxycinnamic acid conjugates biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN
OXIDATIVEPENT-PWY	pentose phosphate pathway (oxidative branch) I	Energy-Metabolism|Pentose-Phosphate-Cycle
PWY-6448	hordatine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN
PWY-7796	pentose phosphate pathway (oxidative branch) II	Energy-Metabolism|Pentose-Phosphate-Cycle
PWY-7673	fusarin C biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
NONOXIPENT-PWY	pentose phosphate pathway (non-oxidative branch)	Energy-Metabolism|Pentose-Phosphate-Cycle
PWY-5885	wax esters biosynthesis II	Biosynthesis|Lipid-Biosynthesis
PENTOSE-P-PWY	pentose phosphate pathway	Super-Pathways
PWY-5884	wax esters biosynthesis I	Biosynthesis|Lipid-Biosynthesis
P21-PWY	pentose phosphate pathway (partial)	Energy-Metabolism|Pentose-Phosphate-Cycle
PWY-5148	acyl-CoA hydrolysis	Biosynthesis|Lipid-Biosynthesis
PWY-241	C4 photosynthetic carbon assimilation cycle, NADP-ME type	Energy-Metabolism|Photosynthesis
PWY-5143	long-chain fatty acid activation	Biosynthesis|Lipid-Biosynthesis
PWY-181	photorespiration	Energy-Metabolism|Photosynthesis
PWY-5142	acyl-ACP thioesterase pathway	Biosynthesis|Lipid-Biosynthesis
PWY-7117	C4 photosynthetic carbon assimilation cycle, PEPCK type	Energy-Metabolism|Photosynthesis
PWY-401	galactolipid biosynthesis I	Biosynthesis|Lipid-Biosynthesis
PWY-7115	C4 photosynthetic carbon assimilation cycle, NAD-ME type	Energy-Metabolism|Photosynthesis
PWY-321	cutin biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|PLANT-CELL-STRUCTURE|EPIDERMAL-STRUCTURE
PHOTOALL-PWY	oxygenic photosynthesis	Super-Pathways
PWY-282	cuticular wax biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|PLANT-CELL-STRUCTURE|EPIDERMAL-STRUCTURE
NAGLIPASYN-PWY	lipid IVA biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|Lipopolysaccharide-Biosynthesis
LIPA-CORESYN-PWY	Lipid A-core biosynthesis	Biosynthesis|Lipid-Biosynthesis
KDO-LIPASYN-PWY	(Kdo)2-lipid A biosynthesis I	Biosynthesis|Lipid-Biosynthesis
PWY-7588	ursodeoxycholate biosynthesis (bacteria)	Biosynthesis|Lipid-Biosynthesis|Sterol-Biosynthesis
PWY-5046	2-oxoisovalerate decarboxylation to isobutanoyl-CoA	Energy-Metabolism|Respiration
PWY-7455	allopregnanolone biosynthesis	Biosynthesis|Lipid-Biosynthesis|Sterol-Biosynthesis
PWY-5084	2-oxoglutarate decarboxylation to succinyl-CoA	Energy-Metabolism|Respiration
PWY-5666	&alpha;-solanine/&alpha;-chaconine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY-7254	TCA cycle VII (acetate-producers)	Energy-Metabolism|TCA-VARIANTS
PWY-7555	&alpha;-cyclopiazonate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
TCA	TCA cycle I (prokaryotic)	Energy-Metabolism|TCA-VARIANTS
PWY-5883	ephedrine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
P105-PWY	TCA cycle IV (2-oxoglutarate decarboxylase)	Energy-Metabolism|TCA-VARIANTS
PWY-7603	stephacidin A biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY-6969	TCA cycle V (2-oxoglutarate:ferredoxin oxidoreductase)	Energy-Metabolism|TCA-VARIANTS
PWY-7607	meleagrin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
REDCITCYC	TCA cycle VIII (helicobacter)	Energy-Metabolism|TCA-VARIANTS
PWY-7138	noscapine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY-5913	partial TCA cycle (obligate autotrophs)	Energy-Metabolism|TCA-VARIANTS
PWY-7612	chaetoglobosin A biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY66-398	TCA cycle III (animals)	Energy-Metabolism|TCA-VARIANTS
PWY-7526	fumitremorgin A biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
TCA-GLYOX-BYPASS	superpathway of glyoxylate bypass and TCA	Energy-Metabolism|TCA-VARIANTS
PWY-7705	4-methoxyviridicatin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY-5690	TCA cycle II (plants and fungi)	Energy-Metabolism|TCA-VARIANTS
PWY-5315	N-methyl-&Delta;<sup>1</sup>-pyrrolinium cation biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY-6389	(S)-acetoin biosynthesis	Energy-Metabolism|Fermentation|Acetoin-Biosynthesis
PWY-7549	asperlicin E biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY-5939	(R)-acetoin biosynthesis II	Energy-Metabolism|Fermentation|Acetoin-Biosynthesis
PWY-5748	&gamma;-coniciene and coniine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY-5938	(R)-acetoin biosynthesis I	Energy-Metabolism|Fermentation|Acetoin-Biosynthesis
PWY-7600	peramine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY3O-440	acetoin biosynthesis III	Energy-Metabolism|Fermentation|Acetoin-Biosynthesis
PWY-6027	capsiconiate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY-5676	acetyl-CoA fermentation to butanoate II	Super-Pathways
PWY-7605	roquefortine C biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
P125-PWY	superpathway of (R,R)-butanediol biosynthesis	Super-Pathways
PWY-7135	emetine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY-6391	meso-butanediol biosynthesis I	Energy-Metabolism|Fermentation|Butanediol-Biosynthesis
PWY-7525	fumitremorgin C biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY-6390	(S,S)-butanediol biosynthesis	Energy-Metabolism|Fermentation|Butanediol-Biosynthesis
PWY-7704	viridicatin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY-6396	superpathway of 2,3-butanediol biosynthesis	Super-Pathways
PWY-7543	5-N-acetylardeemin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY-5951	(R,R)-butanediol biosynthesis	Energy-Metabolism|Fermentation|Butanediol-Biosynthesis
PWY-5710	capsaicin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY-6392	meso-butanediol biosynthesis II	Energy-Metabolism|Fermentation|Butanediol-Biosynthesis
PWY-7558	&alpha;-cyclopiazonate detoxification	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
P142-PWY	pyruvate fermentation to acetate I	Super-Pathways
PWY-5958	acridone alkaloid biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY-6587	pyruvate fermentation to ethanol III	Energy-Metabolism|Fermentation|Pyruvate-Degradation
PWY-7604	notoamide C and D biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY-5100	pyruvate fermentation to acetate and lactate II	Super-Pathways
PWY-6923	ricinine degradation	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY-6594	superpathway of Clostridium acetobutylicum solventogenic fermentation	Super-Pathways
PWY-7608	neoxaline biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY-5482	pyruvate fermentation to acetate II	Energy-Metabolism|Fermentation|Pyruvate-Degradation
PWY-7265	lampranthin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY-6883	pyruvate fermentation to butanol II (engineered)	Energy-Metabolism|Fermentation|Pyruvate-Degradation
PWY-7660	tryptoquialanine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY-5486	pyruvate fermentation to ethanol II	Energy-Metabolism|Fermentation|Pyruvate-Degradation
PWY-7542	fumiquinazoline D biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY4LZ-257	superpathway of fermentation (Chlamydomonas reinhardtii)	Super-Pathways
PWY-7826	Amaryllidacea alkaloids biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN
PWY-5538	pyruvate fermentation to acetate VI	Super-Pathways
PWY-5978	kanosamine biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|SUGAR-DERIVS|Kanosamine-Biosynthesis
P108-PWY	pyruvate fermentation to propanoate I	Energy-Metabolism|Fermentation|Pyruvate-Degradation
PWY8J2-22	kanosamine biosynthesis II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|SUGAR-DERIVS|Kanosamine-Biosynthesis
PWY-6583	pyruvate fermentation to butanol I	Energy-Metabolism|Fermentation|Pyruvate-Degradation
PWY-5418	phenol degradation I (aerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Phenolic-Compounds-Degradation|Phenol-Degradation
PWY-5096	pyruvate fermentation to acetate and alanine	Super-Pathways
PWY-5270	morphine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|ISOQUINOLINE-ALKALOIDS
PWY-6590	superpathway of Clostridium acetobutylicum acidogenic fermentation	Super-Pathways
PWY-5472	bisbenzylisoquinoline alkaloid biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|ISOQUINOLINE-ALKALOIDS
PWY-5481	pyruvate fermentation to lactate	Energy-Metabolism|Fermentation|Pyruvate-Degradation
PWY-7363	papaverine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|ISOQUINOLINE-ALKALOIDS
PWY-6863	pyruvate fermentation to hexanol (engineered)	Energy-Metabolism|Fermentation|Pyruvate-Degradation
PWY-3901	berberine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|ISOQUINOLINE-ALKALOIDS
PWY-5485	pyruvate fermentation to acetate IV	Super-Pathways
PWY-5470	palmatine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|ISOQUINOLINE-ALKALOIDS
PWY-7351	pyruvate fermentation to opines	Energy-Metabolism|Fermentation|Pyruvate-Degradation
PWY-6337	dehydroscoulerine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|ISOQUINOLINE-ALKALOIDS
PWY-5537	pyruvate fermentation to acetate V	Super-Pathways
PWY-7522	(R)-canadine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|ISOQUINOLINE-ALKALOIDS
CENTFERM-PWY	pyruvate fermentation to butanoate	Energy-Metabolism|Fermentation|Pyruvate-Degradation
PWY-5287	sanguinarine and macarpine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|ISOQUINOLINE-ALKALOIDS
PWY-5768	pyruvate fermentation to acetate VIII	Energy-Metabolism|Fermentation|Pyruvate-Degradation
PHENOLDEG-PWY	phenol degradation II (anaerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Phenolic-Compounds-Degradation|Phenol-Degradation
P41-PWY	pyruvate fermentation to acetate and lactate I	Super-Pathways
PWY-7507	chelerythrine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|ISOQUINOLINE-ALKALOIDS
PWY-6588	pyruvate fermentation to acetone	Energy-Metabolism|Fermentation|Pyruvate-Degradation
PWY-5876	magnoflorine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|ISOQUINOLINE-ALKALOIDS
PWY-5480	pyruvate fermentation to ethanol I	Energy-Metabolism|Fermentation|Pyruvate-Degradation
PWY-7155	7-dehydroporiferasterol biosynthesis	Biosynthesis|Lipid-Biosynthesis|Sterol-Biosynthesis
PWY-6604	superpathway of Clostridium acetobutylicum acidogenic and solventogenic fermentation	Super-Pathways
PWY-6132	lanosterol biosynthesis	Biosynthesis|Lipid-Biosynthesis|Sterol-Biosynthesis
PWY-5483	pyruvate fermentation to acetate III	Super-Pathways
PWY-6074	zymosterol biosynthesis	Biosynthesis|Lipid-Biosynthesis|Sterol-Biosynthesis
PWY-7111	pyruvate fermentation to isobutanol (engineered)	Energy-Metabolism|Fermentation|Pyruvate-Degradation
PWY-6061	bile acid biosynthesis, neutral pathway	Biosynthesis|Lipid-Biosynthesis|Sterol-Biosynthesis
PWY-5494	pyruvate fermentation to propanoate II (acrylate pathway)	Energy-Metabolism|Fermentation|Pyruvate-Degradation
PWY-6036	cardenolide glucosides biosynthesis	Biosynthesis|Lipid-Biosynthesis|Sterol-Biosynthesis
PWY-5600	pyruvate fermentation to acetate VII	Super-Pathways
PWY-6032	cardenolide biosynthesis	Biosynthesis|Lipid-Biosynthesis|Sterol-Biosynthesis
PWY-7164	chlorophyll a degradation III	Degradation|COFACTOR-DEGRADATION|Heme-Degradation|Chlorophyll-A-Degradation
PWY-7404	ceramide phosphoethanolamine biosynthesis	Biosynthesis|Lipid-Biosynthesis|Sphingolipid-Biosynthesis
PWY-6927	chlorophyll a degradation II	Degradation|COFACTOR-DEGRADATION|Heme-Degradation|Chlorophyll-A-Degradation
PWY-5366	palmitoleate biosynthesis II (plants and bacteria)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|Palmitoleate-Biosynthesis
PWY-5098	chlorophyll a degradation I	Degradation|COFACTOR-DEGRADATION|Heme-Degradation|Chlorophyll-A-Degradation
PWY-6282	palmitoleate biosynthesis I (from (5Z)-dodec-5-enoate)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|Palmitoleate-Biosynthesis
PWY-5519	D-arabinose degradation III	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|D-Arabinose-Degradation
PWY-7589	palmitoleate biosynthesis III (cyanobacteria)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|Palmitoleate-Biosynthesis
DARABCATK12-PWY	D-arabinose degradation I	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|D-Arabinose-Degradation
PWY-5995	linoleate biosynthesis I (plants)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|Linoleate-Biosynthesis
DARABCAT-PWY	D-arabinose degradation II	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|D-Arabinose-Degradation
PWY-6001	linoleate biosynthesis II (animals)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|Linoleate-Biosynthesis
LACTOSECAT-PWY	lactose and galactose degradation I	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|LACTOSE-DEG
PWY-7593	linoleate biosynthesis III (cyanobacteria)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|Linoleate-Biosynthesis
PWY-6693	D-galactose degradation IV	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|GALACTOSE-DEGRADATION
PWY-5971	palmitate biosynthesis II (bacteria and plants)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Palmitate-Biosynthesis
GALDEG-PWY	D-galactose degradation II	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|GALACTOSE-DEGRADATION
PWY-5994	palmitate biosynthesis I (animals and fungi)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Palmitate-Biosynthesis
PWY-6317	D-galactose degradation I (Leloir pathway)	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|GALACTOSE-DEGRADATION
PWY-5679	clavulanate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY66-422	D-galactose degradation V (Leloir pathway)	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|GALACTOSE-DEGRADATION
PWY-6346	staurosporine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-3821	D-galactose degradation III	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
PWY-7569	arginomycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-7568	L-galactose degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|GALACTOSE-DEGRADATION
PWY-7716	penicillin G and penicillin V biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-5515	L-arabinose degradation II	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|L-Arabinose-Degradation
PWY-5770	phenazine-1-carboxylate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
ARABCAT-PWY	L-arabinose degradation I	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|L-Arabinose-Degradation
PWY-6679	jadomycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-7295	L-arabinose degradation IV	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|L-Arabinose-Degradation
PWY-7045	mithramycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-5517	L-arabinose degradation III	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|L-Arabinose-Degradation
PWY-7626	bacilysin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-6713	L-rhamnose degradation II	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|L-rhamnose-Degradation
PWY-7734	quinoxaline-2-carboxylate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
RHAMCAT-PWY	L-rhamnose degradation I	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|L-rhamnose-Degradation
PWY-5887	albaflavenone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-6714	L-rhamnose degradation III	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|L-rhamnose-Degradation
PWY-6721	sangivamycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
LACTOSEUTIL-PWY	lactose degradation II	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|LACTOSE-DEG
PWY-7352	daunorubicin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
BGALACT-PWY	lactose degradation III	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|LACTOSE-DEG
PWY-7659	viridicatumtoxin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY-3801	sucrose degradation II (sucrose synthase)	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|SUCROSE-DEG
PWY-7769	phosalacine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-7345	superpathway of anaerobic sucrose degradation	Super-Pathways
PWY-5940	streptomycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
SUCUTIL-PWY	sucrose degradation I (sucrose phosphotransferase)	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|SUCROSE-DEG
PWY-6919	neopentalenoketolactone and pentalenate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-621	sucrose degradation III (sucrose invertase)	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|SUCROSE-DEG
PWY-7419	FR-900098 and FR-33289 antibiotics biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
SUCROSEUTIL2-PWY	sucrose degradation VII (sucrose 3-dehydrogenase)	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|SUCROSE-DEG
PWY-7670	fusaridione A biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY-5384	sucrose degradation IV (sucrose phosphorylase)	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|SUCROSE-DEG
PWY-5631	deacetylcephalosporin C biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY66-373	sucrose degradation V (sucrose &alpha;-glucosidase)	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|SUCROSE-DEG
PWY-6003	gramicidin S biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-2723	trehalose degradation V	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|Trehalose-Degradation
PWY-7015	ribostamycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
TREDEGLOW-PWY	trehalose degradation I (low osmolarity)	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|Trehalose-Degradation
PWY-7485	tetracenomycin C biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-2722	trehalose degradation IV	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|Trehalose-Degradation
PWY-7693	guadinomine B biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY0-1466	trehalose degradation VI (periplasmic)	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|Trehalose-Degradation
PWY-7812	tetracycline and oxytetracycline biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-2721	trehalose degradation III	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|Trehalose-Degradation
PWY-6345	K-252 biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY0-1182	trehalose degradation II (trehalase)	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|Trehalose-Degradation
PWY-7019	butirosin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-6760	xylose degradation III	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|Xylose-Degradation
PWY-7564	bacimethrin and bacimethrin pyrophosphate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-5516	xylose degradation II	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|Xylose-Degradation
PWY-7706	dapdiamides biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
XYLCAT-PWY	xylose degradation I	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|Xylose-Degradation
PWY-5757	fosfomycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-7294	xylose degradation IV	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation|Xylose-Degradation
PWY-6666	pyocyanin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-5797	indole-3-acetate degradation VI	Degradation|HORMONE-DEG|PLANT-HORMONE-DEG|AUXINS-DEGRADATION
PWY-7025	gentamicin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-1961	indole-3-acetate degradation I	Degradation|HORMONE-DEG|PLANT-HORMONE-DEG|AUXINS-DEGRADATION
PWY-7733	3-hydroxyquinaldate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-2021	indole-3-acetate degradation IV	Degradation|HORMONE-DEG|PLANT-HORMONE-DEG|AUXINS-DEGRADATION
PWY-5818	validamycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-5788	indole-3-acetate degradation V	Degradation|HORMONE-DEG|PLANT-HORMONE-DEG|AUXINS-DEGRADATION
PWY-6720	toyocamycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-1981	indole-3-acetate degradation III	Degradation|HORMONE-DEG|PLANT-HORMONE-DEG|AUXINS-DEGRADATION
PWY-7655	dechlorogriseofulvin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY-5784	indole-3-acetate conjugate biosynthesis II	Degradation|HORMONE-DEG|PLANT-HORMONE-DEG|AUXINS-DEGRADATION
PWY-7737	thiocoraline biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-5811	indole-3-acetate degradation VII	Degradation|HORMONE-DEG|PLANT-HORMONE-DEG|AUXINS-DEGRADATION
PWY-5930	terpentecin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-1962	indole-3-acetate degradation II	Degradation|HORMONE-DEG|PLANT-HORMONE-DEG|AUXINS-DEGRADATION
PWY-6915	pentalenolactone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-7640	neophaseic acid biosynthesis	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|Abscisic-Acid-Biosynthesis|Abscisic-Acid-Derivative-Biosynthesis
PWY-7355	doxorubicin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-5272	abscisic acid glucose ester metabolism	Activation-Inactivation-Interconversion|Inactivation
PWY-7669	equisetin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY-5271	phaseic acid biosynthesis	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|Abscisic-Acid-Biosynthesis|Abscisic-Acid-Derivative-Biosynthesis
PWY-7797	nocardicin A biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-7642	7-hydroxyabscisate biosynthesis	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|Abscisic-Acid-Biosynthesis|Abscisic-Acid-Derivative-Biosynthesis
PWY-5630	penicillin K biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-2841	cytokinins degradation	Metabolic-Clusters
PWY-7483	elloramycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-6477	gibberellin inactivation II (methylation)	Metabolic-Clusters
PWY-7690	holomycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-102	gibberellin inactivation I (2&beta;-hydroxylation)	Metabolic-Clusters
PWY-7811	6-methylpretetramide biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-6494	gibberellin inactivation III (epoxidation)	Metabolic-Clusters
PWY-5633	cephamycin C biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-7181	pyrimidine deoxyribonucleosides degradation	Degradation|NUCLEO-DEG|Pyrimidine-Degradation|Pyrimidine-Deoxyribonucleosides-Deg
PWY-6324	rebeccamycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-7206	pyrimidine deoxyribonucleotides dephosphorylation	Degradation|NUCLEO-DEG|Pyrimidine-Degradation|Pyrimidine-Deoxyribonucleosides-Deg
PWY-7018	paromomycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-6220	jasmonoyl-amino acid conjugates biosynthesis I	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|Jasmonates-Biosynthesis
PWY-7547	prodigiosin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-6430	thymine degradation	Degradation|NUCLEO-DEG|Pyrimidine-Degradation|Pyrimidine-Base-Degradation
PWY-7702	sch210971 and sch210972 biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY-7209	superpathway of pyrimidine ribonucleosides degradation	Super-Pathways
PWY1A0-6325	actinorhodin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-6233	jasmonoyl-amino acid conjugates biosynthesis II	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|Jasmonates-Biosynthesis
PWY-5737	(5R)-carbapenem carboxylate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-7195	pyrimidine ribonucleosides salvage III	Degradation|NUCLEO-DEG|Pyrimidine-Degradation|Pyrimidine-Ribonucleosides-Degradation
PWY-6511	3-methylarginine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-6556	pyrimidine ribonucleosides salvage II	Degradation|NUCLEO-DEG|Pyrimidine-Degradation|Pyrimidine-Ribonucleosides-Degradation
PWY-7570	blasticidin S biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY0-1295	pyrimidine ribonucleosides degradation	Degradation|NUCLEO-DEG|Pyrimidine-Degradation|Pyrimidine-Ribonucleosides-Degradation
PWY-7718	actinomycin D biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
SALVADEHYPOX-PWY	adenosine nucleotides degradation II	Degradation|NUCLEO-DEG|Purine-Degradation|Adenosine-Nucleotides-Degradation
PWY-5776	2-hydroxyphenazine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-6617	adenosine nucleotides degradation III	Degradation|NUCLEO-DEG|Purine-Degradation|Adenosine-Nucleotides-Degradation
PWY-6682	dehydrophos biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-6596	adenosine nucleotides degradation I	Degradation|NUCLEO-DEG|Purine-Degradation|Adenosine-Nucleotides-Degradation
PWY-7653	griseofulvin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY-6595	superpathway of guanosine nucleotides degradation (plants)	Super-Pathways
PWY-7735	echinomycin and triostin A biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-6608	guanosine nucleotides degradation III	Degradation|NUCLEO-DEG|Purine-Degradation|Guanosine-Nucleotides-Degradation
PWY-5929	puromycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-6607	guanosine nucleotides degradation I	Degradation|NUCLEO-DEG|Purine-Degradation|Guanosine-Nucleotides-Degradation
PWY-6722	candicidin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-6606	guanosine nucleotides degradation II	Degradation|NUCLEO-DEG|Purine-Degradation|Guanosine-Nucleotides-Degradation
PWY-7354	aclacinomycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-6235	hydroxyjasmonate sulfate biosynthesis	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|Jasmonates-Biosynthesis
PWY-7665	aureobasidin A biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-6297	tuberonate glucoside biosynthesis	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|Jasmonates-Biosynthesis
PWY-7770	indolmycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY0-501	lipoate biosynthesis and incorporation I	Biosynthesis|Cofactor-Biosynthesis|Lipoate-Biosynthesis
PWY-5629	isopenicillin N biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-6987	lipoate biosynthesis and incorporation III (Bacillus)	Biosynthesis|Cofactor-Biosynthesis|Lipoate-Biosynthesis
PWY-5984	rifamycin B biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY0-522	lipoate salvage I	Biosynthesis|Cofactor-Biosynthesis|Lipoate-Biosynthesis
PWY-6955	lincomycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-7382	lipoate biosynthesis and incorporation (yeast)	Biosynthesis|Cofactor-Biosynthesis|Lipoate-Biosynthesis
PWY-7441	polymyxin A biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY0-1275	lipoate biosynthesis and incorporation II	Biosynthesis|Cofactor-Biosynthesis|Lipoate-Biosynthesis
PWY-7671	saframycin A biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-6395	superpathway of seleno-compound metabolism	Super-Pathways
PWY-7810	chlorotetracycline biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-6984	lipoate salvage II	Biosynthesis|Cofactor-Biosynthesis|Lipoate-Biosynthesis
PWY-5632	cephalosporin C biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-7274	D-cycloserine biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|Other-Amino-Acid-Biosynthesis
PWY-6322	phosphinothricin tripeptide biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-6116	mannosylfructose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|Oligosaccharides-Biosynthesis
PWY-7016	neomycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-5337	stachyose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|Oligosaccharides-Biosynthesis
PWY-7510	rhizocticin A and B biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-6524	lychnose and isolychnose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|Oligosaccharides-Biosynthesis
PWY-7694	zwittermicin A biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-5025	indole-3-acetate biosynthesis IV (bacteria)	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|Auxin-Biosynthesis
PWY-7821	tunicamycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis
PWY-5342	ajugose biosynthesis I (galactinol-dependent)	Biosynthesis|Carbohydrates-Biosynthesis|Oligosaccharides-Biosynthesis
PWY-3161	indole-3-acetate biosynthesis III (bacteria)	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|Auxin-Biosynthesis
PWY-6525	stellariose and mediose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|Oligosaccharides-Biosynthesis
PWY-5343	ajugose biosynthesis II (galactinol-independent)	Biosynthesis|Carbohydrates-Biosynthesis|Oligosaccharides-Biosynthesis
PWY-5905	hypusine biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|Amino-Acids-Modification
PWY-581	indole-3-acetate biosynthesis II	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|Auxin-Biosynthesis
PWY-6482	diphthamide biosynthesis (archaea)	Biosynthesis|Amino-Acid-Biosynthesis|Amino-Acids-Modification
PWY-7546	diphthamide biosynthesis (eukaryotes)	Biosynthesis|Amino-Acid-Biosynthesis|Amino-Acids-Modification
ALANINE-SYN2-PWY	L-alanine biosynthesis II	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|ALANINE-SYN
PWY-5026	indole-3-acetate biosynthesis V (bacteria and fungi)	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|Auxin-Biosynthesis
PWY0-1061	superpathway of L-alanine biosynthesis	Super-Pathways
ALANINE-VALINESYN-PWY	L-alanine biosynthesis I	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|ALANINE-SYN
PWY0-1021	L-alanine biosynthesis III	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|ALANINE-SYN
PWYDQC-4	indole-3-acetate biosynthesis I	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|Auxin-Biosynthesis
ARGSYN-PWY	L-arginine biosynthesis I (via L-ornithine)	Super-Pathways
PWY-7701	4-hydroxy-4-methyl-L-glutamate biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|Other-Amino-Acid-Biosynthesis
PWY-5154	L-arginine biosynthesis III (via N-acetyl-L-citrulline)	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|ARGININE-SYN
PWY-6154	autoinducer AI-2 biosynthesis II (Vibrio)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Autoinducer-Biosynthesis
ARGSYNBSUB-PWY	L-arginine biosynthesis II (acetyl cycle)	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|ARGININE-SYN
PWY-7649	3-hydroxy-L-homotyrosine biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|Other-Amino-Acid-Biosynthesis
PWY-7400	L-arginine biosynthesis IV (archaebacteria)	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|ARGININE-SYN
PWY-6153	autoinducer AI-2 biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Autoinducer-Biosynthesis
ASPARAGINESYN-PWY	L-asparagine biosynthesis II	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|ASPARAGINE-SYN
PWY-7648	4-methyl-proline biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|Other-Amino-Acid-Biosynthesis
PWY490-4	L-asparagine biosynthesis III (tRNA-dependent)	Biosynthesis|Aminoacyl-tRNAs-Charging
PWY-7621	autoinducer CAI-1 biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Autoinducer-Biosynthesis
ASPARAGINE-BIOSYNTHESIS	L-asparagine biosynthesis I	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|ASPARAGINE-SYN
PWY-7275	L-homophenylalanine biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|Other-Amino-Acid-Biosynthesis
PWY0-1325	superpathway of L-asparagine biosynthesis	Super-Pathways
PWY-6157	autoinducer AI-1 biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Autoinducer-Biosynthesis
ASPARTATESYN-PWY	L-aspartate biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|ASPARTATE-SYN
PWY-6936	seleno-amino acid biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|Other-Amino-Acid-Biosynthesis
PWY-6292	superpathway of L-cysteine biosynthesis (mammalian)	Super-Pathways
PWY-7405	aurachin RE biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis|Aurachin-Biosynthesis
PWY-7289	L-cysteine biosynthesis V	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|CYSTEINE-SYN
PWY-5957	L-nicotianamine biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|Other-Amino-Acid-Biosynthesis
PWY-6308	L-cysteine biosynthesis II (tRNA-dependent)	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|CYSTEINE-SYN
PWY-7407	aurachin A, B, C and D biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis|Aurachin-Biosynthesis
CYSTSYN-PWY	L-cysteine biosynthesis I	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|CYSTEINE-SYN
PWY-5344	L-homocysteine biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|Other-Amino-Acid-Biosynthesis
PWY-6293	L-cysteine biosynthesis IV (fungi)	Super-Pathways
PWY-1186	L-homomethionine biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|Other-Amino-Acid-Biosynthesis
PWY-801	L-homocysteine and L-cysteine interconversion	Activation-Inactivation-Interconversion|Interconversion
HOMOSERSYN-PWY	L-homoserine biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|Other-Amino-Acid-Biosynthesis
GLUTSYN-PWY	L-glutamate biosynthesis I	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|GLUTAMATE-SYN
PWY-5170	melamine degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|s-Triazine-Degradation
PWY-5505	L-glutamate and L-glutamine biosynthesis	Super-Pathways
GLUTAMATE-SYN2-PWY	L-glutamate biosynthesis II	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|GLUTAMATE-SYN
PWY-5171	N-cyclopropylmelamine degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|s-Triazine-Degradation
SHIKIMATEDEG-PWY	shikimate degradation I	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Shikimate-Degradation
PWY-4341	L-glutamate biosynthesis V	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|GLUTAMATE-SYN
GLUGLNSYN-PWY	L-glutamate biosynthesis IV	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|GLUTAMATE-SYN
PWY-5726	deethylsimazine degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|s-Triazine-Degradation
GLUTSYNIII-PWY	L-glutamate biosynthesis III	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|GLUTAMATE-SYN
GLNSYN-PWY	L-glutamine biosynthesis I	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|GLUTAMINE-SYN
P141-PWY	atrazine degradation I (aerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|s-Triazine-Degradation|Atrazine-Degradation
PWY-6419	shikimate degradation II	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Shikimate-Degradation
PWY-6549	L-glutamine biosynthesis III	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|GLUTAMINE-SYN
GLYSYN-ALA-PWY	glycine biosynthesis III	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|GLYCINE-SYN
PWY-5731	atrazine degradation III	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|s-Triazine-Degradation|Atrazine-Degradation
GLYCINE-SYN2-PWY	glycine biosynthesis II	Super-Pathways
GLYSYN-THR-PWY	glycine biosynthesis IV	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|GLYCINE-SYN
PWY-5727	atrazine degradation II	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|s-Triazine-Degradation|Atrazine-Degradation
GLYSYN-PWY	glycine biosynthesis I	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|GLYCINE-SYN
PWY-5029	imidazole-lactate degradation	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|HISTIDINE-SYN
HISTSYN-PWY	L-histidine biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|HISTIDINE-SYN
PWY-3001	superpathway of L-isoleucine biosynthesis I	Super-Pathways
PWY-5104	L-isoleucine biosynthesis IV	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|ISOLEUCINE-SYN
ILEUSYN-PWY	L-isoleucine biosynthesis I (from threonine)	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|ISOLEUCINE-SYN
PWY-5103	L-isoleucine biosynthesis III	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|ISOLEUCINE-SYN
PWY-5101	L-isoleucine biosynthesis II	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|ISOLEUCINE-SYN
PWY-5108	L-isoleucine biosynthesis V	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|ISOLEUCINE-SYN
LEUSYN-PWY	L-leucine biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|LEUCINE-SYN
PWY-2941	L-lysine biosynthesis II	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|LYSINE-SYN
PWY-5097	L-lysine biosynthesis VI	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|LYSINE-SYN
LYSINE-AMINOAD-PWY	L-lysine biosynthesis IV	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|LYSINE-SYN
PWY-3081	L-lysine biosynthesis V	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|LYSINE-SYN
DAPLYSINESYN-PWY	L-lysine biosynthesis I	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|LYSINE-SYN
PWY-2942	L-lysine biosynthesis III	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|LYSINE-SYN
PWY-6628	superpathway of L-phenylalanine biosynthesis	Super-Pathways
PWY-3462	L-phenylalanine biosynthesis II	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|PHENYLALANINE-SYN
PHESYN	L-phenylalanine biosynthesis I	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|PHENYLALANINE-SYN
PWY-7432	L-phenylalanine biosynthesis III (cytosolic, plants)	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|PHENYLALANINE-SYN
PWY-4981	L-proline biosynthesis II (from arginine)	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|PROLINE-SYN
PWY-4281	L-proline biosynthesis IV	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|PROLINE-SYN
PWY-3341	L-proline biosynthesis III	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|PROLINE-SYN
PROSYN-PWY	L-proline biosynthesis I	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|PROLINE-SYN
PWY-6994	L-pyrrolysine biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|Pyrrolysine-Biosynthesis
SERSYN-PWY	L-serine biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|SERINE-BIOSYNTHESIS
PWY0-901	L-selenocysteine biosynthesis I (bacteria)	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|Selenocysteine-Biosynthesis
PWY-6281	L-selenocysteine biosynthesis II (archaea and eukaryotes)	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|Selenocysteine-Biosynthesis
THRESYN-PWY	superpathway of L-threonine biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|THREONINE-BIOSYNTHESIS
HOMOSER-THRESYN-PWY	L-threonine biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|THREONINE-BIOSYNTHESIS
TRPSYN-PWY	L-tryptophan biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|TRYPTOPHAN-BIOSYNTHESIS
PWY-6629	superpathway of L-tryptophan biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|TRYPTOPHAN-BIOSYNTHESIS
PWY-6134	L-tyrosine biosynthesis IV	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|TYROSINE-SYN
PWY-6120	L-tyrosine biosynthesis III	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|TYROSINE-SYN
TYRSYN	L-tyrosine biosynthesis I	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|TYROSINE-SYN
PWY-3461	L-tyrosine biosynthesis II	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|TYROSINE-SYN
PWY-6630	superpathway of L-tyrosine biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|TYROSINE-SYN
VALSYN-PWY	L-valine biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|VALINE-BIOSYNTHESIS
PWY-5345	superpathway of L-methionine biosynthesis (by sulfhydrylation)	Super-Pathways
METSYN-PWY	L-homoserine and L-methionine biosynthesis	Super-Pathways
PWY-702	L-methionine biosynthesis II	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|METHIONINE-SYN|Methionine-De-novo-Biosynthesis
HSERMETANA-PWY	L-methionine biosynthesis III	Super-Pathways
PWY-6598	sciadonate biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis
PWY-5347	superpathway of L-methionine biosynthesis (transsulfuration)	Super-Pathways
HOMOSER-METSYN-PWY	L-methionine biosynthesis I	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|METHIONINE-SYN|Methionine-De-novo-Biosynthesis
PWY-7152	pinolenate and coniferonate biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis
PWY-6151	S-adenosyl-L-methionine cycle I	Super-Pathways
PWY-5041	S-adenosyl-L-methionine cycle II	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|METHIONINE-SYN|Methionine-Salvage|S-adenosyl-L-methionine-cycle
PWY-7595	stearidonate biosynthesis (cyanobacteria)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis
PWY-6057	dimethyl sulfide degradation III (oxidation)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Dimethylsulfide-Degradation
PWY-6047	dimethyl sulfide degradation I	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Dimethylsulfide-Degradation
PWY-7619	juniperonate biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis
PWY-6059	dimethyl sulfide degradation II (oxidation)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Dimethylsulfide-Degradation
PWY-6049	superpathway of dimethylsulfoniopropanoate degradation	Super-Pathways
PWY-7726	(4Z,7Z,10Z,13Z,16Z)-docosa-4,7,10,13,16-pentaenoate biosynthesis (6-desaturase)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis
PWY-6046	dimethylsulfoniopropanoate degradation I (cleavage)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Dimethylsulfoniopropionate-Degradation
PWY-6056	dimethylsulfoniopropanoate degradation II (cleavage)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Dimethylsulfoniopropionate-Degradation
PWY-7728	(4Z,7Z,10Z,13Z,16Z)-docosa-4,7,10,13,16-pentaenoate biosynthesis II (4-desaturase)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis
PWY-6052	dimethylsulfoniopropanoate degradation III (demethylation)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Dimethylsulfoniopropionate-Degradation
PWY0-1534	hydrogen sulfide biosynthesis I	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Hydrogen-Sulfide-Biosynthesis
PWY-7601	arachidonate biosynthesis IV (8-detaturase, lower eukaryotes)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis|Arachidonate-Biosynthesis
PWY66-426	hydrogen sulfide biosynthesis II (mammalian)	Metabolic-Clusters
PWY-6045	methylthiopropanonate degradation II (demethylation)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Methylthiopropionate-Degradation
PWY-7583	arachidonate biosynthesis II (bacteria)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis|Arachidonate-Biosynthesis
PWY-6048	methylthiopropanoate degradation I (cleavage)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Methylthiopropionate-Degradation
DISSULFRED-PWY	sulfate reduction IV (dissimilatory, to hydrogen sufide))	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfate-Reduction
PWY-7725	arachidonate biosynthesis V (8-detaturase, mammals)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis|Arachidonate-Biosynthesis
SO4ASSIM-PWY	sulfate reduction I (assimilatory)	Super-Pathways
PWY-6683	sulfate reduction III (assimilatory)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfate-Reduction
PWY-5353	arachidonate biosynthesis I (6-desaturase, lower eukaryotes)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis|Arachidonate-Biosynthesis
P224-PWY	sulfate reduction V (dissimilatory, to thiosulfate)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfate-Reduction
SULFMETII-PWY	sulfate reduction II (assimilatory)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfate-Reduction
PWY-5274	sulfide oxidation II (sulfide dehydrogenase)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfide-Oxidation
PWY-5335	superpathway of sulfide oxidation (Starkeya novella)	Super-Pathways
PWY-5294	superpathway of sulfide oxidation (Acidithiobacillus ferrooxidans)	Super-Pathways
PWY-6676	superpathway of sulfide oxidation (phototrophic sulfur bacteria)	Super-Pathways
PWY-5278	sulfite oxidation III	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfite-Oxidation
PWY-5279	sulfite oxidation II	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfite-Oxidation
PWY-5326	sulfite oxidation IV	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfite-Oxidation
PWY-1281	sulfoacetaldehyde degradation I	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfoacetaldehyde-Degradation
PWY-6718	sulfoacetaldehyde degradation III	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfoacetaldehyde-Degradation
PWY-5982	sulfoacetaldehyde degradation II	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfoacetaldehyde-Degradation
PWY-6637	sulfolactate degradation II	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfolactate-Degradation
PWY-6616	sulfolactate degradation I	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfolactate-Degradation
PWY-6641	superpathway of sulfolactate degradation	Super-Pathways
PWY-6638	sulfolactate degradation III	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfolactate-Degradation
PWY-5304	superpathway of sulfur oxidation (Acidianus ambivalens)	Super-Pathways
PWY-6675	sulfur oxidation IV (intracellular sulfur)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfur-Oxidation
PWY-5332	sulfur reduction I	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfur-reduction
PWY-5364	sulfur reduction II (via polysulfide)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Sulfur-reduction
PWY-1263	taurine degradation I	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Taurine-degradation
PWY0-981	taurine degradation IV	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Taurine-degradation
PWY-1541	superpathway of taurine degradation	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Taurine-degradation
PWY-1264	taurine degradation II	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Taurine-degradation
TAURINEDEG-PWY	taurine degradation III	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Taurine-degradation
PWY-5359	tetrathionate reductiuon II (to trithionate)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Tetrathionate-Reduction
PWY-5360	superpathway of tetrathionate reduction (Salmonella typhimurium)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Tetrathionate-Reduction
PWY-5350	thiosulfate disproportionation IV (rhodanese)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Thiosulfate-Disproportionation
THIOSULFOX-PWY	thiosulfate oxidation I (to tetrathionate)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Thiosulfate-Oxidation
PWY-5296	thiosulfate oxidation III (multienzyme complex)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Thiosulfate-Oxidation
FASYN-ELONG-PWY	fatty acid elongation -- saturated	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis
PWY-6677	thiosulfate oxidation IV (multienzyme complex)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Thiosulfate-Oxidation
PWY-4381	fatty acid biosynthesis initiation I	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis
PWY-5303	thiosulfate oxidation II (via tetrathionate)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Thiosulfate-Oxidation
PWY-5080	very long chain fatty acid biosynthesis I	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis
PWY-5965	fatty acid biosynthesis initiation III	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis
PWY-5966	fatty acid biosynthesis initiation II	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis
PWY-5970	fatty acids biosynthesis (yeast)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis
PWY-6429	ricinoleate biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis
PWY-6799	fatty acid biosynthesis (plant mitochondria)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis
PWY-7036	very long chain fatty acid biosynthesis II	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis
PWY-7094	fatty acid salvage	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis
PWY-2261	ascorbate glutathione cycle	Detoxification
PWY-7388	octanoyl-[acyl-carrier protein] biosynthesis (mitochondria, yeast)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis
PWY-6502	oxidized GTP and dGTP detoxification	Metabolic-Clusters
PWYG-321	mycolate biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis
PWY-6577	farnesylcysteine salvage pathway	Detoxification
PWY-6646	fluoroacetate degradation	Detoxification
PWY-6745	phytochelatins biosynthesis	Detoxification
PWY-6786	detoxification of reactive carbonyls in chloroplasts	Metabolic-Clusters
PWY-6841	homophytochelatin biosynthesis	Detoxification
PWY-6842	glutathione-mediated detoxification II	Detoxification
PWY-6997	furfural degradation	Detoxification
PWY1G-1	mycothiol-mediated detoxification	Detoxification
PWY-5361	(5Z)-icosenoate biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis
PWY-6931	seleno-amino acid detoxification and volatilization I	Detoxification|Seleno-Amino-Acid-Detoxification
PWY-5362	sapienate biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis
PWY-6935	seleno-amino acid detoxification and volatilization II	Detoxification|Seleno-Amino-Acid-Detoxification
PWY-5367	petroselinate biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis
PWY-6933	seleno-amino acid detoxification and volatilization III	Detoxification|Seleno-Amino-Acid-Detoxification
PWY-5973	cis-vaccenate biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis
PHOSPHONOTASE-PWY	2-aminoethylphosphonate degradation I	Degradation|Noncarbon-Nutrients|Phosphorus-Compounds|2-Aminoethylphosphonate-Degradation
PWY-6014	vernolate biosynthesis I	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Epoxy-Fatty-Acids-Biosynthesis|Vernolate-Biosynthesis
PWY-7447	2-aminoethylphosphonate degradation III	Degradation|Noncarbon-Nutrients|Phosphorus-Compounds|2-Aminoethylphosphonate-Degradation
PWY-7618	ricinoleate biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis
PWY-6832	2-aminoethylphosphonate degradation II	Degradation|Noncarbon-Nutrients|Phosphorus-Compounds|2-Aminoethylphosphonate-Degradation
PWY-7634	vernolate biosynthesis II	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Epoxy-Fatty-Acids-Biosynthesis|Vernolate-Biosynthesis
PWY-7806	glyphosate degradation II	Degradation|Noncarbon-Nutrients|Phosphorus-Compounds|Glyphosate-Degradation
PWY-7663	gondoate biosynthesis (anaerobic)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis
PWY-7804	glyphosate degradation I	Degradation|Noncarbon-Nutrients|Phosphorus-Compounds|Glyphosate-Degradation
PWY-7691	10,13-epoxy-11-methyl-octadecadienoate biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis
PWY-7807	glyphosate degradation III	Degradation|Noncarbon-Nutrients|Phosphorus-Compounds|Glyphosate-Degradation
PWY0-862	(5Z)-dodec-5-enoate biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis
PWY-7399	methylphosphonate degradation II	Degradation|Noncarbon-Nutrients|Phosphorus-Compounds|Methylphosphonate-Degradation
PWY-6000	&gamma;-linolenate biosynthesis II (animals)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|Gamma-linolenate-Biosynthesis
PWY0-1533	methylphosphonate degradation I	Degradation|Noncarbon-Nutrients|Phosphorus-Compounds|Methylphosphonate-Degradation
PWY-5998	&gamma;-linolenate biosynthesis I (plants)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|Gamma-linolenate-Biosynthesis
PWY0-662	PRPP biosynthesis I	Degradation|Noncarbon-Nutrients|Phosphorus-Compounds|PRPP-Biosynthesis
PWY-7594	&gamma;-linolenate biosynthesis III (cyanobacteria)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|Gamma-linolenate-Biosynthesis
PWY0-661	PRPP biosynthesis II	Degradation|Noncarbon-Nutrients|Phosphorus-Compounds|PRPP-Biosynthesis
PWY-5996	oleate biosynthesis II (animals and fungi)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|Oleate-Biosynthesis
AMMASSIM-PWY	ammonia assimilation cycle III	Super-Pathways
PWY-5147	oleate biosynthesis I (plants)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|Oleate-Biosynthesis
PWY-6963	ammonia assimilation cycle I	Super-Pathways
PWY-7664	oleate biosynthesis IV (anaerobic)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|Oleate-Biosynthesis
PWY-3282	superpathway of ammonia assimilation (plants)	Super-Pathways
PWY-7587	oleate biosynthesis III (cyanobacteria)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|Oleate-Biosynthesis
PWY-6964	ammonia assimilation cycle II	Degradation|Noncarbon-Nutrients|NITROGEN-DEG|Ammonia-Assimilation
PWY-5189	tetrapyrrole biosynthesis II (from glycine)	Biosynthesis|Cofactor-Biosynthesis|Tetrapyrrole-Biosynthesis
PWY-7084	nitrifier denitrification	Super-Pathways
PWY-5188	tetrapyrrole biosynthesis I (from glutamate)	Biosynthesis|Cofactor-Biosynthesis|Tetrapyrrole-Biosynthesis
PWY-2242	ammonia oxidation III	Degradation|Noncarbon-Nutrients|NITROGEN-DEG|Ammonia-oxidation
PWY-6748	nitrate reduction VII (denitrification)	Degradation|Noncarbon-Nutrients|NITROGEN-DEG|Nitrate-Reduction
DENITRIFICATION-PWY	nitrate reduction I (denitrification)	Degradation|Noncarbon-Nutrients|NITROGEN-DEG|Nitrate-Reduction
PWY-5674	nitrate reduction IV (dissimilatory)	Degradation|Noncarbon-Nutrients|NITROGEN-DEG|Nitrate-Reduction
PWY-381	nitrate reduction II (assimilatory)	Degradation|Noncarbon-Nutrients|NITROGEN-DEG|Nitrate-Reduction
PWY-5750	itaconate biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis
PWY490-3	nitrate reduction VI (assimilatory)	Degradation|Noncarbon-Nutrients|NITROGEN-DEG|Nitrate-Reduction
PWY-5782	2-keto-L-gulonate biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis
PWY-5675	nitrate reduction V (assimilatory)	Degradation|Noncarbon-Nutrients|NITROGEN-DEG|Nitrate-Reduction
PWY-7576	nitrogen fixation II (flavodoxin)	Degradation|Noncarbon-Nutrients|NITROGEN-DEG|Nitrogen-Fixation
N2FIX-PWY	nitrogen fixation I (ferredoxin)	Degradation|Noncarbon-Nutrients|NITROGEN-DEG|Nitrogen-Fixation
PWY-6821	&kappa;-carrageenan degradation	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Carrageenan-Degradation
PWY-6817	&lambda;-carrageenan degradation	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Carrageenan-Degradation
PWY-6822	&iota;-carrageenan degradation	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Carrageenan-Degradation
PWY-6805	cellulose degradation I (cellulosome)	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Glycan-Pathways
PWY-6788	cellulose degradation II (fungi)	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Cellulose-Degradation
PWY-6784	cellulose and hemicellulose degradation (cellulolosome)	Super-Pathways
PWY-6855	chitin degradation I (archaea)	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Chitin-Degradation
PWY-7822	chitin degradation III (Serratia)	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Chitin-Degradation
PWY-6906	chitin derivatives degradation	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Chitin-Degradation
PWY-6902	chitin degradation II (Vibrio)	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Chitin-Degradation
PWY-6790	L-arabinan degradation	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Glycan-Pathways
PWY-7719	CMP-diacetamido-8-epilegionaminic acid biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|Lipopolysaccharide-Biosynthesis|O-Antigen-Biosynthesis
PWY-7243	pectin degradation I	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Pectin-Degradation
PWY-1269	CMP-3-deoxy-D-<I>manno</I>-octulosonate biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|CMP-Sugar-Biosynthesis
PWY-6717	(1,4)-&beta;-xylan degradation	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Xyloglucan-Degradation
PWY-6144	CMP-N-glycoloylneuraminate biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|CMP-Sugar-Biosynthesis
PWY-842	starch degradation I	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Starch-Degradation
PWY-7674	CMP-8-amino-3,8-dideoxy-D-<I>manno</I>-octulosonate biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|CMP-Sugar-Biosynthesis
PWY-6576	dermatan sulfate degradation (metazoa)	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Glycan-Pathways
PWY-6143	CMP-pseudaminate biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|CMP-Sugar-Biosynthesis
PWY-6791	xyloglucan degradation I (endoglucanase)	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Xyloglucan-Degradation
PWY-7529	CMP-N-acetyl-7-O-acetylneuraminate biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|CMP-Sugar-Biosynthesis
PWY-7586	&beta;-1,4-D-mannosyl-N-acetyl-D-glucosamine degradation	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Glycan-Pathways
PWY-6140	CMP-2-keto-3-deoxy-D-glycero-D-galacto-nononate biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|CMP-Sugar-Biosynthesis
PWY-6573	chondroitin sulfate degradation (metazoa)	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Glycan-Pathways
PWY-5941	glycogen degradation II	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Glycogen-Degradation
GLYCOCAT-PWY	glycogen degradation I	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Glycogen-Degradation
PWY-7662	glycogen degradation III (via anhydrofructose)	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Glycogen-Degradation
PWY-7646	dermatan sulfate degradation I (bacterial)	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Glycosaminoglycan-Degradation
PWY-7645	hyaluronan degradation	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Glycosaminoglycan-Degradation
PWY-6572	chondroitin sulfate degradation I (bacterial)	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Glycosaminoglycan-Degradation
PWY-7644	heparin degradation	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Glycosaminoglycan-Degradation
PWY-7651	heparan sulfate degradation	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Glycosaminoglycan-Degradation
PWY-7246	pectin degradation II	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Pectin-Degradation
PWY-7248	pectin degradation III	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Pectin-Degradation
PWY-6771	rhamnogalacturonan type I degradation II (bacteria)	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Rhamnogalacturonan-Degradation
PWY-6769	rhamnogalacturonan type I degradation I (fungi)	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Rhamnogalacturonan-Degradation
PWY-6724	starch degradation II	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Starch-Degradation
PWY-6737	starch degradation V	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Starch-Degradation
PWY-6735	starch degradation IV	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Starch-Degradation
PWY-6731	starch degradation III	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Starch-Degradation
PWY-6789	(1,3)-&beta;-D-xylan degradation	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Xylan-Degradation
PWY-5317	hyoscyamine and scopolamine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|TROPANE-ALKALOIDS
PWY-6807	xyloglucan degradation II (exoglucanase)	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Xyloglucan-Degradation
PWY-6812	xyloglucan degradation III (cellobiohydrolase)	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG|Xyloglucan-Degradation
PWY0-1299	arginine dependent acid resistance	Detoxification|Acid-Resistance
PWY-5846	colchicine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|TROPANE-ALKALOIDS
PWY-6471	peptidoglycan biosynthesis IV (Enterococcus faecium)	Super-Pathways
PWY-6470	peptidoglycan biosynthesis V (&beta;-lactam resistance)	Super-Pathways
PWY-7096	triclosan resistance	Detoxification|Antibiotic-Resistance
PWY-6828	linezolid resistance	Detoxification|Antibiotic-Resistance
PWY-5843	cocaine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|TROPANE-ALKALOIDS
PWY0-1338	polymyxin resistance	Biosynthesis|Cell-Structure-Biosynthesis|Lipopolysaccharide-Biosynthesis
PWY-4202	arsenate detoxification I (glutaredoxin)	Detoxification|Arsenate-Detoxification
PWY-6421	arsenate detoxification III (mycothiol)	Detoxification|Arsenate-Detoxification
PWY-4621	arsenate detoxification II (glutaredoxin)	Detoxification|Arsenate-Detoxification
PWY-5318	calystegine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|TROPANE-ALKALOIDS
PWY-7142	cyanide detoxification II	Detoxification|Cyanide-Detoxification
P641-PWY	phenylmercury acetate degradation	Detoxification|Mercury-Detoxification
DETOX1-PWY-1	reactive oxygen species degradation	Super-Pathways
DETOX1-PWY	superoxide radicals degradation	Detoxification|REACTIVE-OXYGEN-SPECIES-DEGRADATION
TOLUENE-DEG-DIOL-PWY	toluene degradation to 2-oxopent-4-enoate (<I>via</I> toluene-<I>cis</I>-diol)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|TOLUENE-DEG
PWY-5155	&beta;-alanine biosynthesis III	Biosynthesis|Amino-Acid-Biosynthesis|Other-Amino-Acid-Biosynthesis|Beta-Alanine-Biosynthesis
PWY-3982	uracil degradation I (reductive)	Degradation|NUCLEO-DEG|Pyrimidine-Degradation|Pyrimidine-Base-Degradation|Uracil-Degradation
PWY-3981	&beta;-alanine biosynthesis I	Biosynthesis|Amino-Acid-Biosynthesis|Other-Amino-Acid-Biosynthesis|Beta-Alanine-Biosynthesis
PWY-5760	&beta;-alanine biosynthesis IV	Biosynthesis|Amino-Acid-Biosynthesis|Other-Amino-Acid-Biosynthesis|Beta-Alanine-Biosynthesis
PWY-7700	4-methylphenol degradation to protocatechuate	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|TOLUENE-DEG
PWY-3941	&beta;-alanine biosynthesis II	Biosynthesis|Amino-Acid-Biosynthesis|Other-Amino-Acid-Biosynthesis|Beta-Alanine-Biosynthesis
CITRULBIO-PWY	L-citrulline biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|Other-Amino-Acid-Biosynthesis|Citrulline-Biosynthesis
PWY-5004	superpathway of L-citrulline metabolism	Super-Pathways
PWY-4921	protein citrullination	Biosynthesis|Amino-Acid-Biosynthesis|Other-Amino-Acid-Biosynthesis|Citrulline-Biosynthesis
TOLUENE-DEG-2-OH-PWY	toluene degradation to 2-oxopent-4-enoate I (<I>via</I> o</I>-cresol)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|TOLUENE-DEG
ARGININE-SYN4-PWY	L-ornithine biosynthesis II	Biosynthesis|Amino-Acid-Biosynthesis|Other-Amino-Acid-Biosynthesis|L-Ornithine-Biosynthesis
PWY-6922	L-N<sup>&delta;</sup>-acetylornithine biosynthesis	Biosynthesis|Amino-Acid-Biosynthesis|Other-Amino-Acid-Biosynthesis|L-Ornithine-Biosynthesis
GLUTORN-PWY	L-ornithine biosynthesis I	Biosynthesis|Amino-Acid-Biosynthesis|Other-Amino-Acid-Biosynthesis|L-Ornithine-Biosynthesis
PWY-7177	UTP and CTP dephosphorylation II	Degradation|NUCLEO-DEG|Pyrimidine-Degradation|Pyrimidine-Ribonucleosides-Degradation|UTP-CTP-Dephosphorylation
TOLUENE-DEG-CATECHOL-PWY	toluene degradation to benzoate	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|TOLUENE-DEG
PWY-7185	UTP and CTP dephosphorylation I	Degradation|NUCLEO-DEG|Pyrimidine-Degradation|Pyrimidine-Ribonucleosides-Degradation|UTP-CTP-Dephosphorylation
PWY0-1471	uracil degradation III	Degradation|NUCLEO-DEG|Pyrimidine-Degradation|Pyrimidine-Base-Degradation|Uracil-Degradation
PWY-6426	uracil degradation II (oxidative)	Degradation|NUCLEO-DEG|Pyrimidine-Degradation|Pyrimidine-Base-Degradation|Uracil-Degradation
PWY-81	toluene degradation to benzoyl-CoA (anaerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|TOLUENE-DEG
PWY-6529	chlorate reduction	Energy-Metabolism|Respiration|ANAEROBIC-RESPIRATION
PWY-4601	arsenate reduction (respiratory)	Energy-Metabolism|Respiration|ANAEROBIC-RESPIRATION
PWY-6530	perchlorate reduction	Energy-Metabolism|Respiration|ANAEROBIC-RESPIRATION
TOLUENE-DEG-4-OH-PWY	toluene degradation to 4-methylphenol	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|TOLUENE-DEG
P283-PWY	hydrogen oxidation I (aerobic)	Energy-Metabolism|CHEMOAUTOTROPHIC-ENERGY-METABOLISM|Hydrogen-Oxidation
PWY-1381	fluorene degradation II	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Fluorene-Degradation
PWY-6512	hydrogen oxidation III (anaerobic, NADP)	Energy-Metabolism|CHEMOAUTOTROPHIC-ENERGY-METABOLISM|Hydrogen-Oxidation
PWY-5382	hydrogen oxidation II (aerobic, NAD)	Energy-Metabolism|CHEMOAUTOTROPHIC-ENERGY-METABOLISM|Hydrogen-Oxidation
FLUORENE-DEG-9-ONE-PWY	fluorene degradation I	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Fluorene-Degradation
PWY-7789	erythritol degradation I	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ALCOHOLS-DEG|Erythritol-Degradation
PWY-7703	2,4-xylenol degradation to protocatechuate	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|TOLUENE-DEG
PWY-7788	erythritol degradation II	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ALCOHOLS-DEG|Erythritol-Degradation
GALLATE-DEGRADATION-I-PWY	gallate degradation II	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|GALLATE-DEG
PWY-4781	phytate degradation II	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|CYCLITOLS-DEG|Phytate-Degradation
P3-PWY	gallate degradation III (anaerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|GALLATE-DEG
PWY-4702	phytate degradation I	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|CYCLITOLS-DEG|Phytate-Degradation
GALLATE-DEGRADATION-II-PWY	gallate degradation I	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|GALLATE-DEG
PWY-6989	(-)-camphor degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|TERPENOID-DEG|Camphor-Degradation
TOLUENE-DEG-3-OH-PWY	toluene degradation to 2-oxopent-4-enoate (<I>via 4-methylcatechol)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|TOLUENE-DEG
P601-PWY	(+)-camphor degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|TERPENOID-DEG|Camphor-Degradation
PWY-5927	(4S)-carveol and (4S)-dihydrocarveol degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|TERPENOID-DEG|Carveol-Degradation
PWY-5922	(4R)-carveol and (4R)-dihydrocarveol degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|TERPENOID-DEG|Carveol-Degradation
PWY-6412	ginsenoside degradation II	Degradation|SECONDARY-METABOLITE-DEGRADATION|TERPENOID-DEG|Ginsenoside-Degradation
PWY-7100	spinosyn A biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Insecticides-Biosynthesis
PWY-6411	ginsenoside degradation I	Degradation|SECONDARY-METABOLITE-DEGRADATION|TERPENOID-DEG|Ginsenoside-Degradation
PWY-5133	cohumulone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|TERPENOPHENOLICS-SYN
PWY-5924	limonene degradation II (L-limonene)	Degradation|SECONDARY-METABOLITE-DEGRADATION|TERPENOID-DEG|Limonene-Degradation
PWY-5923	limonene degradation I (D-limonene)	Degradation|SECONDARY-METABOLITE-DEGRADATION|TERPENOID-DEG|Limonene-Degradation
PWY-6526	limonene degradation III (to perillate)	Degradation|SECONDARY-METABOLITE-DEGRADATION|TERPENOID-DEG|Limonene-Degradation
PWY-5808	hyperforin and adhyperforin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|TERPENOPHENOLICS-SYN
PWY-5461	betanidin degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|N-CONTAINING-SECONDARY-CMPD-DEG|ALKALOIDS
PWY-6999	theophylline degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|N-CONTAINING-SECONDARY-CMPD-DEG|ALKALOIDS
PWY-5132	humulone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|TERPENOPHENOLICS-SYN
P23-PWY	reductive TCA cycle I	Degradation|C1-COMPOUNDS|CO2-Fixation|Autotrophic-CO2-Fixation|Reductive-TCA-Cycles
PWY-5392	reductive TCA cycle II	Degradation|C1-COMPOUNDS|CO2-Fixation|Autotrophic-CO2-Fixation|Reductive-TCA-Cycles
PWY-5140	cannabinoid biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|TERPENOPHENOLICS-SYN
P42-PWY	incomplete reductive TCA cycle	Degradation|C1-COMPOUNDS|CO2-Fixation|Autotrophic-CO2-Fixation|Reductive-TCA-Cycles
PWY-7105	olivetol biosynthesis (olivetol synthase by-products synthesis)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|TERPENOPHENOLICS-SYN
PWY-7478	oryzalexin D and E biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHYTOALEXIN-SYN|TERPENOID-PHYTOALEXINS
PWY-7489	oryzalexin A, B, and C biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHYTOALEXIN-SYN|TERPENOID-PHYTOALEXINS
PWY-7477	momilactone A biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHYTOALEXIN-SYN|TERPENOID-PHYTOALEXINS
PWY-7101	5-deoxystrigol biosynthesis	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis
PWY-5827	heliocides biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis
PWY-5828	lacinilene C biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis
PWY-6668	(E,E)-4,8,12-trimethyltrideca-1,3,7,11-tetraene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis
PWY-7087	2-methylisoborneol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis
PWY-2921	capsidiol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHYTOALEXIN-SYN|TERPENOID-PHYTOALEXINS
PWY-2961	sesquiterpenoid phytoalexins biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHYTOALEXIN-SYN|TERPENOID-PHYTOALEXINS
PWY-5425	oleoresin sesquiterpene volatiles biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-5434	(3E)-4,8-dimethylnona-1,3,7-triene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-5733	germacrene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-5773	gossypol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-5950	geosmin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-6128	cis-calamenene related sesquiterpenoids biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-6254	santalene biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-6257	curcumene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-6258	patchoulol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-6265	zerumbone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-6271	eudesmol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-6275	&beta;-caryophyllene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-6278	botrydial biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-6290	&beta;-cubebene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-6291	valencene and 7-epi-&alpha;-selinene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-6294	selinene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-6669	&delta;-guaiene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-6836	santalene biosynthesis II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-6888	zealexin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHYTOALEXIN-SYN|TERPENOID-PHYTOALEXINS
PWY-7711	calonectrin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-7712	deoxynivalenol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-7713	nivalenol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-7714	harzianum A and trichodermin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-7730	T-2 toxin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY4FS-17	abscisic acid biosynthesis shunt	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|Abscisic-Acid-Biosynthesis
PWY-5292	vindoline and vinblastine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TERPENOID-ALKALOIDS
PWY-5290	secologanin and strictosidine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TERPENOID-ALKALOIDS
PWY-5301	ajmaline and sarpagine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TERPENOID-ALKALOIDS
PWY-7532	acetylaszonalenin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|INDOLE-ALKALOIDS
PWY-6493	chanoclavine I aldehyde biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|INDOLE-ALKALOIDS
PWY-7520	terrequinone A biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|QUINONE-SYN
PWY-5877	beta-carboline biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|INDOLE-ALKALOIDS
PWY-7059	fumigaclavine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|INDOLE-ALKALOIDS
PWY-5467	gramine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|INDOLE-ALKALOIDS
PWY-6495	ergotamine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|INDOLE-ALKALOIDS
PWY-5645	4-chloronitrobenzene degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nitroaromatic-Degradation|Nitrobenzene-Degradation
PWY-5640	nitrobenzene degradation II	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nitroaromatic-Degradation|Nitrobenzene-Degradation
PWY-5637	nitrobenzene degradation I	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nitroaromatic-Degradation|Nitrobenzene-Degradation
PWY-6845	nitric oxide biosynthesis (plants)	Biosynthesis|Metabolic-Regulators|Nitric-Oxide-Biosynthesis
PWY-4421	curcumin glucoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY-5882	epoxypseudoisoeugenol-2-methylbutanoate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN
PWY-6673	caffeoylglucarate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN
PWY-6835	6-gingerol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN
PWY-7468	benzoyl-&beta;-D-glucopyranose biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN
PWY-7643	coniferyl alcohol 9-methyl ester biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN
PWY-7682	arabidopyrone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN
PWY-5668	cardiolipin biosynthesis I	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis|Cardiolipin-Biosynthesis
PWY-5269	cardiolipin biosynthesis II	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis|Cardiolipin-Biosynthesis
PWY-7632	hinokiresinol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|NORLIGNANS
PWY0-1545	cardiolipin biosynthesis III	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis|Cardiolipin-Biosynthesis
PWY-4985	mimosine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|NON-PROTEIN-AMINO-ACID-SYN
PWY-5110	trigonelline biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|PYRROLIDINE-ALKALOIDS
PWY-4961	&beta;-pyrazole-1-ylalanine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|NON-PROTEIN-AMINO-ACID-SYN
PWY-6232	chrysoeriol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN
PWY-6325	echinatin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN
PWY-5021	willardiine and isowillardiine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|NON-PROTEIN-AMINO-ACID-SYN
PWY-6515	phloridzin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN
PWY-6787	flavonoid biosynthesis (in equisetum)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN
PWY-1	lathyrine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|NON-PROTEIN-AMINO-ACID-SYN
PWY-7397	naringenin biosynthesis (engineered)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN
PWY1F-FLAVSYN	flavonoid biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN
PWY-5	canavanine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|NON-PROTEIN-AMINO-ACID-SYN
PWY-7627	2,4,6-trinitrophenol and 2,4-dinitrophenol degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nitroaromatic-Degradation|Nitrophenol-Degradation
PWY-5487	4-nitrophenol degradation I	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nitroaromatic-Degradation|Nitrophenol-Degradation|p-Nitrophenol-Degradation
PWY-5488	4-nitrophenol degradation II	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nitroaromatic-Degradation|Nitrophenol-Degradation|p-Nitrophenol-Degradation
PWY-2981	diterpene phytoalexins precursors biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHYTOALEXIN-SYN|TERPENOID-PHYTOALEXINS
PWY-7484	phytocassanes biosynthesis, shared reactions	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHYTOALEXIN-SYN|TERPENOID-PHYTOALEXINS
PWY-5660	taxol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN
PWY-7680	carnosate bioynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN
PWY-5636	2-nitrophenol degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nitroaromatic-Degradation|Nitrophenol-Degradation|o-Nitrophenol-Degradation
PWY-6645	labdane-type diterpenes biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN
PWY-7071	steviol glucoside biosynthesis (rebaudioside A biosynthesis)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN
PWY-7482	cyclooctatin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN
PWY-7572	lolitrem B biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN
PWY-6304	casbene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN
PWY-6691	plaunotol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN
PWY-7070	steviol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN
PWY-7481	oryzalide A biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHYTOALEXIN-SYN|TERPENOID-PHYTOALEXINS
PWY-5107	phytol salvage pathway	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN
PWY-7517	brassicicene C biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN
PWY-7721	methyl phomopsenoate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN
PWY-6659	fusicoccin A biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN
PWY-7065	2&alpha;,7&beta;-dihydroxylation of taxusin	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN
PWY-7667	apicidin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis|CYCLOPEPTIDES
PWY-5034	GA12 biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN|GIBBERELLIN-SYN
PWY-6653	ent -kaurene biosynthesis II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN|GIBBERELLIN-SYN
PWY-7668	apicidin F biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Antibiotic-Biosynthesis|CYCLOPEPTIDES
PWY-5032	ent-kaurene biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN|GIBBERELLIN-SYN
PWY-5035	gibberellin biosynthesis III (early C-13 hydroxylation)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN|GIBBERELLIN-SYN|GIBBERELLINS-BIOSYNTHESIS
PWY-5070	gibberellin biosynthesis I (non C-3, non C-13 hydroxylation)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN|GIBBERELLIN-SYN|GIBBERELLINS-BIOSYNTHESIS
PWY-5047	gibberellin biosynthesis IV (Gibberella fujikuroi)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN|GIBBERELLIN-SYN|GIBBERELLINS-BIOSYNTHESIS
PWY-5036	gibberellin biosynthesis II (early C-3 hydroxylation)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN|GIBBERELLIN-SYN|GIBBERELLINS-BIOSYNTHESIS
PWY-7232	gibberellin biosynthesis V	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN|GIBBERELLIN-SYN|GIBBERELLINS-BIOSYNTHESIS
PWY-3221	dTDP-L-rhamnose biosynthesis II	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis|dTDP-L-Rhamnose-Biosynthesis
PWY-7442	drosopterin and aurodrosopterin biosynthesis	Biosynthesis|Other-biosynthesis
PWY-6730	methylhalides biosynthesis (plants)	Biosynthesis|Other-biosynthesis
PWY-6498	eumelanin biosynthesis	Biosynthesis|Other-biosynthesis
PWY-6481	L-dopachrome biosynthesis	Biosynthesis|Other-biosynthesis
1CMET2-PWY	N<sup>10</sup>-formyl-tetrahydrofolate biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Folate-Biosynthesis
PWY-2161	folate polyglutamylation	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Folate-Biosynthesis
PWY-2161B	glutamate removal from folates	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Folate-Biosynthesis
PWY-6543	4-aminobenzoate biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Folate-Biosynthesis
PWY-6613	tetrahydrofolate salvage from 5,10-methenyltetrahydrofolate	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Folate-Biosynthesis
PWY-6614	tetrahydrofolate biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Folate-Biosynthesis
PWY-3841	folate transformations II	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Folate-Biosynthesis|Folate-Transformations
PWY-2201	folate transformations I	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Folate-Biosynthesis|Folate-Transformations
PWY-7396	butanol and isobutanol biosynthesis (engineered)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Alcohol-Biosynthesis
PWY-5116	sakuranetin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVANONES-SYN
PWY-5793	maysin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONE-SYN
PWY-6024	isovitexin glycosides biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONE-SYN
PWY-7188	C-glycosylflavone biosynthesis II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONE-SYN
PWY-7213	wogonin metabolism	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONE-SYN
HEXPPSYN-PWY	hexaprenyl diphosphate biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Polyprenyl-Biosynthesis
PWY-5122	geranyl diphosphate biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Polyprenyl-Biosynthesis
PWY-7444	luteolin triglucuronide biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONE-SYN
PWY-5783	octaprenyl diphosphate biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Polyprenyl-Biosynthesis
PWY-5785	di-trans,poly-cis-undecaprenyl phosphate biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Polyprenyl-Biosynthesis
PWY-5848	cinchona alkaloids biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|QUINOLINE-ALKALOIDS
PWY-5119	acacetin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONE-SYN
PWY-6326	camptothecin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|QUINOLINE-ALKALOIDS
PWY-5806	all-trans-decaprenyl diphosphate biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Polyprenyl-Biosynthesis
PWY-6015	vitexin and derivatives biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONE-SYN
PWY-5807	heptaprenyl diphosphate biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Polyprenyl-Biosynthesis
PWY-5805	nonaprenyl diphosphate biosynthesis I	Biosynthesis|Cofactor-Biosynthesis|Polyprenyl-Biosynthesis
PWY-5816	all trans undecaprenyl diphosphate biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Polyprenyl-Biosynthesis
PWY-7161	polymethylated quercetin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONE-SYN
PWY-5817	dodecaprenyl diphosphate biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Polyprenyl-Biosynthesis
PWY-5893	tridecaprenyl diphosphate biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Polyprenyl-Biosynthesis
PWY-7212	baicalein metabolism	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONE-SYN
PWY-6383	mono-trans, poly-cis decaprenyl phosphate biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Polyprenyl-Biosynthesis
PWY-6520	nonaprenyl diphosphate biosynthesis II	Biosynthesis|Cofactor-Biosynthesis|Polyprenyl-Biosynthesis
PWY-7325	salvigenin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONE-SYN
PWY-5120	geranylgeranyl diphosphate biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Polyprenyl-Biosynthesis|GGPP-Biosynthesis
PWY-5060	luteolin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONE-SYN
PWY-6010	apigenin glycosides biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONE-SYN
PWY-6239	luteolin glycosides biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONE-SYN
PWY-7189	C-glycosylflavone biosynthesis III	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONE-SYN
PWY-6075	ergosterol biosynthesis I	Biosynthesis|Lipid-Biosynthesis|Sterol-Biosynthesis|Ergosterol-Biosynthesis
PWY-7298	nevadensin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONE-SYN
PWY-7154	ergosterol biosynthesis II	Biosynthesis|Lipid-Biosynthesis|Sterol-Biosynthesis|Ergosterol-Biosynthesis
PWY-6602	C-glycosylflavone biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVANONES-SYN
PWY-5094	naringenin glycoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVANONES-SYN
PWY-5118	ponciretin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVANONES-SYN
PWY-5105	hesperitin glycoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVANONES-SYN
PWY-7377	cob(II)yrinate a,c-diamide biosynthesis I (early cobalt insertion)	Biosynthesis|Cofactor-Biosynthesis|Cobyrinate-diamide-Biosynthesis
PWY-7376	cob(II)yrinate a,c-diamide biosynthesis II (late cobalt incorporation)	Biosynthesis|Cofactor-Biosynthesis|Cobyrinate-diamide-Biosynthesis
PWY-7620	naphthalene degradation (anaerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Naphthalene-Degradation
PWY-5427	naphthalene degradation (aerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Naphthalene-Degradation
PWY-283	benzoate degradation II (aerobic and anaerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Benzoate-Degradation
PWY-2503	benzoate degradation I (aerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Benzoate-Degradation
PWY-7629	yatein biosynthesis II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|LIGNAN-SYN|Yatein-Biosynthesis
PWY-6133	(S)-reticuline biosynthesis II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|ISOQUINOLINE-ALKALOIDS|S-Reticuline-Biosynthesis
PWY-3581	(S)-reticuline biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|ISOQUINOLINE-ALKALOIDS|S-Reticuline-Biosynthesis
PWY-7748	yatein biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|LIGNAN-SYN|Yatein-Biosynthesis
PWY-5203	soybean saponin I biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-6115	avenacin biosynthesis, initial reactions	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-5429	p-xylene degradation to p-toluate	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Xylene-Degradation
PWY-5672	ginsenosides biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-7067	betulinate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-5774	saponin biosynthesis IV	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-7072	hopanoid biosynthesis (bacteria)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-5428	m-xylene degradation to m-toluate	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Xylene-Degradation
PWY-6005	marneral biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-7474	avenacin A-2 biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-6095	dammara-20,24-diene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-6109	mangrove triterpenoid biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-5002	tetrahydroxyxanthone biosynthesis (from 3-hydroxybenzoate)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|XANTHONE-SYN
PWY-7066	glycyrrhetinate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-5759	saponin biosynthesis III	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-7069	oleanolate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-5992	thalianol and derivatives biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-5001	tetrahydroxyxanthone biosynthesis (from benzoate)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|XANTHONE-SYN
PWY-7473	avenacin A-1 biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-6008	baruol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-112	lupeol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-6105	botryococcenes and methylated squalene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-5377	&alpha;-amyrin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-6767	4,4-diapolycopenedioate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-5756	saponin biosynthesis II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-5370	carbon tetrachloride degradation I	Degradation|CHLORINATED-COMPOUNDS-DEG|Carbon-tetrachloride-degradation
PWY-7068	ursolate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-5372	carbon tetrachloride degradation II	Degradation|CHLORINATED-COMPOUNDS-DEG|Carbon-tetrachloride-degradation
PWY-5875	staphyloxanthin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-6007	arabidiol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-7475	des-methyl avenacin A-1 biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-6098	diploterol and cycloartenol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|TRITERPENOID-SYN
PWY-7098	vanillin and vanillate degradation II	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Vanillin-Degradation
PWY-7097	vanillin and vanillate degradation I	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Vanillin-Degradation
PROTOCATECHUATE-ORTHO-CLEAVAGE-PWY	protocatechuate degradation II (ortho-cleavage pathway)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Protocatechuate-Degradation
P184-PWY	protocatechuate degradation I (<I>meta</I>-cleavage pathway)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Protocatechuate-Degradation
PWY-6336	protocatechuate degradation III (<I>para</I>-cleavage pathway)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Protocatechuate-Degradation
PWY-5642	2,4-dinitrotoluene degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nitroaromatic-Degradation|Dinitrotoluene-Degradation
PWY-7681	1-chloro-2-nitrobenzene degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nitroaromatic-Degradation|Dinitrotoluene-Degradation
PWY-5643	2,6-dinitrotoluene degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nitroaromatic-Degradation|Dinitrotoluene-Degradation
PWY-7585	docosahexaenoate biosynthesis II (bacteria)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis|Docosahexaenoate-Biosynthesis
PWY-7053	docosahexaenoate biosynthesis I (lower eukaryotes)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis|Docosahexaenoate-Biosynthesis
PWY-7727	docosahexaenoate biosynthesis IV (4-desaturase, mammals)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis|Docosahexaenoate-Biosynthesis
PWY-7606	docosahexaenoate biosynthesis III (6-desaturase, mammals)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis|Docosahexaenoate-Biosynthesis
PWY-7781	&omega;-sulfo-II-dihydromenaquinone-9 biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis
PWY-5872	ubiquinol-10 biosynthesis (eukaryotic)	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Ubiquinone-Biosynthesis
PWY-7230	ubiquinol-6 biosynthesis from 4-aminobenzoate (eukaryotic)	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Ubiquinone-Biosynthesis
PWY-5407	9-lipoxygenase and 9-allene oxide synthase pathway	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|FATTY-ACID-DERIVATIVE-SYN
PWY3O-19	ubiquinol-6 biosynthesis from 4-hydroxybenzoate (eukaryotic)	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Ubiquinone-Biosynthesis
PWY-5856	ubiquinol-9 biosynthesis (prokaryotic)	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Ubiquinone-Biosynthesis
PWY-5408	9-lipoxygenase and 9-hydroperoxide lyase pathway	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|FATTY-ACID-DERIVATIVE-SYN
PWY-5871	ubiquinol-9 biosynthesis (eukaryotic)	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Ubiquinone-Biosynthesis
PWY-6708	ubiquinol-8 biosynthesis (prokaryotic)	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Ubiquinone-Biosynthesis
PWY-5410	traumatin and (Z)-3-hexen-1-yl acetate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|FATTY-ACID-DERIVATIVE-SYN
PWY-5855	ubiquinol-7 biosynthesis (prokaryotic)	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Ubiquinone-Biosynthesis
PWY-5870	ubiquinol-8 biosynthesis (eukaryotic)	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Ubiquinone-Biosynthesis
PWY-735	jasmonic acid biosynthesis	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|Jasmonates-Biosynthesis
PWY-5873	ubiquinol-7 biosynthesis (eukaryotic)	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Ubiquinone-Biosynthesis
PWY-7233	ubiquinol-6 bypass biosynthesis (eukaryotic)	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Ubiquinone-Biosynthesis
PWY-5409	divinyl ether biosynthesis II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|FATTY-ACID-DERIVATIVE-SYN|Divinyl-Ether-Biosynthesis
PWY-5857	ubiquinol-10 biosynthesis (prokaryotic)	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Ubiquinone-Biosynthesis
PWY0-321	phenylacetate degradation I (aerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Phenolic-Compounds-Degradation|Phenylacetate-Degradation
PWY-1341	phenylacetate degradation II (anaerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Phenolic-Compounds-Degradation|Phenylacetate-Degradation
PWY-6998	CDP-D-arabitol biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|CDP-Sugar-Biosynthesis
PWY-7127	CDP-D-mannitol biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|CDP-Sugar-Biosynthesis
PWY-6378	putrebactin biosynthesis	Biosynthesis|Siderophores-Biosynthesis
PWY-6407	yersiniabactin biosynthesis	Biosynthesis|Siderophores-Biosynthesis
PWY-6574	achromobactin biosynthesis	Biosynthesis|Siderophores-Biosynthesis
PWY-761	rhizobactin 1021 biosynthesis	Biosynthesis|Siderophores-Biosynthesis
PWY-6289	petrobactin biosynthesis	Biosynthesis|Siderophores-Biosynthesis
PWY-6376	desferrioxamine B biosynthesis	Biosynthesis|Siderophores-Biosynthesis
PWY-6381	bisucaberin biosynthesis	Biosynthesis|Siderophores-Biosynthesis
PWY-6409	pyoverdine I biosynthesis	Biosynthesis|Siderophores-Biosynthesis
PWY-7577	ferrichrome biosynthesis	Biosynthesis|Siderophores-Biosynthesis
PWY-5925	hydroxylated mugineic acid phytosiderophore biosynthesis	Biosynthesis|Siderophores-Biosynthesis
PWY-6375	desferrioxamine E biosynthesis	Biosynthesis|Siderophores-Biosynthesis
PWY-6379	alcaligin biosynthesis	Biosynthesis|Siderophores-Biosynthesis
PWY-6408	pyochelin biosynthesis	Biosynthesis|Siderophores-Biosynthesis
AEROBACTINSYN-PWY	aerobactin biosynthesis	Biosynthesis|Siderophores-Biosynthesis
PWY-7571	ferrichrome A biosynthesis	Biosynthesis|Siderophores-Biosynthesis
PWY-5912	2-deoxymugineic acid phytosiderophore biosynthesis	Biosynthesis|Siderophores-Biosynthesis
HCAMHPDEG-PWY	3-phenylpropanoate and 3-(3-hydroxyphenyl)propanoate degradation to 2-oxopent-4-enoate	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Phenolic-Compounds-Degradation
PWY-6080	4-ethylphenol degradation (anaerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Phenolic-Compounds-Degradation
PWY-6690	cinnamate and 3-hydroxycinnamate degradation to 2-oxopent-4-enoate	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Phenolic-Compounds-Degradation
PWY-7046	4-coumarate degradation (anaerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Phenolic-Compounds-Degradation
PWY-7757	bisphenol A degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Phenolic-Compounds-Degradation
PWY-5388	N-glucosylnicotinate metabolism	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|N-CONTAINING-GLUCOSIDE-SYN
PWY-7824	prunasin and amygdalin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|N-CONTAINING-GLUCOSIDE-SYN|CYANOGENIC-GLUCOSIDE-SYN
PWY-3022	linamarin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|N-CONTAINING-GLUCOSIDE-SYN|CYANOGENIC-GLUCOSIDE-SYN
PWY-5422	isopimaric acid biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN|Resin-Acids-Biosynthesis
PWY-7095	3,4-dihydroxymandelonitrile &beta;-D-glucose biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|N-CONTAINING-GLUCOSIDE-SYN|CYANOGENIC-GLUCOSIDE-SYN
PWY-5413	neoabietic acid biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN|Resin-Acids-Biosynthesis
PWY-7088	taxiphyllin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|N-CONTAINING-GLUCOSIDE-SYN|CYANOGENIC-GLUCOSIDE-SYN
PWY-5421	dehydroabietic acid biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN|Resin-Acids-Biosynthesis
PWY-861	dhurrin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|N-CONTAINING-GLUCOSIDE-SYN|CYANOGENIC-GLUCOSIDE-SYN
PWY-5412	levopimaric acid biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN|Resin-Acids-Biosynthesis
PWY-5990	lotaustralin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|N-CONTAINING-GLUCOSIDE-SYN|CYANOGENIC-GLUCOSIDE-SYN
PWY-5411	abietic acid biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN|Resin-Acids-Biosynthesis
PWY-5414	palustric acid biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN|Resin-Acids-Biosynthesis
PWY-361	phenylpropanoid biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|PLANT-CELL-STRUCTURE|SECONDARY-CELL-WALL
PWY-1001	cellulose biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|PLANT-CELL-STRUCTURE|SECONDARY-CELL-WALL
PWY-5800	xylan biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|PLANT-CELL-STRUCTURE|SECONDARY-CELL-WALL
PWY-6302	dihydroconiferyl alcohol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|LIGNIN-SYN
PWY-7461	hydroxycinnamate sugar acid ester biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|CINNAMATE-SYN
PWY-3301	sinapate ester biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|CINNAMATE-SYN
PWY-5049	rosmarinic acid biosynthesis II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|CINNAMATE-SYN
PWY-5859	eugenol and isoeugenol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|CINNAMATE-SYN
PWY-6039	chlorogenic acid biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|CINNAMATE-SYN
PWY-7460	2-O-acetyl-3-O-trans-coutarate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|CINNAMATE-SYN
PWY-2181	free phenylpropanoid acid biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|CINNAMATE-SYN
PWY-5048	rosmarinic acid biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|CINNAMATE-SYN
PWY-5168	ferulate and sinapate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|CINNAMATE-SYN
PWY-5968	cinnamate esters biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|CINNAMATE-SYN
PWY1F-467	phenylpropanoid biosynthesis, initial reactions	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|CINNAMATE-SYN
PWY-5867	t-anethole biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|CINNAMATE-SYN
PWY-6040	chlorogenic acid biosynthesis II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|CINNAMATE-SYN
PWY-6825	phosphatidylcholine biosynthesis V	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis|PhosphatidylcholineBiosynthesis
PWY-7470	phosphatidylcholine biosynthesis VII	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis|PhosphatidylcholineBiosynthesis
PWY4FS-3	phosphatidylcholine biosynthesis III	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis|PhosphatidylcholineBiosynthesis
PWY-7416	phospholipid remodeling (phosphatidylcholine, yeast)	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis|PhosphatidylcholineBiosynthesis
PWY4FS-2	phosphatidylcholine biosynthesis II	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis|PhosphatidylcholineBiosynthesis
PWY-6826	phosphatidylcholine biosynthesis VI	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis|PhosphatidylcholineBiosynthesis
PWY3O-450	phosphatidylcholine biosynthesis I	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis|PhosphatidylcholineBiosynthesis
PWY4FS-4	phosphatidylcholine biosynthesis IV	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis|PhosphatidylcholineBiosynthesis
PWY-5039	theobromine biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|PURINE-ALKALOIDS
PWY-5040	theobromine biosynthesis II (via xanthine)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|PURINE-ALKALOIDS
PWY-5038	caffeine biosynthesis II (via paraxanthine)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|PURINE-ALKALOIDS|Caffeine-Biosynthesis
PWY-5037	caffeine biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|PURINE-ALKALOIDS|Caffeine-Biosynthesis
PWY-3061	menthol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|MONOTERPENOID-SYN
PWY-5423	oleoresin monoterpene volatiles biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|MONOTERPENOID-SYN
PWY-5813	bornyl diphosphate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|MONOTERPENOID-SYN
PWY-5829	geraniol and geranial biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|MONOTERPENOID-SYN
PWY-5835	geranyl acetate biosynthesis	Activation-Inactivation-Interconversion|Interconversion
PWY-5928	(4R)-carvone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|MONOTERPENOID-SYN
PWY-6436	perillyl aldehyde biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|MONOTERPENOID-SYN
PWY-6437	fenchol biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|MONOTERPENOID-SYN
PWY-6445	fenchol biosynthesis II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|MONOTERPENOID-SYN
PWY-6447	trichome monoterpenes biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|MONOTERPENOID-SYN
PWY-6449	fenchone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|MONOTERPENOID-SYN
PWY-6451	3-carene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|MONOTERPENOID-SYN
PWY-7182	linalool biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|MONOTERPENOID-SYN
PWY-7443	(4S)-carvone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|MONOTERPENOID-SYN
PWY-7697	geranyl &beta;-primeveroside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|MONOTERPENOID-SYN
PWY-6914	sophoraflavanone G biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|PRENYLFLAVONOID-SYN
PWY-5135	xanthohumol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|PRENYLFLAVONOID-SYN
PWY-6852	senecionine N-oxide biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|PYRROLIZIDINE-ALKALOIDS
P2-PWY	citrate lyase activation	Biosynthesis|Cofactor-Biosynthesis
P241-PWY	coenzyme B biosynthesis	Biosynthesis|Cofactor-Biosynthesis
PWY-5196	factor 430 biosynthesis	Biosynthesis|Cofactor-Biosynthesis
PWY-5198	factor 420 biosynthesis	Biosynthesis|Cofactor-Biosynthesis
PWY-5199	factor 420 polyglutamylation	Biosynthesis|Cofactor-Biosynthesis
PWY-5207	coenzyme B/coenzyme M regeneration	Biosynthesis|Cofactor-Biosynthesis
PWY-5254	methanofuran biosynthesis	Biosynthesis|Cofactor-Biosynthesis
PWY-5499	vitamin B6 degradation	Biosynthesis|Cofactor-Biosynthesis
PWY-5796	malonate decarboxylase activation	Biosynthesis|Cofactor-Biosynthesis
PWY-5915	phycoerythrobilin biosynthesis I	Biosynthesis|Cofactor-Biosynthesis
PWY-5917	phycocyanobilin biosynthesis	Biosynthesis|Cofactor-Biosynthesis
PWY-4203	volatile benzenoid biosynthesis I (ester formation)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|BENZENOID-SYN|Volatile-Benzenoids-Biosynthesis
PWY-5963	thio-molybdenum cofactor biosynthesis	Biosynthesis|Cofactor-Biosynthesis
PWY-5964	guanylyl molybdenum cofactor biosynthesis	Biosynthesis|Cofactor-Biosynthesis
PWY-6012	acyl carrier protein metabolism	Activation-Inactivation-Interconversion|Interconversion
PWY-6012-1	acyl carrier protein activation	Activation-Inactivation-Interconversion|Activation
PWY-6148	tetrahydromethanopterin biosynthesis	Biosynthesis|Cofactor-Biosynthesis
PWY-6420	pyrroloquinoline quinone biosynthesis	Biosynthesis|Cofactor-Biosynthesis
PWY-6476	cytidylyl molybdenum cofactor biosynthesis	Biosynthesis|Cofactor-Biosynthesis
PWY-7064	3-amino-3-phenylpropanoyl-CoA biosynthesis	Biosynthesis|Cofactor-Biosynthesis
PWY-7170	phytochromobilin biosynthesis	Biosynthesis|Cofactor-Biosynthesis
PWY-7250	[2Fe-2S] iron-sulfur cluster biosynthesis	Biosynthesis|Cofactor-Biosynthesis
PWY-7578	phycoviolobilin biosynthesis	Biosynthesis|Cofactor-Biosynthesis
PWY-7579	phycourobilin biosynthesis	Biosynthesis|Cofactor-Biosynthesis
PWY-7580	phycoerythrobilin biosynthesis II	Biosynthesis|Cofactor-Biosynthesis
PWY-7639	bis(guanylyl molybdenum cofactor) biosynthesis	Biosynthesis|Cofactor-Biosynthesis
PWY-7710	FeMo cofactor biosynthesis	Biosynthesis|Cofactor-Biosynthesis
PWY0-1433	tetrahydromonapterin biosynthesis	Biosynthesis|Cofactor-Biosynthesis
SAM-PWY	S-adenosyl-L-methionine biosynthesis	Biosynthesis|Cofactor-Biosynthesis
PWY-6077	anthranilate degradation II (aerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|2-Aminobenzoate-Degradation
2AMINOBENZDEG-PWY	anthranilate degradation III (anaerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|2-Aminobenzoate-Degradation
PWY-6504	anthranilate degradation IV (aerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|2-Aminobenzoate-Degradation
PWY-6079	anthranilate degradation I (aerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|2-Aminobenzoate-Degradation
PWY-6875	retinoate biosynthesis II	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Vitamin-A-Biosynthesis
PWY-7044	5-nitroanthranilate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nitroaromatic-Degradation
PWY-7403	tetramethylpyrazine degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nitroaromatic-Degradation
PWY-7628	2,4-dinitroanisole degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nitroaromatic-Degradation
PWY-6872	retinoate biosynthesis I	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Vitamin-A-Biosynthesis
PWY-2381	4-nitrobenzoate degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nitroaromatic-Degradation|Nitrobenzoate-Degradation
PWY-5648	2-nitrobenzoate degradation II	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nitroaromatic-Degradation|Nitrobenzoate-Degradation|2-Nirobenzoate-Degradation
PWY-6857	retinol biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Vitamin-A-Biosynthesis
PWY-6861	the visual cycle I (vertebrates)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Vitamin-A-Biosynthesis|Visual-Cycle
PWY-7043	11-cis-3-hydroxyretinal biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Vitamin-A-Biosynthesis|Visual-Cycle
PWY-6146	Methanobacterium thermoautotrophicum biosynthetic metabolism	Super-Pathways
PWY-7042	the visual cycle (insects)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Vitamin-A-Biosynthesis|Visual-Cycle
PWY-7041	the visual cycle II (molluscs)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Vitamin-A-Biosynthesis|Visual-Cycle
PPGPPMET-PWY	ppGpp biosynthesis	Biosynthesis|Metabolic-Regulators
PWY-6100	L-carnitine biosynthesis	Biosynthesis|Metabolic-Regulators
PWY-6158	creatine-phosphate biosynthesis	Biosynthesis|Metabolic-Regulators
PWY-6405	Rapoport-Luebering glycolytic shunt	Biosynthesis|Metabolic-Regulators
PWY-7535	lovastatin biosynthesis	Biosynthesis|Metabolic-Regulators
PWY66-420	carnosine biosynthesis	Biosynthesis|Metabolic-Regulators
PWY66-423	fructose 2,6-bisphosphate biosynthesis	Biosynthesis|Metabolic-Regulators
PWY-722	nicotinate degradation I	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nicotinate-Degradation
PWY-5055	nicotinate degradation III	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nicotinate-Degradation
PWY-5033	nicotinate degradation II	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Nicotinate-Degradation
PWY-4161	superpathway of benzoxazinoid glucosides biosynthesis	Super-Pathways
PWY-6943	testosterone and androsterone degradation to androstendione	Degradation|Steroids-Degradation
PWY-6944	androstenedione degradation	Degradation|Steroids-Degradation
PWY-5765	1,3,5-trimethoxybenzene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-6948	sitosterol degradation to androstenedione	Degradation|Steroids-Degradation
PWY-5826	hypoglycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-5935	tuberculosinol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-5954	(1S,5S)-averufin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-5955	versicolorin B biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-5956	sterigmatocystin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-5959	aflatoxins B1 and G1 biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-5960	aflatoxins B2 and G2 biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-5961	superpathway of aflatoxin biosynthesis	Super-Pathways
PWY-7378	aminopropanol phosphate biosynthesis II	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Cobalamin-Biosynthesis|Aminopropanol-Phosphate-Biosynthesis
PWY-5975	furaneol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-5443	aminopropanol phosphate biosynthesis I	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Cobalamin-Biosynthesis|Aminopropanol-Phosphate-Biosynthesis
PWY-5979	3-amino-5-hydroxybenzoate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-6068	indican biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-6268	adenosylcobalamin salvage from cobalamin	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Cobalamin-Biosynthesis|Adenosylcobalamin-Biosynthesis
PWY-6269		adenosylcobalamin salvage from cobinamide II	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Cobalamin-Biosynthesis|Adenosylcobalamin-Biosynthesis|B12-Salvage-From-Cobinamide
PWY-6069	indigo biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
COBALSYN-PWY	adenosylcobalamin salvage from cobinamide I	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Cobalamin-Biosynthesis|Adenosylcobalamin-Biosynthesis|B12-Salvage-From-Cobinamide
PWY-6330	acetaldehyde biosynthesis II	Biosynthesis|Carbohydrates-Biosynthesis|Acetaldehyde-Biosynthesis
PWY-6585	2-methylketone biosynthesis	Metabolic-Clusters
PWY-6644	fluoroacetate and fluorothreonine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-6648	rhamnolipid biosynthesis	Biosynthesis|Lipid-Biosynthesis
PWY-6660	2-heptyl-3-hydroxy-4(1H)-quinolone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-6661	4-hydroxy-2(1H)-quinolone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-6662	superpathway of quinolone and alkylquinolone biosynthesis	Super-Pathways
PWY-1061	homogalacturonan biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|PLANT-CELL-STRUCTURE|PRIMARY-CELL-WALL
PWY-6699	oxalate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-5936	xyloglucan biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|Polysaccharides-Biosynthesis
PWY-5980	xylogalacturonan biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|PLANT-CELL-STRUCTURE|PRIMARY-CELL-WALL
PWY-6703	preQ0 biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-6403	carrageenan biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|Polysaccharides-Biosynthesis
PWY-6655	xanthan biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|Polysaccharides-Biosynthesis
PWY-6795	diacylglyceryl-N,N,N-trimethylhomoserine biosynthesis	Biosynthesis|Lipid-Biosynthesis
PWY-6658	acetan biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|Polysaccharides-Biosynthesis
PWY-6773	1,3-&beta;-D-glucan biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|Polysaccharides-Biosynthesis
PWY-6802	salidroside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-822	fructan biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|Polysaccharides-Biosynthesis
PWY-6808	dTDP-D-forosamine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-6831	pyrrolnitrin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-6839	2-aminoethylphosphonate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-6920	6-gingerol analog biosynthesis (engineered)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY1G-126	mycothiol oxidation	Biosynthesis|Cofactor-Biosynthesis|Reductants
PWY-6926	pyrethrin I biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-6664	di-myo-inositol phosphate biosynthesis	Biosynthesis|Metabolic-Regulators|Organic-solutes-Biosynthesis
PWY-4121	glutathionylspermidine biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Reductants
PWY-6002	lotaustralin degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|N-CONTAINING-SECONDARY-CMPD-DEG|N-CONTAINING-GLUCOSIDE-DEG|CYANOGENIC-GLUCOSIDE-DEG
PWY-7054	N-acetylglutaminylglutamine amide biosynthesis	Biosynthesis|Metabolic-Regulators|Organic-solutes-Biosynthesis
PWY-7091	linustatin bioactivation	Degradation|SECONDARY-METABOLITE-DEGRADATION|N-CONTAINING-SECONDARY-CMPD-DEG|N-CONTAINING-GLUCOSIDE-DEG|CYANOGENIC-GLUCOSIDE-DEG
PWY-7625	phosphatidylinositol biosynthesis II (eukaryotes)	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis|Phosphatidylinositol-Biosynthesis
PWY-5976	dhurrin degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|N-CONTAINING-SECONDARY-CMPD-DEG|N-CONTAINING-GLUCOSIDE-DEG|CYANOGENIC-GLUCOSIDE-DEG
PWY-6580	phosphatidylinositol biosynthesis I (bacteria)	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis|Phosphatidylinositol-Biosynthesis
PWY-7089	taxiphyllin bioactivation	Degradation|SECONDARY-METABOLITE-DEGRADATION|N-CONTAINING-SECONDARY-CMPD-DEG|N-CONTAINING-GLUCOSIDE-DEG|CYANOGENIC-GLUCOSIDE-DEG
PWY-4661	1D-myo-inositol hexakisphosphate biosynthesis III (Spirodela polyrrhiza)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|SUGAR-DERIVS|Cyclitols-Biosynthesis|Phytate-Biosynthesis
PWY-7093	vicianin bioactivation	Degradation|SECONDARY-METABOLITE-DEGRADATION|N-CONTAINING-SECONDARY-CMPD-DEG|N-CONTAINING-GLUCOSIDE-DEG|CYANOGENIC-GLUCOSIDE-DEG
PWY-6372	1D-myo-inositol hexakisphosphate biosynthesis IV (Dictyostelium)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|SUGAR-DERIVS|Cyclitols-Biosynthesis|Phytate-Biosynthesis
PWY-3121	linamarin degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|N-CONTAINING-SECONDARY-CMPD-DEG|N-CONTAINING-GLUCOSIDE-DEG|CYANOGENIC-GLUCOSIDE-DEG
PWY-2301	myo-inositol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|SUGAR-DERIVS|Cyclitols-Biosynthesis|Myo-Inositol-Biosynthesis
PWY-6011	amygdalin and prunasin degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|N-CONTAINING-SECONDARY-CMPD-DEG|N-CONTAINING-GLUCOSIDE-DEG|CYANOGENIC-GLUCOSIDE-DEG
PWY-6949	DIBOA-glucoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-7092	neolinustatin bioactivation	Degradation|SECONDARY-METABOLITE-DEGRADATION|N-CONTAINING-SECONDARY-CMPD-DEG|N-CONTAINING-GLUCOSIDE-DEG|CYANOGENIC-GLUCOSIDE-DEG
PWY-5371	chrysophanol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|QUINONE-SYN
PWY-6684	glucosinolate breakdown (via thiocyanate-forming protein)	Degradation|SECONDARY-METABOLITE-DEGRADATION|N-CONTAINING-SECONDARY-CMPD-DEG|N-CONTAINING-GLUCOSIDE-DEG|GLUCOSINOLATE-DEG
PWY-5795	juglone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|QUINONE-SYN
PWY-5267	glucosinolate breakdown	Degradation|SECONDARY-METABOLITE-DEGRADATION|N-CONTAINING-SECONDARY-CMPD-DEG|N-CONTAINING-GLUCOSIDE-DEG|GLUCOSINOLATE-DEG
PWY-5987	sorgoleone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|QUINONE-SYN
PWYQT-4477	indole glucosinolate breakdown (active in intact plant cell)	Degradation|SECONDARY-METABOLITE-DEGRADATION|N-CONTAINING-SECONDARY-CMPD-DEG|N-CONTAINING-GLUCOSIDE-DEG|GLUCOSINOLATE-DEG
PWY-7518	atromentin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|QUINONE-SYN
PWYQT-4476	indole glucosinolate breakdown (insect chewing induced)	Degradation|SECONDARY-METABOLITE-DEGRADATION|N-CONTAINING-SECONDARY-CMPD-DEG|N-CONTAINING-GLUCOSIDE-DEG|GLUCOSINOLATE-DEG
PWY-5780	hypericin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|QUINONE-SYN
PWY-6553	caffeine degradation II	Degradation|SECONDARY-METABOLITE-DEGRADATION|N-CONTAINING-SECONDARY-CMPD-DEG|ALKALOIDS|CAFFEINE
PWY-5802	alizarin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|QUINONE-SYN
PWY-6552	caffeine degradation I (main, plants)	Degradation|SECONDARY-METABOLITE-DEGRADATION|N-CONTAINING-SECONDARY-CMPD-DEG|ALKALOIDS|CAFFEINE
PWY-7080	asterrate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|QUINONE-SYN
PWY-6633	caffeine degradation V (bacteria, via trimethylurate)	Degradation|SECONDARY-METABOLITE-DEGRADATION|N-CONTAINING-SECONDARY-CMPD-DEG|ALKALOIDS|CAFFEINE
PWY-5701	shikonin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|QUINONE-SYN
PWY-6538	caffeine degradation III (bacteria, via demethylation)	Degradation|SECONDARY-METABOLITE-DEGRADATION|N-CONTAINING-SECONDARY-CMPD-DEG|ALKALOIDS|CAFFEINE
PWY-5801	lawsone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|QUINONE-SYN
PWY-6632	caffeine degradation IV (bacteria, via demethylation and oxidation)	Super-Pathways
PWY-7079	geodin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|QUINONE-SYN
PWY-5397	crocetin biosynthesis	Degradation|SECONDARY-METABOLITE-DEGRADATION|TERPENOID-DEG|Tetraterpenoid-Degradation|Carotenoid-Degradation
PWYQT-4471	glucosinolate biosynthesis from dihomomethionine	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|N-CONTAINING-GLUCOSIDE-SYN|GLUCOSINOLATE-SYN
PWY-695	abscisic acid biosynthesis	Degradation|SECONDARY-METABOLITE-DEGRADATION|TERPENOID-DEG|Tetraterpenoid-Degradation|Carotenoid-Degradation
PWYQT-4474	glucosinolate biosynthesis from pentahomomethionine	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|N-CONTAINING-GLUCOSIDE-SYN|GLUCOSINOLATE-SYN
PWY-5250	methanogenesis from trimethylamine	Energy-Metabolism|Respiration|ANAEROBIC-RESPIRATION|METHANOGENESIS
PWY-1187	glucosinolate biosynthesis from homomethionine	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|N-CONTAINING-GLUCOSIDE-SYN|GLUCOSINOLATE-SYN
METHFORM-PWY	methyl-coenzyme M reduction to methane	Energy-Metabolism|Respiration|ANAEROBIC-RESPIRATION|METHANOGENESIS
PWYQT-4450	aliphatic glucosinolate biosynthesis, side chain elongation cycle	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|N-CONTAINING-GLUCOSIDE-SYN|GLUCOSINOLATE-SYN
PWY-5248	methanogenesis from dimethylamine	Energy-Metabolism|Respiration|ANAEROBIC-RESPIRATION|METHANOGENESIS
PWYQT-4473	glucosinolate biosynthesis from tetrahomomethionine	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|N-CONTAINING-GLUCOSIDE-SYN|GLUCOSINOLATE-SYN
PWY-5259	methanogenesis from methanethiol	Energy-Metabolism|Respiration|ANAEROBIC-RESPIRATION|METHANOGENESIS
PWY-601	glucosinolate biosynthesis from tryptophan	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|N-CONTAINING-GLUCOSIDE-SYN|GLUCOSINOLATE-SYN
PWY-6830	superpathway of methanogenesis	Super-Pathways
PWYQT-4472	glucosinolate biosynthesis from trihomomethionine	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|N-CONTAINING-GLUCOSIDE-SYN|GLUCOSINOLATE-SYN
METHANOGENESIS-PWY	methanogenesis from H<SUB>2</SUB> and CO<SUB>2</SUB>	Energy-Metabolism|Respiration|ANAEROBIC-RESPIRATION|METHANOGENESIS
PWYQT-4475	glucosinolate biosynthesis from hexahomomethionine	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|N-CONTAINING-GLUCOSIDE-SYN|GLUCOSINOLATE-SYN
PWY-5247	methanogenesis from methylamine	Energy-Metabolism|Respiration|ANAEROBIC-RESPIRATION|METHANOGENESIS
PWY-2821	glucosinolate biosynthesis from phenylalanine	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|N-CONTAINING-GLUCOSIDE-SYN|GLUCOSINOLATE-SYN
PWY-5261	methanogenesis from tetramethylammonium	Energy-Metabolism|Respiration|ANAEROBIC-RESPIRATION|METHANOGENESIS
PWY-6950	DIMBOA-glucoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
METH-ACETATE-PWY	methanogenesis from acetate	Energy-Metabolism|Respiration|ANAEROBIC-RESPIRATION|METHANOGENESIS
PWY-5475	pentagalloylglucose biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|TANNIN-SYN
PWY-5209	methyl-coenzyme M oxidation to CO2	Energy-Metabolism|Respiration|ANAEROBIC-RESPIRATION|METHANOGENESIS
PWY-5477	gallotannin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|TANNIN-SYN|GALLOTANNINS
PWY-6455	vancomycin resistance II	Detoxification|Antibiotic-Resistance|Vancomycin-Resistnace
PWY-6199	quercetin sulfate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-6454	vancomycin resistance I	Detoxification|Antibiotic-Resistance|Vancomycin-Resistnace
PWY-7173	quercetin triglucoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7137	quercetin gentiotetraside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-6917	vernolate biosynthesis III	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Epoxy-Fatty-Acids-Biosynthesis|Vernolate-Biosynthesis
PWY-7448	galloylated catechin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7040	violacein biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-3101	flavonol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7074	phenylethanol glycoconjugate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-5321	quercetin glycoside biosynthesis (Arabidopsis)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7075	phenylethyl acetate biosynthesis	Activation-Inactivation-Interconversion|Interconversion
PWY-7157	eupatolitin 3-O-glucoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-5665	vanillin biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|Vanillin-Biosynthesis
PWY-5390	rutin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7076	3,5-dimethoxytoluene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-7166	kaempferide triglycoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7136	&beta; myrcene degradation	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-6064	methylquercetin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7153	grixazone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-7172	flavonol acylglucoside biosynthesis III - quercetin derivatives	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7394	urate degradation to allantoin II	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Urate-Degradation
PWY-7129	quercetin glucoside biosynthesis (Allium)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7236	mycocyclosin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-7192	quercetin diglycoside biosynthesis (pollen-specific)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7463	N-methylanthraniloyl-&beta;-D-glucopyranose biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-7143	kaempferol gentiobioside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7490	patulin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-5320	kaempferol glycoside biosynthesis (Arabidopsis)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-5691	urate degradation to allantoin I	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Urate-Degradation
PWY-7151	polymethylated quercetin glucoside biosynthesis II - quercetagetin series (Chrysosplenium)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7492	paspaline biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-5363	chrysin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7493	paxilline and diprenylpaxilline biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-7163	polymethylated kaempferol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7531	mannojirimycin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-7171	flavonol acylglucoside biosynthesis II - isorhamnetin derivatives	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7533	gliotoxin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-6631	O -methylation of tricetin	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7534	gliotoxin inactivation	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-7191	kaempferol diglycoside biosynthesis (pollen-specific)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7536	2-amino-3-hydroxycyclopent-2-enone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-7140	myricetin gentiobioside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7540	aflatrem biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-7495	gossypetin metabolism	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7548	methylthiolincosamide biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-5059	pinobanksin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7556	(-)-microperfuranone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-7150	polymethylated quercetin glucoside biosynthesis I - quercetin series (Chrysosplenium)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7599	anditomin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-5348	kaempferol triglucoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7687	stipitatate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-7160	polymethylated myricetin biosynthesis (tomato)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7708	lyngbyatoxin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-5391	syringetin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-6190	2,4-dichlorotoluene degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorotoluene-Degradation|Dichlorotoluene-Degradation
PWY-7168	flavonol acylglucoside biosynthesis I - kaempferol derivatives	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|FLAVONOL-SYN
PWY-7751	shinorine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-7752	gadusol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-6192	3,4-dichlorotoluene degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorotoluene-Degradation|Dichlorotoluene-Degradation
PWY8J2-20	pulcherrimin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS
PWY-7765	3-hydroxy-4-methyl-anthranilate biosynthesis II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|3-OH-4-Methyl-Anthranilate-Biosynthesis
PWY-7717	3-hydroxy-4-methyl-anthranilate biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|3-OH-4-Methyl-Anthranilate-Biosynthesis
PWY-6191	2,5-dichlorotoluene degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorotoluene-Degradation|Dichlorotoluene-Degradation
PWY-4181	glutathione amide metabolism	Biosynthesis|Cofactor-Biosynthesis|Reductants
PWY8J2-1	bacillithiol biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Reductants
PWY-6840	homoglutathione biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Reductants
PWY-4081	glutathione-peroxide redox reactions	Biosynthesis|Cofactor-Biosynthesis|Reductants
GLUT-REDOX-PWY	glutathione-glutaredoxin redox reactions	Biosynthesis|Cofactor-Biosynthesis|Reductants
THIOREDOX-PWY	thioredoxin pathway	Biosynthesis|Cofactor-Biosynthesis|Reductants
PWY1G-0	mycothiol biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Reductants
GLUTATHIONESYN-PWY	glutathione biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Reductants
PWY-7256	cyanidin diglucoside biosynthesis (acyl-glucose dependent)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ANTHOCYANIN-SYN
PWY-5139	pelargonidin conjugates biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ANTHOCYANIN-SYN
PWY-5268	salvianin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ANTHOCYANIN-SYN
PWY-7280	ternatin C3 biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ANTHOCYANIN-SYN
PWY-5295	ternatin C5 biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ANTHOCYANIN-SYN
PWY-7452	cyanidin dimalonylglucoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ANTHOCYANIN-SYN
PWY-7253	apigeninidin 5-O-glucoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ANTHOCYANIN-SYN
PWY-5125	anthocyanin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ANTHOCYANIN-SYN
PWY-7260	delphinidin diglucoside biosynthesis (acyl-glucose dependent)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ANTHOCYANIN-SYN
PWY-5160	rose anthocyanin biosynthesis I (via cyanidin 5-O-&beta;-D-glucoside)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ANTHOCYANIN-SYN
PWY-7267	anthocyanin biosynthesis (pelargonidin 3-O-glucoside)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ANTHOCYANIN-SYN
PWY-7450	anthocyanidin modification (Arabidopsis)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ANTHOCYANIN-SYN
PWY-7464	cyanidin 3,7-diglucoside polyacylation biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ANTHOCYANIN-SYN
PWY-7252	luteolinidin 5-O-glucoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ANTHOCYANIN-SYN
PWY-7259	pelargonidin diglucoside biosynthesis (acyl-glucose dependent)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ANTHOCYANIN-SYN
PWY-5153	anthocyanin biosynthesis (delphinidin 3-O-glucoside)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ANTHOCYANIN-SYN
PWY-7262	rose anthocyanin biosynthesis II (via cyanidin 3-O-&beta;-D-glucoside)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ANTHOCYANIN-SYN
PWY-5284	shisonin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ANTHOCYANIN-SYN
PWY-7449	acylated cyanidin galactoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ANTHOCYANIN-SYN
PWY-5307	gentiodelphin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ANTHOCYANIN-SYN
PWY-7458	violdelphin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ANTHOCYANIN-SYN
PWY-3881	mannitol biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|Sugar-Alcohols-Biosynthesis
PWY-5530	sorbitol biosynthesis II	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|Sugar-Alcohols-Biosynthesis|Sorbitol-Biosynthesis
PWY-5054	sorbitol biosynthesis I	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|Sugar-Alcohols-Biosynthesis|Sorbitol-Biosynthesis
QUINATEDEG-PWY	quinate degradation I	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Quinate-Degradation
PWY-6416	quinate degradation II	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Quinate-Degradation
PWY-1422	vitamin E biosynthesis (tocopherols)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis
PWY-5027	phylloquinol biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Phylloquinone-Biosynthesis
PWY-7436	vitamin E biosynthesis (tocotrienols)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis
PWY3O-6	dehydro-D-arabinono-1,4-lactone biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis
PWY-7756	iso-bile acids biosynthesis II	Biosynthesis|Lipid-Biosynthesis|Sterol-Biosynthesis|Iso-bile-Acids-Biosynthesis
PWY-7755	iso-bile acids biosynthesis I	Biosynthesis|Lipid-Biosynthesis|Sterol-Biosynthesis|Iso-bile-Acids-Biosynthesis
PWY-2761	glyceollin biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ISOFLAVONOID-SYN
PWY-2083	isoflavonoid biosynthesis II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ISOFLAVONOID-SYN
GLUCONEO-PWY	gluconeogenesis I	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|Gluconeogenesis
PWY66-399	gluconeogenesis III	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|Gluconeogenesis
PWY-4502	wighteone and luteone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ISOFLAVONOID-SYN
PWY-6913	methylbutenol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|ISOPRENOIDS
PWY-5729	vestitol and sativan biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ISOFLAVONOID-SYN
PWY-7391	isoprene biosynthesis II (engineered)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|ISOPRENOIDS
NONMEVIPP-PWY	methylerythritol phosphate pathway I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|ISOPRENOIDS|Isopentenyl-Diphosphate-Biosynthesis
PWY-2464	maackiain biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ISOFLAVONOID-SYN
PWY-7560	methylerythritol phosphate pathway II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|ISOPRENOIDS|Isopentenyl-Diphosphate-Biosynthesis
PWY-5851	demethylmenaquinol-9 biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Demethylmenaquinone-Biosynthesis
PWY-922	mevalonate pathway I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|ISOPRENOIDS|Isopentenyl-Diphosphate-Biosynthesis
PWY-7633	calycosin 7-O-glucoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ISOFLAVONOID-SYN
PWY-2321	formononetin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ISOFLAVONOID-SYN
PWY-5061	6,7,4-trihydroxyisoflavone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ISOFLAVONOID-SYN
PWY-7165	L-ascorbate biosynthesis VI (engineered pathway)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Ascorbate-Biosynthesis
PWY-2463	medicarpin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ISOFLAVONOID-SYN
PWY4FS-11	L-ascorbate biosynthesis II (L-gulose pathway)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Ascorbate-Biosynthesis
PWY3DJ-35471	L-ascorbate biosynthesis IV	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Ascorbate-Biosynthesis
PWY-5825	dalpatein and dalnigrein biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ISOFLAVONOID-SYN
PWY-5521	L-ascorbate biosynthesis III	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Ascorbate-Biosynthesis
PWY-882	L-ascorbate biosynthesis I (L-galactose pathway)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Ascorbate-Biosynthesis
PWY-7145	genistin gentiobioside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ISOFLAVONOID-SYN
PWY-5740	GDP-L-colitose biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|Lipopolysaccharide-Biosynthesis|O-Antigen-Biosynthesis
PWY-66	GDP-L-fucose biosynthesis I (from GDP-D-mannose)	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|GDP-Sugar-Biosynthesis
PWY-2762	glyceollin biosynthesis II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ISOFLAVONOID-SYN
PWY-7610	GDP-6-deoxy-D-altro-heptose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|GDP-Sugar-Biosynthesis
PWY-5659	GDP-mannose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|GDP-Sugar-Biosynthesis
PWY-2002	isoflavonoid biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ISOFLAVONOID-SYN
PWY-5739	GDP-D-perosamine biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|Lipopolysaccharide-Biosynthesis|O-Antigen-Biosynthesis
PWY-6262	demethylmenaquinol-8 biosynthesis II	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Demethylmenaquinone-Biosynthesis|Demethylmenaquinol-8-Biosynthesis
PWY-3042	phaseollin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ISOFLAVONOID-SYN
PWY-7573	GDP-mycosamine biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|GDP-Sugar-Biosynthesis
PWY-6478	GDP-D-glycero-&alpha;-D-manno-heptose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|GDP-Sugar-Biosynthesis
PWY-5852	demethylmenaquinol-8 biosynthesis I	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Demethylmenaquinone-Biosynthesis|Demethylmenaquinol-8-Biosynthesis
PWY-4681	kievitone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ISOFLAVONOID-SYN
PWY-5738	GDP-6-deoxy-D-talose biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|Lipopolysaccharide-Biosynthesis|O-Antigen-Biosynthesis
PWY-5115	GDP-L-galactose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|GDP-Sugar-Biosynthesis
PWY-6793	demethylmenaquinol-8 biosynthesis III	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Demethylmenaquinone-Biosynthesis|Demethylmenaquinol-8-Biosynthesis
PWY-5821	dalcochinin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ISOFLAVONOID-SYN
PWY-7613	GDP-6-deoxy-D-manno-heptose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|GDP-Sugar-Biosynthesis
PWY-6	GDP-L-fucose biosynthesis II (from L-fucose)	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|GDP-Sugar-Biosynthesis
PWY-7372	demethylmenaquinol-6 biosynthesis II	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Demethylmenaquinone-Biosynthesis|Demethylmenaquinol-6-Biosynthesis
PWY-2467	pisatin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ISOFLAVONOID-SYN
GDPRHAMSYN-PWY	GDP-D-rhamnose biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|Lipopolysaccharide-Biosynthesis|O-Antigen-Biosynthesis
PWY-5661	GDP-glucose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|GDP-Sugar-Biosynthesis
PWY-5853	demethylmenaquinol-6 biosynthesis I	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Demethylmenaquinone-Biosynthesis|Demethylmenaquinol-6-Biosynthesis
PWY-6332	coumestrol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ISOFLAVONOID-SYN
PWY-7301	dTDP-4-O-demethyl-&beta;-L-noviose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis
PWY-7316	dTDP-N-acetylviosamine biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis
PWY-7413	dTDP-6-deoxy-&alpha;-D-allose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis
PWY-7439	dTDP-&beta;-L-evernitrose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis
PWY-6973	dTDP-D-olivose, dTDP-D-oliose and dTDP-D-mycarose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis
PWY-7688	dTDP-D-ravidosamine and dTDP-4-acetyl-D-ravidosamine biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis
CAMALEXIN-SYN	camalexin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHYTOALEXIN-SYN|INDOLE-PHYTOALEXIN-SYN
PWY-7104	dTDP-L-megosamine biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis
PWY-7315	dTDP-N-acetylthomosamine biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis
PWY-7318	dTDP-3-acetamido-3,6-dideoxy-&alpha;-D-glucose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis
PWY-6953	dTDP-3-acetamido-&alpha;-D-fucos biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis
PWY-7657	dTDP-&beta;-L-digitoxose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis
PWY-6976	dTDP-L-mycarose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis
PWY-7312	dTDP-D-&beta;-fucofuranose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis
PWY-7414	dTDP-&alpha;-D-mycaminose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis
PWY-6942	dTDP-D-desosamine biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis
PWY-7440	dTDP-&beta;-L-4-epi-vancosamine biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis
PWY-6974	dTDP-L-olivose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis
PWY-7814	dTDP-L-daunosamine biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis
PWY-7676	Kdo8N transfer to lipid IVA	Biosynthesis|Lipid-Biosynthesis|KDO-Lipid-IV-Transfer
PWY-7675	Kdo transfer to lipid IVA II	Biosynthesis|Lipid-Biosynthesis|KDO-Lipid-IV-Transfer
KDOSYN-PWY	Kdo transfer to lipid IVA I	Biosynthesis|Lipid-Biosynthesis|KDO-Lipid-IV-Transfer
PWY-5830	CDP-ascarylose biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|Lipopolysaccharide-Biosynthesis|O-Antigen-Biosynthesis
PWY-5833	CDP-4-dehydro-3,6-dideoxy-D-glucose biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|Lipopolysaccharide-Biosynthesis|O-Antigen-Biosynthesis
PWY-7335	UDP-N-acetyl-&alpha;-D-mannosaminouronate biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
PWY-7290	Escherichia coli serotype O86 O-antigen biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|Lipopolysaccharide-Biosynthesis|O-Antigen-Biosynthesis
PWY-7530	&beta;-D-galactosaminyl-(1&rarr;3)-N-acetyl-&alpha;-D-galactosamine biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|Lipopolysaccharide-Biosynthesis|O-Antigen-Biosynthesis
DTDPRHAMSYN-PWY	dTDP-L-rhamnose biosynthesis I	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|dTDP-Sugar-Biosynthesis|dTDP-L-Rhamnose-Biosynthesis
PWY-7331	UDP-N-acetyl-&beta;-L-quinovosamine biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
PWY-5832	CDP-paratose biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|Lipopolysaccharide-Biosynthesis|O-Antigen-Biosynthesis
PWY-7334	UDP-N-acetyl-&alpha;-D-quinovosamine biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
PWY-7090	UDP-2,3-diacetamido-2,3-dideoxy-&alpha;-D-mannuronate biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
PWY-7346	UDP-&alpha;-D-glucuronate biosynthesis (from UDP-glucose)	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
UDPNAGSYN-PWY	UDP-N-acetyl-D-glucosamine biosynthesis I	Biosynthesis|Polyamine-Biosynthesis|UDP-NAc-Glucosamine-Biosynthesis
PWY-7330	UDP-N-acetyl-&beta;-L-fucosamine biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
PWY-5831	CDP-abequose biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|Lipopolysaccharide-Biosynthesis|O-Antigen-Biosynthesis
PWY-7333	UDP-N-acetyl-&alpha;-D-fucosamine biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
PWY-5834	CDP-tyvelose biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|Lipopolysaccharide-Biosynthesis|O-Antigen-Biosynthesis
PWY-7336	UDP-N-acetyl-&alpha;-D-galactosaminuronate biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
PLPSAL-PWY	pyridoxal 5-phosphate salvage I	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Vitamin-B6-Biosynthesis
PWY-7204	pyridoxal 5-phosphate salvage II (plants)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Vitamin-B6-Biosynthesis
PWY-6466	pyridoxal 5-phosphate biosynthesis II	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Vitamin-B6-Biosynthesis
PYRIDOXSYN-PWY	pyridoxal 5-phosphate biosynthesis I	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Vitamin-B6-Biosynthesis
PWY-6710	poly-hydroxy fatty acids biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Hydroxy-Fatty-Acids-Biosynthesis
PWY-6433	hydroxylated fatty acid biosynthesis (plants)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Hydroxy-Fatty-Acids-Biosynthesis
PWY-6349	CDP-archaeol biosynthesis	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis
PWY-6352	3-phosphoinositide biosynthesis	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis
PWY-6367	D-myo-inositol-5-phosphate metabolism	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|SUGAR-DERIVS|Cyclitols-Biosynthesis
PWY-7039	phosphatidate metabolism, as a signaling molecule	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis
PWY-7417	phospholipid remodeling (phosphatidate, yeast)	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis
PWY-7501	phosphatidylserine biosynthesis I	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis
PWY-7506	phosphatidylserine biosynthesis II	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis
PWY-762	phospholipid desaturation	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis
PWY-7782	plasmalogen biosynthesis	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis
PWY-1361	benzoyl-CoA degradation I (aerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Benzoyl-CoA-Degradation
P321-PWY	benzoyl-CoA degradation III (anaerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Benzoyl-CoA-Degradation
CENTBENZCOA-PWY	benzoyl-CoA degradation II (anaerobic)	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Benzoyl-CoA-Degradation
PWY-5663	tetrahydrobiopterin biosynthesis I	Biosynthesis|Cofactor-Biosynthesis|Tetrahydrobiopterin-Biosynthesis
PWY-6983	tetrahydrobiopterin biosynthesis III	Biosynthesis|Cofactor-Biosynthesis|Tetrahydrobiopterin-Biosynthesis
PWY-5664	tetrahydrobiopterin biosynthesis II	Biosynthesis|Cofactor-Biosynthesis|Tetrahydrobiopterin-Biosynthesis
PWY-6654	phosphopantothenate biosynthesis III	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Pantothenate-Biosynthesis
PWY-3961	phosphopantothenate biosynthesis II	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Pantothenate-Biosynthesis
PANTO-PWY	phosphopantothenate biosynthesis I	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Pantothenate-Biosynthesis
PWY-7343	UDP-glucose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
PWY-5113	UDP-D-apiose biosynthesis (from UDP-D-glucuronate)	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
PWY-7028	UDP-N,N-diacetylbacillosamine biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
PWY-4	UDP-D-galacturonate biosynthesis II (from D-galacturonate)	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
PWY-4861	UDP-D-galacturonate biosynthesis I (from UDP-D-glucuronate)	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
PWY-63	UDP-L-arabinose biosynthesis I (from UDP-xylose)	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
PWY-7344	UDP-D-galactose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
PWY-82	UDP-L-arabinose biosynthesis II (from L-arabinose)	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
PWY-3261	UDP-L-rhamnose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|UDP-Sugar-Biosynthesis
PWY-5837	1,4-dihydroxy-2-naphthoate biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|DHNA-Biosynthesis
PWY-7304	3-amino-4,7-dihydroxy-coumarin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|COUMARIN-SYN
PWY-6792	scopoletin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|COUMARIN-SYN
PWY-5176	coumarin biosynthesis (via 2-coumarate)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|COUMARIN-SYN
PWY-5365	linear furanocoumarin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|COUMARIN-SYN
PWY-7055	daphnetin modification	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|COUMARIN-SYN
PWY-7058	esculetin modification	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|COUMARIN-SYN
PWY-7398	coumarins biosynthesis (engineered)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|COUMARIN-SYN
PWY-4922	6-methoxymellein biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY-5349	esculetin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|COUMARIN-SYN
PWY-6982	umbelliferone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|COUMARIN-SYN
PWY-6418	4-hydroxycoumarin and dicoumarol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|POLYKETIDE-SYN
PWY-7739	aucuparin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHYTOALEXIN-SYN
PWY-6366	D-myo-inositol (1,4,5,6)-tetrakisphosphate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|SUGAR-DERIVS|Cyclitols-Biosynthesis
PWY-6365	D-myo-inositol (3,4,5,6)-tetrakisphosphate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|SUGAR-DERIVS|Cyclitols-Biosynthesis
PWY-6369	inositol pyrophosphates biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|SUGAR-DERIVS|Cyclitols-Biosynthesis
PWY-6351	D-myo-inositol (1,4,5)-trisphosphate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|SUGAR-DERIVS|Cyclitols-Biosynthesis
PWY-6364	D-myo-inositol (1,3,4)-trisphosphate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|SUGAR-DERIVS|Cyclitols-Biosynthesis
PWY-5338	galactosylcyclitol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|SUGAR-DERIVS|Cyclitols-Biosynthesis
PWY-6363	D-myo-inositol (1,4,5)-trisphosphate degradation	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|SUGAR-DERIVS|Cyclitols-Biosynthesis
P183-PWY	catechol degradation to 2-oxopent-4-enoate I	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Catechol-Degradation
PWY-5339	chalcone 2-O-glucoside biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|CHALCONE-SYN
PWY-5419	catechol degradation to 2-oxopent-4-enoate II	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Catechol-Degradation
PWY-5161	6-deoxychalcone metabolism	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|CHALCONE-SYN
CATECHOL-ORTHO-CLEAVAGE-PWY	catechol degradation to &beta;-ketoadipate	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Catechol-Degradation
PWY-7469	gentisate degradation II	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Gentisate-Degradation
PWY-6223	gentisate degradation I	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Gentisate-Degradation
PWY-7772	&gamma;-resorcylate degradation II	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Gamma-Resorcylate-Degradation
PWY-7773	&gamma;-resorcylate degradation I	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Gamma-Resorcylate-Degradation
PWY-5476	cornusiin E biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|TANNIN-SYN|ELLAGITANNINS
PWY-6054	dimethylsulfoniopropanoate biosynthesis I (Wollastonia)	Biosynthesis|Metabolic-Regulators|Organic-solutes-Biosynthesis|Dimethylsulfoniopropionate-Biosynthesis
PWY-6053	dimethylsulfoniopropanoate biosynthesis III (algae)	Biosynthesis|Metabolic-Regulators|Organic-solutes-Biosynthesis|Dimethylsulfoniopropionate-Biosynthesis
PWY-6055	dimethylsulfoniopropanoate biosynthesis II (Spartina)	Biosynthesis|Metabolic-Regulators|Organic-solutes-Biosynthesis|Dimethylsulfoniopropionate-Biosynthesis
PWY-5662	glucosylglycerate biosynthesis I	Biosynthesis|Metabolic-Regulators|Organic-solutes-Biosynthesis|Glucosylglycerate-Biosynthesis
PWY-6685	glucosylglycerate biosynthesis II	Biosynthesis|Metabolic-Regulators|Organic-solutes-Biosynthesis|Glucosylglycerate-Biosynthesis
PWY-6686	mannosylglucosylglycerate biosynthesis I	Biosynthesis|Metabolic-Regulators|Organic-solutes-Biosynthesis|Mannosylglucosylglycerate-Biosynthesis
PWY-6687	mannosylglucosylglycerate biosynthesis II	Biosynthesis|Metabolic-Regulators|Organic-solutes-Biosynthesis|Mannosylglucosylglycerate-Biosynthesis
PWY-5658	mannosylglycerate biosynthesis II	Biosynthesis|Metabolic-Regulators|Organic-solutes-Biosynthesis|Mannosylglycerate-Biosynthesis
PWY-5656	mannosylglycerate biosynthesis I	Biosynthesis|Metabolic-Regulators|Organic-solutes-Biosynthesis|Mannosylglycerate-Biosynthesis
PWY-5985	trehalose biosynthesis VII	Biosynthesis|Metabolic-Regulators|Organic-solutes-Biosynthesis|Trehalose-biosynthesis
TRESYN-PWY	trehalose biosynthesis I	Biosynthesis|Metabolic-Regulators|Organic-solutes-Biosynthesis|Trehalose-biosynthesis
PWY-5983	trehalose biosynthesis VI	Biosynthesis|Metabolic-Regulators|Organic-solutes-Biosynthesis|Trehalose-biosynthesis
TREHALOSESYN-PWY	trehalose biosynthesis III	Biosynthesis|Metabolic-Regulators|Organic-solutes-Biosynthesis|Trehalose-biosynthesis
PWY-2661	trehalose biosynthesis V	Biosynthesis|Metabolic-Regulators|Organic-solutes-Biosynthesis|Trehalose-biosynthesis
PWY-881	trehalose biosynthesis II	Biosynthesis|Metabolic-Regulators|Organic-solutes-Biosynthesis|Trehalose-biosynthesis
PWY-2622	trehalose biosynthesis IV	Biosynthesis|Metabolic-Regulators|Organic-solutes-Biosynthesis|Trehalose-biosynthesis
PWY-6073	alginate biosynthesis I (algal)	Biosynthesis|Carbohydrates-Biosynthesis|Polysaccharides-Biosynthesis|Alginate-Biosynthesis
PWY-6082	alginate biosynthesis II (bacterial)	Biosynthesis|Carbohydrates-Biosynthesis|Polysaccharides-Biosynthesis|Alginate-Biosynthesis
PWY-5067	glycogen biosynthesis II (from UDP-D-Glucose)	Biosynthesis|Carbohydrates-Biosynthesis|Polysaccharides-Biosynthesis|GLYCOGEN-BIOSYN
GLYCOGENSYNTH-PWY	glycogen biosynthesis I (from ADP-D-Glucose)	Biosynthesis|Carbohydrates-Biosynthesis|Polysaccharides-Biosynthesis|GLYCOGEN-BIOSYN
PWY-622	starch biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|Polysaccharides-Biosynthesis|GLYCOGEN-BIOSYN
PWY-6558	heparan sulfate biosynthesis (late stages)	Biosynthesis|Carbohydrates-Biosynthesis|Polysaccharides-Biosynthesis|Glycosaminoglycans-Biosynthesis
PWY-6567	chondroitin sulfate biosynthesis (late stages)	Biosynthesis|Carbohydrates-Biosynthesis|Polysaccharides-Biosynthesis|Glycosaminoglycans-Biosynthesis
PWY-6557	glycosaminoglycan-protein linkage region biosynthesis	Macromolecule-Modification|Protein-Modification|Protein-Glycosylation
PWY-6566	chondroitin biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|Polysaccharides-Biosynthesis|Glycosaminoglycans-Biosynthesis
PWY-6568	dermatan sulfate biosynthesis (late stages)	Biosynthesis|Carbohydrates-Biosynthesis|Polysaccharides-Biosynthesis|Glycosaminoglycans-Biosynthesis
PWY-5508	adenosylcobalamin biosynthesis from cobyrinate a,c-diamide II	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Cobalamin-Biosynthesis|Adenosylcobalamin-Biosynthesis|B12-Biosynthesis-From-Cobyrinate-Diamide
PWY-5509	adenosylcobalamin biosynthesis from cobyrinate a,c-diamide I	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Cobalamin-Biosynthesis|Adenosylcobalamin-Biosynthesis|B12-Biosynthesis-From-Cobyrinate-Diamide
PWY-7729	5,6-dimethylbenzimidazole biosynthesis II (anaerobic)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Cobalamin-Biosynthesis|DMB-Biosynthesis
PWY-5523	5,6-dimethylbenzimidazole biosynthesis I (aerobic)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Cobalamin-Biosynthesis|DMB-Biosynthesis
PWY0-541	cyclopropane fatty acid (CFA) biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Cyclopropane-fatty-acid-biosyn
PWY-4942	sterculate biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Cyclopropane-fatty-acid-biosyn
PWY-7371	1,4-dihydroxy-6-naphthoate biosynthesis II	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|14-Dihydroxy-6-Naphthoate-Biosynthesis
PWY-7374	1,4-dihydroxy-6-naphthoate biosynthesis I	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|14-Dihydroxy-6-Naphthoate-Biosynthesis
PWY-6083	chlorobenzene degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorobenzene-Degradation
PWY-6099	1,2,4,5-tetrachlorobenzene degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorobenzene-Degradation
PWY-6081	1,3-dichlorobenzene degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorobenzene-Degradation
PWY-6091	1,2,4-trichlorobenzene degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorobenzene-Degradation
14DICHLORBENZDEG-PWY	1,4-dichlorobenzene degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorobenzene-Degradation
PWY-6090	1,2-dichlorobenzene degradation	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Chloroaromatic-Compounds-Degradation|Chlorobenzene-Degradation
PWY-7425	2-chloroacrylate degradation I	Degradation|CHLORINATED-COMPOUNDS-DEG|2-Chloroacrylates-Degradation
PWY-7428	2-chloroacrylate degradation II	Degradation|CHLORINATED-COMPOUNDS-DEG|2-Chloroacrylates-Degradation
PWY-5468	lupanine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|N-CONTAINING-SECONDARY-CMPD-SYN|ALKALOIDS-SYN|QUINOLIZIDINE-ALKALOIDS
PWY-5512	UDP-N-acetyl-D-galactosamine biosynthesis I	Biosynthesis|Polyamine-Biosynthesis|UDP-Nac-Galactosamine-Biosynthesis
PWY-5514	UDP-N-acetyl-D-galactosamine biosynthesis II	Biosynthesis|Polyamine-Biosynthesis|UDP-Nac-Galactosamine-Biosynthesis
PWY-6397	mycolyl-arabinogalactan-peptidoglycan complex biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|Cell-Wall-Biosynthesis
PWY-6626	cytidine-5-diphosphate-glycerol biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|Cell-Wall-Biosynthesis
PWY-7467	2-acetamido-4-amino-2,4,6-trideoxy-&alpha;-D-galactosyl-diphospho-ditrans,octacis-undecaprenol biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|Cell-Wall-Biosynthesis
PWY-7820	teichuronic acid biosynthesis (B. subtilis 168)	Biosynthesis|Cell-Structure-Biosynthesis|Cell-Wall-Biosynthesis
PWY-7817	type I lipoteichoic acid biosynthesis (S. aureus)	Biosynthesis|Cell-Structure-Biosynthesis|Cell-Wall-Biosynthesis|Teichoic-Acids-Biosynthesis
TEICHOICACID-PWY	poly(glycerol phosphate) wall teichoic acid biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|Cell-Wall-Biosynthesis|Teichoic-Acids-Biosynthesis
PWY-7816	poly(ribitol phosphate) wall teichoic acid biosynthesis II (S. aureus)	Biosynthesis|Cell-Structure-Biosynthesis|Cell-Wall-Biosynthesis|Teichoic-Acids-Biosynthesis
PWY-7819	poly(3-O-&beta;-D-glucopyranosyl-N-acetylgalactosamine 1-phosphate) wall teichoic acid biosynthesis	Biosynthesis|Cell-Structure-Biosynthesis|Cell-Wall-Biosynthesis|Teichoic-Acids-Biosynthesis
PWY-7815	poly(ribitol phosphate) wall teichoic acid biosynthesis I (B. subtilis)	Biosynthesis|Cell-Structure-Biosynthesis|Cell-Wall-Biosynthesis|Teichoic-Acids-Biosynthesis
PWY-7818	type IV lipoteichoic acid biosynthesis (S. pneumoniae)	Biosynthesis|Cell-Structure-Biosynthesis|Cell-Wall-Biosynthesis|Teichoic-Acids-Biosynthesis
PWY-7183	pyrimidine nucleobases salvage I	Biosynthesis|Nucleotide-Biosynthesis|PYR-NUC-SYN|Pyrimidine-Nucleotide-Salvage
PWY-7199	pyrimidine deoxyribonucleosides salvage	Biosynthesis|Nucleotide-Biosynthesis|PYR-NUC-SYN|Pyrimidine-Nucleotide-Salvage
PWY-7205	CMP phosphorylation	Biosynthesis|Nucleotide-Biosynthesis|PYR-NUC-SYN|Pyrimidine-Nucleotide-Salvage
PWY-7193	pyrimidine ribonucleosides salvage I	Biosynthesis|Nucleotide-Biosynthesis|PYR-NUC-SYN|Pyrimidine-Nucleotide-Salvage
PWY-6700	queuosine biosynthesis	Biosynthesis|Nucleotide-Biosynthesis|Nucleic-Acid-Processing
PWY-7375	mRNA capping I	Biosynthesis|Nucleotide-Biosynthesis|Nucleic-Acid-Processing
PWY0-1554	5-(carboxymethoxy)uridine biosynthesis	Biosynthesis|Nucleotide-Biosynthesis|Nucleic-Acid-Processing
PWY-7286	7-(3-amino-3-carboxypropyl)-wyosine biosynthesis	Biosynthesis|Nucleotide-Biosynthesis|Nucleic-Acid-Processing
PWY0-1479	tRNA processing	Biosynthesis|Nucleotide-Biosynthesis|Nucleic-Acid-Processing
PWY-6711	archaeosine biosynthesis	Biosynthesis|Nucleotide-Biosynthesis|Nucleic-Acid-Processing
PWY-7285	methylwyosine biosynthesis	Biosynthesis|Nucleotide-Biosynthesis|Nucleic-Acid-Processing
PWY0-1587	N<sup>6</sup>-L-threonylcarbamoyladenosine<sup>37</sup>-modified tRNA biosynthesis	Biosynthesis|Nucleotide-Biosynthesis|Nucleic-Acid-Processing
PWY-7803	tRNA splicing II	Biosynthesis|Nucleotide-Biosynthesis|Nucleic-Acid-Processing|tRNA-Splicing
PWY-6689	tRNA splicing I	Biosynthesis|Nucleotide-Biosynthesis|Nucleic-Acid-Processing|tRNA-Splicing
PWY-2681	trans-zeatin biosynthesis	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|CYTOKININ-BIOSYNTHESIS
PWY-2901	cytokinins 9-N-glucoside biosynthesis	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|CYTOKININ-BIOSYNTHESIS
PWY-2881	cytokinins 7-N-glucoside biosynthesis	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|CYTOKININ-BIOSYNTHESIS
PWY-5967	lupinate biosynthesis	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|CYTOKININ-BIOSYNTHESIS
PWY-2781	cis-zeatin biosynthesis	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|CYTOKININ-BIOSYNTHESIS
PWY-6794	adenosine 5-phosphoramidate biosynthesis	Biosynthesis|Nucleotide-Biosynthesis|PUR-NUC-SYN
PWY-7219	adenosine ribonucleotides de novo biosynthesis	Biosynthesis|Nucleotide-Biosynthesis|PUR-NUC-SYN|Purine-Nucleotide-De-Novo-Biosynthesis|Purine-Ribonuc-De-Novo-Biosynthesis
PWY-7221	guanosine ribonucleotides de novo biosynthesis	Biosynthesis|Nucleotide-Biosynthesis|PUR-NUC-SYN|Purine-Nucleotide-De-Novo-Biosynthesis|Purine-Ribonuc-De-Novo-Biosynthesis
PWY-5194	siroheme biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Porphyrin-Compounds-Biosynthesis
PWY-5297	siroheme amide biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Porphyrin-Compounds-Biosynthesis
PWY-7760	bacteriochlorophyll e biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Porphyrin-Compounds-Biosynthesis|Chlorophyll-Biosynthesis
PWY-7759	bacteriochlorophyll c biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Porphyrin-Compounds-Biosynthesis|Chlorophyll-Biosynthesis
PWY-7758	bacteriochlorophyll d biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Porphyrin-Compounds-Biosynthesis|Chlorophyll-Biosynthesis
PWY-7762	bacteriochlorophyll b biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Porphyrin-Compounds-Biosynthesis|Chlorophyll-Biosynthesis
PWY-5526	bacteriochlorophyll a biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Porphyrin-Compounds-Biosynthesis|Chlorophyll-Biosynthesis
PWY-7159	3,8-divinyl-chlorophyllide a biosynthesis III (aerobic, light independent)	Biosynthesis|Cofactor-Biosynthesis|Porphyrin-Compounds-Biosynthesis|Chlorophyll-Biosynthesis|Chlorophyllide-a-Biosynthesis
PWY-5531	3,8-divinyl-chlorophyllide a biosynthesis II (anaerobic)	Biosynthesis|Cofactor-Biosynthesis|Porphyrin-Compounds-Biosynthesis|Chlorophyll-Biosynthesis|Chlorophyllide-a-Biosynthesis
CHLOROPHYLL-SYN	3,8-divinyl-chlorophyllide a biosynthesis I (aerobic, light-dependent)	Biosynthesis|Cofactor-Biosynthesis|Porphyrin-Compounds-Biosynthesis|Chlorophyll-Biosynthesis|Chlorophyllide-a-Biosynthesis
PWY-6423	hemoglobin degradation	Degradation|Protein-Degradation
ORNDEG-PWY	superpathway of ornithine degradation	Super-Pathways
PWY-5736	isopropylamine degradation	Degradation|AMINE-DEG
PWY-6534	phenylethylamine degradation II	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Phenolic-Compounds-Degradation|Phenylethylamine-Degradation
PWY-7085	triethylamine degradation	Degradation|AMINE-DEG
PWY6666-2	dopamine degradation	Degradation|AMINE-DEG
GLCMANNANAUT-PWY	superpathway of N-acetylglucosamine, N-acetylmannosamine and N-acetylneuraminate degradation	Super-Pathways
PWY-3621	&gamma;-butyrobetaine degradation	Super-Pathways
PWY-6181	histamine degradation	Degradation|AMINE-DEG
PWY-6968	trimethylamine degradation	Degradation|AMINE-DEG
PWY0-1477	ethanolamine utilization	Degradation|AMINE-DEG
2PHENDEG-PWY	phenylethylamine degradation I	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Phenolic-Compounds-Degradation|Phenylethylamine-Degradation
P561-PWY	L-proline betaine degradation	Degradation|AMINE-DEG
PWY-6071	superpathway of phenylethylamine degradation	Super-Pathways
PWY-6962	superpathway of trimethylamine degradation	Super-Pathways
PWY-7431	aromatic biogenic amine degradation (bacteria)	Degradation|AMINE-DEG
PWY-5373	calendate biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis|Conjugated-Fatty-Acid-Biosynthesis
PWY-7590	(7Z,10Z,13Z)-hexadecatrienoate biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis|Conjugated-Fatty-Acid-Biosynthesis
CYCLOHEXANOL-OXIDATION-PWY	cyclohexanol degradation	Degradation|Alcohol-Degradation
PWY-5368	dimorphecolate biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis|Conjugated-Fatty-Acid-Biosynthesis
GLYCOL-GLYOXDEG-PWY	superpathway of glycol metabolism and degradation	Super-Pathways
PWY-6464	polyvinyl alcohol degradation	Degradation|Alcohol-Degradation
PWY-5375	&alpha;-eleostearate biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis|Conjugated-Fatty-Acid-Biosynthesis
PWY-7013	L-1,2-propanediol degradation	Degradation|Alcohol-Degradation
PWY-7616	methanol oxidation to carbon dioxide	Super-Pathways
PWY-5374	punicate biosynthesis	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis|Conjugated-Fatty-Acid-Biosynthesis
PWY0-1280	ethylene glycol degradation	Degradation|Alcohol-Degradation
PWY66-389	phytol degradation	Degradation|Alcohol-Degradation
PWY-7598	&alpha;-linolenate biosynthesis II (cyanobacteria)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|Alpha-Linolenate-Biosynthesis
MGLDLCTANA-PWY	methylglyoxal degradation VI	Detoxification|Methylglyoxal-Detoxification
PWY-5997	&alpha;-linolenate biosynthesis I (plants and red algae)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|Alpha-Linolenate-Biosynthesis
PWY-5456	methylglyoxal degradation VII	Detoxification|Methylglyoxal-Detoxification
PWY-5462	methylglyoxal degradation II	Detoxification|Methylglyoxal-Detoxification
METHGLYUT-PWY	superpathway of methylglyoxal degradation	Super-Pathways
PWY-5453	methylglyoxal degradation III	Detoxification|Methylglyoxal-Detoxification
PWY-5459	methylglyoxal degradation IV	Super-Pathways
PWY-5386	methylglyoxal degradation I	Detoxification|Methylglyoxal-Detoxification
PWY-5458	methylglyoxal degradation V	Detoxification|Methylglyoxal-Detoxification
PWY-7754	bile acids degradation	Degradation|Bile-Acids-Degradation
PWY-6518	glycocholate metabolism (bacteria)	Degradation|Bile-Acids-Degradation
CO2FORM-PWY	methanogenesis from methanol	Energy-Metabolism|Respiration|ANAEROBIC-RESPIRATION|METHANOGENESIS
PWY-1722	formate assimilation into 5,10-methylenetetrahydrofolate	Degradation|C1-COMPOUNDS
PWY-1881	formate oxidation to CO2	Energy-Metabolism|CHEMOAUTOTROPHIC-ENERGY-METABOLISM
PWY-6183	salicylate degradation I	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Salicylate-Degradation
PWY-1882	superpathway of C1 compounds oxidation to CO2	Super-Pathways
PWY-7750	carbon monoxide oxidation to CO2	Degradation|C1-COMPOUNDS
PWY-6640	salicylate degradation IV	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Salicylate-Degradation
PWY-2582	brassinosteroid biosynthesis II	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|Brassinosteroid-Biosynthesis
PWY-6636	salicylate degradation III	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Salicylate-Degradation
PWY-699	brassinosteroid biosynthesis I	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|Brassinosteroid-Biosynthesis
PWY-6224	salicylate degradation II	Degradation|AROMATIC-COMPOUNDS-DEGRADATION|Salicylate-Degradation
ACETATEUTIL-PWY	superpathway of acetate utilization and formation	Super-Pathways
GLUCONSUPER-PWY	D-gluconate degradation	Degradation|CARBOXYLATES-DEG
P441-PWY	superpathway of N-acetylneuraminate degradation	Super-Pathways
PWY-2361	3-oxoadipate degradation	Degradation|CARBOXYLATES-DEG
PWY-6665	pterostilbene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|STILBENE-SYN
PWY-301	cyclohexane-1-carboxylate degradation (anaerobic)	Degradation|CARBOXYLATES-DEG
PWY-5074	mevalonate degradation	Degradation|CARBOXYLATES-DEG
PWY-5162	2-oxopentenoate degradation	Degradation|CARBOXYLATES-DEG
PWY-5177	glutaryl-CoA degradation	Degradation|CARBOXYLATES-DEG
PWY-84	resveratrol biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|STILBENE-SYN
PWY-5652	2-amino-3-carboxymuconate semialdehyde degradation to glutaryl-CoA	Degradation|CARBOXYLATES-DEG
PWY-5654	2-amino-3-carboxymuconate semialdehyde degradation to 2-oxopentenoate	Degradation|CARBOXYLATES-DEG
PWY-5749	itaconate degradation	Degradation|CARBOXYLATES-DEG
PWY-6021	nitrilotriacetate degradation	Degradation|CARBOXYLATES-DEG
PWY-7736	stellatic acid biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESTERTERPENOID-BIOSYNTHESIS
PWY-6038	citrate degradation	Degradation|CARBOXYLATES-DEG
PWY-6373	acrylate degradation	Degradation|CARBOXYLATES-DEG
PWY-7310	D-glucosaminate degradation	Degradation|CARBOXYLATES-DEG
PWY-804	glycolate degradation II	Energy-Metabolism|Fermentation
PWY-7720	ophiobolin F biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESTERTERPENOID-BIOSYNTHESIS
PWY0-1313	acetate conversion to acetyl-CoA	Degradation|CARBOXYLATES-DEG
PWY-6643	coenzyme M biosynthesis II	Biosynthesis|Cofactor-Biosynthesis|Coenzyme-M-Biosynthesis
PWY0-1465	D-malate degradation	Degradation|CARBOXYLATES-DEG
P261-PWY	coenzyme M biosynthesis I	Biosynthesis|Cofactor-Biosynthesis|Coenzyme-M-Biosynthesis
PYRUVDEHYD-PWY	pyruvate decarboxylation to acetyl CoA	Energy-Metabolism|Acetyl-CoA-Biosynthesis
PWY-6540	costunolide biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN|SESQUITERPENE-LACTONE
PWY-7559	glutathione degradation (DUG pathway - yeast)	Degradation|COFACTOR-DEGRADATION
PWY-6370	ascorbate recycling (cytosolic)	Degradation|COFACTOR-DEGRADATION
PWY-5195	artemisinin biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN|SESQUITERPENE-LACTONE
PWY-6028	acetoin degradation	Degradation|Carbohydrates-Degradation
PWY-7779	methyl tert-butyl ether degradation	Degradation|Carbohydrates-Degradation
PWY-7497	3&beta;-hydroxysesquiterpene lactone biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN|SESQUITERPENE-LACTONE
PWY-7780	butane degradation	Degradation|Carbohydrates-Degradation
ACETOACETATE-DEG-PWY	acetoacetate degradation (to acetyl CoA)	Degradation|Fatty-Acid-and-Lipid-Degradation
LIPAS-PWY	triacylglycerol degradation	Degradation|Fatty-Acid-and-Lipid-Degradation
LIPASYN-PWY	phospholipases	Metabolic-Clusters
PWY-4261	glycerol degradation I	Degradation|Alcohol-Degradation|GLYCEROL-DEG
PWY-6111	mitochondrial L-carnitine shuttle	Degradation|Fatty-Acid-and-Lipid-Degradation
PWY-6368	3-phosphoinositide degradation	Degradation|Fatty-Acid-and-Lipid-Degradation
PWY-6483	ceramide degradation	Degradation|Fatty-Acid-and-Lipid-Degradation
PWY-7119	sphingolipid recycling and degradation (yeast)	Degradation|Fatty-Acid-and-Lipid-Degradation
PWY-7366	phosphatidylglycerol degradation	Degradation|Fatty-Acid-and-Lipid-Degradation
PWY-7367	phosphatidylcholine resynthesis via glycerophosphocholine	Degradation|Fatty-Acid-and-Lipid-Degradation
PWY-7783	plasmalogen degradation	Degradation|Fatty-Acid-and-Lipid-Degradation
PWY3DJ-11470	sphingosine and sphingosine-1-phosphate metabolism	Degradation|Fatty-Acid-and-Lipid-Degradation
PWY6666-1	anandamide degradation	Degradation|Fatty-Acid-and-Lipid-Degradation
SOPHOROSYLOXYDOCOSANOATE-DEG-PWY	sophorosyloxydocosanoate deacetylation	Degradation|Fatty-Acid-and-Lipid-Degradation
PWY-6313	serotonin degradation	Degradation|HORMONE-DEG
PWY-6342	noradrenaline and adrenaline degradation	Degradation|HORMONE-DEG
PWY-6688	thyronamine and iodothyronamine metabolism	Degradation|HORMONE-DEG
PWY-5706	alliin metabolism	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|S-CONTAINING-SECONDARY-CMPD-SYN
PWY-5532	nucleoside and nucleotide degradation (archaea)	Degradation|NUCLEO-DEG
PWY-7614	methiin metabolism	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|S-CONTAINING-SECONDARY-CMPD-SYN
PWY0-1297	superpathway of purine deoxyribonucleosides degradation	Super-Pathways
PWY-5389	3-methylthiopropanoate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|S-CONTAINING-SECONDARY-CMPD-SYN
PWY-5708	ethiin metabolism	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|S-CONTAINING-SECONDARY-CMPD-SYN
ASPSYNII-PWY	cyanide detoxification I	Detoxification|Cyanide-Detoxification
PWY-6900	(Z)-butanethial-S-oxide biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|S-CONTAINING-SECONDARY-CMPD-SYN
P161-PWY	acetylene degradation	Degradation|Other-Degradation
P201-PWY	nitroglycerin degradation	Degradation|Other-Degradation
P221-PWY	octane oxidation	Degradation|Other-Degradation
P401-PWY	cyanide degradation	Degradation|Other-Degradation
PWY-5331	taurine biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|S-CONTAINING-SECONDARY-CMPD-SYN
P481-PWY	adamantanone degradation	Degradation|Other-Degradation
P482-PWY	arsonoacetate degradation	Degradation|Other-Degradation
P621-PWY	nylon-6 oligomer degradation	Degradation|Other-Degradation
PWY-4061	glutathione-mediated detoxification I	Degradation|Other-Degradation
PWY-5707	propanethial S-oxide biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|S-CONTAINING-SECONDARY-CMPD-SYN
PWY-5355	nitroethane degradation	Degradation|Other-Degradation
PWY-5744	glyoxylate assimilation	Degradation|C1-COMPOUNDS|CO2-Fixation|Autotrophic-CO2-Fixation
PWY-6377	&alpha;-tocopherol degradation	Degradation|Other-Degradation
PWY-7112	4-hydroxy-2-nonenal detoxification	Degradation|Other-Degradation
PWY-6539	(Z)-phenylmethanethial S-oxide biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|S-CONTAINING-SECONDARY-CMPD-SYN
PWY-723	alkylnitronates degradation	Degradation|Other-Degradation
PWY-7582	mercaptosuccinate degradation	Degradation|Other-Degradation
PWY0-1569	autoinducer AI-2 degradation	Degradation|Other-Degradation
PWY66-241	bupropion degradation	Degradation|Other-Degradation
PWY-6341	guaiacylglycerol-&beta;-guaiacyl ether degradation	Degradation|Polymer-Degradation
PWY-7794	polyethylene terephthalate degradation	Degradation|Polymer-Degradation
PWY0-1546	muropeptide degradation	Degradation|Polymer-Degradation
PWY-6978	plastoquinol-9 biosynthesis II	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Plastoquinone-Biosynthesis
PWY-1581	plastoquinol-9 biosynthesis I	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Plastoquinone-Biosynthesis
PWY-6813	glucuronoarabinoxylan degradation	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG
PWY-6815	porphyran degradation	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG
PWY-6816	agarose degradation	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG
PWY-6827	gellan degradation	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG
PWY-6986	alginate degradation	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG
PWY-7456	mannan degradation	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG
PWY-6427	rot-2-enonate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ISOFLAVONOID-SYN|Rotenoids-Biosynthesis
PWY-7647	ulvan degradation	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG
PWY-7828	cyclobis-(1&rarr;6)-&alpha;-nigerosyl degradation	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG
PWY-862	fructan degradation	Degradation|Carbohydrates-Degradation|POLYSACCHARIDES-DEG
PWY-7308	acrylonitrile degradation I	Degradation|Other-Degradation|Acrylonitrile-Degradation
PWY-6425	rotenoid biosynthesis II	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ISOFLAVONOID-SYN|Rotenoids-Biosynthesis
P344-PWY	superpathway of acrylonitrile degradation	Super-Pathways
PWY-7309	acrylonitrile degradation II	Degradation|Other-Degradation|Acrylonitrile-Degradation
PWY-5336	carbon disulfide oxidation II (aerobic)	Degradation|Other-Degradation|Carbon-disulfide-degradation
PWY-1164	carbon disulfide oxidation I (anaerobic)	Degradation|Other-Degradation|Carbon-disulfide-degradation
PWY-5775	rotenoid biosynthesis I	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|PHENYLPROPANOID-SYN|FLAVONOID-SYN|ISOFLAVONOID-SYN|Rotenoids-Biosynthesis
PWY-7128	nicotine degradation III (VPP pathway)	Degradation|Other-Degradation|Nicotine-Degradation
PWY-6993	nicotine degradation II (pyrrolidine pathway)	Degradation|Other-Degradation|Nicotine-Degradation
PWY66-221	nicotine degradation V	Degradation|Other-Degradation|Nicotine-Degradation
P181-PWY	nicotine degradation I (pyridine pathway)	Degradation|Other-Degradation|Nicotine-Degradation
PWY66-201	nicotine degradation IV	Degradation|Other-Degradation|Nicotine-Degradation
PWY-743	thiocyanate degradation II	Degradation|Other-Degradation|Thiocyanate-Degradation
P581-PWY	thiocyanate degradation I	Degradation|Other-Degradation|Thiocyanate-Degradation
PWY-6692	Fe(II) oxidation	Degradation|Noncarbon-Nutrients|Iron-Metabolism
PWY-5934	iron reduction and absorption	Degradation|Noncarbon-Nutrients|Iron-Metabolism
PWY-6592	manganese oxidation II	Degradation|Noncarbon-Nutrients|Manganese-Oxidation
PWY-6739	pinitol biosynthesis II	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|Sugar-Alcohols-Biosynthesis|Pinitol-Biosynthesis
PWY-6591	manganese oxidation I	Degradation|Noncarbon-Nutrients|Manganese-Oxidation
PWY-6738	pinitol biosynthesis I	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|Sugar-Alcohols-Biosynthesis|Pinitol-Biosynthesis
CYANCAT-PWY	cyanate degradation	Degradation|Noncarbon-Nutrients|NITROGEN-DEG
PWY-5063	phytyl diphosphate biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|DITERPENOID-SYN|Phytyl-Diphosphate-Biosynthesis
PWY-4984	urea cycle	Degradation|Noncarbon-Nutrients|NITROGEN-DEG
PWY-7387	hypotaurine degradation	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism
PWY-7683	nitrite reduction (hemoglobin)	Degradation|Noncarbon-Nutrients|NITROGEN-DEG
P483-PWY	phosphonoacetate degradation	Degradation|Noncarbon-Nutrients|Phosphorus-Compounds
PWY-5491	diethylphosphate degradation	Degradation|Noncarbon-Nutrients|Phosphorus-Compounds
PWY-6348	phosphate acquisition	Degradation|Noncarbon-Nutrients|Phosphorus-Compounds
PWY-7805	aminomethylphosphonate degradation	Degradation|Noncarbon-Nutrients|Phosphorus-Compounds
PWY-6932	selenate reduction	Degradation|Noncarbon-Nutrients|Selenium-Metabolism
PWY-5889	rhodoquinone-9 biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Rhodoquinone-Biosynthesis
ALKANEMONOX-PWY	two-component alkanesulfonate monooxygenase	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism
PWY-5123	trans, trans-farnesyl diphosphate biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Polyprenyl-Biosynthesis|All-Trans-Farnesyl-PP-Biosynthesis
PWY-2601	isethionate degradation	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism
PWY-5888	rhodoquinone-10 biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Quinone-Biosynthesis|Rhodoquinone-Biosynthesis
PWY-5258	methanogenesis from dimethylsulfide	Energy-Metabolism|Respiration|ANAEROBIC-RESPIRATION|METHANOGENESIS
PWY-5260	methanogenesis from methylthiopropanoate	Energy-Metabolism|Respiration|ANAEROBIC-RESPIRATION|METHANOGENESIS
PWY-5340	sulfate activation for sulfonation	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism
PWY-6043	ethanedisulfonate degradation	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism
PWY-6044	methanesulfonate degradation	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism
PWY-6050	dimethyl sulfoxide degradation	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism
PWY-6058	dimethyl sulfone degradation	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism
PWY-6072	superpathway of dimethylsulfone degradation	Super-Pathways
PWY-6321	homotaurine degradation	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism
PWY-6327	tetrathionate oxidation	Energy-Metabolism|Respiration|AEROBIC-RESPIRATION|Tetrathionate-Oxidation
PWY-6593	sulfoacetate degradation	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism
PWY-6634	3-sulfopropanediol degradation	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism
PWY-6642	(R)-cysteate degradation	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism
PWY-6736	sulfur volatiles biosynthesis	Metabolic-Clusters
PWY-7462	3,3-dithiodipropanoate degradation	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism
PWY-7465	3,3-thiodipropanoate degradation	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism
PWY-7793	dimethyl sulfide biosynthesis from methionine	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism
PWY-7027	hentriaconta-3,6,9,12,15,19,22,25,28-nonaene biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|Olefins-Biosynthesis
SER-GLYSYN-PWY	superpathway of L-serine and glycine biosynthesis I	Super-Pathways
PWY-7026	terminal olefins biosynthesis I	Biosynthesis|Carbohydrates-Biosynthesis|Olefins-Biosynthesis|Alkenes-Biosynthesis
ASPASN-PWY	superpathway of L-aspartate and L-asparagine biosynthesis	Super-Pathways
PWY-7035	(Z)-9-tricosene biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|Olefins-Biosynthesis|Alkenes-Biosynthesis
P4-PWY	superpathway of L-lysine, L-threonine and L-methionine biosynthesis I	Super-Pathways
PWY-7029	terminal olefins biosynthesis II	Biosynthesis|Carbohydrates-Biosynthesis|Olefins-Biosynthesis|Alkenes-Biosynthesis
PWY-821	superpathway of sulfur amino acid biosynthesis (Saccharomyces cerevisiae)	Super-Pathways
PWY-7032	alkane biosynthesis I	Biosynthesis|Carbohydrates-Biosynthesis|Alkanes-Biosynthesis
COMPLETE-ARO-PWY	superpathway of aromatic amino acid biosynthesis	Super-Pathways
PWY-6622	heptadecane biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|Alkanes-Biosynthesis
PWY-724	superpathway of L-lysine, L-threonine and L-methionine biosynthesis II	Super-Pathways
PWY-7033	alkane biosynthesis II	Biosynthesis|Carbohydrates-Biosynthesis|Alkanes-Biosynthesis
BRANCHED-CHAIN-AA-SYN-PWY	superpathway of branched chain amino acid biosynthesis	Super-Pathways
PWY-6122	5-aminoimidazole ribonucleotide biosynthesis II	Biosynthesis|Nucleotide-Biosynthesis|PUR-NUC-SYN|AIR-Biosynthesis
PWY-3481	superpathway of L-phenylalanine and L-tyrosine biosynthesis	Super-Pathways
PWY-6121	5-aminoimidazole ribonucleotide biosynthesis I	Biosynthesis|Nucleotide-Biosynthesis|PUR-NUC-SYN|AIR-Biosynthesis
PWY-7380	biotin biosynthesis from 8-amino-7-oxononanoate II	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|BIOTIN-SYN
PWY0-1507	biotin biosynthesis from 8-amino-7-oxononanoate I	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|BIOTIN-SYN
PWY-7234	inosine-5-phosphate biosynthesis III	Biosynthesis|Nucleotide-Biosynthesis|PUR-NUC-SYN|Purine-Nucleotide-De-Novo-Biosynthesis|Purine-Ribonuc-De-Novo-Biosynthesis|IMP-Biosynthesis
PWY-6124	inosine-5-phosphate biosynthesis II	Biosynthesis|Nucleotide-Biosynthesis|PUR-NUC-SYN|Purine-Nucleotide-De-Novo-Biosynthesis|Purine-Ribonuc-De-Novo-Biosynthesis|IMP-Biosynthesis
PWY-6123	inosine-5-phosphate biosynthesis I	Biosynthesis|Nucleotide-Biosynthesis|PUR-NUC-SYN|Purine-Nucleotide-De-Novo-Biosynthesis|Purine-Ribonuc-De-Novo-Biosynthesis|IMP-Biosynthesis
PWY-7187	pyrimidine deoxyribonucleotides de novo biosynthesis II	Biosynthesis|Nucleotide-Biosynthesis|Deoxyribonucleotide-Biosynthesis|Pyrimid-Deoxyribonucleot-De-Novo-Biosyn
PWY-7176	UTP and CTP de novo biosynthesis	Biosynthesis|Nucleotide-Biosynthesis|PYR-NUC-SYN|Pyrimidine-De-Novo-Biosyn|Pyrimid-Ribonucleot-De-Novo-Biosyn
PWY-7357	thiamine formation from pyrithiamine and oxythiamine (yeast)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Thiamine-Biosynthesis
PWY-7369	thiamine triphosphate metabolism	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Thiamine-Biosynthesis
PWY-6896	thiamine salvage I	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Thiamine-Biosynthesis|Thiamin-Salvage
PWY-6899	base-degraded thiamine salvage	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Thiamine-Biosynthesis|Thiamin-Salvage
PWY-6898	thiamine salvage III	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Thiamine-Biosynthesis|Thiamin-Salvage
PWY-7356	thiamine salvage IV (yeast)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Thiamine-Biosynthesis|Thiamin-Salvage
PWY-6910	hydroxymethylpyrimidine salvage	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Thiamine-Biosynthesis|Thiamin-Salvage
PWY-7282	4-amino-2-methyl-5-diphosphomethylpyrimidine biosynthesis (yeast)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Thiamine-Biosynthesis|HMP-PP-Biosynthesis
PWY-6890	|4-amino-2-methyl-5-diphosphomethylpyrimidine biosynthesis|	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Thiamine-Biosynthesis|HMP-PP-Biosynthesis
PWY-6894	thiamine diphosphate biosynthesis I (E. coli)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Thiamine-Biosynthesis|Thiamin-Diphosphate-Biosynthesis
ADENOSYLHOMOCYSCAT-PWY	L-methionine salvage from L-homocysteine	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|METHIONINE-SYN|Methionine-Salvage
PWY-7527	L-methionine salvage cycle III	Super-Pathways
PWY-6893	thiamine diphosphate biosynthesis II (Bacillus)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Thiamine-Biosynthesis|Thiamin-Diphosphate-Biosynthesis
PWY-7270	L-methionine salvage cycle II (plants)	Super-Pathways
PWY-6908	thiamine diphosphate biosynthesis IV (eukaryotes)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Thiamine-Biosynthesis|Thiamin-Diphosphate-Biosynthesis
PWY-I9	L-cysteine biosynthesis VI (from L-methionine)	Super-Pathways
PWY-6907	thiamine diphosphate biosynthesis III (Staphylococcus)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Thiamine-Biosynthesis|Thiamin-Diphosphate-Biosynthesis
PWY-5441	S-methyl-L-methionine cycle	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|METHIONINE-SYN|Methionine-Salvage
PWY-6163	chorismate biosynthesis from 3-dehydroquinate	Biosynthesis|AROMATIC-COMPOUNDS-BIOSYN|Chorismate-Biosynthesis
PWY-7528	L-methionine salvage cycle I (bacteria and plants)	Super-Pathways
PWY-7174	S-methyl-5-thio-&alpha;-D-ribose 1-phosphate degradation II	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|METHIONINE-SYN|Methionine-Salvage|MTR-1P-Degradation
PWY-4361	S-methyl-5-thio-&alpha;-D-ribose 1-phosphate degradation	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|METHIONINE-SYN|Methionine-Salvage|MTR-1P-Degradation
PWY-1822	indole-3-acetate activation I	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|Auxin-Biosynthesis
PWY-6754	S-methyl-5-thioadenosine degradation I	Degradation|NUCLEO-DEG|Methylthioadenosine-Degradation
PWY-6753	S-methyl-5-thioadenosine degradation III	Degradation|NUCLEO-DEG|Methylthioadenosine-Degradation
PWY-4441	DIMBOA-glucoside activation	Activation-Inactivation-Interconversion|Activation
PWY-1921	indole-3-acetate activation II	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|Auxin-Biosynthesis
PWY0-1391	S-methyl-5-thioadenosine degradation IV	Degradation|NUCLEO-DEG|Methylthioadenosine-Degradation
PWY-6756	S-methyl-5-thioadenosine degradation II	Degradation|NUCLEO-DEG|Methylthioadenosine-Degradation
COA-PWY	coenzyme A biosynthesis I	Biosynthesis|Cofactor-Biosynthesis|CoA-Biosynthesis
COA-PWY-1	coenzyme A biosynthesis II (mammalian)	Biosynthesis|Cofactor-Biosynthesis|CoA-Biosynthesis
P164-PWY	purine nucleobases degradation I (anaerobic)	Degradation|NUCLEO-DEG|Purine-Degradation
PWY-5044	purine nucleotides degradation I (plants)	Super-Pathways
PWY-5497	purine nucleobases degradation II (anaerobic)	Degradation|NUCLEO-DEG|Purine-Degradation
PWY-5695	urate biosynthesis/inosine 5-phosphate degradation	Biosynthesis|Polyamine-Biosynthesis
PWY-6019	pseudouridine degradation	Degradation|NUCLEO-DEG|Purine-Degradation
PWY-6353	purine nucleotides degradation II (aerobic)	Super-Pathways
PWY-7179	purine deoxyribonucleosides degradation I	Degradation|NUCLEO-DEG|Purine-Degradation
PWY-7179-1	purine deoxyribonucleosides degradation II	Degradation|NUCLEO-DEG|Purine-Degradation
PWY0-1296	purine ribonucleosides degradation	Degradation|NUCLEO-DEG|Purine-Degradation
ARGSPECAT-PWY	spermine biosynthesis	Biosynthesis|Polyamine-Biosynthesis
P101-PWY	ectoine biosynthesis	Biosynthesis|Polyamine-Biosynthesis
PWY-5907	homospermidine biosynthesis	Biosynthesis|Polyamine-Biosynthesis
PWY-6173	histamine biosynthesis	Biosynthesis|Polyamine-Biosynthesis
PWY0-1298	superpathway of pyrimidine deoxyribonucleosides degradation	Super-Pathways
PWY-6456	serinol biosynthesis	Biosynthesis|Polyamine-Biosynthesis
PWY-6400	melatonin degradation III	Degradation|HORMONE-DEG|Melatonin-Degradation
PWY-6562	norspermidine biosynthesis	Biosynthesis|Polyamine-Biosynthesis
PWY-6399	melatonin degradation II	Degradation|HORMONE-DEG|Melatonin-Degradation
PWY-7297	octopamine biosynthesis	Biosynthesis|Polyamine-Biosynthesis
PWY-6398	melatonin degradation I	Degradation|HORMONE-DEG|Melatonin-Degradation
PWY0-1303	aminopropylcadaverine biosynthesis	Biosynthesis|Polyamine-Biosynthesis
PWY-6402	superpathway of melatonin degradation	Super-Pathways
TRYPANOSYN-PWY	trypanothione biosynthesis	Biosynthesis|Polyamine-Biosynthesis
PWY-6834	spermidine biosynthesis III	Biosynthesis|Polyamine-Biosynthesis|Spermidine-Biosynthesis
PWY-6559	spermidine biosynthesis II	Biosynthesis|Polyamine-Biosynthesis|Spermidine-Biosynthesis
BSUBPOLYAMSYN-PWY	spermidine biosynthesis I	Biosynthesis|Polyamine-Biosynthesis|Spermidine-Biosynthesis
PWY-6546	brassinosteroids inactivation	Metabolic-Clusters
PWY-6261	thyroid hormone metabolism II (via conjugation and/or degradation)	Degradation|HORMONE-DEG|Thyroid-Hormone-Metabolism
PWY-6260	thyroid hormone metabolism I (via deiodination)	Degradation|HORMONE-DEG|Thyroid-Hormone-Metabolism
PWY-5533	acetone degradation II (to acetoacetate)	Degradation|Fatty-Acid-and-Lipid-Degradation|Acetone-Degradation
PWY-6823	molybdenum cofactor biosynthesis	Biosynthesis|Cofactor-Biosynthesis|Molybdenum-Cofactor-Biosynthesis
PWY-5451	acetone degradation I (to methylglyoxal)	Degradation|Fatty-Acid-and-Lipid-Degradation|Acetone-Degradation
PWY-7466	acetone degradation III (to propane-1,2-diol)	Degradation|Fatty-Acid-and-Lipid-Degradation|Acetone-Degradation
PWY-7337	10-cis-heptadecenoyl-CoA degradation (yeast)	Degradation|Fatty-Acid-and-Lipid-Degradation|Fatty-Acid-Degradation
PWY-7340	9-cis, 11-trans-octadecadienoyl-CoA degradation (isomerase-dependent, yeast)	Degradation|Fatty-Acid-and-Lipid-Degradation|Fatty-Acid-Degradation
PWY-2501	fatty acid &alpha;-oxidation I	Degradation|Fatty-Acid-and-Lipid-Degradation|Fatty-Acid-Degradation
PWY66-388	fatty acid &alpha;-oxidation III	Degradation|Fatty-Acid-and-Lipid-Degradation|Fatty-Acid-Degradation
PWY-5137	fatty acid &beta;-oxidation III (unsaturated, odd number)	Degradation|Fatty-Acid-and-Lipid-Degradation|Fatty-Acid-Degradation
PWY-7288	fatty acid &beta;-oxidation (peroxisome, yeast)	Degradation|Fatty-Acid-and-Lipid-Degradation|Fatty-Acid-Degradation
PWY-7307	oleate &beta;-oxidation (reductase-dependent, yeast)	Degradation|Fatty-Acid-and-Lipid-Degradation|Fatty-Acid-Degradation
PWY-7339	10-trans-heptadecenoyl-CoA degradation (MFE-dependent, yeast)	Degradation|Fatty-Acid-and-Lipid-Degradation|Fatty-Acid-Degradation
NADPHOS-DEPHOS-PWY	NAD phosphorylation and dephosphorylation	Biosynthesis|Cofactor-Biosynthesis|NAD-Metabolism
FAO-PWY	fatty acid &beta;-oxidation I	Degradation|Fatty-Acid-and-Lipid-Degradation|Fatty-Acid-Degradation
PWY-5381	pyridine nucleotide cycling (plants)	Biosynthesis|Cofactor-Biosynthesis|NAD-Metabolism
PWY66-387	fatty acid &alpha;-oxidation II	Degradation|Fatty-Acid-and-Lipid-Degradation|Fatty-Acid-Degradation
PWY-7269	NAD/NADP-NADH/NADPH mitochondrial interconversion (yeast)	Biosynthesis|Cofactor-Biosynthesis|NAD-Metabolism
PWY-5136	fatty acid &beta;-oxidation II (peroxisome)	Degradation|Fatty-Acid-and-Lipid-Degradation|Fatty-Acid-Degradation
PWY-5083	NAD/NADH phosphorylation and dephosphorylation	Biosynthesis|Cofactor-Biosynthesis|NAD-Metabolism
PWY-6837	fatty acid beta-oxidation V (unsaturated, odd number, di-isomerase-dependent)	Degradation|Fatty-Acid-and-Lipid-Degradation|Fatty-Acid-Degradation
PWY-7268	NAD/NADP-NADH/NADPH cytosolic interconversion (yeast)	Biosynthesis|Cofactor-Biosynthesis|NAD-Metabolism
PWY-7292	oleate &beta;-oxidation (thioesterase-dependent, yeast)	Degradation|Fatty-Acid-and-Lipid-Degradation|Fatty-Acid-Degradation
NADPHOS-DEPHOS-PWY-1	NAD phosphorylation and transhydrogenation	Biosynthesis|Cofactor-Biosynthesis|NAD-Metabolism
PWY-7338	10-trans-heptadecenoyl-CoA degradation (reductase-dependent, yeast)	Degradation|Fatty-Acid-and-Lipid-Degradation|Fatty-Acid-Degradation
PWY-6387	UDP-N-acetylmuramoyl-pentapeptide biosynthesis I (meso-diaminopimelate containing)	Biosynthesis|Cell-Structure-Biosynthesis|Cell-Wall-Biosynthesis|Acetylmuramoyl-Pentapeptide-Biosynthesis
PWY0-1337	oleate &beta;-oxidation	Degradation|Fatty-Acid-and-Lipid-Degradation|Fatty-Acid-Degradation
PWY-6386	UDP-N-acetylmuramoyl-pentapeptide biosynthesis II (lysine-containing)	Biosynthesis|Cell-Structure-Biosynthesis|Cell-Wall-Biosynthesis|Acetylmuramoyl-Pentapeptide-Biosynthesis
PWY-2724	alkane oxidation	Degradation|Fatty-Acid-and-Lipid-Degradation|Fatty-Acid-Degradation
PWY66-391	fatty acid &beta;-oxidation VI (peroxisome)	Degradation|Fatty-Acid-and-Lipid-Degradation|Fatty-Acid-Degradation
PWY-5138	unsaturated, even numbered fatty acid &beta;-oxidation	Degradation|Fatty-Acid-and-Lipid-Degradation|Fatty-Acid-Degradation
PWY-6168	flavin biosynthesis III (fungi)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Flavin-Biosynthesis
PWY-7291	oleate &beta;-oxidation (isomerase-dependent, yeast)	Degradation|Fatty-Acid-and-Lipid-Degradation|Fatty-Acid-Degradation
PWY-6167	flavin biosynthesis II (archaea)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Flavin-Biosynthesis
PWY-7776	ethene and chloroethene degradation	Degradation|Carbohydrates-Degradation|Alkene-Degradation
RIBOSYN2-PWY	flavin biosynthesis I (bacteria and plants)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Flavin-Biosynthesis
PWY-5534	propene degradation	Degradation|Carbohydrates-Degradation|Alkene-Degradation
PWY66-366	flavin biosynthesis IV (mammalian)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Flavin-Biosynthesis
PWY-7778	2-methylpropene degradation	Degradation|Carbohydrates-Degradation|Alkene-Degradation
PWY-7777	isoprene degradation	Degradation|Carbohydrates-Degradation|Alkene-Degradation
PWY-7775	propane degradation II	Degradation|Carbohydrates-Degradation|Propane-Degradation
PWY-7774	propane degradation I	Degradation|Carbohydrates-Degradation|Propane-Degradation
DHGLUCONATE-PYR-CAT-PWY	glucose degradation (oxidative)	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ACIDS-DEG
FUCCAT-PWY	fucose degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
GLUCOSE1PMETAB-PWY	glucose and glucose-1-phosphate degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
LYXMET-PWY	L-lyxose degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
MALTOSECAT-PWY	maltose degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
PWY-7761	NAD salvage pathway II	Biosynthesis|Cofactor-Biosynthesis|NAD-Metabolism|NAD-SYN
MANNCAT-PWY	D-mannose degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
PYRIDNUCSYN-PWY	NAD biosynthesis I (from aspartate)	Biosynthesis|Cofactor-Biosynthesis|NAD-Metabolism|NAD-SYN
P124-PWY	Bifidobacterium shunt	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
NAD-BIOSYNTHESIS-III	NAD biosynthesis III	Biosynthesis|Cofactor-Biosynthesis|NAD-Metabolism|NAD-SYN
P302-PWY	L-sorbose degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
PWY-5653	NAD biosynthesis from 2-amino-3-carboxymuconate semialdehyde	Biosynthesis|Cofactor-Biosynthesis|NAD-Metabolism|NAD-SYN
PWY-1081	homogalacturonan degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
PYRIDNUCSAL-PWY	NAD salvage pathway I	Biosynthesis|Cofactor-Biosynthesis|NAD-Metabolism|NAD-SYN
PWY-3861	mannitol degradation II	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ALCOHOLS-DEG
PWY3O-4106	NAD salvage pathway IV	Biosynthesis|Cofactor-Biosynthesis|NAD-Metabolism|NAD-SYN
PWY-5257	superpathway of pentose and pentitol degradation	Super-Pathways
PWY-6527	stachyose degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
PWY-6778	laminaribiose degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
PWY-6901	superpathway of glucose and xylose degradation	Super-Pathways
PWY-7077	N-acetyl-D-galactosamine degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
PWY-7130	L-glucose degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
PWY-7180	2-deoxy-&alpha;-D-ribose 1-phosphate degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
PWY-7395	D-galactosamine and N-acetyl-D-galactosamine degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
PWY-7459	kojibiose degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
PWY-7562	3,6-anhydro-&alpha;-L-galactopyranose degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
PWY-7581	N-acetylneuraminate and N-acetylmannosamine degradation II	Degradation|AMINE-DEG|NAN-MANNACs-degradation
PWY0-1300	2-<I>O</I>-&alpha;-mannosyl-D-glycerate degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
PWY0-1301	melibiose degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
PWY-7409	phospholipid remodeling (phosphatidylethanolamine, yeast)	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis|PhosphatidylethanolamineBiosynthesis
PWY0-1309	chitobiose degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
PWY-6273	phosphatidylethanolamine biosynthesis III	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis|PhosphatidylethanolamineBiosynthesis
PWY0-1314	fructose degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
PWY4FS-6	phosphatidylethanolamine biosynthesis II	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis|PhosphatidylethanolamineBiosynthesis
PWY0-1324	N-acetylneuraminate and N-acetylmannosamine degradation I	Degradation|AMINE-DEG|NAN-MANNACs-degradation
PWY-5669	phosphatidylethanolamine biosynthesis I	Biosynthesis|Lipid-Biosynthesis|Phospholipid-Biosynthesis|PhosphatidylethanolamineBiosynthesis
PWY0-44	D-allose degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
PWY-3561	choline biosynthesis III	Biosynthesis|Lipid-Biosynthesis|Choline-Biosynthesis
PWY4FS-12	VTC2 cycle	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
PWY-3542	choline biosynthesis II	Biosynthesis|Lipid-Biosynthesis|Choline-Biosynthesis
PWY4FS-13	extended VTC2 cycle	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
PWY-3385	choline biosynthesis I	Biosynthesis|Lipid-Biosynthesis|Choline-Biosynthesis
RIBOKIN-PWY	ribose degradation	Degradation|Carbohydrates-Degradation|Sugars-And-Polysaccharides-Degradation
PWY-6854	ethylene biosynthesis III (microbes)	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|Ethylene-Biosynthesis
PWY-5874	heme degradation	Degradation|COFACTOR-DEGRADATION|Heme-Degradation
SALVPURINE2-PWY	xanthine and xanthosine salvage	Biosynthesis|Nucleotide-Biosynthesis|PUR-NUC-SYN|Purine-Nucleotides-Salvage
PWY-5130	2-oxobutanoate degradation I	Super-Pathways
PWY-6605	adenine and adenosine salvage II	Biosynthesis|Nucleotide-Biosynthesis|PUR-NUC-SYN|Purine-Nucleotides-Salvage|Adenine-Adenosine-Salvage
2OXOBUTYRATECAT-PWY	2-oxobutanoate degradation II	Degradation|CARBOXYLATES-DEG|2-Oxobutanoate-Degradation
PWY-6611	adenine and adenosine salvage V	Biosynthesis|Nucleotide-Biosynthesis|PUR-NUC-SYN|Purine-Nucleotides-Salvage|Adenine-Adenosine-Salvage
P121-PWY	adenine and adenosine salvage I	Biosynthesis|Nucleotide-Biosynthesis|PUR-NUC-SYN|Purine-Nucleotides-Salvage|Adenine-Adenosine-Salvage
PWY-6853	ethylene biosynthesis II (microbes)	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|Ethylene-Biosynthesis
PWY-6610	adenine salvage	Biosynthesis|Nucleotide-Biosynthesis|PUR-NUC-SYN|Purine-Nucleotides-Salvage|Adenine-Adenosine-Salvage
ETHYL-PWY	ethylene biosynthesis I (plants)	Biosynthesis|HORMONE-SYN|Plant-Hormone-Biosynthesis|Ethylene-Biosynthesis
PWY-6609	adenine and adenosine salvage III	Biosynthesis|Nucleotide-Biosynthesis|PUR-NUC-SYN|Purine-Nucleotides-Salvage|Adenine-Adenosine-Salvage
PWY-6619	adenine and adenosine salvage VI	Biosynthesis|Nucleotide-Biosynthesis|PUR-NUC-SYN|Purine-Nucleotides-Salvage|Adenine-Adenosine-Salvage
PWY3O-355	stearate biosynthesis III (fungi)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Stearate-Biosynthesis
PWY-5989	stearate biosynthesis II (bacteria and plants)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Stearate-Biosynthesis
PWY-7255	ergothioneine biosynthesis I (bacteria)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Ergothioneine-Biosynthesis
GLYCOLYSIS-E-D	superpathway of glycolysis and Entner-Doudoroff	Super-Pathways
PWY-5972	stearate biosynthesis I (animals and fungi)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Stearate-Biosynthesis
GLYCOLYSIS-TCA-GLYOX-BYPASS	superpathway of glycolysis, pyruvate dehydrogenase, TCA, and glyoxylate bypass	Super-Pathways
PWY-7550	ergothioneine biosynthesis II (fungi)	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Ergothioneine-Biosynthesis
GLYOXYLATE-BYPASS	glyoxylate cycle	Energy-Metabolism
PWY-3041	monoterpene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|MONOTERPENOID-SYN
PWY-5197	lactate biosynthesis (archaea)	Energy-Metabolism
PWY-5277	thiosulfate disproportionation I (thiol-dependent)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Thiosulfate-Disproportionation
PWY-5306	superpathway of thiosulfate metabolism (Desulfovibrio sulfodismutans)	Super-Pathways
PWY-5352	thiosulfate disproportionation II (cytochrome)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Thiosulfate-Disproportionation
PWY-6958	icosapentaenoate biosynthesis I (lower eukaryotes)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis|Eicosapentaenoate-Biosynthesis
PWY-5464	superpathway of cytosolic glycolysis (plants), pyruvate dehydrogenase and TCA cycle	Super-Pathways
PWY-561	superpathway of glyoxylate cycle and fatty acid degradation	Super-Pathways
PWY-7602	icosapentaenoate biosynthesis V (8-desaturase, lower eukaryotes)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis|Eicosapentaenoate-Biosynthesis
PWY-5723	Rubisco shunt	Energy-Metabolism
PWY-5725	farnesene biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|SESQUITERPENOID-SYN
PWY-7049	icosapentaenoate biosynthesis II (6-desaturase, mammals)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis|Eicosapentaenoate-Biosynthesis
PWY-5741	ethylmalonyl-CoA pathway	Energy-Metabolism
PWY-6728	methylaspartate cycle	Energy-Metabolism
PWY-7724	icosapentaenoate biosynthesis III (8-desaturase, mammals)	Biosynthesis|Lipid-Biosynthesis|Fatty-acid-biosynthesis|Unsaturated-Fatty-Acids-Biosynthesis|PUFA-Biosynthesis|Eicosapentaenoate-Biosynthesis
PWY-6859	all-trans-farnesol biosynthesis	Energy-Metabolism
PWY-6871	3-methylbutanol biosynthesis (engineered)	Energy-Metabolism
PWY-6873	long chain fatty acid ester synthesis (engineered)	Energy-Metabolism
PWY-6876	isopropanol biosynthesis (engineered)	Energy-Metabolism
PWY-6938	NADH repair	Energy-Metabolism
PWY-7226	guanosine deoxyribonucleotides de novo biosynthesis I	Biosynthesis|Nucleotide-Biosynthesis|Deoxyribonucleotide-Biosynthesis|Purine-Deoxyribonuc-De-Novo-Biosynthesis|Guanosine-Deoxy-Denovo-Biosynthesis
PWY-6992	1,5-anhydrofructose degradation	Energy-Metabolism
PWY-7222	guanosine deoxyribonucleotides de novo biosynthesis II	Biosynthesis|Nucleotide-Biosynthesis|Deoxyribonucleotide-Biosynthesis|Purine-Deoxyribonuc-De-Novo-Biosynthesis|Guanosine-Deoxy-Denovo-Biosynthesis
PWY-7007	methyl ketone biosynthesis (engineered)	Energy-Metabolism
PWY-7102	bisabolene biosynthesis (engineered)	Energy-Metabolism
PWY-7118	chitin degradation to ethanol	Energy-Metabolism
PWY-7126	ethylene biosynthesis IV (engineered)	Energy-Metabolism
PWY-7178	ethylene glycol biosynthesis (engineered)	Energy-Metabolism
PWY-7709	(3R)-linalool biosynthesis	Biosynthesis|SECONDARY-METABOLITE-BIOSYNTHESIS|Terpenoid-Biosynthesis|MONOTERPENOID-SYN|Linalool-Biosynthesis
PWY-7813	thiosulfate disproportionation III (quinone)	Degradation|Noncarbon-Nutrients|Sulfur-Metabolism|Thiosulfate-Disproportionation
PWY0-1517	sedoheptulose bisphosphate bypass	Energy-Metabolism
PWY2OL-4	superpathway of linalool biosynthesis	Super-Pathways
RUMP-PWY	formaldehyde oxidation I	Degradation|C1-COMPOUNDS|Formaldehyde-Oxidation
ANAEROFRUCAT-PWY	homolactic fermentation	Super-Pathways
FERMENTATION-PWY	mixed acid fermentation	Energy-Metabolism|Fermentation
GLUDEG-II-PWY	L-glutamate degradation VII (to butanoate)	Super-Pathways
P122-PWY	heterolactic fermentation	Energy-Metabolism|Fermentation
P162-PWY	L-glutamate degradation V (via hydroxyglutarate)	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|GLUTAMATE-DEG
P163-PWY	L-lysine fermentation to acetate and butanoate	Degradation|Amino-Acid-Degradation|Proteinogenic-Amino-Acids-Degradation|LYSINE-DEG
P461-PWY	hexitol fermentation to lactate, formate, ethanol and acetate	Super-Pathways
PROPFERM-PWY	L-alanine fermentation to propanoate and acetate	Super-Pathways
PWY-5088	L-glutamate degradation VIII (to propanoate)	Super-Pathways
PWY-5109	2-methylbutanoate biosynthesis	Energy-Metabolism|Fermentation
PWY-5677	succinate fermentation to butanoate	Energy-Metabolism|Fermentation
PWY-7227	adenosine deoxyribonucleotides de novo biosynthesis	Biosynthesis|Nucleotide-Biosynthesis|Deoxyribonucleotide-Biosynthesis|Purine-Deoxyribonuc-De-Novo-Biosynthesis|Adenosine-Deoxy-Denovo-Bbiosynthesis
PWY-6130	glycerol degradation III	Degradation|Alcohol-Degradation|GLYCEROL-DEG
PWY-7220	adenosine deoxyribonucleotides de novo biosynthesis II	Biosynthesis|Nucleotide-Biosynthesis|Deoxyribonucleotide-Biosynthesis|Purine-Deoxyribonuc-De-Novo-Biosynthesis|Adenosine-Deoxy-Denovo-Bbiosynthesis
PWY-7383	anaerobic energy metabolism (invertebrates, cytosol)	Energy-Metabolism|Fermentation
PWY66-421	homocarnosine biosynthesis	Biosynthesis|Metabolic-Regulators|Dipeptide-Biosynthesis
PWY-7384	anaerobic energy metabolism (invertebrates, mitochondrial)	Super-Pathways
PWY-6459	peptidoglycan cross-bridge biosynthesis I (S. aureus)	Biosynthesis|Cell-Structure-Biosynthesis|Cell-Wall-Biosynthesis|Peptidoglycan-Cross-Bridge-Biosynthesis
PWY-7385	1,3-propanediol biosynthesis (engineered)	Energy-Metabolism|Fermentation
PWY-6463	peptidoglycan cross-bridge biosynthesis IV (Weissella viridescens)	Biosynthesis|Cell-Structure-Biosynthesis|Cell-Wall-Biosynthesis|Peptidoglycan-Cross-Bridge-Biosynthesis
PWY-7389	superpathway of anaerobic energy metabolism (invertebrates)	Super-Pathways
PWY-6462	peptidoglycan cross-bridge biosynthesis III (Enterococcus faecalis)	Biosynthesis|Cell-Structure-Biosynthesis|Cell-Wall-Biosynthesis|Peptidoglycan-Cross-Bridge-Biosynthesis
PWY-7401	crotonate fermentation (to acetate and cyclohexane carboxylate)	Super-Pathways
PWY-6461	peptidoglycan cross-bridge biosynthesis II (E. faecium)	Biosynthesis|Cell-Structure-Biosynthesis|Cell-Wall-Biosynthesis|Peptidoglycan-Cross-Bridge-Biosynthesis
PWY-7402	benzoate fermentation (to acetate and cyclohexane carboxylate)	Super-Pathways
PWY-7541	1,2-propanediol biosynthesis from lactate (engineered)	Metabolic-Clusters
PWY0-1312	acetate formation from acetyl-CoA I	Degradation|CARBOXYLATES-DEG|Acetate-Formation
PWY-5536	acetate formation from acetyl-CoA III (succinate)	Degradation|CARBOXYLATES-DEG|Acetate-Formation
PWY-5535	acetate formation from acetyl-CoA II	Degradation|CARBOXYLATES-DEG|Acetate-Formation
PWY-6961	L-ascorbate degradation II (bacterial, aerobic)	Degradation|CARBOXYLATES-DEG|Ascorbate-Degradation
PWY-6960	L-ascorbate degradation III	Degradation|CARBOXYLATES-DEG|Ascorbate-Degradation
PWY-6959	L-ascorbate degradation V	Degradation|CARBOXYLATES-DEG|Ascorbate-Degradation
PWY0-301	L-ascorbate degradation I (bacterial, anaerobic)	Degradation|CARBOXYLATES-DEG|Ascorbate-Degradation
PWY0-1241	ADP-L-glycero-&beta;-D-manno-heptose biosynthesis	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|ADP-Sugar-Biosynthesis
PWY-6704	L-ascorbate degradation IV	Degradation|CARBOXYLATES-DEG|Ascorbate-Degradation
GLYOXDEG-PWY	glycolate and glyoxylate degradation II	Degradation|CARBOXYLATES-DEG|Glycolate-Degradation
GLYCOLATEMET-PWY	glycolate and glyoxylate degradation I	Degradation|CARBOXYLATES-DEG|Glycolate-Degradation
PWY-6649	glycolate and glyoxylate degradation III	Degradation|CARBOXYLATES-DEG|Glycolate-Degradation
PWY-7686	L-malate degradation II	Degradation|CARBOXYLATES-DEG|L-Malate-Degradation
PWY-7685	L-malate degradation I	Degradation|CARBOXYLATES-DEG|L-Malate-Degradation
PWY-5794	malonate degradation I (biotin-independent)	Degradation|CARBOXYLATES-DEG|Malonate-Degradation
PWY-6060	malonate degradation II (biotin-dependent)	Degradation|CARBOXYLATES-DEG|Malonate-Degradation
PWY-6694	oxalate degradation I	Degradation|CARBOXYLATES-DEG|Oxalate-Degradation
PWY-6697	oxalate degradation IV	Degradation|CARBOXYLATES-DEG|Oxalate-Degradation
PWY-6696	oxalate degradation III	Degradation|CARBOXYLATES-DEG|Oxalate-Degradation
PWY-6695	oxalate degradation II	Degradation|CARBOXYLATES-DEG|Oxalate-Degradation
PWY-6698	oxalate degradation V	Degradation|CARBOXYLATES-DEG|Oxalate-Degradation
PROPIONMET-PWY	propanoyl CoA degradation I	Degradation|CARBOXYLATES-DEG|Propionate-Degradation
PWY-7574	propanoyl-CoA degradation II	Degradation|CARBOXYLATES-DEG|Propionate-Degradation
PWY0-43	conversion of succinate to propanoate	Degradation|CARBOXYLATES-DEG|SUCC-DEG
PWY-7216	(R)- and (S)-3-hydroxybutanoate biosynthesis (engineered)	Biosynthesis|Storage-Compounds-Biosynthesis
PWY-6995	5-hydroxymethylfurfural degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION
PWY-7052	cyanophycin metabolism	Biosynthesis|Storage-Compounds-Biosynthesis
PWY-7114	tea aroma glycosidic precursor bioactivation	Degradation|SECONDARY-METABOLITE-DEGRADATION
PWY1-3	polyhydroxybutanoate biosynthesis	Biosynthesis|Storage-Compounds-Biosynthesis
PWY0-1527	curcumin degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION
PWY-6657	polyhydroxydecanoate biosynthesis	Biosynthesis|Storage-Compounds-Biosynthesis
PWY-2345	genistein conjugates interconversion	Activation-Inactivation-Interconversion|Interconversion
PWY-6507	4-deoxy-L-threo-hex-4-enopyranuronate degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives
PWY-2861	biochanin A conjugates interconversion	Activation-Inactivation-Interconversion|Interconversion
PWY-6516	superpathway of microbial D-galacturonate and D-glucuronate degradation	Super-Pathways
PWY-7056	daphnin interconversion	Activation-Inactivation-Interconversion|Interconversion
PWY0-1261	anhydromuropeptides recycling	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives
PWY-2343	daidzein conjugates interconversion	Activation-Inactivation-Interconversion|Interconversion
PWY0-521	fructoselysine and psicoselysine degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives
PWY-2701	maackiain conjugates interconversion	Activation-Inactivation-Interconversion|Interconversion
PWY-5945	zeaxanthin, antheraxanthin and violaxanthin interconversion	Activation-Inactivation-Interconversion|Interconversion
PWY-6972	oleandomycin activation/inactivation	Activation-Inactivation-Interconversion|Interconversion
PWY-2561	medicarpin conjugates interconversion	Activation-Inactivation-Interconversion|Interconversion
PWY-2904	formononetin conjugates interconversion	Activation-Inactivation-Interconversion|Interconversion
GALACTCAT-PWY	D-galactonate degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ACIDS-DEG
PWY-5926	afrormosin conjugates interconversion	Activation-Inactivation-Interconversion|Interconversion
IDNCAT-PWY	L-idonate degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ACIDS-DEG
PWY-6303	methyl indole-3-acetate interconversion	Activation-Inactivation-Interconversion|Interconversion
KETOGLUCONMET-PWY	ketogluconate metabolism	Super-Pathways
PWY-7057	cichoriin interconversion	Activation-Inactivation-Interconversion|Interconversion
PWY-7242	D-fructuronate degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ACIDS-DEG
PWY-7516	L-lyxonate degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ACIDS-DEG
PWY-7566	L-gulonate degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ACIDS-DEG
PWY0-1306	L-galactonate degradation	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|SUGAR-ACIDS-DEG
PWY-5493	reductive monocarboxylic acid cycle	Degradation|C1-COMPOUNDS|CO2-Fixation
PWYQT-4429	CO2 fixation into oxaloacetate (anaplerotic)	Degradation|C1-COMPOUNDS|CO2-Fixation
P185-PWY	formaldehyde assimilation III (dihydroxyacetone cycle)	Degradation|C1-COMPOUNDS|Formaldehyde-Assimilation
PWY-1861	formaldehyde assimilation II (RuMP Cycle)	Degradation|C1-COMPOUNDS|Formaldehyde-Assimilation
PWY-1622	formaldehyde assimilation I (serine pathway)	Degradation|C1-COMPOUNDS|Formaldehyde-Assimilation
FORMASS-PWY	formaldehyde oxidation IV (thiol-independent)	Degradation|C1-COMPOUNDS|Formaldehyde-Oxidation
PWY1G-170	formaldehyde oxidation III (mycothiol-dependent)	Degradation|C1-COMPOUNDS|Formaldehyde-Oxidation
PWY-1801	formaldehyde oxidation II (glutathione-dependent)	Degradation|C1-COMPOUNDS|Formaldehyde-Oxidation
PWY-7238	sucrose biosynthesis II	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|Sucrose-Biosynthesis
PWY-1723	formaldehyde oxidation V (H4MPT pathway)	Degradation|C1-COMPOUNDS|Formaldehyde-Oxidation
PWY-7347	sucrose biosynthesis III	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|Sucrose-Biosynthesis
PWY-6523	nitrite-dependent anaerobic methane oxidation	Degradation|C1-COMPOUNDS|Methane-Oxidation
PWY-1641	methane oxidation to methanol I	Degradation|C1-COMPOUNDS|Methane-Oxidation
PWY-6742	methane oxidation to methanol II	Degradation|C1-COMPOUNDS|Methane-Oxidation
PWY-5506	methanol oxidation to formaldehyde IV	Degradation|C1-COMPOUNDS|Methanol-Oxidation
PWY-6966	methanol oxidation to formaldehyde I	Degradation|C1-COMPOUNDS|Methanol-Oxidation
PWY-5921	glutaminyl-tRNA<sup>gln</sup> biosynthesis via transamidation	Biosynthesis|Aminoacyl-tRNAs-Charging
PWY-1701	methanol and methylamine oxidation to formaldehyde	Super-Pathways
PWY-6510	methanol oxidation to formaldehyde II	Degradation|C1-COMPOUNDS|Methanol-Oxidation
PWY-6139	CMP-N-acetylneuraminate biosynthesis II (bacteria)	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|CMP-Sugar-Biosynthesis|CMP-N-Acetylneuraminate-Biosynthesis
PWY-6509	methanol oxidation to formaldehyde III	Degradation|C1-COMPOUNDS|Methanol-Oxidation
PWY-6138	CMP-N-acetylneuraminate biosynthesis I (eukaryotes)	Biosynthesis|Carbohydrates-Biosynthesis|CARBO-BIOSYNTHESIS|SUGAR-NUCLEOTIDES|CMP-Sugar-Biosynthesis|CMP-N-Acetylneuraminate-Biosynthesis
CITRULLINE-DEG-PWY	L-citrulline degradation	Degradation|Amino-Acid-Degradation|MISCELLANEOUS-DEG
ORN-AMINOPENTANOATE-CAT-PWY	L-ornithine degradation I (L-proline biosynthesis)	Biosynthesis|Amino-Acid-Biosynthesis|IND-AMINO-ACID-SYN|PROLINE-SYN
PWY-4021	&beta;-alanine betaine biosynthesis	Degradation|Amino-Acid-Degradation|MISCELLANEOUS-DEG
PWY-6334	L-dopa degradation	Degradation|Amino-Acid-Degradation|MISCELLANEOUS-DEG
PWY-6344	L-ornithine degradation II (Stickland reaction)	Degradation|Amino-Acid-Degradation|MISCELLANEOUS-DEG
PWY0-1317	L-lactaldehyde degradation (aerobic)	Degradation|Aldehyde-Degradation|Lactaldehyde-Degradation
PWY0-1315	L-lactaldehyde degradation (anaerobic)	Degradation|Aldehyde-Degradation|Lactaldehyde-Degradation
PWY-6388	(S,S)-butanediol degradation	Degradation|Alcohol-Degradation|2-3-Butanediol-Degradation
PWY3O-246	(R,R)-butanediol degradation	Degradation|Alcohol-Degradation|2-3-Butanediol-Degradation
PWY66-162	ethanol degradation IV	Degradation|Alcohol-Degradation|Ethanol-Degradation
PWY66-161	ethanol degradation III	Degradation|Alcohol-Degradation|Ethanol-Degradation
ETOH-ACETYLCOA-ANA-PWY	ethanol degradation I	Degradation|Alcohol-Degradation|Ethanol-Degradation
PWY66-21	ethanol degradation II	Degradation|Alcohol-Degradation|Ethanol-Degradation
GOLPDLCAT-PWY	superpathway of glycerol degradation to 1,3-propanediol	Super-Pathways
PWY-6131	glycerol degradation II	Degradation|Alcohol-Degradation|GLYCEROL-DEG
GLYCEROLMETAB-PWY	glycerol degradation V	Degradation|Alcohol-Degradation|GLYCEROL-DEG
PWY0-381	glycerol and glycerophosphodiester degradation	Super-Pathways
PWY-6952	glycerophosphodiester degradation	Degradation|Alcohol-Degradation|GLYCEROL-DEG
PWY-6117	spermine and spermidine degradation I	Degradation|AMINE-DEG|SPERMINE-SPERMIDINE-DEG
PWY-6441	spermine and spermidine degradation III	Degradation|AMINE-DEG|SPERMINE-SPERMIDINE-DEG
PWY-6440	spermine and spermidine degradation II	Degradation|AMINE-DEG|SPERMINE-SPERMIDINE-DEG
GLUDEG-I-PWY	GABA shunt	Super-Pathways
PWY-6535	4-aminobutanoate degradation I	Degradation|AMINE-DEG|4-Aminobutyraye-Degradation
4AMINOBUTMETAB-PWY	superpathway of 4-aminobutanoate degradation	Super-Pathways
PWY-6473	4-aminobutanoate degradation IV	Degradation|AMINE-DEG|4-Aminobutyraye-Degradation
PWY-6537	4-aminobutanoate degradation II	Degradation|AMINE-DEG|4-Aminobutyraye-Degradation
PWY-5022	4-aminobutanoate degradation V	Degradation|AMINE-DEG|4-Aminobutyraye-Degradation
PWY-6536	4-aminobutanoate degradation III	Degradation|AMINE-DEG|4-Aminobutyraye-Degradation
PWY-3721	choline degradation II	Degradation|AMINE-DEG|Choline-Degradation
P542-PWY	choline-O-sulfate degradation	Super-Pathways
PWY-7494	choline degradation IV	Super-Pathways
CHOLINE-BETAINE-ANA-PWY	choline degradation I	Degradation|AMINE-DEG|Choline-Degradation
PWY-7167	choline degradation III	Degradation|AMINE-DEG|Choline-Degradation
PWY-6967	methylamine degradation I	Degradation|AMINE-DEG|Methylamine-Degradation
PWY-6965	methylamine degradation II	Degradation|AMINE-DEG|Methylamine-Degradation
PUTDEG-PWY	putrescine degradation I	Degradation|AMINE-DEG|Putrescine-Degradation
PWY-3	putrescine degradation V	Degradation|AMINE-DEG|Putrescine-Degradation
PWY-2	putrescine degradation IV	Degradation|AMINE-DEG|Putrescine-Degradation
PWY-0	putrescine degradation III	Degradation|AMINE-DEG|Putrescine-Degradation
PWY0-1221	putrescine degradation II	Degradation|AMINE-DEG|Putrescine-Degradation
CARNMET-PWY	L-carnitine degradation I	Degradation|AMINE-DEG|CARN-DEG
PWY-7471	D-carnitine degradation I	Degradation|AMINE-DEG|CARN-DEG
PWY-3641	L-carnitine degradation III	Degradation|AMINE-DEG|CARN-DEG
PWY-3602	L-carnitine degradation II	Degradation|AMINE-DEG|CARN-DEG
PWY-7472	D-carnitine degradation II	Degradation|AMINE-DEG|CARN-DEG
PWY-6620	guanine and guanosine salvage	Biosynthesis|Nucleotide-Biosynthesis|PUR-NUC-SYN|Purine-Nucleotides-Salvage|Guanine-Guanosine-Salvage
PWY-3661	glycine betaine degradation I	Degradation|AMINE-DEG|Glycine-Betaine-Degradation
PWY-6618	guanine and guanosine salvage III	Biosynthesis|Nucleotide-Biosynthesis|PUR-NUC-SYN|Purine-Nucleotides-Salvage|Guanine-Guanosine-Salvage
PWY-3661-1	glycine betaine degradation II (mammalian)	Degradation|AMINE-DEG|Glycine-Betaine-Degradation
PWY-6599	guanine and guanosine salvage II	Biosynthesis|Nucleotide-Biosynthesis|PUR-NUC-SYN|Purine-Nucleotides-Salvage|Guanine-Guanosine-Salvage
PWY-5704	urea degradation II	Degradation|AMINE-DEG|Urea-Degradation
PWY-5703	urea degradation I	Degradation|AMINE-DEG|Urea-Degradation
PWY-5694	allantoin degradation to glyoxylate I	Super-Pathways
PWY-5705	allantoin degradation to glyoxylate III	Super-Pathways
PWY-5692	allantoin degradation to glyoxylate II	Super-Pathways
PWY-5698	allantoin degradation to ureidoglycolate II (ammonia producing)	Degradation|AMINE-DEG|Allantoin-degradation
URDEGR-PWY	superpathway of allantoin degradation in plants	Super-Pathways
ALLANTOINDEG-PWY	superpathway of allantoin degradation in yeast	Super-Pathways
PWY-5697	allantoin degradation to ureidoglycolate I (urea producing)	Degradation|AMINE-DEG|Allantoin-degradation
PWY0-41	allantoin degradation IV (anaerobic)	Super-Pathways
CRNFORCAT-PWY	creatinine degradation I	Degradation|AMINE-DEG|Creatinine-Degradation
PWY-4741	creatinine degradation III	Degradation|AMINE-DEG|Creatinine-Degradation
PWY-4722	creatinine degradation II	Degradation|AMINE-DEG|Creatinine-Degradation
PWY-6517	N-acetylglucosamine degradation II	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|N-Acetylglucosamine-Degradation
GLUAMCAT-PWY	N-acetylglucosamine degradation I	Degradation|SECONDARY-METABOLITE-DEGRADATION|Sugar-Derivatives|N-Acetylglucosamine-Degradation
PWY-3722	glycine betaine biosynthesis II (Gram-positive bacteria)	Biosynthesis|Polyamine-Biosynthesis|Betaine-Biosynthesis
PWY1F-353	glycine betaine biosynthesis III (plants)	Biosynthesis|Polyamine-Biosynthesis|Betaine-Biosynthesis
BETSYN-PWY	glycine betaine biosynthesis I (Gram-negative bacteria)	Biosynthesis|Polyamine-Biosynthesis|Betaine-Biosynthesis
PWY-6004	glycine betaine biosynthesis V (from glycine)	Biosynthesis|Polyamine-Biosynthesis|Betaine-Biosynthesis
P541-PWY	glycine betaine biosynthesis IV (from glycine)	Biosynthesis|Polyamine-Biosynthesis|Betaine-Biosynthesis
PWY-43	putrescine biosynthesis II	Biosynthesis|Polyamine-Biosynthesis|Putrescine-Biosynthesis
PWY-40	putrescine biosynthesis I	Biosynthesis|Polyamine-Biosynthesis|Putrescine-Biosynthesis
PWY-6305	putrescine biosynthesis IV	Biosynthesis|Polyamine-Biosynthesis|Putrescine-Biosynthesis
PWY-46	putrescine biosynthesis III	Biosynthesis|Polyamine-Biosynthesis|Putrescine-Biosynthesis
UDPNACETYLGALSYN-PWY	UDP-N-acetyl-D-glucosamine biosynthesis II	Biosynthesis|Polyamine-Biosynthesis|UDP-NAc-Glucosamine-Biosynthesis
PWY-6892	thiazole biosynthesis I (facultative anaerobic bacteria)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Thiamine-Biosynthesis|Thiazole-Biosynthesis
PWY-6891	thiazole biosynthesis II (aerobic bacteria)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Thiamine-Biosynthesis|Thiazole-Biosynthesis
PWY-6909	thiazole biosynthesis III (eukaryotes)	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|Thiamine-Biosynthesis|Thiazole-Biosynthesis
PWY-7790	UMP biosynthesis II	Biosynthesis|Nucleotide-Biosynthesis|PYR-NUC-SYN|Pyrimidine-De-Novo-Biosyn|Pyrimid-Ribonucleot-De-Novo-Biosyn|UMP-Biosynthesis
PWY-5686	UMP biosynthesis I	Biosynthesis|Nucleotide-Biosynthesis|PYR-NUC-SYN|Pyrimidine-De-Novo-Biosyn|Pyrimid-Ribonucleot-De-Novo-Biosyn|UMP-Biosynthesis
PWY-7791	UMP biosynthesis III	Biosynthesis|Nucleotide-Biosynthesis|PYR-NUC-SYN|Pyrimidine-De-Novo-Biosyn|Pyrimid-Ribonucleot-De-Novo-Biosyn|UMP-Biosynthesis
PWY-7147	8-amino-7-oxononanoate biosynthesis II	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|BIOTIN-SYN|7-Keto-8-aminopelargonate-Biosynthesis
PWY-6578	8-amino-7-oxononanoate biosynthesis III	Biosynthesis|Cofactor-Biosynthesis|Vitamin-Biosynthesis|BIOTIN-SYN|7-Keto-8-aminopelargonate-Biosynthesis' > pwy.hierarchy
}

module load jdk
module load bedtools

#gff version 3
mygff=GFFFILE
mysam=SAMFILE

if [ "$mygff" == "" ];then
	echo "* No gff provided"
	exit
fi

if [ "$mysam" == "" ];then
        echo "* No sam provided"
        exit
fi

sname=$(echo $mysam |rev |cut -d "/" -f 1 |rev |sed "s/.sam//g")
gname=$(echo $mygff |rev |cut -d "/" -f 1 |rev |sed "s/.gff//g")

samtools view -b $mysam > $sname.bam
samtools sort $sname.bam -o $sname.sort.bam

java -Xms2g -Xmx32g -jar /c1/apps/picard-tools/1.96/MarkDuplicates.jar INPUT=$sname.sort.bam OUTPUT=$sname.markdup.bam METRICS_FILE=$sname.markdup.metrics AS=TRUE VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 REMOVE_DUPLICATES=TRUE

prokkagff2bedFunction $mygff > $gname.bed


#Note: When using bedtools 2.24.0 or later, the A and the B files are switched as follows:
bedtools coverage -hist -a $gname.bed -b $sname.markdup.bam > $sname.hist

getcovFunction
python get_coverage_for_genes.py -i $sname.hist > $sname.coverage
rm -f get_coverage_for_genes.py

checkEC=$(grep "eC_number=" $mygff | cut -f9 | cut -f1,2 -d ';' |awk '{if($0~";"){print 1;exit}}')
#this script is only wuth metacyc pathways

echo "checking gff results"
if [ "$checkEC" == "" ];then
	grep "eC_number=" $mygff | cut -f9 | cut -f1,2 -d ';'| sed 's/ID=//g'| sed 's/;eC_number=/\t/g' > $gname.ec
else
	#check is 1, this means that gff have other format than expected
	grep "eC_number=" $mygff | cut -f9 | cut -f1,3 -d ';'| sed 's/ID=//g'| sed 's/;eC_number=/\t/g' > $gname.ec
fi

echo "gene ID to gene name step"
grep -v "#" $mygff | grep -v ">" |awk -F"\t" '{if($3=="CDS")print}'|cut -f9 |awk -F";" '{for(i=1;i<=NF;i++){if($i~"ID="){printf "%s,",$i};if($i~"gene="){printf "%s",$i}}printf "\n"}' |awk -F"," '{if($2!=""){gsub("="," ");print $1, $2}}' |cut -d " " -f2,4 > id2gene.txt

awk '{if(NR==FNR){n[$1]=$2}else{if($1 in n){print $2"\t"n[$1]}}}' $sname.coverage id2gene.txt |awk -F"_|\t" '{if($3==""){print $1, $2}else{print $1, $3}}' |awk '{n[$1]+=$2}END{for(key in n){print key, n[key]}}' > newcoverage_$sname.tsv
rm id2gene.txt

echo "runing minpath"
minpathFunction
ectopwyFunction

python MinPath1.2.py -any $gname.ec -map ec.to.pwy -report $gname.metacyc.minpath > minpath.$gname.log
rm -f MinPath1.2.py 

echo "converting files for krona"
genes2kronaFunction

python genes.to.kronaTable.py -i $gname.ec -m ec.to.pwy -H pwy.hierarchy -n $gname -l <(grep "minpath 1" $gname.metacyc.minpath) -c $sname.coverage -o $gname.krona.metacyc.minpath.tab
rm -f genes.to.kronaTable.py ec.to.pwy pwy.hierarchy

echo "executing krona"
KRONAHOME/bin/ktImportText -o $gname.krona.metacyc.minpath.html $gname.krona.metacyc.minpath.tab

chainReaction "KRONING" $1 $DEBUG
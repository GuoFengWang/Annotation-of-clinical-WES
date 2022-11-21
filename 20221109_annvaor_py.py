# /usr/bin/env python
# coding=utf-8
###################################
#  Author : wang
#  Version : 1.0(2022/4/10)
#  annovar unknown fileter
#  E-mail : 1316014512@qq.com
###################################
import json
import os,sys
import re
import _thread
import threading
class Annovar2vcf():
	def __int__(self):
		global annotate_variation
		global annovar_database
		global convert2annovar
		global Gene_symbol_re
	def AnnovarPath(self,config):
		fileopen = open(config).read()
		return json.loads(fileopen)
	def AnnovarAll(self,vcffile,outdir,sample,annotate_variation,annovar_database,convert2annovar):
		c2a = convert2annovar+' -format vcf4 '+vcffile+'  > '+outdir+'/'+sample+'.avinput'
		c2g = convert2annovar+' -format vcf4 '+vcffile+'   -include -withzyg > '+outdir+'/'+sample+'.gffout'
		thousand_all = annotate_variation + ' -filter -dbtype 1000g2015aug_all -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput ' + annovar_database
		thousand_esa = annotate_variation + ' -filter -dbtype 1000g2015aug_eas -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput ' + annovar_database
		dbsnp        = annotate_variation+' -filter -dbtype avsnp150  -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput '+annovar_database
		ljb26_all    = annotate_variation+' -filter -dbtype ljb26_all  -otherinfo -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput '+annovar_database
		ref          = annotate_variation+' -buildver hg19  -hgvs  -splicing_threshold 15 -otherinfo -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput '+annovar_database
		#esp_all      = annotate_variation+' -filter -dbtype esp6500siv2_all   -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput '+annovar_database
		esp_all      = annotate_variation+' -filter -dbtype esp_6500_20220425   -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput '+annovar_database
		pLI          = annotate_variation + ' -regionanno -dbtype pLI  -colsWanted 4,5 -buildver hg19 -out ' + outdir + '/' + sample + ' ' + outdir + '/' + sample + '.avinput ' + annovar_database
		Decipher     = annotate_variation+' -regionanno -dbtype Decipher -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput '+annovar_database
		exac         = annotate_variation + ' -filter  -dbtype exac03 -otherinfo -buildver hg19 -out ' + outdir + '/' + sample + ' ' + outdir + '/' + sample + '.avinput ' + annovar_database
		second       = annotate_variation + ' -regionanno -dbtype secondary -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput '+annovar_database
		omim         = annotate_variation + ' -regionanno -dbtype omim_20220926 -buildver hg19 -colsWanted 4,5,6,7,8 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput '+annovar_database
		clinvar      = annotate_variation + ' -filter  -dbtype clinvar_20220908 -otherinfo -buildver hg19 -out ' + outdir + '/' + sample + ' ' + outdir + '/' + sample + '.avinput ' + annovar_database
		revel        = annotate_variation + ' -filter -dbtype revel -otherinfo -buildver hg19 -out ' + outdir + '/' + sample + ' ' + outdir + '/' + sample + '.avinput ' + annovar_database
		gene_intro   = annotate_variation + ' -regionanno -dbtype gene_intro -buildver hg19 -colsWanted 4,5,6,7 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput '+annovar_database
		#hpo          = annotate_variation + ' -regionanno -dbtype HPO4752 -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput '+annovar_database
		hpo          = annotate_variation + ' -regionanno -dbtype HPO4752_CN -buildver hg19 -colsWanted 4,5 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput '+annovar_database
		#intervar     = annotate_variation + ' -filter  -dbtype intervar_20180118 -otherinfo -buildver hg19 -out ' + outdir + '/' + sample + ' ' + outdir + '/' + sample + '.avinput ' + annovar_database
		#cosmic       = annotate_variation + ' -filter  -dbtype cosmic70 -otherinfo -buildver hg19 -out ' + outdir + '/' + sample + ' ' + outdir + '/' + sample + '.avinput ' + annovar_database
		#gnomad       = annotate_variation + ' -filter  -dbtype gnomad_genome -otherinfo -buildver hg19 -out ' + outdir + '/' + sample + ' ' + outdir + '/' + sample + '.avinput ' + annovar_database
		gnomad       = annotate_variation + ' -filter  -dbtype gnomad211_genome -otherinfo -buildver hg19 -out ' + outdir + '/' + sample + ' ' + outdir + '/' + sample + '.avinput ' + annovar_database
		#gnomad_exome = annotate_variation + ' -filter  -dbtype gnomad_exome -otherinfo -buildver hg19 -out ' + outdir + '/' + sample + ' ' + outdir + '/' + sample + '.avinput ' + annovar_database
		gnomad_exome = annotate_variation + ' -filter  -dbtype gnomad211_exome -otherinfo -buildver hg19 -out ' + outdir + '/' + sample + ' ' + outdir + '/' + sample + '.avinput ' + annovar_database
		c2ainfo = os.system(c2a)
		alllist = []
		def c2gtd():
			c2ginfo = os.system(c2g)
			alllist.append("c2g")
		def t_alltd():
			t_allinfo = os.system(thousand_all)
			alllist.append("thousand_all")
		def pLItd():
			pLIinfo = os.system(pLI)
			alllist.append("pLI")
		def Deciphertd():
			Decipherinfo = os.system(Decipher)
			alllist.append("Decipher")
		def gene_introtd():
			gene_introinfo = os.system(gene_intro)
			alllist.append("gene_intro")
		def dstd():
			dsinfo = os.system(dbsnp)
			alllist.append("dbsnp")
		def ljb26_alltd():
			ljb26_allinfo = os.system(ljb26_all)
			alllist.append("ljb26_all")
		def reftd():
			refinfo = os.system(ref)
			alllist.append("ref")
		def esptd():
			espinfo = os.system(esp_all)
			alllist.append('esp6500')
		def espeatd():
			espinfo = os.system(esp_ea)
			alllist.append("esp6500")
		def cybtd():
			cybinfo = os.system(cyb)
			alllist.append('cytoBand')
		def t_esatd():
			t_esainfo = os.system(thousand_esa)
			alllist.append('thousand_esa')
		def exactd():
			exacinfo = os.system(exac)
			alllist.append('exac')
		def intervartd():
			intervarinfo = os.system(intervar)
			alllist.append('intervar')
		def secondtd():
			secondinfo = os.system(second)
			alllist.append('second')
		def omimtd():
			omiminfo = os.system(omim)
			alllist.append('omim')
		def clinvartd():
			clinvarinfo = os.system(clinvar)
			alllist.append('clinvar')
		def hpotd():
			hpoinfo = os.system(hpo)
			alllist.append('hpo')
		def cosmictd():
			cosmicinfo = os.system(cosmic)
			alllist.append('cosmic')
		def gnomadtd():
			gnomadinfo = os.system(gnomad)
			alllist.append('gnomad')
		def gnomad_exometd():
			gnomadinfo = os.system(gnomad_exome)
			alllist.append('gnomad_exome')
		def reveltd():
			revelinfo = os.system(revel)
			alllist.append('revel')
		if c2ainfo == 0:
			t1 = threading.Thread(target=c2gtd)
			t2 = threading.Thread(target=t_alltd)
			t3 = threading.Thread(target=dstd)
			t4 = threading.Thread(target=ljb26_alltd)
			t5 = threading.Thread(target=reftd)
			t6 = threading.Thread(target=esptd)
			#t7 = threading.Thread(target=espeatd)
			#t8 = threading.Thread(target=cybtd)
			t7 = threading.Thread(target=t_esatd)
			t8 = threading.Thread(target=exactd)
			#t11 = threading.Thread(target=intervartd)
			t9 = threading.Thread(target=clinvartd)
			#t13 = threading.Thread(target=cosmictd)
			t10 = threading.Thread(target=gnomadtd)
			t11 = threading.Thread(target=gnomad_exometd)
			t12 = threading.Thread(target=secondtd)
			t13 = threading.Thread(target=omimtd)
			t14 = threading.Thread(target=hpotd)
			t15 = threading.Thread(target=reveltd)
			t16 = threading.Thread(target=gene_introtd)
			t17 = threading.Thread(target=pLItd)
			t18 = threading.Thread(target=Deciphertd)
			t1.start()
			t2.start()
			t3.start()
			t4.start()
			t5.start()
			t6.start()
			t7.start()
			t8.start()
			t9.start()
			t10.start()
			t11.start()
			t12.start()
			t13.start()
			t14.start()
			t15.start()
			t16.start()
			t17.start()
			t18.start()
			while True:
				if len(alllist) == 18:
					#print alllist
					return 0
		else:
			print('annovar is erro please check log')
			sys.exit()
	def Gene_symbol_re(self,content):
		link = re.compile("\(.*?\)")
		info = re.sub(link,"",content)
		return info
	def dictavsnp150(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_avsnp150_dropped').readlines()
		snp150dict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			snp150dict[key] = values
		# fileopen.close()
		return snp150dict
	def dict1000gall(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_ALL.sites.2015_08_dropped').readlines()
		all1000gdict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			all1000gdict[key] = values
		# fileopen.close()
		return all1000gdict
	def dict1000geas(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_EAS.sites.2015_08_dropped').readlines()
		eas1000gdict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			eas1000gdict[key] = values
		# fileopen.close()
		return eas1000gdict
	def dictexac(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_exac03_dropped').readlines()
		exacdict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = [linelist[1].split(',')[0],linelist[1].split(',')[3]]
			exacdict[key] = values
		# fileopen.close()
		return exacdict
	def dictesp(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_esp_6500_20220425_dropped').readlines()
		espdict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			#values = linelist[1]
			values = linelist[1].split('=')[1]
			espdict[key] = values
		# fileopen.close()
		return espdict
	def dictespea(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_esp6500siv2_ea_dropped').readlines()
		espeadict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			espeadict[key] = values
		# fileopen.close()
		return espeadict
	def dictcyb(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_cytoBand').readlines()
		cybdict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			cybdict[key] = values
		return cybdict
	def dictsecond(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_secondary').readlines()
		seconddict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			seconddict[key] = values
		return seconddict
	def dictomim(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_omim_20220926').readlines()
		omimdict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			omimdict[key] = values 
		return omimdict
	def dictpLI(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_pLI').readlines()
		pLIdict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			pLIdict[key] = values 
		return pLIdict
	def dictDecipher(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_Decipher').readlines()
		Decipherdict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			Decipherdict[key] = values 
		return Decipherdict
	def dictgene_intro(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_gene_intro').readlines()
		gene_introdict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			gene_introdict[key] = values
		return gene_introdict
	def dicthpo(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_HPO4752_CN').readlines()
		hpodict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1].replace('Name=','')
			hpodict[key] = values
		return hpodict
			
	def dictgnomad(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_gnomad211_genome_dropped').readlines()
		gnomaddict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			gnomaddict[key] = values
		return gnomaddict
	def dictgnomad_exome(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_gnomad211_exome_dropped').readlines()
		gnomad_exomedict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			gnomad_exomedict[key] = values
		return gnomad_exomedict
	def dictclinvar(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_clinvar_20220422_dropped').readlines()
		clinvardict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values1 = linelist[1].split(',')
			#values1.reverse()
			#values11 = ','.join(values1)
			#values = str(values11).replace('_',' ')
			clinvardict[key] = values1
		return clinvardict
	def dictintervar(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_intervar_20180118_dropped').readlines()
		intervardict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1].split(',')[0]
			intervardict[key] = values
		return intervardict
	def dictcosmic(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_cosmic70_dropped').readlines()
		cosmicdict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			cosmicdict[key] = values
		return cosmicdict
	def dictvariant(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.variant_function').readlines()
		variantdict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = [linelist[0],linelist[1]]
			variantdict[key] = values
		# fileopen.close()
		return variantdict
	def dictexonic(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.exonic_variant_function').readlines()
		exonicdict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[3]+'|'+linelist[4]+'|'+linelist[5]+'|'+linelist[7]
			values = [linelist[1].split('\t')[0],linelist[2]]
			exonicdict[key] = values
		# fileopen.close()
		return exonicdict
	def dictrevel(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_revel_dropped').readlines()
		reveldict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1].split(',')[2]
			reveldict[key] = values
		return reveldict
	def dictljb(self,outdir,sample):
		def substitutions(flt):
			try:
				if float(flt) >= 2:
					return "rejected_substitutions"
				else:
					return "substitutions"
			except:
				return '.'
		ljbsift = {'D':'Deleterious','T':'tolerated','.':'.'}
		ljbP_HVAR = {"D":"Probably damaging","P":"possibly damaging","B":"benign",".":"."}
		MutationTaster = {'A':'disease_causing_automatic','D':'disease_causing','N':'polymorphism','P':'polymorphism_automatic','.':'.'}
		fileopen = open(outdir+'/'+sample+'.hg19_ljb26_all_dropped').readlines()
		ljbdict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			v = linelist[1].split(',')
			#ljbdict[key] = [v[0],ljbsift[v[1]],v[4],ljbP_HVAR[v[5]],v[8],MutationTaster[v[9]],v[21],substitutions(v[21])]
			ljbdict[key] = [v[0],v[1],v[8],v[9],v[10],v[11],v[4],v[5]]
		return ljbdict
	def MergeAll(self,outdir,sample,exonicdict,snp150dict,all1000gdict,eas1000gdict,espdict,variantdict,ljbdict,gnomaddict,gnomad_exomedict,seconddict,omimdict,pLIdict,Decipherdict,gene_introdict,clinvardict,hpodict,reveldict,exacdict):
		def Gene_symbol_re(content):
			link = re.compile("\(.*?\)")
			info = re.sub(link,"",content)
			return info
		fileopen9 = open(outdir+'/'+sample+'.gffout').readlines()
		head = '#Chr\tStart\tStop\tRef\tCall\tGene_name\tNMtranscript\tc.BaseChange\tp.AA_Change\tfull_info\tRs\texon_Number\tFunc.refGene\tExonicFunc.refGene\theterozygosity\tref/alt\tSecondary_finding_genes\tPhenotype MIM number\tomim_disease\tomim_ch\tinheritance_pattern\toe_lof_upper\tpLI\t%HI\tintro_gene\tintro_en\tintro_ch\tClinvar_interpretation\tClinvar_status\tExplanation\tInterpretedCondition\tHPO_terms\tHPO_CN\tgnomad_ALL_AF\tgnomad_ALL_East_Asian_AF\tgnomad_exome_ALL\tgnomad_exome_East_Asian_AF\tExAC_AF\tExAC_East_Asian_AF\t1000G_AF\t1000G_East_AF\tESP6500_AF\tREVEL\tSIFT_score\tSIFT_pred\tMutationTaster_score\tMutationTaster_pred\tMutationAssessor_score\tMutationAssessor_pred\tPolyphen2_HVAR_score\tPolyphen2_HVAR_Pred\n'
		filesave = open(outdir+'/'+sample+'.gff','w+')
		filesave_fileter = open(outdir+'/'+'filter_'+sample+'.gff','w+')
		#filesave_stat = open(outdir+'/'+sample+'.stat','w+')
		filesave.write(head)
		filesave_fileter.write(head)
		#print head
		count = 0
		unknown_count = 0
		for x in fileopen9:
			p =  x.split('\t')
			key = p[0]+'|'+p[1]+'|'+p[2]+'|'+p[4]
			pend = p[0:5]
			if key in variantdict: #Gene基因
				pend.append(Gene_symbol_re(variantdict[key][1]))
			else:
				pend.append('-')
			'''
			if key in exonicdict: # NM转录本，c.碱基改变，p.氨基酸改变
				line = exonicdict[key][1] # NM转录本，c.碱基改变，p.氨基酸改变信息所在行
				nm_r = re.compile(r"NM_\d+")
				nm_list = nm_r.findall(line)
				nm = ",".join(nm_list)
				c_r = re.compile(r"c.\d+[A-Z]>[A-Z]")
				c_list = c_r.findall(line)
				c = ",".join(c_list)
				p_r = re.compile(r"p.[A-Z]\d+[A-Z]")
				p_list = p_r.findall(line)
				p_info = ",".join(p_list)
			'''
			if key in variantdict:
				if variantdict[key][0] == 'exonic':
					line = exonicdict[key][1] # NM转录本，c.碱基改变，p.氨基酸改变信息所在行
					nm_r = re.compile(r"NM_\d+")
					nm_list = nm_r.findall(line)
					nm = ",".join(nm_list)
					c_r = re.compile(r"c.\d+[A-Z]>[A-Z]")
					c_list = c_r.findall(line)
					c = ",".join(c_list)
					p_r = re.compile(r"p.[A-Z]\d+[A-Z]")
					p_list = p_r.findall(line)
					p_info = ",".join(p_list)
					pend.append(nm)
					pend.append(c)
					pend.append(p_info)
					pend.append(line)
				else:
					line = variantdict[key][1]
					nm_r = re.compile(r"NM_\d+")
					nm_list = nm_r.findall(line)
					if nm_list:
						nm = ",".join(nm_list)
					else:
						nm = '-'
					c_r_1 = re.compile(r"c.[-+*]?\d+[+-]?\d+[A-Z]>[A-Z]")
					c_r_2 = re.compile(r"c.[-+*]?\d+[a-z]+[A-Z]+")
					c_list1 = c_r_1.findall(line)
					c_list2 = c_r_2.findall(line)
					c_list = c_list1 + c_list2
					if c_list:
						c = ",".join(c_list)
					else:
						c = "-"
					p_info = "-"	
					pend.append(nm)
					pend.append(c)
					pend.append(p_info)
					pend.append(line)
			if key in snp150dict:  #rs
				pend.append(snp150dict[key])
			else:
				pend.append('-')
			if key in exonicdict: # exon Number
				line = exonicdict[key][1] # NM转录本，c.碱基改变，p.氨基酸改变信息所在行
				exon_re = re.compile(r"exon\d+")
				exon_list = exon_re.findall(line)
				exon_number = str(len(set(exon_list)))
				pend.append(exon_number)
			else:
				pend.append('-')
			if key in variantdict: #Func.refGene
				pend.append(variantdict[key][0])
			else:
				pend.append('-')
			if key in exonicdict: # ExonicFunc.refGene
				pend.append(exonicdict[key][0])
			else:
				pend.append('-')
			pend.append(p[5]) # heterozygosity
			qc = p[-1].split(':')
			if len(qc) == 5:
				pend.append(qc[1]) # ref/alt
			else:
				pend.append('-')
			if key in seconddict: # Secondary finding
				pend.append(seconddict[key])
			else:
				pend.append('-')
			if key in omimdict:
				#pend = pend+omimdict[key] # OMIM数据库：Phenotype MIM number,omim_disease,inheritance pattern
				#pend.append(omimdict[key][0])
				#pend.append(omimdict[key][1])
				#pend.append(omimdict[key][2])
				line2 = omimdict[key]
				omim_num = []
				omim_disease = []
				omim_ch = []
				omim_inher = []
				for i in range(len(line2.split(','))):
					record = line2.split(',')[i].split(':')
					omim_num.append(str(record[0].replace('Name=','')))
					omim_disease.append(record[1])
					omim_ch.append(record[2])
					omim_inher.append(record[-1])
				pend.append(':'.join(omim_num))
				pend.append(','.join(omim_disease))
				pend.append(','.join(omim_ch))
				pend.append(','.join(omim_inher))
			else:
				pend = pend+['-', '-', '-','-']
			if key in pLIdict:       # lof  pLI
				line2 =  pLIdict[key]
				lof_list = []
				pLI_list = []
				for i in range(len(line2.split(','))):
					record = line2.split(',')[i].split(':')
					lof_list.append(str(record[0].replace('Name=','')))
					pLI_list.append(record[1])
				pend.append(','.join(lof_list))
				pend.append(','.join(pLI_list))
			else:
				pend = pend+['-','-']
			if key in Decipherdict:
				line2 = Decipherdict[key]
				hi_list = []
				for i in range(len(line2.split(','))):
					record = line2.split(',')[i].split(':')
					hi_list.append(str(record[0].replace('Name=','')))
				pend.append(','.join(hi_list))
			else:
				pend = pend+['-']
				
			if key in gene_introdict: # gene_intro
				line2 = gene_introdict[key]
				intro_num = []
				intro_gene = []
				intro_en = []
				intro_ch = []
				
				for i in range(len(line2.split(','))):
					record = line2.split(',')[i].split(':')
					intro_num.append(str(record[0].replace('Name=','')))
					try:
						intro_gene.append(record[1])
					except IndexError as e:
						print(line2,e)
						intro_gene.append('-')
					try:
						intro_en.append(record[2])
					except IndexError as e:
						print(line2,e)
						intro_en.append('-')
					try:
						intro_ch.append(record[3])
					except IndexError as e:
						print(line2,e)
						intro_ch.append('-')
					
				#pend.append(':'.join(intro_num))
				pend.append(':'.join(intro_gene))
				pend.append(':'.join(intro_en))
				pend.append(':'.join(intro_ch))
				
			else:
				pend = pend+['-']*3

			if key in clinvardict:  # Clinvar数据库：Clinvar_interpretation，Clinvar_status，Explanation，InterpretedCondition
				#pend = pend+clinvardict[key]
				pend.append(clinvardict[key][0])
				pend.append(clinvardict[key][1])
				pend.append(clinvardict[key][2])
				try:
					pend.append(clinvardict[key][3])
				except IndexError:
					print(key, clinvardict[key])
					pend.append('not provided')

			else:
				pend = pend+['-', '-', '-', '-']
			if key in hpodict: #HPO
				#pend.append(hpodict[key])
				pend.append(hpodict[key].split(':')[0])
				pend.append(hpodict[key].split(':')[1].strip())
			else:
				pend.append('-')
				pend.append('-')
			# 人群频率信息
			if key in gnomaddict:
				pend.append(gnomaddict[key].split(',')[0]) #gnomad_ALL_AF
				pend.append(gnomaddict[key].split(',')[8]) #gnomad_ALL_East Asian_AF
			else:
				pend.append('-')
				pend.append('-')
			if key in gnomad_exomedict:
				pend.append(gnomad_exomedict[key].split(',')[0]) #gnomad_exome_AF
				pend.append(gnomad_exomedict[key].split(',')[8])  # gnomad_exome_Asian_AF
			else:
				pend.append('-')
				pend.append('-')
			if key in exacdict:
				pend.append(exacdict[key][0]) #ExAC_AF
				pend.append(exacdict[key][1]) #ExAC_East Asian_AF
			else:
				pend.append('-')
				pend.append('-')
			if key in all1000gdict: # 1000G_AF
				pend.append(all1000gdict[key])
			else:
				pend.append('-')
			if key in eas1000gdict: #1000G_East_AF
				pend.append(eas1000gdict[key])
			else:
				pend.append('-')
			if key in espdict: # ESP6500_AF
				pend.append(espdict[key])
			else:
				pend.append('-')
			if key in reveldict: # REVEL
				pend.append(reveldict[key])
			else:
				pend.append('-')
			if key in ljbdict:
				pend = pend + ljbdict[key] #SIFT_score,SIFT_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_score,MutationAssessor_pred,Polyphen2 HVAR  score,Polyphen2 HVAR Pred
			else:
				pend = pend + ['-','-','-','-','-','-','-','-']
			if "exonic" in pend[19] or  "splicing" in pend[19] or "UTR" in  pend[19]:
				try:
					filesave_fileter.write( "\t".join(pend)+"\n")
				except TypeError:
					print(pend)
			filesave.write("\t".join(pend)+"\n")
			'''
			#filter
			pend.append(p[14])
			#fs
			fsReg=re.compile(".*FS=(\d*.\d*).*")
			fsRes=re.match(fsReg,p[15])
			if fsRes:
				fs=fsRes.group(1)
			else:
				fs='-'
			pend.append(fs)
			#sor
			fsReg=re.compile(".*SOR=(\d*.\d*).*")
			sorRes=re.match(fsReg,p[15])
			if sorRes:
				sor=sorRes.group(1)
			else:
				sor='-'
			pend.append(sor)
			#qd
			fsReg=re.compile(".*QD=(\d*.\d*).*")
			qdRes=re.match(fsReg,p[15])
			if qdRes:
				qd=qdRes.group(1)
			else:
				qd='-'
			pend.append(qd)
			
			filesave.write("\t".join(pend)+"\n")
			'''
		filesave.close()
		filesave_fileter.close()
if  __name__ == "__main__":
	import argparse
	usage='annovar_py.py [--outdir] [--inputvcf] [--samplename] [--config]'
	if len(sys.argv) < 2:
		print(usage)
		print('please cmd "python annovar_py.py -h "')
		sys.exit()
	else:
		p = argparse.ArgumentParser(usage='annovar_py.py [--outdir] [--inputvcf] [--samplename] [--config]', description="annovar vcf to gff")
		p.add_argument('-i','--inputvcf',   type=str, help = 'input   dir eg: ./sample.vcf',required=True)
		p.add_argument('-s','--samplename', type=str, help = 'sample name',required=True)
		p.add_argument('-o','--outdir',     type=str, help ='oss document directory eg: ./yunkeli',required=True)
		p.add_argument('-c','--config', type=str, help ='config file',required=True)
		args = p.parse_args()
		inputvcf = args.inputvcf
		samplename = args.samplename
		outdir = args.outdir
		config = args.config
		project = Annovar2vcf()
		path = project.AnnovarPath(config)
		#print path
		annotate_variation = path["annotate_variation"]
		annovar_database = path["annovar_database"]
		convert2annovar = path["convert2annovar"]
		otherdatabase = path["otherdatabase"]
		stat = project.AnnovarAll(inputvcf,outdir,samplename,annotate_variation,annovar_database,convert2annovar)
		#print stat
		# inhousedict = project.inhouse(otherdatabase)
		snp150dict = project.dictavsnp150(outdir,samplename)
		all1000gdict = project.dict1000gall(outdir,samplename)
		all1000easgdict = project.dict1000geas(outdir,samplename)
		ljbdict = project.dictljb(outdir,samplename)
		#cybdict = project.dictcyb(outdir,samplename)
		variantdict	= project.dictvariant(outdir,samplename)
		exonicdict = project.dictexonic(outdir,samplename)
		espdict = project.dictesp(outdir,samplename)
		#espeadict = project.dictespea(outdir,samplename)
		gnomaddict = project.dictgnomad(outdir,samplename)
		gnomad_exomedict = project.dictgnomad_exome(outdir,samplename)
		seconddict = project.dictsecond(outdir,samplename)
		omimdict = project.dictomim(outdir,samplename)
		pLIdict = project.dictpLI(outdir,samplename)
		Decipherdict = project.dictDecipher(outdir,samplename)
		gene_introdict = project.dictgene_intro(outdir,samplename)
		clinvardict = project.dictclinvar(outdir,samplename)
		hpodict = project.dicthpo(outdir,samplename)
		#intervardict = project.dictintervar(outdir,samplename)
		#cosmicdict = project.dictcosmic(outdir,samplename)
		exacdict = project.dictexac(outdir,samplename)
		reveldict = project.dictrevel(outdir,samplename)
		annotate_m = project.MergeAll(outdir,samplename,exonicdict,snp150dict,all1000gdict,all1000easgdict,espdict,variantdict,ljbdict,gnomaddict,gnomad_exomedict,seconddict,omimdict,pLIdict,Decipherdict,gene_introdict,clinvardict,hpodict,reveldict,exacdict)
		#print inhousedict

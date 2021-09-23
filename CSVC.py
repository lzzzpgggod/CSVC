def mkdir(path):
        import os
        path=path.rstrip("/")
        isExists=os.path.exists(path)
        if not isExists:
                os.makedirs(path)
                #print (path+'succeed')
                return True
        else:
                #print (path+'failed')
                return False

def distinct(drugs):
    #distinct_drugs = drugs
    drugs = set(drugs)
    distinct_drugs = drugs
    drugs = list(drugs)
    print (drugs)
    for i in range(0,len(drugs)):
        for j in range(0,len(drugs)):
            if drugs[i] in dic_metadrug and drugs[j] in dic_metadrug and i > j:
                if dic_metadrug[drugs[i]] == dic_metadrug[drugs[j]]:
                    if drugs[i] in distinct_drugs:
                        distinct_drugs.remove(drugs[i])
    return distinct_drugs
 

def check_pub(gene,dis,drug):
    metadis_id = dic_metadis[dis]
    metagene_id = dic_metagene[gene]
    if drug in dic_metadrug:
        metadrug_id = dic_metadrug[drug]
    else:
        return 0
    if metadis_id in dic_dis2pmid and metagene_id in dic_gene2pmid and metadrug_id in dic_drug2pmid:
        pmid_dis = set(dic_dis2pmid[metadis_id])
        pmid_gene = set(dic_gene2pmid[metagene_id])
        pmid_drug = set(dic_drug2pmid[metadrug_id])
        pmid = pmid_dis & pmid_gene & pmid_drug
        if pmid:
            print(drug,'win')
            return 1
        else:
            return 0
    else:
        return 0

def check_ct(gene,dis,drug):
    metadis_id = dic_metadis[dis]
    metagene_id = dic_metagene[gene]
    if drug in dic_metadrug:
        metadrug_id = dic_metadrug[drug]
    else:
        return 0
    
    
    if metadis_id in dic_dis2nctid and metagene_id in dic_gene2nctid and metadrug_id in dic_drug2nctid:
        pmid_dis = set(dic_dis2nctid[metadis_id])
        pmid_gene = set(dic_gene2nctid[metagene_id])
        pmid_drug = set(dic_drug2nctid[metadrug_id])
        pmid = pmid_dis & pmid_gene & pmid_drug
        #print (dic_dis2nctid[metadis_id],'zzzzzz')
        if pmid:
            print(drug,'win')
            return 1
        else:
            return 0   
    else:
        return 0

def sorting_el(drugs):
    if drugs == '-':
        return drugs
    else:
        dic_drug = {}
        dic_el = {}
        drug = drugs.split(';')
        num = 0
        sorted_drugs = []
        for i in drug:
            #print (i)
            j = i.split('(')
            if j:
                el = '(' + j[1]
                dic_el[num] = el
                dic_drug[num]=j[0]
                num = num +1
        for k in sorted(dic_el,key=dic_el.__getitem__):
            drugel = dic_drug[k] + dic_el[k]
            sorted_drugs.append(drugel)
        drugs = ';'.join(sorted_drugs)
        return drugs

def sorting_drug(line):
    on = line[8]
    off = line[10]
    re = line[-2]
    print (line)
    print ('zzz')
    if on != '-' and on != '':
        print (on)
        on = sorting_el(on)
    if off != '-' and off != '':
        off = sorting_el(off)
    if re != '_' and re != '':
        re = sorting_el(re)
    line[8] = on
    line[10] = off
    line[-2] = re
    return line

def filter_blacklistdrug(dic_line,black_list):
    dic_newline = {}
    for key,value in dic_line.items():
        line = key.split('\t')
        if value == 'mgc':
            if line[4] not in black_list:
                dic_newline[key] = value
        if value == 'civic':
            filter_drug_c = []
            drugs_civic = line[6].split(',')
            for i in drugs_civic:
                if i  not in black_list:
                    filter_drug_c.append(i)
            drugs_civic = ','.join(filter_drug_c)
            ey = key.replace(line[6],drugs_civic)
            dic_newline[key] = value
        if value == 'oncokb':
            drugs_oncokb = line[8].split(', ')
            filter_drug_o = []
            for i in drugs_oncokb:
                if i  not in black_list:
                    filter_drug_c.append(i)
            drugs_civic = ','.join(filter_drug_c)
            ey = key.replace(line[8],drugs_civic)
            dic_newline[key] = value
    return dic_line

def check_el(drug,el,es,type,label):
    if drug in dic_metadrug:
        metadrugid = dic_metadrug[drug]
    else:
        metadrugid = ''
    pgxel = ''
    if metadrugid in dic_druginfo:

        line = dic_druginfo[metadrugid]
        if 'FDA approved' in line :
            if 'CFDA approved' in line:
                phase = '1'
            else:
                phase = '2'
        else:
            if 'CFDA approved' in line:
                phase = '3'
            else:
                phase = '4'
    else:
        phase = '4'

    if type == 'mcg':
        if es != 'Resistance or Non-Response':
            if label == '1':
                if el == '5':
                    pgxel = '1'
                if el == '4':
                    pgxel = '2' + 'A' + phase
                if el == '3':
                    pgxel = '3' + 'A' + phase
            else:
                if el == '5':
                    pgxel = '2' + 'B' + phase
                if el == '4':
                    pgxel = '2' + 'B' + phase
                if el =='3':
                    pgxel = '3' + 'B' + phase
        else:
            if el == '5':
                pgxel = 'R1'
            else:
                pgxel = 'R2'


    if type == 'civic':
        if es != 'Resistance':
            if label == '1':
                if el == 'A':
                    pgxel = '2'+ 'A' + phase
                if el == 'B' or el == 'C':
                    pgxel = '3' + 'A' + phase
                if el == 'D' or el == 'E':
                    pgxel = '4'
            else:
                if el == 'A':
                    pgxel = '2'+ 'B' + phase
                if el == 'B' or el == 'C':
                    pgxel = '3' + 'B' + phase
                if el == 'D' or el == 'E':
                    pgxel = '4'
        else:
            pgxel = 'R2'


    if type == 'oncokb':
        #(drugs[0],key[7],'','oncokb','2')
        if label == '1':
            if el ==  '1' or el =='R1' or el == 'R2' or el == '4':
                pgxel = el
            else:
                pgxel = el + phase
        else:
            if el == '1' or el == '2A':
                pgxel = '2B' + phase
            if el == '3A':
                pgxel = '3B' + phase
            if el == '4':
                pgxel = '4'
    if pgxel != '':
        pgxel = '('+pgxel + ')'
        return pgxel
    else:
        return 0
            
def compare_el(el1,drug1,drug):
    match = re.search(r'(\(\w*\))' , drug1)
    if match:
        el2 = match.group(1)
    
    el3 = el1[1:-1]
    el4 = el2[1:-1]
    if 'R' not in el1 and 'R' not in el2:
        no1 = float(el3[0])
        no2 = float(el4[0])
        if no1<=no2:
            el1 = drug + el1
            return el1
        else:
            el2 = drug + el2
            return el2
    else:
        if 'R' in el1:
            el1 = drug + el1
            return el1
        else:
            el2 = drug + el2
            return el2

#check distree
def check_distree(dis,drug2dis):
    #print (dis,drug2dis)
    #print (dic_dis2tree[dis],dic_dis2tree[drug2dis])
    if dis in dic_dis2tree and drug2dis in dic_dis2tree:
        #
        if dic_dis2tree[drug2dis] in dic_dis2tree[dis] or dic_dis2tree[dis] in dic_dis2tree[drug2dis]:
            #print(dic_dis2tree[drug2dis])
            return 1
        else:
            a_distreeid = dic_dis2tree[dis].split('|')
            a_drug2distree = dic_dis2tree[drug2dis].split('|')
            
            #a_detailofdistreeid
            #a_treeid = drug2distree.split('\|')
            for i in a_distreeid:
                for j in a_drug2distree:
                    a_detailofdistreeid = i
                    a_detailofdrug2distree = j
                    if i in j or j in i :
                        print (a_detailofdistreeid,a_detailofdrug2distree)
                        print (a_distreeid,a_drug2distree)
                        return 1

    else:
        return 2

def check_disease(dis,info,type):
    line = info.split('\t')
    #print(dis)
    drugin = 0
    if type == 'mcg':
        if line[4]:
            dis_meta = dic_metadis[dis]
            if line[9] in dic_metadis:
                info_meta = dic_metadis[line[9]]
            else:
                return 0
            if dis_meta == info_meta:
                onlabeldrug.append(line[4])
                drugin = 1
                dic_acd2line[info] = 'mcg'
            else:
                re = check_distree(dis_meta,info_meta)
                if re == 1:
                    drugin = 1
                    onlabeldrug.append(line[4])
                    dic_acd2line[info] = 'mcg'
                else:
                    #ct = check_ct(line[2],dis,line[4])
                    #pub = check_pub(line[2],dis,line[4])
                    #if ct == 1 or pub == 1:
                    
                    offlabeldrug.append(line[4])
                    if 'Resistance' not in info:
                        drugin = 1
                        dic_offlabel[info] = 'mcg'
            if drugin == 1:
                if 'Resistance' in info:
                    redrug.append(line[4])
                else:
                    sedrug.append(line[4])
                    
            

    if type == 'civic':
        #print (line[3])
        
        dis_meta = dic_metadis[dis]
        if line[3] in dic_metadis:
            info_meta = dic_metadis[line[3]]
        else:
            return 0
        if line[6]:
            drugs = line[6].split(',')
            if dis_meta == info_meta:
                onlabeldrug.append(line[6])
                dic_acd2line[info] = 'civic'
                drugin = 1
            else:
                re = check_distree(dis_meta,info_meta)
                if re == 1:
                    onlabeldrug.append(line[6])
                    dic_acd2line[info] = 'civic'
                    drugin = 1
                else:
                    if line[6]:
                        for i in drugs:
                            #ct = check_ct(line[0],dis,i) 
                            #pub = check_pub(line[0],dis,i)
                            #if ct == 1 or pub == 1:
                            #offlabeldrug.append(line[6])
                            if 'Resistance' not in info:
                                dic_offlabel[info] = 'civic'
                                drugin = 1
                            break
            if drugin == 1:
                for k in drugs:
                    if 'Resistance' in info:
                        redrug.append(k)
                    else:
                        sedrug.append(k)

    if type == 'oncokb':
        #print (line[6])
        if line[8]:
            drugs = line[8].split(', ')
            dis_meta = dic_metadis[dis]
            if line[6] in dic_metadis:
                info_meta = dic_metadis[line[6]]
            else:
                return 0
            if dis_meta == info_meta:
                onlabeldrug.append(line[8])
                dic_acd2line[info] = 'oncokb'
                drugin = 1
            else:
                re = check_distree(dis_meta,info_meta)
                if re == 1:
                    onlabeldrug.append(line[8])
                    dic_acd2line[info] = 'oncokb'
                    drugin = 1
                else:
                    # ,
                    drugs = line[8].split(', ')
                    for i in drugs:
                        #ct = check_ct(line[3],dis,i) 
                        #pub = check_pub(line[3],dis,i)
                        #if ct == 1 or pub == 1:
                        offlabeldrug.append(line[8])
                        if 'Resistance' not in info:
                            dic_offlabel[info] = 'oncokb'
                            drugin = 1
                        break
            if drugin == 1:
                for k in drugs:
                    if 'Resistance' in info:
                        redrug.append(k)
                    else:
                        sedrug.append(k)

def check_mutation(line,db_type):
    #print ('1')
    #print (line)
    info = line.split('\t')
    if db_type == 'mcg':
        if info[0] in dic_var2info:
            #print(1)
            line = line + dic_var2seqinfo[info[0]]
            return line
        
        if '-|-' in info[0]:
            if info[2] in dic_pathogenic_geneinfo:
                line = line + '\t' + '*' + '\t' + dic_pathogenic_geneinfo[info[2]]
                return line


    if db_type == 'civic':
        ref = info[26]
        alt = info[27]
        if ref == '':
            ref = '-'
        if alt == '':
            alt = '-'
        
        var = info[23] + '|' + info[24] + '|' + ref + '|' + alt
        if var in dic_var2info:
            #print (1)
            #print (line)
            line = line + dic_var2seqinfo[var]
            return line 
        if 'MUTATION' in info[2]:
            if info[0] in dic_pathogenic_geneinfo:
                line = line + '\t' + '*' + '\t' + dic_pathogenic_geneinfo[info[0]]
                return line


    if db_type == 'oncokb':
        
        if info[3] in dic_gene2aachange:
            if info[5] in dic_gene2aachange[info[3]]:
                #print (info[3],info[5])
            #candidata
                for i in dic_geneaa2seqinfo:
                    #print (i)
                    if info[3] in i and info[5] in i:
                        line = line + dic_geneaa2seqinfo[i]
                        return line
        if 'Mutations' in info[5]:
            if info[3] in dic_oncokb_oncogene:
                line = line + '\t' + '*' + '\t' + dic_oncokb_oncogene[info[3]]
                return line
    return 0

def is_pathogenic(score):
    if score > 4:
        return 1
    else :
        return 0
#dis
import os
file_names = os.listdir("./testing/PGxVar/")
file_list = [os.path.join("./testing/PGxVar/",files) for files in file_names]
drugofall = set()
for file_index in range(0,len(file_list)):
    if 'PGxVar' in file_list[file_index]:
        import os
        import re
        matchobj = re.search(r'(\d+)K_mutation' ,file_names[file_index])
        print (file_names[file_index])
        dir_name = matchobj.group(1)
        path_output = './testing/re/' + dir_name
        mkdir(path_output)

        dis = 'Breast Cancer'
        t_rowname = ('Gene','Isoform','ID','gene','Outname','MetaDrugID','Variation','PGxDiseaseID','OfficialDiseaseName','PGI_drugid_tmp',)

        oncokb_oncogene = []
        with open('./oncokb_2019july/allAnnotatedVariants.txt','r') as file:
            for line in file:
                if 'Oncogenic' in line:
                    info = line.split('\t')
                    gene_var = info[5] + '\t' + info[3]
                    oncokb_oncogene.append(gene_var)


        dic_var2info = {}
        dic_gene2aachange = {}
        pathogenic = []
        #output = open("out.txt", "w") 
        aa = []
        gene_list = []
        #import vcf
        dic_oncokb_oncogene = {}
        dic_var2seqinfo = {}
        dic_geneaa2seqinfo = {}
        dic_pathogenic_geneinfo = {}
        with open(file_list[file_index],'r') as file:
            for line in file:
                line = line.replace('\n','')
                if 'normal_sample' in line :
                    matchobj = re.search(r'normal_sample=(\w*)' ,line)
                    normal_name = matchobj.group(1)
                if 'tumor_sample' in line :
                    matchobj = re.search(r'tumor_sample=(\w*)' ,line)
                    tumor_name = matchobj.group(1) 
                if '#CHROM' in line:
                    vcf = line.split('\t')
                    if vcf[9] == tumor_name:
                        tumor_col = 9
                        normal_col = 10
                    else:
                        tumor_col = 10
                        normal_col = 9


                if line[0] == '#':
                    continue
                
                vcf = line.split('\t')
                #print (vcf[7])
                #print(vcf[10])
                #filter quality
                #
                if vcf[6] :
                    chr = vcf[0]
                    #ExonicFunc.refGen
                    info = vcf[7].split(';')
                    #filter unchange protein 
                    if 'AAChange.refGene=.;' not in line or 'splicing' in line:
                        #fix the format
                        ref = vcf[3]
                        alt = vcf[4]
                        if ref == '.':
                            ref = '-'
                        if alt == '.':
                            alt = '-'
                        if chr[0] != 'C' and chr[0] != 'c' :
                            var = chr + '|' + vcf[1] + '|' + ref + '|' + alt
                        else:
                            var = chr[3:] + '|' + vcf[1] + '|' + ref + '|' + alt
                        #print (var)
                        dic_var2info[var] = info
                        #gene
                        matchobj = re.search(r'Gene.refGene=(\w*);' , vcf[7])
                        if matchobj:
                            gene = matchobj.group(1)
                        else:
                            gene = '.'
                        
                        gene_list.append(gene)

                        fuc = re.search(r'Func.refGene=(\w*);' , vcf[7])
                        if fuc:
                            funcgene = fuc.group(1)
                        else:
                            funcgene = '.'
                        #sequening information
                        matchobj = re.search(r'CLINSIG=(.+);CLNDBN',vcf[7])
                        if matchobj:
                            clinsig = matchobj.group(1)
                        matchobj = re.search(r'CLNDBN=(.+);CLNACC',vcf[7])
                        if matchobj:
                            clinbn = matchobj.group(1)
                        matchobj = re.search(r'AAChange.refGene=\w+:(\w+:\w+:c.\w+:p.\w+)(;|,)',vcf[7])
                        if matchobj:
                            aac = matchobj.group(1)
                        else:
                            aac = ''
                        #print (aac)
                        if 'splicing' in line:
                            fuc = 'splicing'
                        else:
                            matchobj = re.search(r'ExonicFunc.refGene=(.+);AAChange.refGene',vcf[7])
                            if matchobj:
                                fuc = matchobj.group(1)
                        


                        #find aa change
                        for i in info:
                            if 'AAChange.refGene=' in i:
                                AAChange = i
                                if gene not in dic_gene2aachange:
                                    dic_gene2aachange[gene] = []
                                matchobj = re.search(r':p\.(\w+)' , i)
                                if matchobj:
                                    alter = matchobj.group(1)
                                matchc = re.search(r':c\.(\w+)' , i)
                                if matchc:
                                    calter = matchc.group(1)
                                    #print (gene,alter)

                                dic_gene2aachange[gene].append(alter)
                                geneaa = alter + '\t' + gene
                            
                        #PGxVar
                        if len(vcf) > 10 :
                            score = float(vcf[11])
                        if len(vcf) < 10 :
                            score = float(vcf[8])
                        
                        seq_info = ''
                        if len(vcf) > 10:
                            match = re.search(r'PopFreqMax=(\w+)',vcf[7])
                            if match:
                                popfre = match.group(1)
                            else:
                                popfre = '.'
                            normal = vcf[9].split(':')
                            tumor = vcf[10].split(':')

                            calter = 'c.' + calter
                            alter = 'p.' + alter
                            var2 = var.replace("|","-")
                            seq_info ='\t'+ fuc + '\t' + aac + '\t' + clinbn + '\t' + clinsig + '\t'+ gene +'\t'+ calter + '\t'+ alter +'\t' + str(score) + '\t' + funcgene +'\t' + popfre + '\t' + normal[2] + '\t' + normal[3] + '\t' + tumor[2] + '\t' +tumor [3] + '\t' + var2 
                        
                        if is_pathogenic(score) == 1:
                            dic_pathogenic_geneinfo[gene] = seq_info
                        
                        if geneaa in oncokb_oncogene:
                            dic_oncokb_oncogene[gene] = seq_info

                        dic_geneaa2seqinfo[geneaa] = seq_info
                        dic_var2seqinfo[var] = seq_info
        



        dic_line = {}
        #dic_oncokb2line = {}
        with open('./oncokb_2019july/allActionableVariants.txt','r',encoding='utf-8') as file:
            for line in file:
                onco_type = 'oncokb' 
                info = line.split('\t')
                if info[0] in t_rowname:
                    continue
                
                #print(line)
                cm_info = check_mutation(line,onco_type)
                if cm_info  != 0:
                    dic_line[cm_info] = 'oncokb'
                
                    
        civic_line = 0
        dic_civicline2evs = {}
        with open('./something/evidence_statement.txt','r',encoding='utf-8') as file:
            for line in file:
                line = line.replace('\n','')
                info = line.split('\t')
                dic_civicline2evs[info[1]] = info[2]
        
        with open('./something/evidence_statement_extra.txt','r',encoding='utf-8') as file:
            for line in file:
                line = line.replace('\n','')
                info = line.split('\t')
                dic_civicline2evs[info[0]] = info[1]

        with open('./civic_2019july01/nightly-ClinicalEvidenceSummaries.tsv','r',encoding='utf-8') as file:
            for line in file:
                civic_type = 'civic'
                info = line.split('\t')
                if info[0] in t_rowname:
                    continue
                cm_info = check_mutation(line,civic_type)
                if cm_info  != 0:
                    dic_line[cm_info] = 'civic'
                    

        dic_metagene = {}                                                                                                 
        with open('./mcg/vg2dr_hg19.txt','r',encoding='utf-8') as file:
            for line in file:
                mcg_type = 'mcg'
                info = line.split('\t')
                if info[0] in t_rowname:
                    continue
                #print (info)
                dic_metagene[info[2]] = info[3]
                cm_info = check_mutation(line,mcg_type)
                if cm_info  != 0:
                    dic_line[cm_info] = 'mcg'
                dic_metagene[info[2]] = info[3] 


        with open('./something/MetaGeneIDSymbol.txt','r',encoding='utf-8') as file:
            for line in file:
                info = line.split('\t')
                dic_metagene[info[1]] = info[0]



        dic_dis2tree={}
        with open('./tmp_import/distable.txt','r',encoding='utf-8') as file:
            for line in file:
                line = line.replace('\n','')
                info = line.split('\t')
                if info[0] in t_rowname:
                    continue
                dic_dis2tree[info[2]] = info[1]

        dic_druginfo2 = {}
        dic_druginfo = {}
        dic_metadrug = {}
        with open('./tmp_import/drugnames.txt','r') as file:
            for line in file:
                line = line.replace('\n','')
                info = line.split('\t')
                dic_metadrug[info[1]] = info[2]

        with open('./something/zz.txt','r',encoding='utf-8') as file:
            for line in file:
                info = line.split('\t')
                #print(info[2])
                #2,3,4,5name,us,cn,descrip
                dic_metadrug[info[2]] = info[0]
                haha = info[3]+'\t'+info[4]+'\t'+info[5]
                dic_druginfo2[info[0]] = haha
                dic_druginfo[info[0]] = line

        dic_metadis = {}
        with open('./tmp_import/diseasename.txt','r') as file:
            for line in file:
                line = line.replace('\n','')
                info = line.split('\t')
                dic_metadis[info[1]] = info[0]

        #pub
        dic_dis2pmid = {}
        with open('./tmp_import/metadisease2pubtator1.txt','r') as file:
            for line in file:
                line = line.replace('\n','')
                info = line.split('\t')
                if info[3] not in dic_dis2pmid:
                    dic_dis2pmid[info[3]] = []
                dic_dis2pmid[info[3]].append(info[2])

        dic_drug2pmid = {}
        with open('./tmp_import/metadrug2pubtator1.txt','r',encoding='utf-8') as file:
            for line in file:
                line = line.replace('\n','')
                info = line.split('\t')
                if info[3] not in dic_drug2pmid:
                    dic_drug2pmid[info[3]] = []
                dic_drug2pmid[info[3]].append(info[2])

        dic_gene2pmid = {}
        with open('./tmp_import/metagene2pubtator1.txt','r',encoding='utf-8') as file:
            for line in file:
                line = line.replace('\n','')
                info = line.split('\t')
                if info[3] not in dic_gene2pmid:
                    dic_gene2pmid[info[3]] = []
                dic_gene2pmid[info[3]].append(info[2])

        #clinical trial
        dic_dis2nctid = {}
        with open('./tmp_import/CT2PGxDisease1.txt','r',encoding='utf-8') as file:
            for line in file:
                line = line.replace('\n','')
                info = line.split('\t')
                if info[2] not in dic_dis2nctid:
                    dic_dis2nctid[info[2]] = []
                dic_dis2nctid[info[2]].append(info[1])

        dic_drug2nctid = {}
        with open('./tmp_import/CT2PGxDrug1.txt','r',encoding='utf-8') as file:
            for line in file:
                line = line.replace('\n','')
                info = line.split('\t')
                if info[2] not in dic_drug2nctid:
                    dic_drug2nctid[info[2]] = []
                dic_drug2nctid[info[2]].append(info[1])

        dic_gene2nctid = {}
        with open('./tmp_import/CT2PGxGene1.txt','r',encoding='utf-8') as file:
            for line in file:
                line = line.replace('\n','')
                info = line.split('\t')
                if info[2] not in dic_gene2nctid:
                    dic_gene2nctid[info[2]] = []
                dic_gene2nctid[info[2]].append(info[1])

        black_list = []
        with open('./blacklist/targetdrug_blacklist.txt','r',encoding='utf-8') as file:
            for line in file:
                info = line.split('\t')
                name = info[0].capitalize()
                black_list.append(name)

        with open('./blacklist/drug_blacklist.txt','r',encoding='utf-8') as file:
            for line in file:
                line = line.replace('\n','')
                info = line.split('\t')
                black_list.append(info[0])


        onlabeldrug = []
        redrug = []
        sedrug = []
        offlabeldrug = []
        dic_acd2line = {}
        dic_offlabel = {}

        d = open("all.txt", "w",encoding='utf-8')
        for key,value in dic_line.items():
            check_disease(dis,key,value)
            key = key.replace('\n','\t')
            key = key + '\n'
            d.write(key)
        
        sedrug = set(sedrug)
        redrug = set(redrug)
        
        print (sedrug)
        print (redrug)
        summary = sedrug - redrug
        drug_array = redrug | sedrug
        drugofall = drugofall | drug_array

        summary_name = './testing/re/'  + dir_name + '_summary.txt'
        druginfo_name = './testing/re/'  +   dir_name + '_druginfo.txt'
        varinfo_name = './testing/re/'  +  dir_name + '_varinfo.txt'
        detail_name ='./testing/re/'  + dir_name +'_detail.txt'

        f = open(summary_name, "w",encoding='utf-8')
        dr = open(druginfo_name, 'w',encoding='utf-8')
        var_seq = open(varinfo_name, 'w',encoding='utf-8')
        drugvar_info = open(detail_name, 'w',encoding='utf-8')
        summary_ti = '基因名' + '\t'  + '突变'  + '\t' +  '在病人所患癌症中可能获益的治疗方案'  + '\t' +  '在其它癌症中可能获益的治疗方案'  + '\t' + '可能不获益的治疗方案' + '\n'
        f.write(summary_ti)
        druginfo = '药物名' + '\t'  + '中文名' + '\t'  + '药物批准状态' + '\t'  + '药物具体介绍' + '\n'
        dr.write(druginfo)
        seqinfo = '基因名' + '\t'+ '突变' + '\t' + '蛋白质序列改变' + '\t' + '突变功能' + '\t'  + '人群基因型频率' + '\t' + '肿瘤组织的测序深度和等位基因频率' + '\t'  + '正常对照的测序深度和等位基因频率' + '\t'  + 'CLINSIG' + '\t' + 'CLNDBN' + '\t' + ' PGxVar' + '\n'
        var_seq.write(seqinfo)
        drugvarseqinfo = '基因名' + '\t' + '突变'  +'\t' + '疾病名' + '\t' +  '药物名'  + '\t' +'证据说明' + '\t' + '文献支持' + '\t'  + '来源' + '\n'
        drugvar_info.write(drugvarseqinfo)
        #whitelist

        metadis = dic_metadis[dis]
        redismeta = dic_metadis['Colorectal Cancer']
        rep = check_distree(metadis,redismeta)
        if rep == 1 and 'KRAS' not in gene_list:
            whitesum = 'KRAS' + '\t' + '野生型' + '\t' + 'Regorafenib(1);Panitumumab(1);Cetuximab(1)' + '\t' + '-' + '\t' + '-' + '\n'
            whiteseq = 'KRAS' + '\t' + '野生型' + '\t' + '.' + '\t' + '.' + '\t'  + '.' + '\t' + '.' + '\t'  + '.' + '\t'  + '.' + '\t' + '.' + '\t' + '.' + '\n'
            f.write(whitesum)
            var_seq.write(whiteseq)
            with open('./tmp_import/whitelist_hg19.txt','r',encoding='utf-8') as file:
                for line in file:
                    if line[0] != 'V':
                        drug_array.add(info[4])
                        info = line.split('\t')
                        whitevar = info[2] + '\t' + '野生型' + '\t' + info[4] + '(1)' + '\t' + info[10] + '\t' + '-' + '\t' + info[11] + '\n'
                        drugvar_info.write(whitevar)


        #analysis the dru
        dic_ondrug = {}
        dic_detail = {}
        dic_offdrug = {}
        dic_seqinfo = {}
        #print (summary)
        onlabeldrug = []
        for key,value in dic_acd2line.items():
            
            key = key.replace('\n','')
            k = key
            key = key.split('\t')
            #print (key)
            #seqinfo
            if key[-1] not in dic_seqinfo:
                dic_seqinfo[key[-1]] = 1
                if float(key[-8]) >6:
                    pgxvar = 'Pathogenic'
                if float(key[-8]) <=6 and float(key[-8])>4:
                    pgxvar = 'Likely_pathogenic'
                if float(key[-8]) <=4 and float(key[-8])>3:
                    pgxvar = 'Uncertain_significance'
                if float(key[-8]) <=3 and float(key[-8])>2:
                    pgxvar = 'Likely_benign'
                if float(key[-8]) <=2:
                    pgxvar = 'Benign'
                linee = key[-11] + '\t' + key[-1] +'\t'+ key[-14] + '\t' + key[-15] + '\t' + key[-6] + '\t' + key[-5] +  ';' + key[-4] +  '\t' + key[-3] +  ';' + key[-2] + '\t' + key[-12] + '\t' + key[-13]+ '\t' +pgxvar + '\n'
                var_seq.write(linee)
            if value == 'mcg':
                #print (key[7])
                #print (key)
                el = check_el(key[4],key[6],key[7],'mcg','1')
                #print(key[4],key[6],key[7],'1',el,'MCG')
                if key[4] in black_list:
                    continue
                #fuc15+aac14+clinbn13+ clinsig12 +gene11 + calter10 + alte9r +str(score)8 +funcgene7 + popfre6  + normal[2]5  + normal[3]4 + tumor[2]3  +tumor [3] 2 + var1 
                drug_el = key[4]+el
                    #'基因名' + '\t' + '突变'  +'\t' + '药物名' + '\t'  + '药物反应' + '\t' +'证据说明' + '\t' + '文献支持' + '\t'  + '来源' + '\n'
                
                #print (key[4])

                if key[-17] != '*':
                    if key[4] in summary or 'R' in el:
                        vgline = key[-11] + '\t' + key[-14] + '\t' + key[9] + '\t' + drug_el   + '\t' + key[10] + '\t' + '-' + '\t' + 'MyCancerGenome' + '\n'
                        drugvar_info.write(vgline)
                        onlabeldrug.append(key[4])
                    #'基因名' + '\t'  + '突变'  + '\t' +  '在病人所患癌症中可能获益的药物'  + '\t' +  '在其它癌症中可能获益的药物'  + '\t' + '抵抗的药物' + '\n'
                        print(drug_el)
                        if key[-1] not in dic_ondrug: 
                            sum = [key[-11], '\t' ,key[-1],';',key[-9],';',key[-10],'\t',drug_el,'\t' ,'-','\t','-','\n'] 
                            dic_ondrug[key[-1]] = sum
                        else:
                            same = 0
                            sum_drug = dic_ondrug[key[-1]][8].split(';')
                            for sd in sum_drug:
                                if key[4] in sd:
                                    drug_el = compare_el(el,sd,key[4])
                                    dic_ondrug[key[-1]][8] = dic_ondrug[key[-1]][8].replace(sd,drug_el)
                                else:
                                    drug2 = sd.split('(')
                                    if drug2[0] in dic_metadrug and key[4] in dic_metadrug:
                                        #
                                        if dic_metadrug[drug2[0]] == dic_metadrug[key[4]]:
                                            drug_el = compare_el(el,sd,key[4])
                                            dic_ondrug[key[-1]][8] = dic_ondrug[key[-1]][8].replace(sd,drug_el)
                                        else:
                                            same = same + 1
                                    else:
                                        same = same + 1
                            if same == len(sum_drug):
                                dic_ondrug[key[-1]][8] = dic_ondrug[key[-1]][8] + ';' + drug_el

            if value == 'civic':
                #print (k)
                drugs = key[6].split(',') 
                el = check_el(drugs[0],key[10],key[11],'civic','1')
                #print(drugs,key[10],key[11],'1',el,'ci')
                if key[6] in black_list:
                    continue
                drugs = key[6].split(',') 
                #print(key)
                #print (drugs[0],key[10],key[11],'civic','1')
                #el = check_el(drugs[0],key[10],key[11],'civic','1')
                #print(drugs,key[10],key[11],'1',el,'ci')
                link = el+';'
                drug_el = key[6].replace(',',link)
                drug_el = drug_el + el
                
                
                if key[-17] != '*':
                    #print (k)
                    vgline = key[-11] + '\t' + key[-14] + '\t' + key[3] + '\t' + drug_el   + '\t'  + dic_civicline2evs[key[12]] + '\t' + key[14] + ':' + key[13] + '\t' + 'CIViC' + '\n'
                    drugvar_info.write(vgline)
                    for drug in drugs:
                        if drug in summary or 'R' in el :
                            onlabeldrug.append(drug)
                            drug1_el = drug + el
                            #print(drug1_el,'zzz')
                            if key[-1] not in dic_ondrug:  
                                sum = [key[-11], '\t' ,key[-1],';',key[-9],';',key[-10],'\t',drug_el,'\t' ,'-','\t','-','\n']
                                #print (sum)
                                dic_ondrug[key[-1]] = sum
                            else:
                                same = 0
                                sum_drug = dic_ondrug[key[-1]][8].split(';')
                                for sd in sum_drug:
                                    if drug in sd :
                                        drug_el = compare_el(el,sd,drug)
                                        dic_ondrug[key[-1]][8] = dic_ondrug[key[-1]][8].replace(sd,drug_el)
                                    else:
                                        drug2 = sd.split('(')
                                        if drug in dic_metadrug and drug2[0] in dic_metadrug:
                                            # 
                                            if dic_metadrug[drug2[0]] == dic_metadrug[drug]:
                                                drug_el = compare_el(el,sd,drug)
                                                dic_ondrug[key[-1]][8] = dic_ondrug[key[-1]][8].replace(sd,drug_el)
                                            else:
                                                same = same + 1
                                        else:
                                            same = same +1
                                print (same,len(sum_drug),sum_drug,'cic')
                                if same == len(sum_drug):
                                    dic_ondrug[key[-1]][8] = dic_ondrug[key[-1]][8] + ';' + drug1_el

                    

            if value == 'oncokb':
                if key[8] in black_list:
                    continue
                drugs = key[8].split(', ')
                el = check_el(drugs[0],key[7],'','oncokb','1')
                #print(drugs,key[7],'','1',el,'oncokb')
                link = el+';'
                drug_el = key[8].replace(', ',link)
                drug_el = drug_el + el
                
                
                if key[-17] != '*':
                    vgline = key[-11] + '\t' + key[-14] + '\t' + key[6] + '\t' + drug_el   + '\t'  + '-' + '\t' + 'PubMed:' + key[9] + '\t' + 'OncoKB' + '\n'
                    drugvar_info.write(vgline)
                    for drug in drugs:
                        if drug in summary or 'R' in el:
                            onlabeldrug.append(drug)
                            drug1_el = drug + el
                            if key[-1] not in dic_ondrug:  
                                
                                sum = [key[-11], '\t' ,key[-1],';',key[-9],';',key[-10],'\t',drug_el,'\t' ,'-','\t','-','\n']
                                dic_ondrug[key[-1]] = sum
                            else:
                                same = 0
                                sum_drug = dic_ondrug[key[-1]][8].split(';')
                                for sd in sum_drug:
                                    if drug in sd :
                                        drug_el = compare_el(el,sd,drug)
                                        dic_ondrug[key[-1]][8] = dic_ondrug[key[-1]][8].replace(sd,drug_el)
                                    else:
                                        drug2 = sd.split('(')
                                        if drug in dic_metadrug and drug2[0] in dic_metadrug:
                                            # 
                                            if dic_metadrug[drug2[0]] == dic_metadrug[drug]:
                                                drug_el = compare_el(el,sd,drug)
                                                dic_ondrug[key[-1]][8] = dic_ondrug[key[-1]][8].replace(sd,drug_el)
                                            else:
                                                same = same + 1
                                        else:
                                            same = same +1
                                if same ==len(sum_drug):
                                    dic_ondrug[key[-1]][8] = dic_ondrug[key[-1]][8] + ';' + drug1_el
                                                

            

        print(onlabeldrug,'zzzz')
        for key,value in dic_offlabel.items():
            key = key.replace("\n","")
            k = key + '\n'
            #d.write(k)
            key = key.split('\t')
            #seqinfo
            #print (key)
            if key[-1] not in dic_seqinfo:
                dic_seqinfo[key[-1]] = 1
                if float(key[-8]) >=6:
                    pgxvar = 'Pathogenic'
                if float(key[-8]) < 6 and float(key[-8])>4:
                    pgxvar = 'Likely_pathogenic'
                if float(key[-8]) <=4 and float(key[-8])>3:
                    pgxvar = 'Uncertain_significance'
                if float(key[-8]) <=3 and float(key[-8])>2:
                    pgxvar = 'Likely_benign'
                if float(key[-8]) <=2:
                    pgxvar = 'Benign'
                linee = key[-11] + '\t' + key[-1] +'\t'+ key[-14] + '\t' + key[-15] + '\t' + key[-6] + '\t' + key[-5] +  ';' + key[-4] +  '\t' + key[-3] +  ';' + key[-2] + '\t' + key[-12] + '\t' + key[-13]+ '\t' +pgxvar + '\n'
                var_seq.write(linee)
            if value == 'mcg':
                #print (key[7])
                if key[4] in black_list:
                    continue
                el = check_el(key[4],key[6],key[7],'mcg','2')
                
                #gene11 + calter10 + alte9r +str(score)8 +  funcgene7 + popfre6  + normal[2]5  + normal[3]4 + tumor[2]3  +tumor [3] 2 + var1 
                drug_el = key[4]+el
                
                if key[-17] != '*':
                    #print (el)
                    if key[4] in summary  and key[4] not in onlabeldrug:
                        vgline = key[-11] + '\t' + key[-14] + '\t' + key[9] + '\t' + drug_el   + '\t' + key[10] + '\t' + '-' + '\t' + 'MyCancerGenome' + '\n'
                        drugvar_info.write(vgline)
                        if key[-1] not in dic_offdrug:  
                            sum = [key[-11], '\t' ,key[-1],';',key[-9],';',key[-10],'\t','-','\t' ,drug_el,'\t','-','\n']
                            dic_offdrug[key[-1]] = sum
                        else:
                            same = 0
                            sum_drug = dic_offdrug[key[-1]][10].split(';')
                            for sd in sum_drug:
                                if key[4] in sd:
                                    drug_el = compare_el(el,sd,key[4])
                                    dic_offdrug[key[-1]][10] = dic_offdrug[key[-1]][10].replace(sd,drug_el)
                                else:
                                    drug2 = sd.split('(')
                                    if drug2[0] in dic_metadrug and key[4] in dic_metadrug:
                                        #
                                        if dic_metadrug[drug2[0]] == dic_metadrug[key[4]]:
                                            drug_el = compare_el(el,sd,key[4])
                                            dic_offdrug[key[-1]][10] = dic_offdrug[key[-1]][10].replace(sd,drug_el)
                                        else:
                                            same = same + 1
                                    else:
                                        same = same + 1
                            if same == len(sum_drug):
                                dic_offdrug[key[-1]][10] = dic_offdrug[key[-1]][10] + ';' + drug_el

            if value == 'civic':
                if key[6] in black_list:
                    continue
                drugs = key[6].split(',')
                #print(drugs[0],key[10],key[11],'civic','2')
                el = check_el(drugs[0],key[10],key[11],'civic','2')
                link = el+';'
                drug_el = key[6].replace(',',link)
                drug_el = drug_el + el
                
                
                if key[-17] != '*':
                    vgline = key[-11] + '\t' + key[-14] + '\t' + key[3] + '\t' + drug_el   + '\t'  + dic_civicline2evs[key[12]] + '\t' + key[14] + ':' + key[13] + '\t' + 'CIViC' + '\n'
                    drugvar_info.write(vgline)
                    for drug in drugs:
                        #print (drug)
                        if drug in summary and drug not in onlabeldrug:
                            
                            drug1_el = drug + el
                            
                            if key[-1] not in dic_offdrug:  
                                sum = [key[-11], '\t' ,key[-1],';',key[-9],';',key[-10],'\t','-','\t' ,drug1_el,'\t','-','\n']
                                dic_offdrug[key[-1]] = sum
                            else:
                                same = 0
                                sum_drug = dic_offdrug[key[-1]][10].split(';')
                                for sd in sum_drug:
                                    if drug in sd :
                                        drug_el = compare_el(el,sd,drug)
                                        dic_offdrug[key[-1]][10] = dic_offdrug[key[-1]][10].replace(sd,drug_el)
                                    else:
                                        drug2 = sd.split('(')
                                        if drug in dic_metadrug and drug2[0] in dic_metadrug:
                                            # 
                                            if dic_metadrug[drug2[0]] == dic_metadrug[drug]:
                                                drug_el = compare_el(el,sd,drug)
                                                dic_offdrug[key[-1]][10] = dic_offdrug[key[-1]][10].replace(sd,drug_el)
                                            else:
                                                same = same + 1
                                        else:
                                            same = same +1
                                #print (same,len(sum_drug),sum_drug,'cic')
                                if same == len(sum_drug):
                                    dic_offdrug[key[-1]][10] = dic_offdrug[key[-1]][10] + ';' + drug1_el
                    

            if value == 'oncokb':
                if key[8] in black_list:
                    continue
                drugs = key[8].split(', ')
                if 'R' in key[7]:
                    continue
                el = check_el(drugs[0],key[7],'','oncokb','2')
                
                link = el + ';'
                drug_el = key[8].replace(', ',link)
                drug_el = drug_el + el
                
                
                if key[-17] != '*':
                    vgline = key[-11] + '\t' + key[-14] + '\t' + key[6] + '\t' + drug_el   + '\t'  + '-' + '\t' + 'PubMed:' + key[9] + '\t' + 'OncoKB' + '\n'
                    drugvar_info.write(vgline)
                    for drug in drugs:
                        #print (drug)
                        if drug in summary and drug not in onlabeldrug:
                            
                            drug1_el = drug + el
                            
                            if key[-1] not in dic_offdrug:  
                                sum = [key[-11], '\t' ,key[-1],';',key[-9],';',key[-10],'\t','-','\t' ,drug1_el,'\t','-','\n']
                                dic_offdrug[key[-1]] = sum
                            else:
                                same = 0
                                sum_drug = dic_offdrug[key[-1]][10].split(';')
                                for sd in sum_drug:
                                    if drug in sd :
                                        drug_el = compare_el(el,sd,drug)
                                        dic_offdrug[key[-1]][10] = dic_offdrug[key[-1]][10].replace(sd,drug_el)
                                    else:
                                        drug2 = sd.split('(')
                                        if drug in dic_metadrug and drug2[0] in dic_metadrug:
                                            # 
                                            if dic_metadrug[drug2[0]] == dic_metadrug[drug]:
                                                drug_el = compare_el(el,sd,drug)
                                                dic_offdrug[key[-1]][10] = dic_offdrug[key[-1]][10].replace(sd,drug_el)
                                            else:
                                                same = same + 1
                                        else:
                                            same = same +1
                                if same ==len(sum_drug):
                                    dic_offdrug[key[-1]][10] = dic_offdrug[key[-1]][10] + ';' + drug1_el
            
            


#predict
        for key,value in dic_acd2line.items():
            key = key.replace('\n','')
            key = key.split('\t')
            if value == 'mcg':
                #print (key[7])
                if key[4] in black_list:
                    continue
                el = check_el(key[4],key[6],key[7],'mcg','1')
                #fuc15+aac14+clinbn13+ clinsig12 +gene11 + calter10 + alte9r +str(score)8 +funcgene7 + popfre6  + normal[2]5  + normal[3]4 + tumor[2]3  +tumor [3] 2 + var1 
                drug_el = key[4]+el
                #print (el)
                if key[4] in summary or 'R' in el:
                    if key[-17] == '*':
                        print(key[4])
                        drug_elx = drug_el.replace('(','*(') 
                        vgline = key[-11] + '\t' + key[-14] + '\t' + key[9] + '\t' + drug_elx   + '\t' + key[10] + '\t' + '-' + '\t' + 'MyCancerGenome' + '\n'
                        drugvar_info.write(vgline)
                        #'基因名' + '\t'  + '突变'  + '\t' +  '在病人所患癌症中可能获益的药物'  + '\t' +  '在其它癌症中可能获益的药物'  + '\t' + '抵抗的药物' + '\n'
                        if key[-1] not in dic_ondrug:
                            drug_elx = drug_el.replace('(','*(')
                            sum = [key[-11], '\t' ,key[-1],';',key[-9],';',key[-10],'\t',drug_elx,'\t' ,'-','\t','-','\n']    
                            dic_ondrug[key[-1]] = sum
                        else:
                            if key[4] not in dic_ondrug[key[-1]][8]:
                                drug_elx = drug_el.replace('(','*(')
                                dic_ondrug[key[-1]][8] = dic_ondrug[key[-1]][8] + ';' + drug_elx
                            else:
                                drugx = key[4] + '*'
                                if drugx in dic_ondrug[key[-1]][8]:
                                    sum_drug = dic_ondrug[key[-1]][8].split(';')
                                    for sd in sum_drug:
                                        if key[4] in sd:
                                            sd1= sd.replace('*','')
                                            drug3_el = compare_el(el,sd1,key[4])
                                            drug_elx = drug3_el.replace('(','*(')
                                            dic_ondrug[key[-1]][8] = dic_ondrug[key[-1]][8].replace(sd,drug_elx)

            if value == 'civic':
                if key[6] in black_list:
                    continue
                drugs = key[6].split(',') 
                #print (drugs[0],key[10],key[11],'civic','1')
                el = check_el(drugs[0],key[10],key[11],'civic','1')
                link = el+';'
                drug_el = key[6].replace(',',link)
                drug_el = drug_el + el
                if key[-17] == '*':
                    drug2_elx = drug_el.replace('(','*(')
                    vgline = key[-11] + '\t' + key[-14] + '\t' + key[3] + '\t' + drug2_elx   + '\t'  + dic_civicline2evs[key[12]] + '\t' + key[14] + ':' + key[13] + '\t' + 'CIViC' + '\n'
                    drugvar_info.write(vgline)
                    for drug in drugs:
                        if drug in summary or 'R' in el:
                            drug1_el = drug + el
                            
                            if key[-1] not in dic_ondrug:
                                drug1_elx = drug1_el.replace('(','*(')
                                sum = [key[-11], '\t' ,key[-1],';',key[-9],';',key[-10],'\t',drug1_elx,'\t' ,'-','\t','-','\n']
                                dic_ondrug[key[-1]] = sum
                            else:
                                if drug not in dic_ondrug[key[-1]][8]:
                                    drug1_elx = drug1_el.replace('(','*(')
                                    dic_ondrug[key[-1]][8] = dic_ondrug[key[-1]][8] + ';' + drug1_elx
                                else:
                                    drugx = drug + '*'
                                    if drugx in dic_ondrug[key[-1]][8]:
                                        #sum_drug = dic_offdrug[key[-1]][10].replace('*','')
                                        sum_drug = dic_ondrug[key[-1]][8].split(';')
                                        for sd in sum_drug:
                                            if drug in sd:
                                                sd1 = sd.replace('*','')
                                                drug3_el = compare_el(el,sd1,drug)
                                                drug1_elx = drug3_el.replace('(','*(')
                                                dic_ondrug[key[-1]][8] = dic_ondrug[key[-1]][8].replace(sd,drug1_elx)
                    

            if value == 'oncokb':
                if key[8] in black_list:
                    continue
                drugs = key[8].split(', ')
                el = check_el(drugs[0],key[7],'','oncokb','1')
                link = el+';'
                drug_el = key[8].replace(', ',link)
                drug_el = drug_el + el

                if key[-17] == '*':
                    drug2_elx = drug_el.replace('(','*(')
                    vgline = key[-11] + '\t' + key[-14] + '\t' + key[6] + '\t' + drug2_elx   + '\t'  + '-' + '\t' + 'PubMed:' + key[9] + '\t' + 'OncoKB' + '\n'
                    drugvar_info.write(vgline)
                    for drug in drugs:
                        if drug in summary or 'R' in el:
                            drug1_el = drug + el
                            
                            if key[-1] not in dic_ondrug:
                                drug1_elx = drug1_el.replace('(','*(')
                                sum = [key[-11], '\t' ,key[-1],';',key[-9],';',key[-10],'\t',drug1_elx,'\t' ,'-','\t','-','\n']
                                dic_ondrug[key[-1]] = sum
                            else:
                                if drug not in dic_ondrug[key[-1]][8]:
                                    drug1_elx = drug1_el.replace('(','*(')
                                    dic_ondrug[key[-1]][8] = dic_ondrug[key[-1]][8] + ';' + drug1_elx
                                else:
                                    drugx = drug + '*'
                                    if drugx in dic_ondrug[key[-1]][8]:
                                        #sum_drug = dic_offdrug[key[-1]][10].replace('*','')
                                        sum_drug = dic_ondrug[key[-1]][8].split(';')
                                        for sd in sum_drug:
                                            if drug in sd:
                                                sd1 = sd.replace('*','')
                                                drug3_el = compare_el(el,sd1,drug)
                                                drug1_elx = drug3_el.replace('(','*(')
                                                dic_ondrug[key[-1]][8] = dic_offdrug[key[-1]][8].replace(sd,drug1_elx)

        for key,value in dic_offlabel.items():
            key = key.replace("\n","")
            key = key.split('\t')
            if value == 'mcg':
                #print (key[7])
                if key[4] in black_list:
                    continue
                el = check_el(key[4],key[6],key[7],'mcg','2')
                #gene11 + calter10 + alte9r +str(score)8 +  funcgene7 + popfre6  + normal[2]5  + normal[3]4 + tumor[2]3  +tumor [3] 2 + var1 
                drug_el = key[4] + el
                if key[4] in summary  and key[4] not in onlabeldrug:  
                    if key[-17] == '*':
                        drug_elx = drug_el.replace('(','*(') 
                        vgline = key[-11] + '\t' + key[-14] + '\t' + key[9] + '\t' + drug_elx  + '\t' + key[10] + '\t' + '-' + '\t' + 'MyCancerGenome' + '\n'
                        drugvar_info.write(vgline)
                        if key[-1] not in dic_offdrug:
                            print (key[4])
                            drug_elx = drug_el.replace('(','*(')  
                            sum = [key[-11], '\t' ,key[-1],';',key[-9],';',key[-10],'\t','-','\t' ,drug_elx,'\t','-','\n']
                            dic_offdrug[key[-1]] = sum
                        else:
                            if key[4] not in dic_offdrug[key[-1]][10]:
                                drug_elx = drug_el.replace('(','*(')
                                dic_offdrug[key[-1]][10] = dic_offdrug[key[-1]][10] + ';' + drug_elx
                            else:
                                drugx = key[4] + '*'
                                if drugx in dic_offdrug[key[-1]][10]:
                                    sum_drug = dic_offdrug[key[-1]][10].split(';')
                                    for sd in sum_drug:
                                        if key[4] in sd:
                                            sd1= sd.replace('*','')
                                            drug3_el = compare_el(el,sd1,key[4])
                                            drug_elx = drug3_el.replace('(','*(')
                                            dic_offdrug[key[-1]][10] = dic_offdrug[key[-1]][10].replace(sd,drug_elx)

            if value == 'civic':
                if key[6] in black_list:
                    continue
                drugs = key[6].split(',')
                el = check_el(drugs[0],key[10],key[11],'civic','2')
                link = el+';'
                drug_el = key[6].replace(',',link)
                drug_el = drug_el + el
            
                if key[-17] == '*':
                    drug2_elx = drug_el.replace('(','*(')
                    vgline = key[-11] + '\t' + key[-14] + '\t' + key[3] + '\t' + drug2_elx   + '\t'  + dic_civicline2evs[key[12]] + '\t' + key[14] + ':' + key[13] + '\t' + 'CIViC' + '\n'
                    drugvar_info.write(vgline)
                    for drug in drugs:
                        if drug in summary  and drug not in onlabeldrug:
                            drug1_el = drug + el
                            print (drug)
                            if key[-1] not in dic_offdrug:
                                drug1_elx = drug1_el.replace('(','*(')  
                                sum = [key[-11], '\t' ,key[-1],';',key[-9],';',key[-10],'\t','-','\t' ,drug1_elx,'\t','-','\n']
                                dic_offdrug[key[-1]] = sum
                            else:
                                if drug not in dic_offdrug[key[-1]][10]:
                                    drug1_elx = drug1_el.replace('(','*(')
                                    dic_offdrug[key[-1]][10] = dic_offdrug[key[-1]][10] + ';' + drug1_elx
                                else:
                                    drugx = drug + '*'
                                    if drugx in dic_offdrug[key[-1]][10]:
                                        #sum_drug = dic_offdrug[key[-1]][10].replace('*','')
                                        sum_drug = dic_offdrug[key[-1]][10].split(';')
                                        for sd in sum_drug:
                                            if drug in sd:
                                                sd1 = sd.replace('*','')
                                                drug3_el = compare_el(el,sd1,drug)
                                                drug1_elx = drug3_el.replace('(','*(')
                                                dic_offdrug[key[-1]][10] = dic_offdrug[key[-1]][10].replace(sd,drug1_elx)
                    

            if value == 'oncokb':
                if key[8] in black_list:
                    continue
                drugs = key[8].split(', ')
                if 'R' in key[7]:
                    continue
                el = check_el(drugs[0],key[7],'','oncokb','2')
                print(drugs[0],key[7],'','2',el,'oncokb')
                print ('zzz')
                link = el+';'
                drug_el = key[8].replace(', ',link)
                drug_el = drug_el + el
            
                if key[-17] == '*':
                    drug2_elx = drug_el.replace('(','*(')
                    vgline = key[-11] + '\t' + key[-14] + '\t' + key[6] + '\t' + drug2_elx   + '\t'  + '-' + '\t' + 'PubMed:' + key[9] + '\t' + 'OncoKB' + '\n'
                    drugvar_info.write(vgline)
                    for drug in drugs:
                        if drug in summary  and drug not in onlabeldrug:
                            drug1_el = drug + el
                            print (drug)
                            if key[-1] not in dic_offdrug:
                                drug1_elx = drug1_el.replace('(','*(')  
                                sum = [key[-11], '\t' ,key[-1],';',key[-9],';',key[-10],'\t','-','\t' ,drug1_elx,'\t','-','\n']
                                dic_offdrug[key[-1]] = sum
                            else:
                                if drug not in dic_offdrug[key[-1]][10]:
                                    drug1_elx = drug1_el.replace('(','*(')
                                    dic_offdrug[key[-1]][10] = dic_offdrug[key[-1]][10] + ';' + drug1_elx
                                else:
                                    drugx = drug + '*'
                                    if drugx in dic_offdrug[key[-1]][10]:
                                        #sum_drug = dic_offdrug[key[-1]][10].replace('*','')
                                        sum_drug = dic_offdrug[key[-1]][10].split(';')
                                        for sd in sum_drug:
                                            if drug in sd:
                                                sd1 = sd.replace('*','')
                                                drug3_el = compare_el(el,sd1,drug)
                                                drug1_elx = drug3_el.replace('(','*(')
                                                dic_offdrug[key[-1]][10] = dic_offdrug[key[-1]][10].replace(sd,drug1_elx)




        


        for i in drug_array:
            if i in dic_metadrug:
                id_metadrug = dic_metadrug[i]
            else:
                continue
            if id_metadrug in dic_druginfo:
                lines = dic_druginfo[id_metadrug]
                #print(lines)
                lines = lines.replace('\n','')
                lines = lines.split('\t')
                if lines[6] and lines[7]:
                    line = lines[2] + '\t' + lines[6] + '\t'  + lines[-2] + ';'  +lines[-1] + '\t' + lines[7] + '\n'
                    #print(line)
                else:
                    if lines[7]:
                        line = lines[2] + '\t' + '-'+ '\t'  + lines[-2] + ';'  +lines[-1] + '\t' + lines[7] + '\n'
                    else:
                        if lines[6]:
                            line = lines[2] + '\t' + lines[6] + '\t'  + lines[-2] + ';'  +lines[-1] + '\t' + '-' + '\n'
                        else:
                            line = lines[2] + '\t' + '-' + '\t'  + lines[-2] + ';'  +lines[-1] + '\t' + '-' + '\n'
            #print(lines[0],lines[1],lines[2],lines[3],lines[4],lines[5],lines[6],lines[7])
                #print(line)
                dr.write(line)
        #something wrong
        for i in dic_ondrug:
            print (dic_ondrug[i])
            if i in dic_offdrug:
                #print (i,'z')
                print (dic_offdrug[i])
                dic_ondrug[i][10] = dic_offdrug[i][10]
                del dic_offdrug[i]
            re = []
            if '(R'in dic_ondrug[i][8]:
                
                drugs = dic_ondrug[i][8].split(';')
                se = []
                for j in drugs:
                    if '(R' in j:
                        re.append(j)
                    else:
                        se.append(j)
                sedrug = ';'.join(se)
                dic_ondrug[i][8] = sedrug
            else:
                
                dic_ondrug[i][-2] = '-'
            redrug = ';'.join(re)
            if re:
                dic_ondrug[i][-2] = redrug

            dic_ondrug[i] = sorting_drug(dic_ondrug[i])
            if dic_ondrug[i][8] == '':
                dic_ondrug[i][8] = '-'
            line = "".join(dic_ondrug[i])

            if line != '':
                f.write(line)

        for i in dic_offdrug:
            #print (dic_offdrug[i])
            #print('z')
            dic_offdrug[i][-1] = '-'
            dic_offdrug[i] = sorting_drug(dic_offdrug[i])
            
            line = "".join(dic_offdrug[i])

            if line != '':
                line = line
                f.write(line)

        f.close()
        dr.close()
        var_seq.close()
        drugvar_info.close()
for he in drugofall:
    print (he)

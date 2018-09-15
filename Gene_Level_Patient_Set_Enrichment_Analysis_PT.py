#Gene Set Enrichment Analysis #is in python3
from __future__ import division
import csv
import os, string, re, numpy as np, statsmodels.stats.multitest, scipy as sp, scipy.stats as stats, pandas as pd #sudo pip install scipy
import sys, logging, json
from functools import reduce
from multiprocessing import pool
from math import ceil
import requests
from numpy import in1d
from pandas import read_table, DataFrame
from gseapy.utils import unique, DEFAULT_LIBRARY
from requests.packages.urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
from scipy import interpolate
import timeit
tic=timeit.default_timer()

rootFolder = os.getcwd() #python script needs to be in the same directory (the current one) as the files of interest
overall_unique_variants = [] #to store the genes that are a variant in at least one patient
matrix_index = open("Percentage_Tumour.csv", "r")
matrix_index = matrix_index.readlines()
for root, dirs, files in os.walk(rootFolder):  
	for name in files: 
		patient_gene_variants = []   #iniate an empty list to store the patient specific unique gene
		if name.endswith(".csv"):  #only the patient data files
			if name != "Percentage_Tumour.csv":
				patient_variants = open(name, "r") #read each patient file
				for line in patient_variants:
					line = line.split(',') #split each line by a comma
					patient_gene_variants.append(line[0]) #get the gene name
				patient_gene_variants = sorted(list(set(patient_gene_variants))) 
				if "Gene" in patient_gene_variants: #remove the word Gene which is the column name
					index_gene = patient_gene_variants.index("Gene")
					del patient_gene_variants[index_gene]
				overall_unique_variants.extend(patient_gene_variants) #1574
overall_unique_variants = sorted(list(set(overall_unique_variants))) #512 variant genes

patient_counts = []
for root, dirs, files in os.walk(rootFolder):  
	for name in files: 
		patient_gene_variants = []   #iniate an empty list to store the patient specific unique genes
		patient_ind_counts = []
		if name.endswith(".csv"):  #only the patient data files
			if name != "Percentage_Tumour.csv":
				patient_variants = open(name, "r") #read each patient file
				for line in patient_variants:
					line = line.split(',') #split each line by a comma
					patient_gene_variants.append(line[0]) #get the gene name
				patient_gene_variants = sorted(list(set(patient_gene_variants)))
				if "Gene" in patient_gene_variants: #remove the word Gene which is the column name
					index_gene = patient_gene_variants.index("Gene")
					del patient_gene_variants[index_gene] 
				for gene in overall_unique_variants:
					if gene in patient_gene_variants:
						patient_ind_counts.append(1)
					else:
						patient_ind_counts.append(0)
				patient = re.sub('\.csv$', '', name) #the name of the file is the patients ID
				patient_ind_counts.insert(0, patient)
				patient_counts.append(patient_ind_counts)
overall_unique_variants.insert(0, "PATIENT")
overall_unique_variants.insert(1, "PERCENTAGE_TUMOUR")
matrix_indexes = []


for count_data in patient_counts:
	count_data = list(np.array(count_data).reshape(-1,))
	count_ID = count_data[0]
	for line in matrix_index:
		line = line.split(',')
		patient_ID = line[0]
		index_value = line[1]
		if patient_ID == count_ID:
			matrix_indexes.append(index_value)

patient_counts = np.array(patient_counts)
patient_counts = np.insert(patient_counts, 1, matrix_indexes, axis=1) #add matrix index values
patient_counts = patient_counts.tolist()
patient_counts.sort(key=lambda ele: float(ele[1]))#s: s[1]) #order according to the matrix index
patient_counts = np.array(patient_counts)
patient_counts = patient_counts.astype('S140')
patient_counts = np.insert(patient_counts, 0, overall_unique_variants, axis=0)
patient_counts_2 = np.matrix(patient_counts)
with open("PT_genes_ordered.csv", "a", newline='') as my_file:
	wr = csv.writer(my_file, dialect='excel')
	for row in range(0, len(patient_counts_2)):
		row = patient_counts_2[row]
		row = list(np.array(row).reshape(-1,))
		row = [el.decode('UTF-8') for el in row] 
		wr = csv.writer(my_file, dialect='excel')
		wr.writerow(row)
patient_counts = np.delete(patient_counts, 1, axis=1)

#gene list:
gene_list = list(patient_counts[:,0])
gene_list.pop(0)
gene_list = [el.decode('UTF-8') for el in gene_list] 

#correl_vector
correlation_vector = np.repeat(1, len(gene_list))

#gene_set
gene_Set_matrix = []
overall_unique_variants.pop(0)
overall_unique_variants.pop(0)
patient_counts = np.delete(patient_counts, 0, 1)
for column in patient_counts.T:
	gene_matrix = []
	gene_matrix.append(column[0])
	gene_matrix.append('description')
	for i in range(1, len(column)):
		pres_abs = column[i]
		if int(pres_abs) == 1:
			gene_matrix.append(gene_list[i-1]) #add to line in file then a tab space
	gene_Set_matrix.append(gene_matrix) #newline
with open('genes.gmt', 'w') as f:
	for line in gene_Set_matrix:
		line = list(np.array(line).reshape(-1,))
		#line = [el.decode('UTF-8') for el in line]
		line = '	'.join(line)
		f.write(line + "\n")
with open('genes.gmt') as gmt:
    gmt.read()

#gene_list = list(patient_counts[:,0])
#gene_list.pop(0)
#gene_list.pop(0)
#gene_list.pop(0)
#overall_unique_variants.pop(0)
#description = ['na']*len(overall_unique_variants)
#gene_set = np.column_stack((gene_list, description, overall_unique_variants))
#gene_set = np.column_stack((overall_unique_variants, description))
#np.savetxt("Genes.gmt", gene_set, delimiter="	", newline = "\n", fmt="%s")
#gene_set = np.column_stack((gene_Set_matrix.T[0], gene_Set_matrix.T[1], gene_Set_matrix.T[3]))
#gene_Set_matrix_decod = []
#for line in gene_Set_matrix:
#	line = list(np.array(line).reshape(-1,))
#	line = [el.decode('UTF-8') for el in line] 
#	gene_Set_matrix_decod.append(line)
#gene_Set_matrix_decod = np.array(gene_Set_matrix_decod)
#gene_Set_matrix = np.asmatrix(gene_Set_matrix)
#np.savetxt("Genes.gmt", gene_Set_matrix, delimiter="", newline = "\n", fmt="%s")

def gsea_gmt_parser(gmt, min_size = 1, max_size = 1000, gene_list=None):
    """Parse gene_sets.gmt(gene set database) file or download from enrichr server.
    :param gmt: the gene_sets.gmt file of GSEA input or an enrichr library name.
                checkout full enrichr library name here: http://amp.pharm.mssm.edu/Enrichr/#stats
    :param min_size: Minimum allowed number of genes from gene set also the data set. Default: 3.
    :param max_size: Maximum allowed number of genes from gene set also the data set. Default: 5000.
    :param gene_list: Used for filtering gene set. Only used this argument for :func:`call` method.
    :return: Return a new filtered gene set database dictionary.
    **DO NOT** filter gene sets, when use :func:`replot`. Because ``GSEA`` Desktop have already
    do this for you.
    """
    if gmt.lower().endswith(".gmt"):
        logging.info("User Defined gene sets is given.......continue..........")
        with open(gmt) as genesets:
             genesets_dict = { line.strip().split("\t")[0]: line.strip().split("\t")[2:]
                              for line in genesets.readlines()}
    # filtering dict
    if sys.version_info[0] >= 3 :
        genesets_filter =  {k: v for k, v in genesets_dict.items() if len(v) >= min_size and len(v) <= max_size}
    elif sys.version_info[0] == 2:
        genesets_filter =  {k: v for k, v in genesets_dict.iteritems() if len(v) >= min_size and len(v) <= max_size}
    else:
        logging.error("System failure. Please Provide correct input files")
        sys.exit(1)
    if gene_list is not None:
        subsets = sorted(genesets_filter.keys())
        for subset in subsets:
            tag_indicator = in1d(gene_list, genesets_filter.get(subset), assume_unique=True)
            tag_len = sum(tag_indicator)
            if tag_len <= min_size or tag_len >= max_size:
                del genesets_filter[subset]
            else:
                continue
    # some_dict = {key: value for key, value in some_dict.items() if value != value_to_remove}
    # use np.intersect1d() may be faster???
    filsets_num = len(genesets_dict) - len(genesets_filter)
    logging.info("%04d gene_sets have been filtered out when max_size=%s and min_size=%s"%(filsets_num, max_size, min_size))

    if filsets_num == len(genesets_dict):
        logging.error("No gene sets passed throught filtering condition!!!, try new paramters again!\n" +\
                         "Note: Gene names for gseapy is case sensitive." )
        sys.exit(1)
    else:
    	return genesets_filter

results = gsea_gmt_parser(gmt='genes.gmt', min_size = 2, max_size = 28, gene_list=None) #only includes those with at least two at at most 28 patients as those isoforms in a single patient are not informative nor are those in all the patients

def enrichment_score(gene_list, correl_vector, gene_set, weighted_score_type=0, 
                     nperm=30000, rs=np.random.RandomState(), single=False, scale=False): #nperm equal to permutation number: 1000?
    """This is the most important function of GSEApy. It has the same algorithm with GSEA and ssGSEA.
    :param gene_list:       The ordered gene list gene_name_list, rank_metric.index.values
    :param gene_set:        gene_sets in gmt file, please use gsea_gmt_parser to get gene_set.
    :param weighted_score_type:  It's same with gsea's weighted_score method. weighting by the correlation
                            is a very reasonable choice that allows significant gene sets with less than perfect coherence.
                            options: 0(classic),1,1.5,2. default:1. if one is interested in penalizing sets for lack of
                            coherence or to discover sets with any type of nonrandom distribution of tags, a value p < 1
                            might be appropriate. On the other hand, if one uses sets with large number of genes and only
                            a small subset of those is expected to be coherent, then one could consider using p > 1.
                            Our recommendation is to use p = 1 and use other settings only if you are very experienced
                            with the method and its behavior.
    :param correl_vector:   A vector with the correlations (e.g. signal to noise scores) corresponding to the genes in
                            the gene list. Or rankings, rank_metric.values #how reliable is each bit of information
    :param nperm:           Only used this parameter when computing esnull for statistical testing. set the esnull value
                            equal to the permutation number.
    :param rs:              Random state for initialize gene list shuffling. Default: np.random.RandomState(seed=None)
    :return:
     ES: Enrichment score (real number between -1 and +1)
     ESNULL: Enrichment score calculated from random permutation.
     Hits_Indices: index of a gene in gene_list, if gene included in gene_set.
     RES: Numerical vector containing the running enrichment score for all locations in the gene list .
    """
    N = len(gene_list)
    # Test whether each element of a 1-D array is also present in a second array
    # It's more intuitived here than orginal enrichment_score source code.
    # use .astype to covert bool to intergers
    #
    #for gene_patients in gene_set:
    #    print(gene_patients)
    #for item in gene_set.items():
    #    print(type(item))
    #k = gene, v = set of patients
    dict_values = list(gene_set.values()) #Modified: as np.in1d can't serach a dictionary
    tag_indicator = np.in1d(gene_list, dict_values, assume_unique=True).astype(int)  # notice that the sign is 0 (no tag) or 1 (tag), tests  if each patient is in the dictionary
    if weighted_score_type == 0:
        correl_vector = np.repeat(1, N)
    else:
        correl_vector = np.abs(correl_vector)**weighted_score_type
    # get indices of tag_indicator
    hit_ind = np.flatnonzero(tag_indicator).tolist() #list of the indexes of all the values of the array that are not zero
    # if used for compute esnull, set esnull equal to permutation number, e.g. 1000
    # else just compute enrichment scores
    # set axis to 1, because we have 2 dimentional array
    axis = 1
    tag_indicator = np.tile(tag_indicator, (nperm+1,1)) #repeat tag indicator nperm
    correl_vector = np.tile(correl_vector,(nperm+1,1))
    # gene list permutation
    for i in range(nperm): rs.shuffle(tag_indicator[i])
    # np.apply_along_axis(rs.shuffle, 1, tag_indicator)
    Nhint = tag_indicator.sum(axis=axis, keepdims=True)
    sum_correl_tag = np.sum(correl_vector*tag_indicator, axis=axis, keepdims=True)
    # compute ES score, the code below is identical to gsea enrichment_score method.
    norm_tag =  1.0/sum_correl_tag
    no_tag_indicator = 1 - tag_indicator #problem is here when tag indicator is a set of 1's i.e. all present 1- all 1's = 0
    Nmiss =  N - Nhint
    norm_no_tag = 1.0/Nmiss
    RES = np.cumsum(tag_indicator * correl_vector * norm_tag - no_tag_indicator * norm_no_tag, axis=axis)
    if scale: RES = RES / N
    if single:
        es_vec = RES.sum(axis=axis)
    else:
        max_ES, min_ES =  RES.max(axis=axis), RES.min(axis=axis)
        es_vec = np.where(np.abs(max_ES) > np.abs(min_ES), max_ES, min_ES)
    # extract values
    es, esnull, RES = es_vec[-1], es_vec[:-1], RES[-1,:]
    return es, esnull, hit_ind, RES

es_values = []
es_null_values = []
hit_ind_values = []
Res_values = []
geneys = []
for k, v in results.items(): #510 genes are 1 < patients < 29
    my_dictt = {k : list(v)}
    geneys.append(k)
    es, esnull, hit_ind, RES = enrichment_score(gene_list=gene_list, correl_vector=correlation_vector, gene_set=my_dictt, weighted_score_type=1, nperm=70000, rs=np.random.RandomState(), single=False, scale=False)
    es_values.append(es)
    es_null_values.append(esnull)
    hit_ind_values.append(hit_ind)
    Res_values.append(RES)

def gsea_pval(es, esnull):
    """Compute nominal p-value.
    From article (PNAS):
    estimate nominal p-value for S from esnull by using the positive
    or negative portion of the distribution corresponding to the sign
    of the observed ES(S).
    """
    # to speed up, using numpy function to compute pval in parallel.
    es = np.array(es)
    esnull = np.array(esnull)
    # try:
    condlist = [ es < 0, es >=0]
    choicelist = [np.sum(esnull < es.reshape(len(es),1), axis=1)/ np.sum(esnull < 0, axis=1),
                  np.sum(esnull >= es.reshape(len(es),1), axis=1)/ np.sum(esnull >= 0, axis=1)]
    pval = np.select(condlist, choicelist)
    return pval
    # except:
    #    return np.repeat(1.0 ,len(es))

def gsea_significance(enrichment_scores, enrichment_nulls):
    """Compute nominal pvals, normalized ES, and FDR q value.
        For a given NES(S) = NES* >= 0. The FDR is the ratio of the percentage of all (S,pi) with
        NES(S,pi) >= 0, whose NES(S,pi) >= NES*, divided by the percentage of
        observed S wih NES(S) >= 0, whose NES(S) >= NES*, and similarly if NES(S) = NES* <= 0.
    """
    # For a zero by zero division (undetermined, results in a NaN),
    # np.seterr(divide='ignore', invalid='ignore')
    import warnings
    warnings.simplefilter("ignore")
    logging.debug("Start to compute pvals..................................")
    # compute pvals.
    enrichmentPVals = gsea_pval(enrichment_scores, enrichment_nulls).tolist()
    # new normalize enrichment score calculating method. this could speed up significantly.
    esnull_meanPos = []
    esnull_meanNeg = []
    es = np.array(enrichment_scores)
    esnull = np.array(enrichment_nulls)
    for i in range(len(enrichment_scores)):
        enrNull = esnull[i]
        meanPos = enrNull[enrNull >= 0].mean()
        esnull_meanPos.append(meanPos)
        meanNeg = enrNull[enrNull < 0 ].mean()
        esnull_meanNeg.append(meanNeg)
    pos = np.array(esnull_meanPos).reshape(len(es), 1)
    neg = np.array(esnull_meanNeg).reshape(len(es), 1)
    # compute normalized enrichment score and normalized esnull
    logging.debug("Compute normalized enrichment score and normalized esnull")
    try:
        condlist1 = [ es >= 0, es < 0]
        choicelist1 = [ es/esnull_meanPos, -es/esnull_meanNeg ]
        nEnrichmentScores = np.select(condlist1, choicelist1).tolist()
        condlist2 = [ esnull >= 0, esnull < 0]
        choicelist2 = [ esnull/pos, -esnull/neg ]
        nEnrichmentNulls = np.select(condlist2, choicelist2)

    except:  #return if according nes, nesnull is uncalculable
        nEnrichmentScores = np.repeat(0.0, es.size).tolist()
        nEnrichmentNulls = np.repeat(0.0 , es.size).reshape(esnull.shape)
    logging.debug("start to compute fdrs..................................")
    # FDR null distribution histogram
    # create a histogram of all NES(S,pi) over all S and pi
    # Use this null distribution to compute an FDR q value,
    # vals = reduce(lambda x,y: x+y, nEnrichmentNulls, [])
    # nvals = np.array(sorted(vals))
    # or
    nvals = np.sort(nEnrichmentNulls.flatten())
    nnes = np.array(sorted(nEnrichmentScores))
    fdrs = []
    # FDR computation
    for i in range(len(enrichment_scores)):
        nes = nEnrichmentScores[i]
        if nes >= 0:
            allPos = int(len(nvals) - np.searchsorted(nvals, 0, side="left"))
            allHigherAndPos = int(len(nvals) - np.searchsorted(nvals, nes, side="left"))
            nesPos = len(nnes) - int(np.searchsorted(nnes, 0, side="left"))
            nesHigherAndPos = len(nnes) - int(np.searchsorted(nnes, nes, side="left"))
        else:
            allPos = int(np.searchsorted(nvals, 0, side="left"))
            allHigherAndPos = int(np.searchsorted(nvals, nes, side="right"))
            nesPos = int(np.searchsorted(nnes, 0, side="left"))
            nesHigherAndPos = int(np.searchsorted(nnes, nes, side="right"))
        try:
            pi_norm = allHigherAndPos/float(allPos)
            pi_obs = nesHigherAndPos/float(nesPos)

            fdr = pi_norm/pi_obs if pi_norm/pi_obs < 1.0  else 1.0
            fdrs.append(fdr)
        except:
            fdrs.append(1000000000.0)

    logging.debug("Statistical testing finished.............................")
    return enrichment_scores, nEnrichmentScores, enrichmentPVals, fdrs

enrichment_scores, nEnrichmentScores, enrichmentPVals, fdrs = gsea_significance(es_values, es_null_values)

def estimate(pv, m=None, verbose=False, lowmem=False, pi0=None): #from: https://github.com/nfusi/qvalue
    """
    Estimates q-values from p-values
    Args
    =====
    m: number of tests. If not specified m = pv.size
    verbose: print verbose messages? (default False)
    lowmem: use memory-efficient in-place algorithm
    pi0: if None, it's estimated as suggested in Storey and Tibshirani, 2003.
         For most GWAS this is not necessary, since pi0 is extremely likely to be
         1
    """
    assert(pv.min() >= 0 and pv.max() <= 1), "p-values should be between 0 and 1"
    original_shape = pv.shape
    pv = pv.ravel()  # flattens the array in place, more efficient than flatten()
    if m is None:
        m = float(len(pv))
    else:
        # the user has supplied an m
        m *= 1.0
    # if the number of hypotheses is small, just set pi0 to 1
    if len(pv) < 100 and pi0 is None:
        pi0 = 1.0
    elif pi0 is not None:
        pi0 = pi0
    else:
        # evaluate pi0 for different lambdas
        pi0 = []
        lam = sp.arange(0, 0.90, 0.01)
        counts = sp.array([(pv > i).sum() for i in sp.arange(0, 0.9, 0.01)])
        for l in range(len(lam)):
            pi0.append(counts[l]/(m*(1-lam[l])))
        pi0 = sp.array(pi0)
        # fit natural cubic spline
        tck = interpolate.splrep(lam, pi0, k=3)
        pi0 = interpolate.splev(lam[-1], tck)
        if verbose:
            print("qvalues pi0=%.3f, estimated proportion of null features " % pi0)

        if pi0 > 1:
            if verbose:
                print("got pi0 > 1 (%.3f) while estimating qvalues, setting it to 1" % pi0)
            pi0 = 1.0
    assert(pi0 >= 0 and pi0 <= 1), "pi0 is not between 0 and 1: %f" % pi0
    if lowmem:
        # low memory version, only uses 1 pv and 1 qv matrices
        qv = sp.zeros((len(pv),))
        last_pv = pv.argmax()
        qv[last_pv] = (pi0*pv[last_pv]*m)/float(m)
        pv[last_pv] = -sp.inf
        prev_qv = last_pv
        for i in range(int(len(pv))-2, -1, -1): #range in python 3 is xrange in python3
            cur_max = pv.argmax()
            qv_i = (pi0*m*pv[cur_max]/float(i+1))
            pv[cur_max] = -sp.inf
            qv_i1 = prev_qv
            qv[cur_max] = min(qv_i, qv_i1)
            prev_qv = qv[cur_max]
    else:
        p_ordered = sp.argsort(pv)
        pv = pv[p_ordered]
        qv = pi0 * m/len(pv) * pv
        qv[-1] = min(qv[-1], 1.0)
        for i in range(len(pv)-2, -1, -1):
            qv[i] = min(pi0*m*pv[i]/(i+1.0), qv[i+1])
        # reorder qvalues
        qv_temp = qv.copy()
        qv = sp.zeros_like(qv)
        qv[p_ordered] = qv_temp
    # reshape qvalues
    qv = qv.reshape(original_shape)
    return qv

enrichmentPVals_array = np.array(enrichmentPVals)
qv = estimate(enrichmentPVals_array, m=None, verbose=False, lowmem=False, pi0=None)

count = 0
for qv_val in qv:
    if qv_val < 0.05:
        count += 1
print("This patient set enrichment analysis has resulted in", count, "significant genes after correcting for multiple hypothesis testing.")

qv = qv.tolist()
genes = geneys
genes.insert(0, "Gene")
enrichment_scores.insert(0, "enrichment_score")
nEnrichmentScores.insert(0, "nEnrichmentScores")
enrichmentPVals.insert(0, "enrichmentPVals")
qv.insert(0, "qv")
my_results = np.column_stack((genes, enrichment_scores, nEnrichmentScores, enrichmentPVals, qv))
my_results = np.matrix(my_results)
np.savetxt("Enrichment_Analysis_Gene_Results_TP.csv", my_results, delimiter=",", newline = "\n", fmt="%s")

#Get the genes in every patient i.e. those that are likely cancer associated
with open('genes.gmt') as genesets:
    genesets_dict = { line.strip().split("\t")[0]: line.strip().split("\t")[2:] for line in genesets.readlines()}
    genesets_filter_whole =  {k: v for k, v in genesets_dict.items() if len(v) == 29}
    all_present_genes = genesets_filter_whole.keys()
    with open('Genes_Present_In_All_Samples_TP.csv', 'w') as filey:
        filey.write("Genes_In_All_Samples" + "\n")
        for key in all_present_genes:
            filey.write(key + "\n")
    genesets_filter_whole_wo_control =  {k: v for k, v in genesets_dict.items() if len(v) == 22}
    all_present_genes_ex_control = genesets_filter_whole_wo_control.keys()
    with open('Genes_Present_In_Cancer_Samples_TP.csv', 'w') as filey:
        filey.write("Genes_In_All_Cancer_Samples" + "\n") 
        for key in all_present_genes_ex_control:
            my_values = genesets_filter_whole_wo_control[key]
            if ("G94", "G93", "G88", "G72", "G58", "G10", "G42") not in my_values == True:
                filey.write(key + "\n")

header = ['Gene', 'Enrichment_Score','normalised_Enrichemnt_Score','p_value', 'q_value']
with open("Significant_gene_TP.csv", "a", newline='') as my_file:
	wr = csv.writer(my_file, dialect='excel')
	wr.writerow(header)
	for row in range(1, len(my_results)):
		row = my_results[row]
		row = list(np.array(row).reshape(-1,))
		if float(row[4]) < 0.05:
			wr = csv.writer(my_file, dialect='excel')
			wr.writerow(row)
  
toc=timeit.default_timer()
print(toc - tic)


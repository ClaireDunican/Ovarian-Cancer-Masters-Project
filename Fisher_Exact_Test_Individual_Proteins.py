#Make new versions of the files with unique genes!
import os, string, re, numpy as np, statsmodels.stats.multitest, scipy as sp, scipy.stats as stats, pandas as pd #sudo pip install scipy
from scipy import interpolate
rootFolder = os.getcwd() #python script needs to be in the same directory (the current one) as the files of interest
overall_unique_variants = [] #to store the genes that are a variant in at least one patient
pre_variant_genes = []
good_variant_genes = []
bad_variant_genes = []
control_variant_genes = []
metadata = open("Patient_Metadata.csv", "r")
metaline = metadata.readlines()
good_total = 0
bad_total = 0
pre_total = 0
control_total = 0
for root, dirs, files in os.walk(rootFolder):  
	for name in files: 
		patient_gene_variants = []   #iniate an empty list to store the patient specific unique genes
		if name.endswith(".csv"):  #only the patient data files
			if name != "Patient_Metadata.csv":
				patient_variants = open(name, "r") #read each patient file
				for line in patient_variants:
					line = line.split(',') #split each line by a comma
					patient_gene_variants.append(line[0]) #get the gene number
				patient_gene_variants = sorted(list(set(patient_gene_variants))) #converting to a set removes the non-unique genes then convert back to a list and sort alphabetically
				patient = re.sub('\.csv$', '', name) #the name of the file is the patients ID
				for meta_line in metaline: #get the matching patient line in the metadata file and get the patients state:
					if patient in meta_line:
						meta_line = meta_line.split(',')
						biopsy_type = meta_line[4]
				patient_gene_variants.append(biopsy_type)
				if "Gene" in patient_gene_variants: #remove the word Gene which is the column name
					index_gene = patient_gene_variants.index("Gene")
					del patient_gene_variants[index_gene]
				if "pre" in patient_gene_variants:
					pre_total += 1 
					index_pre = patient_gene_variants.index("pre")
					del patient_gene_variants[index_pre]
					pre_variant_genes.extend(patient_gene_variants)
				elif "good" in patient_gene_variants:
					good_total += 1
					index_good = patient_gene_variants.index("good")
					del patient_gene_variants[index_good]
					good_variant_genes.extend(patient_gene_variants)
				elif "bad" in patient_gene_variants:
					bad_total += 1
					index_bad = patient_gene_variants.index("bad")
					del patient_gene_variants[index_bad]
					bad_variant_genes.extend(patient_gene_variants)
				elif "control" in patient_gene_variants:
					control_total += 1
					index_control = patient_gene_variants.index("control")
					del patient_gene_variants[index_control]
					control_variant_genes.extend(patient_gene_variants) #add each list for each individual on to the end
				overall_unique_variants.extend(patient_gene_variants) #1574

overall_unique_variants = sorted(list(set(overall_unique_variants))) #512 variant genes
#print(len(overall_unique_variants))
good_total_counts = []
bad_total_counts = []
pre_total_counts = []
control_total_counts = []
pre_total_counts.append("Pre") #column names for the patient subsets
bad_total_counts.append("Bad")
good_total_counts.append("Good")
control_total_counts.append("Control")
for gene in overall_unique_variants: #count how many times each gene appears in each subset of patients
	control_total_counts.append(control_variant_genes.count(gene))
	good_total_counts.append(good_variant_genes.count(gene)) #counts the number of times the gene is in the list
	bad_total_counts.append(bad_variant_genes.count(gene))
	pre_total_counts.append(pre_variant_genes.count(gene))

overall_unique_variants.insert(0, "Gene") #column name for the first column
matrix_gene_counts = np.column_stack((overall_unique_variants, good_total_counts, bad_total_counts, pre_total_counts, control_total_counts))
matrix_gene_counts = np.matrix(matrix_gene_counts)

everything = matrix_gene_counts.tolist()
remove_rows = []
for sublisty in range(1, len(everything)):
	count = everything[sublisty]
	if int(count[1]) == 1:
		if int(count[2]) == 0:
			if int(count[4]) == 0:
				remove_rows.append(everything.index(count))
	elif int(count[2]) == 1:
		if int(count[1]) == 0:
			if int(count[4]) == 0:
				remove_rows.append(everything.index(count))
	elif int(count[4]) == 1:
		if int(count[1]) == 0:
			if int(count[2]) == 0:
				remove_rows.append(everything.index(count))
	elif int(count[3]) == 1:
		if int(count[1]) == 0:
			if int(count[2]) == 0:
				if int(count[4]) == 0:
					remove_rows.append(everything.index(count))

matrix_gene_counts = [i for j, i in enumerate(everything) if j not in remove_rows] #251 left
matrix_gene_counts = np.array(matrix_gene_counts)
matrix_gene_counts = np.matrix(matrix_gene_counts)

np.savetxt("Variant_Protein_Count_Table.csv", matrix_gene_counts, delimiter=",", newline = "\n", fmt="%s") #save to a file

#Fisher Exact Tests:

matrix_gene_counts = np.delete(matrix_gene_counts, (0), axis=0) #remove the rownames

p_values_good_bad = []
odds_good_bad = []
p_val_good = []
odds_good = []
p_val_bad = []
odds_bad = []
p_val_pre = []
odds_pre = []
genes_listy = []
for row in matrix_gene_counts:
	row = list(np.array(row).reshape(-1,))
	genes_listy.append(row[0])
	absent_good = good_total - int(row[1]) #total_good_indiduals - individuals with variant gene present
	absent_bad = bad_total - int(row[2]) 
	absent_pre = pre_total - int(row[3])
	absent_control = control_total - int(row[4])
	oddsratio, pvalue = stats.fisher_exact([[int(row[1]), int(row[2])], [absent_good, absent_bad]])#, alternative='less') #0 = gene name, 1=good, 2=bad, 3=pre-chemo, 4=control, , ‘less’, ‘greater’},
	p_values_good_bad.append(pvalue)
	odds_good_bad.append(oddsratio)
	oddsratio, pvalue = stats.fisher_exact([[int(row[1]), int(row[4])], [absent_good, absent_control]])#, alternative='less') #0 = gene name, 1=good, 2=bad, 3=pre-chemo, 4=control
	p_val_good.append(pvalue)
	odds_good.append(oddsratio)
	oddsratio, pvalue = stats.fisher_exact([[int(row[2]), int(row[4])], [absent_bad, absent_control]])#, alternative='less') #0 = gene name, 1=good, 2=bad, 3=pre-chemo, 4=control
	p_val_bad.append(pvalue)
	odds_bad.append(oddsratio)
	oddsratio, pvalue = stats.fisher_exact([[int(row[3]), int(row[4])], [absent_pre, absent_control]])#, alternative='greater') #0 = gene name, 1=good, 2=bad, 3=pre-chemo, 4=control
	p_val_pre.append(pvalue)
	odds_pre.append(oddsratio)

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
        for i in range(int(len(pv))-2, -1, -1):
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

good_p_val_q = np.array(p_val_good)
good_p_val_q = estimate(good_p_val_q, m=None, verbose=False, lowmem=False, pi0=None) #qv = qvalue.estimate(p_val_good)
Good_q_val = good_p_val_q.tolist()
Good_q_val.insert(0, "Good_q_value")

bad_p_val_q = np.array(p_val_bad)
bad_p_val_q = estimate(bad_p_val_q, m=None, verbose=False, lowmem=False, pi0=None) #qv = qvalue.estimate(p_val_bad)
Bad_q_val = bad_p_val_q.tolist()
Bad_q_val.insert(0, "Bad_q_value")

good_bad_p_val_q = np.array(p_values_good_bad)
good_bad_p_val_q = estimate(good_bad_p_val_q, m=None, verbose=False, lowmem=False, pi0=None) #qv = qvalue.estimate(p_val_good_bad)
Good_Bad_q_val = good_bad_p_val_q.tolist()
Good_Bad_q_val.insert(0, "Good_Bad_q_value")


p_val_pre_q = np.array(p_val_pre)
p_val_pre_q = estimate(p_val_pre_q, m=None, verbose=False, lowmem=False, pi0=None)
pre_q_val = p_val_pre_q.tolist()
pre_q_val.insert(0, "Pre_Control_q_value")

genes_listy.insert(0, "Gene")
odds_good.insert(0, "Good_Odds")
odds_pre.insert(0, "Pre_Odds")
p_val_pre.insert(0, "Pre_Control_p_value")
p_val_good.insert(0, "Good_p_value")
odds_bad.insert(0, "Bad_Odds")
p_val_bad.insert(0, "Bad_p_value")
odds_good_bad.insert(0, "Bad_vs_Good_Odds")
p_values_good_bad.insert(0, "Bad_vs_Good_p_value")		
matrix_gene_results_gene = np.column_stack((genes_listy, p_val_good, Good_q_val, p_val_bad, Bad_q_val, p_values_good_bad, Good_Bad_q_val, p_val_pre, pre_q_val, odds_good, odds_bad, odds_good_bad, odds_pre))
matrix_gene_results_gene = np.matrix(matrix_gene_results_gene)
np.savetxt("Gene_Results_Table.csv", matrix_gene_results_gene, delimiter=",", newline = "\n", fmt="%s") #save to a file




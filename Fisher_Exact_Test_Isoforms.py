import os, string, re, numpy as np, statsmodels.stats.multitest, scipy as sp, scipy.stats as stats, pandas as pd #sudo pip install scipy
from scipy import interpolate
rootFolder = os.getcwd() #python script needs to be in the same directory (the current one) as the files of interest
overall_unique_isoforms = [] #to store the genes that are a variant in at least one patient
pre_isoforms_genes = []
good_isoforms_genes = []
bad_isoforms_genes = []
control_isoforms_genes = []
metadata = open("Patient_Metadata.csv", "r")
metaline = metadata.readlines()
good_total = 0
bad_total = 0
pre_total = 0
control_total = 0
for root, dirs, files in os.walk(rootFolder):  
	for name in files: 
		patient_gene_variants = []   #iniate an empty list to store the patient specific unique gene
		if name.endswith(".csv"):  #only the patient data files
			if name != "Patient_Metadata.csv":
				patient_variants = open(name, "r") #read each file
				for line in patient_variants:
					line = line.split(',') #split each line by a comma
					patient_gene_variants.append("G: " + line[0] + ", V: " + line[2]) #get each isoforms name and gene
				patient_gene_variants = list(set(patient_gene_variants)) #removes the same isoform being included more than once
				patient = re.sub('\.csv$', '', name) #add the type of response to the end of the list for that patient
				for meta_line in metaline: #get the patients state
					if patient in meta_line:
						meta_line = meta_line.split(',')
						biopsy_type = meta_line[4]
				patient_gene_variants.append(biopsy_type)
				if "G: Gene, V: ID" in patient_gene_variants: #remove the word Gene which is the column name
					index_gene = patient_gene_variants.index("G: Gene, V: ID")
					del patient_gene_variants[index_gene]
				if "pre" in patient_gene_variants:
					pre_total += 1
					index_pre = patient_gene_variants.index("pre")
					del patient_gene_variants[index_pre]
					pre_isoforms_genes.extend(patient_gene_variants)
				elif "good" in patient_gene_variants:
					good_total += 1
					index_good = patient_gene_variants.index("good")
					del patient_gene_variants[index_good]
					good_isoforms_genes.extend(patient_gene_variants)
				elif "bad" in patient_gene_variants:
					bad_total += 1
					index_bad = patient_gene_variants.index("bad")
					del patient_gene_variants[index_bad]
					bad_isoforms_genes.extend(patient_gene_variants)
				elif "control" in patient_gene_variants:
					control_total += 1
					index_control = patient_gene_variants.index("control")
					del patient_gene_variants[index_control]
					control_isoforms_genes.extend(patient_gene_variants)
				overall_unique_isoforms.extend(patient_gene_variants)

overall_unique_isoforms = list(set(overall_unique_isoforms)) #list of unique isoforms: 8161
good_total_counts = []
bad_total_counts = []
pre_total_counts = []
control_total_counts = []
pre_total_counts.append("Pre") #column names the patient subset columns
bad_total_counts.append("Bad")
good_total_counts.append("Good")
control_total_counts.append("Control")
corresponding_genes = []
just_ID = []
for gene in overall_unique_isoforms: #for each unique gene: #check that the unique isoforms arn't just subnumbers and that they total properly
	control_total_counts.append(control_isoforms_genes.count(gene))
	good_total_counts.append(good_isoforms_genes.count(gene))
	bad_total_counts.append(bad_isoforms_genes.count(gene))
	pre_total_counts.append(pre_isoforms_genes.count(gene))
	gene = list(np.array(gene).reshape(-1,))
	gene = str(gene)
	ID = gene[(gene.find('V: ') + 3):-2]
	just_ID.append(ID)
	gene = gene[5:gene.find(',')]
	gene = str(gene)
	corresponding_genes.append(gene)

just_ID.insert(0, "Isoform_ID") #column name for the first column
corresponding_genes.insert(0, "Gene") #column name for the first column
matrix_gene_counts = np.column_stack((corresponding_genes, just_ID, good_total_counts, bad_total_counts, pre_total_counts, control_total_counts))
matrix_gene_counts = np.matrix(matrix_gene_counts)

everything = matrix_gene_counts.tolist()
remove_rows = []
for sublisty in range(1, len(everything)):
	count = everything[sublisty]
	if int(count[2]) == 1:
		if int(count[3]) == 0:
			if int(count[5]) == 0:
				#if int(count[4]) == 0:
				remove_rows.append(everything.index(count))
	elif int(count[3]) == 1:
		if int(count[2]) == 0:
			if int(count[5]) == 0:
				#if int(count[4]) == 0:
				remove_rows.append(everything.index(count))
	elif int(count[5]) == 1:
		if int(count[3]) == 0:
			if int(count[2]) == 0:
				#if int(count[4]) == 0:
				remove_rows.append(everything.index(count))
	elif int(count[4]) == 1:
		if int(count[3]) == 0:
			if int(count[2]) == 0:
				if int(count[5]) == 0:
					remove_rows.append(everything.index(count))


matrix_gene_counts = [i for j, i in enumerate(everything) if j not in remove_rows]
matrix_gene_counts = np.array(matrix_gene_counts)
matrix_gene_counts = np.matrix(matrix_gene_counts)

np.savetxt("Isoform_Count_Table.csv", matrix_gene_counts, delimiter=",", newline = "\n", fmt="%s") #save to a file

#Fisher Exact Tests:

matrix_gene_counts = np.delete(matrix_gene_counts, (0), axis=0) #remove the colnames

#get each of the unique genes variants:
p_val_iso_good = []
odds_iso_good = []
p_val_iso_bad = []
odds_iso_bad = []
p_val_iso_good_v_bad = []
odds_iso_good_bad = []
p_val_pre = []
odds_pre = []
gene_listy = []
isoform_listy = []

for row in matrix_gene_counts: #gets all the isoforms for that gene
	#matrix_gene_counts = matrix_gene_counts.astype(int)
	row = list(np.array(row).reshape(-1,)) #isoform of interest row
	gene_listy.append(row[0])
	isoform_listy.append(row[1])
	ref_good = good_total - int(row[2]) #total number of patients in group - total number of isoform of interest
	ref_bad = bad_total - int(row[3])
	ref_pre = pre_total - int(row[4]) #corresponding_genes, just_ID, good_total_counts, bad_total_counts, pre_total_counts, control_total_counts
	ref_control = control_total - int(row[5])
	oddsratio, pvalue = stats.fisher_exact([[ref_good, ref_control], [int(row[2]), int(row[5])]])#, alternative='less')
	p_val_iso_good.append(pvalue)
	odds_iso_good.append(oddsratio)
	oddsratio, pvalue = stats.fisher_exact([[ref_bad, ref_control], [int(row[3]), int(row[5])]])#, alternative='less')
	p_val_iso_bad.append(pvalue)
	odds_iso_bad.append(oddsratio)
	oddsratio, pvalue = stats.fisher_exact([[ref_bad, ref_good], [int(row[3]), int(row[2])]])#, alternative='less')
	p_val_iso_good_v_bad.append(pvalue)
	odds_iso_good_bad.append(oddsratio)
	oddsratio, pvalue = stats.fisher_exact([[int(row[4]), int(row[5])], [ref_pre, ref_control]], alternative='greater')
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
        for i in range(int(len(pv))-2, -1, -1): #for python 2 range is xrange
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

good_p_val_q = np.array(p_val_iso_good)
good_p_val_q = estimate(good_p_val_q, m=None, verbose=False, lowmem=False, pi0=None) #qv = qvalue.estimate(p_val_good)
Good_q_val = good_p_val_q.tolist()
Good_q_val.insert(0, "Good_q_value")

bad_p_val_q = np.array(p_val_iso_bad)
bad_p_val_q = estimate(bad_p_val_q, m=None, verbose=False, lowmem=False, pi0=None) #qv = qvalue.estimate(p_val_good)
Bad_q_val = bad_p_val_q.tolist()
Bad_q_val.insert(0, "Bad_q_value")

good_bad_p_val_q = np.array(p_val_iso_good_v_bad)
good_bad_p_val_q = estimate(good_bad_p_val_q, m=None, verbose=False, lowmem=False, pi0=None) #qv = qvalue.estimate(p_val_good)
Good_Bad_q_val = good_bad_p_val_q.tolist()
Good_Bad_q_val.insert(0, "Good_Bad_q_value")

p_val_pre_q = np.array(p_val_pre)
p_val_pre_q = estimate(p_val_pre_q, m=None, verbose=False, lowmem=False, pi0=None)
pre_q_val = p_val_pre_q.tolist()
pre_q_val.insert(0, "Pre_Control_q_value")

p_val_iso_good.insert(0, "Good_p_value")
odds_iso_good.insert(0, "Good_Odds")
p_val_pre.insert(0, "Pre_Control_p_value")
odds_pre.insert(0, "Pre_Odds")	
p_val_iso_bad.insert(0, "Bad_p_value")
odds_iso_bad.insert(0, "Bad_Odds")	
p_val_iso_good_v_bad.insert(0, "Good_v_bad_p_value")	
odds_iso_good_bad.insert(0, "Good_Bad_Odds")	
gene_listy.insert(0, "Gene")
isoform_listy.insert(0, "Isoform")
matrix_gene_results_iso = np.column_stack((gene_listy, isoform_listy, p_val_iso_good, Good_q_val, p_val_iso_bad, Bad_q_val, p_val_iso_good_v_bad, Good_Bad_q_val, p_val_pre, pre_q_val, odds_iso_good, odds_iso_bad, odds_iso_good_bad, odds_pre))
matrix_gene_results_iso = np.matrix(matrix_gene_results_iso)
np.savetxt("Results_Table_Isoforms.csv", matrix_gene_results_iso, delimiter=",", newline = "\n", fmt="%s") #save to a file



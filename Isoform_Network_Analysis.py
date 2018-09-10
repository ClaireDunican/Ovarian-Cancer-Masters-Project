#Gene Set Enrichment Analysis at the isoform level
import matplotlib.pyplot as plt; plt.rcdefaults()
import pandas as pd
import os, string, re, numpy as np, statsmodels.stats.multitest, scipy as sp, scipy.stats as stats, pandas as pd #sudo pip install scipy

rootFolder = os.getcwd() #python script needs to be in the same directory (the current one) as the files of interest
overall_unique_variants = [] #to store the genes that are a variant in at least one patient
matrix_index = open("Matrix_index_values.csv", "r")
matrix_index = matrix_index.readlines()
for root, dirs, files in os.walk(rootFolder):  
	for name in files: 
		patient_gene_variants = []   #iniate an empty list to store the patient specific unique genes
		if name.endswith(".csv"):  #only the patient data files
			if name != "Matrix_index_values.csv":
				patient_variants = open(name, "r") #read each patient file
				for line in patient_variants:
					line = line.split(',') #split each line by a comma
					patient_gene_variants.append(line[2]) #get the gene name
				patient_gene_variants = sorted(list(set(patient_gene_variants))) #converting to a set removes the non-unique genes then convert back to a list and sort alphabeticall
				if "ID" in patient_gene_variants: #remove the word Gene which is the column name
					index_gene = patient_gene_variants.index("ID")
					del patient_gene_variants[index_gene]
				overall_unique_variants.extend(patient_gene_variants)
overall_unique_variants = sorted(list(set(overall_unique_variants))) #512 variant genes

patient_counts = []
for root, dirs, files in os.walk(rootFolder):  
	for name in files: 
		patient_gene_variants = []   #iniate an empty list to store the patient specific unique genes
		patient_ind_counts = []
		if name.endswith(".csv"):  #only the patient data files
			if name != "Matrix_index_values.csv":
				patient_variants = open(name, "r") #read each patient file
				for line in patient_variants:
					line = line.split(',') #split each line by a comma
					patient_gene_variants.append(line[2]) #get the gene name
				patient_gene_variants = sorted(list(set(patient_gene_variants))) #converting to a set removes the non-unique genes then convert back to a list and sort alphabeticall
				if "ID" in patient_gene_variants: #remove the word Gene which is the column name
					index_gene = patient_gene_variants.index("ID")
					del patient_gene_variants[index_gene] 
				for gene in overall_unique_variants:
					if gene in patient_gene_variants:
						patient_ind_counts.append(1)
					else:
						patient_ind_counts.append(0)
				patient = re.sub('\.csv$', '', name) #the name of the file is the patients ID
				patient_ind_counts.insert(0, patient)
				patient_counts.append(patient_ind_counts)
overall_unique_variants.insert(0, "Patient")

patient_counts = np.array(patient_counts)
patient_counts = patient_counts.astype('S140')
patient_counts = np.insert(patient_counts, 0, overall_unique_variants, axis=0)
decoded_patient_counts = []
for line in patient_counts:
	line = list(np.array(line).reshape(-1,))
	line = [el.decode('UTF-8') for el in line] 
	decoded_patient_counts.append(line)
patient_counts = np.array(decoded_patient_counts)
np.savetxt("Isoforms_For_Pairwise_Comparisons.csv", patient_counts, delimiter=",", newline = "\n", fmt="%s") #save to a file

patient_counts_2 = np.delete(patient_counts, 0, axis=1)
patient_counts_2 = np.delete(patient_counts_2, 0, axis=0)
patient_counts_2 = patient_counts_2.astype(int)
totals = np.sum(patient_counts_2, axis=1)
totals_per_gene = np.sum(patient_counts_2, axis=0)
patient_counts_inv = patient_counts.T
gene_list = list(patient_counts_inv[:,0])
gene_list.pop(0)

patients_1 = []
patients_2 = []
sim_score = []

patient_counts = patient_counts.tolist()
patient_counts_columns = patient_counts.pop(0)
patient_counts = np.array(patient_counts)
for patient_1 in patient_counts:
	patient_1 = list(np.array(patient_1).reshape(-1,))
	for patient_2 in patient_counts:
		patient_2 = list(np.array(patient_2).reshape(-1,))
		if patient_1[0] != patient_2[0]:
			total_mismatches = 0
			total_zero_matches = 0
			total_one_matches = 0
			patients_1.append(patient_1[0])
			patients_2.append(patient_2[0])
			for i in range(1, len(patient_1)):
				if patient_1[i] == patient_2[i]:
					if int(patient_1[i]) == 0:
						total_zero_matches+=1
					else:
						total_one_matches+=1
				else:	
					total_mismatches+=1	
			sim_score.append(float(total_one_matches)/float((len(patient_1) - 1) - total_zero_matches))

simularity_stuff = np.column_stack((patients_1, patients_2, sim_score))
simularity_stuff = np.matrix(simularity_stuff)
np.savetxt("Similarity_Statistics_Isoform_Level_all.csv", simularity_stuff, delimiter=",", newline = "\n", fmt="%s") #save to a file
relevant_rows = []
compliments = []
for row_num in range(0, len(simularity_stuff)):
	row = simularity_stuff[row_num]
	row = list(np.array(row).reshape(-1,))
	row_reverse = [row[1], row[0], row[2]]
	compliments.append(row_reverse)
	if row not in compliments:
		if float(row[2]) >= np.percentile(sim_score, [90]):
			relevant_rows.append(row)
relevant_rows = np.array(relevant_rows)
np.savetxt("Similarity_Statistics_Isoform_Level.csv", relevant_rows, delimiter=",", newline = "\n", fmt="%s") #save to a file

#Similarity graphs:
import networkx as nx #https://networkx.github.io/documentation/networkx-1.10/reference/introduction.html
def sumColumn(matrix):
    return np.sum(matrix, axis=1)

def sumRow(matrix):
    return np.sum(matrix, axis=0)

def square(list):
    return [i ** 2 for i in list]

from networkx.utils import dict_to_numpy_array
from networkx.algorithms.assortativity.pairs import node_degree_xy, \
    node_attribute_xy
_author_ = ' '.join(['Aric Hagberg <aric.hagberg@gmail.com>'])
_all_ = ['attribute_mixing_matrix',
           'attribute_mixing_dict',
           'degree_mixing_matrix',
           'degree_mixing_dict',
           'numeric_mixing_matrix',
           'mixing_dict']

def attribute_mixing_dict(G, attribute, nodes=None, normalized=False):
    xy_iter = node_attribute_xy(G, attribute, nodes)
    return mixing_dict(xy_iter, normalized=normalized)

def attribute_mixing_matrix(G, attribute, nodes=None, mapping=None, normalized=True):
    """Return mixing matrix for attribute.
    Parameters
    ----------
    G : graph
       NetworkX graph object.
    attribute : string
       Node attribute key.
    nodes: list or iterable (optional)
        Use only nodes in container to build the matrix. The default is
        all nodes.
    mapping : dictionary, optional
       Mapping from node attribute to integer index in matrix.
       If not specified, an arbitrary ordering will be used.
    normalized : bool (default=True)
       Return counts if False or probabilities if True.
    Returns
    -------
    m: numpy array
       Counts or joint probability of occurrence of attribute pairs.
    """
    d = attribute_mixing_dict(G, attribute, nodes)
    a = dict_to_numpy_array(d, mapping=mapping)
    if normalized:
        a = a / a.sum()
    return a

def mixing_dict(xy, normalized=False):
    """Return a dictionary representation of mixing matrix.
    Parameters
    ----------
    xy : list or container of two-tuples
       Pairs of (x,y) items.
    attribute : string
       Node attribute key

    normalized : bool (default=False)
       Return counts if False or probabilities if True.
    Returns
    -------
    d: dictionary
       Counts or Joint probability of occurrence of values in xy.
    """
    d = {}
    psum = 0.0
    for x, y in xy:
        if x not in d:
            d[x] = {}
        if y not in d:
            d[y] = {}
        v = d[x].get(y, 0)
        d[x][y] = v + 1
        psum += 1

    if normalized:
        for k, jdict in d.items():
            for j in jdict:
                jdict[j] /= psum
    return d

def exped_standard_dev(Graph):
	matrix = attribute_mixing_matrix(Graph, 'type', nodes=None, mapping=None, normalized=True)
	bi = sumColumn(matrix)
	ai = sumRow(matrix)
	bi_2 = square(bi)
	ai_2 = square(ai)
	AB = sum([a*b for a,b in zip(ai,bi)])
	AB_2 = AB**2
	A2B = sum([a*b for a,b in zip(ai_2,bi)])
	AB2 = sum([a*b for a,b in zip(ai,bi_2)])
	AB_min = 1 - AB
	M_div = 1/Graph.number_of_edges()
	stand_Dev_top = AB + AB_2 - A2B - AB2 
	stand_dev_al = stand_Dev_top/AB_min
	stand_dev_final = stand_dev_al*M_div
	return stand_dev_final

#all 4 individual isoforms
G=nx.MultiGraph()
nodes = []
for row in relevant_rows:
	row = list(np.array(row).reshape(-1,))
	node_1 = str(row[0])
	nodes.append(node_1)

nodes = list(set(nodes))
G.add_nodes_from(nodes)

Patient_States = {'G94':'control', 'G93':'control', 'G92':'Good', 'G88':'control', 'G85':'Good', 'G82':'pre', 'G81':'Good', 'G79':'Bad', 'G76':'pre', 'G75':'Bad', 'G72':'control', 'G69':'pre', 'G67':'Good', 'G58':'control', 'G57':'control', 'G54':'pre', 'G43':'pre', 'G42':'control', 'G36':'pre', 'G33a':'Bad', 'G30':'pre', 'G29a':'Good', 'G17':'Bad', 'G15':'Bad', 'G138':'pre', 'G124':'pre', 'G11':'Bad', 'G102':'Good', 'G10':'control'}

for row in relevant_rows:
	row = list(np.array(row).reshape(-1,))
	node_1 = str(row[0])
	nodes.append(node_1)
	node_2 = str(row[1])
	edgey = float(row[2])*10
	G.add_edge(node_1, node_2, length=int(edgey))

carac = pd.DataFrame({ 'ID':['G94', 'G93', 'G92','G88','G85', 'G82', 'G81', 'G79', 'G76', 'G75', 'G72', 'G69', 'G67', 'G58', 'G57', 'G54', 'G43', 'G42', 'G36', 'G33a', 'G30', 'G29a', 'G17', 'G15', 'G138', 'G124', 'G11', 'G102', 'G10'], 'myvalue':['control', 'control', 'Good', 'control', 'Good', 'pre', 'Good', 'Bad', 'pre', 'Bad', 'control', 'pre', 'Good', 'control', 'control', 'pre', 'pre', 'control', 'pre', 'Bad', 'pre', 'Good', 'Bad', 'Bad', 'pre', 'pre', 'Bad', 'Good', 'control'] })
carac= carac.set_index('ID')
carac=carac.reindex(G.nodes())

carac['myvalue']=pd.Categorical(carac['myvalue'])
carac['myvalue'].cat.codes

labels = nx.get_edge_attributes(G,'length')
nx.set_node_attributes(G, Patient_States, 'type')
nx.draw(G, node_color=carac['myvalue'].cat.codes, cmap=plt.cm.Set1, with_labels=True, font_size = 9, node_size=700, edge_labels=labels)
plt.savefig("All_subsets_isoforms.png") # save as png
plt.show() # display

print("The attribute assortivity coefficient for all subsets: ", nx.attribute_assortativity_coefficient(G,'type')) #-0.04137648504711185
print("The expected standard deviation for this is: ", exped_standard_dev(G))

#similarity graph for cancer vs control
F=nx.MultiGraph()
nodes = []

for row in relevant_rows:
	row = list(np.array(row).reshape(-1,))
	node_1 = str(row[0])
	nodes.append(node_1)
nodes = list(set(nodes))
F.add_nodes_from(nodes)

Patient_States = {'G94':'control', 'G93':'control', 'G92':'cancer', 'G88':'control', 'G85':'cancer', 'G82':'cancer', 'G81':'cancer', 'G79':'cancer', 'G76':'cancer', 'G75':'cancer', 'G72':'control', 'G69':'cancer', 'G67':'cancer', 'G58':'control', 'G57':'control', 'G54':'cancer', 'G43':'cancer', 'G42':'control', 'G36':'cancer', 'G33a':'cancer', 'G30':'cancer', 'G29a':'cancer', 'G17':'cancer', 'G15':'cancer', 'G138':'cancer', 'G124':'cancer', 'G11':'cancer', 'G102':'cancer', 'G10':'control'}

for row in relevant_rows:
	row = list(np.array(row).reshape(-1,))
	node_1 = str(row[0])
	nodes.append(node_1)
	node_2 = str(row[1])
	edgey = float(row[2])*10
	F.add_edge(node_1, node_2, weight=int(edgey))

carac = pd.DataFrame({ 'ID':['G94', 'G93', 'G92','G88','G85', 'G82', 'G81', 'G79', 'G76', 'G75', 'G72', 'G69', 'G67', 'G58', 'G57', 'G54', 'G43', 'G42', 'G36', 'G33a', 'G30', 'G29a', 'G17', 'G15', 'G138', 'G124', 'G11', 'G102', 'G10'], 'myvalue':['control', 'control', 'cancer', 'control', 'cancer', 'cancer', 'cancer', 'cancer', 'cancer', 'cancer', 'control', 'cancer', 'cancer', 'control', 'control', 'cancer', 'cancer', 'control', 'cancer', 'cancer', 'cancer', 'cancer', 'cancer', 'cancer', 'cancer', 'cancer', 'cancer', 'cancer', 'control'] })
carac= carac.set_index('ID')
carac=carac.reindex(F.nodes())

carac['myvalue']=pd.Categorical(carac['myvalue'])
carac['myvalue'].cat.codes
nx.set_node_attributes(F, Patient_States, 'type')
nx.draw(F, node_color=carac['myvalue'].cat.codes, cmap=plt.cm.Set1, with_labels=True, font_size = 9, node_size=700)
plt.savefig("cancer_vs_control_isoforms.png") #save as png
plt.show() #display

print("The attribute assortivity coefficient for cancer vs. control: ", nx.attribute_assortativity_coefficient(F,'type'))
print("The expected standard deviation for this is: ", exped_standard_dev(F))

#Pre vs not pre
I=nx.MultiGraph()
nodes = []

for row in relevant_rows:
	row = list(np.array(row).reshape(-1,))
	node_1 = str(row[0])
	nodes.append(node_1)
nodes = list(set(nodes))
I.add_nodes_from(nodes)

Patient_States = {'G94':'not', 'G93':'not', 'G92':'not', 'G88':'not', 'G85':'not', 'G82':'pre', 'G81':'not', 'G79':'not', 'G76':'pre', 'G75':'not', 'G72':'not', 'G69':'pre', 'G67':'not', 'G58':'not', 'G57':'not', 'G54':'pre', 'G43':'pre', 'G42':'not', 'G36':'pre', 'G33a':'not', 'G30':'pre', 'G29a':'not', 'G17':'not', 'G15':'not', 'G138':'pre', 'G124':'pre', 'G11':'not', 'G102':'not', 'G10':'not'}


for row in relevant_rows:
	row = list(np.array(row).reshape(-1,))
	node_1 = str(row[0])
	nodes.append(node_1)
	node_2 = str(row[1])
	edgey = float(row[2])*10
	I.add_edge(node_1, node_2, weight=int(edgey))

carac = pd.DataFrame({ 'ID':['G94', 'G93', 'G92','G88','G85', 'G82', 'G81', 'G79', 'G76', 'G75', 'G72', 'G69', 'G67', 'G58', 'G57', 'G54', 'G43', 'G42', 'G36', 'G33a', 'G30', 'G29a', 'G17', 'G15', 'G138', 'G124', 'G11', 'G102', 'G10'], 'myvalue':['not', 'not', 'not', 'not', 'not', 'pre', 'not', 'not', 'pre', 'not', 'not', 'pre', 'not', 'not', 'not', 'pre', 'pre', 'not', 'pre', 'not', 'pre', 'not', 'not', 'not', 'pre', 'pre', 'not', 'not', 'not'] })
carac= carac.set_index('ID')
carac=carac.reindex(G.nodes())

carac['myvalue']=pd.Categorical(carac['myvalue'])
carac['myvalue'].cat.codes
nx.set_node_attributes(I, Patient_States, 'type')
nx.draw(I, node_color=carac['myvalue'].cat.codes, cmap=plt.cm.Set1, with_labels=True, font_size = 9, node_size=700)
plt.savefig("Pre-chemo_vs_everything_else_isoforms.png") #save as png
plt.show() #display


print("The attribute assortivity coefficient for pre-chemo vs. none are: ", nx.attribute_assortativity_coefficient(I,'type')) 
print("The expected standard deviation for this is: ", exped_standard_dev(I))

#Bar Graphs
int_totals_per_gene = []
for item in totals_per_gene:
     int_totals_per_gene.append(int(item))
everything = np.column_stack((gene_list, int_totals_per_gene))
everything = everything.tolist()

everything = sorted(everything, key=lambda ele: float(ele[1]), reverse=True)

everything = np.array(everything)
everything = everything.astype('S140')
everything = np.matrix(everything)

one = 0
two = 0
three = 0
four = 0
five = 0
six = 0
seven = 0
eight = 0
nine = 0	
ten = 0
eleven = 0
twelve = 0
thirteen = 0
fourteen = 0
fifth = 0
sixteen = 0
seventeen = 0
eighteen = 0
nineteen = 0
twenty = 0	
twenty_one = 0	
twenty_two = 0	
twnety_three = 0
twenty_four = 0
twenty_five = 0	
twenty_six = 0
twnety_seven = 0
twenty_eight = 0	
twnety_nine = 0

for row in everything:
	row = list(np.array(row).reshape(-1,))
	if int(row[1]) == 1:
		one +=1
	elif int(row[1]) == 2:
		two +=1
	elif int(row[1]) == 3:
		three +=1
	elif int(row[1]) == 4:
		four +=1
	elif int(row[1]) == 5:
		five +=1
	elif int(row[1]) == 6:
		six +=1
	elif int(row[1]) == 7:
		seven +=1
	elif int(row[1]) == 8:
		eight +=1
	elif int(row[1]) == 9:
		nine +=1
	elif int(row[1]) == 10:
		ten +=1
	elif int(row[1]) == 11:
		eleven +=1
	elif int(row[1]) == 12:
		twelve +=1
	elif int(row[1]) == 13:
		thirteen +=1
	elif int(row[1]) == 14:
		fourteen +=1
	elif int(row[1]) == 15:
		fifth +=1
	elif int(row[1]) == 16:
		sixteen +=1
	elif int(row[1]) == 17:
		seventeen +=1
	elif int(row[1]) == 18:
		eighteen +=1
	elif int(row[1]) == 19:
		nineteen +=1
	elif int(row[1]) == 20:
		twenty +=1	
	elif int(row[1]) == 21:
		twenty_one +=1	
	elif int(row[1]) == 22:
		twenty_two +=1	
	elif int(row[1]) == 23:
		twnety_three +=1
	elif int(row[1]) == 24:
		twenty_four +=1	
	elif int(row[1]) == 25:
		twenty_five +=1	
	elif int(row[1]) == 26:
		twenty_six +=1	
	elif int(row[1]) == 27:
		twnety_seven +=1
	elif int(row[1]) == 28:
		twenty_eight +=1	
	elif int(row[1]) == 29:
		twnety_nine +=1

print(one, two, three, four, five, six, seven, eight, nine, ten, eleven, twelve, thirteen, fourteen, fifth, sixteen, seventeen, eighteen, nineteen, twenty, twenty_one, twenty_two, twnety_three, twenty_four, twenty_five, twenty_six, twnety_seven, twenty_eight, twnety_nine)



my_Set = everything[:,1].tolist()
new_list = []
for item in my_Set:
	for sub_item in item:
		sub_item = int(sub_item)
		new_list.append(sub_item)



import matplotlib.pyplot as plt

plt.hist(new_list, bins=29, align="left")
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29])
plt.yticks(range(0, 6000, 250))
plt.title("B. Isoforms")
plt.xlabel("Number of Isoforms")
plt.ylabel("Number of Proteins")
plt.show()


titles = ["Isoforms", "Frequency"]
everything = np.insert(everything, 0, titles, axis = 0)
decoded_everything = []
for line in everything:
	line = list(np.array(line).reshape(-1,))
	line = [el.decode('UTF-8') for el in line]
	decoded_everything.append(line)
decoded_everything = np.array(decoded_everything)
np.savetxt("Total_Isoform_Counts.csv", decoded_everything, delimiter=",", newline = "\n", fmt="%s")
new_gene_list = []
for suby in decoded_everything:
    suby = list(np.array(suby).reshape(-1,))
    new_gene_list.append(str(suby[0]))






#df=pd.read_csv('Total_Isoform_Counts.csv')
#days_values=df.set_index("Isoforms").loc[new_gene_list].plot(kind="bar", legend=False)
#plt.xlabel('Isoforms')
#plt.ylabel('Frequency')
#plt.title('The total count of each isoform in the patient subset')
#plt.show() #overall not broken down into subcategories

from IPython import get_ipython
ipy = get_ipython()
if ipy is not None:
    ipy.run_line_magic('matplotlib', 'inline')
sales=pd.read_csv("Total_Isoform_Counts.csv")
reduced_sales = sales[['Isoforms','Frequency']]
ax = reduced_sales.plot.bar(x='Isoforms', y='Frequency', rot=0)
fig = ax.get_figure()
fig.savefig("isoforms_ordered.png", dpi=200)

#y_pos = np.arange(len(gene_list)) 
#plt.bar(y_pos, int_totals_per_gene, align='center', alpha=0.5)
#plt.xticks(y_pos, gene_list)
#plt.ylabel('Frequency')
#plt.title('Isoform distribution')
plt.show()




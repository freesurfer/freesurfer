import csv
import os
import math
from scipy import stats
import numpy as np
import nibabel as nib
import os.path
import scipy.spatial
import glob
import subprocess
import pandas as pd
from nipy.modalities.fmri.glm import GeneralLinearModel

def getCorrespondingClusters(correspondance,order=True):
		corr=dict()
		distances=dict()
		indeces=dict()
		ind=0
		with open(correspondance, 'r') as csvfile:
				corReader= csv.reader(csvfile, delimiter=',', quotechar='|')	
				header=True
				for row in corReader:
					if not header and len(row)>1:
						if order:
							corr[row[1]]=row[0]
							indeces[ind]=row[1]
							distances[row[1]]=float(row[2])
						else:
							corr[row[0]]=row[1]
							indeces[ind]=row[0]
							distances[row[0]]=float(row[2])
						ind+=1
					else:
						header=False
		return corr, distances, indeces

def averageCorrespondingClusters(correspondences, imagesFolder, outputFolder,  clusterIndeces):
		averages=dict()
		norm=dict()
		for s_i, correspondance in enumerate(correspondences):
				try:
					for clusterIndex in clusterIndeces:
						corr, distances, indeces =  getCorrespondingClusters(correspondance, True)
						#print(correspondance)
						image=imagesFolder[s_i]+""+indeces[clusterIndex]+".nii.gz"
						im =nib.load(image)
						b= im.get_data()
						b = b/b.max()
						b = np.ceil(b) 
						if clusterIndex in averages:
							b += averages[clusterIndex].get_data() 
							averages[clusterIndex]=nib.Nifti1Image(b, averages[clusterIndex].get_affine())
							norm[clusterIndex]+=1
						else:
							averages[clusterIndex] = nib.Nifti1Image(b, im.get_affine())    
							norm[clusterIndex]=1
				except Exception as e:
					print(str(e))

		for clusterIndex in clusterIndeces:
				directory=outputFolder
				if not os.path.exists(directory):
					os.makedirs(directory)
				data=averages[clusterIndex].get_data()/norm[clusterIndex]
				nib.Nifti1Image(data, averages[clusterIndex].get_affine()).to_filename(directory+"/"+str(clusterIndex)+'.nii.gz')
				print ("saving",directory+"/"+str(clusterIndex)+'.nii.gz')

def readTree(numNodes, histogramFile,header=True):
    almostFoundAllClusters=False
    foundAllClusters=False
    nodes_childs=dict()
    whos_dad=dict()
    
    with open(histogramFile, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')	
        clusters=set()

        for row in reader:
            if not header:
                if not foundAllClusters:
                    try:
                        if row[0] in clusters:
                            clusters.remove(row[0])
                        clusters.add(row[1])
                        if len(clusters)==numNodes :
                            foundAllClusters=True
                            for i in clusters:
                                nodes_childs[i]=[]
                            
                    except:
                        #print ("header")
                        None
                else:
                    if row[0] in whos_dad :
                        dad = whos_dad[row[0]]
                        if row[0] in nodes_childs[dad]:
                            nodes_childs[dad].remove(row[0])
                        nodes_childs[dad].append(row[1])
                        whos_dad[row[1]]=dad	
                    else:
                        nodes_childs[row[0]].append(row[1])
                        whos_dad[row[1]]=row[0]
            else:
                header=False
            
    return nodes_childs, whos_dad

#dti_measures=["FA","ADC","RD","AD"]
def groupAnalysis( headers, cols , groups_classification, clustersToAnalyze, subjects_dir, target_subject):
	with open(groups_classification) as f:
		groups_cat = dict(filter(None, csv.reader(f, delimiter=',')))
	significant_childs=set()	
	for c_i, clusterNum in enumerate(clustersToAnalyze):
		
		childs, dads = readTree(clusterNum, f"{subjects_dir}/{target_subject}/HierarchicalHistory.csv")
		#print(childs, dads)
		order_nodes= pd.read_csv(f"{subjects_dir}/{target_subject}/measures/{target_subject}_{target_subject}_c{clusterNum}.csv",delimiter=",", header=0,usecols=[0])
		#print(order_nodes["Cluster"][0])
		ys=[]
		for a in headers:
			ys.append([])
			
		for i in range(clusterNum):
			for j in range(len(headers)):
				ys[j].append([])
		X=[]
		
		for s in groups_cat.keys():
			measures=f"{subjects_dir}/{s}/measures/{target_subject}_{s}_c{clusterNum}.csv"
			data = pd.read_csv(measures,delimiter=",", header=0, names=headers, usecols=cols)
			#print(measures)
			if len(data[headers[0]])>= clusterNum:
					for i,h in enumerate(headers):
						for j in range(clusterNum):
							ys[i][j].append(data[h][j])

					if int(groups_cat[s])>0:
						X=np.append(X,[1, 0])
					else:
						X=np.append(X,[0, 1])

		X= np.array(X).reshape(len(ys[0][0]),2)
		for i, m  in enumerate(headers):
			Y=ys[i][0]
			for j in range(1,len(ys[i])):
				#print(np.shape(ys[i][j]))
				Y = np.vstack((Y,[ys[i][j]]))
			
			Y=np.array(Y).transpose()
			#print(np.shape(Y))
			cval = np.hstack((-1, 1))
			model = GeneralLinearModel(X)
			model.fit(Y)
			p_vals = model.contrast(cval).p_value() # z-transformed statistics
			#print( p_vals)
			for index, p in enumerate(p_vals):
				if p*len(ys[0]) <0.05:
					if len(childs[str(order_nodes["Cluster"][index])]) ==0:
						significant_childs.add(str(order_nodes["Cluster"][index]))
					else:
						for c in childs[str(order_nodes["Cluster"][index])]:
							significant_childs.add(c)
					
				if clusterNum == clustersToAnalyze[-1] and p<0.05 and str(order_nodes["Cluster"][index]) in significant_childs:
					print(m,index,p, p,str(order_nodes["Cluster"][index]))
					
				 
			cval = np.hstack((1, -1))
			model = GeneralLinearModel(X)
			model.fit(Y)
			p_vals = model.contrast(cval).p_value() # z-transformed statistics
			#print( p_vals)
			for index, p in enumerate(p_vals):
				if p*len(ys[0]) <0.05:
					#print(index,p, p*len(ys[0])*4,str(order_nodes["Cluster"][index])) 
					if len(childs[str(order_nodes["Cluster"][index])]) ==0:
						significant_childs.add(str(order_nodes["Cluster"][index]))
					else:
						for c in childs[str(order_nodes["Cluster"][index])]:
							significant_childs.add(c)
		    #plt.show()
				if clusterNum == clustersToAnalyze[-1] and p<0.05 and str(order_nodes["Cluster"][index]) in significant_childs:
					print(m, index,p, p,str(order_nodes["Cluster"][index]))
					
#groupAnalysis(headers=["meanFA","meanADC","meanRD","meanAD"],cols=[2, 6,10,14], groups_classification="/space/vault/7/users/vsiless/lilla/classification.csv", clustersToAnalyze=[50,100, 150,200],target_subject="INF007",subjects_dir="/space/vault/7/users/vsiless/lilla/AnatomiCuts/babybold/")
groupAnalysis(headers=["meanFA"],cols=[2], groups_classification="/space/vault/7/users/vsiless/lilla/classification.csv", clustersToAnalyze=[50,100, 150,200],target_subject="INF007",subjects_dir="/space/vault/7/users/vsiless/lilla/AnatomiCuts/babybold/")
groupAnalysis(headers=["meanADC"],cols=[6], groups_classification="/space/vault/7/users/vsiless/lilla/classification.csv", clustersToAnalyze=[50,100, 150,200],target_subject="INF007",subjects_dir="/space/vault/7/users/vsiless/lilla/AnatomiCuts/babybold/")
groupAnalysis(headers=["meanRD"],cols=[10], groups_classification="/space/vault/7/users/vsiless/lilla/classification.csv", clustersToAnalyze=[50,100, 150,200],target_subject="INF007",subjects_dir="/space/vault/7/users/vsiless/lilla/AnatomiCuts/babybold/")
groupAnalysis(headers=["meanAD"],cols=[14], groups_classification="/space/vault/7/users/vsiless/lilla/classification.csv", clustersToAnalyze=[50,100, 150,200],target_subject="INF007",subjects_dir="/space/vault/7/users/vsiless/lilla/AnatomiCuts/babybold/")

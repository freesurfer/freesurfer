#!/usr/bin/env python
# Andrew Zhang
# Professor Siless
# July 2019
# plots.py

#Now be able to take in multiple csv files, calculate the average of each column for each file, then plot the averages of each column for each file
#Number of plots in violin = Number of files
#Number of violins = number of columns

# Find the matching file names in each subject's directory, take the average, then plot

import sys
import os
import matplotlib.pyplot as plt
import csv
import pandas as pd

# Check for correct inputs
if len(sys.argv) != 4:
    print('Usage: ./executable directory labels structure')
    exit()

# Assign each infant file into a group

# Take in the names and group numbers from the text csv file
group_df = pd.read_csv(sys.argv[2])
group_col = group_df.columns

data = []
for c in group_col:
    data.append(group_df[c].tolist())

#print(data[1][0])

#print(len(os.listdir(sys.argv[1])))

# Find the highest number group needed
max_group = 1
for i in range(len(data[0])):
    #print (data[1][i], max_group)
    if data[1][i] > max_group:
        max_group = data[1][i]

# Create a dictionary of lists: a list for each group number
inf_to_group = {k: [] for k in range(1, max_group + 1)}

# Iterate thru the text csv file and sort the groups
for i in range(len(data[0])):
    for j in range(len(os.listdir(sys.argv[1]))):
        print(data[0][i], os.listdir(sys.argv[1])[j])
        if data[0][i] == os.listdir(sys.argv[1])[j]:
            print(data[1][i], data[0][i])
            inf_to_group[data[1][i]].append(data[0][i])
    
print(inf_to_group)
print(inf_to_group[2])


# Find the files (left/right) of the structure named in the input

# Collect the names of the csv files based on order of groups
directory = sys.argv[1]
file_names = {k: [] for k in range(1, max_group + 1)}

# Case if key has no value?
for i in range(1, len(inf_to_group) + 1):
    for j in range(len(inf_to_group[i])):
        print (inf_to_group[i])
        directory = inf_to_group[i][j]
        print(directory)

# Place this loop into previous one
for subdir, dirs, files in os.walk(directory):
    for file in files:
        filepath = subdir + os.sep + file

        if filepath.endswith(sys.argv[3] + '.csv'):
            print (filepath)
            file_names.append(filepath)

print (file_names)

# Calculate the means of each file

# 4 means total: G1L, G1R, G2L, G2R


'''
# Read in one file to get the dimensions
sample_df = pd.read_csv(sys.argv[1])

width = len(sample_df.columns) - 1
height = len(sys.argv) - 1

# Create a 2D array to hold the mean values
data = [[0 for x in range(height)] for y in range(width)] 

# One loop to iterate thru infants
# Another to iterate thru files
# Another to hold the means

# Loop through the csv files to extract all the mean values from each of its column
for i in range(1, height + 1):
    # Read in the csv file
    df = pd.read_csv(sys.argv[i])
    
    means = df.mean(axis = 0)

    print (means[1])

'''

'''
    # Read in data by columns. This assumes that the first column has irrelevant data
    j = 0
    columns = df.columns
    data = []
    for c in columns:
        if j != 0:
            data.append(df[c].tolist())
            j+=1
        else:
            j+=1

    # Calculate the means and insert into the 2D array
    for k in range(width):
        means[k][i-1] = sum(data[k]) / len(data[k])
'''
'''
# Creates the basis of the plots
# Increase ncols and add [#] to axes to add more plots
#fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(20, 20))

xticklabels=list(sample_df.columns.values.tolist()[1:])


for i in range(width):
    plt.figure(xticklabels[i])
    plt.violinplot(means[i], showmeans=True)
    plt.xlabel('Type of Value')
    plt.ylabel('Observed Values')
'''

'''
# Violin plot
#for i in range(width):
axes.violinplot(means, showmeans=True, showmedians=False)

# Box plot
# axes[1].boxplot(means)

# Adding horizontal grid lines and general labels
# Loop if more than one plot is desired

i = 0
#for ax in axes:
axes.yaxis.grid(True)
axes.set_xticks([y + 1 for y in range(len(means))])
#ax.set_xticks([y for y in range(2)])
axes.set_xlabel('Types of Values')
axes.set_ylabel('Observed Values')
axes.set_title('CSV Comparison')
#i+=1
  
# Add x tick labels
plt.setp(axes, xticks=[y + 1 for y in range(len(means))], 
         xticklabels=list(sample_df.columns.values.tolist()[1:]))
'''
# Print out the plots
#plt.show()

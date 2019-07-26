#!/usr/bin/env python
# Andrew Zhang
# Professor Siless
# July 2019
# plots.py

#Now be able to take in multiple csv files, calculate the average of each column for each file, then plot the averages of each column for each file
#Number of plots in violin = Number of files
#Number of violins = number of columns

import sys
import matplotlib.pyplot as plt
import csv
import pandas as pd

# Read in one file to get the dimensions
sample_df = pd.read_csv(sys.argv[1])

width = len(sample_df.columns) - 1
height = len(sys.argv) - 1

# Create a 2D array to hold the mean values
means = [[0 for x in range(height)] for y in range(width)] 

# Loop through the csv files to extract all the mean values from each of its column
for i in range(1, height + 1):
    # Read in the csv file
    df = pd.read_csv(sys.argv[i])
    
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


# Creates the basis of the plots
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(20, 20))

# Violin plot
axes[0].violinplot(means,
                   showmeans=True, 
                   showmedians=False)
# Box plot
axes[1].boxplot(means)

# Adding horizontal grid lines and general labels
for ax in axes:
    ax.yaxis.grid(True)
    ax.set_xticks([y + 1 for y in range(len(means))])
    ax.set_xlabel('Type of Value')
    ax.set_ylabel('Observed values')
    ax.set_title('Plots of All CSV Values')

# Add x tick labels
plt.setp(axes, xticks=[y + 1 for y in range(len(means))], 
         xticklabels=list(sample_df.columns.values.tolist()[1:]))

# Print out the plots
plt.show()

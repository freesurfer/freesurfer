#!/usr/bin/env python
# Andrew Zhang
# Professor Siless
# July 2019
# plots.py

# This program takes in a directory of subjects, a text csv file that organizes the subjects into
# groups, and the name of a structure that the user wants data for. The program first organizes 
# the subjects into the numbered groups. For each subject, the corresponding files for the named
# structure are extracted and organized by whether they're in the left or right hemisphere. The
# files are then plotted onto violin plots that directly compare each type of mean value and show
# a plot for each group number and both sides of each group. 

# Import libraries
import sys
import os
import matplotlib.pyplot as plt
import csv
import pandas as pd

# Check for correct inputs
if len(sys.argv) != 4:
    print('Usage: ./executable directory labels structure')
    exit()


# Take in the names and group numbers from the text csv file
group_df = pd.read_csv(sys.argv[2])
group_col = group_df.columns

data = []
for c in group_col:
    data.append(group_df[c].tolist())

# Find the highest numerical group in the text csv
max_group = 1
for i in range(len(data[0])):
    if data[1][i] > max_group:
        max_group = data[1][i]


# Dictionary to order the infants in accordance to their group number
inf_to_group = {k: [] for k in range(1, max_group + 1)}

# Iterate through the text csv file and sort the groups
for i in range(len(data[0])):
    for j in range(len(os.listdir(sys.argv[1]))):
        if data[0][i] == os.listdir(sys.argv[1])[j]:
            inf_to_group[data[1][i]].append(data[0][i])


# Collect the names of the csv files based on order of groups
file_names = {k: [] for k in range(1, max_group + 1)}

main_directory = sys.argv[1]

# Loop through the groups
for i in range(1, len(inf_to_group) + 1):
    
    # Loop through each subject within a group
    for j in range(len(inf_to_group[i])):
        subject = inf_to_group[i][j]

        # Iterate through the main directory from the command line
        for subdir, dirs, files in os.walk(main_directory): 
            
            # Get the name of each file
            for file in files:
                filepath = subdir + os.sep + file

                # Filter out the files based on the input from command line
                if filepath.endswith(sys.argv[3] + '.csv') and subject in filepath:
                    file_names[i].append(filepath)


# Take in one sample csv file to find amount of plots needed
sample_df = pd.read_csv(file_names[1][0])

num_plots = len(sample_df.columns) - 1

# Set dimensions of plots. Height concerns both sides, so must be doubled what the largest group is
height = 2 * max_group


# Collect the value type names. This assumes the first column name is irrelevant
xticklabels=list(sample_df.columns.values.tolist()[1:])

# Create labels for the violin plots
violin_labels = ['' for x in range(2*max_group)]

# Each group has 2 labels: for both left and right
i = 2
while i < 2 * max_group + 1:
    violin_labels[i-2] += 'Group ' + str(int(i/2)) + ' Left'
    violin_labels[i-1] += 'Group ' + str(int(i/2)) + ' Right'
    i+=2


# Create a list of lists, two for each group
for k in range(num_plots):
    # Odd indexes will hold left side data, even ones hold right side
    single_plot = [[] for x in range(height)]


    # Loop through each group 
    for i in range(1, len(file_names) + 1):
        
        # Loop though each file within each group
        for j in range(len(file_names[i])):

            # Find files of the left side
            if 'left' in file_names[i][j].lower() or 'lh' in file_names[i][j].lower():
                
                df = pd.read_csv(file_names[i][j])
    
                means = df.mean(axis = 0)

                #print (xticklabels[k])
                #print (violin_labels[2 * i - 2], means[k])

                # Isolate the desired mean value for this specific plot and append it
                single_plot[2 * i - 2].append(means[k])


            # Same idea as code for mapping the left side values
            elif 'right' in file_names[i][j].lower() or 'rh' in file_names[i][j].lower():

                df = pd.read_csv(file_names[i][j])

                means = df.mean(axis = 0)

                #print (xticklabels[k])
                #print (violin_labels[2 * i - 1], means[k])

                # On a separate row and with a separate column iterator, fill in data
                single_plot[2 * i - 1].append(means[k])


    #Violin plots
    plt.figure(xticklabels[k]) 
    plt.violinplot(single_plot, showmeans=True)

    # Plot logistics
    plt.title(xticklabels[k])
    plt.xlabel('Group')
    plt.ylabel('Observed Values')

    # Label the x axis more specifically
    plt.xticks([i+1 for i in range(2 * max_group)], violin_labels)


# Print out all the plots
plt.show()


# Code to read in one csv file and isolate the different types of values
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

# Code to print out both a violin plot and box plot with all data represented in both
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
# plt.show()

# count tandem copies 
import pandas as pd
import numpy as np
import os

# getting a list of pan gene that has NA in the matrix 
# information include pan gene ID, genome it is coming from, and genome where it was missing
infile = 'conflict_matrix.csv'
pan_matrix = open(infile, "r")

# convert list to string 
def listToString(s): 
    # initialize an empty string
    str1 = " " 
    # return string  
    return (str1.join(s))

# Function to print the intersection
def findIntersection(intervals,N):
    # First interval
    l = intervals[0][0]
    r = intervals[0][1]
    # Check rest of the intervals
    # and find the intersection
    for i in range(1,N):
        # If no intersection exists
        if (intervals[i][0] > r or intervals[i][1] < l):
            return("no intersection")
        # Else update the intersection
        else:
            l = max(l, intervals[i][0])
            r = min(r, intervals[i][1])
            return(str(l)+".."+str(r))

# get header info
header = pan_matrix.readline()
header = header.strip()
header = header.split(",")
# key for pan genome name 

tandem_list = []
for line in pan_matrix:
    pan_gene = line.strip()
    pan_gene = line.split(',')
    for idx, gene_id in enumerate(pan_gene): 
        # split gene_id column by ";", if it can be separated, do this below
        if ";" not in gene_id:
            no_action = pan_gene[0].split(';')[0] + " " +header[idx] + " "+gene_id
            tandem_list.append(no_action)
        if ";" in gene_id:
            split_cluster= gene_id.split(';')
            subs = "intact" 
            count_intact = gene_id.count(subs)    #split TE -id get intact, truncated, and null information
            if count_intact == 1: # when only in intact is present, print this intact coordinates to present the TE ID
                intact = [i for i in split_cluster if subs in i]
                singe_intact = pan_gene[0].split(';')[0] + " " +listToString(intact).split("_")[0] + " " +listToString(intact)
                tandem_list.append(singe_intact)
                # above is working 
            if count_intact >1: # when there are more than 1 intact conflict exist, seek for coordinates overlap
                conflict_merge =[]
                intact_2 = [i for i in split_cluster if subs in i]
                for idx in range(len(intact_2)):
                    format_coord = eval(intact_2[idx].split("_")[1].split(":")[1].split("..")[0] + ","+intact_2[idx].split("_")[1].split(":")[1].split("..")[1])
                    conflict_merge.append(format_coord)
                N =len(conflict_merge) # finding overlaps. if there are overlaps, result will print out as range, if not, it will be -1
                if findIntersection(conflict_merge, N) == "no intersection":
                    converted_list = [str(element) for element in intact_2]
                    joined_string = ";".join(converted_list)
                    unmerged= pan_gene[0].split(';')[0]+ " "+ gene_id.split("_")[0] + " "+ joined_string
                    tandem_list.append(unmerged)
                else:
                    merged_coordinate= pan_gene[0].split(';')[0]+ " "+ gene_id.split("_")[0] + " "+ gene_id.split(":")[0] +":" + findIntersection(conflict_merge, N) + "_merge_intact"
                    tandem_list.append(merged_coordinate)
            if count_intact ==0: # if the conflict pair doesn't have intact TE, move into truncated TEs 
                subs2= "truncated"
                count_truncated = gene_id.count(subs2)
                if count_truncated == 1: # when the conflict pair only has one truncated TE, discard the rest and appened it
                    truncated = [i for i in split_cluster if subs2 in i]
                    singe_truncated = pan_gene[0].split(';')[0] + " " +listToString(truncated).split("_")[0] + " " +listToString(truncated)
                    tandem_list.append(singe_truncated)
                if count_truncated >1:
                    conflict_merge =[]
                    truncated_2 = [i for i in split_cluster if subs2 in i]
                    for idx in range(len(truncated_2)):
                        format_coord = eval(truncated_2[idx].split("_")[1].split(":")[1].split("..")[0] + ","+truncated_2[idx].split("_")[1].split(":")[1].split("..")[1])
                        conflict_merge.append(format_coord)
                    N =len(conflict_merge) # finding overlaps. if there are overlaps, result will print out as range, if not, it will be -1
                    if findIntersection(conflict_merge, N) == "no intersection":
                        converted_list = [str(element) for element in truncated_2]
                        joined_string = ";".join(converted_list)
                        unmerged= pan_gene[0].split(';')[0]+ " "+ gene_id.split("_")[0] + " "+ joined_string
                        tandem_list.append(unmerged)
                    else:
                        merged_coordinate= pan_gene[0].split(';')[0]+ " "+ gene_id.split("_")[0] + " "+ gene_id.split(":")[0] +":" + findIntersection(conflict_merge, N) + "_merge_truncated"
                        tandem_list.append(merged_coordinate)         
                if count_truncated == 0:
                    conflict_merge =[]
                    subs3= "null"
                    null2 = [i for i in split_cluster if subs3 in i]
                    for idx in range(len(null2)):
                        format_coord = eval(null2[idx].split("_")[1].split(":")[1].split("..")[0] + ","+null2[idx].split("_")[1].split(":")[1].split("..")[1])
                        conflict_merge.append(format_coord)
                    N =len(conflict_merge) # finding overlaps. if there are overlaps, result will print out as range, if not, it will be failed
                    if findIntersection(conflict_merge, N) == "no intersection":
                        converted_list = [str(element) for element in null2]
                        joined_string = ";".join(converted_list)
                        unmerged= pan_gene[0].split(';')[0]+ " "+ gene_id.split("_")[0] + " "+ joined_string
                        tandem_list.append(unmerged)
                    else:
                        merged_coordinate= pan_gene[0].split(';')[0]+ " "+ gene_id.split("_")[0] + " "+ gene_id.split(":")[0] +":" + findIntersection(conflict_merge, N) + "_merge_null"
                        tandem_list.append(merged_coordinate)

completeName = os.path.join("tandem_list_for_reshape.txt")
sourceFile = open(completeName, 'w')
print(*tandem_list,file=sourceFile,sep = "\n")
sourceFile.close()

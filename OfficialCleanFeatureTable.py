#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 21 15:55:57 2025

@author: kaylapoirier
"""

import pandas as pd #version 2.2.2
import numpy as np #version 1.26.4
from scipy.stats import mannwhitneyu #version 1.15.2
import os

#Before running, make sure to have a folder specifically for this project and add this .py file in it. Also, add MElaS.Rmd to this folder. In this case, these new files will go into my home directory. 
#takes about 5 minutes to run

#Produce MElaSJR0_1000, MElaSJR0_100, MElaSCL0_1000, MElaSCL0_100 folders. In each of these folders, there should be nested folders of all organs of interest
organs=['Cecum', 'Distal_Colon', 'Duodenum', 'Esophagus' , 'Heart_A', 'Heart_B', 'Heart_C', 'Heart_D', 'Ileum', 'Jejunum', 'Proximal_Colon', 'Upper_Stomach'] #list of organs of interest
folder_type=['JR0_1000', 'JR0_100', 'CL0_1000', 'CL0_100'] #folders that will be made 
for f in folder_type:
    
    for organ in organs:
        nested_directory=f'/Users/kaylapoirier/MElaS11.21/MElaS{f}/{organ}' #location to where folders should be added

        try:
            os.makedirs(nested_directory)
            print(f"Nested directories '{nested_directory}' created successfully.")
        except FileExistsError:
            print(f"One or more directories in '{nested_directory}' already exist.")
        except PermissionError:
            print(f"Permission denied: Unable to create '{nested_directory}'.")
        except Exception as e:
            print(f"An error occurred: {e}")


#clean original dataset    
org_dataset=pd.read_csv('/Users/kaylapoirier/Desktop/DisTol_mergeddata.csv') #load original data    
org_dataset=org_dataset.dropna() #drop NAs
org_dataset= org_dataset[(org_dataset['ATTRIBUTE_infected_parasites'] != 10000)] #remove 10,000 parasite load
#org_dataset= org_dataset[(org_dataset['ATTRIBUTE_endpoint'] != '24w')] #only work with 0, 3 and 12w
org_dataset= org_dataset[(org_dataset['ATTRIBUTE_organ'] != 'Lower_Stomach_2')] #unreliable data
org_dataset= org_dataset[(org_dataset['ATTRIBUTE_organ'] != 'Lower_Stomach')] #unreliable data
org_dataset= org_dataset[(org_dataset['ATTRIBUTE_groupset'] == 'Sample')] #keep only sample and remove QC, blanks, etc data
org_dataset['ATTRIBUTE_endpoint'] = org_dataset['ATTRIBUTE_endpoint'].astype('category')

#make dataframe with specific categories

#JR with 0 vs 1000 paraiste load
JR0_1000=org_dataset[(org_dataset['ATTRIBUTE_parasite_species'] != 'CL_Luc')] #keep JR
JR0_1000.loc[JR0_1000['ATTRIBUTE_infected_parasites'] == 0, 'ATTRIBUTE_infected_parasites'] = "uninfected" #change 0 to uninfected
JR0_1000.loc[JR0_1000['ATTRIBUTE_infected_parasites'] == 1000, 'ATTRIBUTE_infected_parasites'] = "infected" #change 1000 to infected 
JR0_1000=JR0_1000[(JR0_1000['ATTRIBUTE_infected_parasites'] != 100)] #remove 100

#JR with 0 vs 100 paraiste load
JR0_100=org_dataset[(org_dataset['ATTRIBUTE_parasite_species'] != 'CL_Luc')] #keep JR
JR0_100.loc[JR0_100['ATTRIBUTE_infected_parasites'] == 0, 'ATTRIBUTE_infected_parasites'] = "uninfected" #change 0 to uninfected
JR0_100.loc[JR0_100['ATTRIBUTE_infected_parasites'] == 100, 'ATTRIBUTE_infected_parasites'] = "infected" #change 100 to infected
JR0_100=JR0_100[(JR0_100['ATTRIBUTE_infected_parasites'] != 1000)] #remove 1000

#CL with 0 vs 1000 paraiste load
CL0_1000=org_dataset[(org_dataset['ATTRIBUTE_parasite_species'] != 'JR_Luc')] #keep CL
CL0_1000.loc[CL0_1000['ATTRIBUTE_infected_parasites'] == 0, 'ATTRIBUTE_infected_parasites'] = "uninfected" #change 0 to uninfected
CL0_1000.loc[CL0_1000['ATTRIBUTE_infected_parasites'] == 1000, 'ATTRIBUTE_infected_parasites'] = "infected" #change 1000 to infected
CL0_1000=CL0_1000[(CL0_1000['ATTRIBUTE_infected_parasites'] != 100)] #remove 100

#CL with 0 vs 100 paraiste load
CL0_100=org_dataset[(org_dataset['ATTRIBUTE_parasite_species'] != 'JR_Luc')] #keep CL
CL0_100.loc[CL0_100['ATTRIBUTE_infected_parasites'] == 0, 'ATTRIBUTE_infected_parasites'] = "uninfected" #change 0 to uninfected
CL0_100.loc[CL0_100['ATTRIBUTE_infected_parasites'] == 100, 'ATTRIBUTE_infected_parasites'] = "infected" #change 100 to infected
CL0_100=CL0_100[(CL0_100['ATTRIBUTE_infected_parasites'] != 1000)] #remove 1000


positions = org_dataset['ATTRIBUTE_organ'].unique() #make an array of all organ types
infection_time = org_dataset['ATTRIBUTE_endpoint'].unique()  #make an arrayr of all endpoints

datasets={'JR0_1000': JR0_1000, 'JR0_100': JR0_100, 'CL0_1000': CL0_1000, 'CL0_100': CL0_100}

#Calculate fold change of median infected/uninfected at time 0w, 3w and 12w of each organ and produce a file for each organ and time point and save it to individual folders by organ
for name, data in datasets.items():
    
    for organ in positions:
        for time in infection_time:

            filtered_df = data[(data['ATTRIBUTE_organ'] == organ) & (data['ATTRIBUTE_endpoint']==time)] #make dataset by specific organ and timepoint
            filtered_df=filtered_df.drop(['sampleid', 'ATTRIBUTE_endpoint', 'ATTRIBUTE_organ', 'ATTRIBUTE_groupset', 'ATTRIBUTE_parasite_species' ], axis=1) #remove ID, sample type, postions and endpoint, groupset and species
            
            infected_df=filtered_df[filtered_df["ATTRIBUTE_infected_parasites"]=='infected'] #makes dataframe of just infected with specific organ and post day infection
            uninfected_df=filtered_df[filtered_df["ATTRIBUTE_infected_parasites"]=='uninfected'] #makes dataframe of just uninfected with specific organ and post day infection
         
            median_df=filtered_df.groupby("ATTRIBUTE_infected_parasites").median(numeric_only=True) #caluclate median of infected and uninfected of all metabolites
            ratio_df=(median_df.loc["infected"]+0.0000005)/(median_df.loc["uninfected"]+0.0000005) #log(median infected+0.0000005/median uninfected+0.0000005 by each time eg infected 3w/uninfected 3w) smallest metabolite values is 0.00000506
            ratio_df=pd.DataFrame(ratio_df, columns=['FC'])
            
            ratio_df.index.name = 'ID'

            filename = f'/Users/kaylapoirier/MElaS11.21/MElaS{name}/{organ}/newratio_{time}.csv' #make a new file called newratio_3w.csv and newratio_12w.csv for each organ and save it to individual folders by organ
            ratio_df.to_csv(filename)
    
    
            print(f'{name}_{organ}_{time}.csv created')
            
    #Calculate logFC of 3w/0w and logFC 12w/3w and make a file called newratio_3v0.txt and newratio_12v3.txt containing logFC and p-values
    for organ in positions:
    
        df_12=pd.read_csv(f'/Users/kaylapoirier/MElaS11.21/MElaS{name}/{organ}/newratio_12w.csv', index_col=0) #make dataframe of each ratio by week
        df_3=pd.read_csv(f'/Users/kaylapoirier/MElaS11.21/MElaS{name}/{organ}/newratio_3w.csv', index_col=0)
        df_0=pd.read_csv(f'/Users/kaylapoirier/MElaS11.21/MElaS{name}/{organ}/newratio_0w.csv', index_col=0)
        
        ratio3_0= df_3['FC'] / df_0['FC'] #FC of 3w/ 0w
        ratio12_3= df_12['FC'] / df_3['FC'] #FC 12w/FC 3w
        
        ratio12_3=np.log(ratio12_3) #log FC12v3
        ratio12_3=pd.DataFrame(ratio12_3, columns=['FC']) #put FC12v3 into dataframe
        ratio12_3.index.name = 'ID'
        filename12v3 = f'/Users/kaylapoirier/MElaS11.21/MElaS{name}/{organ}/newratio_12v3.txt'
        ratio12_3.to_csv(filename12v3, sep='\t') #save FC12v3 into a txt file into desingated organ file
        
        ratio3_0=np.log(ratio3_0) #log FC3v0
        ratio3_0=pd.DataFrame(ratio3_0, columns=['FC']) #put FC3v12 into dataframe
        ratio3_0.index.name = 'ID'
        filename3v0 = f'/Users/kaylapoirier/MElaS11.21/MElaS{name}/{organ}/newratio_3v0.txt'
        ratio3_0.to_csv(filename3v0, sep='\t') #save FC24v12 into a txt file into desingated organ file
        
        #Calcualte p-values comparing infected samples to 3w vs  0w and 12w vs 3w
        filtered_df = data[(data['ATTRIBUTE_organ'] == organ)] #filter original dataset by organ
        
        infected_df_12 = filtered_df[(filtered_df["ATTRIBUTE_infected_parasites"]=='infected') & (filtered_df['ATTRIBUTE_endpoint']=='12w')] #filter infected and 12w
        #infected_df_0 = filtered_df[(filtered_df["ATTRIBUTE_infected_parasites"]=='infected') & (filtered_df['ATTRIBUTE_endpoint']=='0w')] #filter infected and 0w
        infected_df_3 = filtered_df[(filtered_df["ATTRIBUTE_infected_parasites"]=='infected') & (filtered_df['ATTRIBUTE_endpoint']=='3w')] #filter infected and 3w
    
        p_values_12_3={}
        #p_values_3_0={}
        for column in data.select_dtypes(include=[np.number]).columns:
            infected_vals_12=infected_df_12[column].dropna() #remove NAs for infected
            #infected_vals_0=infected_df_0[column].dropna()
            infected_vals_3=infected_df_3[column].dropna()
            
            stats, p_val_12_3=mannwhitneyu(infected_vals_12, infected_vals_3) #p-values comparing infected of 12 w and 3 w of each metabolite in each organ using Mann Whitney test
            p_values_12_3[column]=p_val_12_3 #save to dictionary
            
            #stats, p_val_3_0=mannwhitneyu(infected_vals_3, infected_vals_0) #p-values comparing infected of 24 w and 12 w of each metabolite in each and organ
            #p_values_3_0[column]=p_val_3_0 #save to dictionary
            
            
        pval_df_12_3 = pd.DataFrame.from_dict(p_values_12_3, orient='index', columns=['adj.P.Val'])
        result_df_12_3 = ratio12_3.join(pval_df_12_3) #join p values of 12v3w with logFC of each organ
        result_df_12_3.to_csv(filename12v3, sep='\t')
        
        #pval_df_3_0 = pd.DataFrame.from_dict(p_values_3_0, orient='index', columns=['adj.P.Val'])
        #result_df_3_0 = ratio3_0.join(pval_df_3_0) #join p values of 24v12w with logFC of each organ
        #result_df_3_0.to_csv(filename3v0, sep='\t')
    
    
        
        
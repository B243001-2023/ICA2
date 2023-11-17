#packages for activate the shell command
import subprocess
import os

subprocess.call('export PATH=${HOME}/edirect:${PATH}', shell = True)
print("setting up the environment of edirect")

while 1:
    try:
        while 1:
            #get the input taxonomy from the user
            taxon = ''
            while 1:
                #loop for get a legal taxonomy name
                taxon = input('please pick a taxonomy\n')
                if taxon != ''  and '\n' not in taxon:
                    #testing the taxon
                    print('transforming the taxon into id (please do input anything now)')
                    taxid = subprocess.getoutput('esearch -db taxonomy -query "' + taxon + '"|efetch')
                    #print(taxid)
                if taxid != '' and ('ERROR' not in taxid) and ('FAILURE' not in taxid):
                    break
                print('please input a legal taxonomy name!')
                
            taxid = taxid.split('\n')
            taxid = taxid[::2]
            #get the taxid from the efetch format
            break

        #print(taxid)
        print('Here are the candidates of taxonomy')
        for i in range(0, len(taxid)):
            print(taxid[i])
        while 1:
            #loop for get a legal index of the taxids
            i = input('pick the taxonomy you want (using the number before the name)\n')
            i = int(i)
            if i in list(range(1,len(taxid)+1)):
                #print(taxid[i-1][3:]) used it to check the output
                taxid = subprocess.getoutput('esearch -db taxonomy -query "' + taxid[i-1][3:] + '"|efetch -format taxid')
                print('got the taxid, which is ' + str(taxid) + '\n')
                break
            else:
                print('please input a legal index')
        
        #init the fasta file, make sure there is no information inside
        subprocess.run('echo -n "" > a.fasta', shell = True)
        
        while 1:
            family = input('please pick a protein family\n')
            keywords = '"complete ' + family+' AND txid'+ taxid + '[ORGN] NOT (predicted OR hypothetical OR LOW QUALITY)"'
            #the keywords will be the most important part of the 
            if subprocess.getoutput('esearch -db protein -query ' + keywords + '|efetch -format fasta'):
                break
            else:
                print("You didn't pick a legal name of protein family, please correct it!")
        print('your searching command is:')
        print('esearch -db protein -query ' + keywords + '|efetch -format fasta > a.fasta')
        print('searching the sequences from the protein database\n')
        subprocess.run('esearch -db protein -query ' + keywords + '|efetch -format fasta >> a.fasta', shell=True)
        temp = open('a.fasta').read()
        #if we can get a fasta with something inside, it means we have sucessfully finished the esearch step
        if temp:
            break
    except:
        print('There is something WRONG in your input or network, please do this again/n')


fasta = open('a.fasta')
#read the fasta file into the loop, preparing for next step

seq_fine_names=[]
seq_dict = dict()
seq_name = ''
#init some variable for the following loop

for line in fasta.readlines():
    if (line[0] == '>'):
        seq_name = line[0:-1]
        seq_dict[seq_name] = ''
    else:
        seq_dict[seq_name] += line[0:-1]

    if (line[0] == '>') and (('LOW' in line) or ('PREDICTED' in line) or (family not in line)):
        seq_fine_names += [line]

print('\nAll the ' + str(len(seq_dict)) +'seqs which are collected')
print('Here is an example')
print(list(seq_dict.items())[0])

#import the packages
import numpy as np
import pandas as pd

#%% function list

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import toytree
import toyplot.png

#hmm sequence comparison
def needle(self, k):
    for i in list(self.keys()):
        print(i)
    #picking a seq as the standard one
    standard_seq = ''
    while 1:
        standard_seq = input('pick a seq as standard seq (just copy and paste a name you want)\n')
        if standard_seq in list(self.keys()):
            break
        else:
            print('check your input! you need copy and paste the seq name without any extra char')
    #picking the candidate seqs
    candidate_seq = []
    while 1:
        candidate_seq_1 = input('pick a seq as candidate, input 0 to finish\n')
        print(str(len(candidate_seq)) + ' seqs have been picked')
        if candidate_seq_1 in list(self.keys()):
            candidate_seq += [candidate_seq_1]
        elif candidate_seq_1 == '0' and len(candidate_seq) > 0:
            print('Your candidate_seqs are:')
            print(candidate_seq)
            break
        else:
            print('Please check your input, input a legal seq name or at least get a seq!')
    open('standard.fasta', 'w').write(standard_seq+'\n'+self[standard_seq]+'\n')
    for i in candidate_seq:
        open('candidate.fasta', 'w').write(standard_seq+'\n'+self[standard_seq]+'\n')
    subprocess.run('needle -asequence standard.fasta -bsequence candidate.fasta -outfile needle_'+k+'.output', shell = True)
    print('needle_'+str(k)+'.output has been generated')
    return 0

#using the clustalo analysis
def clustalo_analysis():
    try:
        if os.path.exists('tree.png'):
            os.remove('tree.png') #init environment
        subprocess.run('clustalo -i a.fasta -o a.out --outfmt=phy -v --force --guidetree-out=b', shell = True)
        tre = open('b').read()
        tre = tre.replace('\n','')
        tre = toytree.tree(tre)
        tre_plot = tre.draw(width=1600, height=1600)
        toyplot.png.render(tre_plot[0], 'tree.png')
        image_path = 'tree.png'
        img = mpimg.imread(image_path)
        plt.imshow(img)
        plt.axis('off')
        plt.show()
        return(0)
    except:
        return(1)

def plotcon():
    try:
        subprocess.run('plotcon -sequence a.fasta -scorefile EBLOSUM62 -winsize 160 -graph png', shell = True)
        img = mpimg.imread('plotcon.1.png')
        plt.imshow(img)
        plt.show()
        return(0)
    except:
        return(1)

def patmatmotif(self):
    keys_list = list(self.keys())
    values_list = list(self.values())
    print('anything')
    os.makedirs('seqs', exist_ok=True)
    os.makedirs('motif_report', exist_ok=True)
    for i in range(len(self)):
        a=str(i)
        print(keys_list[i])
        open('seqs/'+a+'.fasta', 'w').write(keys_list[i]+'\n')
        open('seqs/'+a+'.fasta', 'a').write(values_list[i])
        subprocess.run('patmatmotifs -sequence seqs/'+a+'.fasta -outfile motif_report/'+a+'.report -full', shell = True)
    return(0)

#main function area
i = 1
while i!=0:
    print('What are you gonna do?')
    print('0: quit')
    print('1: needle, picking a standard seq and some candidate seqs for a alignment')
    print('2: plotcon, the similarity of seqs depends on the position')
    print('3: patmatmotif, which will find the motif in the seqs')
    print('4: clustalo_analysis, which will generate a tree file and a all to all alignment result')

    i = int(input('\nPick a function, program will do the analysis for you.\n'))
    k=1
    if i == 1:
        k += needle(seq_dict,k)
    elif i == 2:
        a = plotcon()
        if a==1:
            print('sorry, something wrong with this function')
    elif i == 3:
        patmatmotif(seq_dict)
    elif i == 4:
        clustalo_analysis()
    elif i == 0:
        print('Thanks for using me')
        break
    else:
        print('please input a legal number of function')





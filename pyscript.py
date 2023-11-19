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
                    check = input('do you really want use this taxonomy? input y to confirm, other keys to cancel: ')
                    if check == 'y':
                        break
                    else:
                        print('going back')
                        continue
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
                taxid = subprocess.getoutput('esearch -db taxonomy -query "' + taxon + '"|efetch -format taxid')
                taxid = taxid.split('\n')[i-1]
                print('got the taxid, which is ' + str(taxid))
                check = input('do you really want use this taxid? input y to confirm, other keys to cancel: ')
                if check == 'y':
                    break
                else:
                    print('going back')
                    continue
            #found a code breaker
            else:
                print('please input a legal index')
        
        #init the fasta file, make sure there is no information inside
        subprocess.run('echo -n "" > a.fasta', shell = True)
        
        while 1:
            family = input('please pick a protein family\n')
            keywords = '"complete ' + family+' AND txid'+ taxid + '[ORGN] NOT (predicted OR hypothetical OR LOW QUALITY)"'
            #the keywords will be the most important part of the search command line
            result = subprocess.getoutput('esearch -db protein -query ' + keywords + '|efetch -format fasta')
            if 'FAILURE' not in result and 'ERROR' not in result and result:
                check = input('do you really want use this protein family? input y to confirm, other keys to cancel: ')
                if check == 'y':
                    break
                else:
                    print('going back')
                    continue
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

if len(seq_dict) > 1000:
    i = input('your protein family have too many seqs, do you wanna abandon some? (1000 seqs will be kept) [y/n]')
    if i == 'y':
        open('a.fasta', 'w').write('')
        for i in range(0, 1000):
            open('a.fasta', 'a').write(list(seq_dict.keys())[i] + '\n')
            open('a.fasta', 'a').write(list(seq_dict.values())[i] + '\n')
        seq_dict = dict(list(seq_dict.items())[0:1000])
        print('only 1000 seqs have been kept')
    else:
        print('nothing changed, as you wish')

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
        print(str(len(candidate_seq)+1) + ' seqs have been picked')
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
    subprocess.run('needle -asequence standard.fasta -bsequence candidate.fasta -outfile needle_'+str(k)+'.output', shell = True)
    print('needle_'+str(k)+'.output has been generated')
    return('needle_'+str(k)+'.output')

#using the clustalo analysis
def clustalo_analysis():
    try:
        if os.path.exists('tree.png'):
            os.remove('tree.png') #init environment
        subprocess.run('clustalo -i a.fasta -o Clustalo.out --outfmt=phy -v --force --guidetree-out=b', shell = True)
        tre = open('b').read()
        tre = tre.replace('\n','')
        tre = toytree.tree(tre)
        tre_plot = tre.draw()
        toyplot.png.render(tre_plot[0], 'tree.png')
        image_path = 'tree.png'
        img = mpimg.imread(image_path)
        plt.imshow(img)
        plt.axis('off')
        plt.show()
        return(['tree.png', 'Clustalo.out'])
    except:
        return(False)

def plotcon():
    try:
        subprocess.run('plotcon -sequence a.fasta -scorefile EBLOSUM62 -winsize 160 -graph png', shell = True)
        img = mpimg.imread('plotcon.1.png')
        plt.imshow(img)
        plt.show()
        return('plotcon.1.png')
    except:
        return(1)

import shutil
import re
def patmatmotif(self):
    if os.path.exists('motif_report'):
        shutil.rmtree('motif_report')
        #init the environment
    keys_list = list(self.keys())
    values_list = list(self.values())
    print('anything')
    os.makedirs('seqs', exist_ok=True)
    os.makedirs('motif_report', exist_ok=True)
    for i in range(len(self)):
        #make a output for every single file
        a=str(i)
        print(keys_list[i])
        open('seqs/'+a+'.fasta', 'w').write(keys_list[i]+'\n')
        open('seqs/'+a+'.fasta', 'a').write(values_list[i])
        subprocess.run('patmatmotifs -sequence seqs/'+a+'.fasta -outfile motif_report/'+a+'.report -full', shell = True)
    
    open('Motifs_report', 'w').write('')

    for i in range(len(self)):
        report = open('motif_report/'+a+'.report', 'r').read()
        seq = re.search('Sequence:.*', report).group().split(' ')[1]
        open('Motifs_report', 'a').write(seq)
        pattern = re.compile(r'Motif = .*')
        matches = pattern.finditer(report)
        motifs = []
        for match in matches:
            motifs += [match.group().split(' = ')[1]]
        open('Motifs_report', 'a').write('\t'+str(len(motifs)))
        for motif in motifs:
            open('Motifs_report', 'a').write('\t'+motif)
        open('Motifs_report', 'a').write('\n')
    return('Motifs_report')

def cons():
    subprocess.run('cons -sequence a.fasta -outseq cons_output', shell = True)
    return('cons_output')

def garnier():
    subprocess.run('garnier -sequence a.fasta -outfile garnier_output', shell=True)
    return('garnier_output')

def pepstats():
    subprocess.run('pepstats -sequence a.fasta -outfile pepstats_output', shell=True)
    return('pepstats_output')

#%%main function area
output_files = set()
output_files.add(plotcon())
output_files.add(patmatmotif(seq_dict))
for i in clustalo_analysis():
    output_files.add(i)

i = 1
while i!=0:
    print('What are you gonna do?')
    print('0: quit')
    print('1: plotcon (already finished automatically), the similarity of seqs depends on the position')
    print('2: patmatmotif (already finished automatically), which will find the motif in the seqs')
    print('3: clustalo_analysis (already finished automatically), which will generate a tree file and an all to all alignment result')
    print('4: needle, picking a standard seq and some candidate seqs for a alignment')
    print('5: cons, which can generate a consensus alignment file')
    print('6: garnier, predicting protein secondary structure using GOR method')
    print('7: pepstats, generate a file containing all physical information of seqs')

    i = int(input('\nPick a function, program will do the analysis for you.\n'))
    k=1
    if i == 4:
        output_files.add(needle(seq_dict,k))
        k+=1
    elif i == 1:
        a = plotcon()
        if a==False:
            print('sorry, something wrong with this function')
        else:
            output_files.add(a)
    elif i == 2:
        output_files.add(patmatmotif(seq_dict))
    elif i == 3:
        for i in clustalo_analysis():
            output_files.add(i)
    elif i == 5:
        output_files.add(cons())
    elif i == 6:
        output_files.add(garnier())
    elif i == 7:
        output_files.add(pepstats())
    elif i == 0:
        if os.path.exists('RESULTS'):
            shutil.rmtree('RESULTS')
        os.makedirs('RESULTS')
        for i in output_files:
            print(i)
            os.system('cp '+i+' RESULTS')
        print('Thanks for using me')
        break
    else:
        print('please input a legal number of function')





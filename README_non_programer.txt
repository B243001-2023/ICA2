Hello, here is the README file of the ICA2 python file.

At first, the ICA2 python was named as pyscript.py. After you unziped the file, you will see it in the directory of the unzip path.

'''
First step:
'''
#The python version is python 3.8.10.
To use the python file, please type the python path of your server and add the term 'pyscript.py' with a space.
#Here is an example:
python3 pyscript.py

'''
Second step:
'''
After you booted the pyfile sucessfully, the program will ask you to input some information of your analysis target protein family

#Here is an example (in the [] is the commit, when you running the program, they should NOT appear):

please pick a taxonomy
aves	[here is your first input, pick a taxonomy name]
transforming the taxon into id (please do input anything now)
do you really want use this taxonomy? input y to confirm, other keys to cancel: y [you inputed the y to confirm]
Here are the candidates of taxonomy
1. Aves
pick the taxonomy you want (using the number before the name)
1	[picked the index 1, which is the Aves]
got the taxid, which is 8782	[got the id sucessfully]
do you really want use this taxid? input y to confirm, other keys to cancel: y	[you input the y to confirm again]
please pick a protein family
glucose-6-phosphatase [this is the family you picked]
do you really want use this protein family? input y to confirm, other keys to cancel: y [you input the y to confirm again]

#after you finished the step 1 and 2, the program will do the basic analysis for you, if you want do them again, you can do it in the folowing step
#NOTICE: the picture windows need you to close them manufactually. Or the program will stop until you close it.


'''
Third step:
'''
#the program will ask you to finish some analysis in the achievable functions, which will make the output and collect them automaticly.
#Here is the example

What are you gonna do?
0: quit
1: plotcon (already finished automaticly), the similarity of seqs depends on the position
2: patmatmotif (already finished automaticly), which will find the motif in the seqs
3: clustalo_analysis (already finished automaticly), which will generate a tree file and an all to all alignment result
4: needle, picking a standard seq and some candidate seqs for a alignment
5: cons, which can generate a consensus alignment file
6: garnier, predicting protein secondary structure using GOR method
7: pepstats, generate a file containing all physical information of seqs

Pick a function, program will do the analysis for you.
0	[this is your input, you picked the function 0, so the program will quit and you finished your analysis]

[the different function will not affect each other, they can be done sparately, you can only pick one function or some funcitons as you want]

#All the output of different functions will be collected into the directory named RESULT
#NOTICE do please use the 0 function to quit the program, or you wil get nothing in the RESULT file.

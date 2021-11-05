#!/usr/bin/python3

protein = "MSRSLLLRFLLFLLLLPPLP"
aa = "R"
aa_count = protein.count(aa)
protein_length = len(protein)
percentage = aa_count * 100 / protein_length
print(percentage)

def get_aa_percentage(protein, aa):
    aa_count = protein.upper().count(aa.upper())
    protein_length = len(protein)
    percentage = (aa_count / protein_length) * 100
    return percentage

get_aa_percentage( "MSRSLLLRFLLFLLLLPPLP", "R")
get_aa_percentage( "MSRSLLLRFLLFLLLLPPLP", "r")

assert get_aa_percentage("MSRSLLLRFLLFLLLLPPLP", "M") == 5
assert get_aa_percentage("MSRSLLLRFLLFLLLLPPLP", "r") == 10
assert get_aa_percentage("msrslllrfllfllllpplp", "L") == 50
assert get_aa_percentage("MSRSLLLRFLLFLLLLPPLP", "Y") == 0

protein = "MSRSLLLRFLLFLLLLPPLP"
aa_list = ['M', 'L', 'F']

counter = 0

for aa in aa_list: 
    print("counting number of " + aa)
    aa_count = protein.upper().count(aa.upper())
    counter = counter + aa_count 
    print("running total is " + str(counter))
    percentage = (counter / len(protein)) * 100
    print("final percentage is " + str(percentage))

def get_aa_percentage_v2(protein, aa_list): 
    protein_length = len(protein) 
    counter = 0 
    for aa in aa_list: 
        aa_count = protein.upper().count(aa.upper()) 
        counter = counter + aa_count 
    percentage = (counter / protein_length) * 100
    return percentage 

get_aa_percentage_v2("MSRSLLLRFLLFLLLLPPLP", ['M', 'L', 'F'])
get_aa_percentage_v2("MSRSLLLRFLLFLLLLPPLP", 'mlf')

assert get_aa_percentage_v2("MSRSLLLRFLLFLLLLPPLP", ["M"]) == 5
assert get_aa_percentage_v2("MSRSLLLRFLLFLLLLPPLP", ['F', 'S', 'L']) == 70

def get_aa_percentage_v2(protein, aa_list='aILMFWYV'): 
    protein_length = len(protein) 
    counter = 0 
    for aa in aa_list: 
        aa_count = protein.upper().count(aa.upper()) 
        counter = counter + aa_count 
    percentage = (counter / protein_length) * 100
    return percentage

assert get_aa_percentage_v2("MSRSLLLRFLLFLLLLPPLP", ["M"]) == 5
assert get_aa_percentage_v2("MSRSLLLRFLLFLLLLPPLP", ['F', 'S', 'L']) == 70
assert get_aa_percentage_v2("MSRSLLLRFLLFLLLLPPLP",'FSL') == 70
assert get_aa_percentage_v2("MSRSLLLRFLLFLLLLPPLP") == 65

get_aa_percentage_v2("MSRSLLLRFLLFLLLLPPLP","l1234")
get_aa_percentage_v2("MSRSLLLRFLLFLLLLPPLP","lab1234")
get_aa_percentage_v2("MSRSLLLRFLLFLLLLPPLP",set(["l","1","2","l","3","4","l"]))

def count_undetermined_1(dna):
    total_good_bases = 0
    for base in ['A', 'T', 'G', 'C']:
        total_good_bases = total_good_bases + dna.upper().count(base)
    return len(dna) - total_good_bases

count_undetermined_1('atucgtgractanctgactg')

def count_undetermined_2(dna):
    total_undetermined = 0
    for base in dna.upper():
        if base not in ['A', 'T', 'G', 'C']:
            total_undetermined = total_undetermined + 1
    return total_undetermined

count_undetermined_2('atucgtgractanctgactg')

def count_undetermined_3(dna):
    total_undetermined = len(dna.upper().replace('A',"").replace('T',"").replace('G',"").replace('C',""))
    return total_undetermined

count_undetermined_3('atucgtgractanctgactg')

def count_undetermined_4(dna):
    total_undetermined = len(dna.upper().replace('A',"").replace('T',"").replace('G',"").replace('C',""))
    prop_undetermined = total_undetermined / len(dna)
    return prop_undetermined

count_undetermined_4('atucgtgractanctgactg')

def count_undetermined_5(dna,threshold=0.1):
    total_undetermined = len(dna.upper().replace('A',"").replace('T',"").replace('G',"").replace('C',""))
    prop_undetermined = total_undetermined / len(dna)
    return prop_undetermined >= threshold

count_undetermined_5('atucgtgractanctgactg', 0.2)


# Write a FUNCTION that, given any DNA sequence, will print
# all the k-mers (of a chosen size e.g. 4-mers) that occur more than some
# chosen number of times.
# Written by s2258074 on 05 Nov 2021
def find_my_kmers(dna,ksize=2,minkfreq=3) :
   if ksize > len(dna) :
     return "Sorry, your kmer length is longer than your DNA (" + str(len(dna)) +" bases)." 
   if ksize < 2 or ksize > 50 :
     return "Sorry, inappropriate kmer length, only 2 to 50 accepted here."
   print("Processing sequence of length",len(dna),"for kmers longer than",ksize,"and frequency greater than",minkfreq)
   kmersfound = []
   kmerstarts = list(range(0,len(dna)))
   for base in kmerstarts:
       if (base+ksize) < len(dna)+1:
           seqout = (dna)[base:base+ksize]
           kmersfound = kmersfound + [seqout] 
   nrset = list(set(kmersfound)) 

   returnstuff = []
   for kfreqfind in nrset:
       if kmersfound.count(kfreqfind) > minkfreq :
           returnstuff.append(kfreqfind.upper()+": "+str(kmersfound.count(kfreqfind)))
   return print(returnstuff)

find_my_kmers("ttagatcctgaacgtgaacgcacggatttagatcctgaacgtgaacgcacggat",2,2)
find_my_kmers("ttagatcctgaacgtgaacgcacggatttagatcctgaacgtgaacgcacggat")
find_my_kmers("ttagatcctgaacgtgaacgcacggatttagatcctgaacgtgaacgcacggat",55)
find_my_kmers("ttagatcctgaacgtgaacgcacggatttagatcctgaacgtgaacgcacggat",51)
find_my_kmers("ttagatcctgaacgtgaacgcacggatttagatcctgaacgtgaacgcacggat",50)

# Write a Python script/programme that, given any DNA sequence, will print
# all the k-mers (of a chosen size e.g. 4-mers) that occur more than some
# chosen number of times.
# User interaction to supply all details


import os
def new_find_my_kmers(dna,ksize=2,minkfreq=2) :
   kmersfound = []
   kmerstarts = list(range(0,len(dna)))
   for base in kmerstarts:
       if (base+ksize) < len(dna)+1:
           seqout = (dna)[base:base+ksize]
           kmersfound = kmersfound + [seqout] 
   nrset = list(set(kmersfound))
   returnstuff = []
   for kfreqfind in nrset:
       if kmersfound.count(kfreqfind) > minkfreq :
           returnstuff.append(kfreqfind.upper()+": "+str(kmersfound.count(kfreqfind)))
   return returnstuff

inputdna = input ("What is your sequence?\n\t").upper()
if inputdna :
  inputksize = int ( input ("What kmer size shall I use?\n\t") or 2)
  if (inputksize < 2 or inputksize >= len(inputdna) or inputksize > 50) :
     inputksize = 2
     print("Inappropriate value chosen, resetting to 2\n\t")
  inputminkfreq = int ( input ("What minimum frequency shall I use?\n\t") or 2)
  print("Thanks!  Processing:\n"+inputdna+"\n for a kmersize of "
   +str(inputksize)+",\n reporting frequencies greater than "
   +str(inputminkfreq)+"\n")
  outputstuff = new_find_my_kmers(dna=inputdna,ksize=inputksize,minkfreq=inputminkfreq)
  outputstuff.sort()
  myfilename="kmerouts"+"_KMER"+str(inputksize)+"_MIN"+str(inputminkfreq)+".txt"
  outputfilepipe = open(myfilename,"w")

  if len(outputstuff) == 0 :
    print("No kmers met the criteria, so no outputs to file!\n")
    outputfilepipe.close()
  else:
    print("Results:\n")
    print(outputstuff)


    outputfilepipe.write("### Kmer analysis\n#SQ "+str(inputdna)+
     "\n#KMER "+str(inputksize)+
     "\n#MIN "+str(inputminkfreq)+"\n")
    for freqseq in outputstuff :
        outputfilepipe.write(freqseq+"\n")
    outputfilepipe.write("\n")

    outputfilepipe.close()
    print("\n\nContents of the output:\n")
    syscmd="cat " + myfilename
    os.system(syscmd)

else :
  print ("Sorry, really can\'t do any of this without a sequence!\n")


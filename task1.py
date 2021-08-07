# task1 compare CpG islands
# read fasta file contained chromosome 20 1:1000000 as chr20
f=open('chr20_part.fasta')
chr20=''
f.readline() # ignore the first line
for line in f:
    chr20 += line.upper().rstrip()
f.close()

# define the starting position for finding CpG islands: skip the N bases at the beginning of the sequence
start = min(chr20.find('A'),chr20.find('T'),chr20.find('C'),chr20.find('G')) - 50

# finding CpG island
start_pos = 0 # default the CpG islands' start position
CpG={} # default the dictionary storing the coordinates of CpG islands
for i in range(start, len(chr20)-100+1):
    seq = chr20[i:i+100]
    no_C = seq.count('C')
    no_G = seq.count('G')
    if no_C == 0 or no_G == 0:
        start_pos = 0
    if no_C != 0 and no_G != 0: 
        GC_per =(no_C + no_G)/100
        Obs_Exp_CpG = 100 * seq.count('CG')/no_C/no_G
        if GC_per > 0.5 and Obs_Exp_CpG > 0.6:
            if start_pos == 0:
                start_pos = i
            if i+100-start_pos+1 >= 200:
                end_pos = i+100-1
                CpG[start_pos] = end_pos
        else:
            start_pos = 0

# create pandas dataframe from CpG dictionary and cnvert it to bed format
import pandas as pd
import pybedtools
CpG_ser = pd.Series(CpG) # CpG series
CpG_df = pd.DataFrame({'name': 'chr20','start_pos':CpG_ser.index, 'end_pos':CpG_ser.values}) # CpG dataframe
CpG_bf = pybedtools.BedTool.from_dataframe(CpG_df) # CpG bed format
print (CpG_bf) # in command line save to CpG_py.bed
import time
import pandas as pd
import numpy as np
import pyBigWig as bw
import matplotlib.pyplot as plt 
import random
import tqdm 

start_time = time.time()

# peaks from bed files
bed = pd.read_csv("/storage2/asirix/S1/WG_on_ ATAC/ATAC_track.bed", sep="\t", header=None, usecols=[0,1,2] ,
                  names=["chrom","st","end"])

# center of each peak
assert np.all(bed.end - bed.st > 0)
bed["mids"] = bed.st + (bed.end - bed.st) // 2

# set of distances
distances = np.concatenate([
np.arange(-3000,-1000,20), 
np.arange(-1000,1000),
np.arange(1000,3000,20)
])
num_shuffle = 100

# compute coverage
b = bw.open("/storage2/asirix/S1/WG_on_ ATAC/Offset1_K562WG.coverage.bw")

# check consistency bigwig and bed-file

# 1. check overlap of chromosomes

bed_chroms_all = bed.chrom.unique()
bed_chroms_clear = [element for element in bed_chroms_all if len(element)<=5 and element != 'chrM' and element != 'chrY' and element != 'chrX' ]

bigWig_chroms_dict = b.chroms()
bigWig_chroms_all = bigWig_chroms_dict.keys()
bigWig_chroms_clear=[element for element in bigWig_chroms_all if len(element)<=5 and element != 'chrM' and element != 'chrY' and element != 'chrX']

assert len(bed_chroms_clear)==len(bigWig_chroms_clear) 

assert len(set(bed_chroms_clear + bigWig_chroms_clear)) == len(bed_chroms_clear)

# 2. check chrm sizes
for chrm in bed_chroms_clear:
    bed_length = bed.query("chrom == @chrm")["end"].max()
    assert bed_length <= bigWig_chroms_dict[chrm]

coverages = []

#creating_array_with_intervals_values_for_each_chrom
chrom_arrays = {}
for chrom in bigWig_chroms_clear:
  chrom_arrays[chrom] = np.zeros(shape=b.chroms()[chrom], dtype=np.float32)
for chrom in tqdm.tqdm(bigWig_chroms_clear):
  intervals = b.intervals(chrom)
  for start,end,value in intervals:
    chrom_arrays[chrom][start:end] = value

def compute_average_coverage(row, dist, bigwigfile, radius=1): 
        pos = row["mids"] + dist

        if pos > 0 and pos < bigWig_chroms_dict[row["chrom"]]:

          coverage = sum(chrom_arrays[row["chrom"]][pos:pos+radius])

          if coverage is not None and np.isfinite(coverage):
            return coverage
          else:
            coverage=0
            return coverage 
        else: # what happens is this 'if' condition is not TRUE?
          print(pos)
          print(bigWig_chroms_dict[row["chrom"]])
          raise Exception("wrong pos for row "+str(row) + " and distance "+str(dist) )

bed = bed.query("chrom in @bigWig_chroms_clear")
for dist in tqdm.tqdm(distances):
  coverage = bed.apply(compute_average_coverage, dist=dist, bigwigfile = b, axis="columns")
    # join 
  coverages.append(coverage.values.mean())
save_values=pd.DataFrame({'distances': distances, 'values': coverages})
save_values.to_csv('/storage2/asirix/S1/WG_on_ ATAC/values.txt', sep='\t')

plt.plot(distances, coverages)

# shuffle

def calc_random_control(bed=bed, bigWig_chroms_clear=bigWig_chroms_clear):
 chroms = bigWig_chroms_clear 
 random_control = pd.DataFrame(columns=["chrom", "st", "end"])
 for chr in chroms:  
    starts = bed.query('chrom==@chr')["st"].values
    ends = bed.query('chrom==@chr')["end"].values
     
    region_lengths = ends - starts
    inter_region_lengths = starts[1:] - ends[:-1] # well done!
    assert inter_region_lengths.min() > 0 # good point!

    random.shuffle(region_lengths)
    random.shuffle(inter_region_lengths)
    
    new_starts = np.insert(np.cumsum(region_lengths[:-1] + inter_region_lengths)+starts.min(), 0, starts.min())

    # new_starts.sort() # actually, should it be sorted alread?
    # i suggest to comment this line and see whether it will
    # pass through the assert below

    assert (new_starts[1:] - new_starts[:-1]).min() >= 0
    
    new_ends = new_starts +  region_lengths

    ###
    chroms = [chr] * len(starts)

    assert len(set(chroms))==1
 

    random_for_chrom = pd.DataFrame({'chrom': chroms, 'st': new_starts, 'end': new_ends})
    random_control = pd.concat([random_control, random_for_chrom], axis=0, ignore_index=True)
 random_control["mids"] = random_control.st + (random_control.end - random_control.st) // 2

    # this should be better done once, outside the loop
    # you don't whant to recompute the mids each time you add new chromsome
    # just a small step toward fight with the global Earth warming =)
    #random_control["mids"] = random_control.st + (random_control.end - random_control.st) // 2

 return random_control



coverages_many_cont=np.zeros(len(distances))
coverages_cont_var=[]

#shuffling several times
for i in tqdm.tqdm(range(num_shuffle)): # why not for i in range(num_shuffle) ?
 coverages_cont_aver=[]

 control=calc_random_control(bed=bed, bigWig_chroms_clear=bigWig_chroms_clear)

 for dist in distances:
  coverage_cont = control.apply(compute_average_coverage, dist=dist, bigwigfile = b, axis="columns")
  coverages_cont_aver.append(coverage_cont.mean())
 

 coverages_cont_var.append(coverages_cont_aver) 
 coverages_many_cont=coverages_many_cont+coverages_cont_aver

coverages_cont_var=np.reshape(coverages_cont_var, (num_shuffle,len(distances)))

coverages_many_cont_var=np.var(coverages_cont_var, axis=0)
coverages_many_cont=coverages_many_cont/num_shuffle

random_values=pd.DataFrame({'distances': distances, 'mean': coverages_many_cont, 'var': coverages_many_cont_var})
random_values.to_csv('/storage2/asirix/S1/WG_on_ ATAC/random_values.txt', sep='\t')

coverages_3sigma_up=coverages_many_cont+(coverages_many_cont_var**0.5)*3
coverages_3sigma_down=coverages_many_cont-(coverages_many_cont_var**0.5)*3

plt.plot(distances, coverages_many_cont, color='r')
plt.plot(distances, coverages_3sigma_up, '--', color='r')
plt.plot(distances, coverages_3sigma_down,'--', color='r')
plt.fill_between(distances,coverages_3sigma_up,coverages_3sigma_down,color='r',alpha=0.2)
plt.savefig('/storage2/asirix/S1/WG_on_ ATAC/WG_on_ATAC.png')
print("--- %s seconds ---" % (time.time() - start_time))
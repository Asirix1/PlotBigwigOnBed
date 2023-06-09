import time
import pandas as pd
import numpy as np
import pyBigWig as bw
import matplotlib.pyplot as plt 
import random
import tqdm 
import os
import pybedtools as pb

start_time = time.time()

i=os.environ['i'] #bed name 
g=os.environ['g'] #bigwig name

# peaks from bed files
bed = pd.read_csv('/mnt/storage/home/aapopov/S1_final_step_5/S1_final/tracks/'+i, sep="\t", header=None, usecols=[0,1,2],
                  names=["chrom","st","end"])

blacklist= pd.read_csv('/mnt/storage/home/aapopov/S1_final_step_5/S1_final/hg38-blacklist.v2.bed',sep='\t',header=None, usecols=[0,1,2], 
                       names=["chrom","st","end"])

# set of distances
distances = np.arange(-3000,3000,2)
num_shuffle = 100

# compute coverage
b = bw.open('/mnt/storage/home/aapopov/S1_final_step_5/S1_final/bigwigs/' + g)

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

#creating array with intervals values for each chrom
chrom_arrays = {}
for chrom in bigWig_chroms_clear:
  chrom_arrays[chrom] = np.zeros(shape=b.chroms()[chrom], dtype=np.float32)
for chrom in tqdm.tqdm(bigWig_chroms_clear):
  intervals = b.intervals(chrom)
  for start,end,value in intervals:
    chrom_arrays[chrom][start:end] = value

# center of each peak and deleting blacklisted regions

assert np.all(bed.end - bed.st >= 0)
bed["mids"] = bed.st + (bed.end - bed.st) // 2
bed=pb.BedTool.from_dataframe(bed)
blacklist=pb.BedTool.from_dataframe(blacklist)
bed=bed.intersect(blacklist, v=True)
bed = bed.to_dataframe()
bed.rename(columns={'start':'st','name': 'mids'}, inplace=True)

#calculations
coverage=[]

for chrom in tqdm.tqdm(bigWig_chroms_clear):
 bed_chrom = bed.query("chrom == @chrom")
 mids_array=np.array(bed_chrom ["mids"])

 mids_array=mids_array.reshape(len(mids_array),-1)
 coordinates=mids_array+distances
 assert np.all(coordinates) > 0
 assert np.all(coordinates) < bigWig_chroms_dict[chrom]
 chrom_coverage=[chrom_arrays[chrom][item] for item in coordinates]
 assert np.array(chrom_coverage).size == coordinates.size
 coverage.extend(chrom_coverage)
 
coverage=np.array(coverage)
coverage.reshape(-1,len(distances))
coverage_mean=np.mean(coverage,axis=0)

save_values=pd.DataFrame({'distances': distances, 'values': coverage_mean})
save_values.to_csv('/mnt/storage/home/aapopov/S1_final_step_5/S1_final/values/'+g+'_on_'+i+'.txt', sep='\t')

plt.plot(distances, coverage_mean, lw=0.5)

# shuffle

def calc_random_control(bed=bed, bigWig_chroms_clear=bigWig_chroms_clear):
 chroms = bigWig_chroms_clear 
 random_control = pd.DataFrame(columns=["chrom", "st", "end"])
 for chr in chroms:  
    starts = bed.query('chrom==@chr')["st"].values
    ends = bed.query('chrom==@chr')["end"].values
     
    region_lengths = ends - starts
    inter_region_lengths = starts[1:] - ends[:-1] 
    assert inter_region_lengths.min() > 0 

    random.shuffle(region_lengths)
    random.shuffle(inter_region_lengths)
    
    new_starts = np.insert(np.cumsum(region_lengths[:-1] + inter_region_lengths)+starts.min(), 0, starts.min())


    assert (new_starts[1:] - new_starts[:-1]).min() >= 0
    
    new_ends = new_starts +  region_lengths

    chroms = [chr] * len(starts)

    assert len(set(chroms))==1
 

    random_for_chrom = pd.DataFrame({'chrom': chroms, 'st': new_starts, 'end': new_ends})
    random_control = pd.concat([random_control, random_for_chrom], axis=0, ignore_index=True)
 random_control["mids"] = random_control.st + (random_control.end - random_control.st) // 2


 return random_control

coverages_cont_aver=[]

#shuffling several times
for t in tqdm.tqdm(range(num_shuffle)): 
 random_coverage=[]
 control=calc_random_control(bed=bed, bigWig_chroms_clear=bigWig_chroms_clear)
 control=pb.BedTool.from_dataframe(control)
 control=control.intersect(blacklist, v=True)
 control = control.to_dataframe()
 control.rename(columns={'start':'st','name': 'mids'}, inplace=True)

 for chrom in bigWig_chroms_clear:
  control_chrom=control.query("chrom == @chrom")
  mids_array=np.array(control_chrom["mids"],dtype=np.int32)
  mids_array=mids_array.reshape(len(mids_array),-1)
  coordinates=mids_array+distances
  assert np.all(coordinates) > 0
  assert np.all(coordinates) < bigWig_chroms_dict[chrom]
  random_chrom_coverage=[chrom_arrays[chrom][item] for item in coordinates]
  assert np.array(random_chrom_coverage).size == coordinates.size
  random_coverage.extend(random_chrom_coverage)
 
 random_coverage=np.array(random_coverage)
 random_coverage.reshape(-1,len(distances))
 random_coverage_mean=np.mean(random_coverage,axis=0)
 coverages_cont_aver.append(random_coverage_mean)
control_mean=np.average(coverages_cont_aver,axis=0)
control_std=np.std(coverages_cont_aver,axis=0)

random_values=pd.DataFrame({'distances': distances, 'mean': control_mean, 'std': control_std})
random_values.to_csv('/mnt/storage/home/aapopov/S1_final_step_5/S1_final/random_values/'+g+'_on_'+i+'.txt', sep='\t')

plt.plot(distances, control_mean, color='r',lw=0.5)
plt.plot(distances, control_mean+3*control_std, '--', color='r',lw=0.5)
plt.plot(distances, control_mean-3*control_std,'--', color='r',lw=0.5)
plt.fill_between(distances,control_mean+3*control_std,control_mean-3*control_std,color='r',alpha=0.2)
plt.xlabel("Distance to center (bp)")
plt.ylabel("Average signal")
plt.axvline(linestyle ='--',color='black')
for pos in ['right', 'top', 'bottom', 'left']:
    plt.gca().spines[pos].set_visible(False)
plt.savefig('/mnt/storage/home/aapopov/S1_final_step_5/S1_final/figures/'+g+'_on_'+i+'.png')
print("--- %s seconds ---" % (time.time() - start_time))

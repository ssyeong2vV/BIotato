import pandas as pd
import os
import sys
#---------------------------------------------------------------------------------------#
path = os.getcwd()
BATCH = {}
with open(f"{path}/batchconfig.txt", "r") as read :
  for line in read:
    line = line.strip('\n').strip().split("=")
    key = line[0]
    value = line[1]
    BATCH[key] = value
#---------------------------------------------------------------------------------------#
Sample = pd.read_csv("SampleSheet.txt", sep='\t', header=None)
total_sample_list = list(Sample.loc[:,0])
total_sample_path_list= list(Sample.loc[:,1])
Sample_size = list(Sample.loc[:,3])
Sample_size = list(map(lambda x: str(x).replace('MB',''), Sample_size))
Sample_size = list(map(lambda x: str(x).strip(), Sample_size))
sorted_size = sorted(Sample_size, reverse = True)
idx = list(map(lambda x: Sample_size.index(x), sorted_size))
#---------------------------------------------------------------------------------------#
for n in range(Sample.shape[0]):
    folder_name = str(Sample.loc[n,0]).strip()
    if os.path.isdir(folder_name):
        pass
    else :
        os.mkdir(f"{path}/{folder_name}")
#---------------------------------------------------------------------------------------#
if BATCH['Run_type'] == 'WGBS':
    Code = '/labmed/00.Code/BI.Study/00.Pipeline/WGBS.py'
elif BATCH['Run_type'] == 'RNA':
    Code = '/labmed/00.Code/BI.Study/00.Pipeline/RNA.py'
elif BATCH['Run_type'] == 'TARGET':
    Code = '/labmed/01.ALL/03.python/lsy/00.CODE/GATK_Pipeline.py'
elif BATCH['Run_type'] == 'WGS':
    Code = '/labmed/01.ALL/03.python/lsy/00.CODE/GATK_Pipeline.py'
#---------------------------------------------------------------------------------------#
if BATCH['Node'] == 'node01' and int(BATCH['CPU']) > 64:
        raise ValueError("\033[91mValueError: Total CPU is less than 128\033[0m")
elif BATCH['Node'] == 'node02' and int(BATCH['CPU']) > 28:
        raise ValueError("\033[91mValueError: Total CPU is less than 56\033[0m")
elif BATCH['Node'] == 'node03' and int(BATCH['CPU']) > 16:
        raise ValueError("\033[91mValueError: Total CPU is less than 32\033[0m")
elif BATCH['Node'] == 'node04' and int(BATCH['CPU']) > 28:
        raise ValueError("\033[91mValueError: Total CPU is less than 28\033[0m")
#---------------------------------------------------------------------------------------#
Cpu = int(BATCH['CPU'])
Allocated_CPU = int(Cpu / Sample.shape[0])
if Allocated_CPU < 1:
        raise ValueError("\033[91m" + "ValueError: Allocated CPU is less than 1" + "\033[0m")

if Sample.shape[0] == 1 :
    CPU = [Allocated_CPU] * Sample.shape[0]
    if BATCH['Node'] != 'node04' :
        CPU = list(map(lambda x: x*2, CPU))
elif Sample.shape[0] > 1 :
        CPU = [Allocated_CPU] * Sample.shape[0]
        extra_cpu = int(Cpu % Sample.shape[0])
        i = 0
        for i in range(extra_cpu) :
            CPU[idx[i]]+=1
            i+= 1
        if BATCH['Node'] != 'node04' :
            CPU = list(map(lambda x: x*2, CPU))
#---------------------------------------------------------------------------------------#
for n in range(Sample.shape[0]):
    folder_name = str(Sample.loc[n,0]).strip()
    Cpu = CPU[n]
    if os.path.isdir(folder_name):
        os.chdir(folder_name)
        with open("Run.sh", "w") as note:
          note.write("#!/bin/bash" + '\n' +
                    f"#SBATCH -J " + BATCH['Run_type'] + '.' + folder_name + "\n" +
                    f"#SBATCH -o Log.%j.out" + '\n' +
                    f"#SBATCH --time=UNLIMITED" + '\n' +
                    f"#SBATCH --nodelist={BATCH['Node']}" +'\n' +
                    f"#SBATCH -n {Cpu}" + '\n' +
                    '\n' +
                    f"python3 {Code}" + '\n')

        with open(f"{folder_name}_batchconfig.txt", "w") as write_batch_config:
          for keys in BATCH.keys():
            if keys != 'CPU' :
                write_batch_config.write(f"{keys}={BATCH[keys]}\n")
            elif keys == 'CPU' :
                write_batch_config.write(f"{keys}={Cpu}\n")
          write_batch_config.write(f"sample={total_sample_list}\n")
          write_batch_config.write(f"sample_dir={total_sample_path_list}\n")

        os.chdir("..")

#---------------------------------------------------------------------------------------#
with open("Total_run.sh", "w") as note:
    for n in range(Sample.shape[0]):
        folder_name = str(Sample.loc[n,0]).strip()
        note.write(f"cd {path}/{folder_name}; sbatch Run.sh" + '\n')
#---------------------------------------------------------------------------------------#


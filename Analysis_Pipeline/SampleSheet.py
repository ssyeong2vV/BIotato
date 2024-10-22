#!/home/lsy/anaconda3/envs/NGS/bin/python

import sys
import os
import datetime

LIST = sys.argv[1::2]

nowDate = datetime.datetime.now()
filename = nowDate.strftime("%Y-%m-%d_%H%M_") + "Datalist.txt"

with open('SampleSheet.txt', 'w') as note1:
    with open(filename, 'w') as datefile:
        datefile.write(nowDate.strftime("%Y-%m-%d %H:%M")+"\n")
        for sample in LIST:
            if '_R1' in sample:
                Name = sample.split('/')[-1].split('_R1')[0]
                if not os.path.isdir(f'{Name}'):
                    os.mkdir(f'{Name}')
                Size = os.path.getsize(sample) / (1024 ** 2)
                read1 = sample
                read2 = read1.replace('_R1', '_R2')
                note1.write(f'{Name}\t{read1}\t{read2}\t{Size: .2f} MB\n')

                datefile.write(f"{read1}\n{read2}\n")

                with open(f'{Name}/SampleSheet.txt', 'w') as note2:
                    r1 = sample
                    r2 = sample.replace('_R1', '_R2')
                    note2.write(f'{Name}\t{r1}\t{r2}\t')

            elif '_1' in sample:
                Name = sample.split('/')[-1].split('_1')[0]
                if not os.path.isdir(f'{Name}'):
                    os.mkdir(f'{Name}')
                Size = os.path.getsize(sample) / (1024 ** 2)
                read1 = sample
                read2 = read1.replace('_1', '_2')
                note1.write(f'{Name}\t{read1}\t{read2}\t{Size: .2f} MB\n')

                datefile.write(f"{read1}\n{read2}\n")

                with open(f'{Name}/SampleSheet.txt', 'w') as note2:
                    r1 = sample
                    r2 = sample.replace('_1', '_2')
                    note2.write(f'{Name}\t{r1}\t{r2}\t')

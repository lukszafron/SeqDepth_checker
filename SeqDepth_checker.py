#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Łukasz Szafron
@email: lukszafron@gmail.com
"""
appname = "SeqDepth_checker"
app_version = "1.0"

import pandas as pd
import re, os, sys, getopt, multiprocessing, subprocess, functools
from io import StringIO

def usage():
    print(
        "\nWelcome to the " + appname + " app.\n\n"
        "The following options are available:\n\n"
        "\t-b, --bedfile: A BED file with the gene regions the enrichment of which should be evaluated.\n"
        "\t-d, --data_dir: A path to a folder with BAM files to be tested.\n"
        "\t-o, --output_dir: A path to a folder where the resultant Excel file is saved (default: the data_dir).\n"
        "\t-t, --threads: Number of CPU threads to be used (default: 1).\n"
        "\t-s, --bam_suffix: A character string separating a sample name from the file extension (.bam).\n"
        "\t-S, --sd_times: The number of standard devations to be subtracted from each mean sequencing read coverage depth value to get a minimal read coverage depth value (default: 0).\n"
        "\t-T, --threshold: If the minimal sequencing read coverage depth value for a gene region is lower than the given threshold, the enrichment of this region is considered as failed (default: 5).\n"
        "\t-h, --help: prints this help message.\n"
        "\t-v, --version: prints the version of this program.\n"
        )
# The next two lines are for debugging purposes only and should be commented in the final program.
# option_list = ["-b", '/workspace/lukasz/NGS-all-in-one/BED_FILES/KAPA_HyperPETE_Hot_Spot_Panel_hg38_Capture_Targets.bed', "-d", '/workspace/lukasz/NGS-all-in-one/RUNS/TEST1/MAPPINGS_TRIMMED/GRCh38.chr.only.genome/HISAT2', "-o", '/workspace/lukasz/NGS-all-in-one/RUNS/TEST1/MAPPINGS_TRIMMED/GRCh38.chr.only.genome/HISAT2/COMPARISON', "-s", "_sorted.no_dups"]
# opts, args = getopt.getopt(option_list, "b:d:o:t:s:S:T:hv", ["bedfile=", "data_dir=", "output_dir=", "threads=","bam_suffix=", "sd_times=", "threshold=", "help","version"])
opts, args = getopt.getopt(sys.argv[1:], "b:d:o:t:s:S:T:hv", ["bedfile=", "data_dir=", "output_dir=", "threads=","bam_suffix=", "sd_times=", "threshold=", "help","version"])

try:
        opts,args
        if len(opts) == 0:
                usage()
                sys.exit()
        for o, a in opts:
            if o in ("-b", "--bedfile"):
                bedfile = a
            elif o in ("-d", "--data_dir"):
                directory = a
            elif o in ("-o", "--output_dir"):
                output_dir = a
            elif o in ("-t", "--threads"):
                threads = a
            elif o in ("-s", "--bam_suffix"):
                bam_suffix = a
            elif o in ("-S", "--sd_times"):
                sd_times = a
            elif o in ("-T", "--threshold"):
                threshold = a
            elif o in ("-h", "--help"):
                usage()
                sys.exit()
            elif o in ("-v", "--version"):
                print(appname, " version: ", app_version, sep = "")
                sys.exit()
            else:
                assert False, "Unhandled option: "+o

except getopt.GetoptError as err:
    # print help information and exit:
    print("\n"+str(err)) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)
    
try:
    bedfile
except:
    raise(FileExistsError("The BED file was not provided."))
try:
    directory
except:
    raise(Exception("The path to a directory with NGS BAM files was not provided.")) 
try:
    output_dir
except:
    output_dir = directory
try:
    threads
except:
    threads = 1
finally:
    threads = int(threads)
try:
    bam_suffix
except:
    raise(Exception("The BAM file suffix was not provided."))
try:
    sd_times
except:
    sd_times = 0
finally:
    sd_times = int(sd_times)
try:
    threshold
except:
    threshold = 5
finally:
    threshold = int(threshold)

with open(file = bedfile, mode = "rt") as bed:
    bed_data = pd.read_csv(bed, sep='\t', header=None)
if bed_data.shape[1] > 4:
    bed_data = bed_data.loc[:,0:3]
elif bed_data.shape[1] < 4:
    bed_data[3] = "NA"

bed_data.columns = ["Chromosome","Start","End","Description"]
bed_data_list = bed_data.to_dict(orient="records")
rownames = ["_".join((v.get("Chromosome"), str(v.get("Start")), str(v.get("End")), v.get("Description"))).split(";")[0] for v in bed_data_list]
regions = [v.get("Chromosome") + ":" + str(v.get("Start")) + "-" + str(v.get("End")) for v in bed_data_list]

bam_pattern = ''.join((bam_suffix, ".bam"))
bamfiles = [f for f in os.listdir(directory) if re.search(string=f, pattern="^.*"+bam_pattern+"$", flags=re.I)]
samples = [re.sub(string=s, pattern=bam_pattern, repl="", flags=re.I) for s in bamfiles]
samples.sort()
snames = [(index,sname) for index,sname in enumerate(samples)]

class SAMPLE:
    def __init__(self, sname):
        print("Processing sample {} ({} out of {})...\n".format(sname[1], sname[0]+1, len(snames)), flush=True)
        self.sname = sname[1]
        self.bamfile = os.path.join(directory,''.join((sname[1], bam_pattern)))
    def depth_stats(self, region):
        self.res = subprocess.run(["samtools", "depth", "-r", region, self.bamfile], capture_output= True, check= True, text = True)
        if len(self.res.stdout) > 0:
            depth_df = pd.read_csv(StringIO(self.res.stdout), sep="\t", header=None)
            depth_mean,depth_sd = (depth_df.loc[:,2].mean(), depth_df.loc[:,2].std())
            return((depth_mean, depth_sd))
        else:
            return((None, None))
    def exec_stats(self):
       self.means_sds = [self.depth_stats(region) for region in regions]
       self.df_stats = pd.DataFrame(self.means_sds, columns=(self.sname+"_means", self.sname+"_sds"))
       self.df_min = pd.DataFrame(self.df_stats.iloc[:,0] - (self.df_stats.iloc[:,1]*sd_times))
       self.df_min = self.df_min.rename(columns={0:self.sname+"_min"})
       
def exec_func(sname):
    sample = SAMPLE(sname)
    sample.exec_stats()
    return([sample.df_stats, sample.df_min])

# protect the entry point
if __name__ == '__main__':    
    with multiprocessing.Pool(threads) as pool:
        all_stats = pool.map(exec_func, snames)

stats = [s for s,m in all_stats]
mins = [m for s,m in all_stats]

df_means_sds = functools.reduce(lambda x,y:pd.merge(x,y, left_index=True, right_index=True, how="inner", sort=False), stats)
df_means_sds.index = rownames
df_min = functools.reduce(lambda x,y:pd.merge(x,y, left_index=True, right_index=True, how="inner", sort=False), mins)
df_min.index = rownames
df_min_bool = (df_min < threshold) | (df_min.isna())
df_failed_indexes = df_min_bool.loc[df_min_bool.sum(axis=1) > 0].index
df_failed_samples = [', '.join(df_min_bool.loc[index].loc[df_min_bool.loc[index]].index).replace("_min", "") for index in df_failed_indexes]
df_failed_all = pd.DataFrame(df_failed_samples, df_failed_indexes).rename(columns = {0:"Samples"})
min_value = df_min.min(axis=1).to_frame().sort_values(by = 0).iloc[0].values[0]
min_name = df_min.min(axis=1).to_frame().sort_values(by = 0).iloc[0].name
min_sample = df_min.min(axis=0).to_frame().sort_values(by = 0).iloc[0].name.rstrip("_min")

excel_file_name = '_'.join(["Seq_depths", os.path.basename(bedfile), "BAMs", bam_suffix, "threshold", str(threshold), "SDs", str(sd_times), ".xlsx"])
with pd.ExcelWriter(os.path.join(output_dir, excel_file_name)) as writer:
    df_means_sds.to_excel(writer, sheet_name="Means_and_SD_values", index=True, na_rep="NA")
    df_min.to_excel(writer, sheet_name="Minimal_values", index=True, na_rep="NA")
    df_failed_all.to_excel(writer, sheet_name="Failed_regions_and_samples", index=True, na_rep="NA")

print("The analysis is complete.\n")
print("Analysis parameters: sequencing read coverage depth threshold: {}, minimal read coverage depth values were calculated by subtracting {} sd values from mean read coverage depth values.".format(threshold, sd_times))
print("The number of poorly enriched gene regions: {}/{}.".format(len(df_failed_indexes), len(regions)))
print("The most poorly enriched region, name: {}, sample: {}, value: {:.3f}.".format(min_name, min_sample, min_value))

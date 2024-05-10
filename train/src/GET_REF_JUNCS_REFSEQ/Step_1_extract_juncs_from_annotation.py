import os 
import re
import json 
import torch

project_root = '/ccb/cybertron/khchao/splam-analysis-results/' 
output_dir = f'{project_root}train/results/ALL_RefSeq/REF_junctions/'

def chr_name_convert():
    f_chrs = open(f'{project_root}/Dataset/Refseq_2_UCSU_chromosome_names.tsv', "r")
    lines = f_chrs.read().splitlines()
    chrs = {}
    for line in lines:
        eles = line.split("\t")
        chrs[eles[0]] = eles[1]
    return chrs

def main():
    JUNC_COUNTER = 0
    extracted_junc_type = "all"
    os.makedirs(output_dir, exist_ok=True)
    fw = open(f'{output_dir}ref_d_a_{extracted_junc_type}.bed', 'w')
    chrs = chr_name_convert()
    with open(f'{project_root}/Dataset/refseq_GCF_000001405.40_GRCh38.p14_genomic.gff', 'r') as f:
        lists = f.read().splitlines() 
        transcript_id = ""
        prev_transcript_id = ""
        chr = ""
        prev_chr = ""
        strand = "."
        prev_strand = "."
        starts = []
        ends = []
        is_protein_entries = False
        for ele in lists:
            line = ele.split("\t")
            if len(line) < 8:
                continue
            if (line[2] == "gene"):
                if extracted_junc_type == "all":
                    match_prot = re.search(r"gene_biotype=(protein_coding|lncRNA|ncRNA|nc_RNA)", line[8])
                else:
                    match_prot = re.search(r"gene_biotype=protein_coding", line[8])
                if match_prot is not None:
                    is_protein_entries = True
                else:
                    is_protein_entries = False
            if is_protein_entries:
                if (line[2] == "exon"):
                    match = re.search(r"transcript_id=\w+", line[8])
                    if match is not None:
                        transcript_id = match.group()[14:]
                        chr = line[0]
                        strand = line[6]
                        exon_start = int(line[3])
                        exon_end = int(line[4])

                        if prev_transcript_id != transcript_id:
                            if prev_strand == '-':
                                starts.reverse()    
                                ends.reverse()
                            starts = starts[1:]
                            ends = ends[:-1]

                            for idx in range(len(starts)):
                                JUNC_COUNTER += 1
                                fw.write(chrs[prev_chr] + "\t" + str(ends[idx]) + "\t" + str(starts[idx]) + "\t" + "JUNC" + "\t1\t" + prev_strand + "\n")                                
                            starts.clear()
                            ends.clear()
                        starts.append(exon_start)
                        ends.append(exon_end)
                        prev_transcript_id = transcript_id
                        prev_chr = chr
                        prev_strand = strand        

        ###########################
        # Process the last group!
        ###########################
        if prev_strand == '-':
            starts.reverse()    
            ends.reverse()        
        starts = starts[1:]
        ends = ends[:-1]
        for idx in range(len(starts)):
            JUNC_COUNTER += 1
            fw.write(chrs[prev_chr] + "\t" + str(ends[idx]) + "\t" + str(starts[idx]) + "\t" + "JUNC" + "\t1\t" + prev_strand + "\n")
        starts.clear()
        ends.clear()
    fw.close()

if __name__ == "__main__":
    main()
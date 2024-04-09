import pandas as pd
import os 

def get_hg38_chrom_size():
    f_chrs = open("./hg38_chrom_size.tsv", "r")
    lines = f_chrs.read().splitlines()
    chrs = {}
    for line in lines:
        eles = line.split("\t")
        chrs[eles[0]] = int(eles[1])
    return chrs

SEQ_LEN="800"
QUATER_SEQ_LEN = int(SEQ_LEN) // 4

project_root = '/ccb/cybertron/khchao/splam-analysis-results/' 
output_dir = f'{project_root}/train/results/tool_benchmark/splam/'
def main():
    chrs = get_hg38_chrom_size()
    chrs_test = ["chr1", "chr9"]
    #################################
    # For 'd_a.bed': 0-based, 1-based
    # For 'donor.bed': 0-based, 0-based
    # For 'acceptor.bed': 0-based, 0-based
    #################################
    os.makedirs(output_dir, exist_ok=True)
    fw_donor = open(f"{output_dir}/splam.juncs.donor.bed", "w")
    fw_acceptor = open(f"{output_dir}/splam.juncs.acceptor.bed", "w")
    d_a_bed = f"{output_dir}/d_a.sample1.bed"
    fw_da = open(d_a_bed, "w")
    JUNCS = set()
    with open(f'{output_dir}/splam.juncs.bed', "r") as f:
        lines = f.read().splitlines()
        for line in lines:
            eles = line.split("\t")
            if len(eles) < 6:
                continue
            chroms = eles[0]
            junc_name = eles[3]
            score = eles[4]
            strand = eles[5]
            if strand != "+" and strand != "-":
                continue
            if chroms in chrs_test:
                if (strand == "+"):
                    donor = int(eles[1])
                    acceptor = int(eles[2])
                    splice_junc_len = acceptor - donor
                elif (strand == "-"):
                    acceptor = int(eles[1])
                    donor = int(eles[2])
                    splice_junc_len = donor - acceptor
                flanking_size = QUATER_SEQ_LEN
                if splice_junc_len < QUATER_SEQ_LEN:
                    flanking_size = splice_junc_len
                if (strand == "+"):
                    donor_s = donor - QUATER_SEQ_LEN
                    donor_e = donor + flanking_size
                    acceptor_s = acceptor - flanking_size
                    acceptor_e = acceptor + QUATER_SEQ_LEN
                elif (strand == "-"):
                    donor_s = donor - flanking_size
                    donor_e = donor + QUATER_SEQ_LEN
                    acceptor_s = acceptor - QUATER_SEQ_LEN
                    acceptor_e = acceptor + flanking_size
                if donor_e >= chrs[chroms] or acceptor_e >= chrs[chroms]:
                    continue
                if donor_s < 0 or acceptor_s < 0:
                    continue
                new_junc = (chroms, str(donor_s), str(donor_e), str(acceptor_s), str(acceptor_e), strand)
                if new_junc in JUNCS:
                    continue
                else:
                    JUNCS.add(new_junc)
                    fw_donor.write(chroms + "\t" + str(donor_s) + "\t" + str(donor_e) + "\t" + junc_name+"_donor" + "\t" + score + "\t" + strand + "\n")
                    fw_acceptor.write(chroms + "\t" + str(acceptor_s) + "\t" + str(acceptor_e) + "\t" + junc_name+"_acceptor" + "\t" + score + "\t" + strand + "\n")
                    if (strand == "+"):
                        fw_da.write(chroms + "\t" + str(donor) + "\t" + str(acceptor+1) + "\tJUNC\t" + score + "\t" + strand + "\n")
                    elif (strand == "-"):
                        fw_da.write(chroms + "\t" + str(acceptor) + "\t" + str(donor+1) + "\tJUNC\t" + score + "\t" + strand + "\n")
    fw_donor.close()
    fw_acceptor.close()
    fw_da.close()

if __name__ == "__main__":
    main()
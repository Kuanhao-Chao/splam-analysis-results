import os

project_root = '/ccb/cybertron/khchao/splam-analysis-results/' 
output_dir = f'{project_root}/train/results/tool_benchmark/spliceai/'

def main():
    SEQ_LEN = "800"
    HALF_SEQ_LEN = int(SEQ_LEN) // 2
    QUATER_SEQ_LEN = int(SEQ_LEN) // 4
    print("QUATER_SEQ_LEN: ", QUATER_SEQ_LEN)
    output_files = [output_dir]
    for output_file in output_files:
        print(">> output_file")
        fr = open(output_file+"spliceai.juncs.seq.fa", "r")
        fw_noN = open(output_file+"spliceai.noN.juncs.seq.fa", "w")
        
        fw_N = open(output_file+"spliceai.N.juncs.seq.fa", "w")
        lines = fr.read().splitlines()
        line_num = len(lines)

        donors_noN = {}
        donors_N = {}
        acceptors_noN = {}
        acceptors_N = {}

        for idx in range(line_num):
            if idx % 2 == 0:
                chr_name = lines[idx]
                strand = lines[idx][-2]
            else:
                seq = lines[idx]
                length = len(seq)

                for method in ["noN", "N"]:
                    if method == "noN":
                        x = seq
                    else:
                        x = 'N'*(5000) + seq[5000:-5000] + 'N'*(5000)
                    x = x.upper()                    
                    eles = lines[idx-1].split(":")
                    chr = eles[0]
                    strand = eles[-1][-2]
                    s_e = eles[1].split("-")
                    start = str(int(s_e[0])+5200)
                    end = str(int(s_e[1].split("(")[0])-5200)
                    # [:-3]
                    # print(">> start - end: ", start, " - ", end)

                    if method == "noN":
                        fw_noN.write(chr+";"+start+";"+end+";"+strand+"\n")
                        fw_noN.write(x + "\n")
                    else:
                        fw_N.write(chr+";"+start+";"+end+";"+strand+"\n")
                        fw_N.write(x + "\n")
                    
                    donor_dimer = x[QUATER_SEQ_LEN + 5000:QUATER_SEQ_LEN+2 + 5000]
                    acceptor_dimer = x[len(x) - QUATER_SEQ_LEN-5000-2:len(x)-QUATER_SEQ_LEN-5000]

                    if method == "noN":
                        if donor_dimer not in donors_noN.keys():
                            donors_noN[donor_dimer] = 1
                        else:
                            donors_noN[donor_dimer] += 1
                        if acceptor_dimer not in acceptors_noN.keys():
                            acceptors_noN[acceptor_dimer] = 1
                        else:
                            acceptors_noN[acceptor_dimer] += 1
                    
                    else:
                        if donor_dimer not in donors_N.keys():
                            donors_N[donor_dimer] = 1
                        else:
                            donors_N[donor_dimer] += 1
                        if acceptor_dimer not in acceptors_N.keys():
                            acceptors_N[acceptor_dimer] = 1
                        else:
                            acceptors_N[acceptor_dimer] += 1
        with open(f'{output_dir}motif.txt', "w") as f:
            for key, value in donors_noN.items():
                f.write("Donor: "+key+" ("+str(value)+")\n")
                print("\tDonor   : ", key, " (", value, ")")
            for key, value in acceptors_noN.items():
                f.write("Acceptor: "+key+" ("+str(value)+")\n")
                print("\tAcceptor: ", key, " (", value, ")")
        fw_noN.close()
        fw_N.close()

if __name__ == "__main__":
    main()
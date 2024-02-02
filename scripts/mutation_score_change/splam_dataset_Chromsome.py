import torch
from torch.utils.data import Dataset, DataLoader, random_split
from pathlib import Path
import json
import random
import os
import math

from splam_utils import *
# MODEL_VERSION
SEQ_LEN = "800"

def split_seq_name(seq):
    return seq[1:]

class myDataset(Dataset):
    def __init__(self, type, segment_len=800, shuffle=True, input_seq_dir=None):
        print("!!shuffle: ", shuffle, input_seq_dir)
        self.segment_len = segment_len
        self.data = []
        #################################
        ## Processing 'POSITIVE' samples
        #################################
        pidx = 0
        with open(input_seq_dir, "r") as f:
            print("Processing ", input_seq_dir)
            lines = f.read().splitlines()
            seq_name = ""
            seq = ""
            for line in lines:
                # print(line)
                if pidx % 2 == 0:
                    seq_name = split_seq_name(line)
                elif pidx % 2 == 1:
                    seq = line
                    # print(seq)
                    X, Y = create_datapoints(seq, '+')
                    X = torch.Tensor(np.array(X))
                    Y = torch.Tensor(np.array(Y)[0])
                    if X.size()[0] != 800:
                        print("seq_name: ", seq_name)
                        print(X.size())
                        print(Y.size())
                    self.data.append([X, Y, seq_name])
                pidx += 1
                if pidx %10000 == 0:
                    print("pidx: ", pidx)
                if type == "train":
                    pass
                elif type == "test":
                    if pidx >= 3000:
                        break

            print("pidx: ", pidx)

            CONSTANT_SIZE = pidx
            CONSTANT_SIZE_NEG = CONSTANT_SIZE*3
            print("\033[1m[INFO] CONSTANT_SIZE     : ", CONSTANT_SIZE, "\033[0m")
            print("\033[1m[INFO] CONSTANT_SIZE_NEG : ", CONSTANT_SIZE_NEG, "\033[0m")

        # if type == "train" or type == "test" or (type == "eval" and eval_select=="neg_1"):
        #     #################################
        #     ## Processing 'NEGATIVE_1' samples
        #     #################################
        #     n1idx = 0
        #     with open(neg_1_f, "r") as f:
        #         print("Processing ", neg_1_f)
        #         lines = f.read().splitlines()
        #         seq_name = ""
        #         seq = ""
        #         for line in lines:
        #             # print(line)
        #             if n1idx % 2 == 0:
        #                 seq_name = split_seq_name(line)
        #             elif n1idx % 2 == 1:
        #                 seq = line
        #                 X, Y = create_datapoints(seq, '-')
        #                 X = torch.Tensor(np.array(X))
        #                 Y = torch.Tensor(np.array(Y)[0])
        #                 if X.size()[0] != 800:
        #                     print("seq_name: ", seq_name)
        #                     print(X.size())
        #                     print(Y.size())
        #                 self.data.append([X, Y, seq_name])
        #             n1idx += 1
        #             if n1idx %10000 == 0:
        #                 print("n1idx: ", n1idx)

        #             if type == "train" or type == "test":
        #                 if n1idx >= CONSTANT_SIZE_NEG:
        #                     break
        #     print("n1idx: ", n1idx)

        # if type == "train" or type == "test" or (type == "eval" and eval_select=="neg_random"):
        #     #################################
        #     ## Processing 'NEGATIVE_1' samples
        #     #################################
        #     nridx = 0
        #     with open(neg_random_f, "r") as f:
        #         print("Processing ", neg_random_f)
        #         lines = f.read().splitlines()
        #         seq_name = ""
        #         seq = ""
        #         for line in lines:
        #             # print(line)
        #             if nridx % 2 == 0:
        #                 seq_name = split_seq_name(line)
        #             elif nridx % 2 == 1:
        #                 seq = line
        #                 # print(seq)
        #                 X, Y = create_datapoints(seq, '-')
        #                 X = torch.Tensor(np.array(X))
        #                 Y = torch.Tensor(np.array(Y)[0])
        #                 if X.size()[0] != 800:
        #                     print("seq_name: ", seq_name)
        #                     print(X.size())
        #                     print(Y.size())
        #                 self.data.append([X, Y, seq_name])
        #             nridx += 1
        #             if nridx %10000 == 0:
        #                 print("nridx: ", nridx)

        #             if type == "train" or type == "test":
        #                 if nridx >= CONSTANT_SIZE_NEG:
        #                     break
        #     print("nridx: ", nridx)

        # if type == "eval" and eval_select=="pos_ALTS":
        #     #################################
        #     ## Processing 'POSITIVE_REFSEQ_PROTEIN_ISOFORMS' samples
        #     #################################
        #     pp_alts_idx = 0
        #     with open(pos_ALTS_f, "r") as f:
        #         print("Processing ", pos_ALTS_f)
        #         lines = f.read().splitlines()
        #         seq_name = ""
        #         seq = ""
        #         for line in lines:
        #             # print(line)
        #             if pp_alts_idx % 2 == 0:
        #                 seq_name = split_seq_name(line)
        #             elif pp_alts_idx % 2 == 1:
        #                 seq = line
        #                 # print(seq)
        #                 X, Y = create_datapoints(seq, '+')
        #                 X = torch.Tensor(np.array(X))
        #                 Y = torch.Tensor(np.array(Y)[0])
        #                 # print("X.size(): ", X)
        #                 # print("Y.size(): ", Y)
        #                 if X.size()[0] != 800:
        #                     print("seq_name: ", seq_name)
        #                     print(X.size())
        #                     print(Y.size())
        #                 self.data.append([X, Y, seq_name])
        #             pp_alts_idx += 1
        #             if pp_alts_idx %10000 == 0:
        #                 print("pidx: ", pp_alts_idx)
        #                 # print(seq_name)
        #     print("pp_alts_idx: ", pp_alts_idx)

        # if type == "eval" and eval_select=="pos_MANE":
        #     #################################
        #     ## Processing 'POSITIVE_REFSEQ_PROTEIN_ISOFORMS' samples
        #     #################################
        #     pp_MANE_idx = 0
        #     with open(pos_MANE_f, "r") as f:
        #         print("Processing ", pos_MANE_f)
        #         lines = f.read().splitlines()
        #         seq_name = ""
        #         seq = ""
        #         for line in lines:
        #             # print(line)
        #             if pp_MANE_idx % 2 == 0:
        #                 seq_name = split_seq_name(line)
        #             elif pp_MANE_idx % 2 == 1:
        #                 seq = line
        #                 # print(seq)
        #                 X, Y = create_datapoints(seq, '+')
        #                 X = torch.Tensor(np.array(X))
        #                 Y = torch.Tensor(np.array(Y)[0])
        #                 # print("X.size(): ", X)
        #                 # print("Y.size(): ", Y)
        #                 if X.size()[0] != 800:
        #                     print("seq_name: ", seq_name)
        #                     print(X.size())
        #                     print(Y.size())
        #                 self.data.append([X, Y, seq_name])
        #             pp_MANE_idx += 1
        #             if pp_MANE_idx %10000 == 0:
        #                 print("pidx: ", pp_MANE_idx)
        #     print("pp_MANE_idx: ", pp_MANE_idx)
        
        #################################
        ## Shuffle the data 
        #################################
        if shuffle: random.shuffle(self.data)

    def __len__(self):
        return len(self.data)
 
    def __getitem__(self, index):
        # Load preprocessed mel-spectrogram.
        feature = self.data[index][0]
        label = self.data[index][1]
        seq_name = self.data[index][2]
        feature = torch.flatten(feature, start_dim=1)
        return feature, label, seq_name

def get_eval_dataloader(batch_size, TARGET, n_workers, shuffle, input_seq_dir):
    #######################################
    # predicting splice / non-splice
    #######################################
    print("get_eval_dataloader shuffle: ", shuffle)
    testset = myDataset("eval", int(SEQ_LEN), shuffle, input_seq_dir)
    test_loader = DataLoader(
        testset,
        batch_size = batch_size,
        shuffle = shuffle,
        drop_last = False,
        pin_memory = True,
    )
    print("[INFO] Loading dataset (shuffle: " + str(shuffle) + "): ", "./INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/test.pt")
    # torch.save(test_loader, "../src_tools_evaluation/splam_result/splam_dataloader.pt")
    return test_loader


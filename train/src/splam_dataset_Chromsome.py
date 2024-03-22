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

project_root = '/ccb/cybertron/khchao/splam-analysis-results/'
def split_seq_name(seq):
    return seq[1:]


class myDatasetEval(Dataset):
    def __init__(self, type, segment_len=800, shuffle=True, eval_select=None, idx=0):
        print("!!shuffle: ", shuffle, eval_select)
        self.segment_len = segment_len
        self.data = []
        neg_f = f'{project_root}train/results/Negs/Neg_{idx}/INPUTS/800bp/input_neg_{idx}.fa'
        CONSTANT_SIZE_NEG = 10000
        if type == "eval":
            #################################
            ## Processing 'NEGATIVE_1' samples
            #################################
            nidx = 0
            with open(neg_f, "r") as f:
                print("Processing ", neg_f)
                lines = f.read().splitlines()
                seq_name = ""
                seq = ""
                for line in lines:
                    if nidx % 2 == 0:
                        seq_name = split_seq_name(line)
                    elif nidx % 2 == 1:
                        seq = line
                        X, Y = create_datapoints(seq, '-')
                        X = torch.Tensor(np.array(X))
                        Y = torch.Tensor(np.array(Y)[0])
                        if X.size()[0] != 800:
                            print("seq_name: ", seq_name)
                            print(X.size())
                            print(Y.size())
                        self.data.append([X, Y, seq_name])
                    nidx += 1
                    if nidx %10000 == 0:
                        print("nidx: ", nidx)
                    if nidx >= CONSTANT_SIZE_NEG:
                        break
            print("nidx: ", nidx)
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


class myDatasetTrain(Dataset):
    def __init__(self, type, segment_len=800, shuffle=True, eval_select=None, idx=0):
        print("!!shuffle: ", shuffle, eval_select)
        self.segment_len = segment_len
        self.data = []
        pos_f = ""
        if type == "train" or type == "test":
            pos_f = f'{project_root}/train/results/MANE/INPUTS/800bp/input_pos.fa'
            neg_1_f = f'{project_root}/train/results/Neg_1/INPUTS/800bp/input_pos.fa'
            neg_random_f = f'{project_root}/train/results/Neg_1/INPUTS/800bp/input_pos.fa'
        # elif type == "eval":
        #     pos_MANE_f = "./INPUTS/"+SEQ_LEN+"bp/input_pos_MANE/test_pos_MANE.shuffle.fa"
        #     pos_ALTS_f = "./INPUTS/"+SEQ_LEN+"bp/input_pos_ALTS/test_pos_ALTS.shuffle.fa"
        #     neg_1_f = "./INPUTS/"+SEQ_LEN+"bp/input_neg_1/test_neg_1.shuffle.fa"
        #     neg_random_f = "./INPUTS/"+SEQ_LEN+"bp/input_neg_random/test_neg_random.shuffle.fa"
        #     neg_f = f'/ccb/cybertron/khchao/splam-analysis-results/train/results/Negs/Neg_{idx}/INPUTS/800bp/input_neg_{idx}.fa'

        CONSTANT_SIZE_NEG = 10000
        
        #################################
        ## Processing 'POSITIVE' samples
        #################################
        if type == "train" or type == "test" or (type == "eval" and eval_select=="pos"):
            pidx = 0
            with open(pos_f, "r") as f:
                print("Processing ", pos_f)
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

        if type == "train" or type == "test" or (type == "eval" and eval_select=="neg_1"):
            #################################
            ## Processing 'NEGATIVE_1' samples
            #################################
            n1idx = 0
            with open(neg_1_f, "r") as f:
                print("Processing ", neg_1_f)
                lines = f.read().splitlines()
                seq_name = ""
                seq = ""
                for line in lines:
                    # print(line)
                    if n1idx % 2 == 0:
                        seq_name = split_seq_name(line)
                    elif n1idx % 2 == 1:
                        seq = line
                        X, Y = create_datapoints(seq, '-')
                        X = torch.Tensor(np.array(X))
                        Y = torch.Tensor(np.array(Y)[0])
                        if X.size()[0] != 800:
                            print("seq_name: ", seq_name)
                            print(X.size())
                            print(Y.size())
                        self.data.append([X, Y, seq_name])
                    n1idx += 1
                    if n1idx %10000 == 0:
                        print("n1idx: ", n1idx)

                    if type == "train" or type == "test":
                        if n1idx >= CONSTANT_SIZE_NEG:
                            break
            print("n1idx: ", n1idx)

        if type == "train" or type == "test" or (type == "eval" and eval_select=="neg_random"):
            #################################
            ## Processing 'NEGATIVE_1' samples
            #################################
            nridx = 0
            with open(neg_random_f, "r") as f:
                print("Processing ", neg_random_f)
                lines = f.read().splitlines()
                seq_name = ""
                seq = ""
                for line in lines:
                    # print(line)
                    if nridx % 2 == 0:
                        seq_name = split_seq_name(line)
                    elif nridx % 2 == 1:
                        seq = line
                        # print(seq)
                        X, Y = create_datapoints(seq, '-')
                        X = torch.Tensor(np.array(X))
                        Y = torch.Tensor(np.array(Y)[0])
                        if X.size()[0] != 800:
                            print("seq_name: ", seq_name)
                            print(X.size())
                            print(Y.size())
                        self.data.append([X, Y, seq_name])
                    nridx += 1
                    if nridx %10000 == 0:
                        print("nridx: ", nridx)

                    if type == "train" or type == "test":
                        if nridx >= CONSTANT_SIZE_NEG:
                            break
            print("nridx: ", nridx)

        if type == "eval" and eval_select=="pos_ALTS":
            #################################
            ## Processing 'POSITIVE_REFSEQ_PROTEIN_ISOFORMS' samples
            #################################
            pp_alts_idx = 0
            with open(pos_ALTS_f, "r") as f:
                print("Processing ", pos_ALTS_f)
                lines = f.read().splitlines()
                seq_name = ""
                seq = ""
                for line in lines:
                    # print(line)
                    if pp_alts_idx % 2 == 0:
                        seq_name = split_seq_name(line)
                    elif pp_alts_idx % 2 == 1:
                        seq = line
                        # print(seq)
                        X, Y = create_datapoints(seq, '+')
                        X = torch.Tensor(np.array(X))
                        Y = torch.Tensor(np.array(Y)[0])
                        # print("X.size(): ", X)
                        # print("Y.size(): ", Y)
                        if X.size()[0] != 800:
                            print("seq_name: ", seq_name)
                            print(X.size())
                            print(Y.size())
                        self.data.append([X, Y, seq_name])
                    pp_alts_idx += 1
                    if pp_alts_idx %10000 == 0:
                        print("pidx: ", pp_alts_idx)
                        # print(seq_name)
            print("pp_alts_idx: ", pp_alts_idx)

        if type == "eval" and eval_select=="pos_MANE":
            #################################
            ## Processing 'POSITIVE_REFSEQ_PROTEIN_ISOFORMS' samples
            #################################
            pp_MANE_idx = 0
            with open(pos_MANE_f, "r") as f:
                print("Processing ", pos_MANE_f)
                lines = f.read().splitlines()
                seq_name = ""
                seq = ""
                for line in lines:
                    # print(line)
                    if pp_MANE_idx % 2 == 0:
                        seq_name = split_seq_name(line)
                    elif pp_MANE_idx % 2 == 1:
                        seq = line
                        # print(seq)
                        X, Y = create_datapoints(seq, '+')
                        X = torch.Tensor(np.array(X))
                        Y = torch.Tensor(np.array(Y)[0])
                        # print("X.size(): ", X)
                        # print("Y.size(): ", Y)
                        if X.size()[0] != 800:
                            print("seq_name: ", seq_name)
                            print(X.size())
                            print(Y.size())
                        self.data.append([X, Y, seq_name])
                    pp_MANE_idx += 1
                    if pp_MANE_idx %10000 == 0:
                        print("pidx: ", pp_MANE_idx)
            print("pp_MANE_idx: ", pp_MANE_idx)
        
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


def get_train_dataloader(batch_size, TARGET, n_workers):
    """Generate dataloader"""
    # trainset_origin = myDataset("train", int(SEQ_LEN))
    # trainset, valset = torch.utils.data.random_split(trainset_origin, [0.9, 0.1])
    trainset = myDatasetTrain("train", int(SEQ_LEN))
    testset = myDatasetTrain("test", int(SEQ_LEN))
    train_loader = DataLoader(
        trainset,
        batch_size=batch_size,
        shuffle=True,
        drop_last=False,
        pin_memory=True,
    )
    # val_loader = DataLoader(
    #     valset,
    #     batch_size=batch_size,
    #     shuffle=True,
    #     drop_last=False,
    #     pin_memory=True,
    # )
    test_loader = DataLoader(
        testset,
        batch_size = batch_size,
        drop_last=False,
        pin_memory=True,
    )

    #######################################
    # predicting splice / non-splice
    #######################################
    return train_loader, test_loader
    # torch.save(train_loader, "./INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/train.pt")
    # torch.save(val_loader, "./INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/val.pt")
    # torch.save(test_loader, "./INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/test.pt")


def get_test_dataloader(batch_size, TARGET, n_workers, shuffle):
    #######################################
    # predicting splice / non-splice
    #######################################
    testset = myDatasetEval("test", int(SEQ_LEN), shuffle)
    test_loader = DataLoader(
        testset,
        batch_size = batch_size,
        shuffle = shuffle,
        drop_last = False,
        pin_memory = True,
    )
    print("[INFO] Loading dataset (shuffle: " + str(shuffle) + "): ", "./INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/test.pt")
    # test_loader = torch.load("./INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/test.pt")
    return test_loader

def get_eval_dataloader(batch_size, TARGET, n_workers, shuffle, eval_select, idx):
    #######################################
    # predicting splice / non-splice
    #######################################
    testset = myDatasetEval("eval", int(SEQ_LEN), shuffle, eval_select, idx)
    test_loader = DataLoader(
        testset,
        batch_size = batch_size,
        shuffle = shuffle,
        drop_last = False,
        pin_memory = True,
    )
    print(f'[INFO] Loading dataset (shuffle: /ccb/cybertron/khchao/splam-analysis-results/src_tools_evaluation/splam_result/{TARGET}/neg_{idx}/neg_{idx}.splam_dataloader.pt')
    torch.save(test_loader, f'/ccb/cybertron/khchao/splam-analysis-results/src_tools_evaluation/splam_result/{TARGET}/neg_{idx}/neg_{idx}.splam_dataloader.pt')
    return test_loader


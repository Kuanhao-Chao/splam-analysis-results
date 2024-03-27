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
    def __init__(self, type, segment_len=800, shuffle=True, eval_select=None, test_f=""):
        print("!!shuffle: ", shuffle, eval_select)
        self.segment_len = segment_len
        self.data = []
        CONSTANT_SIZE_NEG = 10000
        if type == "eval":
            #################################
            ## Processing 'NEGATIVE_1' samples
            #################################
            nidx = 0
            with open(test_f, "r") as f:
                print("Processing ", test_f)
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
    def __init__(self, process_type, segment_len=800, shuffle=True, eval_select=None, idx=0):
        self.segment_len = segment_len
        self.data = []
        if process_type == "train":
            pos_MANE_f = f'{project_root}/train/results/train_test_dataset/input_pos_mane/train_pos_mane.shuffle.fa'
            pos_ALTS_f = f'{project_root}/train/results/train_test_dataset/input_pos_alts/train_pos_alts.shuffle.fa'
            neg_1_f = f'{project_root}/train/results/train_test_dataset/input_neg_1/train_neg_1.shuffle.fa'
            neg_random_f = f'{project_root}/train/results/train_test_dataset/input_neg_random/train_neg_random.shuffle.fa'
        elif process_type == "test":
            pos_MANE_f = f'{project_root}/train/results/train_test_dataset/input_pos_mane/test_pos_mane.shuffle.fa'
            pos_ALTS_f = f'{project_root}/train/results/train_test_dataset/input_pos_alts/test_pos_alts.shuffle.fa'
            neg_1_f = f'{project_root}/train/results/train_test_dataset/input_neg_1/test_neg_1.shuffle.fa'
            neg_random_f = f'{project_root}/train/results/train_test_dataset/input_neg_random/test_neg_random.shuffle.fa'
        CONSTANT_SIZE = 1000
        #################################
        ## Processing 'Positive_MANE' samples
        #################################
        pp_MANE_idx = 0
        with open(pos_MANE_f, "r") as f:
            print("\nProcessing ", pos_MANE_f)
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
                    print("\tpp_MANE_idx: ", pp_MANE_idx)
                if pp_MANE_idx >= CONSTANT_SIZE:
                    break
        print("\tpp_MANE_idx: ", pp_MANE_idx)
        #################################
        ## Processing 'Positive_ALTS' samples
        #################################
        pp_alts_idx = 0
        with open(pos_ALTS_f, "r") as f:
            print("\nProcessing ", pos_ALTS_f)
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
                    print("\tpp_alts_idx: ", pp_alts_idx)
                if pp_alts_idx >= CONSTANT_SIZE:
                    break
        print("\tpp_alts_idx: ", pp_alts_idx)
        #################################
        ## Processing 'NEGATIVE_1' samples
        #################################
        n1idx = 0
        with open(neg_1_f, "r") as f:
            print("\nProcessing ", neg_1_f)
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
                    print("\tn1idx: ", n1idx)
                if n1idx >= CONSTANT_SIZE:
                    break
        print("\tn1idx: ", n1idx)
        #################################
        ## Processing 'NEGATIVE_Random' samples
        #################################
        nridx = 0
        with open(neg_random_f, "r") as f:
            print("\nProcessing ", neg_random_f)
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
                    print("\tnridx: ", nridx)
                if nridx >= CONSTANT_SIZE:
                    break
        print("\tnridx: ", nridx)
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
    trainset_origin = myDatasetTrain("train", int(SEQ_LEN))
    trainset, valset = torch.utils.data.random_split(trainset_origin, [0.9, 0.1])
    # trainset = myDatasetTrain("train", int(SEQ_LEN))
    testset = myDatasetTrain("test", int(SEQ_LEN))
    train_loader = DataLoader(
        trainset,
        batch_size=batch_size,
        shuffle=True,
        drop_last=False,
        pin_memory=True,
    )
    val_loader = DataLoader(
        valset,
        batch_size=batch_size,
        shuffle=True,
        drop_last=False,
        pin_memory=True,
    )
    test_loader = DataLoader(
        testset,
        batch_size = batch_size,
        drop_last=False,
        pin_memory=True,
    )
    #######################################
    # predicting splice / non-splice
    #######################################
    return train_loader, val_loader, test_loader
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


def get_eval_dataloader(batch_size, TARGET, n_workers, shuffle, eval_select, test_f):
    #######################################
    # predicting splice / non-splice
    #######################################
    testset = myDatasetEval("eval", int(SEQ_LEN), shuffle, eval_select, test_f)
    test_loader = DataLoader(
        testset,
        batch_size = batch_size,
        shuffle = shuffle,
        drop_last = False,
        pin_memory = True,
    )
    dataLoader = os.path.join(os.path.dirname(test_f), "splam_dataloader.pt")
    print(f'[INFO] Loading dataset (shuffle: {test_f}')
    torch.save(test_loader, f'{dataLoader}')
    return test_loader
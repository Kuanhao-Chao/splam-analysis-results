import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import torch.nn as nn
# from TEST_dataset import *
from splam_dataset_Chromsome import *
from SPLAM import *
# from splam_utils import *
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
from tqdm import tqdm
import warnings
from sklearn.metrics import precision_recall_curve, roc_curve
import pickle
import platform

warnings.filterwarnings("ignore")

JUNC_THRESHOLD = 0.5


def parse_junction(name):
    print("name: ", name)
    res = name.split(":")
    strand = name[-2]
    chr_name = res[0]
    if strand == "+":
        start = int(res[1].split("-")[0])+200
        end = int(res[2].split("-")[1].split('(')[0])-200
    elif strand == "-":
        start = int(res[2].split("-")[0])+200
        end = int(res[1].split("-")[1].split('(')[0])-200
    print("chr_name: ", chr_name, start, end, strand)
    return (chr_name, start, end, strand)


def main(MODEL_VERSION):
    # Global variable definition
    EPOCH_NUM = 20
    BATCH_SIZE = 100
    N_WORKERS = 1
    # Selecting device
    device_str = None
    if torch.cuda.is_available():
        device_str = "cuda"
    else:
        if platform.system() == "Darwin":
            device_str = "mps"
        else:
            device_str = "cpu"
    device = torch.device(device_str)
    print(f"\033[1m[Info]: Use {device} now!\033[0m")

    project_root = '/ccb/cybertron/khchao/splam-analysis-results/'
    MODEL_OUTPUT_BASE = f'{project_root}/results/albation_study/MODEL_TEST/{MODEL_VERSION}/'
    os.makedirs(MODEL_OUTPUT_BASE, exist_ok=True)
    
    criterion = nn.BCELoss()
    BATCH_SIZE = 100
    junc_counter = 0
    target = "test_juncs"
    os.makedirs(MODEL_OUTPUT_BASE+target, exist_ok=True)
    d_score_tsv_f = MODEL_OUTPUT_BASE+target+"/splam_all_seq.score.d."+target+".tsv"
    a_score_tsv_f = MODEL_OUTPUT_BASE+target+"/splam_all_seq.score.a."+target+".tsv"
    d_score_fw = open(d_score_tsv_f, "a")
    a_score_fw = open(a_score_tsv_f, "a")
    train_loader, val_loader, test_loader = get_train_dataloader(BATCH_SIZE, MODEL_VERSION, N_WORKERS)
    print(f"[Info]: Finish loading data!", flush = True)
    print("valid_iterator: ", len(test_loader))
    LOG_OUTPUT_TEST_BASE = MODEL_OUTPUT_BASE + "/" + target + "/LOG/"
    os.makedirs(LOG_OUTPUT_TEST_BASE, exist_ok=True)
    ############################
    # Log for testing
    ############################
    OUT_SCORE = LOG_OUTPUT_TEST_BASE + "_junction_score.bed"
    test_log_loss = LOG_OUTPUT_TEST_BASE + "test_loss.txt"
    test_log_A_acc = LOG_OUTPUT_TEST_BASE + "test_A_accuracy.txt"
    test_log_A_auprc = LOG_OUTPUT_TEST_BASE + "test_A_auprc.txt"
    test_log_A_threshold_precision = LOG_OUTPUT_TEST_BASE + "test_A_threshold_precision.txt"
    test_log_A_threshold_recall = LOG_OUTPUT_TEST_BASE + "test_A_threshold_recall.txt"
    test_log_D_acc = LOG_OUTPUT_TEST_BASE + "test_D_accuracy.txt"
    test_log_D_auprc = LOG_OUTPUT_TEST_BASE + "test_D_auprc.txt"
    test_log_D_threshold_precision = LOG_OUTPUT_TEST_BASE + "test_D_threshold_precision.txt"
    test_log_D_threshold_recall = LOG_OUTPUT_TEST_BASE + "test_D_threshold_recall.txt"
    test_log_J_threshold_precision = LOG_OUTPUT_TEST_BASE + "test_J_threshold_precision.txt"
    test_log_J_threshold_recall = LOG_OUTPUT_TEST_BASE + "test_J_threshold_recall.txt"
    fw_test_log_loss = open(test_log_loss, 'w')
    fw_test_log_A_acc = open(test_log_A_acc, 'w')
    fw_test_log_A_auprc = open(test_log_A_auprc, 'w')
    fw_test_log_A_threshold_precision = open(test_log_A_threshold_precision, 'w')
    fw_test_log_A_threshold_recall = open(test_log_A_threshold_recall, 'w')

    fw_test_log_D_acc = open(test_log_D_acc, 'w')
    fw_test_log_D_auprc = open(test_log_D_auprc, 'w')
    fw_test_log_D_threshold_precision = open(test_log_D_threshold_precision, 'w')
    fw_test_log_D_threshold_recall = open(test_log_D_threshold_recall, 'w')

    fw_test_log_J_threshold_precision = open(test_log_J_threshold_precision, 'w')
    fw_test_log_J_threshold_recall = open(test_log_J_threshold_recall, 'w')

    for model_idx in range(0, 15):
        MODEL = f'/ccb/cybertron/khchao/splam-analysis-results/results/albation_study/MODEL/subset_10000/{MODEL_VERSION}/splam_{model_idx}.pt'
        model = torch.load(MODEL)
        model = model.to(device)
        print("########################################")
        print(" Model: ", model)
        print("########################################")
        print(f"[Info]: Finish loading model!",flush = True)
        
        epoch_loss = 0
        epoch_acc = 0
        epoch_donor_acc = 0
        epoch_acceptor_acc = 0
        print("**********************")
        print("** Testing Dataset **")
        print("**********************")
        pbar = tqdm(total=len(test_loader), ncols=0, desc="Test", unit=" step")
        A_G_TP = 1e-6
        A_G_FN = 1e-6
        A_G_FP = 1e-6
        A_G_TN = 1e-6
        D_G_TP = 1e-6
        D_G_FN = 1e-6
        D_G_FP = 1e-6
        D_G_TN = 1e-6
        J_G_TP = 1e-6
        J_G_FN = 1e-6
        J_G_FP = 1e-6
        J_G_TN = 1e-6
        #######################################
        # Important => setting model into evaluation mode
        #######################################
        model.eval()
        for batch_idx, data in enumerate(test_loader):
            DNAs, labels, chrs = data
            DNAs = DNAs.to(torch.float32).to(device)
            labels = labels.to(torch.float32).to(device)
            DNAs = torch.permute(DNAs, (0, 2, 1))
            labels = torch.permute(labels, (0, 2, 1))
            loss, yp = model_fn(DNAs, labels, model, criterion)

            #######################################
            # predicting all bp.
            #######################################
            is_expr = (labels.sum(axis=(1,2)) >= 1)
            Acceptor_YL = labels[is_expr, 1, :].flatten().to('cpu').detach().numpy()
            Acceptor_YP = yp[is_expr, 1, :].flatten().to('cpu').detach().numpy()
            Donor_YL = labels[is_expr, 2, :].flatten().to('cpu').detach().numpy()
            Donor_YP = yp[is_expr, 2, :].flatten().to('cpu').detach().numpy()

            A_YL = labels[is_expr, 1, :].to('cpu').detach().numpy()
            A_YP = yp[is_expr, 1, :].to('cpu').detach().numpy()
            D_YL = labels[is_expr, 2, :].to('cpu').detach().numpy()
            D_YP = yp[is_expr, 2, :].to('cpu').detach().numpy()
            np.savetxt(d_score_fw, D_YP, delimiter=" ")
            np.savetxt(a_score_fw, A_YP, delimiter=" ")

            donor_labels, donor_scores, acceptor_labels, acceptor_scores = get_donor_acceptor_scores(D_YL, A_YL, D_YP, A_YP)

            # Junction statistics
            J_G_TP, J_G_FN, J_G_FP, J_G_TN, J_TP, J_FN, J_FP, J_TN = print_junc_statistics(D_YL, A_YL, D_YP, A_YP, JUNC_THRESHOLD, J_G_TP, J_G_FN, J_G_FP, J_G_TN)
            # Top-k statistics
            A_accuracy, A_auprc = print_top_1_statistics(Acceptor_YL, Acceptor_YP)
            D_accuracy, D_auprc = print_top_1_statistics(Donor_YL, Donor_YP)
            # Donor and Acceptor statistics
            A_G_TP, A_G_FN, A_G_FP, A_G_TN, A_TP, A_FN, A_FP, A_TN = print_threshold_statistics(Acceptor_YL, Acceptor_YP, JUNC_THRESHOLD, A_G_TP, A_G_FN, A_G_FP, A_G_TN)
            D_G_TP, D_G_FN, D_G_FP, D_G_TN, D_TP, D_FN, D_FP, D_TN = print_threshold_statistics(Donor_YL, Donor_YP, JUNC_THRESHOLD, D_G_TP, D_G_FN, D_G_FP, D_G_TN)
            batch_loss = loss.item()
            epoch_loss += loss.item()
            epoch_donor_acc += D_accuracy
            epoch_acceptor_acc += A_accuracy

            pbar.update(1)
            pbar.set_postfix(
                batch_id=batch_idx,
                idx_test=len(test_loader)*BATCH_SIZE,
                loss=f"{batch_loss:.6f}",
                A_accuracy=f"{A_accuracy:.6f}",
                D_accuracy=f"{D_accuracy:.6f}",
                A_auprc = f"{A_auprc:.6f}",
                D_auprc = f"{D_auprc:.6f}",
                A_Precision=f"{A_TP/(A_TP+A_FP+1e-6):.6f}",
                A_Recall=f"{A_TP/(A_TP+A_FN+1e-6):.6f}",
                D_Precision=f"{D_TP/(D_TP+D_FP+1e-6):.6f}",
                D_Recall=f"{D_TP/(D_TP+D_FN+1e-6):.6f}",
                J_Precision=f"{J_TP/(J_TP+J_FP+1e-6):.6f}",
                J_Recall=f"{J_TP/(J_TP+J_FN+1e-6):.6f}"
            )
            fw_test_log_loss.write(str(batch_loss)+ "\n")
            fw_test_log_A_acc.write(str(A_accuracy)+ "\n")
            fw_test_log_A_auprc.write(str(A_auprc)+ "\n")
            fw_test_log_A_threshold_precision.write(f"{A_TP/(A_TP+A_FP+1e-6):.6f}\n")
            fw_test_log_A_threshold_recall.write(f"{A_TP/(A_TP+A_FN+1e-6):.6f}\n")
            fw_test_log_D_acc.write(str(D_accuracy)+ "\n")
            fw_test_log_D_auprc.write(str(D_auprc)+ "\n")
            fw_test_log_D_threshold_precision.write(f"{D_TP/(D_TP+D_FP+1e-6):.6f}\n")
            fw_test_log_D_threshold_recall.write(f"{D_TP/(D_TP+D_FN+1e-6):.6f}\n")
            fw_test_log_J_threshold_precision.write(f"{J_TP/(J_TP+J_FP+1e-6):.6f}\n")
            fw_test_log_J_threshold_recall.write(f"{J_TP/(J_TP+J_FN+1e-6):.6f}\n")
        pbar.close()
        print(f'Epoch {batch_idx+0:03}: | Loss: {epoch_loss/len(test_loader):.5f} | Donor top-k Acc: {epoch_donor_acc/len(test_loader):.3f} | Acceptor top-k Acc: {epoch_acceptor_acc/len(test_loader):.3f}')

        print(f'Junction Precision: {J_G_TP/(J_G_TP+J_G_FP):.5f} | Junction Recall: {J_G_TP/(J_G_TP+J_G_FN):.5f} | TP: {J_G_TP} | FN: {J_G_FN} | FP: {J_G_FP} | TN: {J_G_TN}')
        print(f'Donor Precision   : {D_G_TP/(D_G_TP+D_G_FP):.5f} | Donor Recall   : {D_G_TP/(D_G_TP+D_G_FN):.5f} | TP: {D_G_TP} | FN: {D_G_FN} | FP: {D_G_FP} | TN: {D_G_TN}')
        print(f'Acceptor Precision: {A_G_TP/(A_G_TP+A_G_FP):.5f} | Acceptor Recall: {A_G_TP/(A_G_TP+A_G_FN):.5f} | TP: {A_G_TP} | FN: {A_G_FN} | FP: {A_G_FP} | TN: {A_G_TN}')
        print("\n\n")

    d_score_fw.close()
    a_score_fw.close()
    fw_test_log_loss.close()
    fw_test_log_A_acc.close()
    fw_test_log_A_auprc.close()
    fw_test_log_A_threshold_precision.close()
    fw_test_log_A_threshold_recall.close()
    fw_test_log_D_acc.close()
    fw_test_log_D_auprc.close()
    fw_test_log_D_threshold_precision.close()
    fw_test_log_D_threshold_recall.close()
    fw_test_log_J_threshold_precision.close()
    fw_test_log_J_threshold_recall.close()


def plot_pr_curve(true_y, y_prob, label):
    """
    plots the roc curve based of the probabilities
    """
    fpr, tpr, thresholds = precision_recall_curve(true_y, y_prob)
    plt.plot(fpr, tpr, label=label)
    plt.legend()
    plt.xlabel('Recall')
    plt.ylabel('Precision')


def plot_roc_curve(true_y, y_prob, label):
    """
    plots the roc curve based of the probabilities
    """
    fpr, tpr, thresholds = roc_curve(true_y, y_prob)
    plt.plot(fpr, tpr, label=label)
    plt.legend()
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')


if __name__ == "__main__":

    # rsg = 5
    # for rsb in range(1, 4):
    #     MODEL_VERSION = f'rsg_{rsg}__rsb_{rsb}/'
    #     main(MODEL_VERSION)

    rsb = 4
    for rsg in range(5, 6):
        MODEL_VERSION = f'rsg_{rsg}__rsb_{rsb}/'
        main(MODEL_VERSION)

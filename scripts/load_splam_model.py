import torch
import torch.nn as nn
# from SPLAM import *
import numpy as np
import warnings
import platform

warnings.filterwarnings("ignore")

def main():
    #############################
    # Global variable definition
    #############################
    EPOCH_NUM = 20
    BATCH_SIZE = 100
    N_WORKERS = 1

    #############################
    # Selecting device
    #############################
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
    MODEL = "../model/splam.pt"
    print(">> Using model: ", MODEL)
    model = torch.load(MODEL)

    #############################
    # Model Initialization
    #############################
    print(f"[Info]: Finish loading model!",flush = True)
    print("model: ", model)

if __name__ == "__main__":
    main()

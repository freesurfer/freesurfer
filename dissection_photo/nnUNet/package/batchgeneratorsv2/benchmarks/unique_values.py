import torch
import numpy as np
from time import time
import pandas as pd

def unique_torch(tensor):
    return torch.unique(tensor)

def unique_npy(tensor):
    return np.unique(tensor.numpy())

def unique_pandas(tensor):
    np.sort(pd.unique(tensor.numpy().ravel()))

def unique_bincount(tensor):
    return torch.where(torch.bincount(tensor.ravel()) > 0)[0]


if __name__ == '__main__':
    torch.set_num_threads(1)
    shape = (64, 64, 64)
    labels = 200

    times = []
    for _ in range(10):
        seg = torch.round(torch.rand(shape) * 20, decimals=0).to(torch.uint8)
        st = time()
        unique = unique_torch(seg)
        times.append(time() - st)
    print('unique_torch', np.median(times))

    times = []
    for _ in range(10):
        seg = torch.round(torch.rand(shape) * 20, decimals=0).to(torch.uint8)
        st = time()
        unique = unique_npy(seg)
        times.append(time() - st)
    print('unique_npy', np.median(times))

    times = []
    for _ in range(10):
        seg = torch.round(torch.rand(shape) * 20, decimals=0).to(torch.uint8)
        st = time()
        unique = unique_pandas(seg)
        times.append(time() - st)
    print('unique_pandas', np.median(times))

    times = []
    for _ in range(10):
        seg = torch.round(torch.rand(shape) * 20, decimals=0).to(torch.uint8)
        st = time()
        unique = unique_bincount(seg)
        times.append(time() - st)
    print('unique_bincount', np.median(times))


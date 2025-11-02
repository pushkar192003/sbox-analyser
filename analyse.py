
"""
Author: pushkar raj
THeoratical analysis tool for 8*8 

"""
from collections import Counter
import numpy as np
import csv
import os
import math


SBOX_FLAT = None


CSV_FILENAME = "SM4_SBOX.csv"

import csv

def load_sm4_sbox(path):
    sbox = []
    with open(path, newline='') as f:
        reader = csv.reader(f)
        rows = list(reader)
  
        for row in rows[1:]:  
           
            values = row[1:]
            sbox.extend([int(v, 16) for v in values])
    return sbox

SBOX_FLAT = load_sm4_sbox("SM4_SBOX.csv")

def is_bijective(sbox):
    return len(sbox) == 256 and len(set(sbox)) == 256

def fixed_points(sbox):
    return [i for i,v in enumerate(sbox) if i == v]

def compute_inverse_sbox(sbox):
    inv = [None]*256
    for i,v in enumerate(sbox):
        if inv[v] is not None:
     
            raise ValueError(f"sbox not bijective: value {v:02X} seen at least twice (first at {inv[v]:02X}, now at {i:02X})")
        inv[v] = i
    if any(x is None for x in inv):
        raise ValueError("sbox not bijective: some outputs missing")
    return inv

def walsh_transform_boolean(f):

    W = np.zeros(256, dtype=int)
 
    xs = np.arange(256, dtype=np.uint16)
    for a in range(256):
        ax = (a & xs)
        parity = np.unpackbits(ax.astype(np.uint8).view(np.uint8)).reshape(-1,2)  
       
        tot = 0
        for x in range(256):
            bit = (bin(a & x).count("1") & 1) ^ f[x]
            tot += -1 if bit else 1
        W[a] = tot
    return W

def walsh_transform_boolean_fast(f):
  
    xs = np.arange(256, dtype=np.uint16)
    W = np.zeros(256, dtype=int)
    f_arr = np.array(f, dtype=np.uint8)
    for a in range(256):

        ax = a & xs
        parity = np.array([bin(val).count("1") & 1 for val in ax], dtype=np.uint8)
        expr = (f_arr ^ parity)
    
        tot = expr
        tot = np.where(expr == 0, 1, -1).sum()
        W[a] = int(tot)
    return W

def nonlinearity_per_bit(sbox):

    nl = []
    for bit in range(8):
        f = [ (v >> bit) & 1 for v in sbox ]
        W = walsh_transform_boolean_fast(f)
        maxW = np.max(np.abs(W))
        nl_bit = 128 - (maxW // 2)
        nl.append(int(nl_bit))
    return nl

def build_ddt(sbox):
  
    ddt = [[0]*256 for _ in range(256)]
    for dx in range(256):
        for x in range(256):
            y = sbox[x] ^ sbox[x ^ dx]
            ddt[dx][y] += 1
    return ddt

def differential_uniformity(ddt):
  
    maxval = 0
    for dx in range(1,256):
        for y in range(256):
            if ddt[dx][y] > maxval:
                maxval = ddt[dx][y]
    return maxval

def avalanche_measure(sbox):

    per_input = []
    for i in range(8):
        total_changes = 0
        for x in range(256):
            y1 = sbox[x]
            y2 = sbox[x ^ (1 << i)]
            total_changes += bin(y1 ^ y2).count("1")
        avg = total_changes / 256.0
        per_input.append(avg)
    overall_avg = sum(per_input) / 8.0
    return per_input, overall_avg


def run_all_checks(sbox_flat):
    if len(sbox_flat) != 256:
        raise ValueError("S-box must have 256 entries")

    print("Running S-box checks...\n")
    
    bij = is_bijective(sbox_flat)
    print("Bijective:", bij)
    if not bij:
        counts = Counter(sbox_flat)
        dup = [f"{v:02X}" for v,c in counts.items() if c>1]
        missing = [f"{i:02X}" for i in range(256) if i not in counts]
        print("  Duplicates (values present >1):", dup)
        print("  Missing outputs:", missing)
    fp = fixed_points(sbox_flat)
    print("Fixed points (indices where S[x]==x):", [f"{i:02X}" for i in fp] if fp else "None")


    try:
        inv = compute_inverse_sbox(sbox_flat)
        inv_ok = all(inv[sbox_flat[x]] == x for x in range(256))
        print("Inverse computed. Decryptibility (Inv[S[x]]==x):", "Yes" if inv_ok else "No")
    except Exception as e:
        print("Inverse computation failed:", e)
        inv = None

    print("\nComputing nonlinearity per output bit (this may take a few seconds)...")
    nl = nonlinearity_per_bit(sbox_flat)
    print("Nonlinearity per output bit (bit 0..7):", nl)
    print("Average nonlinearity:", sum(nl)/8.0)

    print("\nBuilding DDT (Difference Distribution Table) ...")
    ddt = build_ddt(sbox_flat)
    du = differential_uniformity(ddt)
    print("Differential uniformity (max entry for dx != 0):", du)

    per_input, overall = avalanche_measure(sbox_flat)
    print("\nAvalanche effect (average number of output bits changed when flipping each input bit):")
    for i,avg in enumerate(per_input):
        print(f"  Input bit {i}: avg output bits changed = {avg:.4f}")
    print(f"Overall average (over all input bits): {overall:.4f}")

    return {
        "bijective": bij,
        "fixed_points": fp,
        "inverse": inv,
        "nonlinearity_per_bit": nl,
        "du": du,
        "ddt": ddt,
        "avalanche_per_input": per_input,
        "avalanche_overall": overall
    }

if __name__ == "__main__":
   
    results = run_all_checks(SBOX_FLAT)

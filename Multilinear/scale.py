import pandas as pd
from enum import Enum
import os
import csv
import shutil
import numpy as np
import math 
import matplotlib.pyplot as plt
import re



df = pd.read_csv('instancedata.csv', sep=';')
opt_dict = {}

for row_ind in range(len(df)):
    row = df.iloc[row_ind]
    if row["name"].find("autocorr") != -1:
        entry = {}
        entry["instance"] = row["name"]
        entry["opt"] = float(row["primalbound"])
        opt_dict[entry["instance"]] = entry

    

bench_path = os.getcwd() + "/benchmark"
scale_file  = open(os.getcwd() + "/scales.txt", "w")
scaledbench_path = os.getcwd() + "/scaledbenchmark"
benchs = os.listdir(bench_path)

def parse(lines):
    newline = ""
    line = " ".join(line.strip() for line in lines) 
    head_index = line.index("objvar")
    head_index += 6 
    head = line[0: head_index]
    tail_index = line.index(" <=")
    tail = line[tail_index:]
    body = line[head_index: tail_index]
    body_es = re.split('(\-|\+)', body)
    if body_es[0] == "":
        body_es=body_es[1:]

    vars = set()
    uni_terms = []
    bi_terms = []
    tri_terms = []
    quad_terms = []

    sum_coef = 0
    uni_sum = 0
    bi_sum = 0
    tri_sum = 0
    quad_sum = 0
    for i, e in enumerate(body_es):
        if e == "+" or e == "-":
            term = body_es[i+1]
        else:
            continue
        spterm = term.split(" ")
        spterm = [s for s in spterm if s != ""]
        coef = int(float(spterm[0]))
        if len(spterm) - 1  == 1:
            var1 = spterm[1]
            vars.add(var1)
            uni_terms.append((e, coef, var1))
            uni_sum += 1
        elif len(spterm) - 1  == 2:
            var1 = spterm[1]
            var2 = spterm[2]
            vars.add(var1)
            vars.add(var2)
            bi_terms.append((e, coef, var1, var2))
            bi_sum += 1
        elif len(spterm) - 1  == 3:
            var1 = spterm[1]
            var2 = spterm[2]
            var3 = spterm[3]
            vars.add(var1)
            vars.add(var2)
            vars.add(var3)
            tri_terms.append((e, coef, var1, var2, var3))
            tri_sum += 1
        elif len(spterm) - 1  == 4:
            var1 = spterm[1]
            var2 = spterm[2]
            var3 = spterm[3]
            var4 = spterm[4]
            vars.add(var1)
            vars.add(var2)
            vars.add(var3)
            vars.add(var4)
            quad_terms.append((e, coef, var1, var2, var3, var4))
            quad_sum += 1
    newbody = head
    assert(quad_sum != 0)
    if quad_sum != 0:
        sum_coef = quad_sum
    sum_coef_ = sum_coef 
    sum_coef = float(sum_coef)
    for term in uni_terms:
        newbody += " " + term[0] + " " + str(round(term[1]/ sum_coef, 12)) + " " + term[2]
    for term in bi_terms:
        newbody += " " + term[0] + " " + str(round(term[1]/ sum_coef, 12)) + " " + term[2] + " " + term[3]
    for term in tri_terms:
        newbody += " " + term[0] + " " + str(round(term[1]/ sum_coef, 12)) + " " + term[2] + " " + term[3] + " " + term[4]
    for term in quad_terms:
        newbody += " " + term[0] + " " + str(round(term[1]/ sum_coef, 12)) + " " + term[2] + " " + term[3] + " " + term[4] + " " + term[5]
    newbody += tail
    return newbody, sum_coef_

lst = []
for bench in benchs:
    file_p = bench_path + "/" + bench
    newfile_p = scaledbench_path + "/" + "scaled" + bench
    data = open(file_p ,"r")
    newdata = open(newfile_p, "w")
    lines = data.readlines()
    start_i = -1
    end_i = -1
    all_i = len(lines)
    for i, line in enumerate(lines):
        #print(i, all_i, line)
        if "Subject To" in line:
            start_i = i+1
        if start_i != -1 and end_i == -1 and "Bounds" in line:
            end_i = i - 1
    print(start_i, end_i)
    newbody, scale = parse(lines[start_i: end_i])
    maxlen = 110
    lst_bd = []
    len_bd = len(newbody)
    pt = 0
    while pt < len_bd:
        j = 0
        while (pt + maxlen + j) < len_bd - 1 and newbody[pt+maxlen + j] != " ":
            j += 1
        next_pt = min(pt + maxlen + j, len_bd)
        lst_bd.append(newbody[pt:next_pt])
        lst_bd.append("\n")
        pt = next_pt
    lst_bd.append("\n")
    newdata.writelines(lines[0:start_i])
    newdata.writelines(lst_bd)
    newdata.writelines(lines[end_i:])
    data.close()
    scale_ = bench[0:-4] + " " + str(scale) + "\n"
    lst.append(scale_)

scale_file.writelines(lst)
scale_file.close()
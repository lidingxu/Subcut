import pandas as pd
from enum import Enum
import os
import csv
import shutil
import numpy as np
import math 
import matplotlib.pyplot as plt



def extract_scip(entry_, file_path):
    ls = open(file_path).readlines()
    #print(ls)
    #print(file)
    stat_keys = ["Total Time",  "nodes", "intersub", "Gap", "interlattice", "solving"]
    stat_dict = {}
    entry = entry_
    for l in ls:
        for stat_key in stat_keys:
            if  l.split(":")[0].strip() == stat_key:
                stat_dict[stat_key] = l
    #print(stat_dict["Dual Bound"].split()[2])
    entry["cut_time"] = float(stat_dict["intersub"].split()[2])  + float(stat_dict["interlattice"].split()[2]) 
    entry["nodes"] = int(stat_dict["nodes"].split()[2])
    entry["Gap"] = 1 - float( stat_dict["Gap"].split()[2]) / 100
    entry["napplied"] =  int(stat_dict["intersub"].split()[11])  + int(stat_dict["interlattice"].split()[11]) 
    entry["solving"] = float(stat_dict["solving"].split()[2])
    #entry["rel_gap"] = float(if stat_dict["Gap"].split()[2] == "infinite") #float(stat_dict["Gap"].split()[2])
    return entry


def parse_name(name):
    s = ""
    for c in name:
        if c == "_":
            s += "\_"
        else:
            s +=c
    return s





opt_file = open("opt.txt", "r")
ls = opt_file.readlines()

opt_dict = {}

pclasses = ["g05", "pw"]

for line in ls:
    line = line.split(" ")
    #print(line)
    entry = {}
    entry["instance"] = line[0]
    entry["opt"] = float(line[2])
    line = line[0].split("_")
    if line[0] == "g05":
        entry["pclass"] = "g05"
        line = line[1].split(".")
        entry["subclass"] = line[0]
    else:
        entry["pclass"] = "pw"
        line = line[0].split(".")
        entry["subclass"] = line[0][2:4]
    #entry = extract_scip(entry, continuous_log_path + "/"+log)
    opt_dict[entry["instance"]] = entry



log_path = os.getcwd() + "/bblogs_60"
logs = os.listdir(log_path)

entries = []


instances = set()


pclassmap = {"g05": ["g05_60", "g05_80", "g05_100", "pw09_100"], "pw": ["pw09_100"]}

def getclass(instance_main):
    if instance_main in pclassmap["g05"]:
        return "g05"
    elif instance_main in pclassmap["pw"]:
        return "pw"

for log in logs:
    log_ = log
    log = log[0:-4]
    setting = log.split("_")[2]
    activate = setting[-1] 
    setting = setting[0:-1] 
    log = log.split("_")[0] +"_" + log.split("_")[1]
    instance = log.split(".")[0] +"." + log.split(".")[1]
    entry={}
    instances.add(instance)
    entry["instance"] = instance
    entry["activate"] = activate
    entry["setting"] = setting
    entry["pclass"] = getclass(instance.split(".")[0])
    entry = extract_scip(entry, log_path + "/"+log_)
    entries.append(entry)


#print(entries)
    

def Stat(aname, sname, pname):
    return {"activate": aname, "setting": sname, "pclass": pname,  "solved": 0, "Gap": 0, "nodes": 0.0, "total": 0, "cut_time": 0.0,  "napplied": 0, "relative": 1., "relative_lst": [], "nodes_lst": [], "cut_time_lst": [], "napplied_lst": [], "Gap_lst": []} 

display_keys = ["Gap", "relative", "nodes", "cut_time",  "napplied", "total"]



activates = ["a"]
settings = ["default", "icuts", "icutsl"]

classstats = {}

defaultdual = {}

for activate in activates:
    for pclass in pclasses:
        for setting in settings:
            classstats[(activate, setting, pclass)] = []
            for entry in entries:
                if entry["activate"] == activate and entry["pclass"] == pclass and entry["setting"] == setting:
                    #print("1")
                    instance = entry["instance"]
                    if setting == "default":
                        defaultdual[(activate, instance)] = entry["Gap"]
                    entry["relative"] = entry["Gap"] /  defaultdual[(activate, instance)]
                    classstats[(activate, setting, pclass)].append(entry)


def add(stat, entry):
    stat["total"] += 1
    stat["cut_time_lst"].append(entry["cut_time"])
    stat["napplied_lst"].append(entry["napplied"])
    stat["Gap_lst"].append(entry["Gap"]) 
    stat["nodes_lst"].append(entry["nodes"])
    stat["relative_lst"].append(entry["relative"])


def SGM(lst, total, bias):
    return np.exp(np.sum([np.log(ele + bias) for ele in lst ]) / total) - bias


def avgStat(stat):
    stat["Gap"] =  SGM(stat["Gap_lst"], stat["total"], 1)
    stat["nodes"] = SGM(stat["nodes_lst"], stat["total"], 1)
    stat["cut_time"] = SGM(stat["cut_time_lst"] , stat["total"], 1)
    stat["napplied"] = SGM(stat["napplied_lst"] , stat["total"], 1)
    stat["relative"] = SGM(stat["relative_lst"] , stat["total"], 1)



def printStat(setting, pclass, stat):
    #print(setting, pclass, stat)
    s = [ str(round(stat[display_key], 3) if display_key is "Gap" or display_key is "relative" else round(stat[display_key], 2)) + " & " for display_key in display_keys]
    print(setting, pclass, "".join(s))



activates = "a"
stats = {}
allstat = {}
for activate in activates:
    for pclass in pclasses:
        for setting in settings:
            allstat[(activate, setting, pclass)] =  Stat(activate, setting, pclass)
            for entry in classstats[(activate, setting, pclass)]:
                add(allstat[(activate, setting, pclass)], entry)
            avgStat(allstat[(activate, setting, pclass)])
            printStat(setting, pclass, allstat[(activate, setting, pclass)])



pairs = (("default", "icuts"), ("default", "icutsl"), ("icutsl","icuts"))

names = {"default": "Default", "icutsl": "Split cut", "icuts": "Submodular cut"}
data= {}
dataname = {}
for i in range(0,2):
    for j in range(0,3):
        data[(pclasses[i], pairs[j][0])] = []
        data[(pclasses[i], pairs[j][1])] = []
        dataname[(pclasses[i], pairs[j][1])] = []
        dataname[(pclasses[i], pairs[j][0])] = []

settings = ["default", "icuts", "icutsl"]
for instance in instances:
    for pclass in pclasses:
        for setting in settings:
            for entry in entries:
                if entry["instance"] == instance and entry["setting"] == setting and entry["pclass"] == pclass:
                    data[(pclass, setting)].append(entry["Gap"])
                    dataname[(pclass, setting)].append(entry["instance"])
                    break


fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 3))  # define the figure and subplots

for i in range(0,1):
    for j in range(0,3):
        pair = pairs[j]
        wins = [0,0]
        num = min(len( data[(pclasses[i], pair[0])]),  len(data[(pclasses[i], pair[1])]))
        maxd = 0
        for k in range(num):
            maxd = max(maxd, data[(pclasses[i], pair[0])][k], data[(pclasses[i], pair[1])][k])
            if data[(pclasses[i], pair[0])][k] > data[(pclasses[i], pair[1])][k]:
                wins[0] += 1
            else:
                wins[1] += 1
        axes[j].scatter(data[(pclasses[i], pair[1])], data[(pclasses[i], pair[0])], color = 'blue', marker = '+')
        axes[j].plot([0,maxd], [0, maxd], color = 'green')
        axes[j].set_xlabel(names[pair[1]] +" wins " + str(wins[1]))
        axes[j].set_ylabel(names[pair[0]] + " wins " +  str(wins[0]))

fig.tight_layout()
plt.savefig('scatter_qubo_bc.pdf') 
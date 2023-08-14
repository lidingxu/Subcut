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
    stat_keys = ["Total Time",  "Dual Bound", "First LP value", "intersub", "interlattice", "solving"]
    stat_dict = {}
    entry = entry_
    for l in ls:
        for stat_key in stat_keys:
            if  l.split(":")[0].strip() == stat_key:
                stat_dict[stat_key] = l
    #print(stat_dict["Dual Bound"].split()[2])
    if "Total Time" not in stat_dict  or "First LP value" not in stat_dict:
        entry["total_time"] = float("NAN")
        entry["dualbound"] = float("NAN")
        entry["FirstLP"] = float("NAN")
        entry["noLP"] = True
        entry["affected"] = False
        entry["intermis"] = 0
        return entry
    entry["total_time"] = float(stat_dict["Total Time"].split()[3])
    #print(stat_dict["Total Time"].split()[3])
    #print(stat_dict["intermis"].split())
    #entry["primal_bound"] = float(stat_dict["Primal Bound"].split()[3])
    entry["dualbound"] = float("NAN") if stat_dict["Dual Bound"].split()[3] == "-" else float(stat_dict["Dual Bound"].split()[3])
    entry["noLP"] = True if stat_dict["First LP value"].split()[4] == "-" else False
    entry["FirstLP"] =   float("NAN")  if stat_dict["First LP value"].split()[4] == "-" else float(stat_dict["First LP value"].split()[4])
    entry["affected"] =   int(stat_dict["intersub"].split()[8]) > 0 
    #print(stat_dict["intersub"].split()[8])
    entry["ncuts"] =  int(stat_dict["intersub"].split()[8])  + int(stat_dict["interlattice"].split()[8]) 
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




bench_path = os.getcwd() + "/benchmark"
benchs = os.listdir(bench_path)

opt_file = open("opt.txt", "r")
ls = opt_file.readlines()

opt_dict = {}

pclasses = ["g05", "pw"]
subpclasses = {"g05": ["60", "80", "100"], "pw": ["01", "05", "09"]}


for line in ls:
    line = line.split(" ")
    #print(line)
    entry = {}
    entry["instance"] = line[0]
    lst_data = open(bench_path + "/"+ entry["instance"] + ".cut", "r").readlines()
    sumweight = 0
    for d in lst_data[1:]:
        sumweight += float(d.split(" ")[2])
    #print(sumweight)
    entry["opt"] = float(line[2]) / sumweight
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


log_path = os.getcwd() + "/rootlogs2"
logs = os.listdir(log_path)

entries = []
instances = set()


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
    entry["pclass"] = opt_dict[instance]["pclass"]
    entry["subclass"] = opt_dict[instance]["subclass"]
    entry = extract_scip(entry, log_path + "/"+log_)
    entries.append(entry)


def Stat(aname, sname, pname, subname):
    return {"activate": aname, "setting": sname, "pclass": pname, "subclass": subname, "solved": 0, "closed": 0, "total_time": 0.0, "total": 0, "ncuts": 0,  "napplied": 0, "relative": 1., "relative_lst": [], "total_time_lst": [], "ncuts_lst": [], "napplied_lst": [], "closed_lst": []} 

display_keys = [ "closed", "relative", "total_time",  "napplied"]



activates = ["d", "a"]
settings = ["default", "icuts", "icutsl"]

subpclasses = {"g05": ["60", "80", "100"], "pw": ["01", "05", "09"]}

classstats = {}

defaultclosed = {}

for activate in activates:
    for setting in settings:
        for pclass in pclasses:
            for subpclass in subpclasses[pclass]:
                classstats[(activate, setting, pclass, subpclass)] = []
                #print(pclass,  subpclass)
                for entry in entries:
                    #print(entry)
                    #print(entry["pclass"], pclass, entry["subclass"], subpclass)
                    if entry["activate"] == activate and entry["pclass"] == pclass and entry["subclass"] == subpclass and entry["setting"] == setting:
                        #print("1")
                        instance = entry["instance"]
                        p_bd = opt_dict[instance]["opt"]
                        entry["closed"] = abs(entry["dualbound"] - entry["FirstLP"]) / abs(p_bd + entry["FirstLP"]) 
                        if setting == "default":
                            defaultclosed[(activate, instance)] = entry["closed"]
                        #print(defaultclosed, setting, pclass, subpclass)
                        entry["relative"] = ( entry["closed"] + 1e-9) / ( defaultclosed[(activate, instance)] + 1e-9)
                        #print(entry["relative"])
                        classstats[(activate, setting, pclass, subpclass)].append(entry)



def add(stat, entry):
    #print(stat, entry, entry["closed"] is float("nan"))
    #print(entry)
    #stat["solved"] += entry["issolved"] 
    stat["total"] += 1
    stat["ncuts_lst"].append(entry["ncuts"])
    stat["napplied_lst"].append(entry["napplied"])
    stat["closed_lst"].append(entry["closed"]) 
    stat["total_time_lst"].append(entry["total_time"])
    #print(entry["relative"])
    stat["relative_lst"].append(entry["relative"])


def SGM(lst, total, bias):
    return np.exp(np.sum([np.log(ele + bias) for ele in lst ]) / total) - bias


def avgStat(stat):
    stat["closed"] =  SGM(stat["closed_lst"], stat["total"], 1)
    stat["total_time"] = SGM(stat["total_time_lst"], stat["total"], 1)
    stat["ncuts"] = SGM(stat["ncuts_lst"] , stat["total"], 1)
    stat["napplied"] = SGM(stat["napplied_lst"] , stat["total"], 1)
    stat["relative"] = SGM(stat["relative_lst"] , stat["total"], 1)



def printStat(activate, setting, pclass, stat):
    #print(setting, pclass, stat)
    s = [ str(round(stat[display_key], 3) if display_key is "relative" or display_key is "closed" else round(stat[display_key], 2)) + " & " for display_key in display_keys]
    print(activate, setting, pclass, "".join(s))

stats = {}
allstat = {}
for activate in activates:
    for setting in settings:
        allstat[(activate,setting)] = Stat(activate, setting, "all", "all")
        for pclass in pclasses:
            for subpclass in subpclasses[pclass]:
                stats[(activate, setting, pclass, subpclass)] = Stat(activate, setting, pclass, subpclass)
                for entry in classstats[(activate, setting, pclass, subpclass)]:
                    add(stats[(activate, setting, pclass, subpclass)], entry)   
                    add(allstat[activate, setting], entry)             
                avgStat(stats[(activate, setting, pclass, subpclass)])
                #stats[(setting, pclass, subpclass)]["relative"] = (stats[(setting, pclass, subpclass)]["closed"] + 1e-9) /  (stats[("default", pclass, subpclass)]["closed"] + 1e-9)
                #printStat(activate, setting, pclass, stats[(activate, setting, pclass, subpclass)])    
        avgStat(allstat[(activate, setting)])
        #allstat[setting]["relative"] = (allstat[setting]["closed"] + 1e-9) /  (allstat["default"]["closed"] + 1e-9)       
        printStat(activate, setting, "all", allstat[(activate, setting)])



pairs = (("default", "icuts"), ("default", "icutsl"), ("icutsl","icuts"))

names = {"default": "Default", "icutsl": "Split cut", "icuts": "Submodular cut"}
data= {}
dataname = {}
for i in range(0,2):
    for j in range(0,3):
        data[(activates[i], pairs[j][0])] = []
        data[(activates[i], pairs[j][1])] = []
        dataname[(activates[i], pairs[j][1])] = []
        dataname[(activates[i], pairs[j][0])] = []

settings = ["default", "icuts", "icutsl"]
for instance in instances:
    for activate in activates:
        for setting in settings:
            for entry in entries:
                if entry["instance"] == instance and entry["setting"] == setting and entry["activate"] == activate:
                    data[(activate, setting)].append(entry["closed"])
                    dataname[(activate, setting)].append(entry["instance"])
                    break


fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(10, 6))  # define the figure and subplots

for i in range(0,2):
    for j in range(0,3):
        pair = pairs[j]
        wins = [0,0]
        print(len( data[(activates[i], pair[0])]),  len(data[(activates[i], pair[1])]))
        num = min(len( data[(activates[i], pair[0])]),  len(data[(activates[i], pair[1])]))
        maxd = 0
        for k in range(num):
            maxd = max(maxd, data[(activates[i], pair[0])][k], data[(activates[i], pair[1])][k])
            if data[(activates[i], pair[0])][k] > data[(activates[i], pair[1])][k]:
                print(data[(activates[i], pair[0])][k], data[(activates[i], pair[1])][k], dataname[(activates[i], pair[1])][k])
                wins[0] += 1
            else:
                wins[1] += 1
        axes[i,j].scatter(data[(activates[i], pair[1])], data[(activates[i], pair[0])], color = 'blue', marker = '+')
        axes[i,j].plot([0,maxd], [0, maxd], color = 'green')
        axes[i,j].set_xlabel(names[pair[1]] +" wins " + str(wins[1]))
        axes[i,j].set_ylabel(names[pair[0]] + " wins " +  str(wins[0]))

fig.tight_layout()
plt.savefig('scatter_qubo.pdf') 
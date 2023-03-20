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
    stat_keys = ["Total Time", "Primal Bound", "Dual Bound", "First LP value", "intersub", "interlattice", "solving"]
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
    entry["primal_bound"] = float(stat_dict["Primal Bound"].split()[3])
    entry["dualbound"] = float("NAN") if stat_dict["Dual Bound"].split()[3] == "-" else float(stat_dict["Dual Bound"].split()[3])
    entry["noLP"] = True if stat_dict["First LP value"].split()[4] == "-" else False
    entry["FirstLP"] =   float("NAN")  if stat_dict["First LP value"].split()[4] == "-" else float(stat_dict["First LP value"].split()[4])
    entry["affected"] =   int(stat_dict["intersub"].split()[8]) > 0 
    entry["ncuts"] =  int(stat_dict["intersub"].split()[8])  + int(stat_dict["interlattice"].split()[8]) 
    entry["solving"] = float(stat_dict["solving"].split()[2])
    #entry["rel_gap"] = float(if stat_dict["Gap"].split()[2] == "infinite") #float(stat_dict["Gap"].split()[2])
    return entry




def Stat(name):
    return {"setting": name, "solved": 0, "closed": 0, "total_time": 0.0, "total": 0, "ncuts": 0} 


#print(opt_dict)


log_path = os.getcwd() + "/logs"
logs = os.listdir(log_path)

entries = []
opt_dict = {}

for log in logs:
    log_ = log
    log = log[0:-4]

    instance = log.split(".")[0]
    setting = log.split("_")[-1]
    setting = setting[0:-1] if setting[0] == "d" else setting[0:6]
    activate = "standalone" if log[-1] == "s" else "embed" 
    insclass = "block" if log[0] == "b" else "normal" 


    entry={}
    entry["instance"] = instance
    entry["activate"] = activate
    entry["setting"] = setting
    entry["class"] =  insclass
    #print(instance, insclass, setting)
    entry["relative"] = 1
    #print(entry, "\n")
    entry = extract_scip(entry, log_path + "/"+log_)
    entries.append(entry)

    #print(entries)

    if not instance in opt_dict:
        entry_ = {}
        entry_["instance"] = instance
        entry_["opt"] = entry["primal_bound"]
        entry_["class"] =  insclass
        opt_dict[instance] = entry_
    elif entry["primal_bound"] < opt_dict[instance]["opt"]:
        opt_dict[instance]["opt"] = entry["primal_bound"]



def Stat(aname, sname, pname):
    return {"activate": aname,"setting": sname, "pclass": pname, "solved": 0, "closed": 0, "total_time": 0.0, "total": 0, "ncuts": 0, "relative": 1., "relative_lst": [], "total_time_lst": [], "ncuts_lst": [], "closed_lst": []} 

display_keys = ["activate","setting", "pclass", "closed", "relative", "total_time", "ncuts", "total"]


details = ""

activates = ["standalone", "embed"]
settings = ["default", "icutss", "icutsl"]
pclasses = ['block', 'normal']

classstats = {}

defaultclosed = {}




for activate in activates:
    for setting in settings:
        for pclass in pclasses:    
            classstats[(activate, setting, pclass)] = []
            #print(pclass,  subpclass)
            for entry in entries:
                if entry["activate"] == activate and entry["class"] == pclass and entry["setting"] == setting:
                    #print("1")
                    instance = entry["instance"]
                    p_bd = opt_dict[instance]["opt"]
                    entry["closed"] = abs(entry["dualbound"] - entry["FirstLP"]) / abs(p_bd - entry["FirstLP"]) 
                    #print(entry["closed"] , ": ",  entry["FirstLP"]," ", entry["dualbound"], " ",p_bd, "\n" )
                    if setting == "default":
                        defaultclosed[(activate, instance)] = entry["closed"]
                    if (activate, instance) not in defaultclosed:
                        print((activate, instance), " ", defaultclosed)
                    entry["relative"] = ( entry["closed"] + 1e-9) / ( defaultclosed[(activate, instance)] + 1e-9)
                    classstats[(activate, setting, pclass)].append(entry)
            #print(classstats[(activate, setting, pclass)], "\n")
        #print(classstats[(setting, pclass)])

def add(stat, entry):
    #print(stat, entry, entry["closed"] is float("nan"))
    #print(entry)
    #stat["solved"] += entry["issolved"] 
    stat["total"] += 1
    stat["ncuts_lst"].append(entry["ncuts"])
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
    stat["relative"] = SGM(stat["relative_lst"] , stat["total"], 1)


def printStat(activate, setting, pclass, stat):
    #print(setting, pclass, stat)
    print(activate, setting, pclass, [(display_key, stat[display_key]) for display_key in display_keys])


stats = {}
allstat = {}

for activate in activates:
    for setting in settings:
        allstat[(activate, setting)] = Stat(activate, setting, "all")
        for pclass in pclasses:
                stats[(activate, setting, pclass)] = Stat(activate, setting, pclass)
                for entry in classstats[(activate, setting, pclass)]:
                    add(stats[(activate, setting, pclass)], entry)   
                    add(allstat[(activate, setting)], entry)             
                avgStat(stats[(activate, setting, pclass)])
                #stats[(setting, pclass)]["relative"] = (stats[(setting, pclass)]["closed"] + 1e-9) /  (stats[(pclass, "default")]["closed"] + 1e-9)
                printStat(activate, setting, pclass, stats[(activate, setting, pclass)])  
        #print(allstat[(activate, setting)])
        avgStat(allstat[(activate, setting)])
        #allstat[setting]["relative"] = (allstat[setting]["closed"] + 1e-9) /  (allstat["default"]["closed"] + 1e-9)       
        printStat(activate, setting, "all", allstat[(activate, setting)])


def mystr(num_str):
    return str(round(num_str, 2))

# display
for pclass in pclasses:
    for activate in activates:
        display_str = ""
        for setting in settings:
            if setting == "default":
                display_str += mystr(stats[(activate, setting, pclass)]["closed"]) + " & " + mystr(stats[(activate, setting, pclass)]["total_time"]) + " & "
            else:
                display_str += mystr(stats[(activate, setting, pclass)]["closed"]) + " & " + mystr(stats[(activate, setting, pclass)]["relative"])  + " & " + mystr(stats[(activate, setting, pclass)]["total_time"]) + " & "  + mystr(stats[(activate, setting, pclass)]["ncuts"]) + " & "
        print(pclass, " ", activate, " ", display_str)

for activate in activates:
    display_str = ""
    for setting in settings:
        if setting == "default":
            display_str += mystr(allstat[(activate, setting)]["closed"]) + " & " + mystr(allstat[(activate, setting)]["total_time"]) + " & "
        else:
            display_str += mystr(allstat[(activate, setting)]["closed"]) + " & " + mystr(allstat[(activate, setting)]["relative"])  + " & " + mystr(allstat[(activate, setting)]["total_time"]) + " & "  + mystr(allstat[(activate, setting)]["ncuts"]) + " & "
    print(activate, " ", display_str)
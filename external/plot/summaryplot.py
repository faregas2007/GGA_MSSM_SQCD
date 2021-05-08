import matplotlib.pyplot as plt
import numpy as np
import csv


def finddatatgb(bottop, tgb, MA, who):
    global everything
    returnlists = {"tgb":[] ,"ReAvg":[], "ReErr":[], "ImAvg":[], "ImErr":[]}
    tgbcompare = [f"{x:.1f}" for x in tgb]
    for line in everything:
        if line[0] == bottop and f"{float(line[1]):.1f}" in tgbcompare and line[2] == str(MA) and line[5] == who:
            try:
                returnlists["tgb"].append(float(line[1]))
                returnlists["ReAvg"].append(float(line[6]))
                returnlists["ReErr"].append(float(line[7]))
                returnlists["ImAvg"].append(float(line[10]))
                returnlists["ImErr"].append(float(line[11]))
            except:
                print(f"Error for {line[1]}")
    return returnlists

def finddataMA(bottop, tgb, MA, who):
    global everything
    returnlists = {"MA":[] ,"ReAvg":[], "ReErr":[], "ImAvg":[], "ImErr":[]}
    for line in everything:
        if line[0] == bottop and f"{float(line[1]):.1f}" == f"{tgb:.1f}" and int(line[2]) in MA and line[5] == who:
            try:
                returnlists["MA"].append(float(line[2]))
                returnlists["ReAvg"].append(float(line[6]))
                returnlists["ReErr"].append(float(line[7]))
                returnlists["ImAvg"].append(float(line[10]))
                returnlists["ImErr"].append(float(line[11]))
            except:
                print(f"Error for {line[2]}")
    return returnlists


fig, ax = plt.subplots(nrows = 2, ncols = 3, figsize = (18, 12),)

everything = []

try:
    for line in csv.reader(open("KIT/summaryfile.csv")):
        everything.append(line)
except:
    print("File KIT/summaryfile.csv not found")
try:
    for line in csv.reader(open("LUKAS/summaryfile.csv")):
        everything.append(line)
except:
    print("File LUKAS/summaryfile.csv not found")

try:
    for line in csv.reader(open("Michael/summaryfile.csv")):
        everything.append(line)
except:
    print("File Michael/summaryfile.csv not found")


# tgb = 5 run
sq =  ["stop", "sbot"]
people = ["MICHAEL","LUKAS", "KIT"]
colors = {"MICHAEL": 'lime', "LUKAS": 'red', "KIT": 'blue'}

for bottop in range(2):
    ax[bottop, 0].set_title(sq[bottop] + r", $\tan \beta = 5$")
    for who in people:
        data = finddataMA(sq[bottop], 5, range(1000, 2801, 100), who)
        # print(data)
        ax[ bottop,0].errorbar(data["MA"], data["ReAvg"],yerr = data["ReErr"], label = f"{who}, Real", ls = '--', color = colors[who])
        ax[ bottop,0].errorbar(data["MA"], data["ImAvg"],yerr = data["ImErr"], label = f"{who}, Imag", ls = 'dotted', color = colors[who])
        
        ax[bottop, 0].set_xlabel(r"$M_A$ in GeV")
        ax[bottop, 0].set_ylabel(r"$C^\mathrm{NLO}$")
for bottop in range(2):
    ax[bottop, 1].set_title(sq[bottop] + r", $\tan \beta = 40$")
    for who in people:
        data = finddataMA(sq[bottop],40, range(1000, 2801, 100), who)
        ax[bottop, 1].errorbar(data["MA"], data["ReAvg"],yerr = data["ReErr"], label = f"{who}, Real", ls = '--', color = colors[who])
        ax[bottop, 1].errorbar(data["MA"], data["ImAvg"],yerr = data["ImErr"], label = f"{who}, Imag", ls = 'dotted', color = colors[who])
        # ax[bottop, 1].legend()
        ax[bottop, 1].set_xlabel(r"$M_A$ in GeV")
        ax[bottop, 1].set_ylabel(r"$C^\mathrm{NLO}$")

for bottop in range(2):
    ax[bottop, 2].set_title(sq[bottop] + r", $M_A = 750$")
    for who in people:
        data = finddatatgb(sq[bottop], np.arange(2,10,0.5), 750, who)
        ax[bottop, 2].errorbar(data["tgb"], data["ReAvg"],yerr = data["ReErr"], label = f"{who}, Real", ls = '--', color = colors[who])
        ax[bottop, 2].errorbar(data["tgb"], data["ImAvg"],yerr = data["ImErr"], label = f"{who}, Imag", ls = 'dotted', color = colors[who])
        # ax[bottop, 2].legend()
        ax[bottop, 2].set_xlabel(r"$\tan \beta$")
        ax[bottop, 2].set_ylabel(r"$C^\mathrm{NLO}$")

ax[0,0].legend()



plt.savefig("summaryplots.pdf")

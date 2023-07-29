#!/usr/bin/env python
import glob, os
from sys import platform, float_info

if platform=="win32":
    file_sep = "\\"
else:
    file_sep = "/"

numerical_acc = False
# numerical_acc = True
alldata = {}
allstumps = []
f = "metabolomics_data_export"
# folders = ["metabolomics_data_export", "PRS_data_export", "immunomics_data_export/blood", "immunomics_data_export/infect", "Phenomics_data_export/charlson", "Phenomics_data_export/icd10"]
folders = ["metabolomics_data_export"]
#data_folder = os.path.join("html", "data")
data_folder = os.path.join("data_raw")

for f_short in folders:
    f = os.path.join(data_folder, f_short)
    for fn in sorted(glob.glob("{}/*.csv".format(f))):
        if "_info.csv" in fn:
            continue

        if "Interaction" in fn:     # worry about interaction later
            continue

        # print(f, fn)
        # continue

        # stump = fn.split(file_sep)[-1].replace(".csv", "").replace(".","")
        stump = fn.split(file_sep)[-1].replace(".csv", "").replace(".","").split("_")[-1]
        print(stump)
        allstumps.append(stump)
        with open(fn) as fin:
            header = None
            for l in fin:
                arr = l.strip().split(",")
                if header == None:
                    header = arr
                    continue
                try:
                    (x, y) = arr[0].split("_")
                except:
                    print(arr[0])
                for i in range(1, len(arr)):
                    key = "###".join((x,y,header[i]))
                    if not key in alldata:
                        alldata[key] = []
                    #alldata[key].append(str(round(float(arr[i]), 3)))
                    alldata[key].append( str( float(arr[i]) ) )

    print(allstumps)
    tokeep = set()
    allmet = set()
    for k in alldata.keys():
        (x, y, met) = k.split("###")
        allmet.add(met)
        if float(alldata[k][1]) < 0.05:     # adjPVal because sorted(glob())
            tokeep.add(met)

    print(len(tokeep))  # this will keep any that met threshold at any grid location!
    print(len(allmet))

    outdir = os.path.join(f, "long")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
        
    # fname = "alldata_v1.csv" if numerical_acc else "alldata.csv"
    # with open(os.path.join(outdir, fname), "w") as fout:
    #     fout.write("metabolite,x,y,%s\n" % ",".join(allstumps))
    #     for k in sorted(alldata.keys()):
    #         (x, y, met) = k.split("###")
    #         if met not in tokeep:
    #             continue
    #         if numerical_acc and float(alldata[k][1])==0:
    #             alldata[k][1] = str(float_info.min)
    #         fout.write("%s,%s,%s,%s\n" % (met, x, y, ",".join(alldata[k])))

    # make location summary files 29k and individual metabolite files for quick loading
    locdir = os.path.join(f, "location_summary")
    if not os.path.isdir(locdir):
        os.makedirs(locdir)
    longdir = os.path.join(f, "long_v2")
    if not os.path.isdir(longdir):
        os.makedirs(longdir)

    for k in sorted(alldata.keys()):
        (slice_num, coord_num, met) = k.split("###")
        metfile = os.path.join(longdir, "{}.csv".format(met))
        if not os.path.isfile(metfile):  # write header
            with open(metfile, 'w') as fout:
                fout.write("slice_num,coord_num,{}\n".format(",".join(allstumps)))
            fout.close()
        with open(metfile, 'a') as fout:
            vals = [slice_num, coord_num] +alldata[k]
            fout.write("{}\n".format(",".join(vals)))
        fout.close()

        locfile = os.path.join(locdir, "location_summary_{}_{}.csv".format(slice_num, coord_num))
        if not os.path.isfile(locfile):
            with open(locfile, "w") as fout:
                fout.write("met,{}\n".format(",".join(allstumps)))
            fout.close()
        with open(locfile, "a") as fout:
            vals = [met] +alldata[k]
            fout.write("{}\n".format(",".join(vals)))
        fout.close()

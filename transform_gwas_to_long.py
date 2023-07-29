import os, json
from sys import platform, float_info

if platform=="win32":
    file_sep = "\\"
else:
    file_sep = "/"

# DATA_FOLDER = "data"
DATA_FOLDER = "."
NUM_SLICES = 128
NUM_COORDS = 256
NUM_COORDS_VALID = 246
NUM_SLICES_VALID = 119
SIG_THRESHOLD_BASE = 0.05
SIG_THRESHOLD = 5*(10**-8) /NUM_SLICES_VALID /NUM_COORDS_VALID  # Vicki approximation of locus wide
SNIP_COUNTS = [842491, 913413, 774528, 796789, 706769, 733026, 639645, 600401, 468652, 553430, 539166, 524442, 398578,
               353816, 309546, 335728, 293820, 309206, 252330, 240886, 148455, 149315, 354574]
NUM_SNIPS = sum(SNIP_COUNTS)
SIG_THRESHOLD_STRICT = 1/(NUM_SNIPS*NUM_SLICES_VALID*NUM_COORDS_VALID)
KEY_SIG_COORDS = "sig_coords"
KEY_SIG_COORDS_STRICT = "sig_coords_strict"

CHROMOSOMES = list(range(1, 23)) + ["X"]
SLICE_NUMS = range(1, NUM_SLICES+1)
COORD_NUMS = range(1, NUM_COORDS+1)


def safe_snipID(snipID_raw):
    # TODO - fix for all windows and linux chars https://stackoverflow.com/questions/1976007/what-characters-are-forbidden-in-windows-and-linux-directory-names
    if ":" in snipID_raw:
        snipID = snipID_raw.split(":")[-1]
        return snipID
    else:
        return snipID_raw
    return snipID_raw


def process_chromosome(chromosomes=[22], slice_nums=[]):
    chr_dict = {}   # all chromosomes - combine at end instead!
    if len(slice_nums)==0:
        slice_nums2 = SLICE_NUMS
    else:
        slice_nums2 = slice_nums

    for chdx, chromosome in enumerate(chromosomes):     # for each chromosome or 1 target chromosome
        print("processing chrom", chdx, chromosome)

        # chr_folder = os.path.join(DATA_FOLDER, "chr{}".format(str(chromosome)))
        chr_folder = os.path.join(DATA_FOLDER, "chr{}".format(str(chromosome)), "sigSNPs")
        if not os.path.isdir(chr_folder):
            continue

        out_folder = os.path.join(chr_folder, "long")
        if not os.path.isdir(out_folder):
            os.makedirs(out_folder)

        chrom_long_file = os.path.join(chr_folder, "coord_summary_{}_long.csv".format(chromosome))
        if not os.path.isfile(chrom_long_file):
            with open(chrom_long_file, 'w') as fout:
                fout.write("snip,x,y,beta,pvalue,bonf\n")
            fout.close()

        snip_info_file = os.path.join(chr_folder, "snip_info_{}.csv".format(chromosome))
        if not os.path.isfile(snip_info_file):
            with open(snip_info_file, 'w') as fout:
                fout.write("snip,pos,A1,snip_index\n")
            fout.close()
        elif len(slice_nums)>0:
            snip_info_counter = 0
            with open(snip_info_file, "r") as fin:
                for l in fin:
                    l_toks = l.rstrip().split(",")
                    snip_info_counter +=1
                    if snip_info_counter==1:
                        snip_info_header = l_toks
                    else:
                        snipID,pos,A1,snip_index = l_toks
                        snip_info_dict[snipID] = 1
            fin.close()

        chr_snip_dict = {}
        snip_info_dict = {}
        for sdx, slice_num in enumerate(slice_nums2):    # for each slice, open and process slice file
            print("processing chrom,slice", chdx, chromosome, slice_num)

            fname = "slice{}_result.txt".format(slice_num)
            fpath = os.path.join(chr_folder, fname)

            if not os.path.isfile(fpath):
                print("slice data not found: ", slice_num, fpath)
                continue

            counter = 0
            with open(fpath, 'r') as fin:
                for l in fin:
                    counter+=1
                    l_toks = l.rstrip().split(",")

                    if counter==1:
                        header = l_toks
                    else:
                        cur_dict = dict(zip(header, l_toks))    # line dict
                        pos = cur_dict["POS"]
                        snipID_raw = cur_dict["ID"]
                        snipID = safe_snipID(snipID_raw)
                        snip_index = counter-1  # track which subfolder for netlify

                        if snipID not in snip_info_dict:
                            # snip_info_dict[snipID] = {"POS":pos, "chrom":chromosome, "A1":cur_dict["A1"]}
                            snip_info_dict[snipID] = 1  # so it only writes once
                            with open(snip_info_file, "a") as fout_info:
                                vals = [snipID, pos, cur_dict["A1"], snip_index]
                                fout_info.write("{}\n".format(",".join([str(x) for x in vals])))
                            fout_info.close()

                        if snipID not in chr_snip_dict:     # init
                            chr_snip_dict[snipID] = {KEY_SIG_COORDS:[], KEY_SIG_COORDS_STRICT:[]}

                        out_sub_folder = os.path.join(out_folder, "{}".format(int(snip_index//25000)))
                        if not os.path.isdir(out_sub_folder):
                            os.makedirs(out_sub_folder)
                        snip_out_file = os.path.join(out_sub_folder,  "{}_{}.txt".format(chromosome, snipID))
                        if not os.path.isfile(snip_out_file):   # write header if not exist
                            with open(snip_out_file, 'w') as fout_snip_file:
                                fout_snip_file.write("x,y,beta,pval\n")
                            fout_snip_file.close()

                        for cdx, coord_num in enumerate(COORD_NUMS):
                            coord_beta_key = "{}_{}_BETA".format(slice_num, coord_num)
                            coord_pval_key = "{}_{}_P".format(slice_num, coord_num)
                            if coord_beta_key in cur_dict and coord_pval_key in cur_dict:
                                coord_beta = cur_dict[coord_beta_key]
                                coord_pval = cur_dict[coord_pval_key]

                                with open(snip_out_file, "a") as fout_snip_file2:
                                    vals = [slice_num, coord_num, coord_beta, coord_pval]
                                    fout_snip_file2.write("{}\n".format(",".join([str(x) for x in vals])))
                                fout_snip_file2.close()

                                # if float(coord_pval)<SIG_THRESHOLD:
                                #     chr_snip_dict[snipID][KEY_SIG_COORDS].append((slice_num, coord_num))
                                # if float(coord_pval)<SIG_THRESHOLD_STRICT:
                                #     chr_snip_dict[snipID][KEY_SIG_COORDS_STRICT].append((slice_num, coord_num))
                                if float(coord_pval) < SIG_THRESHOLD:
                                    is_bonf = float(coord_pval)<SIG_THRESHOLD_STRICT
                                    with open(chrom_long_file, "a") as fout_chrome_file:
                                        vals = [snipID, slice_num, coord_num, coord_beta, coord_pval, int(is_bonf)]
                                        fout_chrome_file.write("{}\n".format(",".join([str(x) for x in vals])))
                                    fout_chrome_file.close()
            fin.close()
            print("processed chrom,slice", chdx, chromosome, slice_num)

        print("processed chrom", chdx, chromosome)
        # with open(os.path.join(out_folder, "snip_info.json"), "w") as fout:
        #     json.dump(snip_info_dict, fout)
        # fout.close()
        #
        # with open(os.path.join(out_folder, "chr_snip.json"), "w") as fout:
        #     json.dump(chr_snip_dict, fout)
        # fout.close()
    return


def process_chromosome_memory(chromosomes=[22], slice_nums=[]):
    chr_dict = {}   # all chromosomes - combine at end instead!
    if len(slice_nums)==0:
        slice_nums2 = SLICE_NUMS
    else:
        slice_nums2 = slice_nums

    for chdx, chromosome in enumerate(chromosomes):     # for each chromosome or 1 target chromosome
        print("processing chrom", chdx, chromosome)

        # chr_folder = os.path.join(DATA_FOLDER, "chr{}".format(str(chromosome)))
        chr_folder = os.path.join(DATA_FOLDER, "chr{}".format(str(chromosome)), "sigSNPs")
        if not os.path.isdir(chr_folder):
            continue

        out_folder = os.path.join(chr_folder, "long")
        if not os.path.isdir(out_folder):
            os.makedirs(out_folder)

        chrom_long_file = os.path.join(chr_folder, "coord_summary_{}_long.csv".format(chromosome))
        if not os.path.isfile(chrom_long_file):
            with open(chrom_long_file, 'w') as fout:
                fout.write("snip,chr,x,y,beta,pvalue,bonf\n")
            fout.close()

        snip_info_file = os.path.join(chr_folder, "snip_info_{}.csv".format(chromosome))
        if not os.path.isfile(snip_info_file):
            with open(snip_info_file, 'w') as fout:
                fout.write("snip,chr,pos,A1,snip_index\n")
            fout.close()

        snip_info_dict = {}
        snip_data_dict = {}     # collect across slices for each SNP
        for sdx, slice_num in enumerate(slice_nums2):    # for each slice, open and process slice file
            slice_coord_dict = {}  # collect SNP for each slice_num, coord_num and combine across chromosomes
            print("processing chrom,slice", chdx, chromosome, slice_num)

            fname = "slice{}_result.txt".format(slice_num)
            fpath = os.path.join(chr_folder, fname)

            if not os.path.isfile(fpath):
                print("slice data not found: ", slice_num, fpath)
                continue

            counter = 0
            with open(fpath, 'r') as fin:
                for l in fin:
                    counter+=1
                    l_toks = l.rstrip().split(",")

                    if counter==1:
                        header = l_toks
                    else:
                        cur_dict = dict(zip(header, l_toks))    # line dict
                        pos = cur_dict["POS"]
                        snipID_raw = cur_dict["ID"]
                        snipID = safe_snipID(snipID_raw)
                        snip_index = counter-1  # track which subfolder for netlify

                        if snipID not in snip_info_dict:
                            snip_info_dict[snipID] = {"SNP":snipID_raw, "POS":pos, "chr":chromosome, "A1":cur_dict["A1"]}

                        for cdx, coord_num in enumerate(COORD_NUMS):
                            coord_beta_key = "{}_{}_BETA".format(slice_num, coord_num)
                            coord_pval_key = "{}_{}_P".format(slice_num, coord_num)
                            if coord_beta_key in cur_dict and coord_pval_key in cur_dict:
                                coord_beta = cur_dict[coord_beta_key]
                                coord_pval = cur_dict[coord_pval_key]

                                is_bonf = float(coord_pval) < SIG_THRESHOLD_STRICT
                                vals = [snipID_raw, slice_num, coord_num, coord_beta, coord_pval, is_bonf]
                                if snipID not in snip_data_dict:
                                    snip_data_dict[snipID] = []
                                snip_data_dict[snipID].append(vals)

                                if coord_num not in slice_coord_dict:
                                    slice_coord_dict[coord_num] = []
                                slice_coord_dict[coord_num].append(vals)
                        # end for across COORD_NUMS
                    # end if
                # end readline
            fin.close()

            # location files
            loc_file = os.path.join(out_folder, "slice{}.txt".format(slice_num))
            with open(loc_file, 'w') as fout_loc:
                fout_loc.write("slice_num,x,y,beta,pval\n")
                for coord_num, coord_data in slice_coord_dict.items():
                    for d in coord_data:
                        vals = d
                        fout_loc.write("{}\n".format(",".join([str(x) for x in vals ])))
                fout_loc.close()

            print("processed chrom,slice", chdx, chromosome, slice_num)

        # out_sub_folder = os.path.join(out_folder, "{}".format(int(snip_index // 25000)))
        out_sub_folder = out_folder
        if not os.path.isdir(out_sub_folder):
            os.makedirs(out_sub_folder)
        # SNP files
        for snipID_raw, snipData in snip_data_dict.items():
            snipID = safe_snipID(snipID_raw)
            snipFile = os.path.join(out_sub_folder, "{}.csv".format(snipID))
            with open(snipFile, "w") as fout_snp:
                fout_snp.write("snip,slice_num,coord_num,beta,pval,is_bonf\n")
                for d in snipData:
                    fout_snp.write("{}\n".format(",".join([str(x) for x in d ])))
            fout_snp.close()

        print("processed chrom", chdx, chromosome)

        # output snip_info_dict
        snip_info_file = os.path.join(out_folder, "snip_info.csv".format(snipID))
        with open(snip_info_file, "w") as fout_info:
            fout_info.write("snip,snip_raw,chr,pos,A1\n")
            for snipID, info_data in snip_info_dict.items():
                vals = [snipID, info_data["SNP"], info_data["chr"], info_data["POS"], info_data["A1"]]
                fout_info.write("{}\n".format(",".join([str(x) for x in vals ])))
        fout_info.close()
    return


def convert_to_flat_csv(chr_num):   # generates summary file for chromosome and snip_info
    chr_long_folder = os.path.join(DATA_FOLDER, "chr{}".format(str(chr_num)), "long")
    chr_path = os.path.join(chr_long_folder, "chr_snip.json")

    chr_snip_dict = json.loads(open(chr_path).read())
    snip_keys = list(chr_snip_dict.keys())
    num_keys = len(snip_keys)

    keys_num_sig_loc = [(key, len(item[KEY_SIG_COORDS])) for key, item in chr_snip_dict.items()]
    keys_num_sig_loc_sorted = sorted(keys_num_sig_loc, key=lambda x: x[1], reverse=True)

    # # sanity check
    # print(chr_snip_dict[keys_num_sig_loc_sorted[0][0]])
    # print(list(set(chr_snip_dict[keys_num_sig_loc_sorted[0][0]][KEY_SIG_COORDS])))
    # # how many really significant
    # num_significant_bonf = []
    # for x in keys_num_sig_loc_sorted:
    #     if x[1]>0:
    #         num_significant_bonf.append(x)

    csv_path = chr_path.replace(".json", ".csv")
    with open(csv_path, "w") as fout:
        fout.write("snip,chr,num_sig,num_sig_strict\n")   # header

        for key_count in keys_num_sig_loc_sorted:
            snipID, num_sig = key_count
            vals = [snipID, chr_num, num_sig,
                    len(chr_snip_dict[snipID][KEY_SIG_COORDS_STRICT] if KEY_SIG_COORDS_STRICT in chr_snip_dict[snipID] else [])]
            fout.write("{}\n".format(",".join([str(x) for x in vals])))
    fout.close()

    snip_info_json_path = os.path.join(chr_long_folder, "snip_info.json")
    snip_info_csv_path = snip_info_json_path.replace(".json", ".csv")
    snip_info_dict = json.loads(open(snip_info_json_path).read())
    with open(snip_info_csv_path, "w") as fout:
        fout.write("snip,chr,pos,A1\n")     # header

        for snipID, snipData in snip_info_dict.items():
            vals = [snipID, snipData['chrom'], snipData['POS'], snipData['A1']]
            fout.write("{}\n".format(",".join([str(x) for x in vals])))
    fout.close()

    return


def summary_by_location():
    loc_summary_folder = os.path.join(DATA_FOLDER, "location_summary")
    if not os.path.isdir(loc_summary_folder):
        os.makedirs(loc_summary_folder)

    # for chdx, chr_num in enumerate(CHROMOSOMES):
    for chdx, chr_num in enumerate([22]):
        chr_long_folder = os.path.join(DATA_FOLDER, "chr{}".format(str(chr_num)), "long")
        chr_path = os.path.join(chr_long_folder, "chr_snip.json")

        if not os.path.isfile(chr_path):
            print("missing chromosome summary file:", chr_num, chr_path)
            continue
        chr_snip_dict = json.loads(open(chr_path).read())
        # snip_keys = list(chr_snip_dict.keys())
        # num_keys = len(snip_keys)

        snip_info_path = os.path.join(chr_long_folder, "snip_info.json")
        snip_info_dict = json.loads(open(snip_info_path).read())

        keys_num_sig_loc = [(key, len(item[KEY_SIG_COORDS])) for key, item in chr_snip_dict.items()]
        keys_num_sig_loc_sorted = sorted(keys_num_sig_loc, key=lambda x: x[1], reverse=True)
        snips_significant_bonf = []
        for x in keys_num_sig_loc_sorted:
            if x[1]>0:
                snips_significant_bonf.append(x[0])

        for snipID_raw in snips_significant_bonf:
            snipID = safe_snipID(snipID_raw)
            snip_path = os.path.join(chr_long_folder, "{}_{}.txt".format(chr_num, snipID))
            if not os.path.isfile(snip_path):
                print("missing snip file:", chr_num, snip_path)
                continue
            snip_dict = parse_snip_file(snip_path)

            for sdx, slice_num in enumerate(SLICE_NUMS):
                for cdx, coord_num in enumerate(COORD_NUMS):
                    summary_file_path = os.path.join(loc_summary_folder, "location_summary_{}_{}.csv".format(slice_num, coord_num))
                    if not os.path.isfile(summary_file_path):
                        with open(summary_file_path, "w") as fout:
                            fout.write("snip,chr,pos,beta,pvalue\n")
                        fout.close()

                    xy = "{},{}".format(slice_num, coord_num)
                    if xy not in snip_dict:
                        print("no data for coordinate:", slice_num, coord_num, snipID)
                        continue

                    cur_data = snip_dict[xy]
                    snipPOS = snip_info_dict[snipID]["POS"]
                    with open(summary_file_path, "a") as fout:
                        if cur_data[1] < SIG_THRESHOLD:
                            vals = [snipID, chr_num, snipPOS] +cur_data
                            fout.write("{}\n".format(",".join([str(x) for x in vals])))
        fout.close()
    return


def parse_snip_file(snip_path):
    snip_dict = {}

    counter = 0
    with open(snip_path, "r") as fin:
        for l in fin:
            counter+=1
            l_toks = l.rstrip().split(",")
            if counter==1:
                header = l_toks
            else:
                cur_dict = dict(zip(header, l_toks))
                x = cur_dict["x"]
                y = cur_dict["y"]
                xy = "{},{}".format(x,y)
                snip_dict[xy] = [float(cur_dict["beta"]), float(cur_dict["pval"])]
    fin.close()

    return snip_dict


if __name__ == "__main__":
    # for parallel
    import sys
    print(sys.argv)
    if len(sys.argv) > 1:
        chr_num = sys.argv[1]
    else:  # default values
        chr_num = 22

    chromomsomes = [chr_num]
    # process_chromosome(chromomsomes)
    process_chromosome_memory(chromomsomes)
    #process_chromosome_memory(CHROMOSOMES[:-1])
    # convert_to_flat_csv(chr_num)
    # summary_by_location()

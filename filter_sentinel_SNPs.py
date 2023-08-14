import os, json
from sys import platform, float_info

from transform_gwas_to_long import file_sep, DATA_FOLDER, NUM_COORDS, NUM_COORDS_VALID, NUM_SLICES, NUM_SLICES_VALID, \
    SNIP_COUNTS, NUM_SNIPS, CHROMOSOMES, SLICE_NUMS, COORD_NUMS, safe_snipID


def parse_sentinelSNPs():
    loci_summary_file = os.path.join(DATA_FOLDER, "gwSigLociSummary.csv")
    loci_out_json = os.path.join(DATA_FOLDER, "gwSigLociSummary.json")

    if os.path.isfile(loci_out_json):
        loci_dict = json.loads(open(loci_out_json).read())
        return loci_dict

    loci_dict = {}
    counter = 0
    with open(loci_summary_file, "r") as fin:
        for l in fin:
            counter +=1
            l_toks = l.rstrip().split(",")

            if counter==1:
                header= l_toks
            else:
                cur_dict = dict(zip(header, l_toks))
                locusID = cur_dict["locusID"]
                loci_dict[locusID] = cur_dict
    fin.close()

    with open(loci_out_json, "w") as fout:
        json.dump(loci_dict, fout)
    fout.close()
    return loci_dict


def parse_lociSNPs():
    # loci_dict = parse_sentinelSNPs()

    locSNPs_file = os.path.join(DATA_FOLDER, "gwSigLociSNPs.csv")
    locSNPIDs_out = os.path.join(DATA_FOLDER, "gwSigLociSNP_IDs.csv")

    counter = 0
    with open(locSNPIDs_out, "w") as fout:
        with open(locSNPs_file, "r") as fin:
            for l in fin:
                counter +=1
                l_toks = l.rstrip().split(",")

                if counter==1:
                    header = l_toks
                else:
                    cur_dict = dict(zip(header, l_toks))
                    # locusID = cur_dict["locusID"]
                    # sentinelSNP = cur_dict["sentinelSNPID"]
                    snipID = cur_dict["locusSNPID"]
                    fout.write("{}\n".format(snipID))
        fin.close()
    fout.close()
    return


def filter_SNPs(chromosomes=CHROMOSOMES):
    locSNPIDs_out = os.path.join(DATA_FOLDER, "gwSigLociSNP_IDs.csv")
    locSNPIDs = []
    with open(locSNPIDs_out, "r") as fin:
        for l in fin:
            locSNPIDs.append(l.rstrip())
    fin.close()
    locSNPIDs_set = set(locSNPIDs)

    for cdx, chrom in enumerate(chromosomes):
        chr_folder = os.path.join(DATA_FOLDER, "chr{}".format(chrom))
        outdir = os.path.join(chr_folder, 'sigSNPs')
        if not os.path.isdir(outdir):
            os.makedirs(outdir)

        for sdx, slice_num in enumerate(SLICE_NUMS):
            print("processing:", cdx, chrom, sdx, slice_num)

            fname = "slice{}_result.txt".format(slice_num)
            fpath = os.path.join(chr_folder, fname)
            outpath = os.path.join(outdir, fname)

            if not os.path.isfile(fpath):
                print("slice data not found: ", slice_num, fpath)
                continue

            with open(outpath, "w") as fout:
                counter = 0
                with open(fpath, 'r') as fin:
                    for l in fin:
                        counter+=1
                        l_toks = l.rstrip().split(",")

                        if counter==1:
                            header = l_toks
                            fout.write("{}\n".format(l))
                        else:
                            cur_dict = dict(zip(header, l_toks))    # line dict
                            snipID_raw = cur_dict["ID"]
                            # pos = cur_dict["POS"]
                            # snipID = safe_snipID(snipID_raw)
                            # snip_index = counter-1  # track which subfolder for netlify
                            if snipID_raw in locSNPIDs_set:
                                fout.write("{}".format(l))
                fin.close()
            fout.close()
            print("processed:", cdx, chrom, sdx, slice_num)
    return


def make_combined_SNPs_info_file():
    data_folder_html = os.path.join("html", "data", "gwas")
    data_folder_raw = os.path.join("data_raw", "gwas")
    # locSNPIDs_file = os.path.join(data_folder_html, "gwSigLociSNP_IDs.csv")
    # locSNPIDs = []
    # with open(locSNPIDs_file, "r") as fin:
    #     for l in fin:
    #         locSNPIDs.append(l.rstrip())
    # fin.close()
    # locSNPIDs_set = set(locSNPIDs)  # NB - non-unique by just rs_name!

    sentinelSNPs_file = os.path.join(data_folder_html, "gwSigLociSummary.csv")
    sentinelDict = {}
    with open(sentinelSNPs_file, "r") as fin:
        counter = 0
        for l in fin:
            counter +=1
            l_toks = l.rstrip().split(",")
            if counter==1:
                header=l_toks
            else:
                cur_dict = dict(zip(header, l_toks))
                curID = cur_dict["ID"]
                sentinelDict[curID] = cur_dict
    fin.close()

    locSNP_file = os.path.join(data_folder_html, "gwSigLociSNPs.csv")
    SNP_IDs = []
    chr_SNP_keys = []
    sigSNPs_dict = {}
    with open(locSNP_file, "r") as fin:
        counter = 0
        for l in fin:
            counter +=1
            l_toks = l.rstrip().split(",")
            if counter==1:
                header=l_toks  # locusID,sentinelSNPID,locusSNPID,r2withSentinel
            else:
                cur_dict = dict(zip(header, l_toks))
                curID = cur_dict["locusSNPID"]
                SNP_IDs.append(curID)
                chr_SNP_keys.append("{}; {}".format(sentinelDict[cur_dict["sentinelSNPID"]]["CHR"], curID))
                sigSNPs_dict[curID] = cur_dict
    fin.close()
    locSNPIDs_set = set(sigSNPs_dict.keys())
    # sanity checks
    print([len(SNP_IDs), len(set(SNP_IDs))])
    print([len(chr_SNP_keys), len(set(chr_SNP_keys))])

    # go through slice1_result.txt for each chromosome
    sigSNPS_outfile = os.path.join(data_folder_html, "info_sigSNPs.csv")
    SNPS_outfile = os.path.join(data_folder_html, "info_SNPs.csv")
    SNPS_outfile = os.path.join(data_folder_html, "info_SNPs_v2.csv")
    with open(sigSNPS_outfile, "w") as fout_sig:
        fout_sig.write("SNP,chr,pos,A1,num_sig,bonf,sentinel\n")
        with open(SNPS_outfile, "w") as fout:
            fout.write("SNP,chr,pos\n")
            # fout.write("SNP,chr\n")    # bare minimum

            for cdx, chromosome in enumerate(CHROMOSOMES):
                fpath = os.path.join(data_folder_raw, "chr{}".format(chromosome), "slice1_result.txt")
                print("processing ", cdx, chromosome, fpath)

                with open(fpath, "r") as fin:
                    counter = 0
                    for l in fin:
                        counter += 1
                        l_toks = l.rstrip().split(",")

                        if counter==1:
                            header = l_toks
                        else:
                            cur_dict = dict(zip(header, l_toks))
                            curID = cur_dict["ID"]
                            if curID in locSNPIDs_set:
                                sentinel_data = sentinelDict[sigSNPs_dict[curID]["sentinelSNPID"]]
                                is_sentinel = int(curID in sentinelDict)
                                out_vals = [curID, chromosome, cur_dict["POS"], cur_dict["A1"],
                                            sentinel_data["nPixelsLocus"], sentinel_data["BonferroniSig"], is_sentinel]
                                fout_sig.write("{}\n".format(",".join([str(x) for x in out_vals])))
                            else:
                                # out_vals = [curID, chromosome]  # bare minimum
                                out_vals = [curID, chromosome, cur_dict["POS"]]  # bare minimum
                                fout.write("{}\n".format(",".join([str(x) for x in out_vals])))
                fin.close()
        fout.close()
    fout_sig.close()

    return


def make_combined_location_files():
    data_folder_html = os.path.join("html", "data", "gwas")
    location_folder = os.path.join(data_folder_html, "location_summary")
    if not os.path.isdir(location_folder):
        os.makedirs(location_folder)

    for cdx, chromosome in enumerate(CHROMOSOMES):
        slice_coord_dict = {}
        chr_folder = os.path.join(data_folder_html, "chr{}".format(chromosome), "sigSNPs", "long")
        for sdx, slice_num in enumerate(range(1,NUM_SLICES_VALID+1)):
            fpath = os.path.join(chr_folder, "slice{}.txt".format(slice_num))
            print("processing: ", cdx, chromosome, sdx, slice_num, fpath)

            with open(fpath, "r") as fin:
                counter = 0
                for l in fin:
                    counter+=1
                    l_toks = l.rstrip().split(",")

                    if counter==1:
                        header = l_toks
                        header = ["SNP","slice_num","coord_num","beta","pval","Bonf"]
                    else:
                        cur_dict = dict(zip(header, l_toks))
                        slice_coord = "{},{}".format(cur_dict["slice_num"],cur_dict["coord_num"])
                        if slice_coord not in slice_coord_dict:
                            slice_coord_dict[slice_coord] = []
                        cur_vals = [cur_dict["SNP"],chromosome,cur_dict["beta"],cur_dict["pval"]]
                        slice_coord_dict[slice_coord].append(cur_vals)
            fin.close()

        # write location_summaries at chromosome level to avoid RAM issues
        for slice_coord, sdata in slice_coord_dict.items():
            stoks = slice_coord.split(",")
            slice_num, coord_num = stoks
            outpath = os.path.join(location_folder, "loc_summary_{}_{}.csv".format(slice_num, coord_num))
            print("writing: ", cdx, chromosome, slice_num, coord_num, outpath)

            if not os.path.isfile(outpath):     # add header if first write
                with open(outpath, "w") as fout:
                    fout.write("SNP,chr,beta,pval\n")
                fout.close()

            with open(outpath, "a") as fout:    # otherwise append
                for ddx, d in enumerate(sdata):
                    fout.write("{}\n".format(",".join([str(x) for x in d])))
            fout.close()
    return


if __name__ == "__main__":
    # parse_sentinelSNPs()
    # parse_lociSNPs()
    # filter_SNPs()
    make_combined_SNPs_info_file()
    # make_combined_location_files()
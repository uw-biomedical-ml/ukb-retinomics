import os, json
from sys import platform, float_info

from .transform_gwas_to_long import file_sep, DATA_FOLDER, NUM_COORDS, NUM_COORDS_VALID, NUM_SLICES, NUM_SLICES_VALID, \
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


if __name__ == "__main__":
    # parse_sentinelSNPs()
    # parse_lociSNPs()
    filter_SNPs()
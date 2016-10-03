#!/usr/bin/env python3

import sys, os, re
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'modules'))
import pop2df as ptd
import pandas as pd
import feather
import argparse
import inspect
import glob
import zipfile

parser = argparse.ArgumentParser(
    formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    fromfile_prefix_chars = '@'
    )

parser.add_argument(
    "--regex", 
    type = str,
    help = "Regular expression to filter files (note: this is applied first)",
    default = "gen_10"
    )

parser.add_argument(
    "--prefix",
    type = str,
    help = "prefix for directories",
    default = None
    )

parser.add_argument(
    "--separate",
    action = "store_true",
    help = "if this flag is signalled, each file will be processed separately. This overrides group_by"
    )

parser.add_argument(
    "--group_by",
    type = str,
    help = "by what elements should the feather files be grouped?",
    default = None
    )

parser.add_argument(
    "-z", "--zip",
    help = "should each feather file be zipped?",
    action = "store_true"
    )

parser.add_argument(
    "-o", "--out",
    help = "directory in which to place the feather output",
    default = "pillow")

parser.add_argument(
    "-s", "--snp",
    help = "Is the data SNP data? Setting this flag will process the genotypes as SNP data, summing the minor allele",
    action = "store_true")


def ruffle(pops, snp, out, fname, zip):
    '''
    Write feather files to a directory
    '''
    df = ptd.pops2df(pops, snp)
    outname = out + "/" + fname + ".feather"
    print("writing file " + outname)
    feather.write_dataframe(df, outname)
    if zip:
        zipup(outname)


def zipup(fname):
    '''
    zip a single file n disk
    '''
    zipf = zipfile.ZipFile(fname + ".zip", "w", zipfile.ZIP_DEFLATED)
    zipf.write(fname)
    zipf.close()


if __name__ == '__main__':

    pars = parser.parse_args()

    dirs = glob.glob(pars.prefix + "*")
    if (len(dirs) < 1):
        print("No directories with the prefix " + pars.prefix)
        sys.exit()
    if not os.path.isdir(pars.out):
        os.mkdir(pars.out)
    for d in dirs:
        print("processing directory " + d)
        pops = [d + "/" + x for x in os.listdir(d)]
        fname  = "dir_" + d
        if (len(pops) < 1):
            print("No files found in " + d)
            continue
        if pars.regex is not None:
            fname += "_" + pars.regex
            pops   = ptd.get_field(pops, pars.regex)
        if pars.separate:
            for p in pops:
                tempfname = fname + os.path.basename(p)
                ruffle([p], pars.snp, pars.out, tempfname, pars.zip)
        elif pars.group_by is not None:
            finder  = re.compile('^.+?(' + pars.group_by + '_[^_]+).*?\.pop$')
            matches = set([re.sub(finder, r'\1', i) for i in pops])
            for g in matches:
                tempfname = fname + "_group_" + g
                print(tempfname)
                pop_group = ptd.get_field(pops, g)
                ruffle(pop_group, pars.snp, pars.out, tempfname, pars.zip)
        else:
            ruffle(pops, pars.snp, pars.out, fname, pars.zip)



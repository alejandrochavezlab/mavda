# -*- coding: utf-8 -*-

"""
--------------------------------------------------------------------------------
Copyright 2023 Benjamin Alexander Albert [Alejandro Chavez Lab at UCSD]
All Rights Reserved
MAVDA Academic License
callhits.py
--------------------------------------------------------------------------------
"""

import os

import argparse

import numpy as np

import pandas as pd

from hitcaller import HitCaller

from platelayout import PlateLayout


def simsubsample(df, frac):
    """
    Simulate a lower sequencing depth by subsampling to
    investigate the potential multiplexing capacity.
    """

    subdf = list()
    for well,welldf in df.groupby("well"):

        welldf.set_index(["gene","barcode"], inplace=True)
        wellsum = welldf.counts.sum()
        choices = list()
        for idx,row in welldf.iterrows():
            choices.extend([idx]*row.counts)
        choices = np.array(choices)
        choices = choices[
            np.random.choice(len(choices), int(frac*wellsum))]
        choices = pd.DataFrame(choices).value_counts()
        welldf["counts"] = 0
        for idx,cnt in choices.to_dict().items():
            welldf.loc[idx,"counts"] = cnt
        subdf.append(welldf.reset_index())

    return pd.concat(subdf, ignore_index=True)


def parseArgs():

    parser = argparse.ArgumentParser(
        "Calculate magratios for all runs in the datadir",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--subsamp",
        type     = int,
        required = False,
        default  = 0,
        help     = (
            "nonzero to conduct an exploratory study in which we subsample"
            "\nNGS reads in order to approximate how much further"
            "\nwe could multiplex the assay. Otherwise, no subsampling."
            "\n(default=0)"))

    return parser.parse_args()


def run(savedir, keepfrac=100):

    datadir = "data"
    
    countsdir = os.path.join(datadir, "counts")
    oddir     = os.path.join(datadir, "od")
    
    for run in sorted(os.listdir(countsdir)):
    
        print('-'*80)
        print("processing run {}".format(run))
    
        rundir = os.path.join(countsdir, run)
    
        od = pd.read_csv(
            os.path.join(oddir, "{}.csv".format(run)))
    
        hitlist = list()
        for countsfile in sorted(os.listdir(rundir)):
            try:
                print("\tprocessing {}".format(countsfile))
                plate = countsfile[:countsfile.rfind('_')]
                plateod = {row.well:row.OD595
                    for _,row in od[od.plate == plate].iterrows()}
                df = pd.read_csv(os.path.join(rundir, countsfile))
                if run == "221126":
                    dmsoWells   = [x for x in df.well if "DMSO" in x]
                    lopControls = {x:["HIV-1","HIV-2"] \
                        for x in df.well if "Lopinavir" in x}
                else:
                    df.well = df.well.apply(
                        lambda x: "{}{:02d}".format(x[0],int(x[1:])))
                    plateLayout = PlateLayout(plate)
                    dmsoWells   = plateLayout.getNegCtrl()
                    lopControls = plateLayout.getLopCtrl()
    
                if keepfrac < 100:
                    df = simsubsample(df, frac=keepfrac/100)
    
                hitcaller = HitCaller(
                    wells       = df.well,
                    genes       = df.gene,
                    barcodes    = df.barcode,
                    counts      = df.counts,
                    od          = plateod,
                    dmsoWells   = dmsoWells,
                    posControls = lopControls)
                hitcaller.preprocess(mincounts=30000*keepfrac/100)
                hitcaller.findDifferentialCounts()
                hits = hitcaller.hits.copy()
                hits.insert(0, "plate", plate)
                hitlist.append(hits)
            except Exception as e:
                print("\tERROR: {}".format(e))
        hitlist = pd.concat(hitlist)
        if not os.path.exists(savedir):
            os.makedirs(savedir)
        hitlist.to_csv(
            os.path.join(savedir,"{}.csv".format(run)),
            index=False)


def main():

    args = parseArgs()

    np.random.seed(1)

    if args.subsamp:
        for keepfrac in range(1,100):
            savedir = os.path.join(
                "subsamp",
                "hits_keepfrac{}".format(keepfrac))
            run(savedir=savedir, keepfrac=keepfrac)
    else:
        run(savedir="hits")


if __name__ == "__main__":
    main()

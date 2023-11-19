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

import pandas as pd

from hitcaller import HitCaller

from platelayout import PlateLayout


def main():

    savedir = "hits"
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

                hitcaller = HitCaller(
                    wells       = df.well,
                    genes       = df.gene,
                    barcodes    = df.barcode,
                    counts      = df.counts,
                    od          = plateod,
                    dmsoWells   = dmsoWells,
                    posControls = lopControls)
                hitcaller.preprocess()
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


if __name__ == "__main__":
    main()

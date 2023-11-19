# -*- coding: utf-8 -*-

"""
--------------------------------------------------------------------------------
Copyright 2023 Benjamin Alexander Albert [Alejandro Chavez Lab at UCSD]
All Rights Reserved
MAVDA Academic License
hitcaller.py
--------------------------------------------------------------------------------
"""

import os

from typing import Tuple, Iterable

import numpy as np

import pandas as pd

from scipy import stats


class HitCaller:

    def __init__(
            self,
            wells     : Iterable,
            genes     : Iterable,
            barcodes  : Iterable,
            counts    : Iterable,
            od        : dict,
            dmsoWells : list = [
                "A01","A12",
                "B01","B12",
                "C01","C12",
                "D12",
                "E01",
                "F01","F12",
                "G01","G12",
                "H01","H12"],
            posControls : dict = {
                "D01":["HIV-1", "HIV-2"],
                "E12":["HIV-1", "HIV-2"]},
            negControlGene : str = "EYFP"):
        """
        Create a HitCaller for a given plate counts table and
        corresponding optical density values.


        Parameters
        ----------

        wells : Iterable[str]
            Iterable of n well names

        genes : Iterable[str]
            Iterable of n gene names

        barcodes : Iterable[str]
            Iterable of n barcodes

        counts : Iterable[int]
            Iterable of n counts

        od : dict[str] -> float
            dict mapping each well to an OD float.
            Each unique well in the wells parameter must
            exist in this od dict keys.

        dmsoWells : list[str]
            list of wells that are only incubated with DMSO.

        posControls : dict[str] -> Iterable[str]
            dict mapping positive control wells to genes.

        negCtrlGene : str
            A gene that should not affect cell growth.
        """
        
        # sanity check that we received equal numbers of:
        # wells, genes, barcodes, and counts
        if len(wells) != len(genes) or \
           len(wells) != len(barcodes) or \
           len(wells) != len(counts):
               raise ValueError(
                    "Unequal number of wells, genes, barcodes, or counts")

        # save controls
        self.dmsoWells      = dmsoWells
        self.posControls    = posControls
        self.negControlGene = negControlGene

        # build counts table
        dmso = [(w in dmsoWells) for w in wells]
        self.df = pd.DataFrame(
            {"well"   :wells,
             "gene"   :genes,
             "barcode":barcodes,
             "dmso"   :dmso,
             "counts" :counts})
        self.df.well    = self.df.well.astype(str)
        self.df.gene    = self.df.gene.astype(str)
        self.df.barcode = self.df.barcode.astype(str)
        self.df.dmso    = self.df.dmso.astype(bool)
        self.df.counts  = self.df.counts.astype(int)

        # sanity check that each well has an OD value
        for well in self.df.well.unique():
            if well not in od.keys():
                raise ValueError(
                    "No OD value found for well {}".format(well))

        # build OD table
        self.od = pd.DataFrame.from_dict(
            data=od,
            columns=["od"],
            orient="index")
        self.od = self.od.astype(float)


        # these values remain none until we execute preprocess
        self.samples  = None
        self.dmsodf   = None
        self.drugdict = None
        self._preprocessed = False


        # these values remain None until we execute testDifferentialCounts
        self.hits    = None
        self.pvalpos = None


    def preprocess(
        self,
        mincounts : int   = 30000,
        nsamp     : int   = 1000,
        odabsmin  : float = 0.40,
        odstdmin  : float = 0.50):
        """

        Preprocess our counts table and OD values to:
        1. high-pass filter wells based on total counts
        2. remove wells that do not pass our OD filter
        3. sample wells for the Monte Carlo method differential counts method
        4. pivot and split the DMSO and drug wells

        mincounts : int
            Minimum number of counts to keep a well.
            See HitCaller._filterCounts

        nsamp : 
            Number of gene barcode replicates to sample per well.
            See Hitcaller._sampleCounts
        """

        # preprocess counts table before calling hits
        self._filterCounts(mincounts)
        self._filterOD(odabsmin=odabsmin, odstdmin=odstdmin)
        self.samples = self._sampleCounts(n=nsamp)
        self.dmsodf, self.drugdict = self._pivotSamples()
        self._preprocessed = True


    def _filterCounts(
        self,
        mincounts : int):
    
        """
        Remove wells from counts table if the well sum is below a threshold.
    
        Parameters
        ----------
    
        mincounts : int
            minimum (inclusive) number of reads to keep the well

        Returns
        -------
        None
        """
    
        for well in self.df.well.unique():
            wellfilter = (self.df.well == well)
            if self.df[wellfilter].counts.sum() < mincounts:
                self.df = self.df[~wellfilter]
    
    
    def _filterOD(
        self,
        odabsmin = float,
        odstdmin = float):

        odmean = self.od.od.mean()
        odstd  = self.od.od.std()

        odmin = min([odabsmin, odmean - (odstdmin*odstd)])

        rmwells = [w for w,row in self.od.iterrows() if row.od < odmin]

        self.df = self.df[self.df.well.apply(lambda w: w not in rmwells)]
    
    
    def _sampleCounts(self, n : int) -> Tuple[dict, pd.DataFrame]:
    
        """
        Generate randomly sampled cell counts from the counts table.

        We perform differential counts analysis as a Monte Carlo method.
        
        Parameters
        ----------
        
        n : int
            Number of times to sample protease counts per plate per well.
        
        Returns
        -------
        
        DataFrame containing sampled protease counts
        """
        
        genes = sorted(self.df.gene.unique())
        
        samples = list()
        for well, welldf in self.df.groupby("well"):

            # all genes must be present in all plate wells
            if not (genes == sorted(welldf.gene.unique())):
                raise ValueError(
                    "plate {} well {} does not contain the correct gene set" \
                        .format(plate, well))

            # dmso values must all be true or false
            if (not welldf.dmso.all()) and welldf.dmso.any():
                raise ValueError(
                    "plate {} well {} contains mixed DMSO values" \
                        .format(plate, well))
                    
            wellsamp = welldf                       \
                .groupby("gene", group_keys=False) \
                .apply(
                    lambda df: df.sample(
                        n=n,
                        replace=True,
                        random_state=1))
            wellsamp.insert(1, "rep", [(x % n) for x in range(len(wellsamp))])
            samples.append(wellsamp)
        return pd.concat(samples)
    

    def _pivotSamples(self) -> Tuple[pd.DataFrame, dict]:
    
        """
        Pivot samples DataFrame to have a gene index and rep columns.
    
        Parameters
        ----------
    
        samples : DataFrame
            The counts table from sampleCounts
    
        Returns
        -------
    
        tuple with elements:
        [0] dmso pivot table with all dmso wells concatenated
        [1] drug dict where each key is a well and each value is a pivot table
        """
    
        dmso = list()
        drug = dict()
    
        for well, welldf in self.samples.groupby("well"):
            pivot = welldf.pivot(
                index   = "gene",
                columns = "rep",
                values  = "counts")
            if welldf.dmso.any():
                dmso.append(pivot)
            else:
                drug[well] = pivot
        dmso = pd.concat(dmso, axis=1)
        return dmso, drug
    
    
    def _testGeneDifferentialCounts(
        self,
        dmsodf   : pd.DataFrame,
        drugdf   : pd.DataFrame,
        gene     : str,
        testreps : int,
        sampreps : int) -> float:
    
        """
        Determine whether a given gene is differentially increased
        from a drug treatment relative to the DMSO control.
    
        To do this, we first normalize the target gene counts
        by EYFP counts for both DMSO and drug wells.
    
        Then, if the mean of the drug well counts is less
        than the DMSO counts, we return a p-value of 1
        because we know that, if there is any differential count,
        it would not be in favor of drug rescue.
    
        Next, we bootstrap independent t-tests for a given number
        of repitions and using a given number of samples per iteration.
        The number of samples drawn from DMSO wells is proportional to
        the number drawn from the drug well, where the multiplicative
        factor is simply the number of DMSO wells. For example,
        if there are 14 DMSO wells, we draw 14 times more samples
        from DMSO than the number of drug well samples.
    
        Lastly, we the geometric mean of the t-test p-values.
    
        Parameters
        ----------
        
        dmsodf : DataFrame
            The DMSO DataFrame returned from _pivotSamples
        
        drugdf : DataFrame
            One of the values in the drug dict returned from _pivotSamples
    
        gene : str
            A gene name found in both the DMSO and drug DataFrames
    
        testreps : int
            Number of iterations to bootstrap the independent t-test
    
        sampreps : int
            Number of samples to draw from the
            drug DataFrame per bootstrap iteration.
            The number of drawsfrom the DMSO DataFrame is equal
            to this number multiplied by the number of DMSO wells.
    
        Returns
        -------
    
        The p-value for the test of differentially increased
        gene counts in the drug well compared to the DMSO wells.
        """
    
        def _sampleCounts(counts, n):
            return counts.sample(
                n=n*testreps,
                replace=True,
                random_state=1) \
            .to_numpy() \
            .reshape(-1, testreps)
    
        dmsocounts = dmsodf.loc[gene] / dmsodf.loc[self.negControlGene]
        drugcounts = drugdf.loc[gene] / drugdf.loc[self.negControlGene]
    
        if (drugcounts.mean() <= dmsocounts.mean()):
            return 1
    
        drugsamps = _sampleCounts(
            counts=drugcounts,
            n=sampreps)
        dmsosamps = _sampleCounts(
            counts=dmsocounts,
            n=sampreps*(len(dmsocounts)//len(drugcounts)))
        return stats.gmean(stats.ttest_ind(drugsamps, dmsosamps)[1])


    @staticmethod
    def _pvalNormalization(
        df : pd.DataFrame,
        geneNormPercentile : float,
        kinases : list) -> pd.DataFrame:
        
        """
        We apply two normalization factors to our p-values:
        gene normalization and kinase normalization.
    
        Gene normailzation applies across wells.

            We divide the p-values for a given gene by the p-value
            at the geneNormPercentile percentile across all wells.
            Most p-values should be high, so this should not affect most values.
            A low p-value at a relatively high percentile would indicate
            erroneous values, so we correct for such error.
    
        Kinase normalization applies to each well individually.

            A drug is assumed to either rescue kinases, non-kinases, or neither.
            If a drug is found to have low p-values in both a kinase well and
            a non-kinase well, then this is likely erroneous.
    
            Given the minimum kinase and non-kinase p-values of a well,
            p_k and p_n, we normalize the well p-values by max(p_k, p_n).
            Therefore, legitimate well p-values are maintained because either
            p_k, p_n, or both are high, and erroneous p-values are increased if
            both p_k and p_n are low.
    
        Parameters
        ----------
        
        df : DataFrame
            The DataFrame returned from testPlateDifferentialCounts.

        geneNormPercentile : float
            The percentile of the gene-wise p-values to use for normalization.
    
        kinases : list
            A list of gene names indicating all kinase columns.
    
        Returns
        -------
    
        DataFrame with with an append "pnorm" column of normalized p-values.
        """

        df["genenorm"]   = None
        df["kinasenorm"] = None

        for gene in df.gene.unique():
            genefilter = (df.gene==gene)
            df.loc[genefilter, "genenorm"] = \
                np.percentile(df[genefilter].p, geneNormPercentile)
    
        for well, welldf in df.groupby("well"):
            kinasefilter = (welldf.gene.apply(lambda x: x in kinases))
            norm = np.max([
                welldf[ kinasefilter].p.min(),
                welldf[~kinasefilter].p.min()])
            df.loc[(df.well==well), "kinasenorm"] = norm

        df["pnorm"] = (df.p / (df.genenorm * df.kinasenorm)).astype(float)
    
        return df
    
    
    def findDifferentialCounts(
        self,
        testreps : int = 1000,
        sampreps : int = 5,
        geneNormPercentile : float = 25,
        kinases  : list = ["YES1","c-Src","FES"]) -> pd.DataFrame:
    
        """
        Determine differentially increased reads for
        each gene and each well in a given plate.

        This can only be run after calling preprocess on this HitCaller.
    
        We populate self.hits as a DataFrame with columns:
            well  (str)
            gene  (str)
            p     (float)
            pnorm (float)

        We also populate self.pvalpos as a float holding
        the best p-value from the positive controls.

        Parameters
        ----------
    
        testreps : int
            Parameter forwarded to _testGeneDifferentialCounts.
    
        sampreps : int
            Parameter forwarded to _testGeneDifferentialCounts.

        geneNormPercentile : float
            Parameter forwarded to _pvalNormalization.

        kinases : list
            Parameter forwarded to _pvalNormalization.

        Returns
        -------
        None
        """

        if not self._preprocessed:
            raise ValueError(
                "Must call preprocess before finding differential counts")

        res = list()
        for well, drugdf in self.drugdict.items():
            for gene in drugdf.index:
                if gene == self.negControlGene:
                    continue
                p = self._testGeneDifferentialCounts(
                    dmsodf=self.dmsodf,
                    drugdf=drugdf,
                    gene=gene,
                    testreps=testreps,
                    sampreps=sampreps)
                res.append([well, gene, p])

        self.hits = pd.DataFrame(res, columns=["well","gene","p"])
        self.hits = HitCaller._pvalNormalization(
            self.hits,
            geneNormPercentile=geneNormPercentile,
            kinases=kinases)
        self.pvalpos = self.hits[self.hits.apply(
            lambda row: any([(row.well==k) and (row.gene in v) \
                for k,v in self.posControls.items()]), axis=1)].pnorm.min()
        self.hits["pvalpos"]  = self.pvalpos
        self.hits["magratio"] = \
            np.log(self.hits.pnorm) / np.log(self.pvalpos)

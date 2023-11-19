## Instructions for Preprocessing Example Reads

In this dir, we provide example data and scripts to convert raw reads into read count files for our hit caller.

1. Create or use the `merged.csv` file containing the indexing pair IDs for each plate - 5 series, 7 series, and internal sequences associated with the A-series sample name from a demultiplexed illumina run.

2. Create or use the `models.csv` file containing the barcode IDs and associated internal sequences.

3. Fill in `setup.csv` with associated plate names and primer plate numbers in the following CSV format:

        plateID,primerplate
        NS1519_022057,1
        NS1519_022058,13
        NS1519_022059,3
        NS1519_022060,4
        NS1519_022061,15
        NS1519_022062,6

4. After importing dependencies in the first code cell of `preprocess.ipynb`, run the second code cell to read in the above files as data frames. A readout of the `setup.csv` should be visible. This cell also contains the definitions of `returnrunindicies`, which retrieves a primer plate subset of the `merged.csv`, and `modelwisedindicies`, which creates a BC-wise index for each well of a primer plate subset.

5. Deposit demultiplexed sequencing files in the `gz` directory.

6. Run the third and final code cell to generate the count tables, which are deposited in the `counttables` directory. They can then be moved to a `counts` folder of the data directory.

7. An output with the provided `setup.csv` and `.fastq.gz` files should produce an output identical to the CSV found in `counttables` dir: `NS1519_022057_PP1.csv`, which has is an example plate with details for rows A, B, G, and H.

8.  The associated optical density (595nm) readouts should be placed in the `data/od` directory (outside of this `preproc` dir) for the entire run in the following CSV format:

        plate, well,OD595
        NCIDIV_4865,A01,0.6682
        NCIDIV_4865,B01,0.6074
        NCIDIV_4865,C01,0.6451
        NCIDIV_4865,D01,0.6494
        NCIDIV_4865,E01,0.6378
        NCIDIV_4865,F01,0.6387
        NCIDIV_4865,G01,0.6356
        NCIDIV_4865,H01,0.6034

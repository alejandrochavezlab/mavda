# Identifying broadly active protease inhibitors through multiplex chemical screens

All code necessary to reproduce the results of our manuscript are available as scripts and Jupyter notebooks in this repository.

Instructions and code for preprocessing the raw sequencing reads are found in the `preproc` dir.

Processed barcode read counts are found in the data `counts` directory, partitioned by run date. Corresponding optical density measurements are found in the data `od` dir.

The results of the `callhits.py` script are written to the `hits` dir, which is then analyzed in the `hitanalysis.ipynb` notebook.

## Installation

### Get the source
```
git clone https://github.com/alejandrochavezlab/mavda.git
```
The repository is fewer than 100Mb, so installation generally takes a few seconds depending on internet speed.

### Environment and Dependencies

Execution is OS agnostic and does not require special hardware.

All code was run with [python](https://www.python.org) (3.10.9) using the following dependencies (versions are parenthesized):
* [numpy](https://numpy.org) (1.24.2)
* [pandas](https://pandas.pydata.org) (1.5.3)
* [scipy](https://scipy.org) (1.10.1)
* [matplotlib](https://matplotlib.org) (3.7.2)
* [seaborn](https://seaborn.pydata.org) (0.12.2)
* [biopython](https://biopython.org) (1.81)

## Usage

To preprocess the raw sequencing reads, please see the instructions and code found in the `preproc` dir.

To run the hit caller, execute the following from within the root repository directory: `python callhits.py`

To analyze hits and generate figures, run the `hitanalysis.ipynb` notebook.

## License

See the [LICENSE](LICENSE) file

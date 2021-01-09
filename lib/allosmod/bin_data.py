"""Bin data from a file."""

from __future__ import print_function, absolute_import, division
import optparse


def bin_data(data_file, datacol, datamin, datamax, nbins):
    data = []
    for line in open(data_file):
        spl = line.split()
        data.append(float(spl[datacol]))

    numrange = [0] * (nbins+1)
    lowcut = [0.] * (nbins+1)
    highcut = [0.] * (nbins+1)
    window = (datamax - datamin) / nbins
    zpart = 0

    lowcut[0] = -1*window+datamin
    highcut[0] = datamin
    for i_range in range(1, nbins + 1):
        numrange[i_range] = 0
        lowcut[i_range] = lowcut[i_range-1] + window
        highcut[i_range] = highcut[i_range-1] + window
        lowcut[1] = lowcut[1]-0.01  # make sure to include datamin point
        for d in data:
            if lowcut[i_range] < d and d <= highcut[i_range]:
                numrange[i_range] += 1
                zpart += 1
        lowcut[1] = lowcut[1]+0.01

    # write final output
    for i_range in range(1, nbins + 1):
        if numrange[i_range] > 0:
            midrange = (lowcut[i_range] + highcut[i_range])/2.
            print("%.5f %.5f %.5f" % (midrange, numrange[i_range] / zpart,
                                      numrange[i_range]))


def parse_args():
    usage = """%prog [opts] <data_file> <column> <min> <max> <nbins>

Bin data from a file, i.e. group a single variable array into a probability
distribution.

Data is read from <data_file> and put into <nbins> bins that span the range
from <min> to <max>. For every non-empty bin, the midrange of the bin, and
the fraction and count of the data points that were placed in that bin is
printed. The data file is expected to contain whitespace-separated columns;
only data from <column> (starting at 0) is read.

In both <min> and <max> the 'm' character is treated as '-'.
"""
    parser = optparse.OptionParser(usage)

    opts, args = parser.parse_args()
    if len(args) != 5:
        parser.error("incorrect number of arguments")
    return(args[0], int(args[1]), float(args[2].replace('m', '-')),
           float(args[3].replace('m', '-')), int(args[4]))


def main():
    data_file, datacol, datamin, datamax, nbins = parse_args()
    bin_data(data_file, datacol, datamin, datamax, nbins)


if __name__ == '__main__':
    main()

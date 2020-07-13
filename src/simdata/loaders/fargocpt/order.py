#!/usr/bin/env python3
# Remove doublicate intervals from a column based output file.

import numpy as np
import argparse
import os


def order(x, bound=np.inf, fullind=True):
    # Return an array of indecis, such that values which are monotonically growing are selected from x.
    #print("ordering {} with bound {}".format(x, bound))
    # go from right to left and detect the first occurence, where the values don't grow monotonically
    # only use values that are smaller than the give bound
    x = x[x < bound]
    #print("X with bound applied {}".format(x))
    up = len(x) - 1
    low = 0

    if len(x) == 0:
        return []
    elif len(x) == 1:
        if fullind:
            return [0]
        else:
            return [(0, 1)]
    else:
        for n in range(up, 0, -1):  # this goes from up, ..., 1
            #print("x[{}] = {}, x[{}] = {}".format(n, x[n], n-1, x[n-1]))
            if x[n] <= x[n - 1]:
                low = n
                break

        #print("found interval ({}, {})".format(low, up))

        if low == 0:  # if we traversed the list without any non-monotonic entry, return a tuple with the full array length
            if fullind:
                return [i for i in range(0, up + 1)]
            else:
                return [(0, up + 1)]
        else:  # else repeat for the part of x left to the non-monotonic entry, and return a list of indices tuples
            if fullind:
                return order(x[:low],
                             bound=x[low]) + [i for i in range(low, up + 1)]
            else:
                return order(x[:low], bound=x[low]) + [(low, up + 1)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("datafile", help="Datafile to order")
    parser.add_argument("timecol",
                        type=int,
                        help="Column in which the time is stored")
    parser.add_argument("-o",
                        "--outfile",
                        help="Path to write the output file")
    parser.add_argument('-v',
                        '--verbose',
                        default=False,
                        action='store_const',
                        const=True,
                        help="verbose output")
    args = parser.parse_args()

    if args.verbose:
        print(
            "Removing overlapping time intervals from {}, time is in column {}"
            .format(args.datafile, args.timecol))

    comment_char = '#'
    time = np.genfromtxt(args.datafile,
                         usecols=args.timecol,
                         comments=comment_char)
    inds = order(time)

    if len(inds) == len(time):
        if args.verbose:
            print("File already in order, nothing to be done.")
    else:
        if args.verbose:
            print("Backup file as {}.bak".format(args.datafile))
        if args.outfile is not None:
            srcfile = args.datafile
            outfile = args.outfile
        else:
            srcfile = args.datafile + ".bak"
            os.rename(args.datafile, srcfile)
            outfile = args.datafile
        # Copy selected lines
        with open(srcfile, 'r') as src:
            with open(outfile, 'w') as dst:
                n = 0
                Ncomments = 0
                nind = 0
                for nl, line in enumerate(src):
                    if line.strip()[0] == comment_char:
                        dst.write(line)
                        Ncomments += 1
                        continue
                    elif nl - Ncomments == inds[nind]:
                        dst.write(line)
                        nind += 1
                        if nind == len(inds):
                            break


if __name__ == '__main__':
    main()

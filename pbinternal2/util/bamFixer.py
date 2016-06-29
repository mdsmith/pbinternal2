
import sys
import pysam
import argparse
import re
from sets import Set


def rectify(bam_fn, out_fn=None):
    # open
    bam_file = pysam.AlignmentFile(bam_fn, 'rb')

    # write to a default name
    if not out_fn:
        out_fn = bam_fn[:-4] + ".pbver1.bam"

    # save the header
    header = bam_file.header

    # The RG field is a list of readgroups, its length corresponds to the
    # number of readgroups in the sam file:
    num_rg = len(header['RG'])

    num_row = 0

    # We shouldn't need to touch the readgroup, but we may if there is only one
    # and its ID is different than the one used by the reads
    set_rg = False
    if num_rg == 1:
        rgs = []
        bam_file.reset()
        for row in bam_file:
            num_row += 1
            for tag in row.tags:
                if tag[0] == 'RG':
                    rgs.append(tag[1])
        # get the number of unique readgroup IDs
        if len(Set(rgs)) == 1:
            set_rg = True
            new_rg = list(Set(rgs))[0]
        else:
            sys.exit("Too many read groups among rows!")

    for rgi in range(num_rg):
        ds = header['RG'][rgi]['DS']
        ds = _rectify_string(ds, 'BASECALLERVERSION', 9)
        ds = _rectify_string(ds, 'BINDINGKIT', 'NA')
        ds = _rectify_string(ds, 'SEQUENCINGKIT', 'NA')
        header['RG'][rgi]['DS'] = ds
        if set_rg:
            header['RG'][rgi]['ID'] = new_rg
    out_file = pysam.AlignmentFile(out_fn, "wb", header=header)
    num_out = 0
    bam_file.reset()
    for row in bam_file:
        num_out += 1
        out_file.write(row)
    out_file.close()
    if num_out != num_row:
        sys.exit("Not all rows were transferred!")
    return 0


def _rectify_string(string, key, value):
    match = re.match(key, string)
    if not match:
        newstring = ';'.join([string, key + "=" + str(value)])
    elif string[match.end() + 1] == ";":
        newstring = ''.join([string[:match.end() + 1],
                             str(value), string[match.end() + 1:]])
    else:
        newstring = string
    return newstring


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("bam_fn", help='A bam file to fix')
    parser.add_argument("--out_fn", help='A bam file to fix',
                        default=None)
    args = parser.parse_args()
    sys.exit(rectify(args.bam_fn, args.out_fn))

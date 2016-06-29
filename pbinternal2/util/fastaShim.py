"""A compatibility layer to make modern (bam, xml) resources work with legacy
fasta based reports)"""
import os
import tempfile
from pbcore.io import openAlignmentFile, DataSet

import logging

log = logging.getLogger(__name__)

TEMP_DIR = tempfile.mkdtemp(suffix="pbinternal2")


def ensure_fasta(in_fn, aligned=False):
    """Ensure a fasta file exists for Gepard to consume. Since this was the
    original format for this entire report, we'll create temporary fasta files
    earlier to minimize changes to the original code paths. There has to be an
    easier way to get from bam to fasta, however.

    Args:
        in_fn: A input file name.

    Returns: If fasta, or an xml containing a fasta, that
             fasta file name will be returned. If bam or xml containing a
             bam, a temporary fasta file will be created and its name
             returned.
    """
    log.debug("Temporary fasta files can be found here: {s}".format(
        s=TEMP_DIR))

    if in_fn.endswith('xml'):
        reader = DataSet(in_fn)
        in_fn = reader.toFofn(uri=False)[0]
        # Fall through to following conversion process

    if in_fn.endswith('bam'):
        reader = openAlignmentFile(in_fn)
        fasta_temp = "%s/%s.fasta" % (TEMP_DIR, os.path.basename(in_fn))
        fasta_fh = open(fasta_temp, "w")
        for row in reader:
            name = row.qName
            sequence = row.read(aligned)
            fasta_fh.write('>{n}\n{s}\n'.format(n=name, s=sequence))
        fasta_fh.close()
        return fasta_temp
    return in_fn


def ensure_fastas(in_fns, aligned=False):
    """Ensure a fasta file exists for Gepard to consume. Since this was the
    original format for this entire report, we'll create temporary fasta files
    earlier to minimize changes to the original code paths. There has to be an
    easier way to get from bam to fasta, however.

    Args:
        in_fns: A list of input file names, be they fastas, bams or xmls

    Returns: A list of fastas.
    """
    log.debug("Temporary fasta files can be found here: {s}".format(
        s=TEMP_DIR))

    stack = in_fns[:]
    fns = []
    while stack:
        in_fn = stack.pop()
        if in_fn.endswith('xml'):
            reader = DataSet(in_fn)
            stack.extend(reader.toFofn(uri=False))
            continue

        if in_fn.endswith('bam'):
            reader = openAlignmentFile(in_fn)
            fasta_temp = "%s/%s.fasta" % (TEMP_DIR, os.path.basename(in_fn))
            fasta_fh = open(fasta_temp, "w")
            for row in reader:
                name = row.qName
                sequence = row.read(aligned)
                fasta_fh.write('>{n}\n{s}\n'.format(n=name, s=sequence))
            fasta_fh.close()
            in_fn = fasta_temp
        fns.append(in_fn)
    return fns

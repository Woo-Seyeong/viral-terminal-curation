#!/usr/bin/env python3
"""
Count PKV2 reads supporting each terminal-loop orientation (flip / flop).

Why the assay is built this way
-------------------------------
The two PKV2 termini are inverted repeats of about 2.65 kb, so the 50 bp
flanking the left loop are the exact reverse complement of the 50 bp flanking
the right loop. Sequence close to a loop therefore cannot say which terminus a
read came from, and a read confined to the repeat cannot be assigned at all.
Every read is consequently classified in two steps:

  1. Terminus and strand are fixed by a 25-bp seed taken from unique sequence
     2.75-3.45 kb inside the genome, beyond the inverted repeat. Six seeds per
     terminus are tried so that a basecalling error in one does not lose the
     read.
  2. Within that read, the 50-bp anchor immediately adjacent to the loop is
     located on the loop side of the unique seed, and the 37-bp window beyond
     it is classified as the loop sequence (flip) or its reverse complement
     (flop) by edit distance.

Because step 1 pins the strand, the call in step 2 reflects the isomer the
molecule carried and not the direction in which the read happened to be
sequenced. Reads are additionally reported as spanning when they continue at
least MIN_FLANK bases past the loop into the opposite arm.

Usage
-----
    python3 loop_orientation_count.py reads.fastq.gz [more.fastq ...]
    samtools fastq -F 0x900 mapped_reads.sorted.bam | python3 loop_orientation_count.py -

-F 0x900 drops secondary and supplementary alignments. Output is a short TSV;
nothing large is written.
"""
import sys, gzip, argparse

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def rc(s): return s.translate(COMP)[::-1]

# --- PKV2 curated genome, 309,183 bp ---------------------------------------
LOOP_L = "TAGCTATACAGCCTGTTTAAAGTGGCAATTTAGCTAC"              # bases 1-37
ANCH_L = "GAAAATGAGATATACTAATGACACCACAAAATTCGCTGACCACGAGTAAA"  # bases 38-87
LOOP_R = "TAGCTATACAGCCTGTTTAAAGTGGCAATTTCGCGAC"              # last 37 bases
ANCH_R = "TTTACTCGTGGTCAGCGAATTTTGTGGTGTCATTAGTATATCTCATTTTC"  # 50 bp before it

# unique 25-mers beyond the inverted repeat; offset = distance from that end
UNIQ_L = ["TCACAGTAGTGTGGTCATTACCCGA", "ATGGGTCATAATTAGTTTTGTTCCT",
          "TGTCATTAGATGTGCTAGGGGTAGG", "CTAGGTTCGGTAGGATCTTCTGGTG",
          "AAATTACTGAATTCCGTTCCATCAC", "TCGATAGTTCCCGCACAGCTGAACG"]
UNIQ_R = ["GTTTATATGGTGTTTTGTATGGCGT", "GTCAAAGGGGGCACCTGCTTCAATG",
          "GCTTCCCAAATTATCTGTGATGTCA", "AGCTGCGTGATGCAGTGGCATGGTT",
          "CACACACACATCGTGACCTTTATGA", "CCCATGAAGCTGCATCGTCGTTCGC"]

ANCH_SEED = 20      # length of the seeds tiled across the loop-proximal anchor
ANCH_STEP = 6
PAD = 12            # slack either side of the loop window, for indels
MAX_EDIT = 9        # maximum edit distance accepted for a loop call
MIN_MARGIN = 4      # required gap between the best and second-best orientation
MIN_FLANK = 200     # bases past the loop needed to count a read as spanning


def edit(pattern, text):
    """Minimum edit distance between `pattern` and any substring of `text`.

    The loop window is extracted with slack either side so that indels cannot
    truncate it, so the comparison must be free to skip unmatched bases at both
    ends of the window; a global distance would be dominated by that slack."""
    prev = [0] * (len(text) + 1)
    for ca in pattern:
        cur = [prev[0] + 1]
        for j, cb in enumerate(text, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return min(prev)


def tile(seq, k, step):
    return [seq[i:i + k] for i in range(0, len(seq) - k + 1, step)]


class Terminus:
    """One terminus. 'outward_left' is True when, with the unique seed read in
    its genome orientation, the loop lies to the left of the anchor."""

    def __init__(self, name, loop, anchor, uniq, outward_left):
        self.name, self.loop, self.anchor = name, loop, anchor
        self.outward_left = outward_left
        self.uniq = [(u, rc(u)) for u in uniq]
        self.aseeds = [(s, i) for i, s in
                       zip(range(0, len(anchor), ANCH_STEP),
                           tile(anchor, ANCH_SEED, ANCH_STEP))]
        self.n = dict(flip=0, flop=0, unclassified=0,
                      no_anchor=0, assigned=0)
        self.span = dict(flip=0, flop=0)

    def _find_anchor(self, seq, upos, leftward):
        """Anchor occurrence closest to the loop, on the loop side of upos."""
        best = None
        for s, off in self.aseeds:
            start = 0
            while True:
                p = seq.find(s, start)
                if p < 0:
                    break
                a0 = p - off                      # implied anchor start
                if leftward and a0 < upos:
                    if best is None or a0 > best:
                        best = a0
                elif (not leftward) and a0 > upos:
                    if best is None or a0 < best:
                        best = a0
                start = p + 1
        return best

    def scan(self, seq):
        for fwd, rev in self.uniq:
            flipped = False
            upos = seq.find(fwd)
            if upos < 0:
                upos = seq.find(rev)
                if upos < 0:
                    continue
                flipped = True
            break
        else:
            return False

        self.n["assigned"] += 1
        loop = self.loop if not flipped else rc(self.loop)
        anchor_len = len(self.anchor)
        leftward = self.outward_left != flipped   # direction of the terminus
        a0 = self._find_anchor(seq, upos, leftward)
        if a0 is None:
            self.n["no_anchor"] += 1
            return True

        if leftward:
            lo, hi = a0 - len(self.loop) - PAD, a0 + PAD
            flank = lo
        else:
            lo, hi = a0 + anchor_len - PAD, a0 + anchor_len + len(self.loop) + PAD
            flank = len(seq) - hi
        if lo < 0 or hi > len(seq):
            self.n["no_anchor"] += 1
            return True

        win = seq[lo:hi]
        f, r = edit(loop, win), edit(rc(loop), win)
        if f <= MAX_EDIT and f + MIN_MARGIN <= r:
            call = "flip"
        elif r <= MAX_EDIT and r + MIN_MARGIN <= f:
            call = "flop"
        else:
            self.n["unclassified"] += 1
            return True
        self.n[call] += 1
        if flank >= MIN_FLANK:
            self.span[call] += 1
        return True


def reads(paths):
    for path in paths:
        fh = sys.stdin if path == "-" else (
            gzip.open(path, "rt") if path.endswith(".gz") else open(path))
        with fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    yield line.strip().upper()


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fastq", nargs="+", help="FASTQ or FASTQ.gz files, or - for stdin")
    ap.add_argument("--min-length", type=int, default=5000,
                    help="skip reads shorter than this (default 5000)")
    a = ap.parse_args()

    termini = [Terminus("left", LOOP_L, ANCH_L, UNIQ_L, True),
               Terminus("right", LOOP_R, ANCH_R, UNIQ_R, False)]
    total = kept = 0
    for seq in reads(a.fastq):
        total += 1
        if len(seq) < a.min_length:
            continue
        kept += 1
        for t in termini:
            if t.scan(seq):
                break

    w = sys.stdout.write
    w(f"# reads read\t{total}\n")
    w(f"# reads >= {a.min_length} bp\t{kept}\n")
    w("terminus\tassigned\tflip\tflop\tunclassified\tno_anchor"
      "\tflip_spanning\tflop_spanning\n")
    for t in termini:
        n, s = t.n, t.span
        w(f"{t.name}\t{n['assigned']}\t{n['flip']}\t{n['flop']}"
          f"\t{n['unclassified']}\t{n['no_anchor']}"
          f"\t{s['flip']}\t{s['flop']}\n")
    w(f"# terminus assigned by a 25-bp exact seed from unique sequence "
      f"2.75-3.45 kb inside the genome (6 seeds per terminus)\n"
      f"# loop window located by {ANCH_SEED}-bp seeds tiled across the 50-bp "
      f"loop-proximal anchor; window +/-{PAD} bp\n"
      f"# loop called at edit distance <= {MAX_EDIT} with a margin of "
      f">= {MIN_MARGIN} over the alternative orientation\n"
      f"# spanning = read continues >= {MIN_FLANK} bp past the loop\n")


if __name__ == "__main__":
    main()

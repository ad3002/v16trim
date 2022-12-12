#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 01.01.2021
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com

import os, sys
import aindex
import argparse
from trseeker.tools.sequence_tools import get_revcomp
from collections import Counter
from trseeker.seqio.fasta_file import sc_iter_fasta_brute


lefts = "CCTACGGGAGGCAGCAG,CCTACGGGAGGCTGCAG,CCTACGGGCGGCAGCAG,CCTACGGGCGGCTGCAG,CCTACGGGTGGCAGCAG,CCTACGGGTGGCTGCAG,CCTACGGGGGGCAGCAG,CCTACGGGGGGCTGCAG"
lefts = lefts.split(",")
rights = "GACTACACGGGTATCTAATCC,GACTACAGGGGTATCTAATCC,GACTACATGGGTATCTAATCC,GACTACGCGGGTATCTAATCC,GACTACGGGGGTATCTAATCC,GACTACGTGGGTATCTAATCC,GACTACTCGGGTATCTAATCC,GACTACTGGGGTATCTAATCC,GACTACTTGGGTATCTAATCC".split(",")
for i in range(len(rights)):
    rights[i] = get_revcomp(rights[i])
k = 23



TO_IUPAC_CODE = {
    ('A','G'): 'R',
    ('G','A'): 'R',

    ('C','T'): 'Y',
    ('T','C'): 'Y',

    ('C','G'): 'S',
    ('G','C'): 'S',

    ('A','T'): 'W',
    ('T','A'): 'W',

    ('G','T'): 'K',
    ('T','G'): 'K',

    ('A','C'): 'M',
    ('C','A'): 'M',
}


def hamming_distance(s1, s2):
    return sum(i != j for i, j in zip(s1,s2) if i != 'N' and j != 'N')

def get_poses(seq, primer):
    '''
    '''
    found = []
    hd_ = 100000
    for i in range(0, len(seq)-len(primer)+1):
        substring = seq[i:i+len(primer)]
        hd = hamming_distance(primer, substring)
        if hd < hd_:
            found = [i]
            hd_ = hd
        elif hd == hd_:
            found.append(i)
    return tuple(found), hd_


def find_primer(seq, pattern):
    '''
    '''
    found = -1
    hd_ = 100000
    k = len(pattern)
    for i in range(0, len(seq)-k+1):
        substring = seq[i:i+k]
        hd = hamming_distance(pattern, substring)
        if hd_ > hd:
            found = i
            hd_ = hd
            if hd_ == 0:
                break
    return found, hd_


def get_shift(left, right):
    '''
    '''
    pos = -1
    hd_ = 100000
    n = len(left)
    kmer = right[:20]
    start = n - ExpParams["left_right_maximal_intersection"]
    for i in range(start, n-20+1):
        left_substring = left[i:i+20]
        right_substring = kmer
        hd = hamming_distance(left_substring, right_substring)
        if hd < hd_:
            pos = i
            hd_ = hd
        if hd_ == 0:
            break

    left_part = left[:pos]
    left_middle = left[pos:]
    right_middle = right[:len(left_middle)]
    right_part = right[len(left_middle):]

    if not left_part:
        return "SHORT", 100000

    try:
        assert len(left_middle) == len(right_middle)
    except:
        print(left, right)
        print(pos, hd_, left_part, "|", left_middle,"|", right_middle,"|", right_part)
        assert len(left_middle) == len(right_middle)

    hd_ = hamming_distance(left_middle, right_middle)

    return pos, hd_


def filter_reads(settings):
    '''
    '''
    sample = settings["sample"]
    file_name = settings["reads_file"]

    left_pattern = settings["left_pattern"]
    right_pattern = settings["right_pattern"]

    total = 0
    left_positive = 0
    too_short = 0
    ok = 0

    mergable = []
    scrap = []
    noleft = []
    noleft_n = 0
    
    left_k = len(left_pattern)

    with open(file_name) as fh:
        for i, line in enumerate(fh):
            if i % 10000 == 0:
                print(sample, "F", i)
            left_found = False
            read = line.strip()
            left, right = read.split("~")
            total += 1
            left_pos = 0

            if left[:left_k] in lefts:
                left_positive += 1
                left_found = True
                left_pos = 0
                left_hd = 0
            else:
                left_pos, left_hd = find_primer(left, left_pattern)
                if left_pos <= ExpParams["left_primer_shift_R"] and left_hd < 5:
                    left_positive += 1
                    left_found = True
                else:
                    noleft.append((left, right, left_pos, left_hd))
                    noleft_n += 1
                    continue
            left = left[left_pos:]
            left_pos = 0

            shift_pos, shift_hd = get_shift(left, right)
            if shift_pos == "SHORT":
                too_short += 1
                scrap.append((left, right, left_pos, shift_pos, shift_hd, "SHORT"))
                continue
            if shift_hd < 10:
                ok += 1
                mergable.append((left, right, left_pos, shift_pos, shift_hd, "OK"))
            else:
                scrap.append((left, right, left_pos, shift_pos, shift_hd, "HD"))
    
    sample2statistics = {
        "total": total,
        "too_short": too_short,
        "left_positive": left_positive,
        "noleft_n": noleft_n,
        "noleft_p": round(100.*noleft_n/total, 3),
        "merged_n": ok,
        "merged_p": round(100.*ok/total, 3),
        "left_positive_p": round(100.*left_positive/total, 3),
        "too_short_p": round(100.*too_short/total, 3), 
        "sample": sample,
    }
    print()
    print(sample, sample2statistics)
    sample2statistics["mergeable"] = mergable
    sample2statistics["scrap"] = scrap
    sample2statistics["noleft"] = noleft
    return sample2statistics


def merge_reads(sample2statistics):
    '''
    '''
    merged = []
    sample = sample2statistics["sample"]
    length_distribution = Counter()
    for left, right, lpos, pos, hd, status in sample2statistics["mergeable"]:
        if len(merged) % 10000 == 0:
            print(sample, "M", len(merged))
        seq = []
        alt = []
        
        left_part = left[:pos]
        left_middle = left[pos:]
        right_middle = right[:len(left_middle)]
        right_part = right[len(left_middle):]

        middle = []
        for i in range(len(left_middle)):
            if left_middle[i] == right_middle[i]:
                middle.append(left_middle[i])
            elif left_middle[i] == 'N' and right_middle[i] != 'N':
                middle.append(right_middle[i])
            elif left_middle[i] != 'N' and right_middle[i] == 'N':
                middle.append(left_middle[i])
            else:
                middle.append('N')

        middle = "".join(middle)

        seq = [left_part, middle, right_part]

        s = "".join(seq)
        length_distribution[len(s)] += 1
        merged.append((s, left_middle, right_middle, pos, len(s), hd))
    sample2statistics["merged"] = merged
    sample2statistics["length_distribution"] = length_distribution
    print()
    return sample2statistics


def get_correct(sample2statistics):
    '''
    '''
    k = 23
    correct = []
    sample = sample2statistics["sample"]
    print(sample, sample2statistics["merged_n"])
    for i,(read,_) in enumerate(sample2statistics["merged"]):
        if i % 10000 == 0:
            print(sample, "C", i)
        cov = min([all2tf[read[j:j+k]] for j in range(len(read)-k+1)])
        if cov > 50:
            correct.append(i)
    sample2statistics["corrected"] = correct
    sample2statistics["correct"] = len(correct)
    sample2statistics["correct_p"] = round(100.*len(correct)/sample2statistics["total"], 2)
    print()
    return sample2statistics


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Filter and merge reads.')
    parser.add_argument('-i', help='File with reads or fastq1,fastq2 pair', required=True)
    parser.add_argument('-t', help='Type fastq or reads', required=False, default="fastq")
    parser.add_argument('-o', help='Output prefix', required=True)
    parser.add_argument('-a', help='Aindex prefix', required=False)
    parser.add_argument('-l', help='Left pattern', required=False, default="CCTACGGGAGGCAGCAG")
    parser.add_argument('-r', help='Right pattern', required=False, default="GGATTAGATACCCGTGTAGTC")
    parser.add_argument('--correction', help='Run correctionAindex prefix', required=False, default=False)

    ExpParams = {
        "left_right_shift_L": 0,
        "left_right_shift_R": 5,
        "left_primer_shift_L": 0,
        "left_primer_shift_R": 10,
        "left_right_maximal_intersection": 50,
    }
    
    args = vars(parser.parse_args())

    settings = {}
    settings["reads_file"] = args["i"]
    settings["sample"] = args["o"]
    aindex_prefix = args["a"]
    prefix = args["o"]
    left_pattern = args["l"]
    right_pattern = args["r"]
    correction = bool(args["correction"])

    settings["left_pattern"] = left_pattern
    settings["right_pattern"] = right_pattern

    if left_pattern == "auto" or right_pattern == "auto":
        if left_pattern == "auto":
            print("Autodetect left primer...")
        else:
            print("Autodetect right primer...") 
        maybe = Counter()
        global_counter = Counter()
        with open(settings["reads_file"]) as fh:
            for i, line in enumerate(fh):
                if i % 10000 == 0:
                    print(settings["sample"], "LAUTODETECT", i)
                if left_pattern == "auto":
                    maybe[line[:13]] += 1
                else:
                    line = line.strip()
                    maybe[line[-13:]] += 1

        # silva_k13_pos_file = "/home/akomissarov/Dropbox/workspace/silva.4sector.k13.pos.v2"
        silva_align = "/mnt/data/prokaryota/16Smock/SILVA_138.1_SSURef_tax_silva_full_align_trunc.fasta"
        for primer, tf in maybe.most_common(20):
            starts = Counter()
            for i,(header, seq) in enumerate(sc_iter_fasta_brute(silva_align)):
                if not "Bacteria" in header:
                    continue
                seq = seq.replace("-","").replace(".","").replace("U", "T")
                s,hd = get_poses(seq, primer)
                starts[(s,hd)] += 1
                global_counter[(s,hd,primer)] += 1
                if i > 1000:
                    break
            print(starts.most_common(10))

        print(global_counter.most_common(10))        
        print("Done.")
        sys.exit(0)

    if correction:
        aindex_settings = {
            "index_prefix": aindex_prefix,
            "aindex_prefix": None,
            "reads_file": None,
            "header_file": None,
            "max_tf": 10000,
        }
        all2tf = aindex.load_aindex(aindex_settings, skip_aindex=False, skip_reads=False)

        sample2statistics = get_correct(merge_reads(filter_reads(settings)))
    else:
        sample2statistics = merge_reads(filter_reads(settings))

    output_file_stats = prefix + ".stats"
    with open(output_file_stats, "w") as fw_stats:
        for key in sample2statistics:
            if key == "scrap":
                output_file = prefix + "_scrap.reads"
                with open(output_file, "w") as fw:
                    for left, right, lpos, pos, hd, status in sample2statistics[key]:
                        fw.write("%s\t%s\t%s\t%s\t%s\n" % (left, right, lpos, pos, hd))
                continue
            if key == "noleft":
                output_file = prefix + "_noleft.reads"
                with open(output_file, "w") as fw:
                    for left, right, pos, hd in sample2statistics[key]:
                        fw.write("%s\t%s\t%s\t%s\n" % (left, right, pos, hd))
                continue
            if key == "merged":
                output_file = prefix + "_merged.reads"
                with open(output_file, "w") as fw:
                    for s, left_middle, right_middle, pos, l, hd in sample2statistics[key]:
                        fw.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (s, left_middle, right_middle, pos, l, hd ))
                continue
            if key == "corrected":
                output_file = prefix + "_correct.rids"
                with open(output_file, "w") as fw:
                    for rid in sample2statistics[key]:
                        fw.write("%s\n" % (rid))
                continue
            if key in "mergeable":
                output_file = prefix + "_mergeable.rids"
                with open(output_file, "w") as fw:
                    for left, right, lpos, pos, hd, status in sample2statistics[key]:
                        fw.write("%s\t%s\t%s\t%s\t%s\n" % (left, right, lpos, pos, hd))
                continue
            print("%s\t%s" % (key, sample2statistics[key]))
            fw_stats.write("%s\t%s\n" % (key, sample2statistics[key]))

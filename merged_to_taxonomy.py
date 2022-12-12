#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 01.01.2021
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com

from collections import *
import argparse
import time
import os, sys
from trseeker.seqio.fasta_file import sc_iter_fasta_brute

AMBIGIOUS_DNA_DICT = {
    "A": "A",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "X": "GATC",
    "N": "GATC",


    "C": "C",
    "G": "G",
    "T": "T",
    
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "B": "CGT",
    
}


DISTANCE_MATRIX = {
    'A': {'A': 0, 'C': 1, 'G': 1, 'T': 1, 'R': 0, 'Y': 1, 'S': 1, 'W': 0, 'K': 1, 'M': 0, 'B': 1, 'D': 0, 'H': 0, 'V': 0, 'N': 0},
    'C': {'A': 1, 'C': 0, 'G': 1, 'T': 1, 'R': 1, 'Y': 0, 'S': 0, 'W': 1, 'K': 1, 'M': 0, 'B': 0, 'D': 1, 'H': 0, 'V': 0, 'N': 0},
    'G': {'A': 1, 'C': 1, 'G': 0, 'T': 1, 'R': 0, 'Y': 1, 'S': 0, 'W': 1, 'K': 0, 'M': 1, 'B': 0, 'D': 0, 'H': 1, 'V': 0, 'N': 0},
    'T': {'A': 1, 'C': 1, 'G': 1, 'T': 0, 'R': 1, 'Y': 0, 'S': 1, 'W': 0, 'K': 0, 'M': 1, 'B': 0, 'D': 0, 'H': 0, 'V': 1, 'N': 0},
    'R': {'A': 0, 'C': 1, 'G': 0, 'T': 1, 'R': 0, 'Y': 1, 'S': 0, 'W': 0, 'K': 0, 'M': 0, 'B': 0, 'D': 0, 'H': 0, 'V': 0, 'N': 0},
    'Y': {'A': 1, 'C': 0, 'G': 1, 'T': 0, 'R': 1, 'Y': 0, 'S': 0, 'W': 0, 'K': 0, 'M': 0, 'B': 0, 'D': 0, 'H': 0, 'V': 0, 'N': 0},
    'S': {'A': 1, 'C': 0, 'G': 0, 'T': 1, 'R': 0, 'Y': 0, 'S': 0, 'W': 1, 'K': 0, 'M': 0, 'B': 0, 'D': 0, 'H': 0, 'V': 0, 'N': 0},
    'W': {'A': 0, 'C': 1, 'G': 1, 'T': 0, 'R': 0, 'Y': 0, 'S': 1, 'W': 0, 'K': 0, 'M': 0, 'B': 0, 'D': 0, 'H': 0, 'V': 0, 'N': 0},
    'K': {'A': 1, 'C': 1, 'G': 0, 'T': 0, 'R': 0, 'Y': 0, 'S': 0, 'W': 0, 'K': 0, 'M': 1, 'B': 0, 'D': 0, 'H': 0, 'V': 0, 'N': 0},
    'M': {'A': 0, 'C': 0, 'G': 1, 'T': 1, 'R': 0, 'Y': 0, 'S': 0, 'W': 0, 'K': 1, 'M': 0, 'B': 0, 'D': 0, 'H': 0, 'V': 0, 'N': 0},
    'M': {'A': 0, 'C': 0, 'G': 1, 'T': 1, 'R': 0, 'Y': 0, 'S': 0, 'W': 0, 'K': 1, 'M': 0, 'B': 0, 'D': 0, 'H': 0, 'V': 0, 'N': 0},
    'B': {'A': 1, 'C': 0, 'G': 0, 'T': 0, 'R': 0, 'Y': 0, 'S': 0, 'W': 0, 'K': 0, 'M': 0, 'B': 0, 'D': 0, 'H': 0, 'V': 0, 'N': 0},
    'D': {'A': 0, 'C': 1, 'G': 0, 'T': 0, 'R': 0, 'Y': 0, 'S': 0, 'W': 0, 'K': 0, 'M': 0, 'B': 0, 'D': 0, 'H': 0, 'V': 0, 'N': 0},
    'H': {'A': 0, 'C': 0, 'G': 1, 'T': 0, 'R': 0, 'Y': 0, 'S': 0, 'W': 0, 'K': 0, 'M': 0, 'B': 0, 'D': 0, 'H': 0, 'V': 0, 'N': 0},
    'V': {'A': 0, 'C': 0, 'G': 0, 'T': 1, 'R': 0, 'Y': 0, 'S': 0, 'W': 0, 'K': 0, 'M': 0, 'B': 0, 'D': 0, 'H': 0, 'V': 0, 'N': 0},
    'N': {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'R': 0, 'Y': 0, 'S': 0, 'W': 0, 'K': 0, 'M': 0, 'B': 0, 'D': 0, 'H': 0, 'V': 0, 'N': 0}
}

def get_reverse_complement(seq):
    '''
    IUPAC reverse complement function
    this function makes reverse complement sequence
    of the IUPAC written sequence
    params: string - which string to edit
    output: string - reverse complement sequence
    '''
    seq_cur = str(seq)
    seq_comp = seq_cur.maketrans("ATCGRYMKSWHBVDN", "TAGCYRKMSWDVBHN")
    return seq_cur.translate(seq_comp)[::-1]


def hamming_distance(exact_string, amb_string):
    """ Get Hamming distance: the number of corresponding symbols that differs in given strings.
    """
    return sum(DISTANCE_MATRIX[i][j] for i, j in zip(amb_string, exact_string))


def hamming_distance_raw(exact_string, amb_string):
    """ Get Hamming distance: the number of corresponding symbols that differs in given strings.
    """
    return sum(i != j for i, j in zip(amb_string, exact_string))



def find_primer(seq, primer, start=None, end=None):
    '''
    '''
    k = len(primer)
    if not start:
        start = 0
    if not end:
        end = len(seq) - k + 1
    pos = []
    hd_ = 100000
    for i in range(start, end):
        substring = seq[i:i+k]
        hd = hamming_distance(primer, substring)
        if hd < hd_:
            pos = [i]
            hd_ = hd
            if hd_ == 0:
                break
        elif hd == hd_:
            pos.append(i)
    return pos, hd_

def save(start_time, prefix, result, statistics):
    print("Saving taxonomy", f"{(time.time()-start_time)/60} min")
    with open(prefix + "_taxonomy.tsv", "w") as fw:
        for d in result:
            fw.write("%s\n" % "\t".join(map(str, d)))

    with open(prefix + "_taxonomy.stats", "w") as fw:    
        fw.write("%s\t%s\n" % ("ok", statistics["ok"]))
        fw.write("%s\t%s\n" % ("failed", statistics["failed"]))
        for (hd,size) in statistics["hds"].most_common():
            fw.write("%s\t%s\n" % ("hd=%s" % hd, size))

    with open(prefix + "_taxonomy.taxons", "w") as fw:
        for (taxon,size) in statistics["taxons"].most_common():
            fw.write("%s\t%s\n" % (taxon, size))

def match_with_kmers(rid_i, read, silva, kmer2sid_good, statistics):    
    found_minhd = 100000
    good_kmer = None
    k = 13
    kmer2check = []
    for i in range(len(read)-k+1):
        kmer = f"{read[i:i+k]}:{i}"
        if kmer in kmer2sid_good:
            kmer2check.append((len(kmer2sid_good[kmer]), kmer))

    kmer2check.sort(key=lambda x: x[0])
    hits = []
    names = [] 
    checked_silva = {}
    good_kmer = None
    current_scope = set()
    hdhits = Counter()
    for tf, kmer in kmer2check:
        # print(hits)
        # print(f"hd {found_minhd} tf {tf} {kmer} scope: {len(current_scope)}")
        # input("?")
        if found_minhd == 0:
            break
        hdhits[found_minhd] += 1
        if hdhits.most_common(1)[0][1] > 50:
            break
        for sid in kmer2sid_good[kmer]:
            if not sid in checked_silva:
                hd = hamming_distance(silva[sid][0][:len(read)], read)
                checked_silva[sid] = hd
            hd = checked_silva[sid]
            if hd < found_minhd:
                found_minhd =  hd
                hits = [silva[sid][2]]
                names = [silva[sid][1]]
                good_kmer = kmer
                current_scope.add(sid)
            elif hd == found_minhd:
                hits.append(silva[sid][2])
                names.append(silva[sid][1])
                current_scope.add(sid)

    hit = [Counter([y[x] for y in hits if len(y) > x ]) for x in range(7)]  
    hit = [x.popitem()[0] for x in hit if len(x) == 1]
    result.append((sample, rid_i, good_kmer, found_minhd, hit, names))
    statistics["ok"] += 1
    statistics["hds"][found_minhd] += 1
    statistics["taxons"][tuple(hit)] += 1
    for name in names:
        statistics["names"][name] += 1
    if not hit:
        statistics["failed"] += 1

    # print(rid_i, read)
    # print(f"hd: {found_minhd} hit: {hit} names: {names} kmer: {good_kmer}")
    # input("?")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Compute taxonomy.')
    parser.add_argument('-i', help='File with merged reads', required=True)
    parser.add_argument('-o', help='Output prefix', required=True)
    parser.add_argument('-l', help='Left pattern : left pos', required=True)
    parser.add_argument('-r', help='right pattern : right_pos', required=True)
    args = vars(parser.parse_args())

    merged_reads_file = args["i"]
    sample = args["o"]
    left_pattern, left_pos = args["l"].split(":")
    right_pattern, right_pos = args["r"].split(":")
    left_pos, right_pos = int(left_pos), int(right_pos)
    prefix = sample

    right_pattern = get_reverse_complement(right_pattern)

    silva_file = "/home/akomissarov/Dropbox/workspace/silva.4sector.data"
    silva_k13_file = "/home/akomissarov/Dropbox/workspace/silva.4sector.k13"
    silva_align = "/mnt/data/prokaryota/16Smock/SILVA_138.1_SSURef_tax_silva_full_align_trunc.fixed.fasta"

    silva_db_file_name = "%s.%s.db" % (left_pattern, right_pattern)
    silva_db_uncall_file_name = "%s.%s.db_uncall" % (left_pattern, right_pattern)
    silva_db_k13_file_name = "%s.%s.db_k13" % (left_pattern, right_pattern)

    start_time = time.time()
    k = 13
    ### check DB for left patter if not build it
    if not os.path.isfile(silva_db_file_name):
        ### build it
        silva = []
        uncalled = []
        print(left_pattern, right_pattern)
        for i,(header, seq) in enumerate(sc_iter_fasta_brute(silva_align)):
            poses_right = None
            hd_right = None
            if i % 25000 == 0:
                print(sample, "BUILDDB", f"total: {i}, uncalled: {len(uncalled)} {round((time.time()-start_time)/60,2)} min")
            if not "Bacteria" in header:
                continue
            poses, hd = find_primer(seq, left_pattern, start=left_pos-100, end=left_pos+100)
            if hd > 3:
                poses, hd = find_primer(seq, left_pattern)
            if hd < 5 and len(poses) == 1:
                pos = poses[0]
                seq = seq[pos:]
                poses_right, hd_right = find_primer(seq, right_pattern)
                if hd_right < 5 and len(poses_right) == 1:
                    right_pos = poses_right[0]
                    seq = seq[:right_pos+len(right_pattern)]
                    silva.append((header, seq, poses, hd, poses_right, hd_right))
                else:
                    uncalled.append((header, seq, poses, hd, poses_right, hd_right))
            else:
                uncalled.append((header, seq, poses, hd, poses_right, hd_right))

        print(f'silva length {len(silva)}, uncalled length: {len(uncalled)}')

        print(sample, "SAVEDB_SILAV")
        with open(silva_db_file_name, "w") as fw:
            for header, seq, poses, hd, poses_right, hd_right in silva:
                poses = ','.join(map(str, poses))
                poses_right = ','.join(map(str, poses_right))
                fw.write(f"{header}\t{seq}\t{poses}\t{hd}\t{poses_right}\t{hd_right}\n")
        print(sample, "SAVEDB_UNCALL")
        with open(silva_db_uncall_file_name, "w") as fw:
            for header, seq, poses, hd, poses_right, hd_right in uncalled:
                poses = ','.join(map(str, poses))
                poses_right = ','.join(map(str, poses_right))
                fw.write(f"{header}\t{seq}\t{poses}\t{hd}\t{poses_right}\t{hd_right}\n")
    
    print(sample, "LOAD SILVA DB:")
    silva = []
    with open(silva_db_file_name) as fh:
        for line in fh:
            header, seq, *_ = line.strip().split('\t')
            silva.append((header, seq))

    print(f"Load SILVA {sample} and {round((time.time()-start_time)/60,2)} min")
    silva_seq2silvaid = {}
    len2silvaids = defaultdict(list)
    silvaid2taxons = defaultdict(list)
    silvaid2seq = {}
    silvaid2name = {}
    length_dist = Counter()
    for i, (header, seq) in enumerate(silva):
        if i % 250000 == 0:
            print(sample, "SILVA", i, f"{round((time.time()-start_time)/60,2)} min")
        sid, *other = header.split()

        taxons = " ".join(other).split(";")
        assert len(taxons) <= 7
        for j in range(7 - len(taxons)):
            taxons.append("")
        assert len(taxons) == 7

        if seq in silva_seq2silvaid:
            silvaid = silva_seq2silvaid[seq]
            silvaid2taxons[silvaid].append(taxons)
            continue
        silva_seq2silvaid[seq] = i
        silvaid2seq[i] = seq
        silvaid2name[i] = sid
        silvaid2taxons[i].append(taxons)
        len2silvaids[len(seq)].append(i)
        length_dist[len(seq)] += 1

    for sid in silvaid2taxons:
        if len(silvaid2taxons[sid]) == 1:
            silvaid2taxons[sid] = silvaid2taxons[sid][0]
        else:
            taxons = []
            for i in range(6):
                taxon = list(set([x[i] for x in silvaid2taxons[sid]]))
                if len(taxon) == 1:
                    taxons.append(taxon[0])
                else:
                    break
            else:
                taxon = list(set([x[6] for x in silvaid2taxons[sid]]))
                if len(taxon) == 1:
                    taxons.append(taxon[0])
                else:
                    taxon = list(set([x[6].split()[0] for x in silvaid2taxons[sid]]))
                    taxon = [x for x in taxon if not x.lower() in ['uncultured', 'unidentified']]
                    if len(taxon) == 1:
                        taxons.append(taxon[0])
                        taxon = list(set([" ".join(x[6].split()[1:]) for x in silvaid2taxons[sid]]))
                        if len(taxon) == 1:
                            taxons.append(taxon[0])
                        else:
                            taxons.append("")
                    else:
                        taxons.append("")
            for j in range(7 - len(taxon)):
                taxon.append("")
            silvaid2taxons[sid] = taxon

    print(f"...loaded {len(silva_seq2silvaid)} uniq seqs")

    if not os.path.isfile(silva_db_k13_file_name):
        k = 13
        print(sample, "BUILD KMER DB:")
        kmer2sid_good = defaultdict(list)
        for i, seq in enumerate(silva_seq2silvaid):
            sid = silva_seq2silvaid[seq]
            L = len(seq)
            if i % 100000 == 0:
                print(sample, "BUILDDB_COUNT", i, len(silva_seq2silvaid), f"{round((time.time()-start_time)/60,2)} min")
            for i in range(len(seq)-k+1):
                kmer2sid_good[f"{seq[i:i+k]}:{L}"].append(sid)

        print(sample, "SAVEDB_COUNT:")
        with open(silva_db_k13_file_name, "w") as fw:
            N = len(kmer2sid_good)
            for i, kmer in enumerate(kmer2sid_good):
                if i % 100000 == 0:
                    print(round(100.*i/N, 3), "%")
                if len(kmer2sid_good[kmer]) < 5:
                    continue
                fw.write("%s\t%s\n" % (kmer, ",".join(map(str,kmer2sid_good[kmer]))))
        print()
    else:
        print("Load kmer2sid_good")
        kmer2sid_good = {}
        with open(silva_db_k13_file_name) as fh:
            for i, line in enumerate(fh):
                if i % 100000 == 0:
                    print(sample, "LOADKMER", i, f"{round((time.time()-start_time)/60,2)} min")
                kmer, hits = line.strip().split("\t")
                hits = list(map(int, hits.split(",")))
                kmer2sid_good[kmer] = hits


    print(f"Load reads {sample} and {round((time.time()-start_time)/60,2)} min")

    merged_reads = defaultdict(int)
    
    with open(merged_reads_file) as fh:
        for line in fh:
            read = line.split()[0]
            merged_reads[read] += 1
    
    reads_length_dist = Counter() 
    for read in merged_reads:
        reads_length_dist[len(read)] += 1
    print(reads_length_dist.most_common(100))

    print(f"...load {i} reads and {len(merged_reads)} uniq reads in {round((time.time()-start_time)/60,2)} min")
    result = []

    statistics = {
        "ok": 0,
        "failed": 0,
        "hds": Counter(),
        "taxons": Counter(),
    }

    print(f"Process dada {sample} of size {len(merged_reads)} and {round((time.time()-start_time)/60,2)} min")

    '''

    L -> (L1 choosing) -> (L2 choosing)

    '''
    for rid_i, read in enumerate(merged_reads):
        if rid_i % 10000 == 0:
            print(sample, rid_i, len(merged_reads), f"{round((time.time()-start_time)/60,2)} min")
            print(statistics["ok"], statistics["failed"])
            print(statistics["hds"].most_common(10))
            print(statistics["taxons"].most_common(10))
            save(start_time, prefix, result, statistics)
        L = len(read)
        found_minhd = 100000
        good_kmer = None
        k = 13
        kmer2check = []
        for i in range(len(read)-k+1):
            kmer = f"{read[i:i+k]}:{L}"
            if kmer in kmer2sid_good:
                kmer2check.append((len(kmer2sid_good[kmer]), kmer))

        kmer2check.sort(key=lambda x: x[0])
        hits = []
        names = [] 
        checked_silva = {}
        good_kmer = None
        current_scope = set()
        hdhits = Counter()
        sids = []
        for tf, kmer in kmer2check:
            # print(hits)
            # print(f"hd {found_minhd} tf {tf} {kmer} scope: {len(current_scope)}")
            # input("?")
            if found_minhd == 0:
                break
            hdhits[found_minhd] += 1
            if hdhits.most_common(1)[0][1] > 50:
                break
            for sid in kmer2sid_good[kmer]:
                if not sid in checked_silva:
                    hd = hamming_distance(silvaid2seq[sid], read)
                    checked_silva[sid] = hd
                hd = checked_silva[sid]
                if hd < found_minhd:
                    found_minhd =  hd
                    hits = [silvaid2taxons[sid]]
                    sids = [sid]
                    good_kmer = kmer
                    current_scope.add(sid)
                elif hd == found_minhd:
                    hits.append(silvaid2taxons[sid])
                    current_scope.add(sid)
                    sids.append(sid)

        '''
        every hit tuple of 7 elements
        '''
        hits = list(set([tuple(x) for x in hits]))
        if len(hits) == 1:
            hit = hits[0]
        else:
            hit = []
            for t in range(7):
                c = Counter([item[t] for item in hits if item[t] != ''])
                if len(c):
                    hit.append(c.most_common()[0][0])
                else:
                    hit.append("")
        result.append((sample, rid_i, merged_reads[read], read, good_kmer, found_minhd, hit))
        statistics["ok"] += 1
        statistics["hds"][found_minhd] += 1
        statistics["taxons"][tuple(hit)] += 1
        if not hit:
            statistics["failed"] += 1



    # for rid_i, read in enumerate(merged_reads):
    #     if rid_i % 10 == 0:
    #         print(sample, rid_i, len(merged_reads), f"{round((time.time()-start_time)/60,2)} min")
    #         # print(statistics["ok"], statistics["failed"])
    #         print(statistics["hds"].most_common(10))
    #         # print(statistics["taxons"].most_common(10))
    #         # print(statistics["names"].most_common(10))
    #         # save(start_time, prefix, result, statistics)

    #     L = len(read)
    #     hd = 10000000
    #     match = None
    #     for silvaid in len2silvaids[L]:
    #         hd_ = hamming_distance_raw(silvaid2seq[silvaid], read)
    #         if hd_ < hd:
    #             hd = hd_
    #             match = silvaid
    #             if hd == 0:
    #                 break
    #     statistics["hds"][hd] += 1
    #     result.append((rid_i, hd, match))

    # with open(f"{sample}.matches", "w") as fw:
    #     for rid_i, hd, match in result:
    #         data = [
    #             str(merged_reads[rid_i]),
    #             str(hd),
    #             str(silvaid2taxons[match]),
    #         ]
    #         fw.write("%s\n" % "\t".join(data))


    save(start_time, prefix, result, statistics)
    print("Done", f"{(time.time()-start_time)/60} min")

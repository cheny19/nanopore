#!/usr/bin/env python
"""
Created on Apr 10, 2015

@author: Chen Yang

This script generates simulated Oxford Nanopore 2D reads.

"""

import sys
import getopt
import random
from time import strftime
import numpy as np

import mixed_models as mm

PYTHON_VERSION = sys.version_info
VERSION = "1.0.0"
PRORAM = "Oxford Nanopore read simulator"
AUTHOR = "Chen Yang (UBC & BCGSC)"
CONTACT = "cheny@bcgsc.ca"

BASES = ['A', 'T', 'C', 'G']


# Check Python version, make sure it's 2.7 or higher
def check_version():
    if PYTHON_VERSION[0] == 2 and PYTHON_VERSION[1] < 7:
        print "You are using Python ", sys.version
        print "Please upgrade Python to 2.7 version or higher."
        sys.exit("Error!")
    return


# Usage information
def usage():
    print >> sys.stderr, "\npython simulator.py [command] <options>"
    print >> sys.stderr, "    [command] circular | linear"
    print >> sys.stderr, "    Do not choose 'circular' when there is more than one sequence in the reference"
    print >> sys.stderr, "    <options>: "
    print >> sys.stderr, "    -r: reference genome in fasta file, specify path and file name"
    print >> sys.stderr, "    -c : Flowcell chemistry, R7 or R7.3"
    print >> sys.stderr, "    -o : The prefix of output file, default = 'simulated'"
    print >> sys.stderr, "    -n : Number of generated reads, default = 20,000 reads"
    print >> sys.stderr, "    -s : Substitution profile, can be omitted if there is no customized profile"
    print >> sys.stderr, "    -p : Error model profile, can be omitted if there is no customized profile"
    print >> sys.stderr, "    -i : Insertion rate, a floating number in the interval [0, 1], default = 0.05"
    print >> sys.stderr, "    -d : Deletion rate, a floating number in the interval [0, 1], default = 0.05"
    print >> sys.stderr, "    -m : Mismatch rate, a floating number in the interval [0, 1], default = 0.1"


def read_ecdf(profile, lanes=1):
    # We need to count the number of zeros. If it's over 10 zeros, l_len/l_ratio need to be changed to higher.
    # Because it's almost impossible that the ratio is much lower than the lowest heuristic value.

    if lanes == 1:
        ecdf_dict = {}
        l_prob = 0.0
        l_len = 0.0
        for line in profile:
            new = line.strip().split('\t')
            length = [float(x) for x in new[0].split('-')]
            prob = float(new[1])
            if prob == l_prob:
                continue
            else:
                if l_prob != 0:
                    ecdf_dict[(l_prob, prob)] = (l_len, length[1])
                else:
                    ecdf_dict[(l_prob, prob)] = (max(l_len, length[1] - 10 * (length[1] - length[0])), length[1])
                l_prob = prob
                l_len = length[1]

    else:
        if lanes == 17:
            ecdf_dict = {(0, 1000): {}, (1000, 2000): {}, (2000, 3000): {}, (3000, 4000): {}, (4000, 5000): {},
                         (5000, 6000): {}, (6000, 7000): {}, (7000, 8000): {}, (8000, 9000): {}, (9000, 10000): {},
                         (10000, 11000): {}, (11000, 12000): {}, (12000, 13000): {}, (13000, 15000): {},
                         (15000, 20000): {}, (20000, 25000): {}, (25000, 50000): {}}
        elif lanes == 14:
            ecdf_dict = {(0, 10): {}, (10, 20): {}, (20, 30): {}, (30, 40): {}, (40, 50): {}, (50, 100): {},
                         (100, 300): {}, (300, 1000): {}, (1000, 2000): {}, (2000, 3000): {}, (3000, 5000): {},
                         (5000, 7500): {}, (7500, 10000): {}, (10000, 50000): {}}
        ecdf_key = sorted(ecdf_dict.keys())
        l_prob = [0.0] * lanes
        l_ratio = [0.0] * lanes

        for line in profile:
            new = line.strip().split('\t')
            ratio = [float(x) for x in new[0].split('-')]
            prob = [float(x) for x in new[1:]]
            for i in xrange(lanes):
                if prob[i] == l_prob[i]:
                    continue
                else:
                    if l_prob[i] != 0:
                        ecdf_dict[ecdf_key[i]][(l_prob[i], prob[i])] = (l_ratio[i], ratio[1])
                    else:
                        ecdf_dict[ecdf_key[i]][(l_prob[i], prob[i])] \
                            = (max(l_ratio[i], ratio[1] - 10 * (ratio[1] - ratio[0])), ratio[1])
                    l_ratio[i] = ratio[1]
                    l_prob[i] = prob[i]
    return ecdf_dict


def read_profile(number, chemistry, substitution_profile, model_profile):
    global unaligned_length, sub_matrix, ref_length
    global match_list, match_ht_list, align_ratio, ht_dict, error_par
    sub_matrix = [[0.0] * 5] * 5

    # Read substitution profile
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read substitution profile\n")
    if substitution_profile == "":
        substitution_profile = chemistry + "_substitution_profile"
    with open(substitution_profile, 'r') as sub_profile:
        for line in sub_profile:
            new = line.split('\t')
            if line[0] == "A":
                for i in xrange(5):
                    sub_matrix[0][i] = float(new[i + 1])
            elif line[0] == "T":
                for i in xrange(5):
                    sub_matrix[1][i] = float(new[i + 1])
            elif line[0] == "C":
                for i in xrange(5):
                    sub_matrix[2][i] = float(new[i + 1])
            elif line[0] == "G":
                for i in xrange(5):
                    sub_matrix[3][i] = float(new[i + 1])
            elif line[0] == "-":
                for i in xrange(5):
                    sub_matrix[4][i] = float(new[i + 1])

    # Read model profile for mismatch, insertion and deletions
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read error distribution profile\n")
    error_par = {}
    if model_profile == "":
        model_profile = chemistry + "_model_profile"
    with open(model_profile, 'r') as mod_profile:
        mod_profile.readline()
        for line in mod_profile:
            new_line = line.strip().split("\t")
            if "mismatch" in line:
                error_par["mis"] = [float(x) for x in new_line[1:]]
            elif "insertion" in line:
                error_par["ins"] = [float(x) for x in new_line[1:]]
            else:
                error_par["del"] = [float(x) for x in new_line[1:]]

    with open(chemistry + "_match.hist", 'r') as m_profile:
        match_hist = [0]
        match_ht_hist = [0]
        for line in m_profile:
            if line == "First and last match\n":
                break
            else:
                new = line.strip().split()
                match_hist.append(int(new[1]) + match_hist[-1])
        for line in m_profile:
            new = line.strip().split()
            match_ht_hist.append(int(new[1]) + match_ht_hist[-1])
    match_hist.remove(0)
    match_ht_hist.remove(0)

    match_list = {}
    match_ht_list = {}
    match_total = match_hist[-1]
    match_ht_total = match_ht_hist[-1]

    last = 0
    for i in xrange(len(match_hist)):
        match_hist[i] = float(match_hist[i]) / match_total
        match_list[(last, match_hist[i])] = i
        last = match_hist[i]
    del match_hist

    last = 0
    for i in xrange(len(match_ht_hist)):
        match_ht_hist[i] = float(match_ht_hist[i]) / match_ht_total
        if match_ht_hist[i] != 0:
            match_ht_list[(last, match_ht_hist[i])] = i
        last = match_ht_hist[i]
    del match_ht_hist

    # Read length of unaligned reads
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read ECDF of unaligned reads\n")
    unaligned_length = []
    with open(chemistry + "_unaligned_length_ecdf", 'r') as u_profile:
        new = u_profile.readline()
        rate = float(new.split('\t')[1])
        number_aligned = int(round(number * rate / (rate + 1)))
        number_unaligned = number - number_aligned
        unaligned_dict = read_ecdf(u_profile)

    for i in xrange(number_unaligned):
        p = random.random()
        for k_p, v_p in unaligned_dict.items():
            if k_p[0] <= p < k_p[1]:
                # consider this small range is linearly distributed:
                unaligned = (p - k_p[0])/(k_p[1] - k_p[0]) * (v_p[1] - v_p[0]) + v_p[0]
                unaligned_length.append(int(round(unaligned)))
                break

    unaligned_dict.clear()
    # Read profile of aligned reads
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read ECDF of aligned reads\n")

    # Read align ratio profile
    with open(chemistry + "_align_ratio", 'r') as a_profile:
        align_ratio = read_ecdf(a_profile, 17)

    # Read head/unaligned region ratio
    with open(chemistry + "_ht_ratio", 'r') as ht_profile:
        ht_dict = read_ecdf(ht_profile, 14)

    # Read length of aligned reads
    with open(chemistry + "_aligned_length_ecdf", 'r') as align_profile:
        aligned_dict = read_ecdf(align_profile)
    ref_length = []
    for i in xrange(number_aligned):
        middle_ref = 0
        while middle_ref < 80:
            p = random.random()
            for k_p, v_p in aligned_dict.items():
                if k_p[0] <= p < k_p[1]:
                    middle_ref = int(round((p - k_p[0])/(k_p[1] - k_p[0]) * (v_p[1] - v_p[0]) + v_p[0]))
                    break
        ref_length.append(middle_ref)

    aligned_dict.clear()


def simulation(ref, out, dna_type):
    global unaligned_length, ref_length, sub_matrix
    global genome_len, seq_dict, seq_len
    global match_list, match_ht_list, align_ratio, ht_dict, error_par

    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read in reference genome\n")
    seq_dict = {}
    seq_len = {}

    # Read in the reference genome
    with open(ref, 'r') as infile:
        for line in infile:
            if line[0] == ">":
                new_line = line.strip()[1:].split()
                chr_name = "-".join(new_line)
            else:
                if chr_name in seq_dict:
                    seq_dict[chr_name] += line.strip()
                else:
                    seq_dict[chr_name] = line.strip()

    if len(seq_dict) > 1 and dna_type == "circular":
        print >> sys.stderr, "Do not choose circular if there is more than one chromosome in the genome"
        sys.exit(1)

    for key in seq_dict.keys():
        seq_len[key] = len(seq_dict[key])
    genome_len = sum(seq_len.values())

    # Start simulation
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Start simulation of random reads\n")
    out_reads = open(out + "_reads.fasta", 'w')
    out_error = open(out + "_error_profile", 'w')
    out_test = open(out + "_matches", 'w')
    out_test2 = open(out + "_start_end_matches", 'w')
    out_error.write("Seq_name\tSeq_pos\terror_type\terror_length\tref_base\tseq_base\n")

    # Simulate random reads
    for i in xrange(len(unaligned_length)):
        unaligned = unaligned_length[i]
        unaligned, error_dict = unaligned_error_list(unaligned, error_par)
        new_read, new_read_name = extract_read(dna_type, unaligned)
        out_reads.write(">" + new_read_name + "_" + str(unaligned) + "-" + str(i) + "\n")
        read_mutated = mutate_read(new_read, new_read_name, out_error, error_dict)
        out_reads.write(read_mutated + "\n")
    del unaligned_length

    middle_length = []
    aligned_length = []
    middle_all_ratio = []
    remainder_length = []
    head_length = []
    tail_length = []

    # Simulate aligned reads
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Start simulation of aligned reads\n")
    for i in xrange(len(ref_length)):
        middle, middle_ref, error_dict = error_list(ref_length[i], match_list, match_ht_list, error_par)

        ref_length[i] = middle_ref
        middle_length.append(middle)

        for k_align in sorted(align_ratio.keys()):
            if k_align[0] <= middle < k_align[1]:
                break

        total = 0
        while total < 400:
            p = random.random()
            for k_r, v_r in align_ratio[k_align].items():
                if k_r[0] <= p < k_r[1]:
                    ratio = (p - k_r[0])/(k_r[1] - k_r[0]) * (v_r[1] - v_r[0]) + v_r[0]
                    total = int(round(middle / ratio))
                    remainder = total - int(round(middle))
                    break
        aligned_length.append(total)
        middle_all_ratio.append(ratio)
        remainder_length.append(remainder)

        for k_ht in sorted(ht_dict.keys()):
            if k_ht[0] <= remainder < k_ht[1]:
                break
        if remainder == 0:
            head = 0
            tail = 0
        else:
            p = random.random()
            for k_h, v_h in ht_dict[k_ht].items():
                if k_h[0] <= p < k_h[1]:
                    ratio = (p - k_h[0])/(k_h[1] - k_h[0]) * (v_h[1] - v_h[0]) + v_h[0]
                    head = int(round(remainder * ratio))
                    tail = remainder - head
                    break
        head_length.append(head)
        tail_length.append(tail)

        # Extract middle region from reference genome
        new_read, new_read_name = extract_read(dna_type, middle_ref)

        # Mutate read
        read_mutated = mutate_read(new_read, new_read_name, out_error, error_dict)

        # Add head and tail region
        for x in xrange(head):
            new_base = random.choice(BASES)
            read_mutated = new_base + read_mutated

        for x in xrange(tail):
            new_base = random.choice(BASES)
            read_mutated = read_mutated + new_base

        p = random.random()
        if p < 0.5:
            read_mutated = reverse_complement(read_mutated)

        out_reads.write(">" + new_read_name + "_" + str(head) + "_" + str(middle) + "_" +
                        str(tail) + "-" + str(i) + "\n")
        out_reads.write(read_mutated + "\n")

    out_reads.close()
    out_error.close()
    out_test.close()
    out_test2.close()

    align_ratio.clear()
    ht_dict.clear()

    o1 = open("head", 'w')
    o2 = open("middle", 'w')
    o3 = open("tail", 'w')
    o4 = open("aligned", 'w')
    o5 = open("ht", 'w')
    o6 = open("ratio", 'w')
    o7 = open("middle_ref", 'w')

    o1.write("\n".join(str(x) for x in head_length))
    o2.write("\n".join(str(x) for x in middle_length))
    o3.write("\n".join(str(x) for x in tail_length))
    o4.write("\n".join(str(x) for x in aligned_length))
    o5.write("\n".join(str(x) for x in remainder_length))
    o6.write("\n".join(str(x) for x in middle_all_ratio))
    o7.write("\n".join(str(x) for x in ref_length))

    o1.close()
    o2.close()
    o3.close()
    o4.close()
    o5.close()
    o6.close()
    o7.close()


def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    seq_list = list(seq)
    reverse_seq_list = reversed([comp.get(base, base) for base in seq_list])
    reverse_seq = ''.join(reverse_seq_list)
    return reverse_seq


def extract_read(dna_type, length):
    global seq_dict, seq_len, genome_len

    # Extract the aligned region from reference
    if dna_type == "circular":
        ref_pos = random.randint(0, genome_len)
        chromosome = seq_dict.keys()[0]
        new_read_name = chromosome + "_" + str(ref_pos)

        if length + ref_pos <= genome_len:
            new_read = seq_dict[chromosome][ref_pos: ref_pos + length]
        else:
            new_read = seq_dict[chromosome][ref_pos:]
            new_read = new_read + seq_dict[chromosome][0: length - genome_len + ref_pos]
    else:
        # Generate a random number within the size of the genome. Suppose chromosomes are connected
        # tail to head one by one in the order of the dictionary. If the start position fits in one
        # chromosome, but the end position does not, then restart generating random number.
        while True:
            new_read = ""
            ref_pos = random.randint(0, genome_len)
            for key in seq_len.keys():
                if ref_pos + length <= seq_len[key]:
                    new_read = seq_dict[key][ref_pos: ref_pos + length]
                    new_read_name = key + "_" + str(ref_pos)
                    break
                elif ref_pos < seq_len[key]:
                    break
                else:
                    ref_pos -= seq_len[key]
            if new_read != "":
                break
    return new_read, new_read_name


def unaligned_error_list(length, error_p):
    e_dict = {}
    error_rate = {(0, 0.4): "match", (0.4, 0.7): "mis", (0.7, 0.85): "ins", (0.85, 1): "del"}
    pos = 0
    last_is_ins = False
    while pos < length:
        p = random.random()
        for k_error in error_rate.keys():
            if k_error[0] <= p < k_error[1]:
                error_type = error_rate[k_error]
                break

        if error_type == "match":
            step = 1

        elif error_type == "mis":
            step = mm.pois_geom(error_p["mis"][0], error_p["mis"][2], error_p["mis"][3])
            e_dict[pos] = ["mis", step]

        elif error_type == "ins":
            step = mm.wei_geom(error_p["ins"][0], error_p["ins"][1], error_p["ins"][2], error_p["ins"][3])
            if last_is_ins:
                e_dict[pos + 0.1][1] += step
            else:
                e_dict[pos + 0.1] = ["ins", step]
                last_is_ins = True

        else:
            step = mm.wei_geom(error_p["del"][0], error_p["del"][1], error_p["del"][2], error_p["del"][3])
            e_dict[pos] = ["del", step]

        if error_type != "ins":
            pos += step
            last_is_ins = False

        if pos > length:
            length = pos

    return length, e_dict


def error_list(m_ref, m_list, m_ht_list, error_p):
    # l_old is the original length, and l_new is used to control the new length after introducing errors
    l_new = m_ref
    pos = 0
    e_dict = {}
    errors = {(0, 0.51027): "mis", (0.51027, 0.72467): "ins", (0.72467, 1): "del"}
    transition_pr = {"mis": {(0, 0.50105): "mis", (0.50105, 0.72018): "ins", (0.72018, 1): "del"},
                     "ins": {(0, 0.52186): "mis", (0.52186, 0.82170): "ins", (0.82170, 1): "del"},
                     "del": {(0, 0.51752): "mis", (0.51752, 0.65589): "ins", (0.65589, 1): "del"}}
    middle_ref = m_ref
    last_error = ""

    # The first match and last match come from m_ht_list
    p = random.random()
    for k in m_ht_list.keys():
        if k[0] <= p < k[1]:
            step = m_ht_list[k]
            break
    pos += step

    p = random.random()
    for k in m_ht_list.keys():
        if k[0] <= p < k[1]:
            last_match = m_ht_list[k]
            break

    # Select an error, then the step size, and then a match and so on so forth.
    while pos < middle_ref - last_match:
        if last_error == "":
            # the first error of a read is randomly selected based on the rate of different errors.
            p = random.random()
            for k in errors.keys():
                if k[0] <= p < k[1]:
                    error = errors[k]
                    break
        else:
            if step != 0:
                # the rest errors are selected based on Markov chain
                p = random.random()
                for k in transition_pr[last_error].keys():
                    if k[0] <= p < k[1]:
                        error = transition_pr[last_error][k]
                        break
            # if there are two consecutive errors, if the first one is mis, the second can be ins or del
            elif last_error == "mis":
                p = random.random()
                if p <= 0.44386:
                    error = "ins"
                else:
                    error = "del"

            # if the first one is ins or del, the following one can only be a mis
            elif last_error in ["ins", "del"]:
                error = "mis"

        if error == "mis":
            step = mm.pois_geom(error_p["mis"][0], error_p["mis"][2], error_p["mis"][3])
        else:
            step = mm.wei_geom(error_p[error][0], error_p[error][1], error_p[error][2], error_p[error][3])
            if error == "ins":
                l_new += step
            else:
                l_new -= step

        if error != "ins":
            e_dict[pos] = [error, step]
            pos += step
            if pos >= middle_ref - last_match:
                l_new += pos + last_match - middle_ref
                middle_ref = pos + last_match
                break
        else:
            e_dict[pos - 0.5] = [error, step]
            if pos == middle_ref - last_match:
                break

        last_error = error

        # Randomly select a match length
        p = random.random()

        for k in m_list.keys():
            if k[0] <= p < k[1]:
                step = m_list[k]
                break

        if pos + step > middle_ref - last_match:
            l_new += pos + step + last_match - middle_ref
            middle_ref = pos + step + last_match

        pos += step
    return l_new, middle_ref, e_dict


def mutate_read(read, read_name, error_log, e_dict):
    # global sub_matrix
    new_read = read

    for key in sorted(e_dict.keys(), reverse=True):
        val = e_dict[key]
        new_bases = ""

        if val[0] == "mis":
            ref_base = read[key: key + val[1]]
            new_read = read[: key]
            for i in xrange(val[1]):
                while True:
                    new_base = random.choice(BASES)
                    if new_base != read[key + i]:
                        new_read += new_base
                        new_bases += new_base
                        break
                    else:
                        continue
            new_read += read[key + val[1]:]

        elif val[0] == "del":
            ref_base = read[key: key + val[1]]
            new_read = read[: key] + read[key + val[1]:]
            new_bases = val[1] * "-"

        elif val[0] == "ins":
            key = int(round(key))
            ref_base = val[1] * "-"
            new_read = read[: key]
            for i in xrange(val[1]):
                new_base = random.choice(BASES)
                new_read += new_base
                new_bases += new_base
            new_read += read[key:]

        read = new_read

        if val[0] != "match":
            error_log.write(read_name + "\t" + str(key) + "\t" + val[0] + "\t" + str(val[1]) +
                            "\t" + ref_base + "\t" + new_bases + "\n")
    return new_read


def main():
    check_version()

    ref = ""
    chemistry = ""
    out = "simulated"
    number = 20000
    substitution_profile = ""
    model_profile = ""
    ins_rate = 0.0676
    del_rate = 0.0901
    mis_rate = 0.1441

    # Parse options and parameters
    if len(sys.argv) < 6:
        usage()
        sys.exit(2)
    else:
        dna_type = sys.argv[1]
        if dna_type not in ["circular", "linear"]:
            usage()
        try:
            opts, args = getopt.getopt(sys.argv[2:], "hr:c:o:n:l:i:d:m:s:p:")
        except getopt.GetoptError:
            usage()
            sys.exit(2)
        for opt, arg in opts:
            if opt == "-r":
                ref = arg
            elif opt == "-c":
                chemistry = arg
                if chemistry not in ("R7", "R7.3"):
                    usage()
                    sys.exit(2)
            elif opt == "-o":
                out = arg
            elif opt == "-n":
                number = int(arg)
            elif opt == "-s":
                substitution_profile = arg
            elif opt == "-p":
                model_profile = arg
            elif opt == "-i":
                ins_rate = float(arg)
            elif opt == "-d":
                del_rate = float(arg)
            elif opt == "-m":
                mis_rate = float(arg)
            elif opt == "-h":
                print "python simulator.py circular|linear -r <reference genome> -c <flowcell chemistry> " \
                      "-o <output prefix> -n <number of simulated reads>"

    # Generate log file
    sys.stdout = open(out + ".log", 'w')
    # Record the command typed to log file
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ' '.join(sys.argv) + '\n')

    # Read in reference genome and generate simulated reads
    read_profile(number, chemistry, substitution_profile, model_profile)

    simulation(ref, out, dna_type)

    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished!")
    sys.stdout.close()

if __name__ == "__main__":
    main()

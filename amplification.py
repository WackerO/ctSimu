import argparse
import copy 
import gzip
import pysam
import random
import sys


alphabet = ["A", "C", "G", "T"]                 #change this list to allow more/less letters for mutations (e.g. add "N" for undefined bases)

#mutates and returns a read with a given error rate and from a given alphabet; mutations are also written to a log file
def mutate_read(read, error_rate, mut_count, log, copies, pre_ampli=False):
    if not pre_ampli and random.random() < error_rate*len(read.sequence) or pre_ampli and random.random() < error_rate:     #only roll once for the whole sequence length, not once for every position in the sequence
        new_read = copy.deepcopy(read)                  #for an error before amplification, the probability is absolute, not relative to the sequence length
        seq, old_letter_s, new_letter_s, index_s, mut_count = mutate_sequence(new_read.sequence, mut_count)
        qual = new_read.quality
        if new_letter_s == "N":
            qual = qual[:index_s] + "!" + qual[index_s+1:]
        else:
            qual = qual[:index_s] + random.choice([",", ":", "?"]) + qual[index_s+1:]
        new_read.sequence = seq
        new_read.quality = qual
        log.write(new_read.name + "\t" + str(copies) + "\t" + str(index_s) + "\t" + old_letter_s + "\t" + new_letter_s + ("\tBefore\n" if pre_ampli else "\tAfter\n"))
        return new_read, mut_count
    else:
        return read, mut_count

#mutates and returns a sequence (String!), helper function for mutate_read
def mutate_sequence(seq, mut_count):
    alph = copy.deepcopy(alphabet)
    mut_count += 1
    index = random.choice(list(range(len(seq))))    #determine the position to be mutated
    if seq[index] in alph:
        alph.remove(seq[index])                     #removes the base that is at the chosen position in the original sequence (to prevent e.g. mutation A -> A)
    old_letter = seq[index]
    new_letter = random.choice(alph)
    seq = seq[:index] + new_letter + seq[index+1:]
    return seq, old_letter, new_letter, index, mut_count


#randomly amplifies a forward and reverse read pair up to multi_maximum times (or does not amplify and only returns a single pair of forward/reverse);
#output: [forward_copy1, forward_copy2, ..., forward_copyN], [reverse_copy1, reverse_copy2, ..., reverse_copyN], mutation_count
def amplify_read(f_read, r_read, preerror_rate, mut_count, logf, logr, multi_maximum):
    f_new_read, r_new_read = copy.deepcopy(f_read), copy.deepcopy(r_read)
    f_list, r_list = [], []
    rand = random.random()
    if rand < 0.25:                                 #decides how many copies of the current read pair are generated
        mult = 1
    elif rand < 0.5:
        mult = 2
    elif rand < 0.75:
        mult = 3
    else:
        mult = random.choice(list(range(4, multi_maximum+1)))
    mult = min(mult, multi_maximum)
    f_new_read, mut_count = mutate_read(f_new_read, preerror_rate, mut_count, logf, mult, True)
    r_new_read, mut_count = mutate_read(r_new_read, preerror_rate, mut_count, logr, mult, True)
    for i in range(mult):
        f = copy.deepcopy(f_new_read)
        f.name = f.name.split(":")[0] + "_" + str(i) + ":" + f.name.split(":")[1]
        r = copy.deepcopy(r_new_read)
        r.name = r.name.split(":")[0] + "_" + str(i) + ":" + r.name.split(":")[1]
        f_list.append(f)
        r_list.append(r)
    return f_list, r_list, mut_count

#mutates and multiplies the reads in an input FASTQ.gz to simulate PCR amplification and sequencing.
def write_fastq(forward, reverse, forward_out, reverse_out, forward_log, reverse_log, multi_maximum, error_rate, preerror_rate):
    with pysam.FastxFile(forward, 'rb') as frw, pysam.FastxFile(reverse, 'rb') as rev, \
    gzip.open(forward_out, 'ab') as outf, gzip.open(reverse_out, 'ab') as outr, \
    open(forward_log, 'a') as logf, open(reverse_log, 'a') as logr:      
        read_count, mut_count = 0, 0
        header = "#Log file for the amplification of " + forward + " and " + reverse +".\n" + "#read name\t" + "number of copies (including this read)\t" + \
                 "0-based index of mutation inside read\t" + "original base\t" + "mutated base\t" + \
                 "Mutation added before amplification? (If so, should persist even after dedup)\n"
        logf.write(header)
        logr.write(header)
        for f, r in zip(frw, rev):
            if f.name.split(":")[0] != r.name.split(":")[0]:
                sys.exit("Different read names", f.name, r.name)
            f_sub_list, r_sub_list, mut_count = amplify_read(f, r, preerror_rate, mut_count, logf, logr, multi_maximum)
            for f_sub in f_sub_list:        #iterate over all copies of the current read
                read_count += 1
                f_new, mut_count = mutate_read(f_sub, error_rate, mut_count, logf, len(f_sub_list))     #each copy has a chance of being mutated
                if (str(f_new) == ""):
                    sys.exit("Empty read, exiting...")
                outf.write((str(f_new)+"\n").encode())
            for r_sub in r_sub_list:
                read_count += 1
                r_new, mut_count = mutate_read(r_sub, error_rate, mut_count, logr, len(r_sub_list))
                if (str(r_new) == ""):
                    sys.exit("Empty read, exiting...")
                outr.write((str(r_new)+"\n").encode())
        print("Introduced %s mutations, %s reads in total." %(mut_count, read_count))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Mutates and multiplies the reads in an input FASTQ.gz to simulate PCR amplification and sequencing.')
    parser.add_argument('-f', '--forward', type=str, help='Input forward FASTQ.gz file to which mutations should be added.')
    parser.add_argument('-r', '--reverse', type=str, help='Input reverse FASTQ.gz file to which mutations should be added.')
    parser.add_argument('-fo', '--forward_out', type=str, help='Output forward FASTQ.gz file to which mutations should be added.')
    parser.add_argument('-ro', '--reverse_out', type=str, help='Output reverse FASTQ.gz file to which mutations should be added.')
    parser.add_argument('-fl', '--forward_log', type=str, default=False, help='Forward log TXT file in which mutations should be documented (if omitted, creates file in same directory as forward input).')
    parser.add_argument('-rl', '--reverse_log', type=str, default=False, help='Reverse log TXT file in which mutations should be documented (if omitted, creates file in same directory as reverse input).')
    parser.add_argument('-m', '--multi', type=int, default=24, help='Maximal duplication number; default: 24 (i.e. a read will be written between 1 and 24 times in the output file.')
    parser.add_argument('-e', '--error', type=float, default=0.0011128, help='Error rate (combined amplification or sequencing error, i.e. can change to A, C, G, T).')
    parser.add_argument('-p', '--preerror', type=float, default=0.001, help='Additional error rate for error before amplification (error will exist in all copies; can change to A, C, G, T).')

    args = parser.parse_args()
    if args.multi <= 0:
        sys.exit("Invalid maximal duplication number, exiting...")
    if (not args.forward.endswith("fastq.gz") and not args.forward.endswith("fq.gz") or 
        not args.reverse.endswith("fastq.gz") and not args.reverse.endswith("fq.gz")):
        sys.exit("Incorrect file format for FASTQ, exiting...")
    forward_log = args.forward_log or ".".join(args.forward.split(".")[:-1]) + "_log.txt"
    reverse_log = args.reverse_log or ".".join(args.reverse.split(".")[:-1]) + "_log.txt"
    write_fastq(args.forward, args.reverse, args.forward_out, args.reverse_out, \
                forward_log, reverse_log, args.multi, args.error, args.preerror)
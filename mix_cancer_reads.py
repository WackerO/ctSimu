import argparse
import copy
import gzip
import pysam
import random
import sys


#parses all mutations from a VCF into a dict: {chromosome_positions: AF, ...}, e.g. {chr1_12096347_C_T: 0.0436, chr3_9916415_T_A: 0.0373, ...}
def parse_mutation_positions(mutation_file):
    with open(mutation_file) as mut_file:
        lines = [l.strip().split("\t") for l in mut_file.readlines() if not l.startswith("#")]
        mutation_dict = {l[0]+"_"+l[1]+"_"+l[3]+"_"+l[4]: l[9].split(":")[1] for l in lines if not "0/1" in l and not "1/1" in l}       #VCF is 1-based
        return mutation_dict

#collects all read names from a golden BAM aligning to one of the VCF variants into a dict, e.g. {chr1_1123155_C_A_0.0475: set(readname1, readname2, readname3, ...), ...}
def get_vcf_positions(bam_file, mutation_file):
    pos_dict = {}
    mutation_dict = parse_mutation_positions(mutation_file)
    with pysam.AlignmentFile(bam_file) as bam:
        for mut in mutation_dict.keys():
            chrom, base = mut.split("_")[0], mut.split("_")[1]
            AF = mutation_dict[mut]
            for pileupcolumn in bam.pileup(chrom, int(base), int(base)+1):                  #https://pysam.readthedocs.io/en/latest/api.html
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:                 #https://stackoverflow.com/questions/39679380/bam-file-getting-all-reads-on-certain-positions-with-pysam
                        new_key = mut+"_"+AF
                        new_key = new_key.strip()
                        if new_key not in pos_dict.keys():
                            pos_dict[new_key] = set()
                        pos_dict[new_key].add(pileupread.alignment.query_name)
    return pos_dict

#gets the positions of all reads from a golden BAM, e.g. {readname1: chr1_125534, readname2: chr1_985141, ...}
def get_read_positions(bam_file):
    pos_dict = {}
    with pysam.AlignmentFile(bam_file) as bam:
        for read in bam.fetch():
            pos_dict[read.query_name] = read.reference_name+"_"+str(read.reference_start+1)
    return pos_dict

#adds mapping position to a read name, e.g. chr1_23124122_readname
def add_mapping(read, name_dict):
    if read.name.split("/")[0] in name_dict.keys():
        read.name = name_dict[read.name.split("/")[0]] + "_" + read.name
    else:
        sys.exit("Could not determine mapping position of read %s." %(read.name))
    return read

#appends UMI sequences to the names of reads for dedup, e.g. readname:ATACAAGG,TTCGAGTC
def add_umis(f, r, umis, umi_pointer1, umi_pointer2):
    combi = umis[umi_pointer1]+","+umis[umi_pointer2]
    f.name, r.name = f.name.rstrip().split("/")[0]+":"+combi+"/"+f.name.rstrip().split("/")[1], r.name.rstrip().split("/")[0]+":"+combi+"/"+r.name.rstrip().split("/")[1]
    return f, r

#changes the UMI pointers as necessary to cycle through all combinations of UMI sequences
def bump_pointers(umis, umi_pointer1, umi_pointer2):
    umi_pointer2 += 1
    if umi_pointer2 >= len(umis):
        umi_pointer2 = 0
        umi_pointer1 += 1
        if umi_pointer1 >= len(umis):
            umi_pointer1 = 0
    return umi_pointer1, umi_pointer2

#writes a true VCF containing for each variant the AF which was actually written into the output FASTQs
def write_vcf(truth, af_healthy_dict, af_cancer_dict):
    with open(truth, "w") as out:
        out.write("##source=mySimuPipeline_mix_cancer_reads\n" +
                  "##INFO=<ID=CONC,Number=1,Type=Float,Description=\"Mixed AF/concentration of variant\">\n" +
                  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")
        for mut in af_cancer_dict.keys():
            if mut not in af_healthy_dict.keys():
                AF = 1
            else:
                AF = af_cancer_dict[mut]/(af_cancer_dict[mut]+af_healthy_dict[mut]) #calculate ratio of variant reads to total reads per variant position, this is the theoretically introduced "true" AF
            line = mut.split("_")[0]+"\t"+mut.split("_")[1] + "\tM\t" + mut.split("_")[2]+"\t"+mut.split("_")[3] + "\t---\tPASS\t---\tGT:AF\t./.:"+str(AF)+"\n"
            out.write(line)

#mixes cancer reads into healthy read FASTQ files and also adds UMIs to the IDs
def write_fastq(forward, reverse, forward_cancer, reverse_cancer, healthy_golden, cancer_golden, forward_out, reverse_out, umifile, mutations, truth):
    umi_pointer1, umi_pointer2 = 0, 0                                       #point to first and second UMI sequence to generate combinations, each 8 BP long in IDT
    umis = []
    healthy_name_dict = get_read_positions(healthy_golden)
    cancer_name_dict = get_read_positions(cancer_golden)
    healthy_relevant_dict = get_vcf_positions(healthy_golden, mutations)    #holds the names of all healthy reads mapped to a variant positions in lists: {chr1_1123155_C_A_0.0475: set(readname1, readname2, readname3, ...), ...}
    mut_relevant_dict = get_vcf_positions(cancer_golden, mutations)         #holds the names of all cancer reads mapped to a variant positions in lists: {chr1_1123155_C_A_0.0475: set(readname1, readname2, readname3, ...), ...}
    count_dict = {}                                                         #saves the necessary number of mutated reads at each variant location: {chr1_1141782: 34, chr3_4223452: 55, ...}
    af_healthy_dict, af_cancer_dict = {}, {}                                #saves the number of variant location reads written as healthy and as cancer (for the true VCF)
    with open(umifile, "r") as u:
        umis = [l.strip() for l in u.readlines()]

    with pysam.FastxFile(forward, 'rb') as frw, pysam.FastxFile(reverse, 'rb') as rev, \
    pysam.FastxFile(forward_cancer, 'rb') as frwt, pysam.FastxFile(reverse_cancer, 'rb') as revt, \
    gzip.open(forward_out, 'ab') as outf, gzip.open(reverse_out, 'ab') as outr:    
    #this loop deals with writing (or skipping) the reads from the healthy dataset
        for f, r in zip(frw, rev):
            write = True                        #boolean switch for whether to write the current healthy read pair or not (if False, don't write healthy, instead write tumor read pair)
            if f.name.split("/")[0] != r.name.split("/")[0]:
                sys.exit("IDs of parallel healthy reads %s and %s do not agree, exiting..." %(f.name, r.name))
            if f.sequence != r.sequence:        #make sure that forward and reverse read do not have an identical sequence
                for key in healthy_relevant_dict.keys():     #check if current read pair is relevant (i.e. one or both reads map to any variant positions)
                    percentage = float(key.split("_")[-1].rstrip())
                    if percentage > 0 and (f.name.split("/")[0] in healthy_relevant_dict[key] or r.name.split("/")[0] in healthy_relevant_dict[key]):    
                        if random.random() < percentage:                        #don't write healthy reads, instead bump
                            write = False                                       #counter for this positions     
                        break
                if write:           #if current healthy reads were not replaced by mutation, write them into the output files
                    new_f, new_r = copy.deepcopy(f), copy.deepcopy(r)
                    new_f, new_r = add_mapping(new_f, healthy_name_dict), add_mapping(new_r, healthy_name_dict)
                    new_f, new_r = add_umis(new_f, new_r, umis, umi_pointer1, umi_pointer2)
                    umi_pointer1, umi_pointer2 = bump_pointers(umis, umi_pointer1, umi_pointer2)
                    outf.write((str(new_f)+"\n").encode())
                    outr.write((str(new_r)+"\n").encode())
                    for key in healthy_relevant_dict.keys():
                        if f.name.split("/")[0] in healthy_relevant_dict[key] or r.name.split("/")[0] in healthy_relevant_dict[key]:
                            if key not in af_healthy_dict.keys():
                                af_healthy_dict[key] = 1
                            else:
                                af_healthy_dict[key] += 1
                else:               #if current healthy reads were replaced by mutation, document this for the respective mutation
                    for key in healthy_relevant_dict.keys():
                        if f.name.split("/")[0] in healthy_relevant_dict[key] or r.name.split("/")[0] in healthy_relevant_dict[key]:
                            if key not in count_dict.keys():                                                    
                                count_dict[key] = 1
                            else:
                                count_dict[key] += 1
                    
        #this loop writes all required reads from the cancer dataset
        for ft, rt in zip(frwt, revt):
            write = False
            if ft.name.split("/")[0] != rt.name.split("/")[0]:
                sys.exit("IDs of parallel cancer reads %s and %s do not agree, exiting..." %(ft.name, rt.name))
            if ft.sequence != rt.sequence:        #make sure that forward and reverse read do not have an identical sequence
                for key in mut_relevant_dict.keys():
                    if (ft.name.split("/")[0] in mut_relevant_dict[key] or rt.name.split("/")[0] in mut_relevant_dict[key]) and key in count_dict.keys():     #check if the current read pair name maps to a variant positions (mut_relevant_dict)
                        write = True                                                                                                                          #and this variant position still needs mutated reads to be written (count_dict)
                        if key not in af_cancer_dict.keys():
                            af_cancer_dict[key] = 1
                        else:
                            af_cancer_dict[key] += 1
                        count_dict[key] -= 1            #at the genomic location of the current read pair, a mutation was introduced, therefore reduce the counter of necessary mutations at this location
                        if count_dict[key] < 1:
                            del count_dict[key]         #current genomic location does not need any more cancer reads (all skipped healthy reads have been replaced with cancer reads), delete it
            if write:
                write = False
                new_ft, new_rt = copy.deepcopy(ft), copy.deepcopy(rt)
                new_ft, new_rt = add_mapping(new_ft, cancer_name_dict), add_mapping(new_rt, cancer_name_dict)
                new_ft, new_rt = add_umis(new_ft, new_rt, umis, umi_pointer1, umi_pointer2)
                umi_pointer1, umi_pointer2 = bump_pointers(umis, umi_pointer1, umi_pointer2)
                outf.write((str(new_ft)+"\n").encode())
                outr.write((str(new_rt)+"\n").encode())
            
            if not count_dict:
                print("All skipped healthy reads have been replaced, breaking loop...")
                break                               #all skipped healthy reads have been replaced; the program is finished
        write_vcf(truth, af_healthy_dict, af_cancer_dict)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Mixes cancer reads into healthy read FASTQ files and also adds UMIs to the IDs.')
    parser.add_argument('-f', '--forward', type=str, help='Healthy FASTQ.GZ forward (R1) file to which the cancer reads should be added.')
    parser.add_argument('-r', '--reverse', type=str, help='Healthy FASTQ.GZ reverse (R2) file to which the cancer reads should be added.')
    parser.add_argument('-fc', '--forward_cancer', type=str, help='cancer FASTQ.GZ forward (R1) file from which the cancer reads should be taken.')
    parser.add_argument('-rc', '--reverse_cancer', type=str, help='cancer FASTQ.GZ reverse (R2) file from which the cancer reads should be taken.')
    parser.add_argument('-fo', '--forward_out', type=str, default=False, help='Output FASTQ.GZ forward (R1) file, if omitted, will save output in folder of -f forward file.')
    parser.add_argument('-ro', '--reverse_out', type=str, default=False, help='Output FASTQ.GZ reverse (R2) file, if omitted, will save output in folder of -r reverse file.')
    parser.add_argument('-hg', '--healthy_golden', type=str, help='Healthy golden BAM file to find relevant reads.')
    parser.add_argument('-cg', '--cancer_golden', type=str, help='cancer golden BAM file to find relevant reads.')
    parser.add_argument('-u', '--umis', type=str, help='TXT file with pool of UMI sequences.')
    parser.add_argument('-m', '--mutations', type=str, help='VCF file containing the relevant mutation positions.')
    parser.add_argument('-t', '--truth', type=str, help='True VCF file in which to save the true AF for each mutation.')

    args = parser.parse_args()
    if any(not f.endswith("fastq.gz") and not f.endswith("fq.gz") for f in [args.forward, args.reverse, args.forward_cancer, args.reverse_cancer]):
        sys.exit("Incorrect file format for FASTQ, exiting...")
    if any(not f.endswith(".bam") for f in [args.healthy_golden, args.cancer_golden]):
        sys.exit("Incorrect file format for golden BAM, exiting...")
    if not args.mutations.endswith(".vcf") or not args.truth.endswith(".vcf"):
        sys.exit("Incorrect file format for VCF, exiting...")
    forward_out = args.forward_out or args.forward.split(".")[0] + "_cancer_added_" + ".".join(args.forward.split(".")[1:])
    if not forward_out.endswith(".gz"):
        forward_out += ".gz"
    reverse_out = args.reverse_out or args.reverse.split(".")[0] + "_cancer_added_" + + ".".join(args.reverse.split(".")[1:])
    if not reverse_out.endswith(".gz"):
        reverse_out += ".gz"
    write_fastq(args.forward, args.reverse, args.forward_cancer, args.reverse_cancer, args.healthy_golden, args.cancer_golden, \
                forward_out, reverse_out, args.umis, args.mutations, args.truth)
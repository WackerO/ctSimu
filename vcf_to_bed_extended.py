import argparse
import sys


#converts VCF to BED and extends it
def write_bed(file, outfile, length, direction):
    with open(file, "r") as f, open(outfile, "w") as o:
        f = f.readlines()
        for l in f:
            if not l.startswith("#"):
                chrom, pos = l.split()[0], l.split()[1]
                start = int(pos)
                if 'u' in direction or 'b' in direction:
                    start -= length
                stop = int(pos)
                if 'd' in direction or 'b' in direction:
                    stop += length
                o.write(chrom + "\t" + str(start) + "\t" + str(stop) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Creates a BED file with additional positions around the mutations in the input VCF.')
    parser.add_argument('-i', '--infile', type=str, help='VCF file containing the mutations.')
    parser.add_argument('-o', '--outfile', type=str, default=False, help='Output BED file, if omitted, will save to folder of infile.')
    parser.add_argument('-l', '--length', type=int, default=0, help='Length in each chosen direction from each variant position, default 0 (start and end will be the same mutation position).')
    parser.add_argument('-d', '--direction', choices=['u', 'up', 'd', 'down', 'b', 'both'], default='b', help='In which directions to extend.')

    args = parser.parse_args()
    infile = args.infile
    if not infile.endswith("vcf") or args.outfile and not args.outfile.endswith("bed"):
        sys.exit("Incorrect file format, exiting...")
    outfile = args.outfile or infile.split(".")[0] + "_cancer" + ".bed"
    print(outfile)
    write_bed(infile, outfile, abs(args.length), args.direction)
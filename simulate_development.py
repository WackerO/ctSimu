import argparse
import os
import random
import sys


time_list = ["before", "therapy1", "therapy2", "therapy3", "after"] #list of timepoints; length has to be same as length of each trend
#the lists in the following dict are ordered! The first entry in each list corresponds to before etc., implementation by index!
trend_dict = {"up": [0.045, 0.048, 0.050, 0.051, 0.055],            #cancer worsens over time, therapy unsuccessful
              "down": [0.048, 0.040, 0.016, 0.001, 0.00001],        #cancer gets better, therapy successful
              "const": [0.045, 0.046, 0.045, 0.046, 0.045],         #cancer stays constant, therapy unsuccessful
              "rel": [0.045, 0.038, 0.021, 0.035, 0.043],           #relapse, cancer gets better, then worse again
              "delay": [0.048, 0.054, 0.026, 0.007, 0.00001]}       #delayed down (abrupt cell death causes initial ctDNA high)

def write_vcfs(infile, outfile, trend, point, variation):
    i = time_list.index(point)          #index of the time point (before: 0, therapy1: 1 etc.)
    average = trend_dict[trend][i]      #average value of the time point in the respective trend (0.045 etc.)
    with open(infile, "r") as f, open(outfile, "w") as out:
        mutations = ["\t".join(l.split("\t")[:6]) for l in f.readlines() if not l.startswith("#")] #get variants from VCF including column index 5 (QUAL)
        out.write("##source=mySimuPipeline_simulate_development\n" +
                  "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Desired/expected AF of variant\">\n" +
                  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" %(point))
        for mut in mutations:
            if "F!" in mut.split()[2]:   #full, the variant is set to an AF of 1 (homozygous)
                line = mut.replace("F!", "")+"\tPASS\t---\tGT:AF\t"
                line += "1/1:"+"1"+"\n"
                AF = 1
            elif "H!" in mut.split()[2]: #half, the variant is set to an AF of 0.5 (heterozygous)
                line = mut.replace("H!", "")+"\tPASS\t---\tGT:AF\t"
                line += "0/1:"+"0.5"+"\n"
                AF = 0.5
            else:                       #the variant is set to a random AF depending on the chosen trend
                line = mut+"\tPASS\t---\tGT:AF\t"
                AF = average+random.uniform(0, average*variation) * random.choice([-1, 1])
                line += "./.:"+str(AF)+"\n"
            if AF > 0:
                out.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Creates VCF files for a development of tumor variant AF of a patient at one timepoint of the following: before therapy, therapy 1, therapy 2, therapy 3, after therapy.')
    parser.add_argument('-i', '--infile', type=str, help='VCF file containing the mutations.')
    parser.add_argument('-o', '--outfile', type=str, default=False, help='Output file, if omitted, will write in same dir as input file.')
    parser.add_argument('-t', '--trend', choices=trend_dict.keys(), help='Type of trend of the patient\'s development; up=worsening, down=improvement, ' +
                                                                                      'const=no overall change (but variation), rel=relapse, i.e. initial improvement, then worsening.')
    parser.add_argument('-p', '--point', choices=time_list, help='Time point in therapy.')
    parser.add_argument('-v', '--variation', type=float, default=0.035, help='Biological variation of the AF around the respective average, default=0.01.')
    
    args = parser.parse_args()
    infile = args.infile
    if not infile.endswith("vcf"):
        sys.exit("Incorrect file format, must be VCF; exiting...")
    outfile = args.outfile or os.path.dirname(os.path.join(infile, args.point+".vcf"))
    if args.variation > 1:
        print("Variation above 1, dividing by 100.")
        variation /= 100
    else:
        variation = args.variation
    write_vcfs(infile, outfile, args.trend, args.point, variation)
import argparse
import os
from read_VCFs import read_vcf_line


#produces a dictionary containing all variants and their allele frequencies from a VCF: {chr1_934851_A_C: 0.00017, ...}
def make_call_dict(file):
    call_dict = {}
    with open(file, "r") as f:
        lines = f.read().splitlines() 
        for l in lines:
            if l.startswith("##source"):  #pick the VarCaller that was used
                source = "lofreq" if "lofreq" in l else "vardict" if "VarDict" in l else "umivar" if "umiVar" in l \
                         else "mySimuPipeline" if "mySimuPipeline_simulate_development" in l or "mySimuPipeline_mix_cancer_reads" in l else False
            if not l.startswith("#"):
                key, value = read_vcf_line(l, source)
                if value > 0:
                    call_dict[key] = value
    return call_dict

#compares the variant calls of the VCFs found in folder and subdirs with the true_file VCF, searches only for file names containing caller
def compare(folder, varcallers):
    varcaller_dict = {varcaller: {} for varcaller in varcallers}
    values = ["file_count", "total_true", "total_called", "tp", "fp", "fn"]     #these are the values that are saved for each variant caller
    for varcaller in varcallers:
        for value in values:
            varcaller_dict[varcaller][value] = 0                               #initialize all values to 0
    #find all variant calling files that have corresponding true files
    for path, subdirs, files in os.walk(folder):
        if "truth.vcf" in files:
            true_file = os.path.join(path, "truth.vcf")
            
            for name in files:
                call_file = os.path.join(path, name)
                if call_file.endswith(".vcf") and os.path.basename(path) != "umivar":
                    varcaller = os.path.splitext(name)[0]
                    if varcaller in varcallers:
                        varcaller_dict[varcaller]["file_count"] += 1
                        true_dict, call_dict =  make_call_dict(true_file), make_call_dict(call_file)
                        varcaller_dict[varcaller]["total_true"] += len(true_dict.keys())
                        varcaller_dict[varcaller]["total_called"] += len(call_dict.keys())
                        
                        for true_var in true_dict.keys():   #check which true variants were found in the variant calling
                            if true_var in call_dict.keys():
                                varcaller_dict[varcaller]["tp"] += 1
                            else:
                                varcaller_dict[varcaller]["fn"] += 1
                        for call_var in call_dict.keys():   #check which false positive variants were called
                            if call_var not in true_dict.keys():
                                varcaller_dict[varcaller]["fp"] += 1

    for varcaller in varcaller_dict.keys():
        print("In %s files, %s called %s variants out of %s true variants. %s calls were true positives (%s %%) and %s were false positives (%s %%), also %s true variants were not called (false negatives, %s %%). Sensitivity: %s, precision: %s\n\n" \
                      %(varcaller_dict[varcaller]["file_count"],
                      varcaller,
                      varcaller_dict[varcaller]["total_called"],
                      varcaller_dict[varcaller]["total_true"],
                      varcaller_dict[varcaller]["tp"],
                      100*varcaller_dict[varcaller]["tp"]/varcaller_dict[varcaller]["total_called"] if varcaller_dict[varcaller]["total_called"] != 0 else "(zero division, none called!)" ,
                      varcaller_dict[varcaller]["fp"],
                      100*varcaller_dict[varcaller]["fp"]/varcaller_dict[varcaller]["total_called"] if varcaller_dict[varcaller]["total_called"] != 0 else "(zero division, none called!)",
                      varcaller_dict[varcaller]["fn"],
                      100*varcaller_dict[varcaller]["fn"]/varcaller_dict[varcaller]["total_called"] if varcaller_dict[varcaller]["total_called"] != 0 else "(zero division, none called!)",
                      varcaller_dict[varcaller]["tp"]/(varcaller_dict[varcaller]["tp"]+varcaller_dict[varcaller]["fn"]),
                      varcaller_dict[varcaller]["tp"]/(varcaller_dict[varcaller]["tp"]+varcaller_dict[varcaller]["fp"])))

#compares VCF files from variant callers to their corresponding true VCF files in a directory
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compares VCF files from variant callers to their corresponding true VCF files in a directory.')
    parser.add_argument('-f', '--folder', type=str, help='Main folder containing the simulated and true VCFs (possibly in subdirectories)')
    parser.add_argument('-v', '--varcallers', type=str, nargs='+', help='List of variant caller names to be searched for in file names')

    args = parser.parse_args()
    compare(args.folder, args.varcallers)



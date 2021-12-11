import argparse
from classify_trend import classify_plot
import os
from plot_AF import get_median
from read_VCFs import read_files
import sys


trends = ["up", "down", "rel", "delay", "const"]                            #names of files for classification have to contain one of these trends
varcallers = ["lofreq", "umivar", "vardict"]                                #names of files for classification have to contain one of these variant callers
timepoints = ["before", "therapy1", "therapy2", "therapy3", "after"]        #names of files for classification have to contain one of these timepoints

#iterates over a number of different parameters for trend classification and finds the best combination
def train(patient_dict):
    ratio_dict = {}                                                         #save results for each combo of parameters to this dict
    for diff_factor_1 in range(10, 15):                                     #these stacked loops try out different combinations of parameters  
        for diff_factor_2 in range(4, 9):
            for count_u in [1, 2, 3]:
                for count_d in [1, 2, 3]:
                    correct_ratio = test(patient_dict, [diff_factor_1, diff_factor_2, count_u, count_d])
                    if correct_ratio not in ratio_dict:
                        ratio_dict[correct_ratio] = []
                    ratio_dict[correct_ratio].append([diff_factor_1, diff_factor_2, count_u, count_d])      #for current accuracy, save current combination to other
    best_val = max(ratio_dict.keys())                                                                       #combinations with the same results        
    print("The best results (%s) were achieved with these combos: %s" %(best_val, ratio_dict[best_val]))

#for one combination of parameters, gets the results for all patient files and calculates the ratio of correctly classified files
def test(main_dict, combo=[14, 4, 1, 1]):
    count_dict = {}                             #saves the number of incorrect predictions for each type (e.g. up-->constant: 12 means that 12 up files were incorrectly classified as constant) -->not necessary for training, only for modifying the classification
    correct, total = 0, 0                       #counters for correctly and total classified files
    for patient in main_dict.keys():
        trend = patient.split("_")[0]
        patient_dict = main_dict[patient]
        for varcaller in patient_dict.keys():
            files = [file for file in patient_dict[varcaller] if file is not None]
            var_dict = read_files(files)[0]
            if var_dict:
                medians = get_median(var_dict)
                prediction_string = classify_plot(medians, combo[0], combo[1], combo[2], combo[3])      #get prediction string, then determine which trend was chosen
                prediction = "up" if "Up" in prediction_string else "down" if "Down" in prediction_string \
                             else "rel" if "Relapse" in prediction_string else "delay" if "Delay" in prediction_string \
                             else "const" if "Constant" in prediction_string else "unknown" if "Unknown" in prediction_string else "HUH?"
                total += 1
                if prediction != trend:
                    error = trend + " --> " + prediction            #for incorrect predictions, formulate the error
                    if error not in count_dict:
                        count_dict[error] = 0
                    count_dict[error] += 1
                else:
                    correct += 1                                    #if correct, bump counter
                    if trend not in count_dict:
                        count_dict[trend] = 0
                    count_dict[trend] += 1
    return correct/total                                        #return ratio

#searches all files in the target folder and determines which files can be used for classification
def make_main_dict(main_folder):
    main_dict = {}                                              #saves eligible files
    for patient_folder in os.scandir(main_folder):
        if patient_folder.is_dir() and any(trend in patient_folder.name for trend in trends):
            patient_dict = {}
            trend = os.path.basename(patient_folder).split("_")[0]
            for sample_folder in os.scandir(patient_folder):
                if sample_folder.is_dir():
                    for sample_file in os.scandir(sample_folder):
                        parent_dir = os.path.basename(os.path.dirname(sample_file))
                        if sample_file.is_file() and sample_file.name.endswith(".vcf") and \
                           any(varcaller in sample_file.name for varcaller in varcallers) and \
                           any(timepoint in parent_dir for timepoint in timepoints):
                            varcaller = os.path.splitext(os.path.basename(sample_file))[0]
                            timepoint = parent_dir.split("_")[-1]
                            if varcaller not in patient_dict.keys():
                                patient_dict[varcaller] = [None]*len(timepoints)
                            patient_dict[varcaller][timepoints.index(timepoint)] = sample_file
            main_dict[patient_folder.name] = patient_dict
    return main_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Trains and or tests the classification of patient disease progress.')
    parser.add_argument('-f', '--folder', type=str, help='Folder containing the subfolders with VCF files.')
    parser.add_argument('-p', '--params', type=int, nargs='+', default=False, help='List of parameters (length 4) to test (e.g. 14 8 3 3); if omitted, will train.')

    args = parser.parse_args()
    main_dict = make_main_dict(args.folder)
    if args.params:
        if len(args.params) == 4:
            print("With parameters %s, %s %% of the predictions were correct." %(args.params, test(main_dict, args.params)*100))
        else:
            sys.exit("Params must be list of int with length 4, exiting...")
    else:
        train(main_dict)

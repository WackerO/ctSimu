import argparse
import math
import sys


#classifies a list of median values for time samples
def classify_plot(medians, diff_factor_1=14, diff_factor_2=4, count_u=1, count_d=1):
    medians = [median for median in medians if not math.isnan(median)]
    slope_string = ""
    for m in range(1, len(medians)):
        previous_median, current_median  = medians[m-1],  medians[m]
        difference = abs(previous_median-current_median)
        if difference >= max(previous_median, current_median)/diff_factor_1:   #is the abs. difference sufficient to declare the development non-constant?
            slope_string += "u" if previous_median < current_median else "d"
        else:
            slope_string += "c"

    start, end = medians[0], medians[-1]
    total_difference = abs(start-end)
    if total_difference >= max(start, end)/diff_factor_2:                      #is the total abs. difference sufficient to declare the development non-constant?
        if start > end:
            if "u" in slope_string:
                return "Delayed down, slope string: " + slope_string
            else:
                return "Down, slope string: " + slope_string
        else:
            if "d" in slope_string:
                return "Relapse, slope string: " + slope_string
            else:
                return "Upwards, slope string: " + slope_string
    elif "c" not in slope_string:
        if "d" in slope_string and slope_string.endswith("u"):
            return "Relapse, slope string: " + slope_string
        else:
            return "Delayed down, slope string: " + slope_string
    elif slope_string.count("u") < count_u and slope_string.count("d") < count_d:   #are too few u and d in the slope string for the trend to be up or down?
        return "Constant, slope string: " + slope_string
    elif "d" not in slope_string and (len(slope_string) > 2 and slope_string.count("d") > 0) or \
                                     (len(slope_string) > 3 and slope_string.count("d") > 1):
        return "Upwards, slope string: " + slope_string
    else:
        return "Unknown, slope string: " + slope_string

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Classifies a number of VCFs for a cancer patient as one of 5 therapy trends.')
    parser.add_argument('-i', '--infiles', type=str, nargs='+', help='VCF files in chronological order.')

    args = parser.parse_args()
    if any(not infile.endswith("vcf") for infile in args.infiles):
        sys.exit("Incorrect file format, must be VCF; exiting...")
    classify_plot(args.infiles)
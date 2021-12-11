import argparse
from classify_trend import classify_plot
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from read_VCFs import read_files, get_median
import sys


#creates and saves the line plot
def plots_af(var_dict, outfile, samples=False, bed=False, filtering=False, legend=False, \
             y_filter=False, y_force=False):
    matplotlib.rcParams.update({'font.size': 17})                           #increase font size for everything
    cmap = plt.get_cmap('winter')                                           #prepare colors and markers for the different lines
    sub_colors = [cmap(i) for i in np.linspace(0, 1, len(var_dict.keys()))] #subscheme for mixed plot (AF and biomarkers)
    markers = [x for x in matplotlib.markers.MarkerStyle.markers.keys()]    #different markers if only AF is plotted (to better distinguish variants)
    cmap = plt.get_cmap('gist_rainbow')
    full_colors = [cmap(i) for i in np.linspace(0, 1, len(var_dict.keys()))]#full scheme if only AF is plotted (to better distinguish variants)

    fig, ax = plt.subplots()                                    #create y axis for the line plots

    if y_filter:
        percentiles = []                                        #make list to save the percentiles into for filtering outliers

    x = np.arange(len(samples))                             #else use sample names to name x ticks
    fig.set_size_inches(10, 8)
    for i, var in enumerate(var_dict.keys()):
        j = i
        while j >= len(markers):
            j -= len(markers)                               #prevent index out of bounds error
        current_marker = markers[j]                         #select marker for current variant
        ax.plot(x, var_dict[var], color=full_colors[i], marker=current_marker, linestyle='-', label=('chr'+var.split('_')[0]+', '+var.split('_')[1]))
        if y_filter:
            percentiles.append(var_dict[var])                   #append value of current variant to list to to rescale y axis
    ax.set_ylabel('Allele frequency')
    ax.set_xlabel('Sample')
    
    if y_filter and not y_force:
        percentiles = [list(x) for x in zip(*percentiles)]      #transpose list to get list of values for each sample
        percentiles = [np.percentile(np.array(p), y_filter) for p in percentiles]   #find desired percentile for each sample
        max_y = max([p for p in percentiles if not np.isnan(p)])*1.1                #y maximum is 110 % of the highest percentile
        if max_y > 0.1:     #if above 0.1, round to make numbers cleaner
            max_y = round(max_y, 2)
        ax.set_ylim(ymax=max_y)
    ax.set_ylim(ymin=0)
    if y_force:
        ax.set_ylim(ymax=y_force)
    plt.xticks(x, samples)
    ax.yaxis.grid()
    
    median_list = get_median(var_dict)
    ax.plot(x, median_list, color='black', marker= 'o', linewidth=4, label='AF median')    #plot median line

    handles, labels = [], []
    for axis in fig.axes:
        axis_handles, axis_labels = axis.get_legend_handles_labels()
        handles.extend(axis_handles)
        labels.extend(axis_labels)

    if len(handles) > 0 and len(labels) > 0:                        #format legend (order and position of legend elements)
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: int(t[0].split(' ')[1]) if t[0].split(',')[0].lstrip('chr').isnumeric() \
                                                                          else 100 if 'X' in t[0].split(',')[0] else 101 if 'Y' in t[0].split(',')[0] else 102))
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: int(t[0].split(',')[0].lstrip('chr')) if t[0].split(',')[0].lstrip('chr').isnumeric() \
                                                                          else 100 if 'X' in t[0].split(',')[0] else 101 if 'Y' in t[0].split(',')[0] else 102))
        if len(labels) < 16 or legend:
            plt.legend(handles, labels, title='Chromosome, position', loc='upper left', bbox_to_anchor=(1.0, 1.0), ncol=max(1, round(len(handles)/15)))

    slope_string = classify_plot(median_list)                       #add classification to plot title
    plt.title('Allele frequency of ctDNA variants in time series;\n %s found %s variants. Classified as %s' %(os.path.splitext(os.path.basename(outfile))[0], len(var_dict.keys()), slope_string))

    if len(x) > 7:                                                  #if many samples, auto-tilt x names to prevent overlapping words
        fig.autofmt_xdate()
    plt.savefig(outfile, bbox_inches='tight')
    plt.clf()

#main, calls the other functions
def evaluate(files, samples, outfile, bed, filtering, maximum, remove, legend, y_filter, y_force):
    var_dict, source = read_files(files, bed, filtering, maximum, remove)
    if not var_dict:
        sys.exit('Empty variant dict, all variants were filtered')
    else:
        plots_af(var_dict, outfile, samples, bed, filtering, legend, y_filter, y_force)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plots the allele frequency of the variants in the input VCFs.')
    parser.add_argument('-i', '--infiles', type=str, nargs='+', help='Input VCFs; order will determine the order of labels in the plot')
    parser.add_argument('-o', '--outfile', type=str, help='PNG file to save the plot to')
    parser.add_argument('-s', '--samples', type=str, default=False, nargs='+', help='Names of the samples in same order as their files for the x labels (optional, will otherwise use file names)')
    parser.add_argument('-b', '--bed', type=str, default=False, help='BED file of the relevant positions (optional)')
    parser.add_argument('-f', '--filtering', action='store_true', default=False, help='If present, considers only entries passing a filter (optional)')
    parser.add_argument('-m', '--maximum', type=float, default=False, help='Max. AF for variants to be considered (optional)')
    parser.add_argument('-r', '--remove', action='store_true', default=False, help='Whether to remove variants that are not found in the first sample t0 (optional)')
    parser.add_argument('-y', '--y_filter', type=int, default=False, help='Percentile for filtering outliers (i.e. 80: Will find the highest 80 % percentile of all samples, then scale AF y axis to 110 % of that value (optional)')
    parser.add_argument('-yf', '--y_force', type=float, default=False, help='Maximal value for the main (AF) y axis to set to, will override y_filter (optional)')
    parser.add_argument('-l', '--legend', action='store_true', default=False, help='Whether to force the inclusion of a legend in the lineplots, otherwise only creates legend for <16 positions (optional)')

    args = parser.parse_args()
    if any(not f.endswith('.vcf') for f in args.infiles):
        sys.exit('Wrong file format for VCF, exiting...')
    if not args.outfile.endswith('.png'):
        sys.exit('Wrong file format for PNG, exiting...')
    if args.bed and not args.bed.endswith('.bed'):
        sys.exit('Wrong file format for BED, exiting...')
    if args.samples and len(args.samples) != len(args.infiles):
        sys.exit('Numbers of input files and sample names are not equal, exiting...')
    samples = args.samples or [os.path.splitext(os.path.split(f)[1])[0] for f in args.infiles]
    evaluate(args.infiles, samples, args.outfile, args.bed, args.filtering, args.maximum, args.remove, args.legend, args.y_filter, args.y_force)
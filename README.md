# ctSimu
## Snakemake pipeline for simulation of cancer patient data (FASTQ->BAM->VCF) with circulating tumor DNA (ctDNA)

This pipeline was written in snakemake, using own python scripts as well as NEAT: https://github.com/zstephens/neat-genreads
It simulates disease progress of patients during therapy for a number of different trends: Up (disease worsens over treatment duration), Down (disease burden lessens), Delayed Down (disease initially worsens, then condition improves), Relapse (disease initially gets better, then worse again), Constant (disease does not change strongly at all). After setting the number of patients to be simulated in the config file, this number of patients is simulated *per trend*, not overall. The current implementation simulates five sample timepoints per patient: before (before treatment start), therapy1, therapy2, therapy3 (these three samples are taken during treatment), after (follow-up sample after treatment ended). The names of these trends *are programmatically relevant*, as simulate_development.py only accepts any of these as valid trend argument; for data evaluation, they are not relevant, however, as the pipeline does not simulate according to a time frame but instead bases simulation of trends on a number of lists of values which give the approximate ctDNA variant allele frequency for each trend and timepoint. This means that new trends can be added and old ones changed in simulate_development.py by modifying the time_list and trend_dict (orders must be kept so that e.g. time_list[0] corresponds to trend_dict[some_trend][0]).

The pipeline  

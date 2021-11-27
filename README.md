# ctSimu
## Snakemake pipeline for simulation of cancer patient data (FASTQ->BAM->VCF) with circulating tumor DNA (ctDNA)

This pipeline was written in snakemake, using own python scripts as well as NEAT, which has to be installed: https://github.com/zstephens/neat-genreads
Additionally, if VarDict is to be used for variant calling, it has to be downloaded from https://github.com/AstraZeneca-NGS/VarDictJava (`git clone --recursive https://github.com/AstraZeneca-NGS/VarDictJava.git`) and the paths to its required scripts have to be configured in the config_ctSimu.yaml (VarDictJava/VarDict/vardict, VarDictJava/VarDict/teststrandbias.R, VarDictJava/VarDict/var2vcf_valid.pl).

The pipeline simulates disease progress of patients during therapy for a number of different trends: Up (disease worsens over treatment duration), Down (disease burden lessens), Delayed Down (disease initially worsens, then condition improves), Relapse (disease initially gets better, then worse again), Constant (disease does not change strongly at all). After setting the number of patients to be simulated in the config file, this number of patients is simulated *per trend*, not overall. All necessary packages can be installed from env_ctSimu.yaml: `conda env create -f env_ctSimu.yaml`

The current implementation simulates five sample timepoints per patient: before (before treatment start), therapy1, therapy2, therapy3 (these three samples are taken during treatment), after (follow-up sample after treatment ended). The names of these trends *are programmatically relevant*, as simulate_development.py only accepts any of these as valid trend argument; for data evaluation, they are not relevant, however, as the pipeline does not simulate according to a time frame but instead bases simulation of trends on a number of lists of values which give the approximate ctDNA variant allele frequency for each trend and timepoint. This means that new trends can be added and old ones changed in simulate_development.py by modifying the time_list and trend_dict (orders must be kept so that e.g. time_list[0] corresponds to trend_dict[some_trend][0]).

The pipeline accepts a number of parameters which can be changed in config_ctSimu.yaml:
-patients: int number of patients to simulate _per trend_
-trends: list of string of trends to simulate, currently any of ['up', 'down', 'const', 'rel', 'delay'] -->addition requires changing simulate_development.py
-timepoints: list of string of timepoints to simulate, currently all of ['before', 'therapy1', 'therapy2', 'therapy3', 'after'] -->needs to be compatible with simulate_development.py
-varcallers: list of string of variant callers to use, currently any of ['vardict', 'lofreq', 'umivar'] -->addition requires writing additional rules for the new variant callers
-template_vcf: VCF file listing the variants to be used for the simulation (do not need to be real cancer variants), the allele frequency is ignored as it is simulated (can be overwritten by adding F! (full, homozygous) or H! (half, heterozygous) to the third column) -->for umiVar, M has to be written in the third column so that all variants are monitored
-variation: float, which initial variation to use for the simulation of allele frequency; will increase in the actual datasets because of RNG function calls; a value of at least 0.02 is recommended

Additonal parameters can be used to modify the NEAT calls. All parameters are explained in the config file as well. A number of parameters are paths to reference genomes or external tools used by the pipeline

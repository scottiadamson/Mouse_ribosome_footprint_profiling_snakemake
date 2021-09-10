### These are two pipelines to support analysis of Ribosome footprint profiling data analysis from mice.

RFP_pipeline.snakemake is a pipeline analyzing data from [NT Ingolia et al., 2012](https://doi.org/10.1038/nprot.2012.086). 

RFP_pipeline_UMI.snakemake is a pipeline analyzing data from [NJ McGlincy & NT Igonlia, 2017](https://doi.org/10.1016/j.ymeth.2017.05.028).  
  * It should be noted that as of right now, this doesn't support multiplexing at the linker level, but does support UMI deduplication with the NI-810 (ATCGT) linker.  

If you are unfamiliar with snakemake, read about installation and documentation [here](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html). Be sure that snakemake, anaconda, and mamba are installed before running this pipeline.  

These also have commands for building genome and transcriptome indicies with hisat2 using mm10 and gencode vM25, and multiple steps rely on the gtf file, so replacing annotation versions may have impacts on multiple steps of the analysis.  


In order to run this you may also need to configure config.yaml files for your cluster setup. The one in this repository is intended for SLURM HPC systems.

This pipeline has a shared references folder for all project and specific project folders for each experiment. Create a shared foler and copy the scripts directory into that folder. 

```
git clone https://github.com/scottiadamson/Mouse_ribosome_footprint_profiling_snakemake.git
mkdir mouse_ribosome_profiling
cp -r scripts mouse_ribosome_profiling/
```

Then create a folder called raw_fastq to put the fastq files in.

```
mkdir mouse_ribosome_profiling/raw_fastq
``` 

Create a project configuration file, Gtpbp1_project.yaml and Gtpbp1_project_UMI.yaml are examples of these files. Properties include project name, sample names (make sure these sample names match the fastq file names i.e. mouse_ribosome_profiling/raw_fastq/sample1.fastq.gz), as well as read length minimum, maximium and reads per codon threshold. These can be altered depending on the periodicity for each read length which can be seen after running the riboWaltz portion of the pipeline, which may be experiment dependent. If you want to consider this, comment out the lines 11-14 for RFP_pipeline.snakemake or RFP_pipeline_UMI.snakemake, then run the pipeline up to the riboWaltz step and adjust the paramters in the project configuration file accordingly. After that you can uncomment lines 11-14 run the pipeline again (it will save the files from the steps you've already run).  

Ensure that the configfile on line 1 of either of the snakemake files corresponds to your project configuration file.  

Then run it:
```
snakemake --snakefile RFP_pipeline_UMI.snakemake --directory mouse_ribosome_profiling --use-conda --jobs 6 --profile ./ --nolock --latency-wait 60
```
or
```
snakemake --snakefile RFP_pipeline.snakemake --directory mouse_ribosome_profiling --use-conda --jobs 6 --profile ./ --nolock --latency-wait 60
```
Be sure config.yaml that has cluster parameters is in ./ and the --jobs corresponds to the number of jobs you want to run simultaneously. Allow about 24 hours to run the traditional pipeline in it's entirety and 32 for the UMI pipeline if --jobs = the number of samples in your experiment.  

#### License
This software is distributed under MIT License


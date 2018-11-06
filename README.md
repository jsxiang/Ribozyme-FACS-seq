# Ribozyme-FACS-seq

This repository contains the scripts used to analyze demultiplexed sequencing file in the form of a .csv file containing sequences in the first column, and read counts in subsequent columns from each indexed run, to generate results from FACS-seq experiments for measuring the fluorescent reporter activity of ribozymes and ribozyme switches in mammalian cells.

Example input file provided is MFS_combinedcounts.csv. Run analyze_MFS.m as a demo. This MatLab script calls facsseq.m, which is a class containing methods to create a facsseq object with computed parameters:

           ttestH: binary vector from unpaired t-test for significance in switching 
           ttestP: p-value from unpaired t-test
             seqs: cell array of m sequences
              mu1: 1xm vector of mean relative fluorescence, mu, for no ligand condition
              mu2: 1xm vector of mean relative fluorescence, mu, for with ligand condition
      adjustedmu2: mu2 values scaled to mu1 values using linear regression of spiked-in control ribozymes 
            ctrl1: mean relative fluorescence, mu, for no ligand condition for the control ribozymes
            ctrl2: mean relative fluorescence, mu, for with ligand condition for the control ribozymes
    adjustedctrl2: mu2 values scaled to mu1 values using linear regression of spiked-in control ribozymes, for the control ribozymes
             val1: mean relative fluorescence, mu, for no ligand condition for the validation candidates
             val2: mean relative fluorescence, mu, for with ligand condition for the validation candidates
     adjustedval2: mu2 values scaled to mu1 values using linear regression of spiked-in control ribozymes, for the validation candidates
      ctrlcounts1: raw read counts for control ribozymes, in the no ligand condition
      ctrlcounts2: raw read counts for control ribozymes, in the with ligand condition
       valcounts1: raw read counts for validation candidates, in the no ligand condition
       valcounts2: raw read counts for validation candidates, in the with ligand condition
         ctrlseqs: cell array of control ribozymes sequences
          valseqs: cell array of control validation candidates
          VYBmus1: 1xm vector of mean relative fluorescence, mu, for no ligand condition, scaled to flow cytometer values using linear regression of control ribozymes
      ctrlVYBmus1: mean relative fluorescence, mu, for no ligand condition for the control ribozymes, scaled to flow cytometer values
       valVYBmus1: mean relative fluorescence, mu, for no ligand condition for the validation candidates, scaled to flow cytometer values
          VYBmus2: 1xm vector of mean relative fluorescence, mu, for with ligand condition, scaled to flow cytometer values using linear regression of control ribozymes
      ctrlVYBmus2: mean relative fluorescence, mu, for with ligand condition for the control ribozymes, scaled to flow cytometer values
       valVYBmus2: mean relative fluorescence, mu, for with ligand condition for the validation candidates, scaled to flow cytometer values
            minus: structure that contain values for two replicates in the no ligand condition
             plus: structure that contain values for two replicates in the with ligand condition
             fold: fold change, or activation ratio, of the mean relative fluorescence 
            fold1: fold change, or activation ratio, of the mean relative fluorescence, for replicate 1
            fold2: fold change, or activation ratio, of the mean relative fluorescence, for replicate 2

To use, analyze_MFS.m needs to be modified to indicate which index columns correspond to which bins in each experiment conditions. Experiment conditions are replicates and different ligand concentrations, and need to be specified. Currently, only with and without ligand concentrations are allowed, and two replicate experiments are required. Validation candidate sequences need to added, control ribozyme sequences are defined in faccseq.m.

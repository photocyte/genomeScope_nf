### Dependencies
 - Nextflow
 - Conda
 - Singularity

### Usage

```
nextflow run photocyte/genomeScope_nf -resume -latest --reads './*paired.fq.gz' --executor pbs --cacheDir ~/scratch/nf_conda_envs
```

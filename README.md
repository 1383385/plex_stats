[![Build Status](https://travis-ci.com/tfursten/plex_stats.svg?branch=master)](https://travis-ci.com/tfursten/plex_stats)

# Plex Stats
This tool provides primer stats for multiplex PCR primers.

## Install

```
conda install --yes --strict-channel-priority --override-channels --channel conda-forge --channel bioconda --channel defaults --channel Tara_Furstenau plex_stats
```
Install in conda environment
```
conda create --yes --strict-channel-priority --override-channels --channel conda-forge --channel bioconda --channel defaults --channel Tara_Furstenau --name <ENV_NAME> plex_stats
```

## Input
The primers should be provided in FASTA format. IUPAC ambiguous DNA values are accepted and all of the possible primers will be generated before running the statistics (A suffix, like "_1", will be added to the degenerate primer name for each expanded version).

## Output
The statistics are output into a csv file, a styled excel file, and a styled html file. 
The table provides:
1. **Primer_Name**
2. **Seq**
3. **Primer_Length**
4. **Tm**: Melting Temperature estimated using NN calculation (salt parameters can be modified on the command line)
5. **Delta_Tm**: The absolute distance from the mean of all of the melting temperatures (values above 5 are highlighted).
6. **GC_Percent**: The GC content of the primer or percent that is made up of G's and C's. (values less than 40 and greater than 60 are highlighted).
7. **Length_Longest_Homopolymer**: Longest run of the same base over 3 bases (values greater than 5 are highlighted).
8. **Percent_Homopolymer**: The overall percent of the primer that contains runs of the same base over 3 bases (values greater than 35% are highlighted)
9. **Mean_Hybrid_Score**: The average score for all pairwise interactions between primers (values are colored as a heatmap).
10. **Mean_Hybrid_Run**: The average longest complimentary region for all pairwise interactions between primers (values are colored as a heatmap).
11. **Mean_3'_Score**: Same as Mean_Hybrid_Score except for just interactions that occur on the 3' end of the primer. This is important for primer dimers and when comparing the primer to itself, it could indicate a potential 3' hairpin.
12. **Mean_Percent_Hybrid**: The average percent of the primer that contains complimentary regions for all pairwise interactions (values are colored as a heatmap).
13. **Median_Hybrid_Score**: Median value for hybrid scores.
14. **Median_Hybrid_Run**: Median value for hybrid runs.
15. **Median_3'_Score**: Median value for 3' scores.
16. **Median_Percent_Hybrid**: Median value for hybrid percent.
13. **Max_Hybrid_Score**: Max value for hybrid scores.
14. **Max_Hybrid_Run**: Max value for hybrid runs.
15. **Max_3'_Score**: Max value for 3' scores.
16. **Max_Percent_Hybrid**: Max value for hybrid percent.

The remaining columns are the values for all of the pairwise interactions between the primers. Scores are highlighted if they are greater than or equal to 25, Runs are highlighted if they are greater than or equal to 6, and Percents are highlighted if they are greater than or equal to 50%. 

## Usage
```
$ plex-stats --help
Usage: plex-stats [OPTIONS] PRIMER_FASTA

Options:
  -o, --outpath DIRECTORY        Output directory.
  -p, --prefix TEXT              Output file name prefix
  -Na, --Na INTEGER RANGE        Millimolar concentration of Na
  -K, --K INTEGER RANGE          Millimolar concentration of K
  -Tris, --Tris INTEGER RANGE    Millimolar concentration of Tris
  -Mg, --Mg INTEGER RANGE        Millimolar concentration of Mg
  -dNTPs, --dNTPs INTEGER RANGE  Millimolar concentration of dNTPs
  --help                         Show this message and exit.
```

## Build
```
bash recipe/conda_build.sh
```










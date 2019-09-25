# TC_hunter

## TC-hunter identifies transgenic insertion sites within host genome

Takes a sorted and indexed bam file mapped to host genome and construct together as input, and output a plot and txt-file with possible breakpoints of insertion sites. 
The blue links correspond to reads with supplementary aligned reads in construct, and the gray links are soft clipped reads (one part of the read map to construct and the other to construct).  


## Install 

```
git clone https://github.com/vborjesson/TC_hunter.git
```

Install required programs and tools using Anaconda
```
conda env create --file TC_hunter/Scripts/Nextflow_env.txt
```

## Run 

```
nextflow TC_hunter.nf 
```

![](Plots/plot1.png)
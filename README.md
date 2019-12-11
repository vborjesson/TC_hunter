# TC_hunter

[](Plots/tc_hunter_out.png)

## TC-hunter identifies transgenic insertion sites within host genome

TC-hunter searches for transgenic insertion sites in a host genome and returns figures and a report to support these findings. 

Theres two programs; TC_hunter and TC_hunter_BWA. 

* TC_hunter.nf
Accepts one or several aligned BAM files (mapped to both host and transgenic sequence) as input. 
TC-hunter then identifies anchors and chimeric reads that maps to both host and transgenig sequence.    

* TC_hunter_BWA.nf 

TC_hunter_BWA accepts raw pair end fastq files (from one or several samples) as inbut and performes BWA MEM alignment before searching for trasgenic insertion site.       

## Software Dependencies

Nextflow v.19.01.0
BWA MEM v.0.7.5a (only if you run TC_hunter_BWA.nf)
Samtools v.1.9
R v.3.5.1
python v.3.5.1
igv v.2.1.7 (You can choose to run this separately)

## Install TC-hunter 

Clone the repository from Github and put it in your path (or add the direct path to config file) 
```
git clone https://github.com/vborjesson/TC_hunter.git

```

Install required programs and tools using Anaconda
```
conda env create --file TC_hunter/Scripts/Nextflow_env.txt
```

## Make Configuration file 

Create a configuration file from template.
```
cp TC_hunter/template/tc_hunter.config /path/to/WorkingDir 
```

Add required information to config file

| First Header  | Second Header |
| ------------- | ------------- |
| Content Cell  | Content Cell  |
| Content Cell  | Content Cell  |

### Run TC_hunter.nf

```
nextflow TC_hunter.nf -c <file.config> [-with-report <report name>]
```

### Run TC_hunter_BWA.nf

Before running, make sure you have a config file with all required information (see Config).  

```
nextflow TC_hunter_BWA.nf -c <file.config> [-with-report <report name>]
```

![](Plots/circlize.png)!


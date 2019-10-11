# How many (and which) resources should I request?

When you submit a job to the SLURM job scheduler on the Xanadu cluster, you must specify how many resources (nodes, cpus, memory) you need and which partition to submit to. But how to decide? This document will give new users some brief guidance. 

## What are the consequences of asking for the wrong resources? 

Why should you care? Why not just copy the same SLURM header from script to script and only worry about it when a job fails? 

Xanadu is a community resource. X users submit Y jobs per month. There can sometimes be a significant wait for jobs to begin running. Submitting jobs with inaccurate resource requests can tie things up unnecessarily and cause problems with job scheduling, wait times and time to completion.

- When you run a job and request too many resources:
	- When the cluster is busy, wait times are unnecessarily extended for everyone. If you only need 1 CPU but request 30, 29 CPUs will sit idly. 
	- You personally will wait longer for your job to run. To fairly allocate resources among users, job priority is lowered for jobs from users consuming many resources. 
- When you run a job and request too few resources: 
	- Your job may grind nearly to a halt and take a very long time to complete. 
	- Your job may fail, requiring you to resubmit it. 
	- The job scheduler may oversubscribe the node you are on, leading to the above problems. 


## How do I decide what to request? 

Ok, now that you're 100% convinced that for the benefit of yourself and the entire Xanadu community you need to try your best to accurately tailor your resource requests for each job, how are you supposed to do that? 

Well, the answer to this _can_ get complicated, but for most common use-cases in bioinformatics, it really isn't. 

To start with, here is a typical slurm header which specifies the details of a job submission:

```bash
#!/bin/bash 
#SBATCH --job-name=MY_JOB
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-user=MY_EMAIL@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
```

There are five lines in there that control the resources requested. You should consider their settings for **every job you submit**. 

```bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
```

### The number of tasks

`-n`: The number of _tasks_ you will run. A single number. These are tasks to be scheduled by SLURM and you invoke them using the SLURM command `srun`. They aren't lines in your script, or anything else like that. This option only applies if you use `srun` more than once within the script you submit using `sbatch`. If you don't use `srun` more than once, or at all, leave this at 1, or leave it out (the default is 1). 

### The number of nodes

`-N`: The number of _nodes_ you are requesting. Xanadu currently has 46 nodes, each with between 36 and 64 CPUs. You can think of nodes as individual computers that are all wired together to make up the cluster, and the CPUs as processing units that can operate independently, or be easily harnessed together for harder tasks. For _most cases_ you will want to specify this as 1. 

To explain: There are two ways that software packages are commonly parallelized. In the simplest, the software can use many CPUs within a single computer. This is how most 'multithreaded' software is written. If you are running something that has an option for how many 'threads' to use, it refers to this kind of parallelization. In the second, parallelization occurs on such a large scale that many nodes must coordinate their computation. This is not common in bioinformatics, and if you are a beginning user, you are unlikely to use software that does this. 

So, if you wish to use typical parallelized (a.k.a. multithreaded) software, you most likely want to set this number to 1. If you don't specify it, then the number of CPUs you request (see below) may be allocated across nodes. As an example, if you wanted to run the short read aligner `bwa`, using 12 cores (a.k.a. threads a.k.a. CPUs) to speed it up, your command might look like this:

```bash
bwa mem -t 12 reference.fa sample.1.fq sample.2.fq -o sample.sam
```
The `-t 12` indicates you want to use 12 cores. If you don't specify `-N 1`, then your job may be allocated 6 cores on node xanadu-04 and 6 cores on xanadu-05. Then `bwa` will only run on one of those nodes and only use 6 cores. The other 6 cores will sit idle because `bwa` is not written so that it can coordinate parallel computations across nodes. 

### The number of cores

`-c`: The number of _cores_ (a.k.a. threads or CPUs) you are requesting. This number is likely to change for every job you submit. If you wanted to run `bwa` with 12 cores as above, you would specify `-c 12`. If you fail to specify `-c 12`, then even though you tell `bwa` to use 12 cores, it won't be able to. If you specify `-c 12` but fail to give `bwa` the option `-t 12` then `bwa` will only use 1 core and 11 will sit idle. 

It is critically important that when using multithreaded software, you tell **both** SLURM and the software how many cores to use. 

How many cores should you request? If you're running a single process, as in this example, set `-c 12`. A general rule of thumb is that you should allocate the maximum number that your script will use at once. If you did something a little more complicated with `bwa`, piping the alignments to two other processes to mark duplicates and sort them, e.g.:

```bash
bwa mem -t 12 reference.fa sample.1.fq sample.2.fq | \
samblaster | \
samtools sort - >sample.sort.sam
```

then because these will all be running at the same time you would want to ask for 12 cores for `bwa` and 1 each for `samblaster` and `samtools`, for a total of `-c 14`. So essentially, you want to ask for the maximum number of cores that you'll be using all at once. 

### The quantity of memory

`--mem`: The quantity of _memory_ requested. The nodes on Xanadu have between 128-256gb for the general partition (see below) and 512gb-1tb for the himem partition. How much memory does your job need? This is a very tricky question, and it can depend on the particulars of not only the software you are using, but also your data. 

Hopefully, the software you are using can give you a sense of how much memory it will need. For example, a genome assembler may tell you that if you have X reads to assemble, then it will require Y memory. This is not always the case, however, and it may not be entirely predictable. 

In cases of high uncertainty, the main approach is to make your best guess and then start the job running. After it gets going, you can enter the compute node and use the command `top` to see how much memory (and CPU) your processes are using. You can do this by first running `squeue`, which will list all your currently running jobs, and the nodes they are running on. Then just do `ssh xanadu-XX` and then `top`. When you ssh into a compute node __you absolutely should not do anything else other than check the status of your running job__. All work done on the cluster needs to be done on resources allocated by SLURM. When you are done checking on your job, be sure to exit. 

There aren't hard limits on memory usage, so your process may start using more memory than you asked for. If it is dramatically more, you should cancel the job using `scancel` and resubmit it requesting a higher allocation. 


### The partition

`--partition`: The _partition_ you are requesting. The nodes on the Xanadu cluster are partitioned into sets with different characteristics. You can see these using `sinfo` which lists sets and their membership. You will notice these sets are overlapping. Most jobs will run on the "general" partition. Though you may specifically need xeon processors (e.g. for variant calling using GATK) or to use nodes with higher memory for memory-intensive tasks such as genome assembly. The most important thing here is not to submit to `himem` unless you're sure you need a very large amount of memory (> 256gb). 

### The quality of service

`--qos`: 

## After the job finishes: did you request the right resources?

After your job finishes running, it's a good idea to check and see if you specified an appropriate resource allocation. There are two ways to get this information. First, if you specify `#SBATCH --mail-user=MY_EMAIL@uconn.edu` in the SLURM header, then the following information will be e-mailed to you:

```bash
Job ID: YOUR_JOBID
Cluster: xanadu
User/Group: YOUR_USERNAME/domain users
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 03:15:37
CPU Efficiency: 99.96% of 03:15:42 core-walltime
Job Wall-clock time: 03:15:42
Memory Utilized: 1.45 GB
Memory Efficiency: 1.45% of 100.00 GB
```

Here you can see this user used their CPU allocation in a highly efficient manner, but requested **way** more memory than the job required, using only 1.45% of the requested memory. 

If you don't use the e-mail option in the SLURM header, you can get this exact same report using the command `seff YOUR_JOBID`. 


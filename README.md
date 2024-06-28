<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-eclipseq_logo_dark.png">
    <img alt="nf-core/eclipseq" src="docs/images/nf-core-eclipseq_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/eclipseq/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/eclipseq/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/eclipseq/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/eclipseq/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/eclipseq/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/eclipseq)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23eclipseq-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/eclipseq)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/eclipseq** is a bioinformatics pipeline that runs a version of [Clipper pipeline](https://www.encodeproject.org/documents/1f171ac6-a36a-41ac-b632-741aeb47aad2/@@download/attachment/eCLIP_analysisSOP_v2.3.pdf). The purpose of this pipeline is to detect RNA Binding protein binding sites.

![Alt text](eclipseq.drawio.png)

Workflow steps:
1.  Extract unique molecular barcodes using umi_tools (?do we need this?)
2.  Trim adapters using cutadapt (?? do we have the adaptors)
3.  QC with FastQC
4.  Alignment to genome (GRCh38 gencode v36) using STAR aligner
5.  Sorting and indexing using samtools
6.  Remove duplicates using umi_tools
7.  Call Peaks using Clipper0. Concatenate if necessary
1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Usage

If you are new to nextflow and nf-core, please create a mamba/conda environment as follows:

Next, make a copy of runEclipseq_template.sh to runEclipseq.sh and edit it appropriately.

Setup the genome: 

1. Download the genome you will use
2. Modify indexGenome.sh
3. Run indexGenome.sh
4. Index the genome with "samtools faidx genomeName.fa"
5. Calculate the chromosome sizes using "cut -f1,2 genomeName.fa.fai > sizes.genome"
6. Edit nextflow.config parameter to use this genome.

Then, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
ID,SAMPLE,REPLICATE,TYPE,FASTQ1,FASTQ2
CONTROL_REP1_SIGNAL,CONTROL,REP1,SIGNAL,s1_R1.fastq.gz,s1_R2.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

Now, you can run the pipeline using:

```bash
sbatch runEclipseq.sh samplesheet.csv outputDirectory
```
## Credits

nf-core/eclipseq was originally written by Ramiro Barrantes Reynolds with the help of Zach Miller.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#eclipseq` channel](https://nfcore.slack.com/channels/eclipseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/eclipseq for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->



An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

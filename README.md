# SNPPhylogenomics

**SNPPhylogenomics** is a reproducible, containerized WDL/Cromwell workflow for SNP-based bacterial phylogenomics. It performs read trimming, quality control, reference-based variant calling, core genome alignment generation, optional recombination filtering, maximum-likelihood phylogeny reconstruction, and tree visualization.

The workflow is designed for bacterial whole-genome sequencing data and is especially useful for genomic epidemiology, outbreak investigation, AMR surveillance, and comparative pathogen genomics.

---


## Workflow overview

<div align="center">
<pre>
Paired-end FASTQ files
⬇
Trimmomatic read trimming
⬇
FastQC quality control
⬇
Snippy variant calling
⬇
snippy-core core genome alignment
⬇
Optional Gubbins recombination filtering
⬇
IQ-TREE2 maximum-likelihood phylogeny
⬇
ETE3 tree visualization
</pre>
</div>

---
## Main outputs

The workflow produces:

Trimmed FASTQ files
FastQC reports
Per-sample Snippy variant calls
core.full.aln
core.aln
core.vcf
Optional Gubbins recombination-filtered alignment
IQ-TREE2 Newick tree
Tree visualization as PNG or SVG
Cleaned Newick tree for downstream tools such as iTOL
Requirements

Install the following:

Java
Cromwell
Docker

Example Cromwell command:

```bash
java -jar ~/cromwell-92.jar run SNPPhylogenomics.wdl --inputs ~/example_inputs.json

```
## Docker images used

| Step                      | Tool        | Docker image                  |
|---------------------------|-------------|-------------------------------|
| Trimming                  | Trimmomatic | staphb/trimmomatic:0.39       |
| QC                        | FastQC      | staphb/fastqc:0.11.9          |
| Variant calling           | Snippy      | staphb/snippy:4.6.0           |
| Recombination filtering   | Gubbins     | staphb/gubbins:3.4.1          |
| Phylogeny                 | IQ-TREE2    | staphb/iqtree2:2.3.4          |
| Tree visualization        | ETE3        | gmboowa/ete3-render:1.18      |

---
Input JSON structure

Example:

'example_inputs.json'

```json
{
  "SNPPhylogenomics.input_reads": [
    "~/sample1_1.fastq.gz",
    "~/sample1_2.fastq.gz",
    "~/sample2_1.fastq.gz",
    "~/sample2_2.fastq.gz"
  ],

  "SNPPhylogenomics.adapters": "~/adapters.fa",
  "SNPPhylogenomics.reference_genome": "~/reference.fasta",
  "SNPPhylogenomics.reference_type": "fasta",

  "SNPPhylogenomics.do_trimming": true,
  "SNPPhylogenomics.do_quality_control": true,
  "SNPPhylogenomics.do_variant_calling": true,
  "SNPPhylogenomics.do_phylogeny": true,

  "SNPPhylogenomics.use_gubbins": true,
  "SNPPhylogenomics.midpoint_root_tree": true,

  "SNPPhylogenomics.iqtree2_model": "GTR+G",
  "SNPPhylogenomics.iqtree2_bootstraps": 1000,

  "SNPPhylogenomics.max_cpus": 8,
  "SNPPhylogenomics.max_memory_gb": 16,

  "SNPPhylogenomics.min_read_length": 50,
  "SNPPhylogenomics.min_mapping_quality": 20,

  "SNPPhylogenomics.tree_width": 4000,
  "SNPPhylogenomics.tree_height": 3000,
  "SNPPhylogenomics.tree_image_format": "png"
}
````
## Important input rules

1. Paired reads must be in order

The workflow expects paired-end reads in strict R1/R2 order:

```bash
[
  "sample1_1.fastq.gz",
  "sample1_2.fastq.gz",
  "sample2_1.fastq.gz",
  "sample2_2.fastq.gz"
]
```
Do not mix the order.

2. Reference genome format

Use:
```bash
"SNPPhylogenomics.reference_type": "fasta"
```
for FASTA references. For phylogeny-only workflows, FASTA is often more stable and sufficient.

## IQ-TREE2 model options

The model is controlled here:

```bash
"SNPPhylogenomics.iqtree2_model": "GTR+G"
```

Common models
## IQ-TREE2 substitution models

| Model    | When to use                              | Notes                                      |
|----------|------------------------------------------|--------------------------------------------|
| GTR+G    | General bacterial SNP phylogeny           | Good default                               |
| GTR+I+G  | If some sites are invariant               | More complex; may improve fit              |
| GTR+F+G  | When empirical base frequencies are needed| Useful for biased base composition         |
| HKY+G    | Simpler datasets                         | Less parameter-rich than GTR               |
| JC       | Very simple testing only                 | Not recommended for final analysis         |
| MFP      | Let IQ-TREE choose best model            | Good for final analysis if runtime allows  |

## Recommended settings

For routine bacterial genomic epidemiology:

```bash
"SNPPhylogenomics.iqtree2_model": "GTR+G"
```
## For model testing:

```bash
"SNPPhylogenomics.iqtree2_model": "MFP"
```
Bootstrap settings

Use:
```bash

"SNPPhylogenomics.iqtree2_bootstraps": 1000
```
IQ-TREE2 ultrafast bootstrap requires at least 1000 replicates.

Recommended:

## Bootstrap settings

| Use case               | Bootstrap value |
|------------------------|-----------------|
| Quick testing          | 1000            |
| Standard publication   | 1000            |
| More rigorous analysis | 2000            |

Tree output format

The image format is controlled by:
```bash
"SNPPhylogenomics.tree_image_format": "png"
```
Supported options:

```bash

"png"
"svg"
"pdf"

```
Recommended

For reports and dashboards:
```bash
"SNPPhylogenomics.tree_image_format": "png"
```
For publication editing:
```bash
"SNPPhylogenomics.tree_image_format": "svg"
```
Recommended PNG resolution

For high-quality PNG output:

```bash

"SNPPhylogenomics.tree_width": 4000,
"SNPPhylogenomics.tree_height": 3000
```

For very high-quality output:

```bash
"SNPPhylogenomics.tree_width": 5000,
"SNPPhylogenomics.tree_height": 3500
```
Running the workflow

```bash
java -jar ~/cromwell-92.jar run SNPPhylogenomics.wdl --inputs ~/example_inputs.json
```
Finding outputs

After completion:

find cromwell-executions -name "core.full.aln"
find cromwell-executions -name "*.treefile"
find cromwell-executions -name "phylogenetic_tree.png"
find cromwell-executions -name "phylogenetic_tree.cleaned.nwk"

Key output files

| File | Description |
|------|-------------|
| core.full.aln | Full core genome alignment from Snippy-core |
| core.aln | Core SNP alignment |
| core.vcf | Combined core variant calls |
| gubbins.filtered_polymorphic_sites.fasta | Recombination-filtered alignment |
| final.treefile | Final Newick tree from IQ-TREE2 |
| phylogenetic_tree.png | Rendered tree image |
| phylogenetic_tree.cleaned.nwk | Cleaned tree used for visualization |

Set:
```bash
"SNPPhylogenomics.iqtree2_bootstraps": 1000
```
Gubbins error: AF_UNIX path too long

This can occur when Cromwell creates deeply nested directories. The workflow runs Gubbins from /tmp to avoid this issue.

Unexpected input error

If Cromwell reports:

Unexpected input provided

Check that every JSON key starts with the correct workflow name:
```bash
"SNPPhylogenomics."
```
not:

"rMAP."
JSON parsing error

Check for missing commas between FASTQ files. Every item must end with a comma except the last item.

Correct:
```bash
"sample1_2.fastq.gz",
"sample2_1.fastq.gz"
```
Incorrect:
```bash
"sample1_2.fastq.gz"
"sample2_1.fastq.gz"
```
Suggested citation text

**SNPPhylogenomics** is a reproducible WDL/Cromwell workflow for bacterial SNP-based phylogenomics. It integrates Trimmomatic, FastQC, Snippy, snippy-core, Gubbins, IQ-TREE2, and ETE3 to generate recombination-aware maximum-likelihood phylogenies from paired-end whole-genome sequencing data.

Recommended method description

Paired-end reads were trimmed using Trimmomatic and assessed using FastQC. Variants were called against a reference genome using Snippy, and a core genome alignment was generated using snippy-core. Recombinant regions were identified and filtered using Gubbins. A maximum-likelihood phylogeny was inferred using IQ-TREE2 under the selected nucleotide substitution model with ultrafast bootstrap support. The final tree was visualized using ETE3.

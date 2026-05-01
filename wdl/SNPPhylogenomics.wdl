version 1.0

workflow SNPPhylogenomics {
  input {
    Array[File]+ input_reads
    File adapters
    File reference_genome

    Boolean do_trimming = true
    Boolean do_quality_control = true
    Boolean do_variant_calling = true
    Boolean do_phylogeny = true
    Boolean use_gubbins = true
    Boolean midpoint_root_tree = true

    String reference_type = "genbank"
    String trimmomatic_quality_encoding = "phred33"
    String iqtree2_model = "GTR+G"
    Int iqtree2_bootstraps = 100
    Int max_cpus = 8
    Int max_memory_gb = 16

    Int min_read_length = 50
    Int min_mapping_quality = 20
    Int tree_width = 1800
    Int tree_height = 1400
    String tree_image_format = "png"
  }

  Int cpu_2 = if max_cpus < 2 then max_cpus else 2
  Int cpu_4 = if max_cpus < 4 then max_cpus else 4
  Int cpu_8 = if max_cpus < 8 then max_cpus else 8

  if (do_trimming) {
    call TRIMMING {
      input:
        input_reads = input_reads,
        adapters = adapters,
        trimmomatic_quality_encoding = trimmomatic_quality_encoding,
        cpu = cpu_4,
        min_length = min_read_length
    }
  }

  Array[File] analysis_reads = select_first([TRIMMING.trimmed_reads, input_reads])

  if (do_quality_control) {
    call QUALITY_CONTROL {
      input:
        input_reads = analysis_reads,
        cpu = cpu_4
    }
  }

  if (do_variant_calling) {
    call SNIPPY_CORE {
      input:
        input_reads = analysis_reads,
        reference_genome = reference_genome,
        reference_type = reference_type,
        cpu = cpu_8,
        memory_gb = max_memory_gb,
        min_quality = min_mapping_quality
    }
  }

  File core_full_alignment_for_phylogeny = select_first([SNIPPY_CORE.core_full_alignment])

  if (do_phylogeny && use_gubbins) {
    call GUBBINS_RECOMBINATION {
      input:
        core_full_alignment = core_full_alignment_for_phylogeny,
        cpu = cpu_8,
        memory_gb = max_memory_gb
    }
  }

  File phylogeny_alignment = select_first([GUBBINS_RECOMBINATION.filtered_alignment, core_full_alignment_for_phylogeny])

  if (do_phylogeny) {
    call IQTREE2_PHYLOGENY {
      input:
        alignment = phylogeny_alignment,
        model = iqtree2_model,
        bootstrap_replicates = iqtree2_bootstraps,
        cpu = cpu_8,
        memory_gb = max_memory_gb,
        midpoint_root_tree = midpoint_root_tree
    }

    call TREE_VISUALIZATION {
      input:
        input_tree = select_first([IQTREE2_PHYLOGENY.final_tree]),
        width = tree_width,
        height = tree_height,
        image_format = tree_image_format,
        title = "SNPPhylogenomics Snippy-core IQ-TREE2 phylogeny"
    }
  }

  output {
    Array[File]? trimmed_reads = TRIMMING.trimmed_reads
    Array[File]? fastqc_reports = QUALITY_CONTROL.fastqc_reports
    Array[File]? snippy_vcfs = SNIPPY_CORE.vcf_files
    File? variant_summary = SNIPPY_CORE.variant_summary
    File? core_full_alignment = SNIPPY_CORE.core_full_alignment
    File? core_snp_alignment = SNIPPY_CORE.core_snp_alignment
    File? snippy_core_tab = SNIPPY_CORE.core_tab
    File? gubbins_filtered_alignment = GUBBINS_RECOMBINATION.filtered_alignment
    File? gubbins_tree = GUBBINS_RECOMBINATION.gubbins_tree
    File? iqtree_newick = IQTREE2_PHYLOGENY.final_tree
    File? iqtree_log = IQTREE2_PHYLOGENY.iqtree_log
    File? tree_image = TREE_VISUALIZATION.tree_image
    File? tree_render_log = TREE_VISUALIZATION.render_log
  }
}

task TRIMMING {
  input {
    Array[File]+ input_reads
    File adapters
    String trimmomatic_quality_encoding = "phred33"
    Int cpu = 4
    Int min_length = 50
  }

  command <<<
    set -euo pipefail
    mkdir -p trimmed

    files=(~{sep=' ' input_reads})
    n=${#files[@]}

    if [ $((n % 2)) -ne 0 ]; then
      echo "ERROR: input_reads must contain paired reads in R1/R2 order." >&2
      exit 1
    fi

    for ((i=0; i<n; i+=2)); do
      R1="${files[$i]}"
      R2="${files[$((i+1))]}"

      sample=$(basename "$R1")
      sample=$(echo "$sample" | sed -E 's/(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$//' | sed -E 's/(_R?1|_1|\.R?1|\.1)(_|$).*//')

      trimmomatic PE \
        -threads ~{cpu} \
        -~{trimmomatic_quality_encoding} \
        "$R1" "$R2" \
        "trimmed/${sample}_R1_paired.fastq.gz" "trimmed/${sample}_R1_unpaired.fastq.gz" \
        "trimmed/${sample}_R2_paired.fastq.gz" "trimmed/${sample}_R2_unpaired.fastq.gz" \
        ILLUMINACLIP:~{adapters}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:~{min_length}
    done
  >>>

  runtime {
    docker: "staphb/trimmomatic:0.39"
    cpu: cpu
    memory: "8 GB"
    disks: "local-disk 100 HDD"
  }

  output {
    Array[File] trimmed_reads = glob("trimmed/*_paired.fastq.gz")
    Array[File] trimming_logs = glob("*.log")
  }
}

task QUALITY_CONTROL {
  input {
    Array[File]+ input_reads
    Int cpu = 4
  }

  command <<<
    set -euo pipefail
    mkdir -p fastqc
    fastqc -t ~{cpu} -o fastqc ~{sep=' ' input_reads}
  >>>

  runtime {
    docker: "staphb/fastqc:0.11.9"
    cpu: cpu
    memory: "8 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    Array[File] fastqc_reports = glob("fastqc/*_fastqc.html")
    Array[File] fastqc_zips = glob("fastqc/*_fastqc.zip")
  }
}

task SNIPPY_CORE {
  input {
    Array[File]+ input_reads
    File reference_genome
    String reference_type = "genbank"
    Int cpu = 8
    Int memory_gb = 16
    Int min_quality = 20
  }

  command <<<
    set -euo pipefail
    mkdir -p snippy_results snippy_core

    files=(~{sep=' ' input_reads})
    n=${#files[@]}

    if [ $((n % 2)) -ne 0 ]; then
      echo "ERROR: input_reads must contain paired reads in R1/R2 order." >&2
      exit 1
    fi

    echo -e "sample\tstatus\tvcf\taligned_fasta" > variant_summary.tsv

    for ((i=0; i<n; i+=2)); do
      R1="${files[$i]}"
      R2="${files[$((i+1))]}"

      sample=$(basename "$R1")
      sample=$(echo "$sample" | sed -E 's/(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$//' | sed -E 's/(_R?1_paired|_R?1|_1_paired|_1|\.R?1|\.1)(_|$).*//')
      outdir="snippy_results/${sample}"

      snippy \
        --cpus ~{cpu} \
        --minqual ~{min_quality} \
        --ref "~{reference_genome}" \
        --R1 "$R1" \
        --R2 "$R2" \
        --outdir "$outdir" \
        --prefix "$sample" \
        --force

      echo -e "${sample}\tsuccess\t${outdir}/${sample}.vcf\t${outdir}/${sample}.aligned.fa" >> variant_summary.tsv
    done

    snippy-core \
      --ref "~{reference_genome}" \
      --prefix snippy_core/core \
      snippy_results/*

    cat > variant_summary.html <<'HTML'
<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>Snippy Variant Calling Summary</title>
<style>
body{font-family:Arial,Helvetica,sans-serif;margin:24px}
table{border-collapse:collapse;width:100%}
th,td{border:1px solid #ddd;padding:8px;text-align:left}
th{background:#0ea5e9;color:white}
</style>
</head>
<body>
<h1>Snippy Variant Calling Summary</h1>
<table>
<thead><tr><th>Sample</th><th>Status</th><th>VCF</th><th>Aligned FASTA</th></tr></thead>
<tbody>
HTML

    awk -F '\t' 'NR>1 {printf "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\n",$1,$2,$3,$4}' variant_summary.tsv >> variant_summary.html

    cat >> variant_summary.html <<'HTML'
</tbody>
</table>
</body>
</html>
HTML
  >>>

  runtime {
    docker: "staphb/snippy:4.6.0"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk 200 HDD"
    timeout: "48 hours"
  }

  output {
    Array[File] vcf_files = glob("snippy_results/*/*.vcf")
    Array[File] aligned_fastas = glob("snippy_results/*/*.aligned.fa")
    File variant_summary = "variant_summary.html"
    File core_full_alignment = "snippy_core/core.full.aln"
    File core_snp_alignment = "snippy_core/core.aln"
    File core_tab = "snippy_core/core.tab"
    File core_vcf = "snippy_core/core.vcf"
  }
}

task GUBBINS_RECOMBINATION {
  input {
    File core_full_alignment
    Int cpu = 8
    Int memory_gb = 16
  }

  command <<<
    set -euo pipefail

    mkdir -p gubbins

    # Avoid Python multiprocessing AF_UNIX socket path-length errors
    # caused by Cromwell's deeply nested execution directories.
    export TMPDIR=/tmp
    export TMP=/tmp
    export TEMP=/tmp

    workdir="/tmp/gubbins_work_${RANDOM}_${RANDOM}"
    mkdir -p "$workdir"

    cp "~{core_full_alignment}" "$workdir/core.full.aln"

    cd "$workdir"

    run_gubbins.py \
      --threads ~{cpu} \
      --prefix gubbins \
      core.full.aln

    cp gubbins.filtered_polymorphic_sites.fasta "$OLDPWD/gubbins/"
    cp gubbins.final_tree.tre "$OLDPWD/gubbins/"
    cp gubbins.recombination_predictions.gff "$OLDPWD/gubbins/"

    cd "$OLDPWD"
    rm -rf "$workdir"
  >>>

  runtime {
    docker: "staphb/gubbins:3.4.1"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk 200 HDD"
    timeout: "48 hours"
  }

  output {
    File filtered_alignment = "gubbins/gubbins.filtered_polymorphic_sites.fasta"
    File gubbins_tree = "gubbins/gubbins.final_tree.tre"
    File recombination_predictions = "gubbins/gubbins.recombination_predictions.gff"
  }
}
task IQTREE2_PHYLOGENY {
  input {
    File alignment
    String model = "GTR+G"
    Int bootstrap_replicates = 100
    Int cpu = 8
    Int memory_gb = 16
    Boolean midpoint_root_tree = true
  }

  command <<<
    set -euo pipefail
    mkdir -p iqtree

    cp "~{alignment}" iqtree/phylogeny_alignment.fasta

    iqtree2 \
      -s iqtree/phylogeny_alignment.fasta \
      -m ~{model} \
      -bb ~{bootstrap_replicates} \
      -nt ~{cpu} \
      -pre iqtree/SNPPhylogenomics_phylogeny

    # The staphb/iqtree2 container does not include python3.
    # Keep the IQ-TREE2 output tree as final.treefile for downstream visualization.
    cp iqtree/SNPPhylogenomics_phylogeny.treefile iqtree/final.treefile
  >>>

  runtime {
    docker: "staphb/iqtree2:2.3.4"
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk 100 HDD"
    timeout: "48 hours"
  }

  output {
    File final_tree = "iqtree/final.treefile"
    File raw_tree = "iqtree/SNPPhylogenomics_phylogeny.treefile"
    File iqtree_log = "iqtree/SNPPhylogenomics_phylogeny.log"
    File iqtree_report = "iqtree/SNPPhylogenomics_phylogeny.iqtree"
  }
}
task TREE_VISUALIZATION {
  input {
    File input_tree
    Int width = 1800
    Int height = 1400
    String image_format = "png"
    String title = "SNPPhylogenomics."
  }

  command <<<
    set -euo pipefail
    mkdir -p tree_visualization
    export QT_QPA_PLATFORM=offscreen
    export MPLBACKEND=Agg

    python3 - <<'PY'
from pathlib import Path
import sys

tree_path = Path("~{input_tree}")
out_img = Path("tree_visualization/phylogenetic_tree.~{image_format}")
out_tree = Path("tree_visualization/phylogenetic_tree.cleaned.nwk")
log = Path("tree_visualization/render.log")

try:
    from ete3 import Tree, TreeStyle, NodeStyle

    t = Tree(str(tree_path), format=1)

    # 1. Remove Reference from final display only
    for leaf in t.get_leaves():
        if leaf.name.lower() == "reference":
            leaf.detach()

    # 2. Convert bootstrap labels from proportions to percentages if needed
    for node in t.traverse():
        if not node.is_leaf():
            try:
                support = float(node.support)
                if 0 < support <= 1:
                    node.support = round(support * 100)
                else:
                    node.support = round(support)
            except Exception:
                pass

    # Save cleaned tree used for figure
    t.write(outfile=str(out_tree), format=1)

    # 3. Style tree
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_support = True

    # 4. Hide confusing large scale bar
    ts.show_scale = False

    # 5. No title added to the rendered image

    ts.margin_left = 20
    ts.margin_right = 20
    ts.margin_top = 20
    ts.margin_bottom = 20

    # Smaller, cleaner nodes
    for node in t.traverse():
        ns = NodeStyle()
        ns["size"] = 4
        node.set_style(ns)

    t.render(str(out_img), w=~{width}, units="px", tree_style=ts)

    log.write_text(
        "Rendered cleaned phylogenetic tree\\n"
        "Reference removed from display\\n"
        "Bootstrap values converted from proportions to percentages when needed\\n"
        "Scale bar hidden\\n"
        "Title intentionally omitted from image\\n"
    )

except Exception as e:
    log.write_text(f"ETE3 rendering failed: {e}\\n")
    try:
        from PIL import Image, ImageDraw
        img = Image.new("RGB", (~{width}, ~{height}), "white")
        d = ImageDraw.Draw(img)
        d.text((30, 30), "Tree rendering failed. Newick tree file was generated successfully.", fill="black")
        d.text((30, 80), str(e), fill="red")
        img.save(out_img)
    except Exception as e2:
        sys.stderr.write(f"Failed to create fallback image: {e2}\\n")
        raise
PY
  >>>

  runtime {
    docker: "gmboowa/ete3-render:1.18"
    cpu: 2
    memory: "8 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File tree_image = "tree_visualization/phylogenetic_tree.~{image_format}"
    File cleaned_tree = "tree_visualization/phylogenetic_tree.cleaned.nwk"
    File render_log = "tree_visualization/render.log"
  }
}

version 1.0 

#IMPORTS
## According to this: https://cromwell.readthedocs.io/en/stable/Imports/ we can import raw from github
## so we can make use of the already written WDLs provided by WARP/VUMC Biostatistics

import "https://raw.githubusercontent.com/shengqh/warp/develop/tasks/vumc_biostatistics/GcpUtils.wdl" as http_GcpUtils
import "https://raw.githubusercontent.com/shengqh/warp/develop/tasks/vumc_biostatistics/Plink2Utils.wdl" as http_Plink2Utils

# WORKFLOW

workflow ancestry_pipeline_vcf_scope_unsupervised {
    input{
        # required inputs original data as array of chromosomes 
        File input_vcf

        # optional inputs for spike in data - required if merging spike in data 
        Boolean external_spike_in = true 
        File? spike_in_pgen
        File? spike_in_pvar
        File? spike_in_psam
        File? spike_in_relatives_exclude 

        String output_prefix

        # optional outputs for exporting 
        String? target_gcp_folder

        #optional inputs for unsupervised scope - required if running unsupervised scope, only use bi-allelic variants
        String scope_plink2_maf_filter = "--maf 0.001 --snps-only --const-fid --max-alleles 2 --set-all-var-ids chr@:#:\\$r:\\$a --allow-extra-chr --chr 1-22"

        String scope_plink2_LD_filter_option = "--indep-pairwise 50000 80 0.1"
        File scope_long_range_ld_file
        Int K = 5
        Int seed = 20260817
    }

    String TargetPgen_output_prefix = output_prefix + ".target"
    call VcfToPgenFilterByMAF as TargetPgen {
        input:
            vcf_file = input_vcf,
            plink2_maf_filter = scope_plink2_maf_filter,
            exclude_range_file = scope_long_range_ld_file,
            output_prefix = TargetPgen_output_prefix
    }

    if (external_spike_in) {
        call http_Plink2Utils.PgenFilter as SpikeinPgen {
            input:
                input_pgen = select_first([spike_in_pgen]),
                input_pvar = select_first([spike_in_pvar]),
                input_psam = select_first([spike_in_psam]),
                exclude_sample_id_file = spike_in_relatives_exclude,
                exclude_range_file = scope_long_range_ld_file,
                plink2_filter_option = scope_plink2_maf_filter,
                output_prefix = output_prefix + ".spikein"
        }

        String MergePgen_output_prefix = output_prefix + ".target_spikein_intersect"
        call MergeSpikeinAndTargetPgen as MergePgen {
            input: 
                target_pgen = TargetPgen.output_pgen,
                target_pvar = TargetPgen.output_pvar,
                target_psam = TargetPgen.output_psam,
                spikein_pgen = SpikeinPgen.output_pgen,
                spikein_pvar = SpikeinPgen.output_pvar,
                spikein_psam = SpikeinPgen.output_psam,
                output_prefix = MergePgen_output_prefix
        }
    }

    String PreparePlinkBed_output_prefix = select_first([MergePgen_output_prefix, TargetPgen_output_prefix]) + ".prepared"
    call PreparePlinkBed {
        input:
            input_pgen = select_first([MergePgen.output_pgen, TargetPgen.output_pgen]),
            input_pvar = select_first([MergePgen.output_pvar, TargetPgen.output_pvar]),
            input_psam = select_first([MergePgen.output_psam, TargetPgen.output_psam]),
            output_prefix = PreparePlinkBed_output_prefix,
            plink2_LD_filter_option = scope_plink2_LD_filter_option
    }

    call RunScopeUnsupervised {    
        input:
            bed_file = PreparePlinkBed.output_bed,
            bim_file = PreparePlinkBed.output_bim,
            fam_file = PreparePlinkBed.output_fam,
            K = K,
            output_string = PreparePlinkBed_output_prefix,
            seed = seed
    }

    if(defined(target_gcp_folder)){
        call http_GcpUtils.MoveOrCopyThreeFiles as CopyFiles{
            input:
                source_file1 = RunScopeUnsupervised.outP,
                source_file2 = RunScopeUnsupervised.outQ,
                source_file3 = RunScopeUnsupervised.outV,
                is_move_file = false,
                target_gcp_folder = select_first([target_gcp_folder])
        }
    }

    output {
        File outP = select_first([CopyFiles.output_file1, RunScopeUnsupervised.outP])   
        File outQ = select_first([CopyFiles.output_file2, RunScopeUnsupervised.outQ])
        File outV = select_first([CopyFiles.output_file3, RunScopeUnsupervised.outV])
    }
}

task VcfToPgenFilterByMAF {
    input {
        File vcf_file

        File? exclude_range_file

        String output_prefix

        String plink2_maf_filter

        Int memory_gb = 20

        String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
    }

    Int disk_size = ceil(size([vcf_file], "GB")  * 2) + 20

    command <<<
  
plink2 \
    --vcf ~{vcf_file} \
    ~{plink2_maf_filter} \
    ~{"--exclude range " + exclude_range_file} \
    --make-pgen \
    --out ~{output_prefix}

>>>

    runtime {
        docker: docker
        preemptible: 1
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_gb + " GiB"
    }

    output {
        File output_pgen = "~{output_prefix}.pgen"
        File output_pvar = "~{output_prefix}.pvar"
        File output_psam = "~{output_prefix}.psam"
    }
}

task MergeSpikeinAndTargetPgen{
    input {
        File target_pgen
        File target_pvar
        File target_psam

        File spikein_pgen
        File spikein_pvar
        File spikein_psam

        String output_prefix

        Int? memory_gb = 20

        String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
    }

    Int disk_size = ceil(size([target_pgen, target_pvar, target_psam, spikein_pgen, spikein_pvar, spikein_psam], "GB")  * 3) + 10

    String target_prefix = basename(target_pgen, ".pgen")
    String spikein_prefix = basename(spikein_pgen, ".pgen")

    String target_snp_list = target_prefix + ".snplist"
    String spikein_snp_list = spikein_prefix + ".snplist"

    String target_prefix_common = target_prefix + "_common"
    String spikein_prefix_common = spikein_prefix + "_common"

    runtime {
        docker: docker
        preemptible: 1
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_gb + " GiB"
    }

    command <<<

    # extract target snp list
plink2 \
    --pgen ~{target_pgen} \
    --pvar ~{target_pvar} \
    --psam ~{target_psam} \
    --rm-dup force-first 'list' \
    --write-snplist \
    --out ~{target_prefix} 

# extract spikein snp list
plink2 \
    --pgen ~{spikein_pgen} \
    --pvar ~{spikein_pvar} \
    --psam ~{spikein_psam} \
    --write-snplist \
    --out ~{spikein_prefix}

# prepare target
plink2 \
    --pgen ~{target_pgen} \
    --pvar ~{target_pvar} \
    --psam ~{target_psam} \
    --extract-intersect ~{target_snp_list} ~{spikein_snp_list} \
    --make-bed \
    --out ~{target_prefix_common}

# prepare spikein
plink2 \
    --pgen ~{spikein_pgen} \
    --pvar ~{spikein_pvar} \
    --psam ~{spikein_psam} \
    --extract-intersect ~{target_snp_list} ~{spikein_snp_list} \
    --make-bed \
    --out ~{spikein_prefix_common}

# merge
plink \
    --bfile ~{target_prefix_common} \
    --bmerge ~{spikein_prefix_common} \
    --make-bed \
    --out ~{output_prefix}

# convert to pgen
plink2 \
    --bfile ~{output_prefix} \
    --make-pgen \
    --out ~{output_prefix}

rm -f ~{target_prefix_common}.bed ~{target_prefix_common}.bim ~{target_prefix_common}.fam
rm -f ~{spikein_prefix_common}.bed ~{spikein_prefix_common}.bim ~{spikein_prefix_common}.fam
rm -f ~{target_prefix}.snplist ~{spikein_prefix}.snplist
rm -f ~{output_prefix}.bed ~{output_prefix}.bim ~{output_prefix}.fam

>>>

    output {
        File output_pgen = output_prefix + ".pgen"
        File output_pvar = output_prefix + ".pvar"
        File output_psam = output_prefix + ".psam"
    }
}

task PreparePlinkBed {
    input {
        File input_pgen
        File input_pvar
        File input_psam

        String output_prefix

        String plink2_LD_filter_option = "--indep-pairwise 50000 80 0.1"

        Int memory_gb = 20
        Float disk_size_multiplier = 4
        Int disk_size_addition = 20
        String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
    }

    Int disk_size = ceil(size([input_pgen, input_pvar, input_psam], "GB")  * disk_size_multiplier) + disk_size_addition

    String input_prefix = basename(input_pgen, ".pgen")

    command <<<

plink2 \
    --pgen ~{input_pgen} \
    --pvar ~{input_pvar} \
    --psam ~{input_psam} \
    ~{plink2_LD_filter_option}

plink2 \
    --pgen ~{input_pgen} \
    --pvar ~{input_pvar} \
    --psam ~{input_psam} \
    --extract plink2.prune.in \
    --make-bed \
    --out ~{output_prefix}

>>>

    runtime {
        docker: docker
        preemptible: 1
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_gb + " GiB"
    }

    output {
        File output_bed = "~{output_prefix}.bed"
        File output_bim = "~{output_prefix}.bim"
        File output_fam = "~{output_prefix}.fam"
    }
}

task RunScopeUnsupervised {
    input {
        File bed_file
        File bim_file
        File fam_file

        Int K
        String output_string
        Int seed = 20260813

        Int memory_gb = 60
        String docker = "blosteinf/scope:0.1"
    }

    String plink_binary_prefix = basename(bed_file, ".bed")
    String relocated_target_bed = plink_binary_prefix + ".bed"
    String relocated_target_bim = plink_binary_prefix + ".bim"
    String relocated_target_fam = plink_binary_prefix + ".fam"

    String unsup_output = output_string + "_unsupervised_" 
    Int disk_size = ceil(size([bed_file, bim_file, fam_file], "GB")  * 2) + 20

    command <<<

ln -s ~{bed_file} ./~{relocated_target_bed}
ln -s ~{bim_file} ./~{relocated_target_bim}
ln -s ~{fam_file} ./~{relocated_target_fam}

scope -g ~{plink_binary_prefix} -k ~{K} -seed ~{seed} -o ~{unsup_output}
awk '{ for (i=1; i<=NF; i++) { a[NR,i] = $i } } NF>p { p = NF } END { for(j=1; j<=p; j++) { str=a[1,j]; for(i=2; i<=NR; i++) { str=str" "a[i,j]; } print str } }' ~{unsup_output}Qhat.txt > transposed_Qhat.txt
cut -f2 ./~{relocated_target_fam} | paste - transposed_Qhat.txt > ~{unsup_output}Qhat.txt

rm -f ./~{relocated_target_bed} ./~{relocated_target_bim} ./~{relocated_target_fam}

    >>>

    runtime {
        docker: docker
        preemptible: 1
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_gb + " GiB"
    }

    output {
        File outP= "${unsup_output}Phat.txt"
        File outQ= "${unsup_output}Qhat.txt"
        File outV= "${unsup_output}V.txt"
    }
}


version 1.0 

#IMPORTS
## According to this: https://cromwell.readthedocs.io/en/stable/Imports/ we can import raw from github
## so we can make use of the already written WDLs provided by WARP/VUMC Biostatistics

import "https://raw.githubusercontent.com/shengqh/warp/develop/tasks/vumc_biostatistics/Plink2Utils.wdl" as http_Plink2Utils
import "https://raw.githubusercontent.com/shengqh/warp/develop/tasks/vumc_biostatistics/GcpUtils.wdl" as http_GcpUtils

# WORKFLOW

workflow ancestry_pipeline_vcf_scope_unsupervised {
    input{
        # required inputs original data as array of chromosomes 
        File input_vcf

        File? g1000_pgen
        File? g1000_pvar
        File? g1000_psam

        String output_prefix

        # optional outputs for exporting 
        String? target_gcp_folder

        #optional inputs for unsupervised scope - required if running unsupervised scope
        String scope_plink2_maf_filter = "--maf 0.01"

        String scope_plink2_LD_filter_option = "--indep-pairwise 50000 80 0.1"
        File scope_long_range_ld_file
        Int K = 4
        Int seed = 1234
    }

    call VcfToPgenFilterByMAF {
        input:
            vcf_file = input_vcf,
            plink2_maf_filter = scope_plink2_maf_filter,
            output_prefix = output_prefix + ".target"
    }

    if defined(g1000_pgen) {
        call http_Plink2Utils.PgenFilter as PgenFilterByMAF{
            input:
                input_pgen = g1000_pgen,
                input_pvar = g1000_pvar,
                input_psam = g1000_psam,
                output_prefix = output_prefix + ".1000g",
                plink2_filter_option = scope_plink2_maf_filter
        }

        call http_Plink2Utils.MergePgenFiles {
            input:
                input_pgen_files = [PgenFilterByMAF.output_pgen, VcfToPgenFilterByMAF.output_pgen],
                input_pvar_files = [PgenFilterByMAF.output_pvar, VcfToPgenFilterByMAF.output_pvar],
                input_psam_files = [PgenFilterByMAF.output_psam, VcfToPgenFilterByMAF.output_psam],
                output_prefix = output_prefix + ".merged"
        }
    }

    call PreparePlinkBed {
        input:
            input_pgen = select_first([MergePgenFiles.output_pgen, VcfToPgenFilterByMAF.output_pgen]),
            input_pvar = select_first([MergePgenFiles.output_pvar, VcfToPgenFilterByMAF.output_pvar]),
            input_psam = select_first([MergePgenFiles.output_psam, VcfToPgenFilterByMAF.output_psam]),
            output_prefix = output_prefix + ".prepared",
            plink2_LD_filter_option = scope_plink2_LD_filter_option,
            long_range_ld_file = scope_long_range_ld_file
    }

    call RunScopeUnsupervised{    
        input:
            bed_file = PreparePlinkBed.output_bed,
            bim_file = PreparePlinkBed.output_bim,
            fam_file = PreparePlinkBed.output_fam,
            K = K,
            output_string = output_prefix,
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

# TASKS

## for scope
task VcfToPgenFilterByMAF {
  input {
    File vcf_file

    String output_prefix

    String plink2_maf_filter = "--maf 0.01"

    Int memory_gb = 20

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  Int disk_size = ceil(size([vcf_file], "GB")  * 2) + 20

  command <<<
  
plink2 \
    --vcf ~{vcf_file} \
    ~{plink2_maf_filter} \
    --snps-only \
    --const-fid \
    --max-alleles 2 \
    --set-all-var-ids chr@:#:\$r:\$a \
    --new-id-max-allele-len 1000 \
    --make-pgen \
    --allow-extra-chr \
    --chr 1-22 \
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

## for scope
task PreparePlinkBed{
  input {
    File input_pgen
    File input_pvar
    File input_psam

    String output_prefix

    String plink2_LD_filter_option = "--indep-pairwise 50000 80 0.1"
    File long_range_ld_file

    Int memory_gb = 20
    Float disk_size_multiplier = 4
    Int disk_size_addition = 20
    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  Int disk_size = ceil(size([input_pgen, input_pvar, input_psam], "GB")  * disk_size_multiplier) + disk_size_addition

  command <<<
        
plink2 \
    --pgen ~{input_pgen} \
    --pvar ~{input_pvar} \
    --psam ~{input_psam} \
    --exclude range ~{long_range_ld_file} \
    --make-pgen \
    --out maf_filtered_longrange

plink2 \
    --pgen maf_filtered_longrange.pgen \
    --pvar maf_filtered_longrange.pvar \
    --psam maf_filtered_longrange.psam \
    ~{plink2_LD_filter_option}

plink2 \
    --pgen maf_filtered_longrange.pgen \
    --pvar maf_filtered_longrange.pvar \
    --psam maf_filtered_longrange.psam \
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

task RunScopeUnsupervised{
    input{
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
    String relocated_bed = plink_binary_prefix + ".bed"
    String relocated_bim = plink_binary_prefix + ".bim"
    String relocated_fam = plink_binary_prefix + ".fam"

    String unsup_output = output_string + "_unsupervised_" 
    Int disk_size = ceil(size([bed_file, bim_file, fam_file], "GB")  * 2) + 20

    command <<<

ln -s ~{bed_file} ./~{relocated_bed}
ln -s ~{bim_file} ./~{relocated_bim}
ln -s ~{fam_file} ./~{relocated_fam}

scope -g ~{plink_binary_prefix} -k ~{K} -seed ~{seed} -o ~{unsup_output}
awk '{ for (i=1; i<=NF; i++) { a[NR,i] = $i } } NF>p { p = NF } END { for(j=1; j<=p; j++) { str=a[1,j]; for(i=2; i<=NR; i++) { str=str" "a[i,j]; } print str } }' ~{unsup_output}Qhat.txt > transposed_Qhat.txt
cut -f2 ./~{relocated_fam} | paste - transposed_Qhat.txt > ~{unsup_output}Qhat.txt

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

version 1.0 

#IMPORTS
## According to this: https://cromwell.readthedocs.io/en/stable/Imports/ we can import raw from github
## so we can make use of the already written WDLs provided by WARP/VUMC Biostatistics

import "https://raw.githubusercontent.com/shengqh/warp/develop/tasks/vumc_biostatistics/GcpUtils.wdl" as http_GcpUtils
import "https://raw.githubusercontent.com/shengqh/warp/develop/pipelines/vumc_biostatistics/genotype/Utils.wdl" as http_GenotypeUtils
import "https://raw.githubusercontent.com/shengqh/warp/develop/pipelines/vumc_biostatistics/agd/AgdUtils.wdl" as http_AgdUtils


# WORKFLOW

workflow ancestry_pipeline_vcf_scope_unsupervised {
    input{
        # required inputs original data as array of chromosomes 
        File input_vcf

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

    call PreparePlink as PreparePlink{
        input:
            vcf_file = input_vcf,
            long_range_ld_file = scope_long_range_ld_file,
            plink2_maf_filter = scope_plink2_maf_filter,
            plink2_LD_filter_option = scope_plink2_LD_filter_option,
            output_prefix = output_prefix
    }

    call ConvertPgenToBed as ConvertPgenToBedForScope{
        input: 
            pgen = PreparePlink.output_pgen, 
            pvar = PreparePlink.output_pvar,
            psam = PreparePlink.output_psam, 
    }

    call RunScopeUnsupervised{    
        input:
            bed_file = ConvertPgenToBedForScope.convert_Pgen_out_bed,
            bim_file = ConvertPgenToBedForScope.convert_Pgen_out_bim,
            fam_file = ConvertPgenToBedForScope.convert_Pgen_out_fam,
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
task ConvertPgenToBed{
    input {
        File pgen 
        File pvar 
        File psam 

        String docker = "hkim298/plink_1.9_2.0:20230116_20230707"

        String? out_prefix

        Int memory_gb = 20

        Int disk_size_multiplier = 4
        Int disk_size_addition = 20
    }

    Int disk_size = ceil(size([pgen, pvar, psam], "GB"))*disk_size_multiplier + disk_size_addition

    String out_string = if defined(out_prefix) then out_prefix else basename(pgen, ".pgen")

    command {
        plink2 \
            --pgen ~{pgen} --pvar ~{pvar} --psam ~{psam} \
            --make-bed \
            --out ~{out_string}
    }

    runtime {
        docker: docker
        preemptible: 1
        disks: "local-disk " + disk_size + " HDD"
        memory: memory_gb + " GiB"
    }

    output {
        File convert_Pgen_out_bed = "${out_string}.bed"
        File convert_Pgen_out_bim = "${out_string}.bim"
        File convert_Pgen_out_fam = "${out_string}.fam"
    }
}

## for scope
task PreparePlink{
  input {
    File vcf_file

    String output_prefix

    String plink2_maf_filter = "--maf 0.01"
    String plink2_LD_filter_option = "--indep-pairwise 50000 80 0.1"
    File long_range_ld_file

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
    --out maf_filtered
      
plink2 \
    --pgen maf_filtered.pgen \
    --pvar maf_filtered.pvar \
    --psam maf_filtered.psam \
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

task RunScopeUnsupervised{
    input{

        File bed_file
        File bim_file
        File fam_file

        Int? K
        String output_string
        Int? seed

        Int memory_gb = 60
        String docker = "blosteinf/scope:0.1"
    }

    String plink_binary_prefix =  basename(bed_file, ".bed")
    String relocated_bed = plink_binary_prefix + ".bed"
    String relocated_bim = plink_binary_prefix + ".bim"
    String relocated_fam = plink_binary_prefix + ".fam"

    String unsup_output = output_string + "_unsupervised_" 
    Int disk_size = ceil(size([bed_file, bim_file, fam_file], "GB")  * 2) + 20

    command <<<

ln ~{bed_file} ./~{relocated_bed}
ln ~{bim_file} ./~{relocated_bim}
ln ~{fam_file} ./~{relocated_fam}
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

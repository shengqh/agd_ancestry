mkdir -p /nobackup/h_cqs/shengq2/temp
cd /nobackup/h_cqs/shengq2/temp
java -Dconfig.file=/data/cqs/softwares/cqsperl/config/wdl/cromwell.local_db.conf \
  -jar /data/cqs/softwares/cromwell/cromwell-90.jar \
  run /nobackup/h_cqs/shengq2/program/agd_ancestry/ancestry_pipeline_vcf_scope_unsupervised.wdl \
  -i /nobackup/h_cqs/shengq2/program/agd_ancestry/ancestry_pipeline_vcf_scope_unsupervised.inputs.json \
  --options /data/cqs/softwares/cqsperl/config/wdl/cromwell.options.json
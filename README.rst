|Build Status| |Coverage Status|

Knight Lab Admin Commands for Qiita. This is an Admin only plugin.

This plugin provides the following commands:

#. Sequence Processing Pipeline (SPP):

This command converts the BCL to FASTQ files, does adaptor and host/human sequence filtering, summarizes results and uploads the resulting files to Qiita as one or more new preparations and artifacts. This command is based on [mg-scripts](https://github.com/biocore/mg-scripts) and (metapool)[https://github.com/biocore/kl-metapool/].

The new preparations that the SPP adds contain some automatically generated extra columns:

* raw_reads_r1r2: this is the total number of reads that the sequencing generated for R1 and R2.

* total_biological_reads_r1r2: the total of reads minus the adapters but before host/human filtering for R1 and R2.

* quality_filtered_reads_r1r2: total of reads that passed host/human filtering for R1 and R2.

* fraction_passing_quality_filter: the fraction of reads that are in quality_filtered_reads_r1r2 from the raw_reads_r1r2.


.. |Build Status| image:: https://github.com/qiita-spots/qp-knight-lab-processing/actions/workflows/qiita-plugin-ci.yml/badge.svg
   :target: https://github.com/qiita-spots/qp-knight-lab-processing/actions/workflows/qiita-plugin-ci.yml
.. |Coverage Status| image:: https://coveralls.io/repos/github/qiita-spots/qp-knight-lab-processing/badge.svg?branch=dev
   :target: https://coveralls.io/github/qiita-spots/qp-knight-lab-processing?branch=master

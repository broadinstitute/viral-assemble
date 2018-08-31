task merge_tarballs {
  Array[File]+  in_tarballs

  String?       out_basename="merged"

  command {
    set -ex -o pipefail

    if [ ${length(in_tarballs)} -gt 1 ]; then
      for i in $in_tarballs; do
          if [[ $i != *.tar.gz ]] && [[ $i != *.tgz ]] ; then
              echo "Input file $i does not end with '.tar.gz' or '.tgz'; quitting..."
              exit 1
          fi
      done
      file_utils.py merge_tarballs ${out_basename}.tar.gz ${sep=' ' in_tarballs}  --loglevel DEBUG
    else
      echo "Skipping merge, only one input file"
      ln -s ${select_first(in_tarballs)} ${out_basename}.tar.gz
    fi
  }

  output {
    File  out_bam = "${out_basename}.tar.gz"
  }

  runtime {
    docker: "quay.io/broadinstitute/viral-ngs"
    memory: "2000 MB"
    cpu: 2
    dx_instance_type: "mem1_ssd2_x4"
  }
}
import "tasks_file_utils.wdl" as file_utils

workflow merge_tarballs {
    call file_utils.merge_tarballs
}

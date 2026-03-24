save_job <- function(job, directory = "./") {
  # takes in a string with commands on separate lines, and writes a job.sh file to the chosen directory
  # may still be a bit buggy, so job.sh files need to be inspected closely
  job <- gsub(pattern = "\n", replacement = "' >> job.sh ;\n ", x = job)
  job <- sub(pattern = ">> job.sh", replacement = "> job.sh", x = job)
  job <- paste(job, "' >> job.sh", sep = "")
  job <- paste("echo '", job, sep = "")
  job <- gsub(pattern = ";\n", replacement = ";\n echo '", x = job)
  job <- gsub(pattern = "job.sh", replacement = paste(directory, "job.sh", sep = ""), x = job)
  system(job)
  system(paste0("cat ", directory, "job.sh | sed -e 's/^[ \t]*//' > ", directory, "temp.sh")) # delete all white space to the left of each line
  system(paste0("mv ", directory, "temp.sh ", directory, "job.sh")) # delete white space left to each line
}
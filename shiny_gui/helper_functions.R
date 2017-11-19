# http://jeromyanglim.tumblr.com/post/33418725712/how-to-source-all-r-files-in-a-directory
source_dir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
       if(trace) cat(nm,":")           
       source(file.path(path, nm), ...)
       if(trace) cat("\n")
    }
}
# Evaluate one named section without depending on source-file line numbers.
source_section <- function(file, start, end, envir = parent.frame()) {
  lines <- readLines(file, warn = FALSE)
  start_marker <- paste0("# ---- ", start, " ----")
  end_marker <- paste0("# ---- ", end, " ----")
  first <- which(lines == start_marker)
  last <- which(lines == end_marker)

  if (length(first) != 1L || length(last) != 1L) {
    stop("Section markers must occur exactly once in ", file, call. = FALSE)
  }
  if (last <= first + 1L) {
    stop("Section markers are out of order or the section is empty in ", file, call. = FALSE)
  }

  section <- lines[seq.int(first + 1L, last - 1L)]
  eval(parse(text = section), envir = envir)
}

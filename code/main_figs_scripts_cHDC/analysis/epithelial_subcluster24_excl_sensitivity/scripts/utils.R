make_logger <- function(tag) {
  function(...) cat(sprintf("[%s] %s\n", tag, sprintf(...)))
}

# make_logger(tag) returns a logger that prefixes every message with [tag],
# keeping the sprintf/cat logic defined in exactly one place.
make_logger <- function(tag) {
  function(...) cat(sprintf("[%s] %s\n", tag, sprintf(...)))
}

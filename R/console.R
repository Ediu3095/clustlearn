#' @importFrom cli console_width
console.log <- function(txt, ...) {
  width <- console_width()

  # Split the text into lines
  tmp <- strsplit(txt, '[\r\n]')[[1]]

  # Process every line
  lines <- NULL
  for (line in tmp) {
    # Get the whitespaces at the beginning of the line
    white <- get_whitespace(line)

    # Remove trailing whitespace from the line
    line <- substring(line, nchar(white) + 1)
    if (is.na(line)) line <- ''

    # Expand the tabs in the whitespace
    white <- gsub('\t', '  ', white)

    # Split the line into several lines which will fit in the console
    parts <- strsplit(
      line,
      paste0("(?<=.{", max(width - nchar(white), 1), "})"),
      perl = TRUE
    )[[1]]

    # Add the whitespace to the beginning of each part
    parts <- paste0(white, parts)

    # Add all parts to the processed lines list
    lines <- c(lines, parts)
  }

  # Make sure to print something...
  if (length(lines) < 1) lines <- ""

  # Print all processed lines
  for (line in lines) {
    message(line, ...)
  }
}

hline <- function() {
  console.log(strrep('_', cli::console_width()))
}

get_whitespace <- function(txt) {
  # Find whitespace at the beggining of the string
  fst_match <- gregexpr("^[ \t]*", txt)[[1]]

  # Extract whitespace
  white.length <- attr(fst_match, "match.length")
  substring(txt, 1, white.length)
}

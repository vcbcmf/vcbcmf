tikzprint <- function(fig, file_name, folder, clean = TRUE, view = FALSE, ...) {
  require(tikzDevice)
  pdf_name <- paste0(file_name, ".pdf")
  tex_name <- paste0(file_name, ".tex")
  tikz(file = tex_name, standAlone = TRUE, ...)
  plot(fig)
  dev.off()
  tools::texi2pdf(tex_name, clean = clean)
  file.remove(tex_name)
  file.rename(pdf_name, paste0(folder, "/", pdf_name))
  return(TRUE)
}
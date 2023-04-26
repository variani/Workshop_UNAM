packages = c(
  # data analysis
  'tidyverse', 'glue', 'data.table',
  # genetic data analysis
  'qqman', 'BEDMatrix',
  # rare-variant association analysis
  'SKAT', 'ACAT'
)

is_installed = function(p) nzchar(system.file(package = p))

for(p in packages) {
  if(!is_installed(p)) {
    install.packages(p)
  }
}

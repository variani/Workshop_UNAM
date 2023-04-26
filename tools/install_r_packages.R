packages = c(
  # data analysis
  'tidyverse', 'glue', 'data.table',
  # genetic data analysis
  'qqman', 'BEDMatrix', 'bigsnpr',
  # rare-variant association analysis
  'SKAT', 'ACAT',
  # misc
  'hexbin'
)

is_installed = function(p) nzchar(system.file(package = p))

for(p in packages) {
  if(!is_installed(p)) {
    install.packages(p)
  }
}

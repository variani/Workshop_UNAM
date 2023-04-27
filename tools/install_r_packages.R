packages = c(
  # main
  'devtools',
  # data analysis
  'tidyverse', 'data.table',
  # genetic data analysis
  'qqman', 'BEDMatrix',
  # rare-variant association analysis
  'SKAT', 'ACAT'
)

is_installed = function(p) nzchar(system.file(package = p))

for(p in packages) {
  if(!is_installed(p)) {
    switch(p,
      'ACAT' = devtools::install_github("yaowuliu/ACAT"),
       install.packages(p))
  }
}

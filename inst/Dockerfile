FROM r-base:3.6.3
LABEL maintainer juank.laria@gmail.com

RUN apt-get update && apt-get install r-cran-dplyr r-cran-rcpparmadillo r-cran-mass r-cran-devtools r-cran-devtools r-cran-ggpubr -y

RUN Rscript -e 'devtools::install_github("jlaria/glasp")'

CMD ["/bin/bash"]

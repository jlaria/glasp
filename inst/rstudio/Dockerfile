FROM rocker/rstudio

EXPOSE 8787
ENV PASSWORD=rstudio

RUN R -e "install.packages(c('rlang', 'tune', 'parsnip', 'yardstick', 'rsample', 'workflows', 'recipes', 'timeDate', 'RcppArmadillo', 'fda.usc', 'aricode', 'roxygen2', 'hardhat', 'usethis', 'pkgbuild', 'reshape2', 'testthat'))"

RUN apt-get update && apt-get install -y libglpk40 libxml2 texlive-fonts-extra
RUN apt-get install -y texinfo

COPY entrypoint.sh /usr/bin/
RUN chmod +x /usr/bin/entrypoint.sh
ENTRYPOINT ["entrypoint.sh"]

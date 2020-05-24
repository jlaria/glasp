# glasp
R package glasp

## Install

To install `glasp` in R use

```
devtools::install_github("jlaria/glasp", dependencies = TRUE)
```

### Using docker

There is a docker image available with R 3.6.3 and all the dependencies of glasp.

``` 
docker pull jlaria/glasp:0.0.1
docker run -it jlaria/glasp:0.0.1
```

Optionally, it can be built using the following `Dockerfile`

```
FROM r-base:3.6.3
LABEL maintainer juank.laria@gmail.com

WORKDIR /root/

RUN apt-get update && apt-get install wget git r-cran-dplyr r-cran-rcpparmadillo r-cran-mass r-cran-devtools r-cran-devtools openjdk-8-jre r-cran-ggpubr -y  

# Install clues (removed from cran these days)
RUN wget https://cran.r-project.org/src/contrib/Archive/clues/clues_0.6.2.2.tar.gz
RUN R CMD INSTALL clues_0.6.2.2.tar.gz

# Install R dependencies and glasp
ADD install_glasp.R .
RUN Rscript install_glasp.R 

# Install sparklyr
ADD install_spark.R .
RUN Rscript install_spark.R 

RUN wget https://cran.r-project.org/src/contrib/mvtnorm_1.1-0.tar.gz
RUN R CMD INSTALL mvtnorm_1.1-0.tar.gz

# Bring the experiments from the paper
ADD https://api.github.com/repos/jlaria/glasp-code/commits  /dev/null
RUN git clone https://github.com/jlaria/glasp-code.git && cd glasp-code

CMD ["/bin/bash"]
```


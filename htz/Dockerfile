FROM ubuntu

RUN apt-get update && apt-get install -y r-base-dev wget sudo


WORKDIR /opt

RUN wget https://cran.r-project.org/src/contrib/segmented_0.5-1.4.tar.gz
RUN R CMD INSTALL segmented_0.5-1.4.tar.gz

RUN wget https://cran.r-project.org/src/contrib/Archive/mixtools/mixtools_1.0.0.tar.gz


RUN R CMD INSTALL mixtools_1.0.0.tar.gz

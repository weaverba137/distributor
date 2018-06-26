# FROM debian:latest
# RUN apt-get update --quiet -y && \
#     apt-get install --quiet -y python-flask astropy bokeh
# WORKDIR /app
# COPY distributor.py /app
# ENTRYPOINT ["python"]
# CMD ["distributor.py"]
FROM continuumio/anaconda3:5.0.0

MAINTAINER Benjamin Alan Weaver <baweaver@lbl.gov>

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN groupadd -r desi --gid=58102 && useradd -r -g desi --uid=34448 bweaver

WORKDIR /opt/app
COPY distributor.py ./
# COPY templates templates/

ENV PATH /opt/conda/bin:$PATH

RUN \
    conda update --all --yes && \
    conda config --add channels astropy-ci-extras && \
    conda config --add channels astropy && \
    conda install --yes astropy flask bokeh

# EXPOSE 8000

# COPY docker-entrypoint.sh /usr/local/bin/
# RUN ln -s usr/local/bin/docker-entrypoint.sh /entrypoint.sh # backwards compat

ENTRYPOINT ["python"]
CMD ["distributor.py"]

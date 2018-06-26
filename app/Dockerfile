FROM debian:latest
RUN apt-get update --quiet -y && \
    apt-get install --quiet -y python-flask astropy bokeh
WORKDIR /app
COPY distributor.py /app
ENTRYPOINT ["python"]
CMD ["distributor.py"]

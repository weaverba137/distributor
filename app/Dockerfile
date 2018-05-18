FROM debian:latest
RUN apt-get update --quiet -y && apt-get install --quiet -y python-flask
WORKDIR /app
COPY app.py /app
ENTRYPOINT ["python"]
CMD ["app.py"] 

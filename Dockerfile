FROM ubuntu:latest

RUN apt-get update

COPY simulate.pl .
COPY Simulate.pm .

RUN chmod +x simulate.pl

CMD ["perl", "simulate.pl"]
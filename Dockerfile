
# sudo docker build -t sequence-simulator:0.1 .
# sudo docker run f95facbc6fae
# sudo docker container run -it 25b380aa6cdd /bin/bash

FROM ubuntu:latest

RUN apt-get update -y
RUN apt-get install -y make gcc libc-dev build-essential
RUN apt-get install -y cpanminus

RUN cpan Math::Random
RUN cpan JSON::Parse

COPY simulate.pl /usr/local/bin/simulate
COPY Simulate.pm /usr/local/bin/Simulate.pm

RUN chmod +x /usr/local/bin/simulate
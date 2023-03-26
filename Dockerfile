FROM python:3.9-slim

# create user with a home directory
ARG NB_USER
ARG NB_UID
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
WORKDIR ${HOME}

COPY . ${HOME}

RUN apt-get update && \
    apt-get install -y git wget clustalo

RUN wget https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb && \
    apt-get install -y ./google-chrome-stable_current_amd64.deb && \
    rm -rf google-chrome-stable_current_amd64.deb

RUN pip install --no-cache --upgrade pip && \
    pip install --no-cache notebook jupyterlab && \
    pip install --no-cache -r requirements.txt

RUN git config --global --add safe.directory '*'

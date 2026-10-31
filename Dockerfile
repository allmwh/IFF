FROM allmwh/iff-base:2026

COPY . ${HOME}

USER root
RUN chown -R ${NB_UID}:${NB_UID} ${HOME}

USER ${NB_USER}
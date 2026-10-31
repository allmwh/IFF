FROM allmwh/iff-base:2026

COPY . ${HOME}

USER root
COPY --chown=${NB_USER}:${NB_USER} . ${HOME}

USER ${NB_USER}
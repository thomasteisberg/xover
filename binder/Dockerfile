FROM python:3.12-slim

RUN pip install uv

# Binder/repo2docker passes these build args
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER=${NB_USER} HOME=/home/${NB_USER}
RUN useradd --create-home --shell /bin/bash --uid ${NB_UID} ${NB_USER}

WORKDIR ${HOME}

# Build context is the repo root when Dockerfile lives in binder/
COPY --chown=${NB_USER}:${NB_USER} pyproject.toml uv.lock ./
COPY --chown=${NB_USER}:${NB_USER} xover/ ./xover/

USER ${NB_USER}
RUN uv sync --frozen --no-dev

ENV PATH="${HOME}/.venv/bin:${PATH}"

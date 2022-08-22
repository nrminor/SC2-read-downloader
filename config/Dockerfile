FROM continuumio/miniconda3

SHELL ["/bin/bash", "-c"]

ADD sra-downloader.yaml /tmp/sra-downloader.yaml

RUN conda env create -f /tmp/sra-downloader.yaml -n sra-downloader python=3.8 \
	&& conda clean --all -f --yes \
	&& conda init bash \
	&& echo "conda activate sra-downloader" > ~/.bashrc

ENV PATH /opt/conda/envs/sra-downloader/bin:$PATH

ENV CONDA_DEFAULT_ENV sra-downloader

CMD ["/bin/bash"]

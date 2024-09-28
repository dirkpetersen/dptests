#!/bin/bash

# Script configuration
SAMTOOLS_VERSION="1.17"
IMAGE_NAME="yourdockerusername/samtools"
BIOCONTAINER_REPO="biocontainers/samtools"
DOCKERHUB_USERNAME="yourdockerusername"
DOCKERHUB_PASSWORD="yourpassword"

# Create a Dockerfile
cat <<EOF > Dockerfile
FROM debian:bullseye-slim

LABEL base.image="debian:bullseye-slim" \
      version="$SAMTOOLS_VERSION" \
      software="samtools" \
      software.version="$SAMTOOLS_VERSION" \
      about.summary="Samtools is a suite of programs for interacting with high-throughput sequencing data" \
      about.home="http://www.htslib.org/" \
      about.documentation="http://www.htslib.org/doc/samtools.html" \
      about.license="https://opensource.org/licenses/BSD-3-Clause" \
      extra.identifiers.biotools="samtools" \
      extra.biotools="samtools" \
      about.tags="Genomics, NGS, BAM, CRAM, VCF, FASTQ"

# Install Samtools
RUN apt-get update && apt-get install -y \
    wget \
    build-essential \
    libncurses5-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    && wget https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VERSION/samtools-$SAMTOOLS_VERSION.tar.bz2 \
    && tar -xjf samtools-$SAMTOOLS_VERSION.tar.bz2 \
    && cd samtools-$SAMTOOLS_VERSION \
    && ./configure \
    && make \
    && make install \
    && cd .. \
    && rm -rf samtools-$SAMTOOLS_VERSION samtools-$SAMTOOLS_VERSION.tar.bz2 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set environment variables
ENV PATH="/usr/local/bin:\$PATH"

# Set the entrypoint to samtools
ENTRYPOINT ["samtools"]
EOF

# Build the Docker image
echo "Building Docker image..."
docker build -t $IMAGE_NAME:$SAMTOOLS_VERSION .

# Log in to DockerHub
echo "Logging in to DockerHub..."
echo "$DOCKERHUB_PASSWORD" | docker login -u "$DOCKERHUB_USERNAME" --password-stdin

# Tag the image for DockerHub
echo "Tagging image for DockerHub..."
docker tag $IMAGE_NAME:$SAMTOOLS_VERSION $DOCKERHUB_USERNAME/samtools:$SAMTOOLS_VERSION

# Push the image to DockerHub
echo "Pushing image to DockerHub..."
docker push $DOCKERHUB_USERNAME/samtools:$SAMTOOLS_VERSION

# Optional: Tag the image for Biocontainers (replace this if you push to another registry)
docker tag $IMAGE_NAME:$SAMTOOLS_VERSION $BIOCONTAINER_REPO:$SAMTOOLS_VERSION

# Push the image to Biocontainers (Optional)
# docker push $BIOCONTAINER_REPO:$SAMTOOLS_VERSION

echo "Docker image for Samtools $SAMTOOLS_VERSION has been built and pushed to DockerHub!"


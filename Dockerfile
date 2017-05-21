# Use an official Python runtime as a base image
FROM ubuntu:16.04

# Set the working directory to /app
WORKDIR /app

# Copy the current directory contents into the container at /app
ADD . /app

# Install any needed packages specified in requirements.txt
RUN apt-get update && apt-get install -y \
    python2.7 \
    python-pip \
    python-dev \
    libmysqlclient-dev \
    r-base

# install R dependencies
RUN echo 'install.packages(c("ggplot2","grid","gridExtra","VennDiagram"), repos="http://cran.us.r-project.org",dependencies=TRUE)' > /tmp/packages.R \
	&& Rscript /tmp/packages.R

# install python dependencies
RUN pip install -r dependencies.txt


# Run app.py when the container launches
CMD ["make"]


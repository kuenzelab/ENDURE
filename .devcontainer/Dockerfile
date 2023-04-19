# From python:3.9.16-bullseye
FROM python:3.9.16-bullseye

# Update aptitude with new repo
#RUN apt-get update

# Install software 
#RUN apt-get install -y git



# Copy Rosetta 
RUN mkdir /app
RUN mkdir /app/ENDURE
COPY lib /app/ENDURE/lib

# Make dir for ENDURE
COPY requirements.txt /tmp/requirements.txt
RUN pip install -r /tmp/requirements.txt --ignore-installed

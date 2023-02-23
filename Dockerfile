# From python:3.9.16-bullseye
FROM python:3.9.16-bullseye

# Update aptitude with new repo
#RUN apt-get update

# Install software 
#RUN apt-get install -y git


# Make dir for ENDURE
RUN mkdir /app

# Clone from https://github.com/kuenzelab/ENDURE.git to /app 
WORKDIR /app
# Currently the clone needs to be done with a personal access token from github the classic ones
RUN git clone https://engelberger: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@github.com/kuenzelab/ENDURE.git
# Set working directory to ENDURE 
WORKDIR /app/ENDURE

# Copy lib/rosetta_linux from local to /app/ENDURE/lib/rosetta_linux 
# This should be changed to a download from the rosetta website with wget and the credentials from the rosetta commons website
COPY lib/rosetta.tar.gz /app/ENDURE/lib/rosetta.tar.gz

# Extract rosetta.tar.gz to /app/ENDURE/lib/rosetta_linux
RUN tar -xvf /app/ENDURE/lib/rosetta.tar.gz -C /app/ENDURE/lib/

# Remove rosetta.tar.gz
RUN rm /app/ENDURE/lib/rosetta.tar.gz

# Copy requirements.txt from local to /app/ENDURE overwriting the existing one
COPY requirements.txt /app/ENDURE/requirements.txt

# Install dependencies
RUN pip install -r requirements.txt  --ignore-installed

# Copy the Welcome.py from local to /app/ENDURE/overwriting the existing one
COPY Welcome.py /app/ENDURE/Welcome.py

# Copy the pages\d_Energy_Heatmap.py from local to /app/ENDURE/pages/overwriting the existing one
COPY pages/d_Energy_Heatmap.py /app/ENDURE/pages/d_Energy_Heatmap.py


# Set working directory to /app/ENDURE
WORKDIR /app/ENDURE

# Delete the /app/ENDURE/pages/testing/
RUN rm -rf /app/ENDURE/pages/testing/

# Expose port 8501
EXPOSE 8501

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

ENTRYPOINT ["python3.9", "-m", "streamlit", "run", "Welcome.py", "--server.port=8501", "--server.address=0.0.0.0"]

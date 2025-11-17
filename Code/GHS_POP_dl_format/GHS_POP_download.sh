#!/bin/bash
#!/bin/bash
#SBATCH --partition=zen4
#SBATCH --time=4:00:00
#SBATCH --mem=4G

for year in {1975..2030..5}
do
    echo "$year"
    wget https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/GHS_POP_E${year}_GLOBE_R2023A_4326_30ss/V1-0/GHS_POP_E${year}_GLOBE_R2023A_4326_30ss_V1_0.zip -P /scratchu/tmandonnet/GHS-POP/
    unzip /scratchu/tmandonnet/GHS-POP/GHS_POP_E${year}_GLOBE_R2023A_4326_30ss_V1_0.zip -d /scratchu/tmandonnet/GHS-POP/
    rm /scratchu/tmandonnet/GHS-POP/GHS_POP_E${year}_GLOBE_R2023A_4326_30ss_V1_0.zip
done
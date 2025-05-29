#!/bin/bash

# Define the base paths
# mac_path="Users/deszoeks/Library/CloudStorage"
from_path_stem="/Users/deszoeks/Library/CloudStorage/OneDrive-OregonStateUniversity/ekamsat"
to_path_stem="/Users/deszoeks/Library/CloudStorage/OneDrive-SharedLibraries-UW/og_ekamsat-st - Documents/Data/radiosonde"

# Define the list of IMD station folders
folders=("ahmedabad" "bhubaneswar" "chennai" "karaikal" "kochi" "kolkata" "machilipatnam" "mangalore" "minicoy" "portblair" "pune" "santacruz" "trivandrum" "visakhapatnam")

# Loop over each folder and run rsync
for folder in "${folders[@]}"; do
    echo "START ${folder}"
    rsync -av "${from_path_stem}/${folder}/" "${to_path_stem}/${folder}/"
    echo "END ${folder}"
done

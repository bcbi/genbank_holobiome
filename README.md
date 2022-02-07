# Genbank Database Pipeline - Version 2 
#### Authors: Vivek Ramanan and Shanti Mechery 

The GenBank DataBase Pipeline has two purposes: to map source site data
from the National Center for Biotechnology Information's GenBank Database
to their preferred name according to the the National Institute of Health's
United Medical Language System and to create a matrix of bacteria species
and host species data, both also from the National Center for Biotechnology
Information's GenBank Database, pertaining to a specific source site. The first
part is the Upstream Process and the second part is the Downstream Process.
This matrix created in the Downstream Process will be a part of a NEXUS file,
which can be used as input for PAUP*, a phylogenetic analysis software that
creates phlyogentic trees. Each section of this guide explains the processes
in further detail. If the most up-to-date input data from the GenBank database
and the UMLS for the pipeline are needed, start with the Upstream Process.
Note that this process only needs to be run once. If all the data is up-to-date,
only the Downstream process needs to be run. The last GenBank download was 
performed in November 2020. 

The data from the GenBank Database can be downloaded using a GenBank Loader to
a MySQL database (details in the github repo BCBI/genbank_loader). In the GenBank
Database, an entry has the following information: Locus, Organism, Source Site 
(isolation_source or tissue_type), and Host. The GenBank Loader imports this data into
a MySQL database named 'genbank' with the following tables: annotations, authors, 
dbxrefs, journals, keywords, and basic. Only the annotations table is needed, with the 
following columns: partitionKey, locus, name, indexedValue, and value. Each row of the
table contains a single piece of type information (name column) with a value (value column)
pertaining to a locus number (locus column). The types of information included in the
name column are organism, source site (isolation_source or tissue_type) and host. 
For the pipeline, use a cleaned version of the host data, which is found in the table 
cleaned_host_annotations. The last cleaned_host_annotations was created in March 2021
and can be updated accordingly using the hostClean scripts included. 

Please refer first to the "GenBank_Write_Up.pdf" for the first version 
of this pipeline. 

Dependencies: 
1. Julia
2. Python3
    - BioPython (Entrez) and NCBI Entrez Account information
    - Pandas and Numpy
    - SKLearn, markov_clustering, and NetworkX
3. R
    - ggTree and treeIO
4. TimeTree.org

** If using Oscar, run ALL of these only using Batch jobs (or interact jobs for testing) **

## Data Download and Processing Process: 

This section processes the source site data (isolation_source or tissue_type) into
CSV files using the UMLS MetaMap program. This process must be run individually
for isolation_source and tissue_type. MetaMap installation can be found at 
https://lhncbc.nlm.nih.gov/ii/tools/MetaMap/documentation/Installation.html and 
can be downloaded with the binary onto Oscar or any other remote server. This section
is built in Julia solely and requires some MySQL commands as well. 

1. Downloading the GenBank data (if needs to be updated): 
    - Refer to https://github.com/bcbi/genbank_loader
    - Original source code: https://bitbucket.org/UVM-BIRD/genbank-loader/src

2. MySQL Queries for source site and bacteria species: 

The following files have been given with their source type from annotations: 
    - isolationLocus.csv: isolation_source
    - tissueLocus.csv: tissue_type
    - allBacteriaSpecies.csv: organism

General Query: 
```
select locus,value from annotations where name="isolation_source"
select locus,value from annotations where name="tissue_type"
select locus,value from annotations where name="organism"
```

HOWEVER, if you are not able to access the MySQL queries directly or more likely, 
because the tables are too large to run directly on the command line, please use
the following command line batch jobs instead. 

ex: IsolationLocus.CSV command
First create the table in mySQL if it does not exist: 
$ create table isolationLocus(locus varchar(20), value longtext)
Then we run the following command outside of mySQL from command line: 
```
mysql --enable-cleantext-plugin -h <DBName> -u <username> -p <password> -D genbank
-e "insert into isolationLocus select locus,value from annotations where name='isolation_source';"
```

This will create the mySQL table in the database properly. Sometimes you will get a timed
out output, but this is likely because the table finished loading and the mySQL server
timed out after that so it should've run correctly. 

**IF YOU ARE UPDATING cleaned_host_annotations: also run hostLocus.csv ("where name='host'")

3. Downloading the mySQL tables:

If you are working off of command line and on a remote server, and cannot download the mySQL
tables directly, you will need to use MySQLDump. 
```
mysqldump --enable-cleartext-plugin -h <DBname> -u <username> -p <password> genbank tissueLocus > tissueLocus.sql
```
** occasionally I have noticed that the space between "-p" and the password can cause errors. 
Try using without a space between the -p flag and the password if this occurs. **

4. Converting .SQL files into .CSV

To convert .SQL files into .CSV, I used a preexisting Python program called mysqldump-to-csv 
that I updated for this project specifically. 
```
python3 mysqldump-to-csvV2.py <file.SQL>
```
Now the data should be in a CSV format and ready to process moving forward. These are the files
you should have: 
    - isolationLocus.csv
    - tissueLocus.csv
    - allBacteriaSpecies.csv
    - (if updating cleaned_host_annotations: hostLocus.csv)
The Locus files are the ones that will go through the upstream process (not allBacteriaSpecies.csv)

## Upstream Process:

1. upstream_1_valueCSV_parser.jl

This file processes CSVs into TXT files for MetaMap to take in as input. 
Input: source site CSV (ex: tissueLocus.csv, isolationLocus.csv, or hostLocus.csv)
Output: appropriate MetaMap input format for text file (sourceSite_MetaMap_input.txt)
```
julia upstream_1_valueCSV_parser.jl sourceSiteLocus.csv sourceSite_MetaMap_input.txt
```
2. Running MetaMap
```
$ ./public_mm/bin/wsdserverctl start
$ ./public_mm/bin/skrmedpostctl start
$ ./public_mm/bin/metamap -I sourceSite_MetaMap_input.txt sourceSite_MetaMap_output.txt
$ ./public_mm/bin/wsdserverctl stop
$ ./public_mm/bin/skrmedpostctl stop
```
** In particular, the isolationLocus.CSV was much larger than the other files for me. I
had to split this file into 5 files to get MetaMap to run on each one individually,
then joined them all together after for the following steps. 

3. upstream_2_metamapParser.jl

Each UMLS concept has a Semantic Type and this script filters for the Semantic Types that we want,
which is specific to body parts and organs. 
** If updating cleaned_host_annotations, please use semanticHost.txt here instead **

Input: sourceSite_MetaMap_output.txt
Output: 1) sourceSite_MetaMap_parsed.csv
            - columns: original, preferred, semanticType
            - original: original source site from Genbank
            - preferred: source site preferred MEtaMap name, all concept names concatenated
            - semanticType: semantic type associated with each MetaMap concept
        2) sourceSite_CUI.csv
            - CUI: concept unique identifiers
```
$ julia upstream_2_metamapParser.jl sourceSite_MetaMap_output.txt semanticTypes.txt sourceSite_MetaMap_parse.csv sourceSite_CUI.csv
```
4. upstream_3_createBacteriaSpecieFile.jl

Now the appropriate bacteria species data need to be filtered out from allBacteriaSpecies.csv
file for each source site. This will be a large file. 
```
$ julia upstream_3_createBacteriaSpecieFile.jl sourceSiteLocus.csv allBacteriaSpecies.csv sourceSite_bacteriaSpecies.csv
```
5. upstream_4_speciesPreferredMapping.jl

Finally, we can merge the bacteria species data for each source site from upstream_3 and the 
metamap preferred names from upstream_2 metamap parsing. 
```
$ julia upstream_4_speciesPreferredMapping.jl sourceSite_MetaMap_parsed.csv sourceSite_bacteriaSpecies.csv sourceSite_bacteriaSpecies_preferred.csv
```
You finished the Upstream Processes! 

## Updating Cleaned_Host_Annotations

If you are updating this file, you should now have the host_bacteriaSpecies_preferred.csv file - this
is after running hostLocus.csv through all the upstream steps and using semanticHost.txt for the processing 
of MetaMap concepts. 

1. hostClean_1_cleaningHosts.py

This script is a first manual run-through of cleaning the host names. It also changes the column names
of the host_bacteriaSpecies_preferred.csv file to the columns of the original cleaned_host_annotations.csv
file and adds a cleaned "updatedScientificName" column that is the most up-to-date host name. 

```
python3 hostClean_1_cleaningHosts.pyy hostData/host_bacteriaSpecies_preferred.csv hostTaxDict.txt
```

Input: 1) host_bacteriaSpecies_preferred.csv
        2) hostTaxDict.txt
            This is a preexisting taxonomic dictionary that I used for the host names. No need to replicate. 
        3) hostSpeciesNames.txt
            A manual creation for common misspellings of the big host categories 
Output: updated_cleanedhost.csv

2. hostClean_2_searchHosts.py

This script is an extra check of the host names in Entrez to find scientific names for all the hosts. This is
not always successful however due to variation and grammatical issues in the way the host name was entered.
I used it because it was valuable in coordinating common names vs. scientific names, however when the host
names are kept as genus rather than species, we did not merge them because of different taxonomical levels.  

```
python3 hostClean_2_searchHosts.py cleaned_hostData.csv
```

Input: updated_cleanedhost.csv
Output: cleaned_hostData.csv

## Downstream Processes: 

1. downstream_1_createHostSpecies.jl

This script maps the data from sourceSite_bacteriaSpecies_preferred.csv to the data in 
cleaned_host_annotations.csv
```
$ julia downstream_1_createHostSpecies.jl cleaned_host_annotations.csv sourceSite_bacteriaSpecies_preferred.csv sourceSite_Host_Bacteria.csv
```
2. downstream_2_mergeAndGroup.py

This script takes the two source site CSVs (isolation_source and tissue_type) and merges them together. 
It also adds a grouping category for the source site, in which body organs of interest (GI in particular), 
are grouped and given group names in a separate column "group". 

Input: 1) isolation_Host_Bacteria.csv
        2) tissue_Host_Bacteria.csv
        3) sourceGroupings.txt
Output: merged_Host_Bacteria.csv
```
$ python3 downstream_2_mergeAndGroup.py isolation_Host_Bacteria.csv tissue_Host_Bacteria.csv sourceGroupings.txt merged_Host_Bacteria.csv 
```
### Cleaning Species Names and Processing Taxonomically

The bacterial species names are quite varied and also do not allow us to look at groupings of bacteria 
(such as family or phyla) rather than individual species. Because of that, we have an additional 
microorganism species cleaning step to add taxonomical data on top of our current dataset.  

1. speciesClean_1_splitSpecies.py

```
$ python3 speciesClean_1_splitSpecies.py merged_Host_Bacteria.csv microSpeciesNames.txt
```
This script simply takes the species names in the merged_Host_Bacteria.csv and sorts them into their 
value counts (most abundant at top and least at bottom). Then, it writes them to the file microSpeciesNames.txt. 
Finally, using the split command, it splits this file into files of 15,000 species in total so that the
next step can be run sequentially. 

2. speciesClean_2_searchSpecies.py

This program runs the Entrez search to find each microorganism's taxonomic information, 
aka it takes a super long time and can often time out from NCBI's side if you run it at 
the wrong time or for too long. I suggest splitting up the microSpeciesNames.txt file
up into around 15,000 entries each (done in prev step) and run on weekdays after 5 pm or weekends. 

Input: a text file with all of the unique species names available
Output: taxonomyDict.txt in the same directory
```
$ python3 speciesClean_2_searchSpecies.py speciesNameFile.txt
```
3. speciesClean_3_updateCSV.py

Now, we update the merged_Host_Bacteria.csv file with the taxonomy information generated by
the speciesClean_1 step and the taxonomyDict.txt file. I recommend running a "wc -l taxonomyDict.txt"
before hand to make sure that the length of that file matches the length of the original species
names file (before splitting) to make sure the file is correct. 

Output: merged_Host_Bacteria_Species.csv 
```
$ python3 speciesClean_3_updateCSV.py merged_Host_Bacteria.csv taxonomyDict.txt 
```
#### Going back to the Downstream Process

3. downstream_3_createMatrix.jl

In downstream 3, we create the matrix file that PAUP will input to create a phylogenetic tree. 
This script has been updated from the previous version to include more command line arguments
to customize the grouping and host requirements before making the PAUP file. 
```
$ julia downstream_3_createMatrix.jl <search term> merged_Host_Bacteria_Spcies.csv <output.nex> <group or specific> <taxonomic resolution> <row number limit>
```
- Search Term: this will be the term that is searched in either the group or species column
    - Group example: GI
    - Specific example: Intestine, Stomach
- Output.nex: name of output file in NEXUS format
- Group or Specific: decides whether you will be searching by Group or by specific source site
    - This is a binary switch - if group is NOT input, it will use the source site column
- Taxonomic resolution: species, genus, family, order, class, phylum 
- Row number limit: this ensures that every host included in the file has AT LEAST the row number limit 
    entries of microorganisms (I used 4 for GI)

Outputs: 
    1) output.nex
    2) hostFiles/<group>_coverage.txt

4. downstream_4_coverage.py

Downstream 4 creates a coverage file from the raw coverage file created by downstream 3 (group_coverage.txt) 
in hostFiles. This coverage file is used by downstream 5 to ensure equal coverage amongst hosts that have
high amounts of data and hosts that have low amounts of data. 
```
$ python3 downstream_4_coverage.py hostFiles/<group>_coverage.txt <output.txt>
``` 
5. downstream_5_dataType.py

Finally, downstream 5 updates the matrix nexus file with a few pieces of information for final visualization.
```
$ python3 downstream_5_dataType.py <outputDS3.nex> hostUpdates.csv <output names> <upper threshold for coverage> <outputDS4.txt>
```
- outputDS3.nex: output from Downstream 3
- hostUpdates.csv: This is a CSV I manually created to ensure which hosts ended up in the final file
    - Columns: Newick_label, updatedScientificName, family, order, class, correctedName
    - Newick_label: the updatedScientificName but with all spaces replaced with underscore so it is one word 
    - correctedName: in case the updatedScientificName is incorrect, I used this to update the name
        or merge two labels together (such as one scientific name and one common name)
- output names: there are multiple output files, so I inputted a phrase that will start all the output files
    - Output files: <group>_updated.nex and hostFiles/<group>_updated.txt
    - First is the nexus file for use, second is a list of all hosts 
    - The second file is necessary for finding out the outgroup in tree creation
- upper threshold for coverage: all hosts ABOVE this number will be coverage updated to this number
    - EX: 20 for GI
- outputDS4.txt: output coverage file from downstream 4

# Visualizaton

1. vis1_speciesComposition.py

This script finds the species composition at the **SPECIES** level. To get the species file, make sure to run
downstream step 3 (need not run downstream steps 4 and 5 because coverage normalization is not required here).
This is an example command:
```
$ julia downstream_3_createMatrix.jl GI merged_Host_Bacteria_Spcies.csv gi_species.nex group species 4
```

Using this output file (gi_species.nex), we can analyze the number of unique species present in each phylum. 
This script will create a percentage based stacked bar chart to visualize the output of host groups of interest.

```
$ python3 vis1_speciesComposition.py gi_species.nex hostUpdates.csv taxonomyDict.txt <outputFile>
```
hostUpdates.csv and taxonomyDict.txt are both located in the Data/ folder.

2. vis2_diversityAnalysis.py

Runs a hill numbers approach of biodiversity estimation across the sample dataset. 
```
python3 vis2_diversityAnalysis.py gi_phyla.nex hostUpdates.csv <outputFile>
```
This can be run on both the phyla and species level files to contrast the difference between 
biodiversity estimation in both. 

3. vis3_networkClustering.py

Runs the Markov Clustering Algorithm on phylum level data to find clustering groups of hosts. 
Also analyzes cluster microbiome composition similarities as well as creates three visualization
files in a directory of choice: MCL cluster network, presence-absence general network, and 
mammalian similarity network. Cosine similarity is the metric for pairwise similarity. 
```
python3 vis3_networkClustering.py gi_phyla.nex hostUpdates.csv <outputDirectory>
```
All network files can be visualized in CytoScape. The basic process I used in Cytoscape is the following:
- Load in the network
- Layout -> Grouped Layout -> Cluster / Taxonomy (cluster for MCL, taxonomy for other two networks)
- Tools -> Analyze Network
- Style -> coloring based on taxonomy and node size based on degree

4. vis4_phyloTree.r - Steps for Phylogenetic Trees

The creation of the microbiome tree requires a few different steps: 
- Using gi hosts file (output from Downstream 5) in timetree.org to find outgroup (I used Octopus_vulgaris)
- Export this tree in newick format to save as true evolutionary tree for hosts
    - there may be hosts that do not have the same naming conventions as hosts in the dataset 
    - I dealt with this by manually including data for these new host names
- Running the phylum presence absence data (gi_phyla.nex) through PAUP neighbor joining with outgroup
- Visualizing both microbiome and evolutionary tree with vis4_phyloTree.r

#### Timetree.org
Input in host file (each line is different host)
** make sure host names (genus and species) are separated by a space and NOT UNDERSCORE **
Note down oldest evolutionary branch and organism (Octupus_vulgaris)
Export tree in newick format in bottom left

#### PAUP Neighbor Joining Commands: 
Run PAUP neighbor joining on the gi_phyla.nex data
```
paup
execute <filename> 
set criterion=distance maxtrees=200 increase=no
hsearch
set root=outgroup
outgroup <name of single host or multiple>
nj brlens=yes treefile=<outputFile>
```
#### Tree Visualization: vis4_phyloTree.r

Visualization of trees using TreeIO and GGTree, using a slightly different hostUpdates CSV file.
I modified this to include hosts with different names from TimeTree.org. 
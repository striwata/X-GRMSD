# X-GRMSD: Exact solution for RMSD minimization with G-RMSD to Determine Molecular Similarity

## Introduction
This repository contains the C++ code to determine similarity of molecules in 3D. It is free software under the terms of the MIT License. Details of the algorithms can be found in our paper:
Tomohiro Nabika, Satoru Iwata, Hiroko Satoh, "Exact solution for minimization of root mean square deviation with G-RMSD to determine molecular similarity", Bulletin of the Chemical Society of Japan, Volume 97, Issue 4, April 2024.
## Compiling
Use cmake to generate desired projects on different platforms.

## Terminology
Data points: points of the source point set to be transformed.

Model points: points of the target point set.

## Running

There are two types of problems: G-RMSD and Substructure Enumeration. 

### G-RMSD 

Run the compiled binary with following parameters: \<METHOD NAME\> \<MODEL FILENAME\> \<DATA FILENAME\> \<IF PERMIT MIRROR ISOMER\>  \<OUTPUT FILENAME\> \<If MULTIPLE MODEL\>, e.g. “./GRMSD MatchFastOpt ./input/model/model_same_random_GRMSD.csv ./input/data/data_same_random_GRMSD.csv 1 ./output/MatchFastOpt_same_random_GRMSD.txt” 1.

\<METHOD NAME\> is the name of METHOD. In this repository, there are AO, TSR, IsometryOpt and MatchFastOpt.
 
\<MODEL FILENAME\> and \<DATA FILENAME\> are the point files of the model and data pointsets respectively. Each point file is in plain text or csv format. It begins with a positive point number N in the first line, followed with N lines of X, Y, Z values of the N points.

\<IF PERMIT MIRROR ISOMER\> is the logical value of whether mirror isomers are permitted. If this value is 1, mirror isomers are considered the same thing. If this value is 0, mirror isomers are considered the different things.

\<OUTPUT FILENAME\> is the output file containing registration results. 

\<If MULTIPLE MODEL\> is the logical value of whether there are multiple models. If this value is 1, there is one model for each data. If this value is 0, there is one model for all data.
 
### Substructure Enumeration

Run the compiled binary with following parameters: \<METHOD NAME\> \<MODEL FILENAME\> \<DATA FILENAME\> \<IF PERMIT MIRROR ISOMER\> \<SUBSTRUCTURE DETERMINATION VALUE\> \<OUTPUT FILENAME\> \<If MULTIPLE MODEL\>, e.g. “./Substructure MatchFPT ./input/model/model_random_Substructure.csv ./input/data/data_random_Substructure.csv 1 0.03 ./output/random_Substructure.txt” 1.

\<METHOD NAME\> is the name of method. In this repository, there are IsometrySearch, MatchFPT, and Improved-MatchFPT.
 
\<MODEL FILENAME\> and \<DATA FILENAME\> are the point files of the model and data pointsets respectively. Each point file is in plain text or csv format. It begins with a positive point number N in the first line, followed with N lines of X, Y, Z values of the N points.

\<IF PERMIT MIRROR ISOMER\> is the logical value of whether mirror isomers are permitted. If this value is 1, mirror isomers are considered the same thing. If this value is 0, mirror isomers are considered the different things.

\<SUBSTRUCTURE DETERMINATION VALUE\> is the value which determines whether data have a model as a substructure. Please read the paper for details.

\<OUTPUT FILENAME\> is the output file containing registration results.

\<If MULTIPLE MODEL\> is the logical value of whether there are multiple models. If this value is 1, there is one model for each data. If this value is 0, there is one model for all data.

# G-RMSD: Exact solution for RMSD minimization between molecules

## Introduction
This repository contains the C++ code to determine similarity of molecules in 3D. It is free software under the terms of the MIT License. Details of the algorithms can be found in our paper:
*This paper has not been published yet
## Compiling
Use cmake to generate desired projects on different platforms.

## Terminology
Data points: points of the source point set to be transformed.

Model points: points of the target point set.

## Running

There are two types of problems: G-RMSD and Substructure Search. 

### G-RMSD 

Run the compiled binary with following parameters: \<METHOD NAME\> \<MODEL FILENAME\> \<DATA FILENAME\> \<IF PERMIT MIRROR ISOMER\>  \<OUTPUT FILENAME\>, e.g. “./GRMSD MatchFastOpt model/random_GRMSD.txt data/random_GRMSD.txt 1 output/random_GRMSD.txt”.

\<METHOD NAME\> is the name of METHOD. In this repository, there are AO, TSR, IsometryOpt and MatchFastOpt.
 
\<MODEL FILENAME\> and \<DATA FILENAME\> are the point files of the model and data pointsets respectively. Each point file is in plain text or csv format. It begins with a positive point number N in the first line, followed with N lines of X, Y, Z values of the N points.

\<IF PERMIT MIRROR ISOMER\> is the logical value of whether mirror isomers are permitted. If this value is 1, mirror isomers are considered the same thing. If this value is 0, mirror isomers are considered the different things.

\<OUTPUT FILENAME\> is the output file containing registration results. 
 
### Substructure Search

Run the compiled binary with following parameters: \<METHOD NAME\> \<MODEL FILENAME\> \<DATA FILENAME\> \<IF PERMIT MIRROR ISOMER\> \<SUBSTRUCTURE DETERMINATION VALUE\> \<OUTPUT FILENAME\>, e.g. “./Substructure MatchFPT model/random_Substructure.txt data/random_Substructure.txt 1 0.03 output/random_Substructure.txt”.

\<METHOD NAME\> is the name of method. In this repository, there are IsometrySearch, MatchFPT, and Improved-MatchFPT.
 
\<MODEL FILENAME\> and \<DATA FILENAME\> are the point files of the model and data pointsets respectively. Each point file is in plain text or csv format. It begins with a positive point number N in the first line, followed with N lines of X, Y, Z values of the N points.

\<IF PERMIT MIRROR ISOMER\> is the logical value of whether mirror isomers are permitted. If this value is 1, mirror isomers are considered the same thing. If this value is 0, mirror isomers are considered the different things.

\<SUBSTRUCTURE DETERMINATION VALUE\> is the value which determines whether data have a model as a substructure. Please read the paper for details.

\<OUTPUT FILENAME\> is the output file containing registration results.

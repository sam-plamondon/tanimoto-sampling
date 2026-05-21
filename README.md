Copyright (c) 2022-2026 Samuel Plamondon

# Overview

This script takes one or more lists of SMILES strings and trims them to a target size using a sampling method that optimizes for having the lowest Tanimoto similarity score between each pair of compounds. This leads to a relatively diverse set of compounds.

# Algorithm

The lists are filtered using a threshold Tanimoto score which is optimized iteratively via a bisection search.

The algorithm works as follows:

1. The list of compounds is sorted randomly.

2. An initial Tanimoto score (0.5) is guessed.

3. The list of compounds is iterated through one by one. The Tanimoto score between the marginal compound and all previous compounds is calculated. If any one of those scores is above 0.5, the marginal compound is thrown out. Thus, any two compounds in the resulting list will have a Tanimoto score of 0.5 or less.

4. The number of compounds in the trimmed list is counted. If this number is within 10% of the target number, the algorithm is exited and the list is output to file. (The 10% error threshold can be changed by changing the value of "err_threshold" on line 102). If not, a new Tanimoto score threshold is guessed. This is done via standard bisection search, with an initial minimum threshold of 0.0 and an initial maximum threshold of 1.0. Thus, if the trimmed list of compounds is too small (i.e. the Tanimoto threshold is too low), we return to step #2 with a guess of 0.75; if the trimmed list is too large (i.e. the Tanimoto threshold is too high), we return to step #2 with a guess of 0.25.

5. Steps #2 to #4 are repeated until the percent error is low enough, or until the maximum and minimum thresholds are within 0.001 of each other.

# Input files

This script takes as an input one or more CSV files, each with at least two columns, one of SMILES strings and one of corresponding compound IDs (these are the only two columns that will be retained in the output files). The order of the columns does not matter, nor does the presence of additional columns. There must be a header row present.

Each input CSV file should contain more entries than the target size.

Input parameters are defined via a text file, "input_tanimoto_sampling.txt", with the following format:

```
Inputs
Header of SMILES column
Header of ID column
Number of lists (n)
Target list size
Input filename (with extension) 1
Input filename (with extension) 2
...
Output filename (with extension) n
Output filename (with extension) 1
Output filename (with extension) 2
...
Output filename (with extension) n
```

For example:

```
Inputs
smiles
compound_id
3
3000
amides.csv
carboxylic_acids.csv
sulfonamides.csv
amides_diverse.csv
carboxylic_acids_diverse.csv
sulfonamides_diverse.csv
```

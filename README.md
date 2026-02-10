# Pruning complex networks in psychiatric symptomatology
**M.Sc. Thesis in Mathematical Statistics | Stockholm University (2023)**

## Advisors
* **Dr. Tom Britton**, Stockholm University
* **Dr. Lars Klintwall**, Karolinska Institutet

## Overview
This repository contains the data and R implementation for my Master’s thesis. The project investigates edge pruning methods for complex (large) networks describing self-reported interactions between psychiatric symptoms for psychiatric patients. The purpose is to simplify the networks while maintaining the most influential interactions to be able to better evaluate the patient's mental health, construct treatment plans targeting specific interactions, and share a more interpretable visual representation with the patient.

* **Pruning Methods:** Edge Betweenness, PageRank, Updated PageRank, and a Connectivity-based "Brute Force" approach.
* **Objective:** Simplify complex networks by pruning redundant edges until the networks are more interpretable while maintaining the most significant structures.
* **Mathematics:** * **Centrality-Based Pruning:** Adapted **Edge Betweenness** and **PageRank** algorithms for directed, weighted graphs to iteratively remove edges with the lowest relative influence.
    * **Connectivity Optimization:** Generalized a "brute force" approach for directed networks that allows for disconnected components, pruning edges that contribute the least to the global connectivity in each step.
    * **Reliability Testing:** Used synthetic networks with added noise to evaluate the methods' ability to recover ground-truth structures, alongside similarity measures for empirical network pairs.
* **Results:** * All methods significantly improved the similarity between synthetic networks and pruned "noisy" networks, particularly when noise consisted of many small-weight edges.
    * **Expert Evaluation:** Psychologists overall preferred the Connectivity-based method for visualization; however, those with direct knowledge of the specific patients preferred **Edge Betweenness** and **PageRank** for maintaining clinical nuance.
    * There is a trade-off between practical usefulness, statistical accuracy, and computational time complexity when selecting an optimal pruning strategy.

## Contents
* `2023_9_report.pdf`: The full thesis text (English).
* `/code`: R script for the pruning methods and data visualization.
* `/data`: Data from human expert evaluations (psychologists) and anonymized questionnaire data on symptom interactions.

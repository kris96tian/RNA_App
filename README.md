# RNA Secondary Structure Predictor

## Visit Website:
- https://rna-secondary.streamlit.app/

<img src="blob:chrome-untrusted://media-app/7cdcde57-db14-41ef-bd44-69f9502a0617" />![image](https://github.com/user-attachments/assets/c512742a-50bc-4a82-9cb0-15f91413680b)


Interactive tool for the prediction and visualization of RNA secondary structures, using dynamic programming. It predicts base pairing and generates a visual representation of the structure based on the input RNA sequence.

## Features

- **Input:** Users provide an RNA sequence composed of the nucleotides `A`, `U`, `G`, and `C` (e.g., `AUGCUAGCUAGC`).
- **Prediction:** The program uses a dynamic programming approach to predict valid base pairs and generate a secondary structure in the dot-bracket notation.
- **Visualization:** A Plotly-based interactive figureof the RNA sequence with its predicted secondary structure.

## Components

- **Dynamic Programming Algorithm:** 
    - A scoring system determines the optimal RNA folding by evaluating base-pairing possibilities using the `can_pair` function (A-U, G-C, and G-U wobble pairs).
    - The traceback process reconstructs the structure by recording paired and unpaired bases.

- **Visualization:** 
    - RNA structure is plotted with the base-pairing information.
    - Base pairs: color-coded for nucleotide pairs (A-U in red, G-C in green, G-U in blue) and arranged into stems.

## Libraries 

- **Streamlit (st)**
- **Pandas**
- **Plotly**
- **NumPy**

## Workflow

1. **Input:** RNA nuleotide-Sequenz (e.g., `GGGGGUCCGUCGAGUGGUCCGGCGGCUUACUGAAGGGGACCGGUGGGAGGCGGACCGGGAGG`).
2. **Prediction:**  "🔍 Predict Structure" button.
3. **Visualization:** Interactive plot of the RNA structure is shown.
4. **Details:** Also shown: the sequence, dot-bracket notation, and detailed list of base pairs in a table.


## Example RNA Sequence

- **Input:** `AUGCUAGCUAGC`
- **Dot-Bracket Structure:** `(((...)))`
- **Base Pairs:** 
  - Position 1-12: A-U
  - Position 2-11: U-A
  - Position 3-10: G-C
  - Position 4-9: C-G
  - Position 5-8: U-A
  - Position 6-7: G-C

## Created by Kristian Alikaj

For more information or updates, visit [My GitHub](https://github.com/kris96tian) or [My Portfolio Website](https://kris96tian.github.io/).

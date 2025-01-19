# RNA Secondary Structure Predictor

## App
[FLASK (faster) RNA Secondary Structure App](https://kr1571an.pythonanywhere.com/)

[Streamlit RNA Secondary Structure App](https://rna-secondary.streamlit.app/?embed_options=dark_theme)

<img src="blob:chrome-untrusted://media-app/7cdcde57-db14-41ef-bd44-69f9502a0617" />![image](https://github.com/user-attachments/assets/c512742a-50bc-4a82-9cb0-15f91413680b)


Interactive tool for the prediction and visualization of RNA secondary structures, using dynamic programming. It predicts base pairing and generates a visual representation of the structure based on the input RNA sequence.
The algorithm I implemented in this code is related to the **Nussinov’s algorithm**.

---

### **Objective**
Optimal secondary structure prediction of RNA sequence by **maximizing the number of base pairs** under specific constraints.

---

### **Definitions**
1. **RNA Sequence**:
   - A string of nucleotides: adenine (**A**), uracil (**U**), guanine (**G**), and cytosine (**C**).
   - Valid base-pair interactions:
     - Canonical: A-U, U-A, G-C, C-G
     - Wobble: G-U, U-G

2. **Dot-Bracket Notation**:
   - The predicted structure is represented as a string of `(`, `)`, and `.`:
     - `.`: Unpaired base.
     - `(`, `)`: Paired bases, indicating a bond.

3. **Dynamic Programming Approach**:
   - DP table `dp[i][j]` where:
     - `i` and `j` represent indices of the RNA sequence.
     - `dp[i][j]` stores the maximum number of base pairs in the subsequence from `i` to `j`.

---

### **Steps**
1. **DP Table initialization**:
   - A 2D matrix `dp` of size \(n \times n\) (where \(n\) is the sequence length) is initialized to 0.
   - A separate `traceback` matrix is used to store the decisions for reconstructing the structure later.

2. **DP Table iterative filling**:
   - Start with short subsequences and expand to larger ones.
   - For each pair of indices `(i, j)` (where \(j > i\)):
     - **Case 1: Unpaired**: The nucleotide at position `j` is left unpaired:
       \[
       dp[i][j] = dp[i][j-1]
       \]
     - **Case 2: Paired**: If `sequence[i]` can pair with `sequence[j]`:
       - Add 1 (for the new pair) to the solution of the inner subsequence \((i+1, j-1)\):
         \[
         dp[i][j] = \max(dp[i][j], 1 + dp[i+1][j-1])
         \]
     - **Case 3: Bifurcation**: Split the subsequence into two parts:
       - Combine solutions from two non-overlapping subsequences \((i, k)\) and \((k+1, j)\):
         \[
         dp[i][j] = \max(dp[i][j], dp[i][k] + dp[k+1][j]) \quad \text{for } k \in (i + \text{MIN\_LOOP}, j - \text{MIN\_LOOP})
         \]

   **Constraints**:
   - Min-loop-size: Ensures there are at least `MIN_LOOP` unpaired bases between paired bases \(i\) and \(j\).
   - Valid base pairing: Ensures following of biological rules.

3. **Traceback**:
   - For reconstructing the optimal solution.
   - Starting from the full sequence `(0, n-1)`, decisions stored in `traceback[i][j]` determine:
     - Whether `i` and `j` are paired.
     - If the subsequence bifurcates into two smaller problems.

4. **Output**:
     - **Dot-bracket notation** for the secondary structure.
     - A list of base pairs (i.e., indices of paired bases).

---

### **Mathematical Recurrence**
The algorithm can be described using the recurrence relation:

![equation](https://latex.codecogs.com/gif.latex?dp%5Bi%5D%5Bj%5D%20%3D%20%5Cmax%5Cbegin%7Bcases%7D%20dp%5Bi%5D%5Bj-1%5D%20%26%20%5Ctext%7B%28Unpaired%29%7D%20%5C%5C%20dp%5Bi&plus;1%5D%5Bj-1%5D%20&plus;%201%20%26%20%5Ctext%7B%28Paired%2C%20if%20%7D%20%5Ctext%7Bcan%5C_pair%28sequence%5Bi%5D%2C%20sequence%5Bj%5D%29%7D%20%5Ctext%7B%29%7D%20%5C%5C%20%5Cmax_%7Bk%3Di&plus;%5Ctext%7BMIN%5C_LOOP%7D%7D%5E%7Bj-%5Ctext%7BMIN%5C_LOOP%7D%7D%20%28dp%5Bi%5D%5Bk%5D%20&plus;%20dp%5Bk&plus;1%5D%5Bj%5D%29%20%26%20%5Ctext%7B%28Bifurcation%29%7D%20%5Cend%7Bcases%7D)
---

### **Time Complexity**
- **Table Filling**:
  - The outer loop iterates over subseq. lengths \(L\) from `MIN_LOOP + 2` to \(n\), and inner loops iterate over \(i\) and \(j\).
  - For each cell `(i, j)`, there is an additional loop over possible bifurcation points \(k\).
  - **Complexity**: \(O(n^3)\).

- **Traceback**:
  - Linear in the number of base pairs.
  - **Complexity**: \(O(n)\).

Overall time complexity is \(O(n^3)\).

---

### **Biological Significance**
Simplified model for RNA secondary structure prediction:
- It **maximizes base pairing** without considering thermodynamic stability or pseudoknots.
- Extensions like Zuker’s algorithm can incorporate free-energy minimization for more accurate predictions.

---

### **Example**
For the sequence `AUGCGUA`:
1. `dp[i][j]` will compute the maximum number of pairs.
2. Traceback reconstructs the structure:
   - Dot-bracket: `.(())..`
   - Base pairs: [(1, 6), (2, 5)]

import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from typing import Tuple, List
import numpy as np

#page config
st.set_page_config(
    page_title="RNA Structure Predictor",
    page_icon="🧬",
    layout="wide",
)

# CSS
st.markdown("""
    <style>
    .stApp {
        background: linear-gradient(135deg, #4e73df22, #1cc88a22);
    }
    .main {
        background-color: transparent;
    }
    .block-container {
        padding-top: 2rem;
    }
    h1 {
        color: #333;
        text-align: center;
        padding-bottom: 2rem;
    }
    .stButton button {
        width: 100%;
        background-color: #4e73df;
        color: white;
        padding: 0.5rem;
        border-radius: 0.5rem;
        transition: all 0.3s;
    }
    .stButton button:hover {
        background-color: #1cc88a;
        transform: scale(1.02);
    }
    </style>
    """, unsafe_allow_html=True)

def can_pair(base1: str, base2: str) -> bool:
    """Check if two bases can form a valid pair"""
    valid_pairs = {
        ('A', 'U'), ('U', 'A'),
        ('G', 'C'), ('C', 'G'),
        ('G', 'U'), ('U', 'G')  # Wobble
    }
    return (base1, base2) in valid_pairs

def predict_rna_structure(sequence: str) -> Tuple[str, List[Tuple[int, int]]]:
    n = len(sequence)
    MIN_LOOP = 3
    
    dp = [[0] * n for _ in range(n)]
    traceback = [[None] * n for _ in range(n)]
    
    for length in range(MIN_LOOP + 2, n + 1):
        for i in range(n - length + 1):
            j = i + length - 1
            
            dp[i][j] = dp[i][j-1]
            traceback[i][j] = ('unpaired', j)
            
            if can_pair(sequence[i], sequence[j]):
                paired_value = (1 + dp[i+1][j-1] if i+1 < j-1 else 1)
                if paired_value > dp[i][j]:
                    dp[i][j] = paired_value
                    traceback[i][j] = ('paired', i, j)
            
            for k in range(i + MIN_LOOP + 1, j - MIN_LOOP):
                if dp[i][k] + dp[k+1][j] > dp[i][j]:
                    dp[i][j] = dp[i][k] + dp[k+1][j]
                    traceback[i][j] = ('bifurcation', k)
    
    structure = ['.' for _ in range(n)]
    base_pairs = []
    
    def traceback_structure(i, j):
        if i >= j: return
        if traceback[i][j] is None: return
            
        move_type = traceback[i][j][0]
        if move_type == 'paired':
            structure[i] = '('
            structure[j] = ')'
            base_pairs.append((i, j))
            if i+1 < j-1:
                traceback_structure(i+1, j-1)
        elif move_type == 'bifurcation':
            k = traceback[i][j][1]
            traceback_structure(i, k)
            traceback_structure(k+1, j)
        elif move_type == 'unpaired':
            traceback_structure(i, j-1)
    
    traceback_structure(0, n-1)
    return ''.join(structure), base_pairs

def visualize_structure(sequence: str, base_pairs: List[Tuple[int, int]]) -> go.Figure:
    n = len(sequence)
    
    x = list(range(n))
    y = [0] * n
    stems = []
    current_stem = []
    
    sorted_pairs = sorted(base_pairs)
    for i, j in sorted_pairs:
        if not current_stem:
            current_stem = [(i, j)]
        else:
            prev_i, prev_j = current_stem[-1]
            if i == prev_i + 1 and j == prev_j - 1:
                current_stem.append((i, j))
            else:
                if len(current_stem) > 0:
                    stems.append(current_stem)
                current_stem = [(i, j)]
    if current_stem:
        stems.append(current_stem)
    for stem_idx, stem in enumerate(stems):
        direction = 1 if stem_idx % 2 == 0 else -1
        for pair_idx, (i, j) in enumerate(stem):
            offset = (pair_idx + 1) * 0.4 * direction
            y[i] = offset
            y[j] = offset
            
            mid_x = (x[i] + x[j]) / 2
            x[i] = mid_x - 0.5
            x[j] = mid_x + 0.5
    
    # base nodes
    node_trace = go.Scatter(
        x=x, y=y,
        mode='markers+text',
        marker=dict(
            size=30,
            color=['#ff7f7f' if b == 'A' else '#7fbf7f' if b == 'U' else 
                   '#7f7fff' if b == 'G' else '#ffbf7f' for b in sequence],
            line=dict(width=1, color='black')
        ),
        text=list(sequence),
        textposition='middle center',
        hoverinfo='text',
        textfont=dict(size=14, color='white')
    )
    
    # backbone connections
    backbone_x = []
    backbone_y = []
    for i in range(n-1):
        backbone_x.extend([x[i], x[i+1], None])
        backbone_y.extend([y[i], y[i+1], None])
    
    backbone_trace = go.Scatter(
        x=backbone_x, y=backbone_y,
        mode='lines',
        line=dict(color='#ccc', width=1),
        hoverinfo='none'
    )
    
    # base pair connections
    pair_traces = []
    for i, j in base_pairs:
        pair_traces.append(go.Scatter(
            x=[x[i], x[j]], y=[y[i], y[j]],
            mode='lines',
            line=dict(
                color='#ff9999' if (sequence[i], sequence[j]) in [('A','U'),('U','A')] else
                      '#99ff99' if (sequence[i], sequence[j]) in [('G','C'),('C','G')] else
                      '#9999ff',
                width=1
            ),
            hoverinfo='none'
        ))
    
    # figure
    fig = go.Figure(data=[backbone_trace, node_trace] + pair_traces)
    fig.update_layout(
        showlegend=False,
        plot_bgcolor='white',
        width=800,
        height=400,
        margin=dict(t=20, b=20, l=20, r=20),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
    )
    
    return fig

# Main
st.title("🧬 RNA Secondary Structure Predictor")
st.markdown(
    "#### Enter RNA sequence (A, U, G, C only)\n"
    "**Example:** `GGGGGUCCGUCGAGUGGUCCGGCGGCUUACUGAAGGGGACCGGUGGGAGGCGGACCGGGAGG`"
)
sequence = st.text_area(
    "",
    height=100,
    placeholder="Example: AUGCUAGCUAGC",
    help="Enter a valid RNA sequence using only A, U, G, and C nucleotides.",
)


if st.button("🔍 Predict Structure"):
    if not sequence:
        st.error("Please enter an RNA sequence")
    elif not all(base in 'AUGC' for base in sequence.upper()):
        st.error("Invalid sequence. Use only A, U, G, C")
    else:
        sequence = sequence.upper()
        with st.spinner("Predicting structure..."):
            dot_bracket, base_pairs = predict_rna_structure(sequence)
            st.success("Structure predicted successfully!")
            st.markdown("<div style='margin-top: 20px; text-align: center;'>", unsafe_allow_html=True)
            st.markdown("### RNA Secondary Structure Visualization")
            fig = visualize_structure(sequence, base_pairs)
            st.plotly_chart(fig, use_container_width=True, height=800)
            st.markdown("### Detailed View")
            st.markdown("**Sequence and Structure:**")
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("**Sequence:**")
                st.code(sequence)
            with col2:
                st.markdown("**Structure:**")
                st.code(dot_bracket)
            st.markdown("**Base Pairs:**")
            pairs_df = pd.DataFrame(
                base_pairs,
                columns=['Position 1', 'Position 2']
            ).assign(
                Pair=lambda x: [f"{sequence[i]}-{sequence[j]}" for i, j in base_pairs]
            )
            st.dataframe(pairs_df, use_container_width=True)

st.markdown("""
        ---
        **Created by Kristian Alikaj**  
        For more, visit [My GitHub](https://github.com/kris96tian) or [My Portfolio Website](https://kris96tian.github.io/)
""")

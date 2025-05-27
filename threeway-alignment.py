from math import log2
import csv
from threeway_align.BLOSUM62 import *

def compute_gap_penalty(indel_rate, P):
    '''
    This function computes the gap penalty for threeway_align using the BLOSUM62 model and the indel rate
    The transition probabilities have been loaded for you and stored in the nested dictionary P (see file threeway_align/BLOSUM62 and threeway_align/blosum62sym.csv ) 
    '''
    q = {}
    total_p = 0
    for aa in AMINOS:
        q[aa] = 0
        for aa2 in AMINOS:
            q[aa] += P[aa][aa2]
            total_p += P[aa][aa2]
    
    for aa in AMINOS:
        q[aa] /= total_p
    
    expected_score = 0
    for aa1 in AMINOS:
        for aa2 in AMINOS:
            p_ij = P[aa1][aa2] / total_p
            q_i_q_j = q[aa1] * q[aa2]
            
            if q_i_q_j > 0 and p_ij > 0:
                log_odds = 2 * log2(p_ij / q_i_q_j)
                expected_score += p_ij * log_odds
    
    scaling_factor = -1.0 - 0.5 * (0.1 / max(indel_rate, 0.001))
    gap_penalty = expected_score * scaling_factor
    gap = max(min(gap_penalty, -0.5), -20)
    return gap

def threeway_align(s1, s2, s3, B, gap, VERBOSE=False):
    '''
    This function computes the threeway sequence alignment between strings s1, s2, and s3 given a substitution matrix and gap penalty
    :param s1 (string): the first sequence
    :param s2 (string): the second sequence
    :param s3 (string): the third sequence
    :param B (char,char -> int) : the substitution matrix (e.g. BLOSUM62)
    :param gap (negative float):  gap penalty
    '''
    # initialize (S[i][j][k] = (score,arrow))
    if VERBOSE:
        from sys import stderr; print("Initializing cube", file=stderr)
    
    m, n, p = len(s1), len(s2), len(s3)
    
    S = {}
    ptr = {}
    
    S[(0,0,0)] = 0
    
    # fill in cube axes
    if VERBOSE:
        print("Filling cube axes", file=stderr)
    
    gap_i = [i * gap for i in range(max(m, n, p) + 1)]
    
    for i in range(1, m+1):
        S[(i,0,0)] = gap_i[i]
        ptr[(i,0,0)] = (i-1, 0, 0)
    
    for j in range(1, n+1): 
        S[(0,j,0)] = gap_i[j]
        ptr[(0,j,0)] = (0, j-1, 0)
        
    for k in range(1, p+1): 
        S[(0,0,k)] = gap_i[k]
        ptr[(0,0,k)] = (0, 0, k-1)
    
    # fill in cube faces
    if VERBOSE:
        print("Filling cube faces", file=stderr)
    
    two_gaps = 2 * gap
    
    for j in range(1, n+1):
        for k in range(1, p+1):
            s2j = s2[j-1]
            s3k = s3[k-1]
            b_s2j_s3k = B[s2j][s3k]
            
            score1 = S[(0,j-1,k-1)] + b_s2j_s3k
            score2 = S[(0,j-1,k)] + gap
            score3 = S[(0,j,k-1)] + gap
            
            best_score = max(score1, score2, score3)
            S[(0,j,k)] = best_score
            
            if best_score == score1:
                ptr[(0,j,k)] = (0, j-1, k-1)
            elif best_score == score2:
                ptr[(0,j,k)] = (0, j-1, k)
            else:
                ptr[(0,j,k)] = (0, j, k-1)
    
    for i in range(1, m+1):
        for k in range(1, p+1):
            s1i = s1[i-1]
            s3k = s3[k-1]
            b_s1i_s3k = B[s1i][s3k]
            
            score1 = S[(i-1,0,k-1)] + b_s1i_s3k
            score2 = S[(i-1,0,k)] + gap
            score3 = S[(i,0,k-1)] + gap
            
            best_score = max(score1, score2, score3)
            S[(i,0,k)] = best_score
            
            if best_score == score1:
                ptr[(i,0,k)] = (i-1, 0, k-1)
            elif best_score == score2:
                ptr[(i,0,k)] = (i-1, 0, k)
            else:
                ptr[(i,0,k)] = (i, 0, k-1)
    
    for i in range(1, m+1):
        for j in range(1, n+1):
            s1i = s1[i-1]
            s2j = s2[j-1]
            b_s1i_s2j = B[s1i][s2j]
            
            score1 = S[(i-1,j-1,0)] + b_s1i_s2j
            score2 = S[(i-1,j,0)] + gap
            score3 = S[(i,j-1,0)] + gap
            
            best_score = max(score1, score2, score3)
            S[(i,j,0)] = best_score
            
            if best_score == score1:
                ptr[(i,j,0)] = (i-1, j-1, 0)
            elif best_score == score2:
                ptr[(i,j,0)] = (i-1, j, 0)
            else:
                ptr[(i,j,0)] = (i, j-1, 0)
    
    # fill in rest of cube
    if VERBOSE:
        print("Filling rest of cube", file=stderr)
    
    use_sparse = max(m, n, p) > 100
    max_keys_to_keep = 5000000
    neg_inf = float('-inf')
    
    for i in range(1, m+1):
        if VERBOSE and i % 10 == 0:
            print(f"Processing i={i}/{m}", file=stderr)
            
        if use_sparse and len(S) > max_keys_to_keep:
            keys_to_keep = {}
            for key in S:
                ki, kj, kk = key
                if ki >= i-1:
                    keys_to_keep[key] = S[key]
            S = keys_to_keep
            
        s1i = s1[i-1]
        
        for j in range(1, n+1):
            s2j = s2[j-1]
            b_s1i_s2j = B[s1i][s2j]
            
            for k in range(1, p+1):
                s3k = s3[k-1]
                b_s1i_s3k = B[s1i][s3k]
                b_s2j_s3k = B[s2j][s3k]
                
                score1 = S.get((i-1,j-1,k-1), neg_inf) + b_s1i_s2j + b_s1i_s3k + b_s2j_s3k
                score2 = S.get((i-1,j-1,k), neg_inf) + b_s1i_s2j + two_gaps
                score3 = S.get((i-1,j,k-1), neg_inf) + b_s1i_s3k + two_gaps
                score4 = S.get((i,j-1,k-1), neg_inf) + b_s2j_s3k + two_gaps
                score5 = S.get((i-1,j,k), neg_inf) + two_gaps
                score6 = S.get((i,j-1,k), neg_inf) + two_gaps
                score7 = S.get((i,j,k-1), neg_inf) + two_gaps
                
                best_score = score1
                best_idx = 0
                
                if score2 > best_score:
                    best_score = score2
                    best_idx = 1
                if score3 > best_score:
                    best_score = score3
                    best_idx = 2
                if score4 > best_score:
                    best_score = score4
                    best_idx = 3
                if score5 > best_score:
                    best_score = score5
                    best_idx = 4
                if score6 > best_score:
                    best_score = score6
                    best_idx = 5
                if score7 > best_score:
                    best_score = score7
                    best_idx = 6
                
                S[(i,j,k)] = best_score
                
                if best_idx == 0:
                    ptr[(i,j,k)] = (i-1, j-1, k-1)
                elif best_idx == 1:
                    ptr[(i,j,k)] = (i-1, j-1, k)
                elif best_idx == 2:
                    ptr[(i,j,k)] = (i-1, j, k-1)
                elif best_idx == 3:
                    ptr[(i,j,k)] = (i, j-1, k-1)
                elif best_idx == 4:
                    ptr[(i,j,k)] = (i-1, j, k)
                elif best_idx == 5:
                    ptr[(i,j,k)] = (i, j-1, k)
                else:
                    ptr[(i,j,k)] = (i, j, k-1)
    
    # backtrack to get alignments
    if VERBOSE:
        print("Backtracking to build alignment", file=stderr)
    
    estimated_aln_length = m + n + p
    aln_s1 = [None] * estimated_aln_length
    aln_s2 = [None] * estimated_aln_length
    aln_s3 = [None] * estimated_aln_length
    
    i, j, k = m, n, p
    idx = 0
    
    while i > 0 or j > 0 or k > 0:
        prev_i, prev_j, prev_k = ptr[(i,j,k)]
        
        if i > prev_i and j > prev_j and k > prev_k:
            aln_s1[idx] = s1[i-1]
            aln_s2[idx] = s2[j-1]
            aln_s3[idx] = s3[k-1]
        elif i > prev_i and j > prev_j:
            aln_s1[idx] = s1[i-1]
            aln_s2[idx] = s2[j-1]
            aln_s3[idx] = '-'
        elif i > prev_i and k > prev_k:
            aln_s1[idx] = s1[i-1]
            aln_s2[idx] = '-'
            aln_s3[idx] = s3[k-1]
        elif j > prev_j and k > prev_k:
            aln_s1[idx] = '-'
            aln_s2[idx] = s2[j-1]
            aln_s3[idx] = s3[k-1]
        elif i > prev_i:
            aln_s1[idx] = s1[i-1]
            aln_s2[idx] = '-'
            aln_s3[idx] = '-'
        elif j > prev_j:
            aln_s1[idx] = '-'
            aln_s2[idx] = s2[j-1]
            aln_s3[idx] = '-'
        elif k > prev_k:
            aln_s1[idx] = '-'
            aln_s2[idx] = '-'
            aln_s3[idx] = s3[k-1]
        
        i, j, k = prev_i, prev_j, prev_k
        idx += 1
    
    aln_s1 = ''.join(aln_s1[:idx][::-1])
    aln_s2 = ''.join(aln_s2[:idx][::-1])
    aln_s3 = ''.join(aln_s3[:idx][::-1])
    
    score = S[(m,n,p)]
    return aln_s1, aln_s2, aln_s3, score

def threeway_align_indel_rate(s1, s2, s3, B, indel_rate, VERBOSE=False):
    """
    Optimized version that computes the appropriate gap penalty based on indel rate
    and then performs the alignment.
    """
    gap = compute_gap_penalty(indel_rate,P)                            
    a1,a2,a3,score = threeway_align(s1, s2, s3, B, gap, VERBOSE=VERBOSE)
    return a1,a2,a3,score 


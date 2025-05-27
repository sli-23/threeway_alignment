from threeway_align.BLOSUM62 import *

def read_FASTA(stream):
    '''
    This function reads a FASTA file from a given stream and returns a dictionary mapping identifiers to sequences
    '''
    seqs = {}; name = None; seq = ''
    for line in stream:
        l = line.strip()
        if len(l) == 0:
            continue
        if l[0] == '>':
            if name is not None:
                assert len(seq) != 0, "Malformed FASTA"
                seqs[name] = seq
            name = l[1:]
            assert name not in seqs, "Duplicate sequence ID: %s" % name
            seq = ''
        else:
            seq += l
    assert name is not None and len(seq) != 0, "Malformed FASTA"
    seqs[name] = seq
    return seqs
    
def sum_pairwise_score(s1,s2,s3,gap,B=B,gap_symbol='-'):
    # compute sum pairwise score of an alignment of 3 sequences
    # assumming that the sequences are already aligned
    assert len(s1) == len(s2) == len(s3), "Sequences are not properly aligned!"
    score = 0
    two_gap = 2*gap
    for c1,c2,c3 in zip(s1,s2,s3):
        a1,a2,a3 = sorted([c1,c2,c3])
        if a1 == gap_symbol:
            if a2 == gap_symbol:
                if a3 in B:
                    score += two_gap
                else:
                    assert a3 == gap_symbol, "Found unrecognized character " + a3
            else:
                assert a2 in B, "Found unrecognized character " + a2
                assert a3 in B, "Found unrecognized character " + a3
                score += two_gap + B[a2][a3]
        else:
            assert a1 in B, "Found unrecognized character " + a1
            assert a2 in B, "Found unrecognized character " + a2
            assert a3 in B, "Found unrecognized character " + a3
            score += B[a1][a2] + B[a1][a3] + B[a2][a3]        
    return score                    

import pandas as pd
import numpy as np
from utils import score, removed
#####SW+blast type aligner/searching algorithm#####

#create kmer indexes of all database sequences
def create_index(database, id=None, k=7):
    index = {}
    lists = enumerate(database) if id==None else zip(id, database)
    for i, sequence in lists:
        for j in range(len(sequence) - k + 1):
            kmer = sequence[j:j+k]
            if kmer not in index:
                index[kmer] = []
            index[kmer].append((i, j))
        
    #dictionary:  {kmer: [ ([genome#, >line, segment, year], sequence pos)*n ]} where n is number of sequences with said kmer
    return index

#smith waterson algorithm matrix creation of query + dbseq's match window
def local_alignment(seq1, seq2, bnmatr, convertd):
    m, n = len(seq1), len(seq2)
    score_matrix = np.zeros((m+1, n+1))
    for i in range(1, m+1):
        for j in range(1, n+1):
            match = score_matrix[i-1,j-1] + score(seq1[i-1], seq2[j-1], bnmatr, convertd)
            delete = score_matrix[i-1,j] - 6
            insert = score_matrix[i,j-1] - 6
            score_matrix[i,j] = max(match, delete, insert)
    #np array where col is genomeseq and rows is query
    return score_matrix

#traceback of built SW matrix, instead finds highest score in last letter of query
def traceback(score_matrix, seq1, seq2, i, j, bnmatr, convertd):
    alignment1, alignment2 = '', ''
    adjustment = 0
    while i > 0 and j > 0:
        score_current = score_matrix[i, j]
        score_diagonal = score_matrix[i-1, j-1]
        score_up = score_matrix[i, j-1]
        score_left = score_matrix[i-1, j]

        if score_current == score_diagonal + score(seq1[i-1], seq2[j-1], bnmatr, convertd):
            alignment1 += seq1[i-1]
            alignment2 += seq2[j-1]
            i -= 1
            j -= 1
        #query insertion, gap score is -6 
        elif score_current == score_left - 6:
            alignment1 += seq1[i-1]
            alignment2 += '-'
            adjustment+=1
            i -= 1
        #database insertion
        else: #score_current == score_up - 6:
            alignment1 += '-'
            alignment2 += seq2[j-1]
            adjustment +=1
            j -= 1   

    #adjustment since for every insertion, the calculated position window is shifted down once (ensures accurate db position)
    return alignment1[::-1], alignment2[::-1], adjustment

#run, determine kmer match db window
def search(query, database, index, k, namedata, bnmatr, convertd):
    query = query.strip()

    results = []  
    #create kmers of the query seq
    for i in range(len(query) - k + 1):
        kmer = query[i:i+k]
        #if match is found in index
        if kmer in index: 
            #check all genome entries of that kmer   
            for seq_id, position in index[kmer]:
                if type(seq_id)!=list:
                    seq_id = [seq_id, seq_id, "unknown"]
                db_seq = database[seq_id[0]-1]
                #pick window around match 
                start = max(0, position - len(query)+k)
                end = min(len(db_seq), position + len(query))
                #create matrix of query and window
                score_matrix = local_alignment(query, db_seq[start:end], bnmatr, convertd)
                
                #reverse score matrix to pick the last occurence of highest match, which is usually the seq instead of insertion
                reva = np.flip(score_matrix, 1)
                temp = np.shape(score_matrix)[1]-1-np.argmax(reva, axis = 1)[-1]
                #coordinate of max position
                max_position = (len(query), temp)
                #max position is therefore the max possible score
                max_score = score_matrix[max_position[0]][max_position[1]]

                #cutoff threshold T, being 4 mismatches (filters trash alignment, saves lots of time)
                if(max_score >= int(len(query)*2-20)):
                    #perform traceback
                    alignment1, alignment2, adj = traceback(score_matrix, query, db_seq[start:end], *max_position, bnmatr, convertd)
                    results.append({
                        'score': max_score,
                        'query_alignment': alignment1,
                        'db_alignment': alignment2,
                        'db_sequence_id': seq_id[0],
                        'db_sequence_name': seq_id[1],
                        'db_sequence_year': seq_id[3],
                        'db_sequence_id_year': f"#{seq_id[0]} ({seq_id[3]})",
                        'db_segment': seq_id[2],
                        'start_position': start + max_position[1] - len(alignment2) +adj, 
                        'end_position': start + max_position[1]+adj,
                        'id': namedata[0],
                        'type': namedata[1]
                    })
    #formatted like above
    return results  

#main function to run algorithm
def blastn(database, query, index, k, bnmatr, convertd, namedata=None):
    if k > len(query):
        raise Exception("Kmer length must be shorter than query length (default k=7)")
    
    namedata = namedata or ["My Primer", "Unknown"]
    results = search(query, database, index, k, namedata, bnmatr, convertd)

    pdata = pd.DataFrame(results)
    
    #no entries
    if pdata.empty:
       
        s = "No sequences found above score threshold for "+namedata[0]
        nomatch = ({'score': -1,
            'query_alignment': "",
            'db_alignment': "",
            'db_sequence_id': s,
            'db_sequence_name': "",
            'db_sequence_year': "",
            'db_sequence_id_year': "",
            'db_segment': "",
            'start_position': -1,
            'end_position': -1,
            'id': namedata[0],
            'type': namedata[1]})
        top = pd.DataFrame(nomatch, index=[0])
        print("No alignment found above threshold score")
         
    else: #remove duplicate alignments, remove insertions 
        sortp = pdata.drop_duplicates(subset=['start_position', 'db_sequence_id', 'db_alignment']) 
        top = sortp[sortp['query_alignment'].map(removed).map(len)==len(query)] 
        top = top[top['db_alignment'].map(removed).map(len)==len(query)]
        print("Alignment found")

        if top.isempty:
            s = "Insertion or deletion in primer region for "+namedata[0]
            nomatch = ({'score': -1,
                        'query_alignment': "",
                        'db_alignment': "",
                        'db_sequence_id': s,
                        'db_sequence_name': "",
                        'db_sequence_year': "",
                        'db_sequence_id_year': "",
                        'db_segment': "",
                        'start_position': -1,
                        'end_position': -1,
                        'id': namedata[0],
                        'type': namedata[1]})
            top = pd.DataFrame(nomatch, index=[0])
            print("Indel present")
    
    return top 
    

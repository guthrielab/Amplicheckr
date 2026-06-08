import pandas as pd
from utils import score,reverse, complement, alignvis, truncate
########Post Algorithm primer rating ########


#for alignment only mode
def alignout(results, bnmatr, bdict, namedata=None):
    namedata = namedata or ["My Primer", "Unknown"]
    top = results
    for index, result in top.iterrows():
        align = alignvis(result['query_alignment'], result['db_alignment'], bnmatr)
        print(f"""\n\tPrimer: {namedata[0]}  (Type: {namedata[1]})
            Score: {result['score']}
            Query Alignment:    {reverse(complement(result['query_alignment'], bdict)) if namedata[1] =='R' else result['query_alignment']}
                                {reverse(align) if namedata[1]=='R' else align}
            Database Alignment: {reverse(complement(result['db_alignment'], bdict)) if namedata[1] =='R' else result['db_alignment']}
            Database Name: {result['db_sequence_name']}
            Database Sequence ID: {result['db_sequence_id']}
            Position: {int(result['start_position']), int(result['end_position'])}\n""")

#post algorithm primer scoring, non align only mode
def primerscore(alignmentset, metadb, bdict, mmmatr, convertd, html):
    metadf = pd.DataFrame(metadb, columns=["db_sequence_id","name","segment", "year"])
    metadf.set_index("db_sequence_id",)
   
    warnings = []
    #seperate full per set into individual primer types (name, id)
    fwd = (alignmentset[alignmentset["type"]=="F"]).copy()
    fwdid = (fwd['id']).iloc[0]
    fwdseq = (fwd['query_alignment']).iloc[0]
    rev = (alignmentset[alignmentset["type"]=="R"]).copy()
    revid = rev['id'].iloc[0]
    revseq = rev['query_alignment'].iloc[0]
    pro = (alignmentset[alignmentset["type"]=="P"]).copy()
    proid = pro['id'].iloc[0]
    proseq = pro[ 'query_alignment'].iloc[0]

    setnum = fwd['set'].iloc[0]

#situation where one primer has no entries as in no matches
    if -1 in alignmentset.score.values:
        if -1 in fwd.score.values:
            warnings.append(f"Forward primer {fwdid} has no alignments above \
threshold score for primer set #{setnum}")
           
        if -1 in rev.score.values:
            warnings.append(f"Reverse primer {revid} has no alignments above \
threshold score for primer set #{setnum}")
           
        if -1 in pro.score.values:
            warnings.append(f"Probe {proid} has no alignments above \
threshold score for primer set #{setnum}")
          
        return "",warnings, "" if not html else {"entry": False,"name": f"Primer set #{setnum}"}
    
    
#situation where a primer binds multiple times to the same sequence
    otfwd = (fwd.loc[fwd.duplicated(subset=["db_sequence_id"])]).drop_duplicates(subset="db_sequence_id")
    otrev = (rev.loc[rev.duplicated(subset=["db_sequence_id"])]).drop_duplicates(subset="db_sequence_id")
    otpro = (pro.loc[pro.duplicated(subset=["db_sequence_id"])]).drop_duplicates(subset="db_sequence_id")
    
    if not otfwd.empty:
        warnings.append(f"Forward primer {fwdid} has multiple binding \
locations on genome sequence(s) {', '.join(truncate(str(i)) for i in list(otfwd.db_sequence_name.values))}")
        fwd.drop_duplicates(subset="db_sequence_id", inplace=True)

    if not otrev.empty:
        warnings.append(f"Reverse primer {revid} has multiple binding \
locations on genome sequence(s) {', '.join(truncate(str(i)) for i in list(otrev.db_sequence_name.values))}")
        rev.drop_duplicates(subset="db_sequence_id", inplace=True)
    if not otpro.empty:
        warnings.append(f"Probe {proid} has multiple binding \
locations on genome sequence(s) {', '.join(truncate(str(i)) for i in list(otpro.db_sequence_name.values))}")
        pro.drop_duplicates(subset="db_sequence_id", inplace=True)

#situation where primers do not bind to a certain sequences
    fullset = set(alignmentset["db_sequence_id"])
    fwdset = set(fwd["db_sequence_id"])
    revset = set(rev["db_sequence_id"])
    proset = set(pro["db_sequence_id"])

    failprimers = (fwdset|revset|proset) - (fwdset&revset&proset)
    faildf = (alignmentset.loc[alignmentset["db_sequence_id"].isin(failprimers)]).drop_duplicates(subset="db_sequence_id")
    
    
    nbfwd = fullset-fwdset
    nbrev = fullset-revset
    nbpro = fullset-proset

        
    if len(nbfwd)!=0:
        nbfwddf = metadf.loc[metadf["db_sequence_id"].isin(nbfwd)].copy()
        nbfwddf["report"] = nbfwddf.apply(lambda x: f"{truncate(x['name'])} (segment {x['segment']}, year {x['year']})", axis=1)
        warnings.append(f"Forward primer {fwdid} failed to bind to \
genome sequence(s) {', '.join(list(nbfwddf['report']))}")
        
    if len(nbrev)!=0:
        nbrevdf = metadf.loc[metadf["db_sequence_id"].isin(nbrev)].copy()
        nbrevdf["report"] = nbrevdf.apply(lambda x: f"{truncate(x['name'])} (segment {x['segment']}, year {x['year']})", axis=1)
        warnings.append(f"Reverse primer {revid} failed to bind to \
genome sequence(s) {', '.join(list(nbrevdf['report']))}")

    if len(nbpro)!=0:
        nbprodf = metadf.loc[metadf["db_sequence_id"].isin(nbpro)].copy()
        nbprodf["report"] = nbprodf.apply(lambda x: f"{truncate(x['name'])} (segment {x['segment']}, year {x['year']})", axis=1)
        warnings.append(f"Probe {proid} failed to bind to \
genome sequence(s) {', '.join(list(nbprodf['report']))}")
    

    fwd.set_index("db_sequence_id",inplace=True,drop=False)
    rev.set_index("db_sequence_id",inplace=True,drop=False)
    pro.set_index("db_sequence_id",inplace=True,drop=False)
  
    fwdmatch = len(fwdset)
    revmatch = len(revset)
    promatch = len(proset)

    #Verify the correct order exists, of fwd < probe < rev
    success = fwdset & revset & proset
    primersuccess = alignmentset.loc[alignmentset['db_sequence_id'].isin(success)]
    verifyorder = [pd.Series(list(success)), 
                   primersuccess.loc[primersuccess["type"]=="F", ["end_position"]].reset_index(drop=True),
                   primersuccess.loc[primersuccess["type"]=="P", ["start_position"]].reset_index(drop=True),
                   primersuccess.loc[primersuccess["type"]=="P", ["end_position"]].reset_index(drop=True), 
                   primersuccess.loc[primersuccess["type"]=="R", ["start_position"]].reset_index(drop=True)]
    a = pd.concat(verifyorder, axis=1, ignore_index=True)
    a.columns = ["db_sequence_id", "fend", "pstar", "pend", "rstar"]
    a.set_index("db_sequence_id", inplace=True)
    adjdf = a.query("fend<pstar and pstar<pend and pend<rstar")

    #update each primers to only include filtered primers/alignments
    fwd = fwd.loc[fwd["db_sequence_id"].isin(adjdf.index.values)]
    rev = rev.loc[rev["db_sequence_id"].isin(adjdf.index.values)]
    pro = pro.loc[pro["db_sequence_id"].isin(adjdf.index.values)]

    #info count data
    segments = pd.concat([metadf["segment"].value_counts(), fwd["db_segment"].value_counts(), faildf["db_segment"].value_counts()], axis=1)
    segments.columns = ["totalsegment", "matchedsegment","invalidsegment"]
    totalgenomein = len(metadf)
    rawgenomesmatched = max(fwdmatch, revmatch, promatch)
    genomesafterorder = len(adjdf)

    results = alignrepr(fwd, fwdseq, fwdid, rev, revseq, revid, pro, proseq, proid, bdict, mmmatr, convertd, html)

    if not html:
        info = f"""Total genomes parsed: {totalgenomein}
Genome segment composition: 
{segments['totalsegment'].to_string()} 

Raw genomes matched: {rawgenomesmatched}
Valid matches: {genomesafterorder} 

Segment composition of matches:
{segments['matchedsegment'].to_string()}
Segment composition of invalidated matches:
{segments['invalidsegment'].to_string()}
"""    
        return info, warnings, results
    
       
    else:
        primerdata = {
            "entry": True,
            "name": f"Primer set #{setnum}",
            "total_genomes_parsed": totalgenomein,
            "genome_segment_composition": segments['totalsegment'].to_dict(),
            "raw_genomes_matched": rawgenomesmatched,
            "valid_matches": genomesafterorder,
            "match_segment_composition": segments['matchedsegment'].to_dict(),
            "invalid_match_segment_composition": segments['invalidsegment'].to_dict(),

            "failures": [],
        
        }
        for i, j in enumerate(warnings): primerdata["failures"].append({"id": "#"+str(i+1), "message": j})
        primerdata["primers"] = results

        return primerdata
    

#create output for filtered primers, one entry per alignment rather than seq to save space
def alignrepr(fwd, fwdseq, fwdid, rev, revseq, revid, pro, proseq, proid, bdict, mmmatr, convertd, html):
    results = []

#fwd
    fwdlist = []
    fwdg = fwd.drop_duplicates(subset=["db_alignment"])
    fwdc = fwd["db_alignment"].value_counts()
    fwdrepr = fwdg["db_alignment"]
    for i in list(fwdrepr):
        fwdtemp = mmscore(fwdseq, complement(i, bdict), -5, len(i), mmmatr, convertd)
        fwdtemp.append(i)
        fwdlist.append(fwdtemp)
    fwddf = pd.DataFrame(fwdlist, columns = ["visual", "rating", "seq"]).set_index("seq")
    finalfwd = pd.concat([fwdg.set_index("db_alignment", drop=False), fwdc, fwddf], axis=1)
    fwdsum = finalfwd["count"].sum()
    finalfwd["report"] = finalfwd.apply(lambda x: primerreport(x, fwdid, "F",fwdsum, fwd, bdict, html), axis=1)
    for i in list(finalfwd["report"]):
        results.append(i)

#rev
    revg = rev.drop_duplicates(subset=["db_alignment"])
    revc = rev["db_alignment"].value_counts()
    revrepr = revg["db_alignment"]
    revlist = []
    for i in list(revrepr):
        revtemp = mmscore(complement(revseq, bdict), i, 0, 5, mmmatr,convertd)
        revtemp[0] = reverse(revtemp[0])
        revtemp.append(i)
        revlist.append(revtemp)
    revdf = pd.DataFrame(revlist, columns = ["visual", "rating", "seq"]).set_index("seq")
    finalrev = pd.concat([revg.set_index("db_alignment", drop=False), revc, revdf], axis=1)
    revsum = finalrev["count"].sum()
    finalrev["report"] = finalrev.apply(lambda x: primerreport(x, revid, "R", revsum, rev, bdict, html), axis=1)
    for i in list(finalrev["report"]):
        results.append(i)

#probe
    prog = pro.drop_duplicates(["db_alignment"])
    proc = pro["db_alignment"].value_counts()
    prorepr = prog["db_alignment"]
    prolist = []
    for i in list(prorepr):
        protemp = mmscore(proseq, complement(i, bdict), 5, -5, mmmatr, convertd)
        protemp.append(i)
        prolist.append(protemp)
    prodf = pd.DataFrame(prolist, columns = ["visual", "rating", "seq"]).set_index("seq")
    finalpro = pd.concat([prog.set_index("db_alignment", drop=False), proc, prodf], axis=1)
    prosum = finalpro["count"].sum()
    finalpro["report"] = finalpro.apply(lambda x: primerreport(x, proid, "P", prosum, pro, bdict, html), axis = 1)
    for i in list(finalpro["report"]):
        results.append(i)
    return results
    
#create the string report for each alignment
def primerreport(row, name, type, sum, df, bdict, html):
    if not html:
        s = f"""\nName: {name} (Type: {type}) 
Alignment Score: {row['score']}
Primer Rating: {row['rating']}
Query Alignment:    {reverse(complement(row['query_alignment'], bdict)) if type =="R" else row['query_alignment']}
                    {row['visual']}
Database Alignment: {reverse(complement(row['db_alignment'], bdict)) if type =="R" else row["db_alignment"]}
Database Sequence Prevalence: {row['count']} ({(row['count']/(sum))*100}%)
Database Matches: {', '.join(truncate(str(i)) for i in ((df.loc[df['db_alignment']==row['db_alignment']])['db_sequence_name'].to_list()))}
"""
    else: 
        s = {
            "name": name, 
            "type": type, 
            "alignment_score": row['score'],
            "primer_rating": row['rating'],
            "query_alignment": reverse(complement(row['query_alignment'], bdict)) if type =="R" else row['query_alignment'],
            "visual": row['visual'],
            "database_alignment": reverse(complement(row['db_alignment'], bdict)) if type =="R" else row["db_alignment"],
            "prevalence": f"{row['count']} ({(row['count']/(sum))*100}%)",
            "matches": ', '.join(truncate(str(i)) for i in ((df.loc[df['db_alignment']==row['db_alignment']])['db_sequence_name'].to_list()))
        }

    return s


#scoring based on mismatch type, 3 mismatches, two purpur,or one mm in critical region = high risk
# anything more than 1 lowest mm in non-critical =moderate risk
def mmscore(seq1, seq2, cstar, cend, mmmatr, convertd):
    vislist = []
    category= ""
    mismatch = 0
    for i, j in zip(seq1[:cstar], seq2[:cstar]):
        ntscore = score(i, j, mmmatr, convertd)
        match ntscore:
            case 2:
                vislist.append("|")
            case -3:
                mismatch += 1
                vislist.append("X")
            case -4:
                mismatch +=1.25
                vislist.append("x")
            case -5:
                mismatch +=1.5
                vislist.append("-")
    for i, j in zip(seq1[cstar:cend], seq2[cstar:cend]):
        ntscore = score(i, j, mmmatr, convertd)
        match ntscore:
            case 2:
                vislist.append("|")
            case -3:
                mismatch += 3
                vislist.append("X")
            case -4:
                mismatch +=3
                vislist.append("x")
            case-5:
                mismatch +=3 
                vislist.append("-")
    for i, j in zip(seq1[cend:], seq2[cend:]):
        ntscore = score(i, j, mmmatr, convertd)
        match ntscore:
            case 2:
                vislist.append("|")
            case -3:
                mismatch += 1
                vislist.append("X")
            case -4:
                mismatch +=1.25
                vislist.append("x")
            case-5:
                mismatch +=1.5
                vislist.append("-")
    if mismatch>=3:
        category = "High Risk"
    elif mismatch>1:
        category = "Moderate Risk"
    else: 
        category = "No Risk"
    
    return ["".join(vislist), category]

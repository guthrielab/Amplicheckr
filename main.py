import argparse
import time
from matrices import *
from parsers import *
from aligner import *
from scorer import *

thumbnail = r"""
_______________________________________________________________________________________________________
 ______                      ___                __              __                    _          __     
/\  _  \                    /\_ \    __        /\ \            /\ \                 /' \       /'__`\   
\ \ \ \ \    ___ ___   _____\//\ \  /\_\    ___\ \ \___      __\ \ \/'\   _ __     /\_, \     /\ \/\ \  
 \ \  __ \ /' __` __`\/\ '__`\\ \ \ \/\ \  /'___\ \  _ `\  /'__`\ \ , <  /\`'__\   \/_/\ \    \ \ \ \ \ 
  \ \ \/\ \/\ \/\ \/\ \ \ \ \ \\_\ \_\ \ \/\ \__/\ \ \ \ \/\  __/\ \ \\`\\ \ \/       \ \ \  __\ \ \_\ \
   \ \_\ \_\ \_\ \_\ \_\ \ ,__//\____\\ \_\ \____\\ \_\ \_\ \____\\ \_\ \_\ \_\        \ \_\/\_\\ \____/
    \/_/\/_/\/_/\/_/\/_/\ \ \/ \/____/ \/_/\/____/ \/_/\/_/\/____/ \/_/\/_/\/_/  _______\/_/\/_/ \/___/ 
                         \ \_\                                                  /\______\               
                          \/_/                                                  \/______/               
_______________________________________________________________________________________________________

"""
description = """
Perform \"local\" alignment via modified Smith-Waterson algorithm for 
primers/sequences along genomes to score and validate amplicon efficiency.
""" 

np.set_printoptions(legacy="1.25")

###### main function, run all, output#####
def amplichekr(primerset, qudb, ntdb, index, k, metadb, bdict, convertd, bnmatr, mmmatr, alignonly=False, html=None):
    try:
        count=1
        alignresults = []
        
        for i in primerset:
            tempin = qudb[qudb["set"]==i]
            tempresult = []
            
            for j, row in tempin.iterrows():
                seq = row["sequence"]
                namedata = (row["id"], row["type"])
                
                if not seq.isalpha():
                    raise Exception(f"Primer {namedata[0]} has non valid characters: {seq}")
                if row["type"]=="R":
                    seq = reverse(complement(seq, bdict))
                print("\nAligning Primer #"+str(count))
                count+=1

                align = blastn(ntdb, seq, index, k,  bnmatr, convertd, namedata)
                
                align['set'] = i
                alignout(align, bnmatr, bdict, namedata) if alignonly else tempresult.append(align)
            if not alignonly:
                setdf = pd.concat(tempresult)
                alignresults.append(setdf)
        
        if not html:
            for i, j in enumerate(alignresults):
                infotitle, warnings, report = primerscore(j, metadb, bdict, mmmatr, convertd, html)
            
                print(f"\nPrimer set #{i+1}:")
                print(infotitle) if len(infotitle)!=0 else print("No amplification")
                print("\nDetails on primer failures:") if  len(j)!=0 else print()
                for i, j in enumerate(warnings): print("#"+str(i+1)+": "+j)
                for i in report: print(i)
        else:
            pass
        

        
    except SystemError:
        print("Missing required rows in primer file")
   

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("primers", 
                        help="(Required) File for primers as csv, must include columns \
                        'name', 'target', 'sequence', 'set' (integers, where a set has exactly 1 each of the 3 types of primers), \
                        and 'type' (either as F,R, or P) in any order")
    parser.add_argument("genomes", 
                        help="(Required) File name for genome as fasta. Acts as directory input when using -f.")
    parser.add_argument("-f", "--files",  
                        help="Use multiple input files for genomes by inputing common suffix")
    parser.add_argument("-k","--kmer", type=int, default=7,
                        help="Specify kmer length to use for indexing and searching (default k=7)")
    parser.add_argument("-v","--virus", default= "influenza",
                        help="Target virus species (influenza/i, hantavirus/h)")
    parser.add_argument("-hi", "--hello", action="store_true",
                        help="hi")
    parser.add_argument("-a", "--align", action="store_true",
                        help="Only perform alignment and skip grading.")
    parser.add_argument("-H", "--html", action="store_true",
                        help="Activate html report generator (silences stdout)")
    parser.add_argument("-o", "--output", 
                        help="Make html report instead of stdout, enter file name")
    args = parser.parse_args()

    primerin = args.primers
    genomein = args.genomes
    mode = args.virus
    k = args.kmer
    html = args.html
    if args.hello:
        print("hi")
    
    qudb, primerset = parsequ(fname=primerin)
    metadb, ntdb = parsedb(fname=genomein, fext=args.files, mode=mode)
   
    
    bnmatr = bnmatrix()
    convertd = convertdict()
    bdict = basedict()
    mmmatr = mmmatrix()

    index = create_index(ntdb, metadb, k)
    
    t0 = time.time_ns()
    
    
    amplichekr(primerset, qudb, ntdb, index, k, metadb,bdict, convertd, bnmatr, mmmatr, args.align, html)

    tf = time.time_ns()
    print (f"finished in {(tf-t0)//1000000}ms")

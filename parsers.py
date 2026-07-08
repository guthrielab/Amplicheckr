import os
import glob
import pandas as pd
from Bio import SeqIO
import re

#####Parsing functions#####
#parse csv file with primer data
def parsequ(path=os.getcwd(), fname=""):
    try:
        qudb= pd.read_csv(os.path.join(path,fname))
        #replace + join target and name to "id"
        qudb["id"] = \
            qudb.apply(lambda row: str(row["target"]).replace(" ", "-")+"_"+str(row["name"]).replace(" ", "-") \
            +"_"+str(row.name), axis=1)
        qudb.drop(["target","name"], axis =1, inplace=True)

        #list of possible primer sets present
        primersets = (qudb["set"].unique())
    except FileNotFoundError or OSError:
        print(f"Primer file error: check path and filename ({os.path.join(path,fname)})")
    except KeyError:
        print("Missing required rows in primer file")
    
    return qudb, primersets

#parse file with full genome sequences
def parsedb(path=os.getcwd(), fname="", fext=None, mode = "influenza"):
    metadb = []
    ntdb = []
    #check if its one file, or multifile mode where fext is the extension
    fext= fname if not fext else os.path.join(fname,"*"+fext)
    count=1
    
    try:
        for file in glob.glob(os.path.join(path, fext)): 
            dict = SeqIO.index(file, "fasta")
            if len(dict)==0:
                raise Exception("Genome file error: data missing or bad format")
            #for every fasta, record genome number, check the > line for segments
            for metadata in dict.keys():
                temp = [count, metadata]
                #[genome#, >line, segment, year]
                if mode=="influenza" or mode=="i":
                    segment = re.findall(r'\b(HA|NA|MP|PA|NP|NS|PB1|PB2)\b', dict[metadata].description)
                    temp.append(segment[0] if len(segment)==1 else "unknown")
                elif mode=="hantavirus" or mode=="h":
                    if "segment S" in dict[metadata].description:
                        temp.append("S")
                    elif "segment M" in dict[metadata].description:
                        temp.append("M")
                    elif "segment L" in dict[metadata].description:
                        temp.append("L")
                    else: temp.append("unknown")
                else: #non specified virus mode
                    temp.append("unknown")
                #find the year
                year = re.findall(r'(?<!\d)\d{4}(?!\d)', dict[metadata].description)
                posyear = [i for i in year if 1990 <= int(i) <=2035]
                temp.append(posyear[0] if len(posyear)==1 else "unknown")
                
                metadb.append(temp)
                count+=1
                #the seq is added to ntdb while name info is added to metadb
                ntdb.append(str(dict[metadata].seq.upper()))
        if len(metadb)==0:
            raise Exception(f"Genome file error: check path and filename ({os.path.join(path,fext)})")
    except FileNotFoundError or OSError:
        print(f"Genome file error: check path and filename ({os.path.join(path,fext)})")
    except KeyError:
        print("Genome file error: check data entries and format")

    if len(metadb)!= len(ntdb):
        raise Exception("Mismatch of entry titles and sequences")
    
    #metadb: [genome#, >line, segment, year]*n,
    #ntdb: [seq]*n
    return metadb, ntdb

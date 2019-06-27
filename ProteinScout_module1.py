import os
import csv
import re
from Bio import SeqIO
from Bio import SearchIO
from Bio import ExPASy
from Bio import SwissProt
from Bio.Blast import NCBIWWW
from datetime import datetime
from shutil import copyfile

databaselist = ("nr", "refseq_protein", "landmark", "swissprot", "pat", "pdb", "env_nr", "tsa_nr")
matrixlist = ("PAM30", "PAM70", "PAM250", "BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM90")
separator = ', '
regex = re.compile(r'[@!#$%^&*<>?/\|}{~:]') #forbidden characters

def spectrosingle(inputline):

    if (len(inputline) > 0) and (not inputline[0] == "\t") :
        resline = inputline.split("\t")[1]
        #resline = re.findall(r'\d\t(.+?)\t\t\t\t|$', inputline)[0]
        return resline
    else:
        return None
    
def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.
    
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
                
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch    
            
def custom_blast(sequence_input):
    """Executes NCBIWWW.qblast query on input file with specified parameters.
    """
    while True:
        bdatabase = input("Please enter the code of the pBLAST database to search against (for full list of codes, type help)\n\n")
        if bdatabase in databaselist:
            break
        elif bdatabase == "help":
            print("Non-redundant protein sequences: nr\nReference proteins: refseq_protein\nModel organisms: landmark\nUniProtKB/SwissProt: swissprot\nPatented protein sequences: pat\nProtein Data Bank proteins: pdb\nMetagenomic proteins: env_nr\nTranscriptome Shotgun Assembly proteins: tsa_nr\n")
        else:
            print("Invalid input.\n")

    bentrez = input("If you wish to restrict the search to a single organism, please enter its binomial name (e.g. Arabidopsis thaliana). If you wish to make a general search, just press Enter.\n\n")
    
    while True:
        hitlist_input = input("Please enter the maximum number of hits for each entry: \n\n")
        try:
            hitlist_val = int(hitlist_input)
            break
        except ValueError:
            print("Invalid input.\n")

    if len(bentrez) == 0:
        print("pBLAST for the file %s launched at %s\nParameters:\tDatabase: %s\t Hitcount: %i\nThis may take a very long time depending on BLAST server load and the length of your query, please wait...\n"
              "To save time, it's recommended to manually start the queries on WegoLoc, NucPred, Localizer, CELLO2GO and GOAnna now.\n\nIMPORTANT: Save the webtool results into the folder WEBTOOL_RESULTS that has been created in your specified directory.\n"
              "If a webtool has multiple result files, please make sure that they are numbered the same as their input files (e.g. for HORVU_input_nucpred00.fasta, result should be named HORVU_result_nucpred00.txt)\nand save them into the folder %s." % (fastaseq, datetime.now(), bdatabase, hitlist_val, finpath))
        output = NCBIWWW.qblast(program="blastp",database="swissprot", sequence=sequence_input.read(), hitlist_size=hitlist_val)
        
    else:
        bentrez = '"%s"[organism]' % (bentrez,)
        print("pBLAST for the file %s launched at %s\nParameters:\tDatabase: %s\t Model organism: %s\t Hitcount: %i\nThis may take a very long time depending on BLAST server load and the length of your query, please wait...\n"
              "To save time, it's recommended to manually start the queries on WegoLoc, NucPred, Localizer, CELLO2GO and GOAnna now.\n\nIMPORTANT: Save the webtool results into the folder WEBTOOL_RESULTS that has been created in your specified directory.\n"
              "If a webtool has multiple result files, please make sure that they are numbered the same as their input files (e.g. for HORVU_input_nucpred00.fasta, result should be named HORVU_result_nucpred00.txt)\nand save them into the folder %s." % (fastaseq, datetime.now(), bdatabase, bentrez, hitlist_val, finpath))
        output = NCBIWWW.qblast(program="blastp",database=bdatabase, sequence=sequence_input.read(), entrez_query=bentrez, hitlist_size=hitlist_val)
        
    return output

def slicer(inputfile, seqcount, prefix):
    """Separates a large fasta file into multiple files of 400, 1000 and 2000 sequences for the purpose of webtool usage"""
    record_slice = SeqIO.parse(inputfile,"fasta")
    for i, batch in enumerate(batch_iterator(record_slice, 400)):
        outname = prefix + "input_wegoloc_{:02d}.fasta"
        output_file = outname.format(i)
        with open(output_file, "w") as handle:
            count = SeqIO.write(batch, handle, "fasta")
        print("Wrote %i records to %s" % (count, output_file))
        
    if seqcount > 1000:
        inputfile.seek(0)
        record_slice = SeqIO.parse(inputfile,"fasta")
        for i, batch in enumerate(batch_iterator(record_slice, 1000)):
            outname = prefix + "input_nucpred_{:02d}.fasta"
            output_file = outname.format(i)
            with open(output_file, "w") as handle:
                count = SeqIO.write(batch, handle, "fasta")
            print("Wrote %i records to %s" % (count, output_file))
    
    if seqcount > 2000:
        inputfile.seek(0)
        record_slice = SeqIO.parse(inputfile,"fasta")
        for i, batch in enumerate(batch_iterator(record_slice, 2000)):
            outname = prefix + "input_localizer_{:02d}.fasta"
            output_file = outname.format(i)
            with open(output_file, "w") as handle:
                count = SeqIO.write(batch, handle, "fasta")
            print("Wrote %i records to %s" % (count, output_file))
        
        
def strappend(sourcestr, lsep='', rsep='', midsep=separator):
    """If given an array, this function takes the individual items and appends them into a single string, with items separated by a substring specified by midsep.
        Optionally, this function can append lsep and rsep specified substrings to the start and end of the strings - this is useful in case the final string needs to be in a specific format (e.g. in brackets)"""
    result = ''
    result += lsep
    for i in sourcestr:
        result += str(i)+midsep
    result += rsep
    return result

def textcutter(sourcestr, startstr, endstr, delimstr='', count=1):
    """Returns any substring that is located inside the sourcestr, between startstr and endstr. Used for extracting snippets of protein info from SwissProt records and output files from webtools.
        If there are multiple desired substrings between the startstr and endstr delimited by a uniform separator, this function can return them in an array if the delimiter and number of substrings are provided."""
    wildcard = '(.+?)' #returns anything
    metastr = wildcard+delimstr
    for i in range(count-1):
        metastr = metastr+wildcard+delimstr
    final = startstr+metastr+endstr
    n = re.findall(final, sourcestr)
    return n
#_________________________________________________________PREP SECTION_________________________________________________________
print("Protein Scout module 1, v. 1.01, by Alois Kozubik, Palacky University Olomouc, 2019\n----------------------------------------------------------------------------------\n")

while True:
    mainpath = input("Please specify a full path to the main work folder for this project (e.g. c:\\Output\\Results\\2019_05_05\\). Character \\ can be typed by pressing Right Alt+Q).\nIf the folder doesn't exist, it will be created. Alternatively, if this script is already located in the designated output folder, just press Enter:\n")
    if mainpath == '':
        mainpath = os.getcwd()
        break
    elif os.path.exists(mainpath):
        os.chdir(mainpath)
        print("Directory changed\n")
        break
    elif mainpath[1] == ':':
        try:
            os.mkdir(mainpath)
            os.chdir(mainpath)
            print("Directory changed\n")
            break
        except:
            print("Error when creating new folder!\n")
    else:
        print("Error - invalid path!\n")
            
while True:
    fastaseq = input("Please enter the full name of the FASTA FILE inside the %s folder containing sequences designated for pBLAST (e.g. protein.fasta), or the full path to a file elsewhere (e.g. c:\\folder1\\folder2\\protein.txt):\n\n" % str(os.getcwd()))
    if os.path.isfile(fastaseq):
        break
    else:
        print("File corrupted or not found!\n")

while True:
    spectseq = input("\nPlease enter the full name of the PEPTIDE SHAKER FILE inside the %s folder containing sequences designated for pBLAST (e.g. peptideshaker.txt), or the full path to a file elsewhere (e.g. c:\\folder1\\folder2\\peptideshaker.txt):\n\n" % str(os.getcwd()))
    if os.path.isfile(fastaseq):
        break
    else:
        print("File corrupted or not found!\n")

while True:
    filterdec = input("Do you wish to filter contaminants from the FASTA and Peptide Shaker files?\n\nY/N:" )          
    if filterdec in 'Yy':
        filterdec = True
        filter = input("Please enter ID prefix of sequences which you want in the result (NOT contaminants, e.g. HORVU).:")
        break
    elif filterdec in 'Nn':
        filterdec = False
        break
    else: 
        print("Invalid input\n.")

while True:
    prepid = input("If you wish to set a file prefix for better identification of files (e.g. Batch_2019_05) to all of the output files from this program, please type it below (strongly recommended if the target folder already contains multiple files - this program may overwrite some existing files).\nFor standard filenames, just press Enter:\n" )          
    if(regex.search(prepid) == None):
        if prepid == '':
            break
        else:
            prepid = prepid + '_'
            break
    else: 
        print("Identifier contains forbidden characters.")
     
respath = mainpath+'\\'+prepid+'WEBTOOL_RESULTS'
try:
    os.mkdir(respath)
except:
    pass

finpath = mainpath+'\\'+prepid+'FINAL'
try:
    os.mkdir(finpath)
except:
    pass


data = {}
with open(fastaseq, 'r') as fastafile:
    if filterdec:
        for line_id in fastafile:
            mainacc = line_id.strip()[1:] # remove newline and ">"
            if mainacc.startswith(filter):  #check if the sequence prefix is on the whitelist
                filterflag = True
            else:
                filterflag = False
            if filterflag:
                line_data = next(fastafile)  # get the sequence from file
                data[mainacc] = line_data.strip()
            else:
                try:
                    discard = next(fastafile) # discard both the ID and the sequence
                except:
                    break
    else:
        for line_id in fastafile:
            mainacc = line_id.strip()[1:] # remove newline and ">"
            line_data = next(fastafile)  # get the sequence from file
            data[mainacc] = line_data.strip()

"""
with open(fastaseq, 'r') as fastafile:
    if filterdec:
        for line_id in fastafile:
            mainacc = line_id.strip()[1:] # remove newline and ">"
            if mainacc.startswith(filter):  #check if the sequence prefix is on the whitelist
                filterflag = True
            else:
                filterflag = False
                fastacontaminants.append(mainacc)
            if filterflag:
                if mainacc not in checkarr:
                    checkarr.append(mainacc)
                    line_data = next(fastafile)  # get the sequence from file
                    data[mainacc] = line_data.strip()
                else:
                    duparr.append(mainacc)
                    try:
                        discard = next(fastafile) # discard both the ID and the sequence
                    except:
                        break # EOF break
            else:
                try:
                    discard = next(fastafile) # discard both the ID and the sequence
                except:
                    break # EOF break
    else:
        for line_id in fastafile:
            if mainacc not in checkarr:
                checkarr.append(mainacc)
                mainacc = line_id.strip()[1:] # remove newline and ">"
                line_data = next(fastafile)  # get the sequence from file
                data[mainacc] = line_data.strip()
            else:
                duparr.append(mainacc)
                try:
                    discard = next(fastafile) # discard both the ID and the sequence
                except:
                    break # EOF break
"""       

spectrocontaminants = []
duparr = []

with open(spectseq, "r") as spectrometry, open(os.path.splitext(finpath)[0]+'\\' + prepid + 'FINAL_shakerfile.txt', mode='w') as final_output: #writing the final spectrometry file with deleted sequences which would cause a mismatch during the final merger of data
    final_output.write(spectrometry.readline()) # loads the header
    for line in spectrometry:
        idno, mainacc, *rest = line.split("\t")
        final_output.write(line)

with open(mainpath+'\\'+prepid+'MAIN.fasta', mode='w') as output_handle, open(os.path.splitext(finpath)[0] + '\\' + prepid + 'FINAL_shakerfile.txt', "r") as spectrometry:
    spectrometry.readline() #read the header
    for line in spectrometry:
        idno, mainacc, *rest = line.split("\t")
        if mainacc in data:
            the_data = data.get(mainacc)
            print("Sorted sequence ID %s." % (mainacc))
            output_handle.write('>'+mainacc+'\n')
            the_data = the_data + '\n'
            output_handle.write(the_data)
        else:
            print("Filtered contaminant sequence ID %s." % (mainacc))
            spectrocontaminants.append(mainacc)

copyfile(mainpath+'\\'+prepid+'MAIN.fasta', mainpath+'\\'+prepid+'MAIN_ProteinCutter_input.txt') #ProteinCutter can't handle files with .fasta extension, so copy the file to a plain .txt

print("Sorting finished! There were \t\t%i contaminant sequences in Peptide Shaker file and %i duplicate sequences\nin the Peptide Shaker file, which are listed below and have been filtered out of the final files.\n" % (len(spectrocontaminants), len(duparr)))
print('Duplicate entries:\n', duparr)
print('Contaminants in Peptide Shaker:\n', spectrocontaminants)
spectrometry.close()

#_________________________________________________________BLAST SECTION_________________________________________________________

mainfastafile = open(mainpath+'\\'+prepid+'MAIN.fasta', "r") #open the sorted FASTA file

fastacount = 0
for line in mainfastafile: #count the number of FASTA sequences in the file
    if line.startswith(">"):
        fastacount += 1
mainfastafile.seek(0) #return to the beginning of file

if fastacount > 400:
    print("The specified file contains %i sequences. To ensure functionality of all necessary webtools, the FASTA file will be split into smaller batches.\n" % (fastacount))
    os.chdir(respath)
    slicer(mainfastafile, fastacount, prepid)
    mainfastafile.seek(0)
    os.chdir(mainpath)
    print("\nBatch files created. If you don't see a batch file with a name of a required webtool, it's because that webtool can handle the original FASTA file, so use it as that webtool's input instead.\n")

while True:
    dec_filter = input("\nDo you wish to set a positivity treshold? If you do so, any BLAST hit with lower than specified positivity percentage will be discarded.\nY/N:")
    if dec_filter in 'Yy':
        filter_input = input("\nPlease enter a numerical value between 0 and 100:\n")
        try:
            val_filter = int(filter_input)
            if val_filter >= 0 and val_filter <= 100:
                break
            else:
                print("Value is not between 0 and 100.\n")
        except ValueError:
            print("Value is not numerical.\n")
    elif dec_filter in 'Nn':
        val_filter = 0
        break
    else:
        print("Invalid input\n")

os.chdir(finpath)

while True:
    dec_default = input("\nDo you wish to use the default pBLAST settings? (Swissprot database, a single top hit for each entry, results filtered for Arabidopsis thaliana)\nY/N:")
    if dec_default in 'Yy':
        entrquery = '"Arabidopsis thaliana"[organism]'
        print("\npBLAST for the file %s\t launched at %s\n\nParameters:\nDatabase: \tSwissprot\t\t Model organism: \t%s\t\t Hitcount: \t1\nThis may take a very long time depending on BLAST server load and the length of your query, please wait...\n"
        "\nTo save time, it's recommended to manually start the queries on WegoLoc, NucPred, Localizer, CELLO2GO and GOAnna now.\n\nIMPORTANT: Save the webtool results in .txt format and with prefix 'FINAL_' into the folder %s that has been created in your specified directory.\n"
        "If a webtool has multiple result files, please make sure that they are numbered the same as their input files (e.g. for HORVU_input_nucpred00.fasta, result should be named HORVU_result_nucpred00.txt)\nand save them into the folder %s." % (mainpath, datetime.now(), entrquery, finpath, respath))
        result_handle = NCBIWWW.qblast(program="blastp",database="swissprot", sequence=mainfastafile.read(), entrez_query=entrquery, hitlist_size=1)
        break
    elif dec_default in 'Nn':
        result_handle = custom_blast(mainfastafile)
        break
    else:
        print("Invalid input.\n")

xml_filename = prepid+'blast_raw.xml'
save_file = open(xml_filename, "w")
save_file.write(result_handle.read())
save_file.close()
mainfastafile.close()
result_handle.close()
print("\npBLAST finished at %s. \nComplete results have been saved as %s to the folder %s\nProceeding with SwissProt annotation...\n" % (datetime.now(), xml_filename, finpath))

#_________________________________________________________SPROT+ANNOTATION SECTION_________________________________________________________
#Saving the data to a .tsv file
tsv_filename = finpath+'\\'+prepid+'FINAL_blast_tab.tsv'
with open(tsv_filename, mode='w', newline='') as BLAST_input:
    output_writer = csv.writer(BLAST_input, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    output_writer.writerow(['Entry', 'SwissProt Entry', 'Gene names', 'Protein names', 'Status', 'Organism', 'e-value', 'Bitscore', '% identity', '% positives', 'Gene ontology(biological process)', 'Gene ontology(cellular component)', 'Gene ontology(molecular function)', 'Subcellular location [CC]']) #csv table header
    qresults = SearchIO.parse(xml_filename, 'blast-xml')
    queryno = 0
    for qresult in qresults:
        try:
            tophit = qresult[0] #check if the sequence has a BLAST result available
            has_results = True
        except:
            has_results = False
        
        if has_results:
            tophsp = qresult[0][0]
            seqlen = tophit.seq_len
            positives = round(((tophsp.pos_num / seqlen) *100), 4)
            if positives < val_filter:
                #hit is below the positivity treshold, and is therefore discarded
                tophit = ''
                genestr = ''
                proteinstr = ''
                swissacc = ''
                datatype = ''
                organism = ''
                tophsp = ''
                evl = ''
                btscr = ''
                identity = ''
                positives = ''
                cmpstr=''
                funstr=''
                prcstr=''
                ccstr=''
                queryno = queryno + 1
                print("%i\t\tQuery %s\t\t\t\tskipped - positivity below required percentage.\n" % (queryno, qresult.id))
            else:
                #hit is above the positivity treshold, and is therefore looked up via get_sprot_raw()
                identity = round(((tophsp.ident_num / seqlen) * 100), 4)
                evl = tophsp.evalue
                btscr = tophsp.bitscore
                swissacc = tophit.accession
                swisshandle = ExPASy.get_sprot_raw(swissacc)
                records = SwissProt.parse(swisshandle)
                for record in records: #there's always just one record here, I'm just not sure if there's any smarter way to access the record than using 'for'
                    datatype = record.data_class
                    organism = record.organism
                    var = str(record.gene_name)
                    genestr = strappend(textcutter(var, r"\S*=", ";")) #select the text between "(any non-whitespace)=" and ";"
                    var = str(record.description)
                    proteinstr = strappend(textcutter(var, r"\S*=", ";"))
                    var = str(record.comments)
                    ccstr = strappend(textcutter(var, "SUBCELLULAR LOCATION: ", ".'"))
                    cmpstr=''
                    funstr=''
                    prcstr=''
                    var = str(record.cross_references)
                    xrefsource = textcutter(var, "'GO', '", "I", "', '", 2) #cut only the GO terms from the cross_references. there are 2 values inside each GO object, so we call the textcutter with the optional delimiter and repeat number = 1
                    for x in xrefsource:
                        #check what the GO term is refering to, then append it to the appropriate variable
                        if x[1][0]=='C':
                            cmpstr=cmpstr+str(x)+separator
                        elif x[1][0]=='F':
                            funstr=funstr+str(x)+separator
                        elif x[1][0]=='P':
                            prcstr=prcstr+str(x)+separator
                        else:
                            continue #GO term is discarded because it doesn't refer to either function, process or cell. component
                queryno = queryno + 1
                print("%i\t\tQuery %s\t\tsuccessfully annotated at %s.\n" % (queryno, qresult.id, datetime.now()))
        else:
            #sequence doesn't have any BLAST result, so there's no data of a similar protein to display here
            tophit = ''
            genestr = ''
            proteinstr = ''
            swissacc = ''
            datatype = ''
            organism = ''
            tophsp = ''
            evl = ''
            btscr = ''
            identity = ''
            positives = ''
            cmpstr=''
            funstr=''
            prcstr=''
            ccstr=''
            queryno = queryno + 1
            print("%i\t\tQuery %s\t\t\t\tskipped - no relevant hit available.\n" % (queryno, qresult.id))
        
        output_writer.writerow([qresult.id, swissacc, genestr, proteinstr, datatype, organism, evl, btscr, identity, positives, prcstr, cmpstr, funstr, ccstr])
        
print("File saved as %s\n" % (tsv_filename))
quitkey = input("Press Enter to close this window.\n")


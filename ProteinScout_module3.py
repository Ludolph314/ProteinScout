import os
import sys
import csv
import re
from operator import add
from functools import reduce
from decimal import Decimal

validnuc = 0
prednuc = 0
combnuc = 0
disuni = 0
dispred = 0

def textcutter(sourcestr, startstr, endstr, delimstr='', count=1):
    """Returns any substring that is located inside the sourcestr, between startstr and endstr. Used for extracting snippets of protein info from SwissProt records and output files from webtools.
        If there are multiple desired substrings between the startstr and endstr delimited by a uniform separator, this function can cut them if the delimiter and number of substrings are provided."""
    wildcard = '(.+?)' #returns anything
    metastr = wildcard+delimstr
    for i in range(count-1):
        metastr = metastr+wildcard+delimstr
    final = startstr+metastr+endstr
    n = re.findall(final, sourcestr)
    return n

def protcutsingle(inputfile, mode = False):
    """Returns an array with values of a single accession from the Protein Cutter record"""
    anyprefix = '<td class="c'
    idprefix = '<td class="c1">&gt;'
    #due to the format of the ProteinCutter HTML code, which splits the lines in the middle of sequence ID, the final thread must be assembled from two arrays of values delimited by two different substrings
    if mode:
        checkline = 'temp'
        lengtharr = []
        while checkline:
            checkline = inputfile.readline()
            try:
                lengtharr.append(int(re.findall(r'</td><td class="c0r">[0-9\.]+<br></td><td class="c0r">[0-9\.]+<br></td><td class="c0r">(.+?)</td><td class="c0r">[-0-9\.]+</td><td class="c0r">[0-9\.]+</td><td class="c0r">[0-9\.]+</td>', checkline)[0]))
            except:
                continue
        inputfile.seek(0)
        return lengtharr
        
    checkline = 'temp'
    while checkline:
        if anyprefix in checkline:
            if checkline and '</tbody>' not in checkline:
                
                pid_input1 = textcutter(checkline, idprefix, '\n')
                
                checkline = inputfile.readline()
                pid_check = checkline.partition("<") #can't use textcutter because the searched text is not delimited from the left
                if not pid_check[2]:  # ProteinCutter HTML code has some rare anomalous lines where after the second part of accession, there is not a html code following, but just a simple \n. This checks for those cases and skips to the next line.
                    pid_check = pid_check[0].partition("\n")
                    pid_input2 = pid_check[0] 
                    pid_input = pid_input1[0] + pid_input2
                    checkline = inputfile.readline()
                else:
                    pid_input2 = pid_check[0] 
                    pid_input = pid_input1[0] + pid_input2

                checkline = inputfile.readline()
                while len(checkline) < 12:
                    checkline = inputfile.readline() #line contains only a part of amino acid sequence, hence it should be skipped
                valarr = re.findall(r'</td><td class="c0r">(.+?)<br></td><td class="c0r">(.+?)<br></td><td class="c0r">(.+?)</td><td class="c0r">(.+?)</td><td class="c0r">(.+?)</td><td class="c0r">(.+?)</td>', checkline)
                mono_input = float(valarr[0][0])
                avg_input = float(valarr[0][1])
                length_input = int(valarr[0][2])
                hdrpath_input = float(valarr[0][3])
                nps_input = float(valarr[0][4])
                isoel_input = float(valarr[0][5])

                output = [pid_input, mono_input, avg_input, length_input, hdrpath_input, nps_input, isoel_input]
                pos = proteincutterfile.tell()
                proteincutterfile.seek(pos-50)
                break
            else:
                checkline = inputfile.readline()
                return None
        else:
            checkline = inputfile.readline()
    return output

def localizersingle(inputfile):
    """Returns an array with values of a single accession from the Localizer record"""
    checkline = 'temp'
    while checkline:
        checkline = inputfile.readline()
        if checkline == '\n':
            continue
        elif re.search(r"#", checkline):
            continue
        elif re.search("\tChloroplast", checkline):
            continue
        else:
            resultid = checkline.partition("\t")[0]
            resultid = resultid.partition(" ")[0]
            resultseq = re.search(r'Y \(([A-Z,]+)\)', checkline)
            if not resultseq:
                resultseq = ' '
            output = [resultid, resultseq[0]]
            return output

def wegolocsingle(inputfile):
    """Returns an array with values of a single accession from the WegoLoc record"""
    checkline = 'temp'
    while checkline:
        checkline = inputfile.readline()
        try:
            eofbool = checkline[0].isdigit()
        except:
            break
        if eofbool:
            #resultid = re.search(r'[0-9 ]+)\)', checkline)
            resfull = checkline.partition("\t ")[2]
            resfull = resfull.partition("\t")
            resid = resfull[0]
            resloc = resfull[2].partition("\t")[0]
            resfull = [resid, resloc]
            return resfull
        else:
            continue
        
def nucpredsingle(inputfile):
    """Returns an array with values of a single accession from the NucPred record"""
    checkline = 'temp'
    while checkline:
        checkline = inputfile.readline()
        if ('Sequence-ID' in checkline) and ('NucPred-score' in checkline):
            continue
        elif len(checkline) == 0 :
            continue
        else:
            resarr = checkline.partition(r" ")
            resid = resarr[0]
            resval = float(resarr[2])
            resarr = [resid, resval]
            return resarr

def cello2gosingle(inputfile):
    """Returns an array with an identifier and another array with predicted cellular localizations"""
    resvalarr = []
    checkline = 'temp'
    while checkline:
        checkline = inputfile.readline()
        checkhead = checkline.partition(r"   ")
        if checkhead[0] == 'SD':
            resid = checkhead[2].partition(r";")[0]
            checkline = inputfile.readline()
            checkline = inputfile.readline()
            checkhead = checkline.partition(r"   ")
            while checkhead[0] == 'CP':
                resvalarr.append(checkhead[2].partition(r";")[0])  #if there are multiple cell locations (max. 8 different, albeit highly improbable), return list of all of them
                checkline = inputfile.readline()
                checkhead = checkline.partition(r"   ")
            if 'Nuclear' in resvalarr:
                resarr = [resid, 'Nuclear']
                return resarr
            else:
                resarr = [resid, '']
                return resarr
        else:
            continue

def shakersingle(inputfile, mode = False):
    checkline = 'temp'
    if mode: #mode for extracting only the validated PSMs into an array (needed for NSAF calculation)
        psmsarr = []
        while checkline:
            checkline = inputfile.readline()
            try:
                psmsarr.append(round(float(checkline.split("\t")[22])))
            except:
                inputfile.seek(0)
                checkline = inputfile.readline()
                return psmsarr
        inputfile.seek(0)
        checkline = inputfile.readline()
        return psmsarr
            
    while checkline:
        checkline = inputfile.readline()
        if (len(checkline) > 0) and (not checkline[0] == "\t") :
            resline = re.findall(r'(.+?)\t(.+?)\t\t\t\t(.+?)\t\d{1,3}.\d{1,2}\t(.+?)\t(.+?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.+?)\t(.*?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\n', checkline)[0]
            return resline
        elif (len(checkline) > 0) and (checkline[0] == '\t') :
            resline = re.findall(r'\t(.+?)\t[A-Za-z ]+?\t[A-Za-z ]+?\t[A-Za-z ]+?\t(.+?)\t[A-Za-z %\[\]]+?\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\n', checkline)[0]
            return resline
        else:
            continue
            
def blastsingle(inputfile):
    checkline = inputfile.readline()
    while checkline:
        checkline = inputfile.readline()
        if len(checkline) > 0:
            resline = re.findall(r'(.+?)[;\n]', checkline)
            return resline
        else:
            continue

def judge(GO, GOperc, GOtresh, nucval, nuctres, wegoval, nlsval, celloval):
    score = 0
    global validnuc, combnuc, prednuc, disuni, dispred
    verdictstr = ''
    
    if GOperc == -1:
        flag = 0 #No uniprot entry available
    elif GOperc >= GOtresh and ('nucleus' in GO or 'Nucleus' in GO) and ('chromosome' in GO or 'Chromosome' in GO):
        flag = 1 #uniprot entry available and indicates nuclear
        verdictstr = 'Possibly chromosomal - '
    elif GOperc >= GOtresh and ('nucleus' in GO or 'Nucleus' in GO):
        flag = 1 #uniprot entry available and indicates nuclear
    else:
        flag = 2 # uniprot entry available and shows discrepancy
    if nucval >= nuctres:
        score = score + 1
    if 'nucleus' in wegoval:
        score = score + 1
    if 'Y' in nlsval:
        score = score + 1
    if 'Nuclear' in celloval:
        score = score + 1
    if flag == 1 and score >= 3 :
        validnuc = validnuc + 1
        verdictstr = verdictstr + 'Validated nuclear'
        return verdictstr
    elif flag == 1 and score == 2 :
        combnuc = combnuc + 1
        verdictstr = verdictstr + 'Combined nuclear'
        return verdictstr
    elif flag == 0 and score >= 3 :
        prednuc = prednuc + 1
        verdictstr = verdictstr + 'Predicted nuclear'
        return verdictstr
    elif flag == 2 and score >= 3 :
        disuni = disuni + 1
        verdictstr = verdictstr + 'Discrepancy - Uniprot'
        return verdictstr
    elif flag == 1 and score == 1 :
        dispred = dispred + 1
        verdictstr = verdictstr + 'Discrepancy - prediciton'
        return verdictstr
    else:
        return verdictstr

#_________________________________________________________MAIN CODE_________________________________________________________

print('Protein Scout module 3, v. 1.01, by Alois Kozubik, Palacky University Olomouc, 2019\n----------------------------------------------------------------------------------\n'
      'USER ADVICE: filenames accepted by this module are always in the following format: "any possible prefix + FINAL_ + all lowercase webtool name + extension (either .txt or .csv)"\n\n'
      'Examples: HORVU0506_FINAL_wegoloc.txt\t\tFINAL_nucpred.txt\t\tFINAL_proteincutter.txt\t\tFINAL_cello2go.txt\t\tFINAL_localizer.txt\t\tHORVU5_FINAL_blast_tab.csv\t\t2019_05_09_FINAL_spectrometry.txt')
while True:
    usrpath = input("\nPlease specify a full path (e.g. c:\\Output\\Results\\FINAL) to the folder FINAL, created with module 1 and containing the files listed below.\n\n- Final webtool output files for WegoLoc, NucPred, ProteinCutter, Cello2GO and Localizer (e.g. FINAL_wegoloc.txt)\n"
                    "- BLAST output file from module A (e.g. FINAL_blast_tab.csv)\n- Peptide Shaker output file with table-separated values (e.g. FINAL_shakerfile.txt). \n\nCharacter \\ can be typed by pressing Right Alt+Q.\n"
                    "Alternatively, if this script is located in the designated output folder, just press Enter\n\nWARNING: The folder must contain ONLY the files specified above and only from one run. Having any extra files with similar naming pattern in the folder may result in ID mismatch.\n"
                    "If the original FASTA file was split, merge the multiple webtool results with module B first and put the FINAL files into a separate folder.\n")
    if usrpath == '':
        usrpath = os.getcwd()
        break
    elif os.path.exists(usrpath):
        # Change the current working Directory    
        os.chdir(usrpath)
        print("Directory changed\n")
        break
    else:
        print("Error when opening the specified folder!\n")

        
while True:
    treshold_nuc = input("Please specify a treshold for NucPred score (decimal number between 0 and 1). Default value is 0.5, which represents 70% specifity and 62% sensitivity. For sensitivity and specifity values on different tresholds, please consult the NucPred website.\n")
    try:
        treshold_nuc = float(treshold_nuc)
        break
    except:
        print("Error - please enter a numerical value with dot as a decimal mark (e.g. 0.48, NOT 0,48)\n")
    if treshold_nuc >= 0 and treshold_nuc <= 1:
        break
    else:
        print("Error - please enter a decimal number between 0 and 1")

#open all the files
fileflags = [True] * 7 #ensures that only one file of each kind loads and saves time on file.endswith() operation
namelist = ["FINAL_blast_tab.tsv", "FINAL_cello2go.txt", "FINAL_localizer.txt", "FINAL_nucpred.txt", "FINAL_proteincutter.txt", "FINAL_shakerfile.txt", "FINAL_wegoloc.txt"]
for file in os.listdir(usrpath):
    if fileflags[0] and file.endswith(namelist[0]):
        blastfile = open(file, mode='r')
        blasttab = csv.reader(blastfile, delimiter = '\t')
        blastheader = next(blasttab)
        print("File %s loaded" % (str(file)))
        blasttab_linecount = 0
        fileflags[0] = False
    elif fileflags[1] and file.endswith(namelist[1]):
        cello2gofile = open(file, mode='r') #output = [resid, resvalarr]
        print("File %s loaded" % (str(file)))#            this is array
        fileflags[1] = False
    elif fileflags[2] and file.endswith(namelist[2]):
        localizerfile = open(file, mode='r')  #output = [resultid, resultseq[0]]
        print("File %s loaded" % (str(file)))
        fileflags[2] = False
    elif fileflags[3] and file.endswith(namelist[3]):
        nucpredfile = open(file, mode='r')  #output = [resid, resval]
        print("File %s loaded" % (str(file)))
        fileflags[3] = False
    elif fileflags[4] and file.endswith(namelist[4]):
        proteincutterfile = open(file, mode='r')
        cutterheader = ['pid_input', 'Mono', 'Avg', 'Length', 'HdrPath (GRAVY)', 'NPS', 'pI (Isoel. point)']
        print("File %s loaded" % (str(file)))
        fileflags[4] = False
    elif fileflags[5] and file.endswith(namelist[5]):
        shakerfile = open(file, mode='r')
        shakerheader = list(shakersingle(shakerfile))
        print("File %s loaded" % (str(file)))
        fileflags[5] = False
    elif fileflags[6] and file.endswith(namelist[6]):
        wegolocfile = open(file, mode='r')  #output = [resid, resloc]
        print("File %s loaded\n\n" % (str(file)))
        fileflags[6] = False
    elif True not in fileflags: # checks if all files have been loaded
        break
    else:
        continue
    
indices = [i for i, x in enumerate(fileflags) if x == True]
if indices:
    missingfiles = []
    for i in indices:
        missingfiles.append(namelist[i])
    print("Error! File(s) %s not found in specified folder!"% (missingfiles))
    input("Press Enter to exit the program.\n")
    sys.exit()

valpsmsarr = shakersingle(shakerfile, True)
protlengtharr = protcutsingle(proteincutterfile, True)

NSAFarr = []
SAFarr = [round(x/y, 10) for x, y in zip(valpsmsarr, protlengtharr)]

nsafcoef = round(sum(SAFarr), 10)
for i in SAFarr:
    NSAFarr.append(i/nsafcoef)

with open('FINAL_annotation.csv', mode='w', newline='') as BLAST_input:
    output_writer = csv.writer(BLAST_input, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    otherheader = ['ID', 'SAF', 'NSAF', 'NucPred Score', 'WegoLoc', 'NLS', 'Cello2GO', 'FINAL']
    finalheader = reduce(add,[otherheader[:1], shakerheader[:2], cutterheader[1:], shakerheader[2:], otherheader[1:3], blastheader[:], otherheader[3:]])
    output_writer.writerow(finalheader) #csv table header
    maincount = 0
    iter = -1
    while True:
        result7 = cello2gosingle(cello2gofile) # check if there are still more entries left ahead. Using result7 because it is small
        iter = iter + 1
        if result7:
            try:
                result1 = list(shakersingle(shakerfile))
                result2 = protcutsingle(proteincutterfile)
                result3 = next(blasttab)
                result4 = nucpredsingle(nucpredfile)
                result5 = wegolocsingle(wegolocfile)
                result6 = localizersingle(localizerfile)
                
                
                if result1[1] == result2[0] == result3[0] == result4[0] == result5[0] == result6[0] == result7[0]: #check for ID mismatch (if all of files were obtained via module 1 or the FASTA file output from module 1, all of the files will have the same order of sequences, which removes the need to search every file for each matching sequence)
                    #spectral count = 

                    treshold_GO = 0
                    
                    GOcheck = result3[1]
                    if len(GOcheck) == 0:
                        verdict = judge(str(result3[11]+result3[13]), -1, treshold_GO, result4[1], treshold_nuc, result5[1], result6[1], result7[1])  #sequence doesn't have an Uniprot entry, so we set the percentage to -1, which also serves as a special flag to sign the absence of entry, which is crucial for making the final verdict
                    else:
                        verdict = judge(str(result3[11]+result3[13]), float(result3[9]), treshold_GO, result4[1], treshold_nuc, result5[1], result6[1], result7[1])
                    
                    valarr = [SAFarr[iter], NSAFarr[iter], verdict]
                    resultarr = reduce(add, [result1[:3], result2[1:], result1[3:], valarr[:2], result3[:], result4[1:], result5[1:], result6[1:], result7[1:], valarr[3:]])
                    resultarr.append(verdict)
                    
                    output_writer.writerow(resultarr)
                    print("Sequence %s\t\tannotated" % (result1[1]))
                    maincount = maincount + 1
                    
                else:
                    print("ID Mismatch at index %s!\n" % (result1[1]))
                    break
                
            except Exception as e:
                print("Critical error during file handling!\n", result1, '\n\n', result2, '\n\n', result3, '\n\n', result4, '\n\n', result5, '\n\n', result6, '\n\n', result7, '\n\n', SAFarr[iter], '\n\n', NSAFarr[iter], '\n\n', e)
                break
        else:
            validperc = (validnuc / maincount) * 100
            predperc = (prednuc / maincount) * 100
            combperc = (combnuc / maincount) * 100
            disuperc = (disuni / maincount) * 100
            dispperc = (dispred / maincount) * 100
            print("Full annotation finished! File saved as %s.\n\n"
                  "SAF sum:\t\t\t\t\t\t%f\n"
                  "Total sequences:\t\t\t\t\t\t%i\n"
                  "Validated nuclear: \t\t\t\t\t\t%i (%% %i)\n"
                  "Predicted nuclear: \t\t\t\t\t\t%i (%% %i)\n"
                  "Combined nuclear: \t\t\t\t\t\t%i (%% %i)\n"
                  "Uniprot discrepancies: \t\t\t\t\t\t%i (%% %i)\n"
                  "Prediction discrepancies:\t\t\t\t\t\t%i (%% %i)\n" % (usrpath+'\FINAL_annotation.tsv', nsafcoef, maincount, validnuc, validperc, prednuc, predperc, combnuc, combperc, disuni, disuperc, dispred, dispperc))
            with open('FINAL_statistics.txt', 'w') as statistics:
                statistics.write("Total sequences: \t\t\t%i\n" % (maincount))
                statistics.write("Validated nuclear: \t\t\t%i (%% %i)\n" % (validnuc, validperc))
                statistics.write("Predicted nuclear: \t\t\t%i (%% %i)\n" % (prednuc, predperc))
                statistics.write("Combined nuclear: \t\t\t%i (%% %i)\n" % (combnuc, combperc))
                statistics.write("Uniprot discrepancies: \t\t%i (%% %i)\n" % (disuni, disuperc))
                statistics.write("Prediction discrepancies: \t%i (%% %i)\n" % (dispred, dispperc))
            break
blastfile.close()
shakerfile.close()
proteincutterfile.close()
nucpredfile.close()
wegolocfile.close()
localizerfile.close()
cello2gofile.close()

quitkey = input("Press Enter to close this window.\n")
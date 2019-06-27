import os
import re
import glob

print('ProteinPacker module 2, v. 1.01, by Alois Kozubik, Palacky University Olomouc, 2019\n----------------------------------------------------------------------------------\n'
      'USER ADVICE - Please note that ALL result files inside the specified folder must have a .txt extension and must be uniformly named as "result_webtool_number.txt"'
      '(with optional ID prefix, if you chose one in the module A), e.g.: \n\nresult_wegoloc_00.txt, result_wegoloc_01.txt for WegoLoc results\nHORVU5_result_nucpred_01.txt, '
      'HORVU5_result_nucpred_02.txt for NucPred (with HORVU5_ as prefix)\nresult_localizer_00.txt, result_localizer_01.txt\n\notherwise the intended functionality of this script is not guaranteed.\n\n')

while True:
    usrpath = input("Please specify full path to an existing folder containing all mergeable webtool result files  (e.g. c:\\Output\\Results\\)\nCharacter \\ can be typed by pressing Right Alt+Q). Alternatively, if the files are located in the same folder as this script, just press Enter:\n")
    if usrpath == '':
        usrpath = os.getcwd()
        break
    elif os.path.exists(usrpath):
        # Change the current working Directory    
        os.chdir(usrpath)
        print("Directory changed\n")
        break
    else:
        print("Can't find the specified directory, please try again.\n\n")

if usrpath[-1:] != '\\':
    usrpath = usrpath+"\\"

extension = '.txt' # allowed extension
toolnames = ['wegoloc', 'localizer', 'nucpred', 'cello2go'] # allowed filenames
for toolname in toolnames:
    try:
        filenames = [] # initialize array of filenames
        # Search the folder for all files with allowed filenames and extensions (e.g. result_wegoloc_01.txt)
        filenames.extend(glob.glob(usrpath + '*' + 'result_' + toolname + '_'+('[0-9]' * 2)+extension)) 
        if not filenames:
            print("No result files for %s found!" % (toolname)) # no file for the webtool was found
        else:
            with open('FINAL_'+toolname+extension, 'w') as outfile: # files found, output file created
                for fname in filenames:
                    with open(fname) as infile: # write all found files into the output in the correct order
                        for line in infile:
                            if line.endswith('\n'):
                                outfile.write(line)
                            else:
                                outfile.write(line+'\n')
                print('Results from %s merged successfully.\n' % (toolname))
    except:
        print("Error - file merging failed!\n")

quitkey = input("Press Enter to close this window.\n")

import os.path
import shutil
import sys
import time
import re

# function to find all positions in a character
def findpos(ch, string1):
    last=[]
    for i in range(len(string1)):
        if(string1[i]==ch):
            last.append(i)
    return last

# prints current directory
print("Current directory is: " + os.getcwd())
# check if directory tosend exists
exists = os.path.isdir('tosend')
# if it does, compress it as a backup with date and time
if exists:
    timestr = time.strftime("%Y%m%d-%H%M%S")
    concat = 'tosend' + timestr
    # Creating the ZIP file and remove directory
    archived = shutil.make_archive(concat.strip(), 'zip', 'tosend')
    shutil.rmtree('tosend', ignore_errors=True)
# create the new directory
os.mkdir("tosend")
# folder path
dir_path = r'.'
# list to store files
resl = []
# Iterate for lyx (1)
nl = len(sys.argv)
if nl == 2:
    resl.append(sys.argv[1])
else:
    for file in os.listdir(dir_path):
    # check only lyx files
        if file.endswith('.lyx'):
            resl.append(file)
print(r"Lyx Files:")
print(resl)
if len(resl) == 1:
    os.system("lyx --export latex " + resl[0])
# Iterate for tex (2)
basefile=resl[0].removesuffix(".lyx")
latexfile=basefile+".tex"
print(r"Tex File:")
print(latexfile)
lists = []
file = open(latexfile, 'r', encoding='latin-1')
print(file.readable())
count = 0
countgraphics = 0
listofgraphics = []
countbib = 0
# Loop for file lines
for line in file:
    count += 1
    lists.append(line)
    if r'includegraphics' in line:
        countgraphics = countgraphics + 1
        line = line.rstrip('\n')
        # whole path
        istart = line.find(r'{')
        iend = line.find(r'}')
        subst = line[istart + 1:iend]
        # just the file
        allpos = findpos(r'/', line)
        filename = ""
        if (len(allpos) != 0):
            filename = line[allpos[len(allpos) - 1]+1:iend]
        else:
            filename = subst
        print("filename={}".format(filename))
        print("substring={}".format(subst))
        if filename in listofgraphics:
            print("File {} is duplicated".format(filename))
            sys.exit()
        listofgraphics.append(filename)
        subst = subst.removesuffix(".eps")
        # copy file to "tosend"
        if os.path.isfile(subst + ".eps"):
            shutil.copy2(subst + ".eps", "tosend")
        subst = subst.removesuffix(".jpg")
        # copy file to "tosend"
        if os.path.isfile(subst + ".jpg"):
            shutil.copy2(subst + ".jpg", "tosend")
        subst = subst.removesuffix(".pdf")            
        # copy file to "tosend"
        if os.path.isfile(subst + ".pdf"):
            shutil.copy2(subst + ".pdf", "tosend")            
        if (len(allpos) > 0):
            front = line[:istart+1]  # up to but not including n
            back = line[allpos[len(allpos) - 1] + 1:]  # n+1 through end of string
            lists[count - 1] = front + back
            print("substring={}".format(lists[count-1]))
    elif r'bibliography{' in line:
        countbib = countbib + 1
        line = line.rstrip('\n')
        # whole path
        istart = line.find(r'{')
        iend = line.find(r'}')
        subst = line[istart + 1:iend] + r".bib"
        print("bib file={}".format(subst))
        if os.path.isfile(subst):
            shutil.copy2(subst, r"tosend")

print("Number of graphics {}".format(countgraphics))
print("Number of bibs {}".format(countbib))
assert (countbib == 1)
file.close()
# now creates a temporary tex file
fileo = open("temporary.tex", 'w', encoding='latin-1')
for line in lists:
    fileo.write(line)
fileo.close()
shutil.move("temporary.tex", "tosend/manuscript.tex")
for file in os.listdir(dir_path):
    # check only lyx files
    if file.endswith('.cls'):
        shutil.copy(file,"tosend/")
for file in os.listdir(dir_path):
    # check only lyx files
    if file.endswith('.sty'):
        shutil.copy(file,"tosend/")

# changes to the new directory
os.chdir(r"tosend")
print("Current directory: " + os.getcwd())
os.system("pdflatex manuscript.tex")
os.system("pdflatex manuscript.tex")
os.system("bibtex manuscript.aux")
os.system("pdflatex manuscript.tex")
os.system("bibtex manuscript.aux")
os.system("pdflatex manuscript.tex")
# now tackles the tex and bbl
lists = []
file = open("manuscript.tex", 'r', encoding='latin-1')
count = 0
# Using for loop
for line in file:
    count += 1
    if r'bibliography{' in line:
        fileb = open("manuscript.bbl", 'r', encoding='latin-1')
        for bibs in fileb:
            lists.append(bibs)
        fileb.close()
    else:
        lists.append(line)
file.close()
file = open("manuscript.tex", 'r', encoding='latin-1')
# extract abstract and keywords
# Regular expression to capture everything between \begin{abstract} and \end{abstract}
# (?s) or re.DOTALL lets '.' match newlines as well.
content = file.read()
pattern = r'\\begin{abstract}(.*?)\\end{abstract}'
match = re.search(pattern, content, flags=re.DOTALL)
if match:
    # The first capture group will be the text inside the environment
    abstract_text = match.group(1).strip()
    # Write the abstract to the output file
    with open("abstract.txt", 'w', encoding='utf-8') as out:
        out.write(abstract_text)
    print(f"Abstract extracted to abstract.txt")
else:
   print("No abstract environment found in the file.")

pattern_env = r'\\begin\{keyword\}(.*?)\\end\{keyword\}'
match_env = re.search(pattern_env, content, flags=re.DOTALL)
if match_env:
    keywords_text = match_env.group(1).strip()
else:
    # If not found, fall back to \keywords{...}
    pattern_macro = r'\\keywords\s*\{(.*?)\}'
    match_macro = re.search(pattern_macro, content, flags=re.DOTALL)
    if match_macro:
        keywords_text = match_macro.group(1).strip()
    else:
        keywords_text = ""
        
if keywords_text:
    with open("keywords.txt", 'w', encoding='utf-8') as out:
        out.write(keywords_text)
        print(f"Keywords extracted to keywords.txt")
else:
    print("No keywords found in environment or macro form.")
        
file.close()

# now inserts again in a temporary
fileo = open("temporary.tex", 'w', encoding='latin-1')
for line in lists:
    fileo.write(line)
fileo.close()
os.system("rm *.bib")
shutil.move("temporary.tex", "manuscript.tex")
os.system("pdflatex manuscript.tex")
os.system("pdflatex manuscript.tex")
os.chdir('..')
archived = shutil.make_archive("tosend", 'zip', 'tosend')
print("Finalized the operation")

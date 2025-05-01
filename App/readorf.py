#ReadORF
import gzip

class ReadORF:
    '''
    ReadORF: class object that handles whole genome sequenced files (eg fasta/fastq) to extract ORF information
    !!! DOES NOT WORK IF YOUR FILE CONTAINS MULTIPLE SEQUENCES, ONLY SUPPORTS ONE LARGE GENOME SEQUENCE !!!
    '''
    #works seamlessly with large files, because it calculates ORF's one by one and so u have to orfsearch for each orf file

    #initializer for our object
    def __init__(self, file):
        #print('ReadORF has been called')
        #INSTANCE ATTRIBUTES:
        self.file=file
        self.count=0

        #AA to Codon
        self.codontable={
            #Codon table, technically codon in translation use U instead of T because it is using rna
            #but for the purpose of this app, im reading the DNA to locate ORF's so I will be using T to detect codons
            'F':['TTT', 'TTC'], #Phe-Phenylalanine
            'L':['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], #Leu-Leucine
            'I':['ATT','ATC','ATA'], #Ile-Isoleucine
            'M':['ATG'], #Met-Methionine START Codon
            'V':['GTT','GTC','GTA','GTG'], #Val-Valine
            'S':['TCT','TCC','TCA','TCG','AGT','AGC'], #Ser-Serine
            'P':['CCT','CCC','CCA','CCG'], #Pro-Proline
            'T':['ACT', 'ACC','ACA','ACG'], #Thr-Threonine
            'A':['GCT','GCC','GCA','GCG'], #Ala-Alanine
            'Y':['TAT','TAC'], #Tyr-Tyrosine
            '*':['TAA','TAG','TGA'], #STOP Codon
            'H':['CAT','CAC'], #His-Histidine
            'Q':['CAA','CAG'], #Gln-Glutamine
            'N':['AAT','AAC'], #Asn-Asparagine
            'K':['AAA','AAG'], #Lys-Lysine
            'D':['GAT','GAC'], #Asp-Aspartic acid
            'E':['GAA','GAG'], #Glu-Glutamic acid
            'C':['TGT','TGC'], #Cys-Cysteine
            'W':['TGG'], #Trp-Tryptophan
            'R':['CGT','CGC','CGA','CGG','AGA','AGG'], #Arg-Arginine
            'G':['GGT','GGC','GGA','GGG'], #Gly-Glycine
        }
        #Codon to AA for direct calls
        self.codon_aa_table = {
            'TCA': 'S',  #Serine
            'TCC': 'S',  #Serine
            'TCG': 'S',  #Serine
            'TCT': 'S',  #Serine
            'TTC': 'F',  #Phenylalanine
            'TTT': 'F',  #Phenylalanine
            'TTA': 'L',  #Leucine
            'TTG': 'L',  #Leucine
            'TAC': 'Y',  #Tyrosine
            'TAT': 'Y',  #Tyrosine
            'TAA': '*',  #Stop *
            'TAG': '*',  #Stop *
            'TGC': 'C',  #Cysteine
            'TGT': 'C',  #Cysteine
            'TGA': '*',  #Stop *
            'TGG': 'W',  #Tryptophan
            'CTA': 'L',  #Leucine
            'CTC': 'L',  #Leucine
            'CTG': 'L',  #Leucine
            'CTT': 'L',  #Leucine
            'CCA': 'P',  #Proline
            'CCC': 'P',  #Proline
            'CCG': 'P',  #Proline
            'CCT': 'P',  #Proline
            'CAC': 'H',  #Histidine
            'CAT': 'H',  #Histidine
            'CAA': 'Q',  #Glutamine
            'CAG': 'Q',  #Glutamine
            'CGA': 'R',  #Arginine
            'CGC': 'R',  #Arginine
            'CGG': 'R',  #Arginine
            'CGT': 'R',  #Arginine
            'ATA': 'I',  #Isoleucine
            'ATC': 'I',  #Isoleucine
            'ATT': 'I',  #Isoleucine
            'ATG': 'M',  #Methionine
            'ACA': 'T',  #Threonine
            'ACC': 'T',  #Threonine
            'ACG': 'T',  #Threonine
            'ACT': 'T',  #Threonine
            'AAC': 'N',  #Asparagine
            'AAT': 'N',  #Asparagine
            'AAA': 'K',  #Lysine
            'AAG': 'K',  #Lysine
            'AGC': 'S',  #Serine
            'AGT': 'S',  #Serine
            'AGA': 'R',  #Arginine
            'AGG': 'R',  #Arginine
            'GTA': 'V',  #Valine
            'GTC': 'V',  #Valine
            'GTG': 'V',  #Valine
            'GTT': 'V',  #Valine
            'GCA': 'A',  #Alanine
            'GCC': 'A',  #Alanine
            'GCG': 'A',  #Alanine
            'GCT': 'A',  #Alanine
            'GAC': 'D',  #Aspartic Acid
            'GAT': 'D',  #Aspartic Acid
            'GAA': 'E',  #Glutamic Acid
            'GAG': 'E',  #Glutamic Acid
            'GGA': 'G',  #Glycine
            'GGC': 'G',  #Glycine
            'GGG': 'G',  #Glycine
            'GGT': 'G',  #Glycine
        }

        #Complement table:
        self.complementtable={
            'A': 'T',
            'C': 'G',
            'T': 'A',
            'G': 'C'
        }

        #File opening/handling initiatlizing:
        self.openwith=''
        self.writemode=''
        if self.file[-2:]=='gz':
            self.openwith=gzip.open
            self.writemode='wb'
        else:
            self.openwith=open
            self.writemode='w'

    def subset(self, out, max):
        #used for TESTING
        #should probably make this static method where I can choose any file to subset?

        with self.openwith(self.file, 'r') as fastq:
            with open(out, self.writemode) as out:

                for lines in fastq:
                    self.count+=1
                    if self.count<=max:
                        out.write(lines)
                    else: break
        self.count=0 #reset self.count to 0 so if called into another class method, it behaves correctly

    def preview(self, file, number):
        #used for TESTING
        #returns a print of each line upto <number> lines
        with self.openwith(file, 'r') as filein:
            self.count=0

            for lines in filein:
                self.count+=1
                if self.count<=number:
                    print(lines)
                else: break

        return

    def header(self):
        #returns header of sequencing file
        with self.openwith(self.file, 'r') as filein:
            header=filein.readline()

        return header

    def complement(self, sequence):
        '''
        Returns the complement sequence of sequence inserted.
        :param sequence:  DNA Nucleotide sequence (str)
        :return: Complemented sequence (str)
        '''
        #could make this static method if i remove complementtable from __init__
        seq=[]

        for base in sequence: #Iterate over all nucleotide bases in the sequence
            if base in self.complementtable: #Find the base in complement table
                seq.append(self.complementtable[base]) #append its complement into seq
            else:
                seq.append(base) #if there is N, append it instead of deleting

        return (''.join(seq)) #return one string of the complement strand

    def ORF1(self, out):
        '''

        Translates a <self.file> fasta file or file with a nucleotide sequence into its amino acids
        (eg: 5' ACT GCT AGT A 3' > 5' (ACT) (GCT) (AGT) (A) 3' TRANSLATE > 5' T A S 3'

        :param out: Out file containing translated amino acid sequence
        :return:
        '''
        with self.openwith(self.file, 'r') as filein, open(out, self.writemode) as fileout:

            header=filein.readline()
            fileout.write(header)
            codonseq=[]

            for lines in filein:
                lines=lines.rstrip() #remove trailing artefacts and whitespace
                aa_seq=[] #Stores the aa_sequence of each new sequence line

                for base in lines:
                    codonseq.append(base) #add a base ntd to codonseq

                    if len(codonseq)==3: #translate codon if its full
                        codonjoin="".join(codonseq) #join up codonseq list into a single string
                        codonseq=[] #empty codonseq to prepare for new codon

                        try:
                            aa_seq.append(self.codon_aa_table[codonjoin]) #direct call in dictionary for amino acid of codon
                        except:
                            aa_seq.append('-') #If the codon sequence has any other letters apart from actg (due to sequencing error or low quality/confidence base call) it will replace it with a '-'
                            #OR if codon sequence can not be found in codon_aa_table['XYZ'] then it appends '-'

                fileout.write(''.join(aa_seq)+'\n') #write line of translated amino acids and add new line for next line + neater formating

        return

    def ORF2(self, out):
        '''

        Translates a <self.file> fasta file or file with a nucleotide sequence into its amino acids with a frame shift of one
        (eg: 5' A CTG CTA GTA 3' > 5' (A) (CTG) (CTA) (GTA) 3' TRANSLATE > 5' L L V 3'

        :param out: Out file containing translated amino acid sequence
        :return:
        '''
        with self.openwith(self.file, 'r') as filein, open(out, self.writemode) as fileout:

            header=filein.readline()
            fileout.write(header)
            codonseq=[]
            linecount=0 #Linecount to prepare first sequence line for frame shift reading

            for lines in filein:
                lines=lines.rstrip() #remove trailing artefacts and whitespace
                aa_seq=[] #Stores the aa_sequence of each new sequence line

                linecount+=1 #Using line counts to frame shift for reading
                if linecount==1:
                    lines=lines[1:]

                for base in lines:
                    codonseq.append(base) #add a base ntd to codonseq

                    if len(codonseq)==3: #translate codon if its full
                        codonjoin="".join(codonseq)  #join up codonseq list into a single string
                        codonseq=[]  #empty codonseq to prepare for new codon

                        try:
                            aa_seq.append(self.codon_aa_table[codonjoin]) #direct call in dictionary for amino acid of codon
                        except:
                            aa_seq.append('-') #If the codon sequence has any other letters apart from actg (due to sequencing error or low quality/confidence base call) it will replace it with a '-'
                            #OR if codon sequence can not be found in codon_aa_table['XYZ'] then it appends '-'

                fileout.write(''.join(aa_seq)+'\n') #new line for neater formating

        return

    def ORF3(self, out):
        '''

        Translates a <self.file> fasta file or file with a nucleotide sequence into its amino acids with a frame shift of two
        (eg: 5' AC TGC TAG TA 3' > 5' (AC) (TGC) (TAG) (TA) 3' TRANSLATE > 5' C * 3'

        :param out: Out file containing translated amino acid sequence
        :return:
        '''
        with self.openwith(self.file, 'r') as filein, open(out, self.writemode) as fileout:

            header=filein.readline()
            fileout.write(header)
            codonseq=[]
            linecount=0 #Linecount to prepare first sequence line for frame shift reading

            for lines in filein:
                lines=lines.rstrip() #remove trailing artefacts and whitespace
                aa_seq=[] #Stores the aa_sequence of each new sequence line

                linecount+=1 #Using line counts to frame shift for reading
                if linecount==1:
                    lines=lines[2:]

                for base in lines:
                    codonseq.append(base) #add a base ntd to codonseq

                    if len(codonseq)==3: #translate codon if its full
                        codonjoin="".join(codonseq)  #join up codonseq list into a single string
                        codonseq=[]  #empty codonseq to prepare for new codon

                        try:
                            aa_seq.append(self.codon_aa_table[codonjoin]) #direct call in dictionary for amino acid of codon
                        except:
                            aa_seq.append('-') #If the codon sequence has any other letters apart from actg (due to sequencing error or low quality/confidence base call) it will replace it with a '-'
                            #OR if codon sequence can not be found in codon_aa_table['XYZ'] then it appends '-'

                fileout.write(''.join(aa_seq)+'\n') #new line for neater formating

        return

    def mORF1(self, out):
        '''

        Complements a <self.file> fasta file or file's nucleotide sequence and then translates it into its amino acids
        (eg: 5' ACT GCT AGT A 3' > COMPLEMENT > 5' TGA CGA TCA T 3' > CODON REVERSE >5' AGT AGC ACT T 3'> TRANSLATE > 5' S S T 3')
        #Normally when you complement it, you return it in the 3' > 5' direction eg 3' T ACT AGC AGT  5' and then TRANSLATE > 3' T S S 5'
        #But to return the 5'>3' sequence into 3'<5', I would have to read the whole sequence in and then reverse it. This will take lot of memory for very large files.
        #So I decided to just reverse codon translate it, and leave the sequence 5' > 3'

        :param out: Out file containing complemented translated amino acid sequence
        :return:
        '''
        with self.openwith(self.file, 'r') as filein, open(out, self.writemode) as fileout:

                header=filein.readline()
                fileout.write(header)
                codonseq=[]

                for lines in filein:
                    lines=lines.rstrip() #remove trailing artefacts and whitespace
                    lines=self.complement(lines) #complement the sequence
                    aa_seq=[] #Stores the aa_sequence of each new lines

                    for base in lines:
                        codonseq.append(base) #add a base ntd to codonseq

                        if len(codonseq)==3: #translate codon if its full
                            codonjoin=(''.join(codonseq)) #join up codonseq list into a single string
                            codonjoin=codonjoin[::-1]  #reverse codon string because im dealing with 3' to 5' utr now
                            codonseq=[] #empty codonseq to prepare for a new codon

                            try:
                                aa_seq.append(
                                    self.codon_aa_table[codonjoin])  #direct call in dictionary for amino acid of codon
                            except:
                                aa_seq.append(
                                    '-')  #If the codon sequence has any other letters apart from actg (due to sequencing error or low quality/confidence base call) it will replace it with a '-'
                                # OR if codon sequence can not be found in codon_aa_table['XYZ'] then it appends '-'

                    fileout.write(''.join(aa_seq)+'\n')  #new line for neater formating

        return

    def mORF2(self, out):
        '''

        Complements a <self.file> fasta file or file's nucleotide sequence and then translates it into its amino acids with a frame shift of one
        (eg: 5' A CTG CTA GTA 3' > COMPLEMENT > 5' T GAC GAT CAT 3' > CODON REVERSE >5' T CAG TAG TAC 3'> TRANSLATE > 5' Q - Y 3')
        #Normally when you complement it, you return it in the 3' > 5' direction eg 3' TAC TAG CAG T  5' and then TRANSLATE > 3' Y - Q 5'
        #But to return the 5'>3' sequence into 3'<5', I would have to read the whole sequence in and then reverse it. This will take lot of memory for very large files.
        #So I decided to just reverse codon translate it, and leave the sequence 5' > 3'

        :param out: Out file containing complemented translated amino acid sequence
        :return:
        '''
        with self.openwith(self.file, 'r') as filein, open(out, self.writemode) as fileout:

                header=filein.readline()
                fileout.write(header)
                codonseq=[]
                linecount = 0  # Linecount to prepare first sequence line for frame shift reading

                for lines in filein:
                    lines=lines.rstrip()
                    lines=self.complement(lines) #complement the sequence
                    aa_seq=[] #Stores the aa_sequence of each new lines

                    linecount += 1  # Using line counts to frame shift for reading
                    if linecount == 1:
                        lines = lines[1:]


                    for base in lines:
                        codonseq.append(base)

                        if len(codonseq)==3: #translate codon if its full
                            codonjoin=(''.join(codonseq)) #join up codonseq list into a single string
                            codonjoin=codonjoin[::-1]  #reverse codon string because im dealing with 3' to 5' utr now
                            codonseq=[] #empty codonseq to prepare for a new codon

                            try:
                                aa_seq.append(
                                    self.codon_aa_table[codonjoin])  #direct call in dictionary for amino acid of codon
                            except:
                                aa_seq.append(
                                    '-')  #If the codon sequence has any other letters apart from actg (due to sequencing error or low quality/confidence base call) it will replace it with a '-'
                                # OR if codon sequence can not be found in codon_aa_table['XYZ'] then it appends '-'


                    fileout.write(''.join(aa_seq)+'\n')  #new line for neater formating

        return

    def mORF3(self, out):
        '''

        Complements a <self.file> fasta file or file's nucleotide sequence and then translates it into its amino acids with a frame shift of one
        (eg: 5' AC TGC TAG TA 3' > COMPLEMENT > 5' TG ACG ATC AT 3' > CODON REVERSE >5' TC AGT AGT AC 3'> TRANSLATE > 5' A L 3')
        #Normally when you complement it, you return it in the 3' > 5' direction eg 3' TA CTA GCA GT  5' and then TRANSLATE > 3' L A 5'
        #But to return the 5'>3' sequence into 3'<5', I would have to read the whole sequence in and then reverse it. This will take lot of memory for very large files.
        #So I decided to just reverse codon translate it, and leave the sequence 5' > 3'

        :param out: Out file containing complemented translated amino acid sequence
        :return:
        '''
        with self.openwith(self.file, 'r') as filein, open(out, self.writemode) as fileout:

            header=filein.readline()
            fileout.write(header)
            codonseq=[]
            linecount=0  #Linecount to prepare first sequence line for frame shift reading

            for lines in filein:
                lines=lines.rstrip()  #remove trailing artefacts and whitespace
                lines = self.complement(lines)  # complement the sequence
                aa_seq=[]  #Stores the aa_sequence of each new lines

                linecount+=1  #Using line counts to frame shift for reading
                if linecount==1:
                    lines=lines[2:]

                for base in lines:
                    codonseq.append(base)  #add a base ntd to codonseq

                    if len(codonseq)==3:  #translate codon if its full
                        codonjoin=(''.join(codonseq))  #join up codonseq list into a single string
                        codonjoin=codonjoin[::-1]  #reverse codon string because im dealing with 3' to 5' utr now
                        codonseq=[]  #empty codonseq to prepare for a new codon

                        try:
                            aa_seq.append(
                                self.codon_aa_table[codonjoin])  #direct call in dictionary for amino acid of codon
                        except:
                            aa_seq.append(
                                '-')  #If the codon sequence has any other letters apart from actg (due to sequencing error or low quality/confidence base call) it will replace it with a '-'
                            # OR if codon sequence can not be found in codon_aa_table['XYZ'] then it appends '-'

                fileout.write(''.join(aa_seq)+'\n')  #new line for neater formating

        return

    def ORFSearch(self, ORFile, RESfile):
        '''

        Analyses a file from <ORF1>, <ORF2> or <ORF3> (which have been translated in 5' > 3') and writes each ORF protein found into <RESfile>.
        Each ORF protein is written with a sequence identifier as well, which follows @[amino acid start index, amino acid sequence length] >[raw file (nucleotide) start index, nucleotide sequence length].
        #ORF protein: A potential protein sequence starting with a start codon (M) and ending with an end codon (*).

        :param ORFile: <ORF1>, <ORF2>, <ORF3>
        :param RESfile: File containing ORF protein sequences
        :return: Dictionary hash map of ORF protein X : [amino acid start index, amino acid sequence length]
        #^I might remove this dictionary as I have no need for it, but it might be nice to have if it doesn't take too much memory? But this might be hard if im storing other frame shifts as well within dictionaries...

        '''
        with self.openwith(ORFile, 'r') as filein, open(RESfile, self.writemode) as fileout:

            header=filein.readline()
            #intialize some variables:
            maybe_orf=0
            ORFcount=0
            ORF=[]
            ORF_coords={
            'co-ords':0,
            'orf protein x': ['start index', 'orf protein length'],
            }

            for sequence in filein: #iterate over sequence lines in our translated ORF files
                sequence=sequence.rstrip() #remove trailing artefacts ('\n' and whitespaces

                for aa in sequence: #iterate over each aminno acid in sequence
                    ORF_coords['co-ords']+=1 # #adding count of each aa

                    if aa=='M' and maybe_orf==0: #simple and gate to start initalizing orf sequence
                        orfcoord=ORF_coords['co-ords'] #save start co-ordinate of this ORF protein

                        ORF.append(aa) #add aa to our ORF seq
                        maybe_orf=1 #constant to detect orf transcribing
                    elif maybe_orf==1: #if orf transcribing has been intiailized, this will continue transcribing the sequence
                        ORF.append(aa)
                    if aa=='*' and ORF: #simple and gate to stop transcribing of ORF seqeunce (if ORF list is NOT empty AND if sequence iterates into a STOP (*) codon:
                        maybe_orf=0 #stop transcribing
                        ORFcount+=1 #increase ORFcount
                        ORFjoin=''.join(ORF) #join ORF aa sequence
                        ORF=[] #empty ORF list to prepare for next ORF
                        orfname='ORF Protein '+str(ORFcount)

                        #Adding co-ordinates to orf_coords dict: #THIS IS MADE REDUNDANT BY orf_id[orfcoord-1, len(ORFjoin)] and being able to WRITE this as a seq id in the ORF protein header lines
                        ORF_coords[orfname]=[(orfcoord-1), len(ORFjoin)] #orfcoord-1 because i'm including the start codon (M) in the total ORF sequence (as well as end codon * but this doesnt matter on the start ORF coord)

                        orf_id=[orfcoord-1, len(ORFjoin)] #saves the co-ordinates for this specific ORF protein X and the length of its orf protein sequence
                        #^ this kinda makes the dictionary redudant for storing ORF protein X information, so I might delete adding co-ords into the dict, but counting ['co-ords'] for total seq length is nice

                        fileout.write('>'+orfname+'\n')#' @'+str(orf_id)+' >['+str(orf_id[0]*3)+', '+str(orf_id[1]*3)+'] \n') #write ORF number as mini header
                        fileout.write(ORFjoin+'\n') #write our aa joined ORF sequence

            #print(ORF_coords.items()) #TESTS

        return (ORF_coords)

    def mORFSearch(self, ORFile, RESfile):
        '''

        Analyses a file from <ORFm1>, <ORFm2> or <ORFm3> (which have been complemented and translated in 5' > 3') and writes each ORF protein found into <RESfile>
        #ORF protein: A potential protein sequence starting with a start codon (M) and ending with an end codon (*).

        :param ORFile: <ORF1>, <ORF2>, <ORF3>
        :param RESfile: File containing ORF protein sequences
        :return:
        '''

    #this algorithm works by using a for loop to iterate over the lines (sequences)
    #and then using a while loop to iterate over the letters (amino acids) of the lines.
    #the while loop uses a counter variable to index each letter within the line,
    #counter is added by 1 at start of every while loop, and continues until it equals length of lines (sequences)
    #it checks for a stop codon, and then starts transcribing the amino acid sequence,
    #until it meets another stop codon. It will then check for a start codon within that transcribed aa sequence
    #if start codon is true then, it joins the list into a string and reverses it*, and then use a find function to
    #locate index of start codon, and save the orf protein from the start codon index to write into the file.
    #if start codon is NOT found, within the *xyzxyz* sequence then the list is cleared and the counter is subtracted by one
    #so that the while loop counts again from the previous stop codon> this is in place because otherwise, the loop would search between *xyzx*yz*xy*zxy* every two stop codons
    #meaning it would find *xyzx* and then skip *yz* and start searching at *xy*. by using while loop to count until it equals len(lines) i can subtract counter by one so it
    #restarts the count at the previous stop codon and it will search between each stop codon found.

    #* it is reversed because this func searches the 3'>5' orf region and when i translated the fastq ntd base into aa sequence for 3'>5', i didnt
    #reverse the whole file/sequence, i simply just translated the codon in its place. to reverse the whole file, i would have to store the big sequence in memory and then reverse it like [::-1], which might take long

        with self.openwith(ORFile, 'r') as filein, open(RESfile, self.writemode) as fileout:

            header=filein.readline()
            #intialize some variables:
            maybe_orf=0
            ORFcount=0
            ORF=[]
            ORF_coords={
            'co-ords':0,
            'orf protein x': ['start index', 'orf protein length'],
            }


            for sequence in filein:
                sequence=sequence.rstrip()
                seqlen=len(sequence)
                counter=0

                while counter!=seqlen:
                    #THIS SKIPS THE LAST PROTEIN, SO I SHOULD CHECK THE LAST PROTEIN BECAUSE IF THE LAST PROTEIN DOESNT FIND ANOTHER STOP CODON, IT WILL NEVER REGISTER IT, AND SO I SHOULD
                    #SEARCH THE LAST ORF LIST TO SEE IF IT CONTAINS * AND M TO MAKE UP THE PROTEIN
                    aa=sequence[counter]
                    ORF_coords['co-ords'] += 1  # #adding count of each aa>

                    if aa=='*' and maybe_orf==1:
                        ORF.append(aa)
                        counter-=1
                        ORF_coords['co-ords']-=1  #since this while loop can iterate over more than the number of letters in the line, i have to subtract by one each time counter subtracts as well
                        maybe_orf=0
                        ORFjoin=''.join(ORF)[::-1] #reverse to start from start codon to end stop codon
                        ORF=[]
                        if 'M' in ORFjoin:
                            ORFcount+=1
                            mloc=ORFjoin.find('M')
                            ORFjoin=ORFjoin[mloc:]
                            orfname='ORF Protein '+str(ORFcount)
                            fileout.write('>'+orfname+'\n')
                            fileout.write(ORFjoin+'\n')

                            ORF_coords[orfname]=[(ORF_coords['co-ords']-mloc+1), len(ORFjoin)] #-mloc because co-ords right now is at stop codon, mloc finds the index of m within the sequence,
                            #so i subtract mloc from co-ords to get start position of start codon M, +1 because i want to include the M start codon
                            #>now i just have to find the real co-ords in relative to the transposed 3'>5'

                    elif maybe_orf==1:
                        ORF.append(aa)

                    elif aa=='*' and maybe_orf==0:
                        ORF.append(aa)
                        maybe_orf=1

                    counter+=1
        #problem with the algorithm above is that it only searches for sequences between * and *
        #so if its the last line and it doesnt end in *, the last remaining sequence will not be processed
        #so to process it, we can print the ORF and see if it contains a valid start codon and end codon ourselves:
            ORFjoin=''.join(ORF)[::-1]
            if 'M' in ORFjoin:
                mloc = ORFjoin.find('M')
                ORFjoin = ORFjoin[mloc:]  # cut string from M index
                orfname='ORF Protein '+str(ORFcount+1)
                fileout.write('>'+orfname+'\n')  # write ORF number as mini header
                fileout.write(ORFjoin+'\n')  # write our aa joined ORF sequence

                ORF_coords[orfname] = [(ORF_coords['co-ords']-mloc), len(ORFjoin)] #adding last orf protein to our orf_coords dict, #dont have to do +1 here because M codon is already included

        #now i have to fix the real co-ordinates, by looping thru the dict, and subtracting orf protein x co-ord [0] by the total length of the sequence ['co-ords']
        #i can double check the real co-ords by getting the real sequence and reversing it and then calculating the orfsearch53 on it
        #print(ORF_coords.items()) #this currently returns the dict orf positions in my orf (3'>5' which isnt reversed)

        #if it was reversed, the co-ords would be different, so here i fix it:
        coords=ORF_coords['co-ords']
        del ORF_coords['co-ords'], ORF_coords['orf protein x']

        ORF_coords_sorted={
            'co-ords': coords,
            'orf protein x': ['start index', 'orf protein length'],
        }

        #this loops over the remaining ORF Protein X's from old dict
        seqlen=len(ORF_coords)
        for x in range(seqlen):

            new=(ORF_coords['ORF Protein '+str(seqlen)]) #getting the values from the old dict for the last (seqlen) ORF Protein X (last because this is the first orf protein in the reversed orf gene)
            new[0]=coords-new[0] #correcting the co-ord by subtracting it by the total length
            ORF_coords_sorted['ORF Protein '+str(x+1)]=new #adding to new dictionary key in new sorted dict
            seqlen-=1 #minus seqlen by 1 so it loops over the next last protein

        ORF_coords['co-ords']=coords

        return (ORF_coords) #ORF_coords_sorted (for the real transposed co_ords) / ORF_coords (for my non transposed co_ords)

    def ntd_recall(self, resfile, orf_dict):
        """
        Searches the raw fasta/fastq sequence file to return an ORF's source nucleotide sequence.
        :param orfile: Raw fasta/fastq sequence file: <self.file>
        :param resfile: Out file to store results in
        :param orf_dict: Dictionary containing orf protein list from <ORF1/2/3>
        :return: NONE
        """
        with self.openwith(self.file, 'r') as filein, open(resfile, 'w') as fileout:
            header=filein.readline()

            counter=0
            orf_prot=1
            ntd_seq=[]

            append=0

            for lines in filein:
                lines=lines.rstrip()
                seqlen=len(lines)
                seqcount=0

                while seqcount!= seqlen:
                    base=lines[seqcount]

                    try: #try, except blocks to manage keyerror when orf protein dictionary calls out of table
                        #setting up checkpoints of orf/ntd seq co-ordinates from dictionary
                        prot=orf_dict['ORF Protein '+str(orf_prot)]
                        ntd_start=prot[0]*3
                        ntd_end=ntd_start+(prot[1]*3)
                    except KeyError:
                        break

                    if counter==ntd_start: #start transcribing nucleotide sequence
                        append=1
                        ntd_seq.append(base)

                    elif counter==ntd_end: #stop transcribing nucleotide sequence
                        append=0

                        #write to file:
                        fileout.write('>ORF Sequence '+str(orf_prot)+'\n')
                        ntd_seqjoin=''.join(ntd_seq)
                        fileout.write(ntd_seqjoin+'\n')

                        orf_prot+=1 #read for next orf protein
                        ntd_seq=[]

                        #minus 1 to counters so it backtracks if another orf protein starts immediately after one ends it will detect it,
                        seqcount-=1
                        counter-=1

                    elif append==1: #continue transcribing nucleotide sequence
                        ntd_seq.append(base)

                    counter+=1
                    seqcount+=1
        del(orf_dict) #save mem space

        return

    def mntd_recall(self, resfile, orf_dict):
        """
        Searches the raw fasta/fastq sequence file to return an ORF's source nucleotide sequence.
        :param orfile: Raw fasta/fastq sequence file: <self.file>
        :param resfile: Out file to store results in
        :param orf_dict: Dictionary containing orf protein list from <ORF1/2/3>
        :return: NONE
        """
        with self.openwith(self.file, 'r') as filein, open(resfile, 'w') as fileout:
            header=filein.readline()

            counter=0
            orf_prot=1
            ntd_seq=[]

            append=0
            coords=orf_dict['co-ords']

            for lines in filein:
                lines=lines.rstrip()
                seqlen=len(lines)
                seqcount=0

                while seqcount!= seqlen:
                    base=lines[seqcount]
                    base=self.complement(base) #complement it because this is the 3'<5'

                    try: #try, except blocks to manage keyerror when orf protein dictionary calls out of table
                        #setting up checkpoints of orf/ntd seq co-ordinates from dictionary
                        #for 3'<5' its different calcs:
                        prot=orf_dict['ORF Protein '+str(orf_prot)]
                        ntd_start=(3*(coords-prot[0]))-(prot[1]*3)
                        ntd_end=ntd_start+(prot[1]*3)
                    except KeyError:
                        break

                    if counter==ntd_start: #start transcribing nucleotide sequence
                        append=1
                        ntd_seq.append(base)

                    elif counter==ntd_end: #stop transcribing nucleotide sequence
                        append=0

                        #write to file:
                        fileout.write('>ORF Sequence '+str(orf_prot)+'\n')
                        ntd_seqjoin=''.join(ntd_seq)
                        ntd_seqjoin=ntd_seqjoin[::-1] #reversing it to match proper mORF protein seq too
                        fileout.write(ntd_seqjoin+'\n')

                        orf_prot+=1 #read for next orf protein
                        ntd_seq=[]

                        #minus 1 to counters so it backtracks if another orf protein starts immediately after one ends it will detect it,
                        seqcount-=1
                        counter-=1

                    elif append==1: #continue transcribing nucleotide sequence
                        ntd_seq.append(base)

                    counter+=1
                    seqcount+=1
        del(orf_dict) #save mem space

        return

class DBTest:
#FOR TESTS
    def __init__(self, file):
        self.file=file


    def add_entry(self):
        from . import Session
        from .database import raw_file

        entry=raw_file(file_name='FILE_NAME_1', header='>HEADER_TEST_1')
        session=Session()

        session.add(entry)
        session.commit()

        return

    def delete_entry(self):
        from . import Session
        from .database import raw_file

        session=Session()
#this deletes one
        #query= session.query(raw_file).filter(raw_file.id>=1).first()
        # session.delete(query)
#this deletes ALL entries
        session.query(raw_file).delete()
        session.commit()

        return


    @staticmethod
    def curr_dir():
        import os
        dir=os.getcwd()
        print(dir)
        abs_path=os.path.abspath(os.path.curdir)
        print(abs_path)
        return dir








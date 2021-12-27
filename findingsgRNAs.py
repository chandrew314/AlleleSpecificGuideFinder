import functools
import argparse
from multiprocessing import Process



#command line input: python new_method.py /location/for/mouse1/fasta [genome name for mouse1] /location/for/mouse1/CRISPOR /location/for/mouse2/fasta [genome name for mouse2] /location/for/mouse2/CRISPOR

# locationformouse1fasta = sys.argv[1]
# genomenameformouse1 = sys.argv[2]
# locationformouse1crispor = sys.argv[3]
#
# locationformouse2fasta = sys.argv[4]
# genomenameformouse2 = sys.argv[5]
# locationformouse2crispor = sys.argv[6]

parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group()
parser.add_argument("-PAM", help='PAM motif to be used in CRISPOR; default is NGG', nargs="?", default="NGG")
group.add_argument("-m", help='this will only output in-PAM and near-PAM guides', action='store_true')
group.add_argument("-control", help='this will only output guides that can target both species', action='store_true')
parser.add_argument('species1FastaInFile', help='location of the FASTA file for one species')
parser.add_argument('species1Genome', help='genome identifier for this species')
parser.add_argument('species1GuideOutFile', help='tab-separated file that will hold all CRISPOR-generated guides for this species')
parser.add_argument('species2FastaInFile', help='location of the FASTA file for another species')
parser.add_argument('species2Genome', help='genome identifier for this species')
parser.add_argument('species2GuideOutFile', help='tab-separated file that will hold all CRISPOR-generated guides for this species')
parser.add_argument('guideOutputFile', help='tab-separated file that will contain the output of this program')
args = parser.parse_args()

if args.PAM is not None:
    print("pam has been set " + args.PAM)

@functools.lru_cache(maxsize=None)
def make_crispor_tables():
    from subprocess import Popen

    crisporcommandlineargumentforspecies1 = ['sudo', 'python2.7', 'crispor.py', '-p', args.PAM, args.species1Genome, args.species1FastaInFile, args.species1GuideOutFile]

    crisporcommandlineargumentforspecies2 = ['sudo', 'python2.7', 'crispor.py', '-p',  args.PAM, args.species2Genome, args.species2FastaInFile, args.species2GuideOutFile]

    crisporcommandlinearguments = [crisporcommandlineargumentforspecies1, crisporcommandlineargumentforspecies2]
    procs = [Popen(i) for i in crisporcommandlinearguments]
    for p in procs:
        p.wait()
    return procs



def reverseComplement(kwargs):
    complementarynucleotides = {'T': 'A', 'G': 'C', 'C': 'G', 'A': 'T'}
    dictreversecomplement = {}
    for exonnumber in kwargs.keys():
        actuallytranssequence = ''
        reversetranssequence = ''
        for listsequences in kwargs[exonnumber]:
            transsequence = listsequences.maketrans(complementarynucleotides)
            actuallytranssequence += listsequences.translate(transsequence)
            reversetranssequence = actuallytranssequence[::-1]
        dictreversecomplement.setdefault(exonnumber, reversetranssequence)
        # dictreversecomplement[exonnumber] = ''.join(dictreversecomplement[exonnumber])
    # print(dictreversecomplement)
    return dictreversecomplement


class Mouse:
    def __init__(self, species, gene, orientation):
        self.species = species
        self.gene = gene
        self.orientation = orientation
        self.exonsequence = self.exonsequence()

    def exonsequence(self):
        from Bio import SeqIO

        if self.species == 0:
            if self.orientation == 0:  # 0 = sense
                seq_dict = {}
                i = 1
                for rec in SeqIO.parse(args.species1FastaInFile,
                                       "fasta"):  # parameter
                    seq_dict['exon ' + str(i)] = str(rec.seq)
                    i += 1
                return seq_dict

            if self.orientation == 1:  # 1 = antisense
                forwardseq_dict = {}
                i = 1
                for rec in SeqIO.parse(args.species1FastaInFile,
                                       "fasta"):  # parameter
                    forwardseq_dict['exon ' + str(i)] = str(rec.seq)
                    i += 1
                seq_dict = reverseComplement(forwardseq_dict)
                return seq_dict
        if self.species == 1:
            if self.orientation == 0:  # 0 = sense
                seq_dict = {}
                i = 1
                for rec in SeqIO.parse(args.species2FastaInFile,
                                       "fasta"):  # parameter
                    seq_dict['exon ' + str(i)] = str(rec.seq)
                    i += 1
                return seq_dict

            if self.orientation == 1:  # 1 = antisense
                forwardseq_dict = {}
                i = 1
                for rec in SeqIO.parse(args.species2FastaInFile,
                                       "fasta"):  # parameter
                    forwardseq_dict['exon ' + str(i)] = str(rec.seq)
                    i += 1
                seq_dict = reverseComplement(forwardseq_dict)
                return seq_dict

@functools.lru_cache(maxsize=None)
def import_crispor_table(mouse):
    import pandas as pd
    pd.set_option('display.max_columns', None)

    species = mouse.species
    gene = mouse.gene
    procs = make_crispor_tables()

    name = "none"
    if species == 0:
        name = args.species1Genome
    elif species == 1:
        name = args.species2Genome

    crisportable_df = 0
    poll1 = procs[0].poll()
    poll2 = procs[1].poll()
    print(args.species1GuideOutFile)
    if species == 0:
        if poll1 == 0:
            crisportable_df = pd.read_csv(args.species1GuideOutFile, sep='\t')
    elif species == 1:
        if poll2 == 0:
            crisportable_df = pd.read_csv(args.species2GuideOutFile, sep='\t')
    print(crisportable_df)
    crisportable_df['species'] = name

    return crisportable_df

@functools.lru_cache(maxsize=None)
def combine_mouse_crispor_tables(mouse1, mouse2):
    import pandas as pd
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    pd.set_option('expand_frame_repr', False)

    mouse1crisportable_df = import_crispor_table(mouse1)
    mouse2crisportable_df = import_crispor_table(mouse2)
    if args.control == False:
        combinedcrisportable_df = pd.concat([mouse1crisportable_df, mouse2crisportable_df])
        uniquecrisportable_df = combinedcrisportable_df.drop_duplicates(['targetSeq'], keep=False)
        print(uniquecrisportable_df)
        return uniquecrisportable_df
    elif args.control == True:
        combinedcrisportable_df = pd.concat([mouse1crisportable_df, mouse2crisportable_df])
        nonuniquecrisportable_df = combinedcrisportable_df[
            combinedcrisportable_df.duplicated(subset=('targetSeq'), keep=False)]
        print(nonuniquecrisportable_df)
        nonuniquecrisportable_df.to_csv(args.guideOutputFile, sep='\t')
        return nonuniquecrisportable_df

def get_exonsequence_as_a_string(mouse):
    dict_mouseexonsequence = mouse.exonsequence
    mouseexonsequence = ''

    for exonnumber in dict_mouseexonsequence.keys():
        mouseexonsequence = mouseexonsequence + dict_mouseexonsequence[exonnumber]

    return mouseexonsequence

@functools.lru_cache(maxsize=None)
def get_alignment_between_two_species(mouse1, mouse2):
    import biotite as biotite
    import biotite.sequence as seq
    from biotite.sequence.align.matrix import SubstitutionMatrix

    mouse1exonsequence = get_exonsequence_as_a_string(mouse1)
    mouse2exonsequence = get_exonsequence_as_a_string(mouse2)

    matrix = SubstitutionMatrix.std_nucleotide_matrix()

    seqmouse1 = seq.NucleotideSequence(mouse1exonsequence)
    seqmouse2 = seq.NucleotideSequence(mouse2exonsequence)
    alignment = seq.align.align_optimal(seqmouse1, seqmouse2, matrix, gap_penalty=(-6, -2), terminal_penalty=False)
    listofalignedsequences = []
    for a in range(0, len(alignment)):
        listofalignedsequences = alignment[a].get_gapped_sequences()

    sequenceidentity = biotite.sequence.align.alignment.get_sequence_identity(alignment[0], mode="not_terminal")
    if sequenceidentity < 0.8:
        print('These two genes are too far diverged.')
    else:
        print('The sequence identity is ' + str(sequenceidentity))

    print(listofalignedsequences)
    #print('huh')
    return listofalignedsequences

@functools.lru_cache(maxsize=None)
def find_index_for_each_guide(mouse1, mouse2):
    import re

    listofalignedsequences = get_alignment_between_two_species(mouse1, mouse2)
    mousealignedexonsequence = listofalignedsequences[0]
    mousecrisportable_df = combine_mouse_crispor_tables(mouse1, mouse2)
    targetseqaslist = mousecrisportable_df['targetSeq'].tolist()
    print(targetseqaslist)
    #print('why')


    listofguidestartindexes = []
    for targetseq in targetseqaslist:
        #print(targetseq)
        searchregex = r'(?=' + targetseq[0] + r'-*' + targetseq[1] + r'-*' + targetseq[2] + r'-*' + targetseq[
                3] + r'-*' + targetseq[4] + r'-*' + targetseq[5] + r'-*' + targetseq[6] + r'-*' + targetseq[7] + r'-*' + \
                          targetseq[8] + r'-*' + targetseq[9] + r'-*' + targetseq[10] + r'-*' + targetseq[11] + r'-*' + \
                          targetseq[12] + r'-*' + targetseq[13] + r'-*' + targetseq[14] + r'-*' + targetseq[
                              15] + r'-*' + targetseq[16] + r'-*' + targetseq[17] + r'-*' + targetseq[18] + r'-*' + \
                          targetseq[19] + r'-*' + targetseq[20] + r'-*' + targetseq[21] + r'-*' + targetseq[22] + r')'
        #print(searchregex)
        for guide in re.finditer(searchregex, mousealignedexonsequence):
            listofguidestartindexes.append(guide.start())
                # listofguidestartindexes = [guide.start() for guide in re.finditer(searchregex, exonsequenceforthisexon)]
                # print(listofguidestartindexes)
    #dictofguidestartindexes[exonnumber] = listofguidestartindexes
    #print(dictofguidestartindexes)

    print(listofguidestartindexes)
    #print('omg')

    lengthofguides = 23
    mouse1guides = []
    for index in listofguidestartindexes:
        necessarytopreserveindexing = mousealignedexonsequence[index:]
        mouse1guides.append(necessarytopreserveindexing.replace('-', '')[:lengthofguides])


    return mouse1guides, listofguidestartindexes, listofalignedsequences


def find_differences_between_guides(mouse1, mouse2):
    import pandas as pd
    listofalignedsequences = find_index_for_each_guide(mouse1, mouse2)[2]
    mouse1guides = find_index_for_each_guide(mouse1, mouse2)[0]
    mouse1listofguidestartindexes = find_index_for_each_guide(mouse1, mouse2)[1]
    mouse2alignedexonsequence = listofalignedsequences[1]

    lengthofguide = 23
    indexofmutationpositionneededtobeconsideredinPAM = -3

    mouse2guides = []
    for index in mouse1listofguidestartindexes:
        necessarytopreserveindexing = mouse2alignedexonsequence[index:]
        mouse2guides.append(necessarytopreserveindexing.replace('-', '')[:lengthofguide])

    print(mouse1guides)
    print(len(mouse1guides))
    print(mouse2guides)
    print(len(mouse2guides))
    #print('hi')
    listsofinformationaboutmutations = []
    inpamornotinpam = []
    extraindexes = '-----------------------'
    for guide in range(0, len(mouse1guides)):

        specificmouse1guide = mouse1guides[guide]
        specificmouse2guide = mouse2guides[guide]
        print(len(specificmouse2guide))
        print(len(specificmouse1guide))
        #print('boo')
        if len(specificmouse1guide) <= len(specificmouse2guide):
            neededindexes = extraindexes[:(lengthofguide - len(specificmouse1guide))]
            specificmouse1guide = specificmouse1guide + neededindexes
            wherethemutationsare = [i for i in range(0, len(specificmouse1guide)) if
                                specificmouse1guide[i] != specificmouse2guide[i]]
        elif len(specificmouse2guide) < len(specificmouse1guide):
            neededindexes = extraindexes[:(lengthofguide - len(specificmouse1guide))]
            specificmouse2guide = specificmouse2guide + neededindexes
            wherethemutationsare = [i for i in range(0, len(specificmouse2guide)) if
                                specificmouse1guide[i] != specificmouse2guide[i]]

        if len(wherethemutationsare) != 0:
            greatestmutationindex = max(wherethemutationsare)
            if greatestmutationindex >= lengthofguide + indexofmutationpositionneededtobeconsideredinPAM:
                inpamornotinpam.append('in-PAM')
                listsofinformationaboutmutations.append(wherethemutationsare)
            else:
                listsofinformationaboutmutations.append(wherethemutationsare)
                inpamornotinpam.append('not in-PAM')
        else:
            listsofinformationaboutmutations.append('No analogous guide found.')
            inpamornotinpam.append('in-PAM')

    return listsofinformationaboutmutations, inpamornotinpam


def make_dataframe_for_targetseq_information(mouse1, mouse2):
    import pandas as pd

    mouse1guides = find_index_for_each_guide(mouse1, mouse2)[0]
    listsofinformationaboutmutations = find_differences_between_guides(mouse1, mouse2)[0]
    inpamornotinpam = find_differences_between_guides(mouse1, mouse2)[1]

    targetseq_df = pd.DataFrame(mouse1guides, columns=['targetSeq_' + str(mouse1.species)])
    targetseq_df['mutation locations'] = listsofinformationaboutmutations
    targetseq_df['in-Pam?'] = inpamornotinpam
    if mouse1.orientation == 0:
        targetseq_df['orientation'] = 'sense'
    else:
        targetseq_df['orientation'] = 'antisense'

    print(targetseq_df)
    return targetseq_df


def length_of_exons(mouse1):
    exonsequence = mouse1.exonsequence

    dictofcumulativeexonlengths = {}
    cumulativeexonlengths = 0

    for exonnumber in exonsequence.keys():
        cumulativeexonlengths += len(exonsequence[exonnumber])
        dictofcumulativeexonlengths[exonnumber] = cumulativeexonlengths

    print(dictofcumulativeexonlengths)
    #print('help')
    return dictofcumulativeexonlengths


def target_cut_percentage(mouse1, mouse2):
    import re

    targetseq_df = make_dataframe_for_targetseq_information(mouse1, mouse2)
    mouse1guides = targetseq_df['targetSeq_' + str(mouse1.species)]

    dictofcumulativeexonlengths = length_of_exons(mouse1)
    print(dictofcumulativeexonlengths)
    #print('bruh')
    lastkeyofdictofcumulativeexonlengths = list(dictofcumulativeexonlengths.keys())[-1]
    distancefromindexwhereendonucleasecuts = 18
    totallengthofexon = dictofcumulativeexonlengths[lastkeyofdictofcumulativeexonlengths]
    print(totallengthofexon)
    lengthofeachguide = 23

    mouse1listofguidestartindexes = find_index_for_each_guide(mouse1, mouse2)[1]

    mouse1dictexonsequence = mouse1.exonsequence
    mouse1exonsequence = ''
    for exonnumber in mouse1dictexonsequence.keys():
        mouse1exonsequence = mouse1exonsequence + mouse1dictexonsequence[exonnumber]

    mouse1listofguidestartindexeswithoutdashes = [mouse1exonsequence.find(guide) for guide in mouse1guides]

    listoftargetcutindexes = []
    if mouse1.orientation == 0:
        for index in mouse1listofguidestartindexeswithoutdashes:
            listoftargetcutindexes.append(int(index + distancefromindexwhereendonucleasecuts))
    if mouse1.orientation == 1:
        for index in mouse1listofguidestartindexeswithoutdashes:
            targetcutindexantisense = int(index) + distancefromindexwhereendonucleasecuts
            targetcutindex = totallengthofexon - targetcutindexantisense
            listoftargetcutindexes.append(targetcutindex)

    listoftargetcutpercentages = [(targetcutindex / totallengthofexon) * 100 for targetcutindex in
                                  listoftargetcutindexes]

    targetseq_df['target cut %'] = listoftargetcutpercentages
    return targetseq_df


def check_with_crispor_table(mouse1, mouse2, mouse1anti, mouse2anti):
    import pandas as pd
    crispor_df = combine_mouse_crispor_tables(mouse1, mouse2)
    combined_df = ''
    if __name__ == "__main__":
        targetseqformouse1sense_df = target_cut_percentage(mouse1, mouse2)
        targetseqformouse2sense_df = target_cut_percentage(mouse2, mouse1)
        targetseqofmouse1antisense_df = target_cut_percentage(mouse1anti, mouse2anti)
        targetseqofmouse2antisense_df = target_cut_percentage(mouse2anti, mouse1anti)
        combined_df = pd.concat([targetseqformouse1sense_df, targetseqformouse2sense_df, targetseqofmouse1antisense_df,
                             targetseqofmouse2antisense_df], axis=0)
    print(combined_df)
    #print('helpp')

    mouse1merged_df = pd.merge(crispor_df, combined_df, left_on='targetSeq', right_on='targetSeq_' + str(mouse1.species),
                         how='inner')
    mouse2merged_df = pd.merge(crispor_df, combined_df, left_on='targetSeq', right_on='targetSeq_' + str(mouse2.species),
                         how='inner')
    merged_df = pd.concat([mouse1merged_df, mouse2merged_df], axis = 0)

    print(merged_df)
    #print('popcorn')
    return merged_df


def refined_crispor_table(mouse1, mouse2, mouse1anti, mouse2anti):
    import pandas as pd
    degeneratePAMnucleotidepositionandundesirablemutationlocations = [20, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3,
                                                                      2, 1, 0]
    sgrnasforaspecificgene_df = check_with_crispor_table(mouse1, mouse2, mouse1anti, mouse2anti)
    semifinal_df = sgrnasforaspecificgene_df.loc[:, ['targetSeq', 'species', 'target cut %', 'cfdSpecScore',
                                         'Doench \'16-Score', 'Moreno-Mateos-Score',
                                         'offtargetCount', 'mutation locations', 'in-Pam?', 'orientation']]  # , 'gene']

    if args.m == True:
        semifinal_df['mutation locations'] = semifinal_df['mutation locations'].astype(str)
        seriesofmutationlocations = pd.Series(semifinal_df['mutation locations'])
        print(seriesofmutationlocations)
        results = seriesofmutationlocations.str.contains(pat=r'\b(1[6-9]|2[1-2]|[A-Za-z]+)\b', regex=True)
        results.name = 'mutations locations boolean'
        print(results)
        print('sad')
        semifinal_df = pd.concat([semifinal_df, results.to_frame()], axis=1)

        final_df = semifinal_df.loc[semifinal_df['mutations locations boolean'] == True]
        superfinal_df = final_df.loc[:, ['targetSeq', 'species', 'target cut %', 'cfdSpecScore',
                                         'Doench \'16-Score', 'Moreno-Mateos-Score',
                                         'offtargetCount', 'mutation locations', 'in-Pam?', 'orientation']]
        print(superfinal_df)
        superfinal_df.to_csv(args.guideOutputFile, sep='\t')
        print('done')
	return superfinal_df
    elif args.m == False:
        semifinal_df.to_csv(args.guideOutputFile, sep='\t')
        print(semifinal_df)
        print('done')
        return(semifinal_df)



mouse1 = Mouse(species=0, gene=6, orientation=0)
mouse2 = Mouse(species=1, gene=6, orientation=0)
mouse1anti = Mouse(species=0, gene=6, orientation=1)
mouse2anti = Mouse(species=1, gene=6, orientation=1)

if __name__ == "__main__":
    if args.control == False:
        gettinginformationabobutguides = Process(target=refined_crispor_table, args=(mouse1, mouse2, mouse1anti, mouse2anti,))
        gettinginformationabobutguides.start()
        makingCRISPORtables = Process(target=make_crispor_tables)
        makingCRISPORtables.start()
    elif args.control == True:
        makingCRISPORtables = Process(target=make_crispor_tables)
        makingCRISPORtables.start()
        makingnonuniqueCRISPORtables = Process(target=combine_mouse_crispor_tables, args=(mouse1, mouse2))
        makingnonuniqueCRISPORtables.start()

exit()

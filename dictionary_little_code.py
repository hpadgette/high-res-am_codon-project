#dictionary for Mixture project
import csv, sys

####### parameters #######
time_Steps = [2, 80, 160, 600, 800]
singleAntigenTest = 1
flag = 1
cocktailFirst = 0
l = 46 #number of residues, 36 variable, 10 conserved

concsingle = 0.01
conccock   = 1.10
concseq0   = 0.25
concseq1   = 0.20
concseq2   = 0.85
GCDur1 = time_Steps[1] # 80
GCDur2 = time_Steps[2] # 130
GCDur3 = time_Steps[4] # 300

####### antigens creation #######
BG505          = "KENLWVTVYYGVPVWKDAETTLFCASDAKAYETEKHNVWATHACVPTDPNPQEIHLENVTEEFNMWKNNMVEQMHTDIISLWDQSLKPCVKLTPLCVTLQCTNVTNNITDDMRGELKNCSFNMTTELRDKKQKVYSLFYRLDVVQINENQGNRSNNSNKEYRLINCNTSAITQACPKVSFEPIPIHYCAPAGFAILKCKDKKFNGTGPCPSVSTVQCTHGIKPVVSTQLLLNGSLAEEEVMIRSENITNNAKNILVQFNTPVQINCTRPNNNTRKSIRIGPGQAFYATGDIIGDIRQAHCNVSKATWNETLGKVVKQLRKHFGNNTIIRFANSSGGDLEVTTHSFNCGGEFFYCNTSGLFNSTWISNTSVQGSNSTGSNDSITLPCRIKQIINMWQRIGQAMYAPPIQGVIRCVSNITGLILTRDGGSTNSTTETFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTRAKRR"

BG505_KR423280 = "AENLWVTVYYGVPVWKDAETTLFCASDAKAYETEKHNVWATHACVPTDPNPQEIHLENVTEEFNMWKNNMVEQMQTDIISLWDQSLKPCVKLTPLCVTLQCTNVTNNITDDMRGELKNCSFNMTTELRDKKQKVYSLFYRLDVVQINENQGNRSNNSNKEYRLISCNTSAITQACPKVSFEPIPIHYCAPAGFAILKCKDKKFNGTGPCPSVSTVQCTHGIKPVVSTQLLLNGSLAEEEVMIRSENITNNAKVILVQFNTPVQINCTRPNNNTRKSIRIGPGQAFYATGDIIGDIRQAHCNVSKATWNETLGKVVKQLRKHEGNNTIIVFANSSGGDLEITTHSFNCGGEFFYCNTSGLFNSTWISNTSVQGSNSTGSNDSITLPCRIKQIINMWQEIGQAMYAPPIQGVIRCVSNITGLILTRDGGNNTNTTEIFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTRAKRR"

BG505_HQ217523 = "AENLWVTVYYGVPVWKDAETTLFCASDAKAYETEKHNVWATHACVPTDPNPQEIHLENVTEEFNMWKNNMVEQMQTDIISLWDQSLKPCVKLTPLCVTLQCTNVTNNITDDMRGELKNCSFNMTTELRDKKQKVYSLFYRLDVVQINENQGNRSNNSNKEYRLISCNTSAITQACPKVSFEPIPIHYCAPAGFAILKCKDKKFNGTGPCPSVSTVQCTHGIKPVVSTQLLLNGSLAEEEVMIRSENITNNAKTILVQFNTPVQINCTRPNNNTRKSIRIGPGQAFYATGDIIGDIRQAHCNVSKATWNETLGKVVKQLRKHFGNNTIISFANSSGGDLEITTHSFNCGGEFFYCNTSGLFNSTWISNTSVQGSNSTGSNDSITLPCRIKQIINMWQEIGQAMYAPPIQGVIRCVSNITGLILTRDGGTNGNTTEIFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTRAKRR"

BG505_EU577271 = "AENLWVTVYYGVPVWKDAETTLFCASDAKAYETEKHNVWATHACVPTDPNPQEIHLENVTEEFNMWKNNMVEQMQTDIISLWDQSLKPCVKLTPLCVTLQCTNVTNNITDDMRGELKNCSFNMTTELRDKKQKVYSLFYRLDVVQINENQGNRSNNSNKEYRLINCNTSAITQACPKVSFEPIPIHYCAPAGFAILKCKDKKFNGTGPCPSVSTVQCTHGIKPVVSTQLLLNGSLAEEEVMIRSENITNNAKTILVQFNTPVQINCTRPNNNTRKSIRIGPGQAFYATGDIIGDIRQAHCNVSKATWNETLGKVVKQLRKHFGNNTIIIFANSSGGDLEITTHSFNCGGEFFYCNTSGLFNSTWISNTSVQGSNSTGSNDSITLPCRIKQIINMWQEIGQAMYAPPIQGVIRCVSNITGLILTRDGGNPNGTTEIFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTRAKRR"

if singleAntigenTest ==1:
        dicAgs = {c:[BG505_HQ217523] for c in list(range(time_Steps[0], time_Steps[4]))}
        dicconc= {c:concsingle for c in list(range(time_Steps[0], time_Steps[4]))}
        dicGCDur = {c:GCDur3 for c in list(range(time_Steps[0],time_Steps[4]))}
        with open('output_singleAntigenDict.csv', 'w') as output:
            w = csv.writer(output)
            w.writerows(dicAgs.items())
else:
        if cocktailFirst==1:
            dicAgs = {c: [Ag0, AgC1, AgC2] for c in list(range(time_Steps[0],time_Steps[1]))}
            dicAgs.update({c:[Ag1] for c in list(range(time_Steps[1],time_Steps[2]))})
            dicAgs.update({c:[Ag2] for c in list(range(time_Steps[2],time_Steps[3]))})

            dicconc= {c:conccock for c in list(range(time_Steps[0], time_Steps[1]))}
            dicconc.update({c:concseq1 for c in list(range(time_Steps[1], time_Steps[2]))})
            dicconc.update({c:concseq2 for c in list(range(time_Steps[2], time_Steps[3]))})

            #dicGCDur = {c:GCDur3 for c in list(range(time_Steps[0],time_Steps[3]))}
            dicGCDur = {c:GCDur1 for c in list(range(time_Steps[0],time_Steps[1]+1))}
            dicGCDur.update({c:GCDur2 for c in list(range(time_Steps[1]+1,time_Steps[2]+1))})
            dicGCDur.update({c:GCDur3 for c in list(range(time_Steps[2]+1,time_Steps[3]))})

            with open('output_cocktailDict.csv', 'w') as output:
                w = csv.writer(output)
                w.writerows(dicAgs.items())

        else:
            if flag == 1:
                dicAgs = {c:[BG505] for c in list(range(time_Steps[0],time_Steps[1]))}
                dicconc = {c:concseq0 for c in list(range(time_Steps[0], time_Steps[1]))}
                dicGCDur = {c:GCDur1 for c in list(range(time_Steps[0],time_Steps[1]+1))}
            elif flag == 2: 
                dicAgs = {c:[BG505] for c in list(range(time_Steps[0],time_Steps[1]))}
                dicconc = {c:concseq1 for c in list(range(time_Steps[0], time_Steps[1]))}
                dicGCDur = {c:GCDur1 for c in list(range(time_Steps[0],time_Steps[1]+1))}       
            elif flag == 3: 
                dicAgs = {c:[BG505] for c in list(range(time_Steps[0],time_Steps[1]))}         
                dicconc = {c:concseq2 for c in list(range(time_Steps[0], time_Steps[1]))}
                dicGCDur = {c:GCDur1 for c in list(range(time_Steps[0],time_Steps[1]+1))}

            with open('output_sequenceDict.csv', 'w') as output:
                w = csv.writer(output)
                w.writerows(dicAgs.items())

with open('output_GC_duration.csv', 'w') as output:
    w = csv.writer(output)
    w.writerows(dicGCDur.items())

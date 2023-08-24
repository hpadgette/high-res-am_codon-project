import csv

founderseq = []
with open('seed.csv') as csvfile:
    readCSV = csv.reader(csvfile)
    for row in readCSV:
        founderseq.append(str(row[0]))
print(founderseq)

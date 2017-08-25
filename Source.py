import math, re, sys, os, urllib.request
from Bio import Entrez, SeqIO
from tkinter import *
from tkinter import ttk

Entrez.email = "jrosa.balle@gmail.com" 
def getTemp(*args):
	try:
		
		cMolar = float(concMolar.get())
		cDNA = float(concDNA.get())
		if (genome.get()):
			ncbiID = str(strand.get())
			string = getSeq(ncbiID)
			temporary = 0
			for i in range (0, len(string), 70):
				temporary += (7.35 * getStrPerBase(string[i:i+72])) + (17.34 * math.log(len(string[i:i+72]))) + (4.96 * math.log(cMolar)) + (0.89 * math.log(cDNA)) - 25.42

			temp.set(round(temporary / (len(string)/70),2) - 0.5)
			# return tempC
		else:
			string = str(strand.get())
			tempC = round(((7.35 * getStrPerBase(string)) + (17.34 * math.log(len(string))) + (4.96 * math.log(cMolar)) + (0.89 * math.log(cDNA)) - 25.42), 2) 
			temp.set(tempC)

			# return tempC
	except ValueError:
		pass

def getStrPerBase(strand):
	return float(getTotalStr(strand)/len(strand))

def getTotalStr(strand): 
	tot = 0
	for i in range(len(strand)):
		dinucleotide = strand[i:i+2]
		if (len(dinucleotide) == 1):
			return tot
		elif (dinucleotide == 'GC'):
			tot += 13
		elif (dinucleotide == 'CC' or dinucleotide == 'GG'):
			tot += 11
		elif (dinucleotide == 'CG' or dinucleotide == 'AC' or dinucleotide == 'GT'):
			tot += 10
		elif (dinucleotide == 'TC' or dinucleotide == 'AG' or dinucleotide == 'CT' or dinucleotide == 'GA'):
			tot += 8
		elif (dinucleotide == 'CA' or dinucleotide == 'AT' or dinucleotide == 'TG'):
			tot += 7
		elif (dinucleotide == 'TT' or dinucleotide == 'AA'):
			tot += 5
		elif (dinucleotide == 'TA'):
			tot += 4
	return float(tot)

def getSeq(ncbiID):
	ncbiID.replace(" ", "")
	handle = Entrez.efetch(db="nucleotide",id=ncbiID,rettype="fasta")
	record = SeqIO.read(handle, "fasta")
	handle.close()
	string = str(record.seq)
	return string

# GUI
root = Tk()
root.title("Melting Sequences")


mainframe = ttk.Frame(root, padding= "4 4 10 10")
mainframe.grid(column = 0, row = 0, sticky=(N,W,E,S))
mainframe.columnconfigure(0, weight = 1)
mainframe.rowconfigure(0, weight = 1)

ncbiID = StringVar()
strand = StringVar()
concMolar = StringVar()
concDNA = StringVar()
genome = BooleanVar()
temp = DoubleVar()

ttk.Label(mainframe, text = "Enter either the sequence or the NCBI ID to look up the genome.").grid(column=1, row=1, sticky= N)

R1 = Radiobutton(mainframe, text="NCBI ID", variable=genome, value=True).grid(column=1, row=2, sticky=W)
R2 = Radiobutton(mainframe, text="Strand sequence", variable=genome, value=False).grid(column=1, row=3, sticky=W)
strand_entry = ttk.Entry(mainframe, width = 14, textvariable = strand)
strand_entry.grid(column = 3, row = 2, sticky=(W,E))

ttk.Label(mainframe, text = "Na+ Concentration (Molar)").grid(column=1, row=4, sticky=W)
concMolar_entry = ttk.Entry(mainframe, width = 14, textvariable = concMolar)
concMolar_entry.grid(column = 3, row = 4, sticky=(W,E))

ttk.Label(mainframe, text = "Total nucleotide strand concentration (g/ml)").grid(column=1, row=5, sticky=W)
concDNA_entry = ttk.Entry(mainframe, width = 14, textvariable = concDNA)
concDNA_entry.grid(column = 3, row = 5, sticky=(W,E))

ttk.Label(mainframe, text = "Predicted temperature (°C)").grid(column=1, row=7, sticky=W)
ttk.Label(mainframe, textvariable=temp).grid(column=3, row=7, sticky = W)
string = str(strand)


ttk.Button(mainframe, text="Calculate", command=getTemp).grid(column=3, row=8, sticky=W)

for child in mainframe.winfo_children(): child.grid_configure(padx=5, pady=5)

strand_entry.focus()
root.bind('<Return>', getTemp)
root.mainloop()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
    	description='Teloreader.py <strain> <pwd> <fastafile> <size_window> <max_size_gap> <min_mean_window> <min_mean_telo> <max_dist_term>')
    parser.add_argument(
        "strain",
        help="strain name",
    )
    parser.add_argument(
        "pwd",
        help="path of the folder where the fasta file is located",
    )
    parser.add_argument(
        "fastafile",
        help="name of the fasta file of reads",
    )
    parser.add_argument(
        "-size_window",
        "-w",
        default=15,
        type=int,
        help="Size of the sliding window for the calculation of the score average (must be smaller than max_size_gap) default=15",
    )
    parser.add_argument(
        "-max_size_gap",
        "-g",
        default=20,
        type=int,
        help="Maximum size of a non telomeric gap in telomere default=20",
    )
    parser.add_argument(
        "-min_mean_window",
        "-mw",
        default=7,
        type=float,
        help="Minimum average score of a telomere window default=20",
    )
    parser.add_argument(
        "-min_mean_telo",
        "-mt",
        default=6.5,
        type=float,
        help="Minimum average score of a telomere default=6.5",
    )
    parser.add_argument(
        "-min_len",
        "-ml",
        default=16,
        type=int,
        help="Minimum length of a telomere default=16",
    )
    parser.add_argument(
        "-max_dist_term",
        "-t",
        default=50,
        help="Distance from the extremity to be considered as terminal",
    )
    parser.add_argument(
        "-out_pwd",
        "-o",
        default='None',
        help="Directory of output result",
    )

    return parser.parse_args()

def Make_reverse(seq):
	'''
	Takes sequence and returns the reverse complement 
	'''
    rev={'A':'T','T':'A','C':'G','G':'C'}
    rseq=''
    for i in range(len(seq)-1,-1,-1):
        rseq+=rev[seq[i]]
    return(rseq)


def Fasta_to_8mertype(fasta,file8mer,typ):
	'''
	Converts a fasta file to an 8-mer score file.

    :param fasta: path to input fasta file
    :param file8mer: path to output 8-mer score file
    :param typ: type of scoring system ('C' for C-rich or 'G' for G-rich)
	'''
    f=open(fasta,'r')
    f8=open(file8mer,'w')
    rest=''
    for l in f :
        if l[0]=='>':
            f8.write(l)
            rest=''
        else:
            f8.write(l)
            if l[-1]!='\n':
                f8.write('\n')
            if rest=='':
                f8.write('       ')
            if len(rest)>1:
                seq=rest+l.strip('\n')
            else:
                seq=l.strip('\n')

            for i in range(len(seq)-7):
                mer=seq[i:i+8].upper()
                if typ=='C' or typ=='c': # if C rich scoring, reverse complement the 8-mer
                    mer=Make_reverse(mer)
                f8.write(str(int(Tab_score.at[mer,'score'])))
            rest=seq[-7:]
            f8.write('\n')
    f.close()
    f8.close()


def Add_Telo(dfTelo,nam,deb,fin,typ,len_read,mean_score):
	'''
	Adds a telomere to a telomere DataFrame.

    :param dfTelo: DataFrame to add telomere to
    :param nam: name of telomere
    :param deb: start position of telomere
    :param fin: end position of telomere
    :param typ: type of telomere ('C' for C-rich or 'G' for G-rich)
    :param len_read: length of read
    :param mean_score: mean 8-mer score of telomere
    :return: updated DataFrame with new telomere
	'''
	ind=len(dfTelo)
	dfTelo.at[ind,'strain']=strain
	dfTelo.at[ind,'name']=nam
	dfTelo.at[ind,'type']=typ
	dfTelo.at[ind,'len']=fin-deb
	dfTelo.at[ind,'start']=deb+1
	dfTelo.at[ind,'end']=fin
	if deb<max_dist_term or fin+max_dist_term > len_read : 
		dfTelo.at[ind,'Loc']='term'
	else:
		dfTelo.at[ind,'Loc']='intern'
	dfTelo.at[ind,'Score_8mer']=mean_score
	dfTelo.at[ind,'reads_len']=len_read
	return(dfTelo)

def Cal_all_mean(score):
	'''
    Calculates mean scores for all windows of given size in a list of scores,
    and returns a list of decisions based on whether the mean score in each
    window is above a given threshold.

    :param score: list of integer scores
    :return: list of decisions (0 if mean score is below threshold, 1 if above)
    '''
	mean_score=[]
	for i in range(len(score)-size_window+mer_size):
		mean_score.append(np.mean([int(k) for k in score[i:i+size_window-mer_size+1]]))
	decision=[]
	for ms in mean_score:
		if ms<min_mean_window:
			decision.append(0)
		else:
			decision.append(1)
	return(decision)

def Eval_One_Read(dfTelo,nam,seq,typ,score):
	'''
    Given a score of a sequencing read, identifies and evaluates telomeric regions.
    Adds the telomeric regions to the data frame dfTelo, with the corresponding read name and type.

    :param dfTelo: Pandas DataFrame to which the telomeric regions will be added
    :param nam: read name
    :param seq: DNA sequence of the read
    :param typ: telomere type, 'G' Grich and 'C' for Crich
    :param score: list of 8mer score (in string format) of the sequence
    :returns dfTelo: the updated Pandas DataFrame with the telomeric regions added
    '''

    #Checks if 7 or 8 is in score, which may indicate the presence of a telomeric repeat.  If this is not the case, the function does not execute the rest of the code.
	if '7' in score or '8' in score :
		telo=0 
		i=0
		start_gap=''
		decision=Cal_all_mean(score)
		if typ == 'G' : # if we search a G-rich telomere, we go through the sequence from the end
			ind_seq=range(len(seq)-1,-1,-1)
			ind_score=range(len(score)-1,-1,-1)
			ind_d=range(len(decision)-1,-1,-1)
		else: # if we search a C-rich telomere, we go through the sequence from the start
			ind_seq=range(len(seq))
			ind_score=range(len(score))
			ind_d=range(len(decision))

		while i < len(score): # if we reach the end of the read
			if telo == 0 :
				if score[ind_score[i]] in  ['7','8'] : # if the quality score is 7 or 8, start a new telomeric region
					telo=1
					start=ind_score[i]
				i+=1
			else : # if in telomeric stretch
				if i >= len(decision): # if we are close to the end of the sequence
					if typ == 'G':
						score_telo = score[:start+1]
						while score_telo[0] not in  ['7','8'] or np.mean([int(k) for k in score_telo]) < min_mean_telo: 
							score_telo=score_telo[1:]
						# removes the first or/and the last nucleotide at the boundaries, when transitioning from a 8 to 7 score
						if len(score_telo)>1 and  score_telo[-1]=='7' and score_telo[-2]=='8':
							start-=1
							score_telo=score_telo[:-1]
						end=start+mer_size
						if len(score_telo)>1 and score_telo[0]=='7' and score_telo[1]=='8':
							score_telo=score_telo[1:]
						start=max(0,end - len(score_telo)- mer_size + 1)
					else : 
						score_telo = score[start:]
						while score_telo[-1] not in  ['7','8'] or np.mean([int(k) for k in score_telo]) < min_mean_telo :
							score_telo=score_telo[:-1]
						# removes the first or/and the last nucleotide at the boundaries, when transitioning from a 8 to 7 score 
						if len(score_telo)>1 and score_telo[0] =='7' and score_telo[1]=='8':
							start+=1
							score_telo=score_telo[1:]
						end=min(start+len(score_telo)-1,len(score)-1)+mer_size
						if len(score_telo)>1 and score_telo[-1]=='7' and score_telo[-2]=='8':
							score_telo=score_telo[:-1]
					if len(score_telo)+mer_size> min_len and score_telo.count('8')>=8:
						dfTelo = Add_Telo(dfTelo,nam,start,end,typ,len(seq),np.mean([int(k) for k in score_telo]))
					i=len(score) 
				elif decision[ind_d[i]] == 1: # if decision is 1, telomeric stretch continues
					i+=1
				else : # enter in a new gap 
					if start_gap != '': # if there is already a gap in the telomeric stretch
						# Verifies if the previous gap is accepted, depending on the mean of telomeric stretch up to the current position before the new gap.
						if typ == 'G':
							if np.mean([int(k) for k in score[max(ind_score[i]-size_window+mer_size,0):start+1]]) < min_mean_telo :
								telo=0
						else :
							if np.mean([int(k) for k in score[start:min(ind_score[i]+size_window-mer_size,len(score))+1]]) < min_mean_telo :
								telo=0
					if telo==1 : # if the previous gap is accepted, verifies if the new gap is tested
						start_gap=ind_score[i]
						# if the gap is too long, stops the telomere stretch.
						if 1 not in decision[min(ind_d[i],ind_d[min(i+max_size_gap+mer_size,len(ind_d)-1)]):max(ind_d[i],ind_d[min(i+max_size_gap+mer_size,len(ind_d)-1)])+1] : 
							telo=0
						else: # telomeric stretch continues
							while decision[ind_d[i]]!=1:
								i+=1
				if telo == 0 :
					if typ == 'G':
						score_telo = score[max(0,start_gap-size_window+mer_size):start+1]
						while score_telo[0] not in  ['7','8'] or np.mean([int(k) for k in score_telo]) < min_mean_telo: 
							score_telo=score_telo[1:]
						# removes the first or/and the last nucleotide at the boundaries, when transitioning from a 8 to 7 score 
						if len(score_telo)>1 and  score_telo[-1]=='7' and score_telo[-2]=='8':
							start-=1
							score_telo=score_telo[:-1]
						if len(score_telo)>1 and score_telo[0]=='7' and score_telo[1]=='8':
							score_telo=score_telo[1:]
						end=start+mer_size
						start=max(0,end - len(score_telo)- mer_size + 1)
						
						i=ind_seq.index(max(start-1,0))

					else : 

						score_telo = score[start:min(start_gap+size_window-mer_size,len(score))+1]
						j=min(len(score_telo)-(start_gap-start),len(score_telo))+1

						while score_telo[-1] not in  ['7','8'] or np.mean([int(k) for k in score_telo]) < min_mean_telo :
							score_telo=score_telo[:-1]


						if len(score_telo)>1 and score_telo[0] =='7' and score_telo[1]=='8':
							start+=1
							score_telo=score_telo[1:]
						if len(score_telo)>1 and score_telo[-1]=='7' and score_telo[-2]=='8':
							score_telo=score_telo[:-1]
						end=min(start+len(score_telo),len(score)-1)+mer_size-1
						i=ind_seq.index(end)
					start_gap=''
					# Check if there are at least three scores of 8 in the telomere
					if len(score_telo)+mer_size> min_len and score_telo.count('8')>=8 :
						dfTelo = Add_Telo(dfTelo,nam,start,end,typ,len(seq),np.mean([int(k) for k in score_telo]))
	return(dfTelo)

def Add_N(dfTelo):
	'''
	Adds name of read in the dfTelo DataFrame

	:param dfTelo: DataFrame to add name
	:return: updated DataFrame with read name of telomere
	'''
	for i in dfTelo.index:
		dfTelo.at[i,'N']=len(dfTelo[dfTelo['name']==dfTelo.at[i,'name']])
	return(dfTelo)

def Make_Fasta_Telo_Complet_dfTelo(csv,fasta,original_fasta):
	'''
	Creates a fasta file with all telomeric regions and adds telomeric lengths in csv. 

	:param csv: DataFrame (dfTelo) with all telomeric region found
	:fasta: path of the new created fasta file with all telomeric regions
	:original_fasta: path of original fasta file with all sequences
	:return :updated DataFrame with read lengths 
	'''
    foriginal=open(original_fasta,'r')
    f=open(fasta,'w')
    seq=''
    ok=0
    for l in foriginal:
        if l[0]=='>':
            if ok==1:
                for i in ind:
                    s=int(csv.at[i,'start'])-1
                    e=int(csv.at[i,'end'])
                    csv.at[i,'reads_len']=len(seq)
                    f.write('>'+nam+' Loc : '+csv.at[i,'Loc']+' Telo start: '+str(csv.at[i,'start'])+' Telo end: '+str(csv.at[i,'end'])+'\n')
                    f.write(seq[s:e]+'\n')
            nam=l[1:-1]
            if nam not in list(csv['name']):
                ok=0
            else:
                ok=1
                ind=list(csv[csv['name']==nam].index)
            seq=''
        elif ok==1:
            seq+=(l[:-1])
            
    if ok==1:
        for i in ind:
            s=int(csv.at[i,'start'])-1
            e=int(csv.at[i,'end'])
            csv.at[i,'reads_len']=len(seq)
            f.write('>'+nam+' Loc : '+csv.at[i,'Loc']+' Telo start: '+str(csv.at[i,'start'])+' Telo end: '+str(csv.at[i,'end'])+'\n')
            f.write(seq[s:e]+'\n')
    f.close()
    foriginal.close()
    return(csv)

def Find_Telo_on_8mer(original_fasta,out_fasta,fileC8mer,fileG8mer,out_csv):
	'''
	This function initializes dfTelo DataFrame, and executes the different functions to obtain output results.

	:param original_fasta: path of original fasta file with all sequences
	:param out_fasta: path of the new created fasta file with all telomeric regions
	:param fileC8mer: path of the file created with the function Fasta_to_8mertype with C-rich 8mer scores of given sequences.
	:param fileG8mer: path of the file created with the function Fasta_to_8mertype with G-rich 8mer scores of given sequences.
	:param out_csv: path of the out Dataframe summarizing all telomeric regions found.
	'''
	dfTelo=pd.DataFrame(columns=['strain','name','N','type','len','start','end','Loc','Score_8mer','reads_len'])
	for file,typ in ((fileC8mer,'C'),(fileG8mer,'G')):
		f = open(file,'r')
		seq=''
		score=''
		for l in f :
			if l[0]=='>':
				if seq != '':
					dfTelo=Eval_One_Read(dfTelo,nam,seq,typ,score)
				nam = l[1:-1]
				seq=''
				score=''
			elif l.split(' ')[-1][0] in ['4','5','6','7','8']:
				score+=l.split(' ')[-1].strip('\n')
			else : 
				seq+=l[:-1]
		dfTelo=Eval_One_Read(dfTelo,nam,seq,typ,score)
	dfTelo=Add_N(dfTelo)
	dfTelo.to_csv(out_csv,sep='\t')
	Make_Fasta_Telo_Complet_dfTelo(dfTelo,out_fasta,original_fasta)


if __name__ == "__main__":
	args = parse_arguments()

	print("strain : ",args.strain, 
	"\n pwd : ",args.pwd,
	"\n fastafile : ",args.fastafile,
	"\n size_window : ",args.size_window,
	"\n max_size_gap : ",args.max_size_gap,
	"\n min_mean_window : ",args.min_mean_window,
	"\n min_mean_telo : ",args.min_mean_telo,
	"\n max_dist_term : ",args.max_dist_term,
	"\n min_len : ",args.min_len,
	"\n out_pwd : ",args.out_pwd,'\n')

	if args.out_pwd == 'None':
		out_pwd=args.pwd
	else:
		out_pwd=args.out_pwd
	strain=args.strain
	min_len=args.min_len
	original_fasta=args.pwd+args.fastafile
	fileC8mer=args.pwd+'C8mer_'+args.fastafile
	fileG8mer=args.pwd+'G8mer_'+args.fastafile
	out_fasta=out_pwd+'OUTMer_'+args.fastafile
	out_csv=out_pwd+'OUTMer_'+args.fastafile.rsplit('.',1)[0]+'.csv'


	size_window = args.size_window 
	max_size_gap = args.max_size_gap 
	min_mean_window = args.min_mean_window 
	min_mean_telo = args.min_mean_telo 
	max_dist_term = args.max_dist_term 
	mer_size=8

	Tab_score=pd.read_csv('Score_for_all_8Mer_min4.tab',sep='\t',index_col='Mer')

	Fasta_to_8mertype(original_fasta,fileC8mer,'C')
	Fasta_to_8mertype(original_fasta,fileG8mer,'G')
	Find_Telo_on_8mer(original_fasta,out_fasta,fileC8mer,fileG8mer,out_csv)






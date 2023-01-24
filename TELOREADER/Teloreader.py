#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
    	description='TeloFind_8mer.py <strain> <pwd> <fastafile> <size_window> <max_size_gap> <min_mean_window> <min_mean_telo> <max_dist_term>')
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
        help="Taille maximal d'un gap non télomérique dans le télomère default=20",
    )
    parser.add_argument(
        "-min_mean_window",
        "-mw",
        default=7,
        type=float,
        help="Maximum size of a non telomeric gap in the telomere default=20",
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
    rev={'A':'T','T':'A','C':'G','G':'C'}
    rseq=''
    for i in range(len(seq)-1,-1,-1):
        rseq+=rev[seq[i]]
    return(rseq)


def Fasta_to_8mertype(fasta,file8mer,typ):
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
                if typ=='C' or typ=='c':
                    mer=Make_reverse(mer)
                f8.write(str(int(Tab_score.at[mer,'score'])))
            rest=seq[-7:]
            f8.write('\n')


def Add_Telo(dfTelo,nam,deb,fin,typ,len_read,mean_score):
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
	if '7' in score or '8' in score :
		telo=0 
		i=0
		deb_gap=''
		decision=Cal_all_mean(score)
		if typ == 'G' : 
			ind_seq=range(len(seq)-1,-1,-1)
			ind_score=range(len(score)-1,-1,-1)
			ind_d=range(len(decision)-1,-1,-1)
		else:
			ind_seq=range(len(seq))
			ind_score=range(len(score))
			ind_d=range(len(decision))

		while i < len(score):
			if telo == 0 :
				if score[ind_score[i]] in  ['7','8'] :
					telo=1
					deb=ind_score[i]
				i+=1
			else :
				if i >= len(decision):
					if typ == 'G':
						score_telo = score[:deb+1]
						while score_telo[0] not in  ['7','8'] or np.mean([int(k) for k in score_telo]) < min_mean_telo: 
							score_telo=score_telo[1:]
						if len(score_telo)>1 and  score_telo[-1]=='7' and score_telo[-2]=='8':
							deb-=1
							score_telo=score_telo[:-1]
						fin=deb+mer_size
						if len(score_telo)>1 and score_telo[0]=='7' and score_telo[1]=='8':
							score_telo=score_telo[1:]
						deb=max(0,fin - len(score_telo)- mer_size + 1)
					else : 
						score_telo = score[deb:]
						while score_telo[-1] not in  ['7','8'] or np.mean([int(k) for k in score_telo]) < min_mean_telo :
							score_telo=score_telo[:-1]




						if len(score_telo)>1 and score_telo[0] =='7' and score_telo[1]=='8':
							deb+=1
							score_telo=score_telo[1:]
						fin=min(deb+len(score_telo)-1,len(score)-1)+mer_size
						if len(score_telo)>1 and score_telo[-1]=='7' and score_telo[-2]=='8':
							score_telo=score_telo[:-1]
					if len(score_telo)+mer_size> min_len and score_telo.count('8')>=8:
						dfTelo = Add_Telo(dfTelo,nam,deb,fin,typ,len(seq),np.mean([int(k) for k in score_telo]))
					i=len(score) 
				elif decision[ind_d[i]] == 1:
					i+=1
				else : 
					if deb_gap != '': 
						if typ == 'G':
							if np.mean([int(k) for k in score[max(ind_score[i]-size_window+mer_size,0):deb+1]]) < min_mean_telo :
								telo=0
						else :
							if np.mean([int(k) for k in score[deb:min(ind_score[i]+size_window-mer_size,len(score))+1]]) < min_mean_telo :
								telo=0
					if telo==1 :
						deb_gap=ind_score[i]
						if 1 not in decision[min(ind_d[i],ind_d[min(i+max_size_gap+mer_size,len(ind_d)-1)]):max(ind_d[i],ind_d[min(i+max_size_gap+mer_size,len(ind_d)-1)])+1] : # le gap est trop long
							telo=0
						else:
							while decision[ind_d[i]]!=1:
								i+=1
				if telo == 0 :
					if typ == 'G':
						score_telo = score[max(0,deb_gap-size_window+mer_size):deb+1]
						while score_telo[0] not in  ['7','8'] or np.mean([int(k) for k in score_telo]) < min_mean_telo: 
							score_telo=score_telo[1:]

						if len(score_telo)>1 and  score_telo[-1]=='7' and score_telo[-2]=='8':
							deb-=1
							score_telo=score_telo[:-1]
						if len(score_telo)>1 and score_telo[0]=='7' and score_telo[1]=='8':
							score_telo=score_telo[1:]
						fin=deb+mer_size
						deb=max(0,fin - len(score_telo)- mer_size + 1)
						
						i=ind_seq.index(max(deb-1,0))

					else : 

						score_telo = score[deb:min(deb_gap+size_window-mer_size,len(score))+1]
						j=min(len(score_telo)-(deb_gap-deb),len(score_telo))+1

						while score_telo[-1] not in  ['7','8'] or np.mean([int(k) for k in score_telo]) < min_mean_telo :
							score_telo=score_telo[:-1]


						if len(score_telo)>1 and score_telo[0] =='7' and score_telo[1]=='8':
							deb+=1
							score_telo=score_telo[1:]
						if len(score_telo)>1 and score_telo[-1]=='7' and score_telo[-2]=='8':
							score_telo=score_telo[:-1]
						fin=min(deb+len(score_telo),len(score)-1)+mer_size-1
						i=ind_seq.index(fin)
					deb_gap=''

					if len(score_telo)+mer_size> min_len and score_telo.count('8')>=8 :
						dfTelo = Add_Telo(dfTelo,nam,deb,fin,typ,len(seq),np.mean([int(k) for k in score_telo]))
	return(dfTelo)

def Add_N(dfTelo):
	for i in dfTelo.index:
		dfTelo.at[i,'N']=len(dfTelo[dfTelo['name']==dfTelo.at[i,'name']])
	return(dfTelo)

def Make_Fasta_Telo_Complet_dfTelo(csv,fasta,original_fasta):
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





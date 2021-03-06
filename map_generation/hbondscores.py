# Hydrogen bonding scoring function for surface points - by base and atom
# returns a tuple (a, b) where a is 0, 1, 2 or 3 hydrogen donation capability;
# b is 0, 1 or 2 depending on H-bond acceptance.


def bondscore(a):
	return (float(score_table[a][0]),float(score_table[a][1]))


score_table = {
'ALAC': (0, 0),
'ALACA': (0, 0),
'ALACB': (0, 0),
'ALAH': (0, 0),
'ALAHA': (0, 0),
'ALAHB1': (0, 0),
'ALAHB2': (0, 0),
'ALAHB3': (0, 0),
'ALAN': (1, 0),
'ALAO': (0, 2),
'ALAOXT':(0, 2),
'ARGC': (0, 0),
'ARGCA': (0, 0),
'ARGCB': (0, 0),
'ARGCD': (0, 0),
'ARGCG': (0, 0),
'ARGCZ': (0, 0),
'ARGH': (0, 0),
'ARGHA': (0, 0),
'ARGHB2': (0, 0),
'ARGHB3': (0, 0),
'ARGHD2': (0, 0),
'ARGHD3': (0, 0),
'ARGHE': (0, 0),
'ARGHG2': (0, 0),
'ARGHG3': (0, 0),
'ARGHH11': (0, 0),
'ARGHH12': (0, 0),
'ARGHH21': (0, 0),
'ARGHH22': (0, 0),
'ARGN': (1, 0),
'ARGNE': (1, 0),
'ARGNH1': (2, 0),
'ARGNH2': (2, 0),
'ARGO': (0, 2),
'ARGOXT':(0, 2),
'ASNC': (0, 0),
'ASNCA': (0, 0),
'ASNCB': (0, 0),
'ASNCG': (0, 0),
'ASNH': (0, 0),
'ASNHA': (0, 0),
'ASNHB2': (0, 0),
'ASNHB3': (0, 0),
'ASNHD21': (0, 0),
'ASNHD22': (0, 0),
'ASNN': (1, 0),
'ASNND2': (2, 0),
'ASNO': (0, 2),
'ASNOXT':(0, 2),
'ASNOD1': (0, 2),
'ASPC': (0, 0),
'ASPCA': (0, 0),
'ASPCB': (0, 0),
'ASPCG': (0, 0),
'ASPH': (0, 0),
'ASPHA': (0, 0),
'ASPHB2': (0, 0),
'ASPHB3': (0, 0),
'ASPN': (1, 0),
'ASPO': (0, 2),
'ASPOXT':(0, 2),
'ASPOD1': (0, 2),
'ASPOD2': (0, 2),
'ASXN' : (1, 0),
'ASXCA' : (0, 0),
'ASXC' : (0, 0),
'ASXO' : (0, 2),
'ASXCB' : (0, 0),
'ASXCG' : (0, 0),
'ASXXD1' : (0, 2),
'ASXXD2' : (1, 1), #average of scores for ASPOD2 and ASNND2
'CYSC': (0, 0),
'CYSCA': (0, 0),
'CYSCB': (0, 0),
'CYSH': (0, 0),
'CYSHA': (0, 0),
'CYSHB2': (0, 0),
'CYSHB3': (0, 0),
'CYSHG': (0, 0),
'CYSN': (1, 0),
'CYSO': (0, 2),
'CYSOXT':(0, 2),
'CYSSG': (0, 0),
'GLNC': (0, 0),
'GLNCA': (0, 0),
'GLNCB': (0, 0),
'GLNCD': (0, 0),
'GLNCG': (0, 0),
'GLNH': (0, 0),
'GLNHA': (0, 0),
'GLNHB2': (0, 0),
'GLNHB3': (0, 0),
'GLNHE21': (0, 0),
'GLNHE22': (0, 0),
'GLNHG2': (0, 0),
'GLNHG3': (0, 0),
'GLNN': (1, 0),
'GLNNE2': (2, 0),
'GLNO': (0, 2),
'GLNOE1': (0, 2),
'GLNOXT': (0, 2),
'GLUC': (0, 0),
'GLUCA': (0, 0),
'GLUCB': (0, 0),
'GLUCD': (0, 0),
'GLUCG': (0, 0),
'GLUH': (0, 0),
'GLUHA': (0, 0),
'GLUHB2': (0, 0),
'GLUHB3': (0, 0),
'GLUHG2': (0, 0),
'GLUHG3': (0, 0),
'GLUN': (1, 0),
'GLUO': (0, 2),
'GLUOXT':(0, 2),
'GLUOE1': (0, 2),
'GLUOE2': (0, 2),
'GLYC': (0, 0),
'GLYCA': (0, 0),
'GLYH': (0, 0),
'GLYHA2': (0, 0),
'GLYHA3': (0, 0),
'GLYN': (1, 0),
'GLYO': (0, 2),
'GLYOXT':(0, 2),
'HISC': (0, 0),
'HISCA': (0, 0),
'HISCB': (0, 0),
'HISCD2': (0, 0),
'HISCE1': (0, 0),
'HISCG': (0, 0),
'HISH': (0, 0),
'HISHA': (0, 0),
'HISHB2': (0, 0),
'HISHB3': (0, 0),
'HISHD1': (0, 0),
'HISHD2': (0, 0),
'HISHE1': (0, 0),
'HISHE2': (0, 0),
'HISN': (1, 0),
'HISND1': (1, 0),
'HISNE2': (0, 1),
'HISO': (0, 2),
'HISOXT':(0, 2),
'ILEC': (0, 0),
'ILECA': (0, 0),
'ILECB': (0, 0),
'ILECD1': (0, 0),
'ILECG1': (0, 0),
'ILECG2': (0, 0),
'ILEH': (0, 0),
'ILEHA': (0, 0),
'ILEHB': (0, 0),
'ILEHD11': (0, 0),
'ILEHD12': (0, 0),
'ILEHD13': (0, 0),
'ILEHG12': (0, 0),
'ILEHG13': (0, 0),
'ILEHG21': (0, 0),
'ILEHG22': (0, 0),
'ILEHG23': (0, 0),
'ILEN': (1, 0),
'ILEO': (0, 2),
'ILEOXT':(0, 2),
'LEUC': (0, 0),
'LEUCA': (0, 0),
'LEUCB': (0, 0),
'LEUCD1': (0, 0),
'LEUCD2': (0, 0),
'LEUCG': (0, 0),
'LEUH': (0, 0),
'LEUHA': (0, 0),
'LEUHB2': (0, 0),
'LEUHB3': (0, 0),
'LEUHD11': (0, 0),
'LEUHD12': (0, 0),
'LEUHD13': (0, 0),
'LEUHD21': (0, 0),
'LEUHD22': (0, 0),
'LEUHD23': (0, 0),
'LEUHG': (0, 0),
'LEUN': (1, 0),
'LEUO': (0, 2),
'LEUOXT':(0, 2),
'LYSC': (0, 0),
'LYSCA': (0, 0),
'LYSCB': (0, 0),
'LYSCD': (0, 0),
'LYSCE': (0, 0),
'LYSCG': (0, 0),
'LYSH': (0, 0),
'LYSHA': (0, 0),
'LYSHB2': (0, 0),
'LYSHB3': (0, 0),
'LYSHD2': (0, 0),
'LYSHD3': (0, 0),
'LYSHE2': (0, 0),
'LYSHE3': (0, 0),
'LYSHG2': (0, 0),
'LYSHG3': (0, 0),
'LYSHZ1': (0, 0),
'LYSHZ2': (0, 0),
'LYSHZ3': (0, 0),
'LYSN': (1, 0),
'LYSNZ': (3, 0),
'LYSO': (0, 2),
'LYSOXT':(0, 2),
'METC': (0, 0),
'METCA': (0, 0),
'METCB': (0, 0),
'METCE': (0, 0),
'METCG': (0, 0),
'METH': (0, 0),
'METHA': (0, 0),
'METHB2': (0, 0),
'METHB3': (0, 0),
'METHE1': (0, 0),
'METHE2': (0, 0),
'METHE3': (0, 0),
'METHG2': (0, 0),
'METHG3': (0, 0),
'METN': (1, 0),
'METO': (0, 2),
'METOXT':(0, 2),
'METSD': (0, 0),
'PHEC': (0, 0),
'PHECA': (0, 0),
'PHECB': (0, 0),
'PHECD1': (0, 0),
'PHECD2': (0, 0),
'PHECE1': (0, 0),
'PHECE2': (0, 0),
'PHECG': (0, 0),
'PHECZ': (0, 0),
'PHEH': (0, 0),
'PHEHA': (0, 0),
'PHEHB2': (0, 0),
'PHEHB3': (0, 0),
'PHEHD1': (0, 0),
'PHEHD2': (0, 0),
'PHEHE1': (0, 0),
'PHEHE2': (0, 0),
'PHEHZ': (0, 0),
'PHEN': (1, 0),
'PHEO': (0, 2),
'PHEOXT':(0, 2),
'PROC': (0, 0),
'PROCA': (0, 0),
'PROCB': (0, 0),
'PROCD': (0, 0),
'PROCG': (0, 0),
'PROHA': (0, 0),
'PROHB2': (0, 0),
'PROHB3': (0, 0),
'PROHD2': (0, 0),
'PROHD3': (0, 0),
'PROHG2': (0, 0),
'PROHG3': (0, 0),
'PRON': (1, 0),
'PROO': (0, 2),
'PROOXT':(0, 2),
'SERC': (0, 0),
'SERCA': (0, 0),
'SERCB': (0, 0),
'SERH': (0, 0),
'SERHA': (0, 0),
'SERHB2': (0, 0),
'SERHB3': (0, 0),
'SERHG': (0, 0),
'SERN': (1, 0),
'SERO': (0, 2),
'SEROXT':(0, 2),
'SEROG': (1, 2),
'THRC': (0, 0),
'THRCA': (0, 0),
'THRCB': (0, 0),
'THRCG2': (0, 0),
'THRH': (0, 0),
'THRHA': (0, 0),
'THRHB': (0, 0),
'THRHG1': (0, 0),
'THRHG21': (0, 0),
'THRHG22': (0, 0),
'THRHG23': (0, 0),
'THRN': (1, 0),
'THRO': (0, 2),
'THROXT':(0, 2),
'THROG1': (1, 2),
'TRPC': (0, 0),
'TRPCA': (0, 0),
'TRPCB': (0, 0),
'TRPCD1': (0, 0),
'TRPCD2': (0, 0),
'TRPCE2': (0, 0),
'TRPCE3': (0, 0),
'TRPCG': (0, 0),
'TRPCH2': (0, 0),
'TRPCZ2': (0, 0),
'TRPCZ3': (0, 0),
'TRPH': (0, 0),
'TRPHA': (0, 0),
'TRPHB2': (0, 0),
'TRPHB3': (0, 0),
'TRPHD1': (0, 0),
'TRPHE1': (0, 0),
'TRPHE3': (0, 0),
'TRPHH2': (0, 0),
'TRPHZ2': (0, 0),
'TRPHZ3': (0, 0),
'TRPN': (1, 0),
'TRPNE1': (1, 0),
'TRPO': (0, 2),
'TRPOXT':(0, 2),
'TYRC': (0, 0),
'TYRCA': (0, 0),
'TYRCB': (0, 0),
'TYRCD1': (0, 0),
'TYRCD2': (0, 0),
'TYRCE1': (0, 0),
'TYRCE2': (0, 0),
'TYRCG': (0, 0),
'TYRCZ': (0, 0),
'TYRH': (0, 0),
'TYRHA': (0, 0),
'TYRHB2': (0, 0),
'TYRHB3': (0, 0),
'TYRHD1': (0, 0),
'TYRHD2': (0, 0),
'TYRHE1': (0, 0),
'TYRHE2': (0, 0),
'TYRHH': (0, 0),
'TYRN': (1, 0),
'TYRO': (0, 2),
'TYROXT':(0, 2),
'TYROH': (1, 1),
'VALC': (0, 0),
'VALCA': (0, 0),
'VALCB': (0, 0),
'VALCG1': (0, 0),
'VALCG2': (0, 0),
'VALH': (0, 0),
'VALHA': (0, 0),
'VALHB': (0, 0),
'VALHG11': (0, 0),
'VALHG12': (0, 0),
'VALHG13': (0, 0),
'VALHG21': (0, 0),
'VALHG22': (0, 0),
'VALHG23': (0, 0),
'VALN': (1, 0),
'VALO': (0, 2),
'VALOXT':(0, 2)
}

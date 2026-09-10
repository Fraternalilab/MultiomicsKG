HPRDPPIpool = set()
HPRDUniprotID = dict()
HPRDmultiUniprotID = dict()


def findHPRDUniprot(finName):
    import csv
    import sys

    csv.field_size_limit(sys.maxsize)
    with open(finName) as csvfile:
        protIDreader = csv.reader(csvfile, delimiter="\t")
        for row in protIDreader:
            if not row[6] == None:
                uniACs = row[6].split(",")
                uniIDs = list()
                for anUniAC in uniACs:
                    tUniID = findUniprotID(anUniAC)
                    if tUniID != "NA":
                        if not (tUniID.startswith(anUniAC)):
                            uniIDs.append(tUniID)
                if len(uniIDs) > 1:
                    # print(uniIDs)
                    for aID in uniIDs:
                        HPRDmultiUniprotID[aID] = row[0]
                elif len(uniIDs) == 1:
                    HPRDUniprotID[row[0]] = uniIDs[0]


# 				else :
# 					HPRDUniprotID[row[0]]='NA'


def findMultiUniIDs(hprdid):
    multiUniIDs = list()
    for k, v in HPRDmultiUniprotID.items():
        if hprdid == v:
            multiUniIDs.append(k)
    return multiUniIDs


def buildHPRDPPI(finName):
    import csv
    import sys

    csv.field_size_limit(sys.maxsize)
    with open(finName) as csvfile:
        protppireader = csv.reader(csvfile, delimiter="\t")
        for row in protppireader:
            if row[1] in HPRDUniprotID.keys():
                prot1_UniID = HPRDUniprotID[row[1]]

                if row[4] in HPRDUniprotID.keys():
                    prot2_UniID = HPRDUniprotID[row[4]]
                    apair = frozenset([prot1_UniID, prot2_UniID])
                    if len(apair) == 2:
                        HPRDPPIpool.add(apair)
                elif row[4] in HPRDmultiUniprotID.values():
                    prot2IDs = findMultiUniIDs(row[4])
                    for aprot2ID in prot2IDs:
                        apair = frozenset([prot1_UniID, aprot2ID])
                        if len(apair) == 2:
                            HPRDPPIpool.add(apair)
            elif row[1] in HPRDmultiUniprotID.values():
                prot1IDs = findMultiUniIDs(row[1])
                if row[4] in HPRDUniprotID.keys():
                    prot2_UniID = HPRDUniprotID[row[4]]
                    for aprot1ID in prot1IDs:
                        apair = frozenset([aprot1ID, prot2_UniID])
                        if len(apair) == 2:
                            HPRDPPIpool.add(apair)
                elif row[4] in HPRDmultiUniprotID.values():
                    prot2IDs = findMultiUniIDs(row[4])
                    for aprot1ID in prot1IDs:
                        for aprot2ID in prot2IDs:
                            apair = frozenset([aprot1ID, aprot2ID])
                            if len(apair) == 2:
                                HPRDPPIpool.add(apair)

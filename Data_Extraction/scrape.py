import json
from collections import Counter

with open("uniprot_proteins.txt", "r") as f:
    proteins = f.read().splitlines()

print("Total Uniprot proteins: ", len(proteins))

my_count = 0
my_proteins_uniprot = []
with open("CD_data.txt") as f:
    entire_file = f.read()
    for protein in proteins:
        if protein in entire_file:
            my_count += 1
            my_proteins_uniprot.append(protein)

print("Total Uniprot proteins identified by me: ", my_count)
print(my_proteins_uniprot)

#--------------glycosid extraction-------------
with open("glycosid_data.txt") as file:
    content = file.read()

seperation = content.split("Clan: ")
limit = len(seperation)

z = 0
uniprot_ids = []
glycosid_clan = []
for protein in my_proteins_uniprot:
    for i in range(limit):
        if protein in seperation[i]:
            print(True)
            uniprot_ids.append(protein)
            l = seperation[i].split(" ")
            okk = l[0].split("\n")
            print(f"Protein : {protein} => Glycosid_Clan : {okk}")
            glycosid_clan.append(okk[0])
            z = z+ 1
print(f"No. of proteins extracted : {z}")
print(uniprot_ids)
print(glycosid_clan)


with open("glycosid_data.txt") as file:
    data = file.read()

family_seperation = data.split("Family ")
j = family_seperation[-1].split("Non-classified enzymes")
family_seperation[-1] = j[0]
print(family_seperation[-2])
family_seperation.append(j[1])
limit = len(family_seperation)

g = 0
glycosid_family = []
for protein in my_proteins_uniprot:
    for i in range(limit):
        if protein in family_seperation[i]:
            l = family_seperation[i].split(" ")
            okk = l[0].split("\n")
            if okk[0] == "":
                print(True)
                glycosid_family.append("Non-Classified Enzymes")
            else:
                glycosid_family.append(okk[0])
            print(f"Protein : {protein} => Glycosid_Family : {okk[0]}")

            g = g + 1
print(glycosid_family)
glycosid_dict = {}
for i in range(len(my_proteins_uniprot)):
    protein = my_proteins_uniprot[i]
    glycosid_dict[protein] = {"clan" : glycosid_clan[i], "family" : glycosid_family[i]}

with open("glycosid_dict.json", "w") as f:
    json.dump(glycosid_dict, f, indent=4)
#
#--------------kinases extraction---------------

with open('kinases_data.txt', 'r') as file:
    splitted_parts = []
    temp_part = []

    for line in file:
        if line.startswith('='):
            splitted_parts.append(temp_part)
            temp_part = []
        temp_part.append(line.strip())

    if temp_part:
        splitted_parts.append(temp_part)


j = 0

kinases_clan = []
kinases_families = []
for protein in my_proteins_uniprot:
    for i in range(1,len(splitted_parts)):
        kinases_family = splitted_parts[i][1]
        for okk in splitted_parts[i]:
            if protein in okk:
                now = okk.split(" ")
                print(f"Protein : {protein} => Kinases_Clan : {now[0]}, Kinases_Family : {kinases_family}")
                kinases_clan.append(now[0])
                kinases_families.append(kinases_family)
                j = j+1
print(j)
print(len(kinases_families))
print(len(kinases_clan))
kinases_dict = {}
for i in range(len(my_proteins_uniprot)):
    protein = my_proteins_uniprot[i]
    kinases_dict[protein] = {"clan" : kinases_clan[i], "family" : kinases_families[i]}

with open("kinases_dict.json", "w") as f:
    json.dump(kinases_dict, f, indent=4)

#------hox_extraction-------
with open('hox_proteins_data.txt', 'r') as file:
    splitted_parts = []
    temp_part = []

    for line in file:
        if not line.startswith(' '):
            splitted_parts.append(temp_part)
            temp_part = []
        temp_part.append(line.strip())

    if temp_part:
        splitted_parts.append(temp_part)

count = 0
hox_families = []
for protein in my_proteins_uniprot:
    for l in splitted_parts:
        if len(l) !=0:
            extract = l[0]
            now = extract.split(" ")
        for j in l:
            if protein in j:
                print(f"Protein : {protein} => hox_family : {now[0]}")
                hox_families.append(now[0])
                count += 1
print(count)

hox_dict = {}
for i in range(len(my_proteins_uniprot)):
    protein = my_proteins_uniprot[i]
    hox_dict[protein] = hox_families[i]

with open("hox_proteins_dict.json", "w") as f:
    json.dump(hox_dict, f, indent=4)

#--------peptides_extraction-----------

with open('peptides_data.txt', 'r') as file:
    splitted_parts = []
    temp_part = []

    for line in file:
        if line.startswith('-'):
            splitted_parts.append(temp_part)
            temp_part = []
        temp_part.append(line.strip())

    if temp_part:
        splitted_parts.append(temp_part)

count = 0
pept = []
for protein in my_proteins_uniprot:
    for i in range(1,len(splitted_parts)):
        peptides = splitted_parts[i][1]
        for j in splitted_parts[i]:
            if protein in j:
                print(f"Protein : {protein} => peptides : {peptides}")
                pept.append(peptides)
                count += 1
print(count)

peptide_dict = {}
for i in range(len(my_proteins_uniprot)):
    protein = my_proteins_uniprot[i]
    peptide_dict[protein] = pept[i]

with open("peptides_dict.json", "w") as f:
    json.dump(peptide_dict, f, indent=4)

#----------------ribome_extraction---------

with open('ribosomal_proteins_data.txt', 'r') as file:
    splitted_parts = []
    temp_part = []

    for line in file:
        if line.startswith('\n'):
            splitted_parts.append(temp_part)
            temp_part = []
        temp_part.append(line.strip())

    if temp_part:
        splitted_parts.append(temp_part)

count = 0
ribo = []
for protein in my_proteins_uniprot:
    for i in splitted_parts:
        for okk in i:
            if protein in okk:
                print(f"Protein : {protein} => ribosomes : {i[1]}")
                ribo.append(i[1])

                count +=1
print(count)

ribo_dict = {}
for i in range(len(my_proteins_uniprot)):
    protein = my_proteins_uniprot[i]
    ribo_dict[protein] = ribo[i]

with open("ribosomal_proteins_dict.json", "w") as f:
    json.dump(ribo_dict, f, indent=4)


#---------CD_Extraction---------------
listing = ['A6NI73', 'A8K4G0', 'O00144', 'O00206', 'O00220', 'O00238', 'O00241', 'O00481', 'O00574', 'O14638', 'O14672', 'O14763', 'O14786', 'O14788', 'O14798', 'O14836', 'O14931', 'O15389', 'O15455', 'O43490', 'O43557', 'O43699', 'O60449', 'O60486', 'O60603', 'O75015', 'O75019', 'O75022', 'O75023', 'O75144', 'O75326', 'O75330', 'O75509', 'O75888', 'O76036', 'O95256', 'O95727', 'O95944', 'O95971', 'P01589', 'P01730', 'P01732', 'P02724', 'P02730', 'P02786', 'P04156', 'P04216', 'P04233', 'P04234', 'P04626', 'P04921', 'P05106', 'P05107', 'P05362', 'P05556', 'P06028', 'P06126', 'P06127', 'P06213', 'P06729', 'P06731', 'P06734', 'P06756', 'P07204', 'P07333', 'P07359', 'P07766', 'P08069', 'P08138', 'P08174', 'P08183', 'P08195', 'P08473', 'P08514', 'P08571', 'P08575', 'P08582', 'P08637', 'P08648', 'P08887', 'P08962', 'P09326', 'P09564', 'P09619', 'P09693', 'P10721', 'P10747', 'P10966', 'P11049', 'P11215', 'P11279', 'P11362', 'P11464', 'P11717', 'P11836', 'P11912', 'P12018', 'P12314', 'P12318', 'P12821', 'P12830', 'P13164', 'P13224', 'P13473', 'P13591', 'P13598', 'P13612', 'P13688', 'P13726', 'P13987', 'P14151', 'P14209', 'P14770', 'P14778', 'P14784', 'P15144', 'P15151', 'P15260', 'P15391', 'P15509', 'P15529', 'P15812', 'P15813', 'P15814', 'P15941', 'P16070', 'P16109', 'P16144', 'P16150', 'P16234', 'P16284', 'P16410', 'P16422', 'P16581', 'P16671', 'P16871', 'P17301', 'P17813', 'P17927', 'P18577', 'P18627', 'P18827', 'P19022', 'P19256', 'P19320', 'P19397', 'P19438', 'P19440', 'P20023', 'P20138', 'P20273', 'P20333', 'P20701', 'P20702', 'P20963', 'P21453', 'P21589', 'P21730', 'P21757', 'P21802', 'P21854', 'P21926', 'P22413', 'P22455', 'P22607', 'P22897', 'P23229', 'P23276', 'P23510', 'P24071', 'P24394', 'P25024', 'P25025', 'P25063', 'P25445', 'P25942', 'P26006', 'P26715', 'P26717', 'P26718', 'P26842', 'P26951', 'P27487', 'P27701', 'P27930', 'P28906', 'P28907', 'P28908', 'P29016', 'P29017', 'P29965', 'P30203', 'P31358', 'P31785', 'P31994', 'P31995', 'P31997', 'P32004', 'P32246', 'P32248', 'P32302', 'P32927', 'P32942', 'P32970', 'P32971', 'P33151', 'P33681', 'P34741', 'P34810', 'P34910', 'P35613', 'P35968', 'P36888', 'P36894', 'P38570', 'P40189', 'P40197', 'P40198', 'P40199', 'P40200', 'P40238', 'P40259', 'P41217', 'P41597', 'P41732', 'P42081', 'P42701', 'P42702', 'P43121', 'P43489', 'P43626', 'P43628', 'P43629', 'P43630', 'P43631', 'P43632', 'P48023', 'P48357', 'P48509', 'P48960', 'P49682', 'P49961', 'P50591', 'P50895', 'P51677', 'P51679', 'P51681', 'P51684', 'P51685', 'P51686', 'P52961', 'P54709', 'P56199', 'P57087', 'P59901', 'P60033', 'P61073', 'P78324', 'P78325', 'P78504', 'P78536', 'P78552', 'Q01113', 'Q01151', 'Q01344', 'Q02094', 'Q02161', 'Q02223', 'Q02763', 'Q03405', 'Q04900', 'Q04912', 'Q07011', 'Q07075', 'Q07108', 'Q07954', 'Q08334', 'Q08345', 'Q08708', 'Q08722', 'Q10588', 'Q10589', 'Q12913', 'Q12918', 'Q13241', 'Q13261', 'Q13291', 'Q13349', 'Q13478', 'Q13651', 'Q13740', 'Q14242', 'Q14627', 'Q14773', 'Q14953', 'Q14954', 'Q15116', 'Q15223', 'Q15399', 'Q15762', 'Q16570', 'Q16832', 'Q496F6', 'Q5QGZ9', 'Q5ZPR3', 'Q6GTX8', 'Q6UXB8', 'Q6UXG3', 'Q6UXN8', 'Q6UXZ3', 'Q6YHK3', 'Q7Z6A9', 'Q86VB7', 'Q8IUN9', 'Q8IX05', 'Q8N149', 'Q8N423', 'Q8N6C8', 'Q8N6Q3', 'Q8N743', 'Q8NHJ6', 'Q8NHL6', 'Q8TDQ0', 'Q8TDQ1', 'Q8WTT0', 'Q8WWI5', 'Q8WWV6', 'Q8WXI8', 'Q92692', 'Q92854', 'Q92956', 'Q93033', 'Q93070', 'Q969P0', 'Q96D42', 'Q96DU3', 'Q96F46', 'Q96LA5', 'Q96LA6', 'Q96P31', 'Q96PJ5', 'Q96RD9', 'Q96RJ3', 'Q99062', 'Q99467', 'Q99706', 'Q9BQ51', 'Q9BXN2', 'Q9BXR5', 'Q9BZW8', 'Q9BZZ2', 'Q9H2X3', 'Q9H5V8', 'Q9HBE5', 'Q9HBG7', 'Q9HCU0', 'Q9NNX6', 'Q9NP84', 'Q9NP99', 'Q9NPF0', 'Q9NPY3', 'Q9NQ25', 'Q9NQS3', 'Q9NR16', 'Q9NR96', 'Q9NR97', 'Q9NYZ4', 'Q9NZQ7', 'Q9P0V8', 'Q9P1W8', 'Q9P2B2', 'Q9UBG0', 'Q9UBN6', 'Q9UGN4', 'Q9UHX3', 'Q9UIB8', 'Q9UJ71', 'Q9ULV1', 'Q9ULW2', 'Q9UM73', 'Q9UMR7', 'Q9UNN8', 'Q9UNQ0', 'Q9UQV4', 'Q9Y275', 'Q9Y286', 'Q9Y2C9', 'Q9Y5U5', 'Q9Y5Y4', 'Q9Y624', 'Q9Y6Q6', 'Q9Y6W8', 'P43627', 'Q6ISS4', 'Q8N109', 'Q8NHK3']
okk = {
    "A6NI73": [
        "CD85f"
    ],
    "A8K4G0": [
        "CD300b"
    ],
    "O00144": [
        "CD349"
    ],
    "O00206": [
        "CD284"
    ],
    "O00220": [
        "CD261"
    ],
    "O00238": [
        "CDw293"
    ],
    "O00241": [
        "CD172b"
    ],
    "O00481": [
        "CD277"
    ],
    "O00574": [
        "CD186"
    ],
    "O14638": [
        "CD203c"
    ],
    "O14672": [
        "CD156c"
    ],
    "O14763": [
        "CD262"
    ],
    "O14786": [
        "CD304"
    ],
    "O14788": [
        "CD254"
    ],
    "O14798": [
        "CD263"
    ],
    "O14836": [
        "CD267"
    ],
    "O14931": [
        "CD337"
    ],
    "O15389": [
        "CD170"
    ],
    "O15455": [
        "CD283"
    ],
    "O43490": [
        "CD133"
    ],
    "O43557": [
        "CD258"
    ],
    "O43699": [
        "CD327"
    ],
    "O60449": [
        "CD205"
    ],
    "O60486": [
        "CD232"
    ],
    "O60603": [
        "CD282"
    ],
    "O75015": [
        "CD16b"
    ],
    "O75019": [
        "CD85i"
    ],
    "O75022": [
        "CD85a"
    ],
    "O75023": [
        "CD85c"
    ],
    "O75144": [
        "CD275"
    ],
    "O75326": [
        "CD108"
    ],
    "O75330": [
        "CD168"
    ],
    "O75509": [
        "CD358"
    ],
    "O75888": [
        "CD256"
    ],
    "O76036": [
        "CD335"
    ],
    "O95256": [
        "CD218b"
    ],
    "O95727": [
        "CD355"
    ],
    "O95944": [
        "CD336"
    ],
    "O95971": [
        "CD160"
    ],
    "P01589": [
        "CD25"
    ],
    "P01730": [
        "CD4"
    ],
    "P01732": [
        "CD8a"
    ],
    "P02724": [
        "CD235a"
    ],
    "P02730": [
        "CD233"
    ],
    "P02786": [
        "CD71"
    ],
    "P04156": [
        "CD230"
    ],
    "P04216": [
        "CD90"
    ],
    "P04233": [
        "CD74"
    ],
    "P04234": [
        "CD3d"
    ],
    "P04626": [
        "CD340"
    ],
    "P04921": [
        "CD236"
    ],
    "P05106": [
        "CD61"
    ],
    "P05107": [
        "CD18"
    ],
    "P05362": [
        "CD54"
    ],
    "P05556": [
        "CD29"
    ],
    "P06028": [
        "CD235b"
    ],
    "P06126": [
        "CD1a"
    ],
    "P06127": [
        "CD5"
    ],
    "P06213": [
        "CD220"
    ],
    "P06729": [
        "CD2"
    ],
    "P06731": [
        "CD66e"
    ],
    "P06734": [
        "CD23"
    ],
    "P06756": [
        "CD51"
    ],
    "P07204": [
        "CD141"
    ],
    "P07333": [
        "CD115"
    ],
    "P07359": [
        "CD42b"
    ],
    "P07766": [
        "CD3e"
    ],
    "P08069": [
        "CD221"
    ],
    "P08138": [
        "CD271"
    ],
    "P08174": [
        "CD55"
    ],
    "P08183": [
        "CD243"
    ],
    "P08195": [
        "CD98"
    ],
    "P08473": [
        "CD10"
    ],
    "P08514": [
        "CD41"
    ],
    "P08571": [
        "CD14"
    ],
    "P08575": [
        "CD45"
    ],
    "P08582": [
        "CD228"
    ],
    "P08637": [
        "CD16a"
    ],
    "P08648": [
        "CD49e"
    ],
    "P08887": [
        "CD126"
    ],
    "P08962": [
        "CD63"
    ],
    "P09326": [
        "CD48"
    ],
    "P09564": [
        "CD7"
    ],
    "P09619": [
        "CD140b"
    ],
    "P09693": [
        "CD3g"
    ],
    "P10721": [
        "CD117"
    ],
    "P10747": [
        "CD28"
    ],
    "P10966": [
        "CD8b"
    ],
    "P11049": [
        "CD37"
    ],
    "P11215": [
        "CD11b"
    ],
    "P11279": [
        "CD107a"
    ],
    "P11362": [
        "CD331"
    ],
    "P11464": [
        "CD66f"
    ],
    "P11717": [
        "CD222"
    ],
    "P11836": [
        "CD20"
    ],
    "P11912": [
        "CD79a"
    ],
    "P12018": [
        "CD179a"
    ],
    "P12314": [
        "CD64"
    ],
    "P12318": [
        "CD32"
    ],
    "P12821": [
        "CD143"
    ],
    "P12830": [
        "CD324"
    ],
    "P13164": [
        "CD225"
    ],
    "P13224": [
        "CD42c"
    ],
    "P13473": [
        "CD107b"
    ],
    "P13591": [
        "CD56"
    ],
    "P13598": [
        "CD102"
    ],
    "P13612": [
        "CD49d"
    ],
    "P13688": [
        "CD66a"
    ],
    "P13726": [
        "CD142"
    ],
    "P13987": [
        "CD59"
    ],
    "P14151": [
        "CD62L"
    ],
    "P14209": [
        "CD99"
    ],
    "P14770": [
        "CD42a"
    ],
    "P14778": [
        "CD121a"
    ],
    "P14784": [
        "CD122"
    ],
    "P15144": [
        "CD13"
    ],
    "P15151": [
        "CD155"
    ],
    "P15260": [
        "CD119"
    ],
    "P15391": [
        "CD19"
    ],
    "P15509": [
        "CD116"
    ],
    "P15529": [
        "CD46"
    ],
    "P15812": [
        "CD1e"
    ],
    "P15813": [
        "CD1d"
    ],
    "P15814": [
        "CD179b"
    ],
    "P15941": [
        "CD227"
    ],
    "P16070": [
        "CD44"
    ],
    "P16109": [
        "CD62P"
    ],
    "P16144": [
        "CD104"
    ],
    "P16150": [
        "CD43"
    ],
    "P16234": [
        "CD140a"
    ],
    "P16284": [
        "CD31"
    ],
    "P16410": [
        "CD152"
    ],
    "P16422": [
        "CD326"
    ],
    "P16581": [
        "CD62E"
    ],
    "P16671": [
        "CD36"
    ],
    "P16871": [
        "CD127"
    ],
    "P17301": [
        "CD49b"
    ],
    "P17813": [
        "CD105"
    ],
    "P17927": [
        "CD35"
    ],
    "P18577": [
        "CD240CE"
    ],
    "P18627": [
        "CD223"
    ],
    "P18827": [
        "CD138"
    ],
    "P19022": [
        "CD325"
    ],
    "P19256": [
        "CD58"
    ],
    "P19320": [
        "CD106"
    ],
    "P19397": [
        "CD53"
    ],
    "P19438": [
        "CD120a"
    ],
    "P19440": [
        "CD224"
    ],
    "P20023": [
        "CD21"
    ],
    "P20138": [
        "CD33"
    ],
    "P20273": [
        "CD22"
    ],
    "P20333": [
        "CD120b"
    ],
    "P20701": [
        "CD11a"
    ],
    "P20702": [
        "CD11c"
    ],
    "P20963": [
        "CD247"
    ],
    "P21453": [
        "CD363"
    ],
    "P21589": [
        "CD73"
    ],
    "P21730": [
        "CD88"
    ],
    "P21757": [
        "CD204"
    ],
    "P21802": [
        "CD332"
    ],
    "P21854": [
        "CD72"
    ],
    "P21926": [
        "CD9"
    ],
    "P22413": [
        "CD203a"
    ],
    "P22455": [
        "CD334"
    ],
    "P22607": [
        "CD333"
    ],
    "P22897": [
        "CD206"
    ],
    "P23229": [
        "CD49f"
    ],
    "P23276": [
        "CD238"
    ],
    "P23510": [
        "CD252"
    ],
    "P24071": [
        "CD89"
    ],
    "P24394": [
        "CD124"
    ],
    "P25024": [
        "CD181"
    ],
    "P25025": [
        "CD182"
    ],
    "P25063": [
        "CD24"
    ],
    "P25445": [
        "CD95"
    ],
    "P25942": [
        "CD40"
    ],
    "P26006": [
        "CD49c"
    ],
    "P26715": [
        "CD159a"
    ],
    "P26717": [
        "CD159c"
    ],
    "P26718": [
        "CD314"
    ],
    "P26842": [
        "CD27"
    ],
    "P26951": [
        "CD123"
    ],
    "P27487": [
        "CD26"
    ],
    "P27701": [
        "CD82"
    ],
    "P27930": [
        "CD121b"
    ],
    "P28906": [
        "CD34"
    ],
    "P28907": [
        "CD38"
    ],
    "P28908": [
        "CD30"
    ],
    "P29016": [
        "CD1b"
    ],
    "P29017": [
        "CD1c"
    ],
    "P29965": [
        "CD154"
    ],
    "P30203": [
        "CD6"
    ],
    "P31358": [
        "CD52"
    ],
    "P31785": [
        "CD132"
    ],
    "P31997": [
        "CD66b"
    ],
    "P32004": [
        "CD171"
    ],
    "P32246": [
        "CD191"
    ],
    "P32248": [
        "CD197"
    ],
    "P32302": [
        "CD185"
    ],
    "P32927": [
        "CD131"
    ],
    "P32942": [
        "CD50"
    ],
    "P32970": [
        "CD70"
    ],
    "P32971": [
        "CD153"
    ],
    "P33151": [
        "CD144"
    ],
    "P33681": [
        "CD80"
    ],
    "P34741": [
        "CD362"
    ],
    "P34810": [
        "CD68"
    ],
    "P34910": [
        "CD361"
    ],
    "P35613": [
        "CD147"
    ],
    "P35968": [
        "CD309"
    ],
    "P36888": [
        "CD135"
    ],
    "P36894": [
        "CD292"
    ],
    "P38570": [
        "CD103"
    ],
    "P40189": [
        "CD130"
    ],
    "P40197": [
        "CD42d"
    ],
    "P40198": [
        "CD66d"
    ],
    "P40199": [
        "CD66c"
    ],
    "P40200": [
        "CD96"
    ],
    "P40238": [
        "CD110"
    ],
    "P40259": [
        "CD79b"
    ],
    "P41217": [
        "CD200"
    ],
    "P41597": [
        "CD192"
    ],
    "P41732": [
        "CD231"
    ],
    "P42081": [
        "CD86"
    ],
    "P42701": [
        "CD212"
    ],
    "P42702": [
        "CD118"
    ],
    "P43121": [
        "CD146"
    ],
    "P43489": [
        "CD134"
    ],
    "P43626": [
        "CD158a"
    ],
    "P43628": [
        "CD158b2"
    ],
    "P43629": [
        "CD158e"
    ],
    "P43630": [
        "CD158k"
    ],
    "P43631": [
        "CD158j"
    ],
    "P43632": [
        "CD158i"
    ],
    "P48023": [
        "CD178"
    ],
    "P48357": [
        "CD295"
    ],
    "P48509": [
        "CD151"
    ],
    "P48960": [
        "CD97"
    ],
    "P49682": [
        "CD183"
    ],
    "P49961": [
        "CD39"
    ],
    "P50591": [
        "CD253"
    ],
    "P50895": [
        "CD239"
    ],
    "P51677": [
        "CD193"
    ],
    "P51679": [
        "CD194"
    ],
    "P51681": [
        "CD195"
    ],
    "P51684": [
        "CD196"
    ],
    "P51685": [
        "CDw198"
    ],
    "P51686": [
        "CDw199"
    ],
    "P52961": [
        "CD296"
    ],
    "P54709": [
        "CD298"
    ],
    "P56199": [
        "CD49a"
    ],
    "P57087": [
        "CD322"
    ],
    "P59901": [
        "CD85g"
    ],
    "P60033": [
        "CD81"
    ],
    "P61073": [
        "CD184"
    ],
    "P78324": [
        "CD172a"
    ],
    "P78325": [
        "CD156a"
    ],
    "P78504": [
        "CD339"
    ],
    "P78536": [
        "CD156b"
    ],
    "P78552": [
        "CD213a1"
    ],
    "Q01113": [
        "CD129"
    ],
    "Q01151": [
        "CD83"
    ],
    "Q01344": [
        "CD125"
    ],
    "Q02094": [
        "CD241"
    ],
    "Q02161": [
        "CD240D"
    ],
    "Q02223": [
        "CD269"
    ],
    "Q02763": [
        "CD202b"
    ],
    "Q03405": [
        "CD87"
    ],
    "Q04900": [
        "CD164"
    ],
    "Q04912": [
        "CD136"
    ],
    "Q07011": [
        "CD137"
    ],
    "Q07075": [
        "CD249"
    ],
    "Q07108": [
        "CD69"
    ],
    "Q07954": [
        "CD91"
    ],
    "Q08334": [
        "CDw210b"
    ],
    "Q08345": [
        "CD167a"
    ],
    "Q08708": [
        "CD300c"
    ],
    "Q08722": [
        "CD47"
    ],
    "Q10588": [
        "CD157"
    ],
    "Q10589": [
        "CD317"
    ],
    "Q12913": [
        "CD148"
    ],
    "Q12918": [
        "CD161"
    ],
    "Q13241": [
        "CD94"
    ],
    "Q13261": [
        "CD215"
    ],
    "Q13291": [
        "CD150"
    ],
    "Q13349": [
        "CD11d"
    ],
    "Q13478": [
        "CD218a"
    ],
    "Q13651": [
        "CD210"
    ],
    "Q13740": [
        "CD166"
    ],
    "Q14242": [
        "CD162"
    ],
    "Q14627": [
        "CD213a2"
    ],
    "Q14773": [
        "CD242"
    ],
    "Q14953": [
        "CD158g"
    ],
    "Q14954": [
        "CD158h"
    ],
    "Q15116": [
        "CD279"
    ],
    "Q15223": [
        "CD111"
    ],
    "Q15399": [
        "CD281"
    ],
    "Q15762": [
        "CD226"
    ],
    "Q16570": [
        "CD234"
    ],
    "Q16832": [
        "CD167b"
    ],
    "Q496F6": [
        "CD300e"
    ],
    "Q5QGZ9": [
        "CD371"
    ],
    "Q5ZPR3": [
        "CD276"
    ],
    "Q6GTX8": [
        "CD305"
    ],
    "Q6UXB8": [
        "CD364"
    ],
    "Q6UXG3": [
        "CD300g"
    ],
    "Q6UXN8": [
        "CD370"
    ],
    "Q6UXZ3": [
        "CD300d"
    ],
    "Q6YHK3": [
        "CD109"
    ],
    "Q7Z6A9": [
        "CD272"
    ],
    "Q86VB7": [
        "CD163"
    ],
    "Q8IUN9": [
        "CD301"
    ],
    "Q8IX05": [
        "CD302"
    ],
    "Q8N149": [
        "CD85h"
    ],
    "Q8N423": [
        "CD85d"
    ],
    "Q8N6C8": [
        "CD85e"
    ],
    "Q8N6Q3": [
        "CD177"
    ],
    "Q8N743": [
        "CD158z"
    ],
    "Q8NHJ6": [
        "CD85k"
    ],
    "Q8NHL6": [
        "CD85j"
    ],
    "Q8TDQ0": [
        "CD366"
    ],
    "Q8TDQ1": [
        "CD300f"
    ],
    "Q8WTT0": [
        "CD303"
    ],
    "Q8WWI5": [
        "CD92"
    ],
    "Q8WWV6": [
        "CD351"
    ],
    "Q8WXI8": [
        "CD368"
    ],
    "Q92692": [
        "CD112"
    ],
    "Q92854": [
        "CD100"
    ],
    "Q92956": [
        "CD270"
    ],
    "Q93033": [
        "CD101"
    ],
    "Q93070": [
        "CD297"
    ],
    "Q969P0": [
        "CD316"
    ],
    "Q96D42": [
        "CD365"
    ],
    "Q96DU3": [
        "CD352"
    ],
    "Q96F46": [
        "CD217"
    ],
    "Q96LA5": [
        "CD307b"
    ],
    "Q96LA6": [
        "CD307a"
    ],
    "Q96P31": [
        "CD307c"
    ],
    "Q96PJ5": [
        "CD307d"
    ],
    "Q96RD9": [
        "CD307e"
    ],
    "Q96RJ3": [
        "CD268"
    ],
    "Q99062": [
        "CD114"
    ],
    "Q99467": [
        "CD180"
    ],
    "Q99706": [
        "CD158d"
    ],
    "Q9BQ51": [
        "CD273"
    ],
    "Q9BXN2": [
        "CD369"
    ],
    "Q9BXR5": [
        "CD290"
    ],
    "Q9BZW8": [
        "CD244"
    ],
    "Q9BZZ2": [
        "CD169"
    ],
    "Q9H2X3": [
        "CD299"
    ],
    "Q9H5V8": [
        "CD318"
    ],
    "Q9HBE5": [
        "CD360"
    ],
    "Q9HBG7": [
        "CD229"
    ],
    "Q9HCU0": [
        "CD248"
    ],
    "Q9NNX6": [
        "CD209"
    ],
    "Q9NP84": [
        "CD266"
    ],
    "Q9NP99": [
        "CD354"
    ],
    "Q9NPF0": [
        "CD320"
    ],
    "Q9NPY3": [
        "CD93"
    ],
    "Q9NQ25": [
        "CD319"
    ],
    "Q9NQS3": [
        "CD113"
    ],
    "Q9NR16": [
        "CD163b"
    ],
    "Q9NR96": [
        "CD289"
    ],
    "Q9NR97": [
        "CD288"
    ],
    "Q9NYZ4": [
        "CD329"
    ],
    "Q9NZQ7": [
        "CD274"
    ],
    "Q9P0V8": [
        "CD353"
    ],
    "Q9P1W8": [
        "CD172g"
    ],
    "Q9P2B2": [
        "CD315"
    ],
    "Q9UBG0": [
        "CD280"
    ],
    "Q9UBN6": [
        "CD264"
    ],
    "Q9UGN4": [
        "CD300a"
    ],
    "Q9UHX3": [
        "CD312"
    ],
    "Q9UIB8": [
        "CD84"
    ],
    "Q9UJ71": [
        "CD207"
    ],
    "Q9ULV1": [
        "CD344"
    ],
    "Q9ULW2": [
        "CD350"
    ],
    "Q9UM73": [
        "CD246"
    ],
    "Q9UMR7": [
        "CD367"
    ],
    "Q9UNN8": [
        "CD201"
    ],
    "Q9UNQ0": [
        "CD338"
    ],
    "Q9UQV4": [
        "CD208"
    ],
    "Q9Y275": [
        "CD257"
    ],
    "Q9Y286": [
        "CD328"
    ],
    "Q9Y2C9": [
        "CD286"
    ],
    "Q9Y5U5": [
        "CD357"
    ],
    "Q9Y5Y4": [
        "CD294"
    ],
    "Q9Y624": [
        "CD321"
    ],
    "Q9Y6Q6": [
        "CD265"
    ],
    "Q9Y6W8": [
        "CD278"
    ],
    "P43627": [
        "CD158b1"
    ],
    "Q6ISS4": [
        "CD306"
    ],
    "Q8N109": [
        "CD158f1"
    ],
    "Q8NHK3": [
        "CD158f2"
    ]
}
missing_ids = [id for id in listing if id not in okk]

# Printing missing IDs
print("Missing IDs:")
print(missing_ids)

import json
import re

with open("uniprot_proteins.txt", "r") as f:
    proteins = f.read().splitlines()

print("Total Uniprot proteins: ", len(proteins))

my_count = 0
my_proteins_uniprot = []
with open("hox_proteins_data.txt") as f:
    entire_file = f.read()
    for protein in proteins:
        if protein in entire_file:
            my_count += 1
            my_proteins_uniprot.append(protein)

print("Total Uniprot proteins identified by me: ", my_count)
print(my_proteins_uniprot)

#Allergens Extraction
# uniprot_ids = ['O43290', 'P02538', 'Q13765', 'Q9BPX6', 'Q9BQE9']
# allergens = ["Hom s 1", "Hom s 5", "Hom s 2", "Hom s 4", "Hom s 3"]
# allergen_dict = dict(zip(uniprot_ids, allergens))
#
# print(allergen_dict)
# with open("allergen_dict.json", "w") as f:
#     json.dump(allergen_dict, f, indent=4)

#Glycogen Extraction
# Extracting corresponding glycosid for each identified UniProt protein
# protein_glycosid_mapping = {}
# with open("glycosid_data.txt") as f:
#     for line in f:
#         parts = line.split("(")
#         if len(parts) > 1:
#             glycosid = parts[0].strip()
#             protein_id = parts[1].split(")")[0].strip()
#             protein_glycosid_mapping[protein_id] = glycosid
#
# # Displaying glycogen ID for each identified UniProt protein
# uniprot_ids = []
# glycosid_data = []
# for protein in my_proteins_uniprot:
#     if protein in protein_glycosid_mapping:
#         glycogen_id = protein_glycosid_mapping[protein]
#         uniprot_ids.append(protein)
#         glycosid_data.append(glycogen_id)
#         print(f"Protein ID: {protein}, Glycogen ID: {glycogen_id}")
#     else:
#         print(f"No glycogen ID found for protein {protein}")
#
# glycosid_dict = dict(zip(uniprot_ids, glycosid_data))
# with open("glycosid_dict.json", "w") as f:
#     json.dump(glycosid_dict, f, indent=4)

#Blood Antigen Proteins
# Open the file and read its content
# with open("blood_anntigen_data.txt", "r") as f:
#     content = f.read()
#
# # Initialize a dictionary to store the antigen values for each UniProt protein ID
# antigen_values = {}
#
# # Iterate through each UniProt protein ID
# for protein_id in my_proteins_uniprot:
#     # Use regular expression to find the antigen values associated with the protein ID
#     pattern = rf"{protein_id}.*(?:\n\s*Antigens:\s*([\s\S]*?))(?=\n{{2,}})"
#     matches = re.findall(pattern, content)
#
#     # Store the antigen values in the dictionary
#     antigen_values[protein_id] = matches
#
# uniprot_ids = []
# antigens_data = []
# # Print the antigen values corresponding to each UniProt protein ID
# for protein_id, antigens in antigen_values.items():
#     if antigens:
#         uniprot_ids.append(protein_id)
#         antigens_data.append(antigens)
#     else:
#         continue
# antigens_dict = dict(zip(uniprot_ids, antigens_data))
# with open("antigens_dict.json", "w") as f:
#     json.dump(antigens_dict, f, indent=4)

#Kinases Extraction
# # Open the file and read its content
# with open("kinases_data.txt", "r") as f:
#     content = f.read()
#
# kinases_values = {}
#
# # Iterate through each UniProt protein ID
# for protein_id in my_proteins_uniprot:
#     # Use regular expression to find the antigen values associated with the protein ID
#     pattern = fr"\b(\w+)\s+(\w+_HUMAN)\s+\({protein_id}\s*\)"
#     matches = re.findall(pattern, content)
#
#     # Store the antigen values in the dictionary
#     kinases_values[protein_id] = matches[0][0]
#
# uniprot_ids = []
# kinases_data = []
# # Print the antigen values corresponding to each UniProt protein ID
# for protein_id, antigens in kinases_values.items():
#     if antigens:
#         uniprot_ids.append(protein_id)
#         kinases_data.append(antigens)
#     else:
#         continue
# kinases_dict = dict(zip(uniprot_ids, kinases_data))
# with open("kinases_dict.json", "w") as f:
#     json.dump(kinases_dict, f, indent=4)

#CD Extraction
# Open the file and read its content
# with open("CD_data.txt", "r") as f:
#     content = f.read()
#
# CD_values = {}
#
# # Iterate through each UniProt protein ID
# for protein_id in my_proteins_uniprot:
#     # Use regular expression to find the antigen values associated with the protein ID
#     pattern = fr"\b(\w+)\s+(\w+_HUMAN)\s+{protein_id}\s+"
#     matches = re.findall(pattern, content)
#
#     # Store the antigen values in the dictionary
#     CD_values[protein_id] = [match[0] for match in matches]
#
#
# uniprot_ids = []
# CD_data = []
# # Print the antigen values corresponding to each UniProt protein ID
# for protein_id, antigens in CD_values.items():
#     if antigens:
#         uniprot_ids.append(protein_id)
#         CD_data.append(antigens)
#     else:
#         continue
# CD_dict = dict(zip(uniprot_ids, CD_data))
# with open("CD_dict.json", "w") as f:
#     json.dump(CD_dict, f, indent=4)

#MT Extraction
# # Open the file and read its content
# with open("MT_data.txt", "r") as f:
#     content = f.read()
#
# MT_values = {}
#
# # Iterate through each UniProt protein ID
# for protein_id in my_proteins_uniprot:
#     # Use regular expression to find the antigen values associated with the protein ID
#     pattern = rf"(\w+)_\w+\s+\({re.escape(protein_id)}\)"
#     matches = re.findall(pattern, content)
#     print(matches)
#     # Store the antigen values in the dictionary
#     MT_values[protein_id] = matches[0]
#
#
# uniprot_ids = []
# MT_data = []
# # Print the antigen values corresponding to each UniProt protein ID
# for protein_id, antigens in MT_values.items():
#     if antigens:
#         uniprot_ids.append(protein_id)
#         MT_data.append(antigens)
#     else:
#         continue
# MT_dict = dict(zip(uniprot_ids, MT_data))
# with open("MT_dict.json", "w") as f:
#     json.dump(MT_dict, f, indent=4)

#peptides_extraction
# # Open the file and read its content
# with open("peptides_data.txt", "r") as f:
#     content = f.read()
#
# peptides_values = {}
#
# # Iterate through each UniProt protein ID
# for protein_id in my_proteins_uniprot:
#     # Use regular expression to find the antigen values associated with the protein ID
#     pattern = rf"(\w+)_\w+\s+\({re.escape(protein_id)}\)"
#     matches = re.findall(pattern, content)
#     print(matches)
#     # Store the antigen values in the dictionary
#     peptides_values[protein_id] = matches[0]
#
#
# uniprot_ids = []
# peptides_data = []
# # Print the antigen values corresponding to each UniProt protein ID
# for protein_id, antigens in peptides_values.items():
#     if antigens:
#         uniprot_ids.append(protein_id)
#         peptides_data.append(antigens)
#     else:
#         continue
# peptides_dict = dict(zip(uniprot_ids, peptides_data))
# with open("peptides_dict.json", "w") as f:
#     json.dump(peptides_dict, f, indent=4)

#Restriction_Enzyme_Extraction ------ 0 uniprot_proteins_identified

#Ribosomal_Proteins_Extraction
# # Open the file and read its content
# with open("ribosomal_proteins_data.txt", "r") as f:
#     content = f.read()
#
# peptides_values = {}
#
# # Iterate through each UniProt protein ID
# for protein_id in my_proteins_uniprot:
#     # Use regular expression to find the antigen values associated with the protein ID
#     pattern = rf"(\w+)_\w+\s+\({re.escape(protein_id)}\)"
#     matches = re.findall(pattern, content)
#     print(matches)
#     # Store the antigen values in the dictionary
#     peptides_values[protein_id] = matches[0]
#
#
# uniprot_ids = []
# peptides_data = []
# # Print the antigen values corresponding to each UniProt protein ID
# for protein_id, antigens in peptides_values.items():
#     if antigens:
#         uniprot_ids.append(protein_id)
#         peptides_data.append(antigens)
#     else:
#         continue
# peptides_dict = dict(zip(uniprot_ids, peptides_data))
# with open("ribosomal_proteins_dict.json", "w") as f:
#     json.dump(peptides_dict, f, indent=4)

#toxins_data ------- 0 uniprot proteins unidentified

#hox_proteins extraction
# Open the file and read its content
with open("hox_proteins_data.txt", "r") as f:
    content = f.read()

peptides_values = {}

# Iterate through each UniProt protein ID
for protein_id in my_proteins_uniprot:
    # Use regular expression to find the antigen values associated with the protein ID
    pattern = rf"{re.escape(protein_id)};\s+(\w+_HUMAN)"
    matches = re.findall(pattern, content)
    print(matches)
    # Store the antigen values in the dictionary
    peptides_values[protein_id] = matches[0]


uniprot_ids = []
peptides_data = []
# Print the antigen values corresponding to each UniProt protein ID
for protein_id, antigens in peptides_values.items():
    if antigens:
        uniprot_ids.append(protein_id)
        peptides_data.append(antigens)
    else:
        continue
peptides_dict = dict(zip(uniprot_ids, peptides_data))
with open("hox_proteins_dict.json", "w") as f:
    json.dump(peptides_dict, f, indent=4)

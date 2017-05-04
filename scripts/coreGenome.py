ID_H = ['JF72','F259','JG30','JF76','F260','F261','F262','F263',
        'JF73','JF74','JF75','WB8','WB10','L185','L186',
        'L184','L183'] # Honeybee-bacterial genomes
ID_B = ['F225','F230','F233','F234','F236','F237',
        'F228','F245','F246','F247'] # Bumblebee-bacterial genomes
ID_O = ['LA14','LA2','LDB','LGAS','LHV','LJP','WANG','JG29'] # Outgroup bacterial genomes

# Loading ortholog table and creating output file for core genes (= genes present in all strain genomes) 
ortho_tab = open("mclOutput",'r')
Core_genes = open("Core_genes.txt",'w')
#==========================================
allStrainsIDs = ID_H + ID_B + ID_O # Get all strain IDs in one list by concatenating the 3 lists above
count=0   # Variable to store number of lines in which all strains are present
for line in ortho_tab:  # Go through each line of the orthoMCL output file 
    allStrainsOnLine=False # if all strains are on a line, will be set to True. Reset to False for each line of ortho_tab (= each gene family)
    for element in allStrainsIDs: # Go through all strain IDs
        if element not in line:   # If a strain is absent from current line
            allStrainsOnLine = False # This line doesn't have all strains
            break                    # Exit from this ortho_line (goes to next line, no need to waste time to test the presence of other strains at this point)
        else:                        # Otherwise
            allStrainsOnLine = True  # set all StrainsOnLine = True
    if allStrainsOnLine:             # if allStrainsOnline (is True),
                                     # IMPORTANT: this condition will be evaluated only if ALL strains where present
                                     # (if one is absent, the break above will preclude this line from being executed)
        count+=1                     # Count this line as a line in which all strains are present
        Core_genes.write(line)       # Store this line in file
print(count)                         # Print the number of Core genes (determined as number of lines from orthoMCL output file containing ALL strains)
Core_genes.close()                   # Close file once you no longer write anything in it.





        

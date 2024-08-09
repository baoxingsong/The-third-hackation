import re
sequence = "ACTGNNNNNNNNNNNNNNNNNNNNACGTGCAACTGNNNNNNNNNNNNNNNNNNNNACGTGCAACTGNNNNNNNNNNNNNNNNNNNNACGTGCAACTGNNNNNNNNNNNNNNNNNNNNACGTGCA"
largestGap = 5
pattern = "N{" + str(largestGap) + ",}"
print(pattern)
sequence = re.sub(pattern, "N"*largestGap, sequence)
print(sequence)

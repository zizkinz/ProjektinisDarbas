import pandas as pd
import matplotlib.pyplot as plt


# df = pd.read_csv("test_data/DP23_sequence.bed")
# L = ["Geeks\n", "for\n", "Geeks\n"]
# Glul
# Location   chr1:153,775,692-153,785,469 (+)
# Using readline()
file1 = open("test_data/DP23_sequence.bed", 'r')
file2 = open("test_data/DP23_sequence_Glul.bed", 'w')
count = 0
# list = []

while True:
    count += 1

    # Get next line from file
    line = file1.readline()
    z = line.strip().split(sep = "	")
    z[0] = int(z[0])
    z[1] = int(z[1])
    z[2] = int(z[2])
    if z[1] > 153775692 and z[1] < 153785469:
        file2.writelines(line)
        print(z)
    else:
        continue
    # if line is empty
    # end of file is reached
    if count == 1e6 or not line:
        break
    # print("Line{}: {}".format(count, line.strip()))
    # list.append(int(z[2]) - int(z[1]))


file1.close()
#

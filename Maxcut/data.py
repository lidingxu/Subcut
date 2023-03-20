import numpy as np

benchmark = "benchmark"

ns = [20, 30, 40]


for n in ns:
    ws = {}
    m = 0
    for i in range(n):
        for j in range(i+1,n):
            ws[i,j] = [i, j, np.random.uniform(0, 1)]
            m += 1

    descriptor = str(n) + " " + str(m)
    file_name = str(n) + "_" + str(m) + ".cut"
    f = open(benchmark + "/" + file_name, "x")
    f.write(descriptor +  '\n')

    for i in range(n):
        for j in range(i+1,n):
            f.write(str(ws[i,j][0]) + " " + str(ws[i,j][1]) + " " + str(ws[i,j][2]) + "\n")
    f.close()


            



